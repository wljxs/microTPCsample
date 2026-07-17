#include "Script/DeconvAnalysisScript.h"
#include "WaveformAnalysis.h"
#include "IAlgorithm.h"
#include "Script/Base/ScriptFactory.h"

#include "AnalysisUtils.h"
#include "Detector/DetectorFactory.h"
#include "Detector/Planar.h"
#include "Event/DataModel.h"
#include "Event/DetectorFrame.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "RtypesCore.h"
#include "Script/Base/RawDataParser.h"
#include "Script/Base/ScriptManager.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TMarker.h"
#include "TGraphErrors.h"
#include "TInterpreter.h"
#include "TMultiGraph.h"
#include "TVector3.h"
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TH1.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <new>
#include <regex>
#include <string>
#include <vector>

double pitchwidth[2] = {0.4, 0.4};
int n_tier = 5;
int voltage = 500;
int drawmaxeventid = 1000; // 为了画图对齐，限制只画前1000个事件的波形和电荷分布图

// 残差协方差逆矩阵（正则化），用于协方差对齐（在 RunDeconvolution 中 n_tier 同步后初始化）
static Eigen::MatrixXd s_tier_cov_inv[2];

double thickness = 5.0;    // mm，增益层厚度，用于计算电子漂移时间
double vdrift = 0.04;      // mm/ns，电子漂移速度，用于计算电子漂移时间
double drift_ratio = 1.0;          // 电子漂移速度与增益层厚度的比值，用于计算电子漂移时间


double z_ref_position = thickness / 5.0*2.5; // mm，参考位置，通常取增益层中间位置，作为对齐的基准点

constexpr int apvsamples = 27;
double delayindex2d[2] = {97/25.0,88/25.0};
// double delayindex2d[2] = {3.,3.};
int delayindex = 0; // 从画出来的平均波形看出来的粗对齐，用于补偿 t0 延迟
int fit_point = std::min(22,
                         apvsamples - 3); // 成形函数拟合点数，不能太多，否则矩阵太大不好求逆，不能太少，否则成形函数不够精细
using namespace std;

void RunDeConvAlign(const std::vector<Event> &events, std::shared_ptr<Detector> detector, int detID, bool is3D, bool useposition = false, double reslimit = 4.0);
void RunDeConvAlignCov(const std::vector<Event> &events, std::shared_ptr<Detector> detector, int detID, bool is3D);
double DeConvChi2Cov_3D(const double *par, const std::vector<Event> &events, std::shared_ptr<Detector> detector, int detID);//打算用反卷积进行对齐的

// 更新 delay.root 缓存（定义在文件末尾）
static void UpdateDelayRootCache(const std::string& cacheFile,
                                  const double fitted_delay[2][256],
                                  const double t0_per_strip[2][256]);
double DeConvChi2Cov_6D(const double *par, const std::vector<Event> &events, std::shared_ptr<Detector> detector, int detID);
double DeConvTrackChi2_3D(const double *par, const std::vector<Event> &events, std::shared_ptr<Detector> detector, int detID,double reslimit);//用microTPC拟合的track进行对齐
double DeConvTrackChi2_6D(const double *par, const std::vector<Event> &events, std::shared_ptr<Detector> detector, int detID,double reslimit);
double DeConvPositionChi2_3D(const double *par, const std::vector<Event> &events, std::shared_ptr<Detector> detector, int detID,double reslimit);//用microTPC拟合的position进行对齐
double DeConvPositionChi2_6D(const double *par, const std::vector<Event> &events, std::shared_ptr<Detector> detector, int detID,double reslimit);
std::map<int,std::pair<double,double>> Positionanalysis(std::vector<Event> events,std::string mode,TDirectory* fout);

void DeconvAnalysisScript::LoadConfig(const json &config) {
  m_message = config.value("message", "Hello from DeconvAnalysisScript!");
  m_eventLimit = config.value("eventLimit", 10);
  m_mode = config.value("mode", "full");//"deconv" / "align" / "full"
  m_usePerChannelT0 = config.value("usePerChannelT0", false);
  m_detectBadStrips = config.value("detectBadStrips", true);
  m_cuts.LoadFromJson(config);
}

void DeconvAnalysisScript::Print() const {
  std::cout << "  DeconvAnalysisScript Configuration:" << std::endl;
  std::cout << "  Message: " << m_message << std::endl;
  std::cout << "  Event Limit: " << m_eventLimit << std::endl;
  std::cout << "  Mode: " << m_mode << std::endl;
  std::cout << "  Use Per-Channel T0: " << (m_usePerChannelT0 ? "yes" : "no") << std::endl;
  std::cout << "  Detect Bad Strips: " << (m_detectBadStrips ? "yes" : "no") << std::endl;
  m_cuts.Print();
}

bool DeconvAnalysisScript::Execute() {
  std::cout << "\n========================================" << std::endl;
  std::cout << " Deconv Analysis Script (mode: " << m_mode << ")" << std::endl;
  std::cout << "========================================" << std::endl;

  if (m_mode == "deconv" || m_mode == "full") {
    if (!RunDeconvolution()) return false;
  }
  if (m_mode == "align" || m_mode == "full") {
    if (!RunAlignment()) return false;
  }
  return true;
}

bool DeconvAnalysisScript::RunDeconvolution() {
  gROOT->SetBatch(kTRUE);    // 关闭图形界面，避免弹出画布
  TH1::AddDirectory(kFALSE); // 直方图不自动添加到当前目录，避免内存泄漏
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);

  auto tStart = chrono::high_resolution_clock::now();
  int m_progressInterval = 1000;
  int m_maxEvents = 60000;
  std::cout << m_maxEvents << " events will be processed at most." << std::endl;
  std::cout << "\n========================================" << std::endl;
  std::cout << " Deconv Analysis Script" << std::endl;
  std::cout << "========================================" << std::endl;

  // ---- 将 m_cuts 中的可配置参数同步到全局变量 ----
  pitchwidth[0] = m_cuts.pitchwidth[0];  pitchwidth[1] = m_cuts.pitchwidth[1];
  n_tier        = m_cuts.n_tier;
  voltage       = m_cuts.voltage;
  thickness     = m_cuts.thickness;
  vdrift        = m_cuts.vdrift;
  drift_ratio   = m_cuts.drift_ratio;
  z_ref_position = m_cuts.z_ref_position;
  delayindex2d[0] = m_cuts.delayindex2d[0];  delayindex2d[1] = m_cuts.delayindex2d[1];
  delayindex    = m_cuts.delayindex;
  fit_point     = m_cuts.fit_point;

  // 在 n_tier 同步后初始化协方差逆矩阵
  s_tier_cov_inv[0] = Eigen::MatrixXd::Identity(n_tier, n_tier);
  s_tier_cov_inv[1] = Eigen::MatrixXd::Identity(n_tier, n_tier);
  // -------------------------------------------------

  auto parser = GetParser();
  if (!parser) {
    cerr << "Error: Parser not set!" << endl;
    return false;
  }



  // ====== Phase 0: 加载track信息 ======
  string trackFile = GetOutputDir() + "TrackInfo.root";
  std::cout << "\nLoading track info..." << std::endl;
  std::cout << "File: " << trackFile << std::endl;

  TFile *f = TFile::Open(trackFile.c_str(), "READ");
  if (!f || f->IsZombie()) {
    cerr << "Error: Cannot open " << trackFile << endl;
    cerr << "Please run Track Analysis first!" << endl;
    return false;
  }

  TTree *trackTree = (TTree *)f->Get("Tracks");
  if (!trackTree) {
    cerr << "Error: No Tracks tree!" << endl;
    f->Close();
    return false;
  }

  Int_t eventID;
  Track *track = nullptr;
  double sigTime;
  trackTree->SetBranchAddress("eventID", &eventID);
  trackTree->SetBranchAddress("track", &track);
  trackTree->SetBranchAddress("t0", &sigTime);

  cout << "\nProcessing Deconv data..." << endl;

  Long64_t nEntries = trackTree->GetEntries();
  if (m_maxEvents > 0 && nEntries > m_maxEvents) {
    nEntries = m_maxEvents;
  }
  std::cout << nEntries << " entries to process." << std::endl;

  // ====== Phase 1: 导入pedestal std ======
  TH1D *h_pedestal[2];
  for(int type=0;type<2;type++){
    h_pedestal[type] = new TH1D(Form("h_pedestal_type%d", type), Form("Pedestal Distribution Type %d;strip number;Pedestal STD", type), 256, 1, 257);
  }

  auto &factory = DetectorFactory::GetInstance();
  const auto &deconvs = factory.GetDetectorsByRole(Detector::Role::DeConv);
  std::cout << deconvs.size() << " DeConv detectors found." << std::endl;
  std::map<int, std::shared_ptr<const std::map<int, std::vector<float>>>> pedStdByDet;
  for (const auto &det : deconvs) {
    const int detID = det->GetID();
    auto pedByType = std::make_shared<std::map<int, std::vector<float>>>();
    for (int type : det->getConfig().readoutPlaneType) {
      if(type<0||type>=2){
        std::cout << "Warning: Invalid readout plane type " << type << " for detector " << detID << ", skipping pedestal loading for this type." << std::endl;
        continue;
      }
      (*pedByType)[type] = parser->GetPedstd(detID, type);
      for (int i = 0; i < (*pedByType)[type].size(); ++i) {
        // std::cout << "Det " << detID << " Type " << type << " Strip " << (i + 1) << " PedStd: " << (*pedByType)[type][i] << std::endl;
          h_pedestal[type]->Fill(i,(*pedByType)[type][i]);
        if (type == 1 && i < 128) {
          (*pedByType)[type][i] = 12;
        }
      }
    }
    pedStdByDet[detID] = pedByType;
  }
  
  // ====== Phase 2: 加入模拟的波形 ======
  TFile *fmean_sample;
  if(voltage==250){
  fmean_sample= TFile::Open("../raw/gain_and_mean_sample250.root", "READ");}
  else{
  fmean_sample= TFile::Open("../raw/gain_and_mean_sample.root", "READ");
  }

  if (!fmean_sample || fmean_sample->IsZombie()) {
    cerr << "Error: Cannot open mean sample file" << endl;
    return false;
  }
  TTree *tree_mean_sample = (TTree *)fmean_sample->Get("tree_mean_sample");
  const char *treeName[2] = {"tree_x", "tree_y"};
  TTree *tree_sample[2] = {nullptr, nullptr};
  std::vector<double> *sigzero[2] = {nullptr, nullptr};
  std::vector<double> *sigtrackmean[2] = {nullptr, nullptr};
  std::vector<double> *sigtrackmean_conv[2] = {nullptr, nullptr};
  TGraph *gr_trackmean_tran[2] = {nullptr, nullptr};
  gr_trackmean_tran[0] = (TGraph *)fmean_sample->Get("gr_trackmean_tran_d0");
  gr_trackmean_tran[1] = (TGraph *)fmean_sample->Get("gr_trackmean_tran_d1");
  TGraph *gr_trackmean_conv_tran[2] = {nullptr, nullptr};
  gr_trackmean_conv_tran[0] = (TGraph *)fmean_sample->Get("gr_trackmean_conv_tran_d0");
  gr_trackmean_conv_tran[1] = (TGraph *)fmean_sample->Get("gr_trackmean_conv_tran_d1");
  if (!gr_trackmean_tran[0] || !gr_trackmean_tran[1] || !gr_trackmean_conv_tran[0] || !gr_trackmean_conv_tran[1]) {
    cerr << "Error: Cannot find gr_trackmean_tran or gr_trackmean_conv_tran in mean sample file" << endl;
  }
  for (int d = 0; d < 2; ++d) {
    tree_sample[d] = (TTree *)fmean_sample->Get(treeName[d]);
    tree_sample[d]->SetBranchAddress("sigzero", &sigzero[d]);
    tree_sample[d]->SetBranchAddress("sigtrackmean", &sigtrackmean[d]);
    tree_sample[d]->SetBranchAddress("sigtrackmean_conv", &sigtrackmean_conv[d]);
    tree_sample[d]->GetEntry(0);
  }

  // sigzero_sample[0]=x信号模板, [1]=y信号模板
  std::vector<double> sigzero_sample[2];
  sigzero_sample[0] = std::vector<double>(sigzero[0]->begin(), sigzero[0]->end());
  sigzero_sample[1] = std::vector<double>(sigzero[1]->begin(), sigzero[1]->end());
  // 拟合用模板：不加延迟；下面的 b_vec 已去除延迟以保持一致。
  // 分 x/y 两个反卷积矩阵
  Eigen::MatrixXd deconvMatrix[2];
  deconvMatrix[0] = AnalysisUtils::DeconvBaseMatrix(sigzero_sample[0], fit_point, n_tier);
  deconvMatrix[1] = AnalysisUtils::DeconvBaseMatrix(sigzero_sample[1], fit_point, n_tier);
  if (deconvMatrix[0].size() == 0 || deconvMatrix[1].size() == 0) {
    cerr << "Error: Failed to initialize deconvolution base matrix!" << endl;
    return false;
  }

  // ====== Phase 3: 检查/生成 per-channel t0 ROOT 缓存文件,导入t0 ======
  string t0cacheFile = GetOutputDir() + "delay.root";
  double t0_per_strip[2][256];
  for (int type = 0; type < 2; type++) {
    for (int strip = 0; strip < 256; strip++) {
      t0_per_strip[type][strip] = delayindex2d[type] * 25.0; // 默认用粗延迟
    }
  }

  if (m_usePerChannelT0) {
    TFile *ft0cache = TFile::Open(t0cacheFile.c_str(), "READ");
    if (!ft0cache || ft0cache->IsZombie()) {
      if (ft0cache) { ft0cache->Close(); delete ft0cache; }
      std::cout << "Phase 0: No t0 cache found, creating " << t0cacheFile << " with defaults..." << std::endl;
      ft0cache = TFile::Open(t0cacheFile.c_str(), "RECREATE");
      const char *t0names[2] = {"t0_strip_arr_x", "t0_strip_arr_y"};
      for (int type = 0; type < 2; type++) {
        TH1D h_t0(t0names[type], t0names[type], 256, 0, 256);
        for (int strip = 0; strip < 256; strip++)
          h_t0.SetBinContent(strip + 1, t0_per_strip[type][strip]);
        h_t0.Write();
      }
      ft0cache->Close();
      delete ft0cache;
      std::cout << "Created " << t0cacheFile << std::endl;
    } else {
      std::cout << "Phase 0: Loading per-channel t0 from " << t0cacheFile << std::endl;
      const char *t0names[2] = {"t0_strip_arr_x", "t0_strip_arr_y"};
      for (int type = 0; type < 2; type++) {
        TH1D *h_t0 = (TH1D *)ft0cache->Get(t0names[type]);
        if (h_t0 && h_t0->GetNbinsX() >= 256) {
          for (int strip = 0; strip < 256; strip++)
            t0_per_strip[type][strip] = h_t0->GetBinContent(strip + 1);
        } else {
          std::cout << "Warning: " << t0names[type] << " not found, using defaults for " << (type==0?"X":"Y") << std::endl;
        }
      }
      std::cout << "Loaded, e.g. X[50]=" << t0_per_strip[0][50]
                << " ns, Y[50]=" << t0_per_strip[1][50] << " ns" << std::endl;
      ft0cache->Close();
      delete ft0cache;
    }
  } else {
    std::cout << "Phase 0: Per-channel t0 disabled, using fixed default delays." << std::endl;
  }
  std::cout << "Per-channel t0 ready. Starting main deconvolution pass..." << std::endl;



  // ====== Phase 4: 建立root文件，保存相关信息 ======
  string deconvFile = GetOutputDir() + "DeConvInfo.root";
  TFile *fDeconv = new TFile(deconvFile.c_str(), "RECREATE");
  fDeconv->cd();
  // strip和cluster的信息
  TTree *tDeconv = new TTree("DeconvTree", "Deconv data");
  TTree *tDeconvori[2]; // 分开记录好坏道事件，0表示正常事件，1表示筛选后的事件
  Int_t dutID;
  vector<StripHitDeconv> StripHitDeconvs[2]; // 同上，分开记录过没过阈事件
  vector<ClusterDeconv> ClusterDeconvs[2]; // 同上，分开记录过没过阈事件
  int clusternum[2][2]; // 记录每个事件的cluster数量，第一维表示好坏道，第二维表示X/Y向
  for (int i = 0; i < 2; i++) {
    tDeconvori[i] = new TTree(Form("DeconvTree_%s", i == 0 ? "Normal" : "Good"), Form("%s events Deconv data", i == 0 ? "Normal" : "Good"));
    tDeconvori[i]->Branch("eventID", &eventID);
    tDeconvori[i]->Branch("dutID", &dutID);
    tDeconvori[i]->Branch("StripHitDeconvs", &StripHitDeconvs[i]);
    tDeconvori[i]->Branch("ClusterDeconvs", &ClusterDeconvs[i]);
    tDeconvori[i]->Branch("track", &track);
    for(int j=0;j<2;j++){
      tDeconvori[i]->Branch(Form("clusternum_%s", j == 0 ? "X" : "Y"), &clusternum[i][j]);
    }
  }
  
  // 统计分析后的信息，一方面区分两个方向，另一方面来看经过没经过筛选的
  TDirectory *dir[2]; // 用来记录经没经过筛选的的事例
  // 加载事件并处理DUT数据,前面那个表示好不好，后面那个表示X/Y向
  TDirectory *dirwaveform[2][2];
  TDirectory *dirwaveform2d[2][2];      // 用来记前1000个事件的X/Y向条号vs采样点的二维图
  TDirectory *dircharge2d[2][2];        // 用来记前1000个事件的X/Y向条号vs电荷层数的二维图
  TDirectory *dirapvdistribution[2][2]; // 用来记前1000个事件的X/Y向每个apv点的电荷分布图
  TDirectory *dirnormdistribution[2];   // 用来记录1000个事件的相关拟合程度vsStripid的图像
  TDirectory *dirtest[2]; // 用来记录一些测试图像，比如拟合的图像
  for (int i = 0; i < 2; i++) {
    dir[i] = fDeconv->mkdir(Form("%s", i == 0 ? "Normal" : "Good"));
    dirtest[i] = fDeconv->mkdir(Form("%s", i == 0 ? "X" : "Y"));
  }
  for(int j=0;j<2;j++){
    dirnormdistribution[j] = dir[0]->mkdir(Form("%s_normdistribution", j == 0 ? "X" : "Y"));
  }
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      dirwaveform[i][j] = dir[i]->mkdir(Form("%s_waveform", j == 0 ? "X" : "Y"));
      dirwaveform2d[i][j] = dir[i]->mkdir(Form("%s_waveform2d", j == 0 ? "X" : "Y"));
      dircharge2d[i][j] = dir[i]->mkdir(Form("%s_charge2d", j == 0 ? "X" : "Y"));
      dirapvdistribution[i][j] = dir[i]->mkdir(Form("%s_apvdistribution", j == 0 ? "X" : "Y"));
    }
  }

  // 在线统计量：不存储原始值，只累加 count/sum/sumsq，避免内存爆炸
  struct OnlineStat {
    long long n = 0;
    double sum = 0.0;
    double sumsq = 0.0;
    void add(double v) { n++; sum += v; sumsq += v * v; }
    double mean() const { return n > 0 ? sum / n : 0.0; }
    double se() const {
      if (n < 2) return 0.0;
      double m = mean();
      double var = (sumsq - n * m * m) / n;
      if (var < 0) var = 0;
      return std::sqrt(var / (n - 1));
    }
  };
  TH1D *hapv[2][2][apvsamples];                          // 每个apv点的采样分布图（堆分配），分有没有筛选，分X/Y向
  OnlineStat stripStats[2][2][256][apvsamples];          // 每条strip每个apv点的在线统计量，分有没有筛选，分X/Y向
  OnlineStat totalStats[2][2][apvsamples];               // 所有strip合计的在线统计量，分有没有筛选，分X/Y向
  TH1D *h_res[2];                                       // 看预期位置和用microTPC位置的残差分布，分X/Y向
  TH1D *h_res_charge[2];                                // 看预期的位置和用电荷中心法的残差分布，分X/Y向
  TH1D *h_test1[2];                                     //当clustersize为2时，看cluster的距离
  for (int i = 0; i < 2; i++) {
    h_test1[i] = new TH1D(Form("h_test1_%s", i == 0 ? "X" : "Y"), Form("Test histogram %s", i == 0 ? "X" : "Y"), 256, 0, 256);
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < apvsamples; k++) {
        hapv[i][j][k] = new TH1D(Form("%s_hapv_%s_apv%d", i == 0 ? "Normal" : "Good", j == 0 ? "X" : "Y", k), Form("APV%d distribution for %s signals on %s strips;ADC", k, i == 0 ? "Normal" : "Good", j == 0 ? "X" : "Y"), 1100, -100, 1000);
      }
    }
    h_res[i] = new TH1D(Form("hres_%s", i == 0 ? "X" : "Y"), Form("Residual distribution on %s strips;Residual [mm]", i == 0 ? "X" : "Y"), 100, -5, 5);
    h_res_charge[i] = new TH1D(Form("hres_charge_%s", i == 0 ? "X" : "Y"),
                               Form("Residual distribution using charge "
                                    "centroid on %s strips;Residual [mm]",
                                    i == 0 ? "X" : "Y"),
                               100, -5, 5);
  }


  // 000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
  // ===================================逐事例分析分割线===========================================

  for (int i = 0; i < nEntries; ++i) {
    if (i % m_progressInterval == 0) std::cout << "\r  Processing: " << i << "/" << nEntries << std::endl;
    
    trackTree->GetEntry(i);
    auto rawHits = parser->LoadEvent(eventID);
    if (rawHits.empty()) {
      std::cout << "Event " << eventID << ": No raw hits found, skipping." << std::endl;
      continue;
    }

    //===================================== begin: 处理每个DUT的事件数据 =====================================
    for (auto &det : deconvs) { // 这里目前默认是一个DUT,后续如果有多个DUT的话需要改成循环处理每个DUT的事件数据
      
      for(int isgood = 0; isgood < 2; isgood++) {
        StripHitDeconvs[isgood].clear();
        ClusterDeconvs[isgood].clear();
        for(int type = 0; type < 2; type++){
          clusternum[isgood][type] = 0;
        }
      }

      const int id = det->GetID();
      const auto &rawData = rawHits[det->GetID()];
      const auto *badStripMask = det->GetBadStripMask();

      vector<vector<short>> vec_ApvValue_temp[2];        // 保存每个条的每个apv点的adc值，用于后续画分布图
      vector<vector<double>> vec_tracksignal_temp[2];      // 保存每个条的每个apv点的拟合信号值，用于后续画分布图
      vector<StripHitDeconv> vec_StripHitDeconv_temp[2]; // 按照类型保存，后续聚类成clusterdeconv时用
      vector<int> vec_StripHitDeconvid_temp[2];          // 同上，但保存的是条号和类型信息，用于后续聚类成clusterdeconv

      if (rawData.size() == 0)  continue;

      //=== phase 4: 处理每个条的波形，提取峰值、过阈时间、半高宽、上升时间等信息 ===
      for (size_t k = 0; k < rawData.size(); ++k) {

        if (rawData[k].type < 0 || rawData[k].type >= 2) {
          std::cerr << "Warning: Invalid strip type " << rawData[k].type << " in event " << eventID << ", skipping." << std::endl;
          continue;
        }

        const auto &pedstdVec = pedStdByDet[id]->at(rawData[k].type);
        if (rawData[k].stripID < 1 || rawData[k].stripID > static_cast<int>(pedstdVec.size())) {
          std::cerr << "Warning: Strip ID " << rawData[k].stripID << " out of range [1, " << pedstdVec.size() << "] in event " << eventID << ", skipping." << std::endl;
          continue;
        }
        double noiseTh = 5.0 * pedstdVec[rawData[k].stripID - 1];

        StripHitDeconv sh;
        sh.isBad = false;
        if (badStripMask && badStripMask->IsBad(rawData[k].type, rawData[k].stripID)) {
          sh.isBad = true;
        }

        sh.ID = rawData[k].stripID;
        sh.type = rawData[k].type;

        const auto &waveform = rawData[k].adc;
        if (!AnalyzeWaveformFeatures(waveform, noiseTh, t0_per_strip[sh.type][sh.ID], &sh))
          continue;

        // 用 per-channel t0 做反卷积（基线），后续 per-type 再做精细调整
        double t0_samples = t0_per_strip[sh.type][sh.ID] / 25.0;
        if (t0_samples + fit_point > static_cast<int>(waveform.size()))
          continue;

        // 将波形转为 double 用于反卷积
        std::vector<double> smooth(waveform.size(), 0.0);
        for (size_t j = 0; j < waveform.size(); ++j) smooth[j] = waveform[j];

        std::vector<double> delayed_signal = AnalysisUtils::DelaySignal(smooth, t0_samples);
        std::vector<double> b_vec = std::vector<double>(delayed_signal.begin(), delayed_signal.begin() + fit_point);
        auto [x_vec, rnorm] = AnalysisUtils::SolveNNLSLawHanson(b_vec, deconvMatrix[sh.type]);
        const double charge = std::accumulate(x_vec.begin(), x_vec.end(), 0.0);
        if (charge <= 0.0) {
          continue;
        }
        if (!std::isfinite(rnorm)) {
          continue;
        }
        sh.chargetier = x_vec;
        sh.rnorm = rnorm / charge;
        sh.rnorm_no = rnorm;
        sh.t0 = t0_per_strip[sh.type][sh.ID]; // per-channel t0 基线 (ns)

        if (sh.rnorm > m_cuts.rnormlim[sh.type]) {
          continue;
        }

        sh.rawIndices = static_cast<int>(k);

        vec_StripHitDeconv_temp[sh.type].push_back(sh);
        vec_StripHitDeconvid_temp[sh.type].push_back(sh.ID);
        vec_ApvValue_temp[sh.type].push_back(waveform);
      }

      //======phase 5: 利用事例的平均信号，对每个条进行事例级t0修正 ======
      for(int type=0;type<2;type++){
        int n_sample_signal;
        if(vec_StripHitDeconv_temp[type].size() == 0)
          continue;
        else{
          n_sample_signal = vec_ApvValue_temp[type][0].size();
        }
        // 各条先用自己的 per-channel t0 延迟，再累加
        std::vector<double> vec_totalsignal;
        vec_totalsignal.resize(n_sample_signal, 0.0);
        for(int m=0;m<vec_ApvValue_temp[type].size();m++){
          auto &sh = vec_StripHitDeconv_temp[type][m];
          std::vector<double> delayed = AnalysisUtils::DelaySignal(vec_ApvValue_temp[type][m], t0_per_strip[type][sh.ID] / 25.0);
          for(int j=0;j<n_sample_signal;j++){
            vec_totalsignal[j] += delayed[j];
          }
        }

        // 搜索全局残差 shift（在 0 附近，因为各条已用自己 t0 对齐）
        double best_score = std::numeric_limits<double>::infinity();
        double best_shift_ns = 0.0;

        for (int m = -12; m <= 12; ++m) {
          std::vector<double> delayed_signal = AnalysisUtils::DelaySignal(vec_totalsignal, m/25.0);
          std::vector<double> b_vec = std::vector<double>(delayed_signal.begin(), delayed_signal.begin() + fit_point);
          auto [x_vec, rnorm] =
              AnalysisUtils::SolveNNLSLawHanson(b_vec, deconvMatrix[type]);
          const double score = rnorm;
          if (score < best_score) {
            best_score = score;
            best_shift_ns = m;
          }
        }
        // 用 per-channel t0 + 全局残差 重做反卷积
        for(int m=0;m<vec_StripHitDeconv_temp[type].size();m++){
          auto &sh = vec_StripHitDeconv_temp[type][m];
          double delay_samples = t0_per_strip[type][sh.ID] / 25.0 + best_shift_ns / 25.0;
          std::vector<double> delayed_signal = AnalysisUtils::DelaySignal(vec_ApvValue_temp[type][m], delay_samples);
          std::vector<double> b_vec = std::vector<double>(delayed_signal.begin(), delayed_signal.begin() + fit_point);
          auto [x_vec, rnorm] =
              AnalysisUtils::SolveNNLSLawHanson(b_vec, deconvMatrix[type]);
          sh.chargetier = x_vec;
          sh.rnorm = rnorm / std::accumulate(x_vec.begin(), x_vec.end(), 0.0);
          sh.t0 =  best_shift_ns;//t0_per_strip[type][sh.ID] + best_shift_ns; // per-channel + global
          vec_tracksignal_temp[type].push_back(b_vec);
        }
      }

      //======phase 6: 聚类，把相邻的条聚成cluster ======
      int hasCluster[2] = {0, 0};      
      vector<ClusterDeconv> vec_ClusterDeconv_temp; // cluster的临时变量
      for (int type = 0; type < 2; type++) {
        if (vec_StripHitDeconvid_temp[type].size() == 0)
          continue;
        for (int m = 0; m < vec_StripHitDeconvid_temp[type].size() - 1; m++) {
          for (int n = m + 1; n < vec_StripHitDeconvid_temp[type].size(); n++) {
            if (vec_StripHitDeconvid_temp[type][m] > vec_StripHitDeconvid_temp[type][n]) {
              swap(vec_StripHitDeconvid_temp[type][m], vec_StripHitDeconvid_temp[type][n]);
              swap(vec_StripHitDeconv_temp[type][m], vec_StripHitDeconv_temp[type][n]);
              swap(vec_tracksignal_temp[type][m], vec_tracksignal_temp[type][n]);
              swap(vec_ApvValue_temp[type][m], vec_ApvValue_temp[type][n]);
            }
          }
        }
      }
      // 聚类（调用独立函数）
      for (int type = 0; type < 2; type++) {
        ClusterStripHits(vec_StripHitDeconv_temp[type], badStripMask, m_cuts, eventID, type,
                         hasCluster[type], vec_ClusterDeconv_temp);
      }


      //======phase 7: 画图，和区分好坏 ======
      TH2D *h_waveform2d[2][2];
      TH2D *h_charge2d[2][2];
      TH1D *h_norm[2];
  
      for (int j = 0; j < 2; j++) {
        for (int k = 0; k < 2; k++) {
          h_waveform2d[j][k] = new TH2D(Form("h_waveform2d_evt%d_%s_%s", eventID, k == 0 ? "X" : "Y", j == 0 ? "Normal" : "Good"), Form("Event %d %s Strip Waveforms;Strip ID;Sample Index", eventID, k == 0 ? "X" : "Y"), 256, 0.5, 256.5, apvsamples, -0.5, apvsamples - 0.5);
          h_charge2d[j][k] = new TH2D(Form("h_charge2d_evt%d_%s_%s", eventID, k == 0 ? "X" : "Y", j == 0 ? "Normal" : "Good"), Form("Event %d %s Strip Charge Tiers;Strip ID;Position[mm]", eventID, k == 0 ? "X" : "Y"), 256, 0.5, 256.5, n_tier, 0, thickness);
          h_waveform2d[j][k]->SetMinimum(1);
        }
        h_norm[j] = new TH1D(Form("h_norm_evt%d_%s", eventID, j == 0 ? "X" : "Y"), Form("Event %d %s Fit Rnorm strip distribution;Rnorm", eventID, j == 0 ? "X" : "Y"), 256, 0.5, 256.5);
        h_norm[j]->SetOption("hist");
      }

      for (auto kclusterdeconv : vec_ClusterDeconv_temp) {
        if (kclusterdeconv.isValid && !kclusterdeconv.isBad) {
          for (int k = 0; k < 2; k++) {
            ClusterDeconvs[k].push_back(kclusterdeconv);      
          }
          clusternum[0][kclusterdeconv.type]++;
          clusternum[1][kclusterdeconv.type]++;
          for (int j = 0; j < kclusterdeconv.stripHitIndices.size(); ++j) {
            int stripindex = kclusterdeconv.stripHitIndices[j];
            auto kstripdeconv = vec_StripHitDeconv_temp[kclusterdeconv.type][stripindex];
            for (int k = 0; k < 2; ++k) {
              StripHitDeconvs[k].push_back(kstripdeconv);
              vector<short> waveform = vec_ApvValue_temp[kclusterdeconv.type][stripindex];
              for (int m = 0; m < waveform.size() && m < apvsamples; ++m) {
                hapv[k][kclusterdeconv.type][m]->Fill(waveform[m]);
                stripStats[k][kclusterdeconv.type][kstripdeconv.ID][m].add(waveform[m]);
                totalStats[k][kclusterdeconv.type][m].add(waveform[m]);
              }
              if (i < drawmaxeventid) {
                for (int m = 0; m < waveform.size() && m < apvsamples; ++m) {
                  h_waveform2d[k][kclusterdeconv.type]->SetBinContent(kstripdeconv.ID, m + 1, waveform[m]);
                }
                for (int m = 0; m < n_tier; ++m) {
                  h_charge2d[k][kclusterdeconv.type]->SetBinContent(kstripdeconv.ID, m + 1, kstripdeconv.chargetier[m]);
                }
                h_norm [kclusterdeconv.type]->Fill(kstripdeconv.ID, kstripdeconv.rnorm);
              }
            }
          }
        } else {
          ClusterDeconvs[0].push_back(kclusterdeconv);
          clusternum[0][kclusterdeconv.type]++;
          for (int j = 0; j < kclusterdeconv.stripHitIndices.size(); ++j) {
            int stripindex = kclusterdeconv.stripHitIndices[j];
            auto kstripdeconv = vec_StripHitDeconv_temp[kclusterdeconv.type][stripindex];
            StripHitDeconvs[0].push_back(kstripdeconv);
            vector<short> waveform = vec_ApvValue_temp[kclusterdeconv.type][stripindex];
            for (int m = 0; m < waveform.size() && m < apvsamples; ++m) {
              hapv[0][kclusterdeconv.type][m]->Fill(waveform[m]);
              stripStats[0][kclusterdeconv.type][kstripdeconv.ID][m].add(waveform[m]);
              totalStats[0][kclusterdeconv.type][m].add(waveform[m]);
            }
            if (i < drawmaxeventid) {
              for (int m = 0; m < waveform.size() && m < apvsamples; ++m) {
                h_waveform2d[0][kclusterdeconv.type]->SetBinContent(kstripdeconv.ID, m + 1, waveform[m]);
              }
              for (int m = 0; m < n_tier; ++m) {
                h_charge2d[0][kclusterdeconv.type]->SetBinContent(kstripdeconv.ID, m + 1, kstripdeconv.chargetier[m]);
              }
              h_norm[kclusterdeconv.type]->Fill(kstripdeconv.ID, kstripdeconv.rnorm);
            }
          }
        }

        if (i < drawmaxeventid) {
          for (int j = 0; j < kclusterdeconv.stripHitIndices.size(); ++j) {
            int stripindex = kclusterdeconv.stripHitIndices[j];
            auto kstripdeconv = vec_StripHitDeconv_temp[kclusterdeconv.type][stripindex];
            std::vector<double> waveform = vec_tracksignal_temp[kclusterdeconv.type][stripindex];
            std::vector<short> waveform_ori = vec_ApvValue_temp[kclusterdeconv.type][stripindex];
            std::vector<double>* deconvsignal = &sigzero_sample[kstripdeconv.type];
            TCanvas *c = new TCanvas(Form("c_evt%d_strip%d", eventID, kstripdeconv.ID), Form("Event %d Strip %d Waveform", eventID, kstripdeconv.ID), 800, 600);
            TGraph *gr_waveform = new TGraph();
            gr_waveform->SetTitle(Form("Event %d Strip %d Waveform %s t0 %.2f chi2 %.2f r2 %.2f scale %.2f rnorm %.3f(no %.3f);Time (ns)  ns;ADC", eventID, kstripdeconv.ID, kstripdeconv.type == 0 ? "X" : "Y", kstripdeconv.t0, kstripdeconv.chi2, kstripdeconv.r2, kstripdeconv.scale, kstripdeconv.rnorm, kstripdeconv.rnorm_no));
            gr_waveform->SetName(Form("gr%swaveform_evt%d_strip%d", kstripdeconv.type == 0 ? "x" : "y", eventID, kstripdeconv.ID));
            gr_waveform->SetLineColor(kBlue);
            gr_waveform->SetLineStyle(1);
            gr_waveform->SetMarkerStyle(20);
            TGraph *gr_waveform_ori = new TGraph();
            gr_waveform_ori->SetTitle(Form("Event %d Strip %d Original Waveform %s;Time (ns);ADC", eventID, kstripdeconv.ID, kstripdeconv.type == 0 ? "X" : "Y"));
            gr_waveform_ori->SetName(Form("gr%sorwaveform_evt%d_strip%d", kstripdeconv.type == 0 ? "x" : "y", eventID, kstripdeconv.ID));
            gr_waveform_ori->SetLineColor(kGreen+2);
            gr_waveform_ori->SetLineStyle(1);
            TGraph *gr_fit = new TGraph();
            gr_fit->SetTitle(Form("Event %d Strip %d Fit %s t0 %.2f ns;Time (ns);ADC", eventID, kstripdeconv.ID, kstripdeconv.type == 0 ? "X" : "Y", kstripdeconv.t0));
            gr_fit->SetName(Form("gr%sfitt_evt%d_strip%d", kstripdeconv.type == 0 ? "x" : "y", eventID, kstripdeconv.ID));
            gr_fit->SetLineColor(kRed);
            gr_fit->SetLineStyle(2);
            gr_fit->SetMarkerStyle(22);
            TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend->AddEntry(gr_waveform, "Raw Waveform", "l");
            legend->AddEntry(gr_fit, "Fit Waveform", "l");
            legend->AddEntry(gr_waveform_ori, "Original Waveform", "l");
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            // 仅用于画图：用加延迟的模板重建以进行视觉对齐。
            std::vector<double> fit_waveform(waveform.size(), 0.0);
            for (int tier = 0; tier < n_tier; ++tier) {
              for (int k = 0; k < deconvsignal->size(); ++k) {
                if (tier + k < fit_waveform.size())
                  fit_waveform[tier + k] += kstripdeconv.chargetier[tier] * deconvsignal->at(k);
              }
            }
            for (size_t k = 0; k < waveform.size(); ++k) {
              gr_waveform->SetPoint(k, k * 25.0, waveform[k]);
              gr_fit->SetPoint(k, k * 25.0, fit_waveform[k]);
            }
            for(size_t k = 0; k < waveform_ori.size(); ++k){
              gr_waveform_ori->SetPoint(k, k * 25.0, waveform_ori[k]);
            }
            gr_waveform->Draw("AL");
            gr_fit->Draw("L SAME");
            gr_waveform_ori->Draw("L SAME");
            legend->Draw();

            if(kstripdeconv.t0 == -20){
              dirtest[kstripdeconv.type]->WriteTObject(gr_waveform);
            }
            if (kclusterdeconv.isValid)
              dirwaveform[1][kstripdeconv.type]->WriteTObject(c);
            dirwaveform[0][kstripdeconv.type]->WriteTObject(c);
            delete gr_waveform;
            delete gr_fit;
            delete gr_waveform_ori;
            delete legend;
            delete c;
          }
        }
      }

          //画有两个cluster的事件，画出两个cluster的距离
        for(int j=0;j<2;j++){
          if(clusternum[1][j] == 2 && clusternum[1][1-j] != 0){
            double distance1 = 0.0, distance2 = 0.0;
            bool found = false;
            for(int k=0;k<ClusterDeconvs[1].size();k++){
              if(ClusterDeconvs[1][k].type == j){
                distance1 = *std::max_element(ClusterDeconvs[1][k].stripIDs.begin(), ClusterDeconvs[1][k].stripIDs.end());
                distance2 = *std::min_element(ClusterDeconvs[1][k].stripIDs.begin(), ClusterDeconvs[1][k].stripIDs.end());
                for(int m=k+1;m<ClusterDeconvs[1].size();m++){
                  if(ClusterDeconvs[1][m].type == j){
                    distance2 -= *std::max_element(ClusterDeconvs[1][m].stripIDs.begin(), ClusterDeconvs[1][m].stripIDs.end());
                    distance1 -= *std::min_element(ClusterDeconvs[1][m].stripIDs.begin(), ClusterDeconvs[1][m].stripIDs.end());
                    h_test1[j]->Fill(fabs(distance2)>fabs(distance1)?fabs(distance1):fabs(distance2));
                    found = true;
                    break;
                  }
                }
                break;
              }
            }
            if (!found) {
              std::cout << "Warning: Event " << eventID << " has 2 clusters in direction " << j << " but could not find both for distance calculation." << std::endl;
            }
         }
        }


      TVector3 predG = det->CalcHitFromTrack((*track));
      TVector3 predL = det->GlobalToLocal(predG);
      TVector3 trackDirGlobal((*track).kx, (*track).ky, 1);
      TVector3 localDir = det->GlobalDirToLocal(trackDirGlobal);//变成局域的斜率
      // std::cout << "Event " << eventID << ": " << trackDirGlobal.X() << " " << trackDirGlobal.Y() << " " << trackDirGlobal.Z() << std::endl;
      // std::cout << localDir.X() << " " << localDir.Y() << " " << localDir.Z() << std::endl;
      double hit[2];
      TGraph *gr_track[2];
      TCanvas *c_compare[2];
      for (int type = 0; type < 2; type++) {
        if (fabs(localDir.Z()) > 1e-9) {
          hit[type] = predL[type]+ localDir[type] / localDir.Z() * z_ref_position; // 考虑入射角的影响
        }
        gr_track[type] = new TGraph(n_tier);
        gr_track[type]->SetLineColor(kRed);
        gr_track[type]->SetLineWidth(2);
        for (int k = 0; k < n_tier; k++) {
          double z_tier = thickness / n_tier * (k + 0.5);
          double predStrip = predL[type] * 2.5;
          double predStrip_microTPC;
          if (fabs(localDir.Z()) > 1e-9) {
            predStrip = (predL[type] + localDir[type] / localDir.Z() * z_tier) * 2.5;
          }
          gr_track[type]->SetPoint(k, predStrip, k);
          // std::cout << "Event " << eventID << ": Tier " << k <<"z_tier:"<< z_tier << ": pred strip = " << predStrip << std::endl;
        }
        c_compare[type] = new TCanvas(Form("c_compare_evt%d_type%d", eventID, type), Form("Event %d %s Comparison", eventID, type == 0 ? "X" : "Y"), 800, 600);
        c_compare[type]->cd();
        if (hasCluster[0] && hasCluster[1]) { // 如果X和Y都有cluster才画对比图
          h_charge2d[1][type]->Draw("COLZ");
          gr_track[type]->Draw("L SAME");
          dircharge2d[1][type]->WriteTObject(c_compare[type]);
          h_charge2d[0][type]->Draw("COLZ");
          gr_track[type]->Draw("L SAME");
          dircharge2d[0][type]->WriteTObject(c_compare[type]);
          double minres_charge = std::numeric_limits<double>::max();
          double minres = std::numeric_limits<double>::max();
          for (auto kclusterdeconv : ClusterDeconvs[1]) { // 选择最靠近hit的cluster
            if (kclusterdeconv.type == type) {
              double res_charge = kclusterdeconv.chargeposition*pitchwidth[type] - hit[type];
              if (fabs(res_charge) < fabs(minres)) {
                minres = kclusterdeconv.microTPCposition*pitchwidth[type] - hit[type];
                minres_charge = res_charge;
              }
            }
          }
          bool isBadCluster = false;// 如果最靠近的cluster是坏条cluster，就不填分布了
          for(auto kclusterdeconv : ClusterDeconvs[0]){ // 选择最靠近hit的cluster
            if (kclusterdeconv.type == type) {
              double res_charge = kclusterdeconv.chargeposition*pitchwidth[type] - hit[type];
              if (fabs(res_charge) < fabs(minres) && kclusterdeconv.isBad){
                isBadCluster = true;
              }
            }
          }
          if (hasCluster[0] == 1 && hasCluster[1] == 1 && !isBadCluster) { // 只有当X和Y都有且最靠近的cluster不是坏条cluster时才填分布
            h_res[type]->Fill(minres);
            h_res_charge[type]->Fill(minres_charge);
          }
        }
        delete gr_track[type];
        delete c_compare[type];
      }

      if (i < drawmaxeventid) {
        for (int j = 0; j < 2; j++) {
          dirwaveform2d[0][j]->WriteTObject(h_waveform2d[0][j]);
          dircharge2d[0][j]->WriteTObject(h_charge2d[0][j]);
          dirnormdistribution[j]->WriteTObject(h_norm[j]);
        }
        if (hasCluster[0] && hasCluster[1]) { // 如果X和Y都有cluster才保存二维图
          for (int j = 0; j < 2; j++) {
            dirwaveform2d[1][j]->WriteTObject(h_waveform2d[1][j]);
          }
        }
      }
      // 填充 DeconvTree: 好事件和所有事件
      dutID = id;
      if (hasCluster[0] && hasCluster[1]) {
        tDeconvori[1]->Fill();
      }
      tDeconvori[0]->Fill();
      for (int j = 0; j < 2; j++) {
        delete h_charge2d[0][j];
        delete h_waveform2d[0][j];
        delete h_charge2d[1][j];
        delete h_waveform2d[1][j];
        delete h_norm[j];
      }
    }
  }
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < apvsamples; k++) {
        dirapvdistribution[i][j]->WriteTObject(hapv[i][j][k]);
      }
    }
  }

  std::vector<double> meansignal[2][2]; // 保存每个条的每个apv点的电荷均值，分有没有筛选，分X/Y向
  std::vector<double> sesignal[2][2];
  std::vector<double> meansignal_onestrip[2][2][256]; // 保存每个条的每个apv点的电荷均值，分有没有筛选，分X/Y向，分每条
  std::vector<double> sesignal_onestrip[2][2][256]; // 保存每个条的每个apv点的电荷标准误，分有没有筛选，分X/Y向，分每条  
  int hitnum[2][2][256]; // 保存每个条的每个apv点的hit数量，分有没有筛选，分X/Y向，分每条
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      meansignal[i][j].resize(apvsamples, 0.0);
      sesignal[i][j].resize(apvsamples, 0.0);
      for(int k=0; k<256; k++){
        meansignal_onestrip[i][j][k].resize(apvsamples, 0.0);
        sesignal_onestrip[i][j][k].resize(apvsamples, 0.0);
       }
    }
  }

  // 直接从在线统计量获取均值和标准误，无需遍历原始数据
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int m = 0; m < 256; m++) {
        hitnum[i][j][m] = static_cast<int>(stripStats[i][j][m][0].n); // 所有 apv 采样点的 n 相同，只取第一个
        for (int k = 0; k < apvsamples; k++) {
          meansignal_onestrip[i][j][m][k] = stripStats[i][j][m][k].mean();
          sesignal_onestrip[i][j][m][k] = stripStats[i][j][m][k].se();
        }
      }
      for (int k = 0; k < apvsamples; k++) {
        meansignal[i][j][k] = totalStats[i][j][k].mean();
        sesignal[i][j][k] = totalStats[i][j][k].se();
      }
    }
  }

     


  TGraph *gr_sigtrackmean_conv[2][2];
  TGraphErrors *gr_apvmean[2][2];
  TGraph *gr_fit_mean[2][2];
  TMultiGraph *mg_apvmean[2][2];
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      gr_apvmean[i][j] = new TGraphErrors();
      gr_apvmean[i][j]->SetTitle(Form("Mean signal of each APV channel for %s signals on %s "
                                      "strips;APV channel;Mean signal (ADC)",
                                      i == 0 ? "Normal" : "Good", j == 0 ? "X" : "Y"));
      gr_apvmean[i][j]->SetName(Form("gr_apvmean_%s_%s", i == 0 ? "normal" : "good", j == 0 ? "x" : "y"));
      gr_fit_mean[i][j] = new TGraph();
      gr_fit_mean[i][j]->SetTitle(Form("Fit of mean signal of each APV channel for %s signals on %s "
                                       "strips;APV channel;Mean signal (ADC)",
                                       i == 0 ? "Normal" : "Good", j == 0 ? "X" : "Y"));
      gr_fit_mean[i][j]->SetName(Form("gr_fit_mean_%s_%s", i == 0 ? "normal" : "good", j == 0 ? "x" : "y"));
      mg_apvmean[i][j] = new TMultiGraph();
      mg_apvmean[i][j]->SetTitle(Form("Mean signal of each APV channel for %s signals on %s "
                                      "strips;APV channel;Mean signal (ADC)",
                                      i == 0 ? "Normal" : "Good", j == 0 ? "X" : "Y"));
      gr_sigtrackmean_conv[i][j] = new TGraph();
      gr_sigtrackmean_conv[i][j]->SetTitle("Ref Signal;APV channel;Mean signal (ADC)");
      gr_sigtrackmean_conv[i][j]->SetName(Form("gr_sigtrackmean_conv_%s_%s", i == 0 ? "normal" : "good", j == 0 ? "x" : "y"));
    }
  }

  
  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
        double dt,search_min,search_max,step;
        dt = 25.0;
        search_min = 0;
        search_max = 125;
        step = 1;
        std::vector<double> b_vec = std::vector<double>(meansignal[i][j].begin(), meansignal[i][j].begin()+delayindex+fit_point);
        std::vector<double> sigtrackmean_sample = std::vector<double>(sigtrackmean_conv[j]->begin(), sigtrackmean_conv[j]->end());
        AnalysisUtils::DelayFitResult fitResult = AnalysisUtils::FitDelayLeastSquaresImpl_TGraph(b_vec, gr_trackmean_conv_tran[j],dt,search_min,search_max,step,false,delayindex,delayindex+fit_point);
        for(int k=0;k<sigtrackmean_conv[j]->size();k++){
          gr_sigtrackmean_conv[i][j]->SetPoint(k,k,sigtrackmean_sample[k]*fitResult.scale);
        }
        std::cout <<  (j == 0 ? "X" : "Y") << " " << (i == 0 ? "Normal" : "Good")<< ": Delay = " << fitResult.delay << ", Chi2 = " << fitResult.chi2 << ", Scale = " << fitResult.scale << std::endl;
    }
  }

  std::vector<double> chargetiermean[2][2]; // 保存每个条的每个电荷层的均值，分有没有筛选，分X/Y向
  std::vector<double> fit_meansignal[2][2]; // 保存每个条的每个apv点的拟合值，分有没有筛选，分X/Y向
  double rnormmean[2][2];
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      fit_meansignal[i][j].resize(apvsamples, 0.0);
      std::vector<double> signalmean_delay = AnalysisUtils::DelaySignal(meansignal[i][j], delayindex2d[j]);
      if (delayindex2d[j] + fit_point > static_cast<int>(meansignal[i][j].size())) {
        chargetiermean[i][j].clear();
        rnormmean[i][j] = 0.0;
        continue;
      }
      auto [x_vec, rnorm_vec] = AnalysisUtils::SolveNNLSLawHanson(std::vector<double>(signalmean_delay.begin() , signalmean_delay.begin() + fit_point), deconvMatrix[j]);
      chargetiermean[i][j] = x_vec;
      rnormmean[i][j] = rnorm_vec/std::accumulate(x_vec.begin(), x_vec.end(), 0.0);
      std::cout << "APV " << (j == 0 ? "X" : "Y") << " " << (i == 0 ? "Normal" : "Good") << ": Rnorm of fit = " << rnorm_vec << ", Rnorm mean = " << rnormmean[i][j] << std::endl;
      for (int tier = 0; tier < n_tier; tier++) {
        std::cout << "APV " << (j == 0 ? "X" : "Y") << " " << (i == 0 ? "Normal" : "Good") << " Strip Tier " << tier << ": Mean charge = " << chargetiermean[i][j][tier] << std::endl;
        for (int k = 0; k < sigzero_sample[j].size(); k++) {
          if (tier + k < fit_meansignal[i][j].size())
            fit_meansignal[i][j][tier + k] += x_vec[tier] * sigzero_sample[j][k];
        }
      }

      for (int k = 0; k < signalmean_delay.size(); k++) {
        gr_apvmean[i][j]->SetPoint(k, k, signalmean_delay[k]);
        gr_fit_mean[i][j]->SetPoint(k, k, fit_meansignal[i][j][k]);
      }
    }
  }


  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      gr_apvmean[i][j]->SetLineColor(kBlue);
      gr_fit_mean[i][j]->SetLineColor(kRed);
      gr_sigtrackmean_conv[i][j]->SetLineColor(kGreen);
      gr_apvmean[i][j]->SetMarkerStyle(20);
      gr_fit_mean[i][j]->SetMarkerStyle(22);
      gr_sigtrackmean_conv[i][j]->SetMarkerStyle(24);
      mg_apvmean[i][j]->Add(gr_apvmean[i][j], "LP");
      mg_apvmean[i][j]->Add(gr_fit_mean[i][j], "LP");
      mg_apvmean[i][j]->Add(gr_sigtrackmean_conv[i][j], "LP");
      TCanvas *c_apvmean = new TCanvas(Form("c_apvmean_%s_%s", i == 0 ? "normal" : "good", j == 0 ? "x" : "y"), Form("Mean signal of each APV channel for %s signals on %s strips", i == 0 ? "Normal" : "Good", j == 0 ? "X" : "Y"), 800, 600);
      TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
      legend->AddEntry(gr_apvmean[i][j], "Mean signal", "lp");
      legend->AddEntry(gr_fit_mean[i][j], "Fit", "lp");
      legend->AddEntry(gr_sigtrackmean_conv[i][j],"min chi2","lp");
      legend->SetBorderSize(0);
      legend->SetFillStyle(0);
      mg_apvmean[i][j]->Draw("ALP");
      legend->Draw("SAME");
      // fDeconv->WriteTObject(mg_apvmean[i][j]);
      fDeconv->WriteTObject(c_apvmean);
      delete c_apvmean;
      delete legend;
    }
  }



  //256个通道的每个apv点的均值和标准误
  TDirectory *dir_apvmean = fDeconv->mkdir("APVMean");
  TDirectory *dir_apvmean_fit[2][2] = {nullptr};
  TH2D *h_apvmean_2d[2][2];
  TH1D *h_delay[2][2];
  TH1D *h_delay_strip[2][2];
  double fitted_delay[2][256]; // 保存 Good 事件的 per-channel 拟合 delay
  for (int type = 0; type < 2; type++)
    for (int strip = 0; strip < 256; strip++)
      fitted_delay[type][strip] = -1.0;

  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
      dir_apvmean_fit[i][j] = fDeconv->mkdir(Form("APVMeanFit_%s_%s", j == 0 ? "X" : "Y", i == 0 ? "Normal" : "Good"));
      h_apvmean_2d[i][j] = new TH2D(Form("h_apvmean_2d_%s_%s", j == 0 ? "X" : "Y", i == 0 ? "Normal" : "Good"), Form("Mean signal of each APV channel for %s signals on %s strips;APV channel;APV sample", i == 0 ? "Normal" : "Good", j == 0 ? "X" : "Y"), 256, -0.5, 255.5, apvsamples, -0.5, apvsamples-0.5);
      h_delay[i][j] = new TH1D(Form("h_delay_%s_%s", j == 0 ? "X" : "Y", i == 0 ? "Normal" : "Good"), Form("Delay of mean signal of each APV channel for %s signals on %s strips;Delay (ns);Entries", i == 0 ? "Normal" : "Good", j == 0 ? "X" : "Y"), 100, 50, 150);
      h_delay_strip[i][j] = new TH1D(Form("h_delay_strip_%s_%s", j == 0 ? "X" : "Y", i == 0 ? "Normal" : "Good"), Form("Delay of mean signal of each APV channel for %s signals on %s strips;Delay (ns);Entries", i == 0 ? "Normal" : "Good", j == 0 ? "X" : "Y"), 256, -0.5, 255.5);

      TDirectory *dir_strip = dir_apvmean->mkdir(Form("%s_%s", j == 0 ? "X" : "Y", i == 0 ? "Normal" : "Good"));
      for(int k=0; k<256; k++){
        if(hitnum[i][j][k] == 0) continue; // 如果这个条没有hit，就不画了
        for(int m=0;m<meansignal_onestrip[i][j][k].size(); m++){
          h_apvmean_2d[i][j]->SetBinContent(k+1, m+1, meansignal_onestrip[i][j][k][m]);
        }
        double dt,search_min,search_max,step;
        dt = 25.0;
        search_min = 0;
        search_max = 125;
        step = 1;
        std::vector<double> b_vec = std::vector<double>(meansignal_onestrip[i][j][k].begin(), meansignal_onestrip[i][j][k].begin()+delayindex+fit_point);
        std::vector<double> sigtrackmean_sample2 = std::vector<double>(sigtrackmean_conv[j]->begin(), sigtrackmean_conv[j]->end());
        AnalysisUtils::DelayFitResult fitResult = AnalysisUtils::FitDelayLeastSquaresImpl_TGraph(b_vec, gr_trackmean_tran[j],dt,search_min,search_max,step,false,delayindex,delayindex+fit_point);
        if (i == 1) fitted_delay[j][k] = fitResult.delay; // 保存 Good 事件拟合结果
        std::vector<double> sigstrip_delay = AnalysisUtils::DelaySignal(meansignal_onestrip[i][j][k], fitResult.delay/25.);

        TMultiGraph *mg_strip_fit = new TMultiGraph();
        mg_strip_fit->SetTitle(Form("Mean signal of APV channel %d for %s signals on %s strips(hitnum=%d) delay: %f;APV channel;Mean signal (ADC)", k, i == 0 ? "Normal" : "Good", j == 0 ? "X" : "Y", hitnum[i][j][k], fitResult.delay));
        mg_strip_fit->SetName(Form("mg_strip_fit_strip%d", k));
        TGraph *gr_strip_fit = new TGraph();
        gr_strip_fit->SetTitle(Form("Fit of mean signal of APV channel %d for %s signals on %s strips(hitnum=%d) delay: %f;APV channel;Mean signal (ADC)", k, i == 0 ? "Normal" : "Good", j == 0 ? "X" : "Y", hitnum[i][j][k], fitResult.delay));
        gr_strip_fit->SetName(Form("gr_fit_mean_strip%d", k));
        gr_strip_fit->SetLineColor(kRed);
        gr_strip_fit->SetLineStyle(1);
        gr_strip_fit->SetMarkerStyle(22);
        for(int m=0; m<sigtrackmean_sample2.size(); m++){
          gr_strip_fit->SetPoint(m, m, sigtrackmean_sample2[m]*fitResult.scale);
        }
        mg_strip_fit->Add(gr_strip_fit, "LP");
        TGraph *gr_sigstrip_delay = new TGraph();
        gr_sigstrip_delay->SetTitle(Form("Delay of mean signal of APV channel %d for %s signals on %s strips(hitnum=%d) delay: %f;APV channel;Mean signal (ADC)", k, i == 0 ? "Normal" : "Good", j == 0 ? "X" : "Y", hitnum[i][j][k], fitResult.delay));
        gr_sigstrip_delay->SetName(Form("gr_sigstrip_delay_strip%d", k));
        gr_sigstrip_delay->SetLineColor(kGreen);
        gr_sigstrip_delay->SetLineStyle(2);
        gr_sigstrip_delay->SetMarkerStyle(24);
        for(int m=0; m<sigstrip_delay.size(); m++){
          gr_sigstrip_delay->SetPoint(m, m, sigstrip_delay[m]);
        }
        mg_strip_fit->Add(gr_sigstrip_delay, "LP");
        h_delay[i][j]->Fill(fitResult.delay);
        h_delay_strip[i][j]->Fill(k, fitResult.delay);
        dir_apvmean_fit[i][j]->WriteTObject(mg_strip_fit);


        TGraphErrors *gr_mean = new TGraphErrors();
        gr_mean->SetTitle(Form("Mean signal of APV channel %d for %s signals on %s strips(hitnum=%d);APV channel;Mean signal (ADC)", k, i == 0 ? "Normal" : "Good", j == 0 ? "X" : "Y", hitnum[i][j][k]));
        gr_mean->SetName(Form("gr_mean_strip%d", k));
        for(int m=0; m<apvsamples; m++){
          gr_mean->SetPoint(m, m, meansignal_onestrip[i][j][k][m]);
          gr_mean->SetPointError(m, 0, sesignal_onestrip[i][j][k][m]);
        }
        dir_strip->WriteTObject(gr_mean);
        delete gr_mean;
        delete mg_strip_fit;
      }
    }
  }

  // 将拟合得到的 per-channel delay 更新回 delay.root (仅在启用逐通道t0时)
  if (m_usePerChannelT0) {
    TFile *ft0update = TFile::Open(t0cacheFile.c_str(), "UPDATE");
    if (ft0update && !ft0update->IsZombie()) {
      std::cout << "Updating per-channel t0 in " << t0cacheFile << " from fitted delays..." << std::endl;
      const char *t0names[2] = {"t0_strip_arr_x", "t0_strip_arr_y"};
      for (int type = 0; type < 2; type++) {
        TH1D h_t0(t0names[type], t0names[type], 256, 0, 256);
        int nUpdated = 0;
        for (int strip = 0; strip < 256; strip++) {
          if (fitted_delay[type][strip] > 0) {
            h_t0.SetBinContent(strip + 1, fitted_delay[type][strip]);
            nUpdated++;
          } else {
            h_t0.SetBinContent(strip + 1, t0_per_strip[type][strip]);
          }
        }
        h_t0.Write("", TObject::kOverwrite);
        std::cout << "  " << (type == 0 ? "X" : "Y") << ": updated " << nUpdated << " / 256 channels" << std::endl;
      }
      ft0update->Close();
      delete ft0update;
      std::cout << "delay.root updated." << std::endl;
    } else {
      if (ft0update) { ft0update->Close(); delete ft0update; }
      std::cerr << "Warning: Could not open " << t0cacheFile << " for update." << std::endl;
    }
  }

  // ====== 坏道检测 ======
  if (m_detectBadStrips) {
    TDirectory *dirBadStrips = fDeconv->mkdir("BadStrips");
    for (int type = 0; type < 2; type++) {
      TH1D *h_badStrip = new TH1D(Form("h_badStrip_type%d", type),
          Form("Data-driven bad strip mask for %s;Strip ID;isBad",
               type == 0 ? "X" : "Y"),
          256, -0.5, 255.5);
      int nBad = 0;
      double meanHitRate[2] = {0.0, 0.0};
      int nStripsWithHits[2] = {0, 0};
      for (int strip = 0; strip < 256; strip++) {
        if (hitnum[1][type][strip] > 0) {
          meanHitRate[type] += hitnum[1][type][strip];
          nStripsWithHits[type]++;
        }
      }
      if (nStripsWithHits[type] > 0)
        meanHitRate[type] /= nStripsWithHits[type];
      std::cout << "[BadStrip] Type " << (type == 0 ? "X" : "Y")
                << " mean hit rate (good events): " << meanHitRate[type]
                << " (from " << nStripsWithHits[type] << " active strips)" << std::endl;

      for (int strip = 0; strip < 256; strip++) {
        bool isBad = false;
        // 判断标准1: 命中数为0（完全没有信号的死道）
        if (hitnum[1][type][strip] == 0) {
          isBad = true;
        }
        // 判断标准2: 命中数远低于平均（< 平均的10%且不是0）
        else if (meanHitRate[type] > 0 &&
                 hitnum[1][type][strip] < meanHitRate[type] * 0.1) {
          isBad = true;
        }
        if (isBad) {
          h_badStrip->Fill(strip, 1.0);
          nBad++;
        }
      }
      std::cout << "[BadStrip] Type " << (type == 0 ? "X" : "Y")
                << ": detected " << nBad << " bad strips out of 256" << std::endl;
      dirBadStrips->WriteTObject(h_badStrip);
      delete h_badStrip;
    }
  }

  // 保存 DeConvInfo.root 并关闭文件
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      dircharge2d[i][j]->Write();
      dirwaveform2d[i][j]->Write();
      dirwaveform[i][j]->Write();
      delete mg_apvmean[i][j];
    }
  }
  fDeconv->cd();  // 目录 Write 后 gDirectory 可能变了，回到文件根
  for (int i = 0; i < 2; i++) {
    tDeconvori[i]->Write();
    dirtest[i]->Write();
  }
  tDeconv->Write();
  fDeconv->Write();
  // 保存诊断直方图并清理
  for (int i = 0; i < 2; i++) {
    fDeconv->WriteTObject(h_res[i]);
    fDeconv->WriteTObject(h_res_charge[i]);
    fDeconv->WriteTObject(h_pedestal[i]);
    fDeconv->WriteTObject(h_test1[i]);
    delete h_res[i];
    delete h_res_charge[i];
    delete h_pedestal[i];
    delete h_test1[i];
    for (int j = 0; j < 2; j++) {
      fDeconv->WriteTObject(h_delay[j][i]);
      fDeconv->WriteTObject(h_apvmean_2d[i][j]);
      fDeconv->WriteTObject(h_delay_strip[j][i]);
      delete h_delay[j][i];
      delete h_apvmean_2d[j][i];
      delete h_delay_strip[j][i];
    }
  }
  fDeconv->Close();

  if (fmean_sample) {
    fmean_sample->Close();
    delete fmean_sample;
  }
  if (f) {
    f->Close();
    delete f;
  }

  std::cout << "\n=== Deconvolution Phase Complete ===" << std::endl;
  std::cout << "Output: " << deconvFile << std::endl;
  return true;
}

// ================================================================
// Phase 2: 对齐和位置分辨分析
// 从 DeConvInfo.root 读取反卷积结果，重建 Event 对象，进行探测器对齐
// ================================================================
bool DeconvAnalysisScript::RunAlignment() {
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);

  std::cout << "\n=== Alignment Phase ===" << std::endl;

  // ---- 将 m_cuts 中的可配置参数同步到全局变量 ----
  pitchwidth[0] = m_cuts.pitchwidth[0];  pitchwidth[1] = m_cuts.pitchwidth[1];
  thickness     = m_cuts.thickness;
  z_ref_position = m_cuts.z_ref_position;
  // -------------------------------------------------

  auto &factory = DetectorFactory::GetInstance();
  const auto &deconvs = factory.GetDetectorsByRole(Detector::Role::DeConv);
  if (deconvs.empty()) {
    std::cerr << "Error: No DeConv detectors found!" << std::endl;
    return false;
  }

  // 1. 读取 TrackInfo.root
  string trackFile = GetOutputDir() + "TrackInfo.root";
  TFile *fTrack = TFile::Open(trackFile.c_str(), "READ");
  if (!fTrack || fTrack->IsZombie()) {
    std::cerr << "Error: Cannot open " << trackFile << std::endl;
    return false;
  }
  TTree *trackTree = (TTree *)fTrack->Get("Tracks");
  if (!trackTree) {
    std::cerr << "Error: No Tracks tree in " << trackFile << std::endl;
    fTrack->Close(); delete fTrack;
    return false;
  }

  // 2. 读取 DeConvInfo.root
  string deconvFile = GetOutputDir() + "DeConvInfo.root";
  TFile *fDeconvIn = TFile::Open(deconvFile.c_str(), "READ");
  if (!fDeconvIn || fDeconvIn->IsZombie()) {
    std::cerr << "Error: Cannot open " << deconvFile << std::endl;
    std::cerr << "Please run deconvolution first (mode=deconv or mode=full)!" << std::endl;
    fTrack->Close(); delete fTrack;
    return false;
  }
  TTree *tDeconvIn = (TTree *)fDeconvIn->Get("DeconvTree_Good"); // Good events
  if (!tDeconvIn) {
    std::cerr << "Error: No DeconvTree_Good in " << deconvFile << std::endl;
    fDeconvIn->Close(); delete fDeconvIn;
    fTrack->Close(); delete fTrack;
    return false;
  }

  // 3. 建立 track 查找表 (eventID → Track)
  std::map<int, std::pair<Track, double>> trackMap; // eventID → (Track, t0)
  Int_t trackEventID;
  Track *trackPtr = nullptr;
  double sigTime;
  trackTree->SetBranchAddress("eventID", &trackEventID);
  trackTree->SetBranchAddress("track", &trackPtr);
  trackTree->SetBranchAddress("t0", &sigTime);
  for (Long64_t i = 0; i < trackTree->GetEntries(); ++i) {
    trackTree->GetEntry(i);
    trackMap[trackEventID] = {*trackPtr, sigTime};
  }
  fTrack->Close(); delete fTrack;

  // 4. 从 DeconvTree 重建 Event 对象
  Int_t evtID;
  Int_t dutID_in;
  vector<StripHitDeconv> *stripHitsIn = nullptr;
  vector<ClusterDeconv> *clustersIn = nullptr;
  tDeconvIn->SetBranchAddress("eventID", &evtID);
  tDeconvIn->SetBranchAddress("dutID", &dutID_in);
  tDeconvIn->SetBranchAddress("StripHitDeconvs", &stripHitsIn);
  tDeconvIn->SetBranchAddress("ClusterDeconvs", &clustersIn);

  std::vector<Event> singleClusterEvents;
  int nReconstructed = 0;
  Long64_t nDeconvEntries = tDeconvIn->GetEntries();
  std::cout << "Reading " << nDeconvEntries << " deconv events..." << std::endl;

  for (Long64_t i = 0; i < nDeconvEntries; ++i) {
    tDeconvIn->GetEntry(i);
    auto trackIt = trackMap.find(evtID);
    if (trackIt == trackMap.end()) continue;

    Event evt{.eventID = int(evtID), .track = trackIt->second.first, .t0 = trackIt->second.second};

    for (auto &det : deconvs) {
      auto detEvt = std::make_shared<DetectorFrame>(*det);
      auto &mutableStrips = detEvt->GetMutableStripHitsDeconv();
      auto &mutableClusters = detEvt->GetMutableClustersDeconv();
      mutableStrips = *stripHitsIn;
      mutableClusters = *clustersIn;

      // 分离坏条cluster
      auto &mutableBad = detEvt->GetMutableClustersDeconvBad();
      for (auto &c : *clustersIn) {
        if (c.isBad) mutableBad.push_back(c);
      }

      evt.detectorFramesMap[det->GetID()] = detEvt;
    }

    nReconstructed++;

    // 筛选单cluster好事件
    int goodX = 0, goodY = 0;
    bool hasLargeCluster = true;
    for (auto &c : *clustersIn) {
      if (c.isValid && !c.isBad) {
        if (c.type == 0) goodX++;
        if (c.type == 1) goodY++;
        if (c.size <= 1) hasLargeCluster = false;
      }
    }
    if (goodX == 1 && goodY == 1 && hasLargeCluster) {
      singleClusterEvents.push_back(evt);
    }
  }
  fDeconvIn->Close(); delete fDeconvIn;

  std::cout << "Reconstructed " << nReconstructed << " events, "
            << singleClusterEvents.size() << " single-cluster events for alignment." << std::endl;

  if (singleClusterEvents.size() < 10) {
    std::cerr << "Error: Too few single-cluster events for alignment!" << std::endl;
    return false;
  }

  // 5. 创建对齐输出文件
  string alignFile = GetOutputDir() + "AlignResult.root";
  TFile *fAlign = new TFile(alignFile.c_str(), "RECREATE");
  TDirectory *dir_afteralign = fAlign->mkdir("AfterAlignment");

  // 6. 对齐: 粗对齐 3D
  std::cout << "\n--- Coarse Alignment (3D: dx, dy, rotZ) ---" << std::endl;
  for (auto &det : deconvs) {
    RunDeConvAlign(singleClusterEvents, det, det->GetID(), true, true, 5);
  }

  // 7. 第一次位置分辨（粗对齐后）
  std::map<int, std::pair<double, double>> res_map_afteralign;
  res_map_afteralign = Positionanalysis(singleClusterEvents, "microTPC", dir_afteralign);
  std::cout << "(X)mean=" << res_map_afteralign[0].first << ", sigma=" << res_map_afteralign[0].second << std::endl;
  std::cout << "(Y)mean=" << res_map_afteralign[1].first << ", sigma=" << res_map_afteralign[1].second << std::endl;

  // 8. 滤波: 剔除残差过大的事件
  std::vector<Event> filtered;
  filtered.reserve(singleClusterEvents.size());
  for (auto &event : singleClusterEvents) {
    bool isGood = true;
    for (const auto &framePair : event.detectorFramesMap) {
      int detID = framePair.first;
      auto &detEvt = framePair.second;
      auto &clusterDeconvs = detEvt->ClustersDeconv();
      auto alignedDet = factory.GetDetector(detID);
      if (!alignedDet) continue;
      const Detector &det = *alignedDet;
      TVector3 predG = det.CalcHitFromTrack(event.track);
      TVector3 predL = det.GlobalToLocal(predG);
      TVector3 trackDirGlobal(event.track.kx, event.track.ky, 1);
      TVector3 localDir = det.GlobalDirToLocal(trackDirGlobal);
      for (auto &cluster : clusterDeconvs) {
        double hit;
        if (fabs(localDir.Z()) > 1e-9)
          hit = predL[cluster.type] + localDir[cluster.type] / localDir.Z() * z_ref_position;
        else
          hit = predL[cluster.type];
        double res = cluster.microTPCposition * pitchwidth[cluster.type] - hit;
        if (res > res_map_afteralign[cluster.type].first + 5 * res_map_afteralign[cluster.type].second ||
            res < res_map_afteralign[cluster.type].first - 5 * res_map_afteralign[cluster.type].second) {
          isGood = false;
          break;
        }
      }
    }
    if (isGood) filtered.push_back(event);
  }
  std::cout << "After filtering: " << filtered.size() << " events." << std::endl;

  // 9. 精对齐 6D + 逐层对齐
  std::cout << "\n--- Fine Alignment (6D + tier-by-tier) ---" << std::endl;
  for (auto &det : deconvs) {
    RunDeConvAlign(filtered, det, det->GetID(), true, true);
    RunDeConvAlign(filtered, det, det->GetID(), true, false);
    RunDeConvAlign(filtered, det, det->GetID(), false, false);
  }

  // 10. 位置分辨分析
  TH1D *h_res_afteralign[2];
  TH1D *h_res_charge_afteralign[2];
  TH1D *h_slope_res_afteralign[2];
  std::vector<TH1D*> h_res_tier[2];
  for (int type = 0; type < 2; type++) {
    h_res_afteralign[type] = new TH1D(Form("h_res_afteralign_type%d", type),
        Form("Residual after alignment for %s;Residual [mm]", type == 0 ? "X" : "Y"), 200, -2, 2);
    h_res_charge_afteralign[type] = new TH1D(Form("h_res_charge_afteralign_type%d", type),
        Form("Charge residual after alignment for %s;Residual [mm]", type == 0 ? "X" : "Y"), 200, -2, 2);
    h_slope_res_afteralign[type] = new TH1D(Form("h_slope_res_afteralign_type%d", type),
        Form("Slope residual for %s;Residual", type == 0 ? "X" : "Y"), 400, -0.2, 0.2);
    h_res_tier[type].resize(n_tier);
    for (int tier = 0; tier < n_tier; tier++) {
      h_res_tier[type][tier] = new TH1D(Form("h_res_tier_t%d_type%d", tier, type),
          Form("%s Tier %d residual;Residual [mm]", type == 0 ? "X" : "Y", tier), 200, -2, 2);
    }
  }

  std::vector<double> sum_tier_res[2];
  std::vector<double> sum_tier_cross[2];
  for (int type = 0; type < 2; type++) {
    sum_tier_res[type].resize(n_tier, 0.0);
    sum_tier_cross[type].resize(n_tier * n_tier, 0.0);
  }
  int n_tier_complete[2] = {0, 0};
  std::vector<double> res_tier_arr[2];
  res_tier_arr[0].resize(n_tier);
  res_tier_arr[1].resize(n_tier);
  int entry_count = 0;
  TTree *tres_afteralign = new TTree("trees_afteralign", "Residuals after alignment");
  tres_afteralign->Branch("entry_count", &entry_count, "entry_count/I");
  tres_afteralign->Branch("res_tier_x", &res_tier_arr[0]);
  tres_afteralign->Branch("res_tier_y", &res_tier_arr[1]);

  for (auto &event : filtered) {
    for (const auto &framePair : event.detectorFramesMap) {
      for (int type = 0; type < 2; type++)
        for (int tier = 0; tier < n_tier; tier++)
          res_tier_arr[type][tier] = std::numeric_limits<double>::infinity();

      int detID = framePair.first;
      auto &detEvt = framePair.second;
      auto &clusterDeconvs = detEvt->ClustersDeconv();
      auto alignedDet = factory.GetDetector(detID);
      if (!alignedDet) continue;
      const Detector &det = *alignedDet;
      TVector3 predG = det.CalcHitFromTrack(event.track);
      TVector3 predL = det.GlobalToLocal(predG);
      TVector3 trackDirGlobal(event.track.kx, event.track.ky, 1);
      TVector3 localDir = det.GlobalDirToLocal(trackDirGlobal);
      bool all_valid2 = true;

      for (auto &cluster : clusterDeconvs) {
        double hit;
        if (fabs(localDir.Z()) > 1e-9)
          hit = predL[cluster.type] + localDir[cluster.type] / localDir.Z() * z_ref_position;
        else
          hit = predL[cluster.type];

        double res = cluster.microTPCposition * pitchwidth[cluster.type] - hit;
        double res_charge = cluster.chargeposition * pitchwidth[cluster.type] - hit;
        h_res_afteralign[cluster.type]->Fill(res);
        h_res_charge_afteralign[cluster.type]->Fill(res_charge);

        if (fabs(localDir.Z()) > 1e-9 && cluster.k > -900) {
          double expectedSlope = localDir[cluster.type] / localDir.Z();
          h_slope_res_afteralign[cluster.type]->Fill(cluster.k * pitchwidth[cluster.type] - expectedSlope);
        }

        bool all_valid = true;
        for (int tier = 0; tier < n_tier; tier++) {
          if (!std::isfinite(cluster.pos[tier])) { all_valid = false; all_valid2 = false; continue; }
          double predTier = (fabs(localDir.Z()) > 1e-9)
              ? predL[cluster.type] + localDir[cluster.type] / localDir.Z() * (thickness / n_tier * (tier + 0.5))
              : predL[cluster.type];
          res_tier_arr[cluster.type][tier] = cluster.pos[tier] * pitchwidth[cluster.type] - predTier;
          h_res_tier[cluster.type][tier]->Fill(res_tier_arr[cluster.type][tier]);
        }
        if (all_valid) {
          for (int i = 0; i < n_tier; i++) {
            sum_tier_res[cluster.type][i] += res_tier_arr[cluster.type][i];
            for (int j = 0; j < n_tier; j++)
              sum_tier_cross[cluster.type][i * n_tier + j] += res_tier_arr[cluster.type][i] * res_tier_arr[cluster.type][j];
          }
          n_tier_complete[cluster.type]++;
        }
      }
      if (all_valid2) { entry_count++; tres_afteralign->Fill(); }
    }
  }

  // 打印分辨率
  res_map_afteralign = Positionanalysis(singleClusterEvents, "microTPC", dir_afteralign);
  std::cout << "\n(X) Final mean=" << res_map_afteralign[0].first
            << ", sigma=" << res_map_afteralign[0].second << std::endl;
  std::cout << "(Y) Final mean=" << res_map_afteralign[1].first
            << ", sigma=" << res_map_afteralign[1].second << std::endl;

  // 逐层相关矩阵
  TDirectory *dir_corr = dir_afteralign->mkdir("TierCorrelation");
  for (int type = 0; type < 2; type++) {
    if (n_tier_complete[type] < 2 * n_tier) continue;
    TH2D *h_corr = new TH2D(Form("h_tier_corr_type%d", type),
        Form("%s %d-tier residual correlation;Tier;Tier", type == 0 ? "X" : "Y", n_tier),
        n_tier, 0, n_tier, n_tier, 0, n_tier);
    double N = n_tier_complete[type];
    for (int i = 0; i < n_tier; i++) {
      for (int j = 0; j < n_tier; j++) {
        double cov = (sum_tier_cross[type][i * n_tier + j] - sum_tier_res[type][i] * sum_tier_res[type][j] / N) / (N - 1);
        double corr = cov / sqrt(
            (sum_tier_cross[type][i * n_tier + i] - sum_tier_res[type][i] * sum_tier_res[type][i] / N) / (N - 1) *
            (sum_tier_cross[type][j * n_tier + j] - sum_tier_res[type][j] * sum_tier_res[type][j] / N) / (N - 1));
        h_corr->SetBinContent(i + 1, j + 1, corr);
      }
    }
    dir_corr->WriteTObject(h_corr);
    delete h_corr;
  }

  // 保存结果
  fAlign->cd();
  tres_afteralign->Write();
  for (int i = 0; i < 2; i++) {
    fAlign->WriteTObject(h_res_afteralign[i]);
    fAlign->WriteTObject(h_res_charge_afteralign[i]);
    fAlign->WriteTObject(h_slope_res_afteralign[i]);
    for (int tier = 0; tier < n_tier; tier++) {
      fAlign->WriteTObject(h_res_tier[i][tier]);
      delete h_res_tier[i][tier];
    }
    delete h_res_afteralign[i];
    delete h_res_charge_afteralign[i];
    delete h_slope_res_afteralign[i];
  }
  dir_afteralign->Write();
  fAlign->Write();
  fAlign->Close();

  std::cout << "\n=== Alignment Phase Complete ===" << std::endl;
  std::cout << "Output: " << alignFile << std::endl;
  return true;
}


// ============ 对齐辅助函数 ============
void RunDeConvAlign(const std::vector<Event> &events, std::shared_ptr<Detector> detector, int detID,
                    bool is3D, bool useposition,double reslimit) { // 和DUTAlign类似的对齐函数，输入是已经处理好的事件数据，和需要对齐的detector对象，detector
                                       // ID，以及是否第一次对齐的标志, useposition 选择用逐层 pos[k] (false) 还是 microTPCposition (true)
  if (events.empty()) {
    std::cout << "No events to align for detector " << detID << std::endl;
    return;
  }

  std::cout << "[DeConv " << detID << "] Aligning (" << (useposition ? "microTPC" : "tier-by-tier")
            << ") ..." << std::endl;

  auto minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  minimizer->SetTolerance(0.005);
  minimizer->SetPrintLevel(1);
  UInt_t nPar = is3D ? 3 : 6; // 第一次对齐只调整x，y和rotz

  // 使用lambda捕获this和参数
  auto chi2Func = [&events, &detector, detID, is3D, useposition,reslimit](const double *par) -> double {
    {
      if (useposition) {
        return is3D ? DeConvPositionChi2_3D(par, events, detector, detID, reslimit) : DeConvPositionChi2_6D(par, events, detector, detID, reslimit);
      } else {
        return is3D ? DeConvTrackChi2_3D(par, events, detector, detID, reslimit) : DeConvTrackChi2_6D(par, events, detector, detID, reslimit);
      }
    }
  };

  ROOT::Math::Functor f(chi2Func, nPar);
  minimizer->SetFunction(f);

  // 使用当前探测器对齐作为初始猜测，避免第二次对齐从(0,0,0)开始覆盖粗对齐结果
  TVector3 curAlignPos = detector->GetAlignPos();
  TVector3 curAlignRot = detector->GetAlignRot();

  // === 诊断：扫描 chi2 vs dx, dy, rotZ，评估约束强度 ===
  {
    std::cout << "\n--- [DeConv " << detID << "] Chi2 scan " << (is3D ? "(coarse)" : "(fine)") << " ---" << std::endl;
    auto eval = [&](double dx_v, double dy_v, double rz_v) -> double {
      double par[6] = {dx_v, dy_v, 0, 0, 0, rz_v};
      if (!is3D) {
        par[2] = curAlignPos.Z(); par[3] = curAlignRot.X();
        par[4] = curAlignRot.Y(); par[5] = rz_v;
      }
      return is3D ? DeConvPositionChi2_3D(par, events, detector, detID,reslimit)
                         : DeConvPositionChi2_6D(par, events, detector, detID,reslimit);
    };
    auto scan = [&](const char *name, double center, double range, int n) {
      std::cout << "  " << name << ":" << std::endl;
      double bestVal = center, bestChi2 = 1e9, c_minus = 0, c_center = 0, c_plus = 0;
      double step = range * 2.0 / n;
      for (int i = 0; i <= n; i++) {
        double val = center - range + i * step;
        double chi2;
        if (name[1] == 'x') chi2 = eval(val, curAlignPos.Y(), curAlignRot.Z());
        else if (name[1] == 'y') chi2 = eval(curAlignPos.X(), val, curAlignRot.Z());
        else chi2 = eval(curAlignPos.X(), curAlignPos.Y(), val);
        std::cout << "    " << name << "=" << std::fixed << std::setprecision(4) << val
                  << "  c2=" << std::setprecision(5) << chi2;
        if (chi2 < bestChi2) { bestChi2 = chi2; bestVal = val; std::cout << " *"; }
        std::cout << std::endl;
        if (fabs(val - (center - step*2)) < step*0.6) c_minus = chi2;
        if (fabs(val - center) < step*0.6) c_center = chi2;
        if (fabs(val - (center + step*2)) < step*0.6) c_plus = chi2;
      }
      double h = step * 2;
      double curv = (h > 0 && c_center > 0) ? (c_plus + c_minus - 2*c_center) / (h*h) : 0;
      double sig = (curv > 1e-9) ? sqrt(2.0 / curv) : -1;
      std::cout << "  -> best " << name << "=" << bestVal << ", curv≈" << curv << ", σ_est≈" << sig << std::endl;
    };
    double R = is3D ? 5.0 : 10.0;
    scan("dx", curAlignPos.X(), R, 20);
    scan("dy", curAlignPos.Y(), R, 20);
    scan("rotZ", curAlignRot.Z(), 0.05, 20);
    std::cout << "--- End scan ---\n" << std::endl;
  }

  if (is3D) {
    // 参数: dx, dy, rotZ (alignment offsets)
    minimizer->SetVariable(0, "dx", curAlignPos.X(), 0.001);
    minimizer->SetVariable(1, "dy", curAlignPos.Y(), 0.001);
    minimizer->SetVariable(2, "rotZ", curAlignRot.Z(), 0.001);
  } else {
    // 参数: dx, dy, dz, rotX, rotY, rotZ (alignment offsets)
    minimizer->SetVariable(0, "dx", curAlignPos.X(), 0.001);
    minimizer->SetVariable(1, "dy", curAlignPos.Y(), 0.001);
    minimizer->SetVariable(2, "dz", curAlignPos.Z(), 0.001);
    minimizer->SetVariable(3, "rotX", curAlignRot.X(), 0.0001);
    minimizer->SetVariable(4, "rotY", curAlignRot.Y(), 0.0001);
    minimizer->SetVariable(5, "rotZ", curAlignRot.Z(), 0.0001);
    minimizer->SetVariableLimits(0,curAlignPos.X()-2,curAlignPos.X()+2);
     //限制dx在合理范围，避免过度调整
    minimizer->SetVariableLimits(1,curAlignPos.Y()-2,curAlignPos.Y()+2); 
    //限制dy在合理范围，避免过度调整
    minimizer->SetVariableLimits(2, curAlignPos.Z() - 20,
                                 curAlignPos.Z() + 20); // 限制dz在合理范围，避免过度调整
    // minimizer->SetFixedVariable(2, "dz", curAlignPos.Z());
    minimizer->SetVariableLimits(3, curAlignRot.X() - 0.1,
                                 curAlignRot.X() + 0.1); // 限制rotX在合理范围，避免过度调整
    minimizer->SetVariableLimits(4, curAlignRot.Y() - 0.1,
                                 curAlignRot.Y() + 0.1); // 限制rotY在合理范围，避免过度调整
    minimizer->SetVariableLimits(5, curAlignRot.Z() - 0.1,
                                 curAlignRot.Z() + 0.1); // 限制rotZ在合理范围，避免过度调整
  }
  minimizer->Minimize();

  // 应用结果
  if (is3D) {
    const double *result = minimizer->X();
    double dx = result[0];
    double dy = result[1];
    double rotZ = result[2];
    if (!std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(rotZ)) {
      std::cout << "[DeConv " << detID << "] Alignment result is NaN/Inf, skip apply." << std::endl;
      delete minimizer;
      return;
    }

    detector->SetAlignment(dx, dy, 0, 0, 0, rotZ);

    TVector3 pos = detector->GetPos();
    TVector3 rot = detector->GetRot();

    std::cout << "DUT " << detID << " coarse alignment: " << std::fixed << std::setprecision(5) << "\"position\": [" << pos.X() << "," << pos.Y() << "," << pos.Z() << "],"
              << "\"rotation\": [" << rot.X() << "," << rot.Y() << "," << rot.Z() << "]" << std::endl;
    std::cout << "DUT " << detID << " coarse alignment corrections: "
              << "dx=" << dx << ", dy=" << dy << ", rotZ=" << rotZ << std::endl;

    std::cout << "[DUT " << detID << "] Coarse alignment chi2: " << minimizer->MinValue() << std::endl;
    delete minimizer;
  } else {
    const double *result = minimizer->X();
    double dx = result[0];
    double dy = result[1];
    double dz = result[2];
    double rotX = result[3];
    double rotY = result[4];
    double rotZ = result[5];
    if (!std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dz) || !std::isfinite(rotX) || !std::isfinite(rotY) || !std::isfinite(rotZ)) {
      std::cout << "[DeConv " << detID << "] Alignment result is NaN/Inf, skip apply." << std::endl;
      delete minimizer;
      return;
    }

    detector->SetAlignment(dx, dy, dz, rotX, rotY, rotZ);

    TVector3 pos = detector->GetPos();
    TVector3 rot = detector->GetRot();

    std::cout << "DUT " << detID << " alignment: " << std::fixed << std::setprecision(5) << "\"position\": [" << pos.X() << "," << pos.Y() << "," << pos.Z() << "],"
              << "\"rotation\": [" << rot.X() << "," << rot.Y() << "," << rot.Z() << "]" << std::endl;
    std::cout << "DUT " << detID << " alignment corrections: "
              << "dx=" << dx << ", dy=" << dy << ", dz=" << dz << ", "
              << "rotX=" << rotX << ", rotY=" << rotY << ", rotZ=" << rotZ << std::endl;

    std::cout << "[DUT " << detID << "] Final chi2: " << minimizer->MinValue() << std::endl;

    delete minimizer;
  }
}

// ====== 协方差对齐 ======
void RunDeConvAlignCov(const std::vector<Event> &events, std::shared_ptr<Detector> detector, int detID,
                       bool is3D) {
  if (events.empty()) return;
  std::cout << "[DeConvCov " << detID << "] Aligning (covariance) ..." << std::endl;

  auto minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  minimizer->SetTolerance(0.005);
  minimizer->SetPrintLevel(1);
  UInt_t nPar = is3D ? 3 : 6;

  auto chi2Func = [&events, &detector, detID, is3D](const double *par) -> double {
    return is3D ? DeConvChi2Cov_3D(par, events, detector, detID)
                       : DeConvChi2Cov_6D(par, events, detector, detID);
  };
  ROOT::Math::Functor f(chi2Func, nPar);
  minimizer->SetFunction(f);

  TVector3 curAlignPos = detector->GetAlignPos();
  TVector3 curAlignRot = detector->GetAlignRot();

  if (is3D) {
    minimizer->SetVariable(0, "dx", curAlignPos.X(), 0.001);
    minimizer->SetVariable(1, "dy", curAlignPos.Y(), 0.001);
    minimizer->SetVariable(2, "rotZ", curAlignRot.Z(), 0.001);
  } else {
    minimizer->SetVariable(0, "dx", curAlignPos.X(), 0.001);
    minimizer->SetVariable(1, "dy", curAlignPos.Y(), 0.001);
    minimizer->SetVariable(2, "dz", curAlignPos.Z(), 0.001);
    minimizer->SetVariable(3, "rotX", curAlignRot.X(), 0.0001);
    minimizer->SetVariable(4, "rotY", curAlignRot.Y(), 0.0001);
    minimizer->SetVariable(5, "rotZ", curAlignRot.Z(), 0.0001);
    minimizer->SetVariableLimits(0,curAlignPos.X()-2,curAlignPos.X()+2);
    minimizer->SetVariableLimits(1,curAlignPos.Y()-2,curAlignPos.Y()+2);
    minimizer->SetVariableLimits(2,curAlignPos.Z()-20,curAlignPos.Z()+20);
    // minimizer->SetFixedVariable(2, "dz", curAlignPos.Z());
    minimizer->SetVariableLimits(3, curAlignRot.X()-0.1, curAlignRot.X()+0.1);
    minimizer->SetVariableLimits(4, curAlignRot.Y()-0.1, curAlignRot.Y()+0.1);
    minimizer->SetVariableLimits(5, curAlignRot.Z()-0.1, curAlignRot.Z()+0.1);
  }
  minimizer->Minimize();

  const double *result = minimizer->X();
  if (is3D) {
    double dx = result[0], dy = result[1], rotZ = result[2];
    if (!std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(rotZ)) { delete minimizer; return; }
    detector->SetAlignment(dx, dy, 0, 0, 0, rotZ);
    std::cout << "[DUT " << detID << "] Cov coarse: dx=" << dx << " dy=" << dy << " rotZ=" << rotZ
              << "  chi2=" << minimizer->MinValue() << std::endl;
  } else {
    double dx = result[0], dy = result[1], dz = result[2], rotX = result[3], rotY = result[4], rotZ = result[5];
    if (!std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dz) ||
        !std::isfinite(rotX) || !std::isfinite(rotY) || !std::isfinite(rotZ)) { delete minimizer; return; }
    detector->SetAlignment(dx, dy, dz, rotX, rotY, rotZ);
    std::cout << "[DUT " << detID << "] Cov fine: dx=" << dx << " dy=" << dy << " rotZ=" << rotZ
              << "  chi2=" << minimizer->MinValue() << std::endl;
  }
  delete minimizer;
}

double DeConvChi2Cov_3D(const double *par, const std::vector<Event> &events,
                         std::shared_ptr<Detector> detector, int detID) {
  detector->SetAlignment(par[0], par[1], 0, 0, 0, par[2]);
  double chi2 = 0;
  int n = 0;
  auto planar = dynamic_cast<Planar*>(detector.get());
  if (!planar) return 1e9;

  for (const auto &evt : events) {
    auto it = evt.detectorFramesMap.find(detID);
    if (it == evt.detectorFramesMap.end()) continue;
    const auto &clusters = it->second->ClustersDeconv();
    auto [globalHits, localHits] = planar->CalcHitsFromTrack(evt.track, 5, thickness);
    if ((int)localHits.size() != 5) continue;

    for (const auto &c : clusters) {
      bool complete = true;
      for (int k = 0; k < 5; k++)
        if (!std::isfinite(c.pos[k])) { complete = false; break; }
      if (!complete) continue;

      Eigen::VectorXd res(5);
      int t = c.type;
      for (int k = 0; k < 5; k++)
        res(k) = c.pos[k] * pitchwidth[t] - ((t == 0) ? localHits[k][0] : localHits[k][1]);
      chi2 += res.transpose() * s_tier_cov_inv[t] * res;
      n++;
    }
  }
  return n > 0 ? chi2 : 1e9;
}

double DeConvChi2Cov_6D(const double *par, const std::vector<Event> &events,
                         std::shared_ptr<Detector> detector, int detID) {
  detector->SetAlignment(par[0], par[1], par[2], par[3], par[4], par[5]);
  double chi2 = 0;
  int n = 0;
  auto planar = dynamic_cast<Planar*>(detector.get());
  if (!planar) return 1e9;

  for (const auto &evt : events) {
    auto it = evt.detectorFramesMap.find(detID);
    if (it == evt.detectorFramesMap.end()) continue;
    const auto &clusters = it->second->ClustersDeconv();
    auto [globalHits, localHits] = planar->CalcHitsFromTrack(evt.track, 5, thickness);
    if ((int)localHits.size() != 5) continue;

    for (const auto &c : clusters) {
      bool complete = true;
      for (int k = 0; k < 5; k++)
        if (!std::isfinite(c.pos[k])) { complete = false; break; }
      if (!complete) continue;

      Eigen::VectorXd res(5);
      int t = c.type;
      for (int k = 0; k < 5; k++)
        res(k) = c.pos[k] * pitchwidth[t] - ((t == 0) ? localHits[k][0] : localHits[k][1]);
      chi2 += res.transpose() * s_tier_cov_inv[t] * res;
      n++;
    }
  }
  return n > 0 ? chi2 : 1e9;
}

double DeConvTrackChi2_3D(const double *par, const std::vector<Event> &events, std::shared_ptr<Detector> detector, int detID,double reslimit) {

  const double dx = par[0];
  const double dy = par[1];
  const double rotZ = par[2];

  detector->SetAlignment(dx, dy, 0, 0, 0, rotZ);

  // 对所有事件求平均 χ²（n_tier层位置）
  double chi2 = 0.0;
  int eventnum = 0;

  auto planar = dynamic_cast<Planar*>(detector.get());
  if (!planar) return 1e9;

  for (const auto &evt : events) {
    auto frameIt = evt.detectorFramesMap.find(detID);
    if (frameIt == evt.detectorFramesMap.end() || frameIt->second->ClustersDeconv().empty()) {
      std::cout << "Event " << evt.eventID << " has no clusters for detector " << detID << std::endl;
      continue;
    }

    const auto &clusterdeconvs = frameIt->second->ClustersDeconv();

    // 直接用 Planar::CalcHitsFromTrack 获取各层的局域预测位置
    auto [globalHits, localHits] = planar->CalcHitsFromTrack(evt.track, n_tier, thickness);
    if ((int)localHits.size() != n_tier)
      continue;

    // 好cluster：每层独立计算残差，按type分别累加残差平方和
    double goodSumRes2[2] = {0.0, 0.0};
    double eventChi2 = 0.0;
    int tierCount[2] = {0, 0};

    for (auto clusterdeconv : clusterdeconvs) {
      int type = clusterdeconv.type;
      for (int k = 0; k < n_tier; k++) {
        if (clusterdeconv.chargetier[k] <= 0 || !std::isfinite(clusterdeconv.pos[k]))
          continue;
        double mpos = clusterdeconv.pos[k] * pitchwidth[type];
        double predictedPos = (type == 0) ? localHits[k][0] : localHits[k][1];
        double res = mpos - predictedPos;
        if(abs(res)>reslimit){continue;}//剔除残差过大的
        tierCount[type]++;
        goodSumRes2[type] += res * res;
      }
      if(tierCount[type] == 0){
        continue;
      }
      goodSumRes2[type] = goodSumRes2[type]/tierCount[type];
      eventChi2 += goodSumRes2[type];
      eventnum++;
    }

    // 坏条cluster检查：比较同type所有有效层的残差平方和
    bool isBadCluster = false;
    double badSumRes2[2] = {0.0, 0.0};
    int badTierCount[2] = {0, 0};
    for (auto clusterdeconv : frameIt->second->ClustersDeconvBad()) {
      int type = clusterdeconv.type;
      for (int k = 0; k < n_tier; k++) {
        if (clusterdeconv.chargetier[k] <= 0 || !std::isfinite(clusterdeconv.pos[k]))
          continue;
        double mpos = clusterdeconv.pos[k] * pitchwidth[type];
        double predictedPos = (type == 0) ? localHits[k][0] : localHits[k][1];
        double res = mpos - predictedPos;
        if(abs(res)>10){continue;}//剔除残差过大的
        badTierCount[type]++;
        badSumRes2[type] += res * res;
      }
      if(badTierCount[type] == 0){
        continue;
      }
      badSumRes2[type] = badSumRes2[type]/badTierCount[type];
      if (badSumRes2[type] > 0 && badSumRes2[type] < goodSumRes2[type]) {
        isBadCluster = true;
        break;
      }
    }

    if (!isBadCluster && tierCount[0] > 0 && tierCount[1] > 0) {
      chi2 += eventChi2;
    }
  }

  // 返回平均 χ²（按总层数归一化）
  return (eventnum > 0) ? chi2 / eventnum : 1e9;
}

double DeConvTrackChi2_6D(const double *par, const std::vector<Event> &events, std::shared_ptr<Detector> detector, int detID,double reslimit) {

  const double dx = par[0];
  const double dy = par[1];
  const double dz = par[2];
  const double rotX = par[3];
  const double rotY = par[4];
  const double rotZ = par[5];

  detector->SetAlignment(dx, dy, dz, rotX, rotY, rotZ);

  // 对所有事件求平均 χ²（5层位置）
  double chi2 = 0.0;
  int eventnum = 0;

  auto planar = dynamic_cast<Planar*>(detector.get());
  if (!planar) return 1e9;

  for (const auto &evt : events) {
    auto frameIt = evt.detectorFramesMap.find(detID);
    if (frameIt == evt.detectorFramesMap.end() || frameIt->second->ClustersDeconv().empty()) {
      continue;
    }

    const auto &clusterdeconvs = frameIt->second->ClustersDeconv();

    // 直接用 Planar::CalcHitsFromTrack 获取各层的局域预测位置
    auto [globalHits, localHits] = planar->CalcHitsFromTrack(evt.track, n_tier, thickness);
    if ((int)localHits.size() != n_tier)
      continue;

    // 好cluster：每层独立计算残差，按type分别累加残差平方和
    double goodSumRes2[2] = {0.0, 0.0};
    double eventChi2 = 0.0;
    int tierCount[2] = {0, 0};

    for (auto clusterdeconv : clusterdeconvs) {
      int type = clusterdeconv.type;
      for (int k = 0; k < n_tier; k++) {
        if (clusterdeconv.chargetier[k] <= 0 || !std::isfinite(clusterdeconv.pos[k]))
          continue;
        double mpos = clusterdeconv.pos[k] * pitchwidth[type];
        double predictedPos = (type == 0) ? localHits[k][0] : localHits[k][1];
        double res = mpos - predictedPos;
        if(abs(res)>reslimit){continue;}//剔除残差过大的
        tierCount[type]++;
        goodSumRes2[type] += res * res;
      }
      if(tierCount[type] == 0){
        continue;
      }
      goodSumRes2[type] = goodSumRes2[type]/tierCount[type];
      eventChi2 += goodSumRes2[type];
      eventnum++;
    }

    // 坏条cluster检查：比较同type所有有效层的残差平方和
    bool isBadCluster = false;
    double badSumRes2[2] = {0.0, 0.0};
    int badTierCount[2] = {0, 0};
    for (auto clusterdeconv : frameIt->second->ClustersDeconvBad()) {
      int type = clusterdeconv.type;
      for (int k = 0; k < n_tier; k++) {
        if (clusterdeconv.chargetier[k] <= 0 || !std::isfinite(clusterdeconv.pos[k]))
          continue;
        double mpos = clusterdeconv.pos[k] * pitchwidth[type];
        double predictedPos = (type == 0) ? localHits[k][0] : localHits[k][1];
        double res = mpos - predictedPos;
        if(abs(res)>10){continue;}//剔除残差过大的
        badTierCount[type]++;
        badSumRes2[type] += res * res;
      }
      if(badTierCount[type] == 0){
        continue;
      }
      badSumRes2[type] = badSumRes2[type]/badTierCount[type];
      if (badSumRes2[type] > 0 && badSumRes2[type] < goodSumRes2[type]) {
        isBadCluster = true;
        break;
      }
    }

    if (!isBadCluster && tierCount[0] > 0 && tierCount[1] > 0) {
      chi2 += eventChi2;
    }
  }

  // 返回平均 χ²（按总层数归一化）
  return (eventnum > 0) ? chi2 / eventnum : 1e9;
}

// ================================================================
// 微 TPC 对齐目标函数（用 microTPCposition，非逐层 pos[k]）
// ================================================================

// 三维对齐：dx, dy, rotZ
double DeConvPositionChi2_3D(const double *par, const std::vector<Event> &events,
                              std::shared_ptr<Detector> detector, int detID,double reslimit) {
  const double dx = par[0];
  const double dy = par[1];
  const double rotZ = par[2];

  detector->SetAlignment(dx, dy, 0, 0, 0, rotZ);

  double chi2 = 0.0;
  int eventnum = 0;

  for (const auto &evt : events) {
    auto frameIt = evt.detectorFramesMap.find(detID);
    if (frameIt == evt.detectorFramesMap.end() ||
        frameIt->second->ClustersDeconv().empty())
      continue;

    const auto &clusterdeconvs = frameIt->second->ClustersDeconv();

    // track 在 z_ref_position 处的局域预测位置
    TVector3 predG = detector->CalcHitFromTrack(evt.track);
    TVector3 predL = detector->GlobalToLocal(predG);
    TVector3 trackDirGlobal(evt.track.kx, evt.track.ky, 1);
    TVector3 localDir = detector->GlobalDirToLocal(trackDirGlobal);
    double predictedPos[2];
    predictedPos[0] = (fabs(localDir.Z()) > 1e-9)
                          ? predL[0] + localDir[0] / localDir.Z() * z_ref_position
                          : predL[0];
    predictedPos[1] = (fabs(localDir.Z()) > 1e-9)
                          ? predL[1] + localDir[1] / localDir.Z() * z_ref_position
                          : predL[1];

    // 好cluster：microTPCposition vs track 预测位置
    double goodSumRes2[2] = {0.0, 0.0};
    double eventChi2 = 0.0;
    int clusterCount[2] = {0, 0};

    for (auto clusterdeconv : clusterdeconvs) {
      int type = clusterdeconv.type;
      if (!std::isfinite(clusterdeconv.chargeposition))
        continue;
      double mpos = clusterdeconv.microTPCposition * pitchwidth[type];
      double res = mpos - predictedPos[type];
      if (fabs(res) > reslimit) continue;
      clusterCount[type]++;
      goodSumRes2[type] += res * res;
    }
    if (clusterCount[0] > 0) {
      goodSumRes2[0] /= clusterCount[0];
      eventChi2 += goodSumRes2[0];
      eventnum++;
    }
    if (clusterCount[1] > 0) {
      goodSumRes2[1] /= clusterCount[1];
      eventChi2 += goodSumRes2[1];
      eventnum++;
    }

    // 坏条cluster检查
    bool isBadCluster = false;
    double badSumRes2[2] = {0.0, 0.0};
    int badClusterCount[2] = {0, 0};
    for (auto clusterdeconv : frameIt->second->ClustersDeconvBad()) {
      int type = clusterdeconv.type;
      if (!std::isfinite(clusterdeconv.microTPCposition))
        continue;
      double mpos = clusterdeconv.microTPCposition * pitchwidth[type];
      double res = mpos - predictedPos[type];
      badClusterCount[type]++;
      badSumRes2[type] += res * res;
    }
    for (int type = 0; type < 2; ++type) {
      if (badClusterCount[type] > 0) {
        badSumRes2[type] /= badClusterCount[type];
        if (badSumRes2[type] > 0 && badSumRes2[type] < goodSumRes2[type]) {
          isBadCluster = true;
          break;
        }
      }
    }

    if (!isBadCluster && clusterCount[0] > 0 && clusterCount[1] > 0) {
      chi2 += eventChi2;
    }
  }

  return (eventnum > 0) ? chi2 : 1e9;
}

// 六维对齐：dx, dy, dz, rotX, rotY, rotZ
double DeConvPositionChi2_6D(const double *par, const std::vector<Event> &events,
                              std::shared_ptr<Detector> detector, int detID,double reslimit) {
  const double dx = par[0];
  const double dy = par[1];
  const double dz = par[2];
  const double rotX = par[3];
  const double rotY = par[4];
  const double rotZ = par[5];

  detector->SetAlignment(dx, dy, dz, rotX, rotY, rotZ);

  double chi2 = 0.0;
  int eventnum = 0;

  for (const auto &evt : events) {
    auto frameIt = evt.detectorFramesMap.find(detID);
    if (frameIt == evt.detectorFramesMap.end() ||
        frameIt->second->ClustersDeconv().empty())
      continue;

    const auto &clusterdeconvs = frameIt->second->ClustersDeconv();

    TVector3 predG = detector->CalcHitFromTrack(evt.track);
    TVector3 predL = detector->GlobalToLocal(predG);
    TVector3 trackDirGlobal(evt.track.kx, evt.track.ky, 1);
    TVector3 localDir = detector->GlobalDirToLocal(trackDirGlobal);
    double predictedPos[2];
    predictedPos[0] = (fabs(localDir.Z()) > 1e-9)
                          ? predL[0] + localDir[0] / localDir.Z() * z_ref_position
                          : predL[0];
    predictedPos[1] = (fabs(localDir.Z()) > 1e-9)
                          ? predL[1] + localDir[1] / localDir.Z() * z_ref_position
                          : predL[1];

    double goodSumRes2[2] = {0.0, 0.0};
    double eventChi2 = 0.0;
    int clusterCount[2] = {0, 0};

    for (auto clusterdeconv : clusterdeconvs) {
      int type = clusterdeconv.type;
      if (!std::isfinite(clusterdeconv.microTPCposition))
        continue;
      double mpos = clusterdeconv.microTPCposition * pitchwidth[type];
      double res = mpos - predictedPos[type];
      if (fabs(res) > reslimit) continue;
      clusterCount[type]++;
      goodSumRes2[type] += res * res;
    }
    if (clusterCount[0] > 0) {
      goodSumRes2[0] /= clusterCount[0];
      eventChi2 += goodSumRes2[0];
      eventnum++;
    }
    if (clusterCount[1] > 0) {
      goodSumRes2[1] /= clusterCount[1];
      eventChi2 += goodSumRes2[1];
      eventnum++;
    }

    bool isBadCluster = false;
    double badSumRes2[2] = {0.0, 0.0};
    int badClusterCount[2] = {0, 0};
    for (auto clusterdeconv : frameIt->second->ClustersDeconvBad()) {
      int type = clusterdeconv.type;
      if (!std::isfinite(clusterdeconv.microTPCposition))
        continue;
      double mpos = clusterdeconv.microTPCposition * pitchwidth[type];
      double res = mpos - predictedPos[type];
      badClusterCount[type]++;
      badSumRes2[type] += res * res;
    }
    for (int type = 0; type < 2; ++type) {
      if (badClusterCount[type] > 0) {
        badSumRes2[type] /= badClusterCount[type];
        if (badSumRes2[type] > 0 && badSumRes2[type] < goodSumRes2[type]) {
          isBadCluster = true;
          break;
        }
      }
    }

    if (!isBadCluster && clusterCount[0] > 0 && clusterCount[1] > 0) {
      chi2 += eventChi2;
    }
  }

  return (eventnum > 0) ? chi2  : 1e9;
}

std::map<int,std::pair<double,double>> Positionanalysis(std::vector<Event> events,std::string mode,TDirectory* fout){
  TF1 *f_res[2];
  TH1D *h_res[2];
  std::map<int,std::pair<double,double>> res_map;
  for(int i = 0;i<2;i++){
    f_res[i] = new TF1(Form("f_res_%s",i==0?"X":"Y"),"gaus",-5,5);
    h_res[i] = new TH1D(Form("h_res_%s",i==0?"X":"Y"),Form("Residual distribution(mode: %s);Residual (mm);Entries",mode.c_str()),1000,-5,5);
  }
   for (auto event : events) {
    for (const auto &framePair : event.detectorFramesMap) {
      int detID = framePair.first;
      auto &detEvt = framePair.second;
      auto &clusterDeconvs = detEvt->ClustersDeconv();
      auto&factory = DetectorFactory::GetInstance();
      // 使用已对齐的 detector，而非 DetectorFrame 中的拷贝
      auto alignedDet = factory.GetDetector(detID);
      if (!alignedDet) continue;
      const Detector &det = *alignedDet;
      TVector3 predG = det.CalcHitFromTrack(event.track);
      TVector3 predL = det.GlobalToLocal(predG);
      TVector3 trackDirGlobal(event.track.kx, event.track.ky, 1);
      TVector3 localDir = det.GlobalDirToLocal(trackDirGlobal);

      for (auto &cluster : clusterDeconvs) {
        double hit,res;
        if (fabs(localDir.Z()) > 1e-9) {
          hit = predL[cluster.type]+ localDir[cluster.type] / localDir.Z() * z_ref_position; // 考虑入射角的影响
        } else {
          hit = predL[cluster.type];
        }
        if(mode == "charge"){
           res = cluster.chargeposition * pitchwidth[cluster.type] - hit;
        }else if(mode == "microTPC"){
           res = cluster.microTPCposition * pitchwidth[cluster.type] - hit;
        }

        h_res[cluster.type]->Fill(res);

      }
    }
  }
  for(int i = 0;i<2;i++){
    h_res[i]->Fit(f_res[i],"Q");
    double mean = f_res[i]->GetParameter(1);
    double sigma = f_res[i]->GetParameter(2);
    std::cout<<"Residual "<<(i==0?"X":"Y")<<" resolution (mode: "<<mode<<"): "<<sigma<<" mm"<<std::endl;
    fout->WriteTObject(h_res[i]);
    res_map[i] = {mean, sigma};
    delete f_res[i];
    delete h_res[i];
  }
  return res_map;
}


// 注册脚本
REGISTER_SCRIPT("DeconvAnalysisScript", DeconvAnalysisScript);
