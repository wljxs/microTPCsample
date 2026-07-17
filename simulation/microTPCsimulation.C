// ============================================================
// microTPCsimulation.C — 基于抽样库的 Track 信号模拟（20° 入射角）,不过每一个都是从每层的中心点发射的
//
//   流程：
//     1. TrackHeed 生成 20° 入射角的 muon track
//     2. 对 track 中每个 cluster 的每个初级电子，漂移到读出平面
//     3. 按漂移终点 (ex1, ey1) 在抽样库中查表获取信号
//     4. 聚合所有电子信号，输出 stripid + 波形 + track 斜率/截距
//
//   依赖：需先运行 analysis/sampleconstruction.C 中的 buildSampleLibrary()
// ============================================================

#include <Rtypes.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TSystem.h>
#include <TTree.h>
#include <TF1.h>
#include <TRandom3.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <limits>
#include <random>
#include <string>

#include "Garfield/ComponentComsol.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/RandomEngine.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/Shaper.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/TrackHeed.hh"

using namespace Garfield;

// ---------- 抽样库元数据（从 samplelibrary.root 的 meta tree 读取） ----------
static double gLib_xmin, gLib_xmax, gLib_ymin, gLib_ymax;
static int    gLib_nbins;
static std::vector<double> gLib_t;  // 共享时间轴（从 meta 读取一次）
static const int kNbins = 10;       // 与 buildSampleLibrary.C 中 kNbins 一致

// ---- 信号通道名称表（14 个，与 microTPC.C 一致）----
static const int N_SIG = 8;
static const char* SIG_NAME[N_SIG] = {
        "sigxtestminus1","sigxtest1", 
    "sigxtestminus2",  "sigzero","sigxtest2",   
     "sigytestminus1", "sigy","sigytest1"
};//前4个信号是x方向落到两个条中间时相邻4个条产生的，中间5个信号是x方向落到条上时相邻5个条产生的，后5个信号是y方向上相邻5个条产生的

// ============================================================
// main — 主模拟程序
// ============================================================
int main(int argc, char *argv[])
{
    // ---------- 参数 ----------
    int runid = (argc > 1) ? std::atoi(argv[1]) : 0;
    const std::string outputFile = (argc > 2)
        ? argv[2] : "../result/nodelta/20test.root";

    std::string detector = "microTPC5mmone";
    constexpr double pitch = 0.04;
    const unsigned int nEvents = (argc > 3)
        ? static_cast<unsigned int>(std::max(1, std::atoi(argv[3]))) : 5000;
    constexpr int N_STRIPS = 20;                 // x/y 方向各 20 条
    constexpr double trackAngle = 20;         // 入射角（度，相对于垂直方向）
    constexpr double trackAngleRad = trackAngle * M_PI / 180.0;
    TRandom3 rnd(runid * 12345);

    // ---------- 加载抽样库 ----------
    const char *libFile = "../analysis/samplelibrarysmall.root";
    TFile *fLib = TFile::Open(libFile, "READ");
    if (!fLib || fLib->IsZombie()) {
        std::cerr << "Error: Cannot open " << libFile << std::endl;
        return 1;
    }

    // —— 读取分箱元数据 + 共享时间轴 ——
    TTree *meta = (TTree *)fLib->Get("meta");
    if (!meta) {
        std::cerr << "Error: No meta tree found in " << libFile << std::endl;
        fLib->Close();
        return 1;
    }
    {
        double xmin, xmax, ymin, ymax;
        int    nbins;
        std::vector<double> *meta_t = nullptr;
        meta->SetBranchAddress("xmin",  &xmin);
        meta->SetBranchAddress("xmax",  &xmax);
        meta->SetBranchAddress("ymin",  &ymin);
        meta->SetBranchAddress("ymax",  &ymax);
        meta->SetBranchAddress("nbins", &nbins);
        meta->SetBranchAddress("t",     &meta_t);
        meta->GetEntry(0);
        gLib_xmin  = xmin;
        gLib_xmax  = xmax;
        gLib_ymin  = ymin;
        gLib_ymax  = ymax;
        gLib_nbins = nbins;
        gLib_t     = *meta_t;

        if (gLib_nbins != kNbins) {
            std::cerr << "Error: sample library nbins=" << gLib_nbins
                      << " does not match compiled kNbins=" << kNbins << std::endl;
            fLib->Close();
            return 1;
        }

        std::cout << "[loadSampleLibrary] Binning: nbins=" << nbins
                  << "  x=[" << xmin << ", " << xmax << "]"
                  << "  y=[" << ymin << ", " << ymax << "]"
                  << "  t.size=" << gLib_t.size() << std::endl;
    }

    // —— 逐一打开各区域的独立树，预绑定分支 ——
    TTree *regionTrees[kNbins][kNbins] = {};
    std::vector<double> *region_sig[kNbins][kNbins][N_SIG] = {};
    Long64_t region_n[kNbins][kNbins] = {};
    {
        int nNonEmpty = 0;
        for (int ix = 0; ix < gLib_nbins; ++ix) {
            for (int iy = 0; iy < gLib_nbins; ++iy) {
                std::string tname = "region_" + std::to_string(ix) + "_" + std::to_string(iy);
                TTree *rt = (TTree *)fLib->Get(tname.c_str());
                if (!rt) continue;

                regionTrees[ix][iy] = rt;
                region_n[ix][iy]   = rt->GetEntries();
                for (int is = 0; is < N_SIG; ++is)
                    rt->SetBranchAddress(SIG_NAME[is], &region_sig[ix][iy][is]);

                ++nNonEmpty;
                if (nNonEmpty <= 10)
                    std::cout << "  region (" << ix << "," << iy << "): "
                              << region_n[ix][iy] << " samples" << std::endl;
            }
        }
        std::cout << "  Total non-empty regions: " << nNonEmpty
                  << " / " << gLib_nbins * gLib_nbins << std::endl;
    }
    // NOTE: fLib 不关闭，各 region tree 需要它保持打开

    // ---------- 气体 ----------
    MediumMagboltz gas;
    gas.LoadGasFile("../" + detector + "/ar_45_cf4_40_co2_15100-100000.gas");
    constexpr double rPenning = 0.20;
    constexpr double lambdaPenning = 0.;
    gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");
    gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");
    gas.SetMaxElectronEnergy(200.);

    // ---------- 电场 ----------
    ComponentComsol fm;
    fm.Initialise(("../" + detector + "/mesh.mphtxt").c_str(),
                  ("../" + detector + "/mplist.txt").c_str(),
                  ("../" + detector + "/field.txt").c_str(), "mm");
    fm.EnablePeriodicityX();
    fm.EnablePeriodicityY();
    fm.SetGas(&gas);

    // ---------- Sensor ----------
    Sensor sensor(&fm);
    sensor.SetArea(0, 0, 0.01, N_STRIPS * pitch, N_STRIPS * pitch, 0.51);


    const unsigned int nTimeBins = 140;
    const double tmin = 0.;
    const double tmax = 700;
    const double tstep = (tmax - tmin) / nTimeBins;
    sensor.SetTimeWindow(tmin, tstep, nTimeBins);
    Shaper shaper(1, 50., 1., "unipolar");
    sensor.SetTransferFunction(shaper);

    TF1 *kernel = new TF1("kernel", "[0]*exp([1])*pow(x/[1]/[2],[1])*exp(-x/[2])", tmin, tmax);
    kernel->SetParameter(0, 1.0);
    kernel->SetParameter(1, 1.0);
    kernel->SetParameter(2, 50.0);
    std::vector<double> kernel_vec(nTimeBins);
    for (int i = 0; i < nTimeBins; ++i) {
        double t = tmin + i * tstep;
        kernel_vec[i] = kernel->Eval(t);
    }
 
    // ---------- TrackHeed（20° 入射角 muon track） ---------- 
    TrackHeed track(&sensor);
    track.SetParticle("muon");
    track.SetEnergy(150.e9);  // 150 GeV muon

    // ---------- Drift ----------
    AvalancheMC drift(&sensor);
    drift.SetTimeWindow(0, 1400);
    drift.SetDistanceSteps(1.e-3);

    // ---------- 输出 Tree ----------
    TFile *fout = new TFile(outputFile.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << outputFile << std::endl;
        return 1;
    }

    int event;
    // 时间轴
    std::vector<double> vec_t;
    for(int i = 0; i < nTimeBins; ++i) vec_t.push_back(tmin + i * tstep);
    
    // 内部信号累加器（固定 N_STRIPS 条 x + N_STRIPS 条 y）
    std::vector<double> vec_sigx[N_STRIPS];
    std::vector<double> vec_sigy[N_STRIPS];
    std::vector<double> vec_sigx_conv[N_STRIPS];
    std::vector<double> vec_sigy_conv[N_STRIPS];
    
    // 每 event 的输出：strip ID 向量 + 波形向量 + track 斜率/截距
    std::vector<int>                  stripIDs;
    std::vector<std::vector<double>>  waveforms;
    std::vector<std::vector<double>>  waveforms_conv;
    std::vector<int>                  types;  // 0: x strip, 1: y strip
    double kx, bx;
    double ky, by;

    TTree *tout = new TTree("tree", "Track Sampled data");
    tout->Branch("event",         &event,          "event/I");
    tout->Branch("t",             &vec_t);
    tout->Branch("types",         &types);
    tout->Branch("stripIDs",     &stripIDs);
    tout->Branch("waveforms",     &waveforms);
    tout->Branch("waveforms_conv",&waveforms_conv);
    tout->Branch("kx",      &kx,       "kx/D");
    tout->Branch("bx",  &bx,   "bx/D");
    tout->Branch("ky",      &ky,       "ky/D");
    tout->Branch("by",  &by,   "by/D");

    // 每个 event 的电子终点位置（tree2，与 tree 一一对应）
    TTree *tpos = new TTree("tree2", "Electron endpoint positions per event");
    std::vector<double> vec_xe1, vec_ye1, vec_ze1, vec_te1, vec_e1;
    std::vector<int>    vec_xe_stripid, vec_ye_stripid;
    tpos->Branch("event",        &event,           "event/I");
    tpos->Branch("xe1",          &vec_xe1);
    tpos->Branch("ye1",          &vec_ye1);
    tpos->Branch("ze1",          &vec_ze1);
    tpos->Branch("te1",          &vec_te1);
    tpos->Branch("e1",           &vec_e1);
    tpos->Branch("xe_stripid",   &vec_xe_stripid);
    tpos->Branch("ye_stripid",   &vec_ye_stripid);

    // 每个 event 的电子初始电离位置（tree3，与 tree 一一对应）
    TTree *tinit = new TTree("tree3", "Electron initial ionization positions per event");
    std::vector<double> vec_xe0, vec_ye0, vec_ze0, vec_te0, vec_ee0;
    std::vector<int>    vec_xe0_stripid, vec_ye0_stripid;
    tinit->Branch("event",        &event,           "event/I");
    tinit->Branch("xe0",          &vec_xe0);
    tinit->Branch("ye0",          &vec_ye0);
    tinit->Branch("ze0",          &vec_ze0);
    tinit->Branch("te0",          &vec_te0);
    tinit->Branch("ee0",          &vec_ee0);
    tinit->Branch("xe0_stripid",  &vec_xe0_stripid);
    tinit->Branch("ye0_stripid",  &vec_ye0_stripid);

    // 统计漂移终点到抽样库的接受率，便于检查 (ix, iy) 是否越界。
    Long64_t nPrimary = 0;
    Long64_t nDriftSuccess = 0;
    Long64_t nReadoutInRange = 0;
    Long64_t nLibraryOutOfRange = 0;
    Long64_t nEmptyRegion = 0;
    Long64_t nSampledElectrons = 0;
    // 性能统计：区分 Garfield 漂移、ROOT 模板读取/叠加和电子学卷积。
    double driftSeconds = 0.;
    double librarySeconds = 0.;
    double convolutionSeconds = 0.;
    const auto eventLoopStart = std::chrono::steady_clock::now();
    // ==================== 事件循环 ====================
    track.DisableDeltaElectronTransport();
    for (unsigned int i = 0; i < nEvents; ++i) {

        // —— 重置信号累加器 ——
        for (int indexx = 0; indexx < N_STRIPS; ++indexx) {
            vec_sigx[indexx].assign(nTimeBins, 0);
            vec_sigy[indexx].assign(nTimeBins, 0);
            vec_sigx_conv[indexx].assign(nTimeBins, 0);
            vec_sigy_conv[indexx].assign(nTimeBins, 0);
        }
        
        if(i % 100 == 0) std::cout << "Event: " << i << "/" << nEvents << std::endl;
        sensor.ClearSignal();

        // —— 清空本 event 的电子位置向量 ——
        vec_xe1.clear(); vec_ye1.clear(); vec_ze1.clear(); vec_te1.clear(); vec_e1.clear();
        vec_xe_stripid.clear(); vec_ye_stripid.clear();
        vec_xe0.clear(); vec_ye0.clear(); vec_ze0.clear(); vec_te0.clear(); vec_ee0.clear();
        vec_xe0_stripid.clear(); vec_ye0_stripid.clear();

        // —— 随机 track 起始位置 ——
        double x0 = 2*pitch + rnd.Uniform(0*pitch, 2 * pitch);
        double y0 = 2*pitch + rnd.Uniform(0*pitch, 2 * pitch);
        double z0 = 0.49;   // 从探测器顶部附近开始
        double t0 = 0.;
        
        // —— track 方向：入射角 20°（相对于垂直 z 方向，在 y-z 平面内） ——
        double dx = 0.;
        double dy = std::sin(trackAngleRad);
        double dz = -std::cos(trackAngleRad);  // 向下
        
        // —— 计算斜率与截距（z 为自变量） ——
        kx     = dx / dz;                     // = -tan(angle)
        bx = x0 - kx * z0;
        ky     = dy/dz;                          // dy=0，无 y 方向偏转
        by = y0 - ky *z0;
        
        event = i + runid * nEvents;

        // —— 生成 track ——
        track.NewTrack(x0, y0, z0, t0, dx, dy, dz);
        
        // —— 获取 clusters ——
        double xc, yc, zc, tc, ec, extra;
        int ne;
        
        // 遍历每个 cluster
        while (track.GetCluster(xc, yc, zc, tc, ne, ec, extra)) {
            // 遍历 cluster 中的每个初级电子
            for (int ie = 0; ie < ne; ++ie) {
                double xe0, ye0, ze0, te0, ee0,dxe, dye, dze;
                track.GetElectron(ie, xe0, ye0, ze0, te0, ee0, dxe, dye, dze);
                ++nPrimary;
                
                // —— 记录初始电离位置（tree3） ——
                vec_xe0.push_back(xe0);
                vec_ye0.push_back(ye0);
                vec_ze0.push_back(ze0);
                vec_te0.push_back(te0);
                vec_ee0.push_back(ee0);
                vec_xe0_stripid.push_back((int)(xe0 / pitch));
                vec_ye0_stripid.push_back((int)(ye0 / pitch));
                
                // 漂移模拟
                const auto driftStart = std::chrono::steady_clock::now();
                drift.DriftElectron(xe0, ye0, ze0, te0, ee0);
                driftSeconds += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - driftStart).count();
                
                int status;
                double xe1, ye1, ze1, te1;
                drift.GetElectronEndpoint(0, xe0, ye0, ze0, te0,
                                          xe1, ye1, ze1, te1, status);
                
                if(status != -1) {
                    // 漂移失败，跳过
                    continue;
                }
                ++nDriftSuccess;
                
                // 计算 strip 编号
                int xe_stripid = (int)(xe1 / pitch);
                int ye_stripid = (int)(ye1 / pitch);
                
                // 边界检查
                if(xe_stripid < 0 || xe_stripid >= N_STRIPS ||
                   ye_stripid < 0 || ye_stripid >= N_STRIPS) continue;
                ++nReadoutInRange;
                
                // 记录每个电子的终点（tree2）
                vec_xe1.push_back(xe1);
                vec_ye1.push_back(ye1);
                vec_ze1.push_back(ze1);
                vec_te1.push_back(te1);
                vec_e1.push_back(ee0);
                vec_xe_stripid.push_back(xe_stripid);
                vec_ye_stripid.push_back(ye_stripid);
                
                // —— 区域映射（与 microTPCsample.C 一致） ——
                double deltax = gLib_xmax - gLib_xmin;
                double binWx = (gLib_xmax - gLib_xmin) / gLib_nbins;
                double binWy = (gLib_ymax - gLib_ymin) / gLib_nbins;
                int ix, iy;
                bool iscenter = false;
                int isleft = 0;
                int stripidx = (int)(xe1 / pitch);
                int signaljudege = (int)(xe1 / (pitch/4)) % 4;
                int stripidy = (int)(ye1 / pitch);
                
                if(signaljudege == 0){
                    ix = (int)((std::fmod(xe1, pitch) + deltax/2) / binWx);
                    isleft = 1;
                }
                else if(signaljudege == 1){
                    ix = (int)((std::fmod(xe1, pitch) - 0.5*pitch + deltax/2) / binWx);
                    iscenter = true;
                }
                else if(signaljudege == 2){
                    ix = (int)((std::fmod(xe1, pitch) - 0.5*pitch + deltax/2) / binWx);
                    iscenter = true;
                }
                else{
                    ix = (int)((std::fmod(xe1, pitch) - pitch + deltax/2) / binWx);
                }
                iy = (int)((std::fmod(ye1, pitch)) / binWy);
                
                // 边界检查
                if(ix < 0 || ix >= gLib_nbins || iy < 0 || iy >= gLib_nbins) {
                    ++nLibraryOutOfRange;
                    continue;
                }
                
                // —— 从抽样库获取信号 ——
                if (!regionTrees[ix][iy] || region_n[ix][iy] == 0) {
                    ++nEmptyRegion;
                    continue;
                }
                
                const auto libraryStart = std::chrono::steady_clock::now();
                Long64_t entry_choose = rnd.Integer(region_n[ix][iy]);
                regionTrees[ix][iy]->GetEntry(entry_choose);
                
                // 检查信号向量非空
                bool sigEmpty = false;
                if (iscenter) {
                    for (int is = 2; is <= 4; ++is) {
                        if (!region_sig[ix][iy][is] || region_sig[ix][iy][is]->empty()) {
                            sigEmpty = true; break;
                        }
                    }
                } else {
                    for (int is = 0; is <= 1; ++is) {
                        if (!region_sig[ix][iy][is] || region_sig[ix][iy][is]->empty()) {
                            sigEmpty = true; break;
                        }
                    }
                }
                for (int is = 5; is <= 7; ++is) {
                    if (!region_sig[ix][iy][is] || region_sig[ix][iy][is]->empty()) {
                        sigEmpty = true; break;
                    }
                }
                if (sigEmpty) continue;
                
                // 浮点索引 + 线性插值，解决 te1 不是 tstep 整数倍时的截断误差
                double fracBegin = te1 / tstep;
                int    beginindex = (int)std::floor(fracBegin);
                double alpha      = fracBegin - beginindex;
                double w0 = 1.0 - alpha, w1 = alpha;

                if (beginindex < 0 || beginindex >= static_cast<int>(nTimeBins)) continue;

                // 所有实际使用的通道必须具有同一段可插值的有效长度。
                int libSize = std::numeric_limits<int>::max();
                const int xFirst = iscenter ? 2 : 0;
                const int xLast = iscenter ? 4 : 1;
                for (int is = xFirst; is <= xLast; ++is) {
                    libSize = std::min(libSize, static_cast<int>(region_sig[ix][iy][is]->size()));
                }
                for (int is = 5; is <= 7; ++is) {
                    libSize = std::min(libSize, static_cast<int>(region_sig[ix][iy][is]->size()));
                }
                if (libSize < 2) continue;
                ++nSampledElectrons;

                for (int k = beginindex; k < nTimeBins; ++k) {
                    int j = k - beginindex;               // 库内索引
                    if (j + 1 >= libSize) break;          // 越界保护：需要 j 和 j+1 两个样本

                    if(iscenter){
                        int signalidx = 2;
                        for(int indexx = stripidx-1; indexx <= stripidx+1; ++indexx) {
                            double val = w0 * region_sig[ix][iy][signalidx]->at(j)
                                    + w1 * region_sig[ix][iy][signalidx]->at(j+1);
                            if (indexx >= 0 && indexx < N_STRIPS) vec_sigx[indexx][k] += val;
                            signalidx++;
                        }
                    }
                    else{
                        int signalidx = 0;
                        for(int indexx = stripidx-isleft; indexx <= stripidx-isleft+1; ++indexx) {
                            double val = w0 * region_sig[ix][iy][signalidx]->at(j)
                                    + w1 * region_sig[ix][iy][signalidx]->at(j+1);
                            if (indexx >= 0 && indexx < N_STRIPS) vec_sigx[indexx][k] += val;
                            signalidx++;
                        }
                    }
                    int signalidy = 5;
                    for(int indexy = stripidy-1; indexy <= stripidy+1; ++indexy) {
                        double val = w0 * region_sig[ix][iy][signalidy]->at(j)
                                + w1 * region_sig[ix][iy][signalidy]->at(j+1);
                        if (indexy >= 0 && indexy < N_STRIPS) vec_sigy[indexy][k] += val;
                        signalidy++;
                    }
                }
                librarySeconds += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - libraryStart).count();
            } // end electron loop
        } // end cluster loop
        
        // —— 卷积 ——
        const auto convolutionStart = std::chrono::steady_clock::now();
        for(int indexx = 0; indexx < N_STRIPS; ++indexx) {
            for (int k = 0; k < nTimeBins; ++k) {
                double convx = 0;
                double convy = 0;
                for (int j = 0; j <= k; ++j) {
                    int idx = k - j;
                    convx += vec_sigx[indexx][idx] * kernel_vec[j];
                    convy += vec_sigy[indexx][idx] * kernel_vec[j];
                }
                vec_sigx_conv[indexx][k] = convx * tstep;
                vec_sigy_conv[indexx][k] = convy * tstep;
            }
        }
        convolutionSeconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - convolutionStart).count();
        
        // —— 转换为输出格式：收集非零信号的 strip ——
        types.clear();
        stripIDs.clear();
        waveforms.clear();
        waveforms_conv.clear();
        for(int idx = 0; idx < N_STRIPS; ++idx) {
            // 检查 x 方向是否有信号
            bool hasSignalX = false;
            for(int k = 0; k < nTimeBins; ++k) {
                if(vec_sigx_conv[idx][k] != 0) { hasSignalX = true; break; }
            }
            if(hasSignalX) {
                stripIDs.push_back(idx);           // x strip: 0..N_STRIPS-1
                types.push_back(0);  // 标记为 x strip
                waveforms.push_back(vec_sigx[idx]);
                waveforms_conv.push_back(vec_sigx_conv[idx]);
            }
        }
        for(int idy = 0; idy < N_STRIPS; ++idy) {
            // 检查 y 方向是否有信号
            bool hasSignalY = false;
            for(int k = 0; k < nTimeBins; ++k) {
                if(vec_sigy_conv[idy][k] != 0) { hasSignalY = true; break; }
            }
            if(hasSignalY) {
                stripIDs.push_back(idy);      // y strip: N_STRIPS..2*N_STRIPS-1
                types.push_back(1);  // 标记为 y strip
                waveforms.push_back(vec_sigy[idy]);
                waveforms_conv.push_back(vec_sigy_conv[idy]);
            }
        }
        
        tpos->Fill();
        tinit->Fill();
        tout->Fill();
    }

    // 写入输出
    tout->Write();
    tpos->Write();
    tinit->Write();
    fout->Close();
    fLib->Close();
    delete kernel;

    const double eventLoopSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - eventLoopStart).count();

    std::cout << "\nSample-library acceptance:" << std::endl;
    std::cout << "  primary electrons: " << nPrimary << std::endl;
    std::cout << "  drift successful: " << nDriftSuccess << std::endl;
    std::cout << "  readout in range: " << nReadoutInRange << std::endl;
    std::cout << "  library (ix, iy) out of range: " << nLibraryOutOfRange << std::endl;
    std::cout << "  empty library region: " << nEmptyRegion << std::endl;
    std::cout << "  sampled electrons: " << nSampledElectrons << std::endl;
    std::cout << "\nTiming (event loop): " << eventLoopSeconds << " s" << std::endl;
    std::cout << "  DriftElectron: " << driftSeconds << " s"
              << " (" << 100. * driftSeconds / eventLoopSeconds << "%)" << std::endl;
    std::cout << "  template read + accumulation: " << librarySeconds << " s"
              << " (" << 100. * librarySeconds / eventLoopSeconds << "%)" << std::endl;
    std::cout << "  waveform convolution: " << convolutionSeconds << " s"
              << " (" << 100. * convolutionSeconds / eventLoopSeconds << "%)" << std::endl;
    std::cout << "Done! Output written to " << outputFile << std::endl;
    return 0;
}
