
#include "AnalysisUtils.h"
#include "TInterpreter.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TSystem.h"
#include "WaveformAnalysis.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TMultiGraph.h"

#include <numeric>
#include <vector>

/**
 * 从模拟ROOT文件的tree3中读取每个event的初始电离电子位置，
 * 将气隙厚度等分为n_layer层，统计每层的电子数量和平均位置。
 * 同时画出前100个事件的电子二维分布图（x vs z 和 y vs z）。
 *
 * @param fin       已打开的模拟ROOT文件（包含tree3）
 * @param dirOut    输出目录，用于保存2D分布图（可为nullptr，则不画图）
 * @param n_layer   气隙分层数（默认5）
 * @param thickness 气隙厚度（默认5.0 mm）
 * @return          输出的TTree，包含每event每层的电子统计信息
 *                  - event/I
 *                  - layer_electronCount : vector<int>   每层电子数
 *                  - layer_electronMeanX : vector<double> 每层电子平均x (mm)
 *                  - layer_electronMeanY : vector<double> 每层电子平均y (mm)
 *                  - layer_electronMeanZ : vector<double> 每层电子平均z (mm)
 */
TTree* AnalyzeElectronLayers(TFile *fin, TDirectory *dirOut = nullptr, int n_layer = 5, double thickness = 5.0) {
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: Invalid input file for AnalyzeElectronLayers" << std::endl;
        return nullptr;
    }

    TTree *tree3 = (TTree*)fin->Get("tree3");
    if (!tree3) {
        std::cerr << "Error: Cannot find tree3 in input file" << std::endl;
        return nullptr;
    }

    // ====== 设置 tree3 分支地址 ======
    int event;
    std::vector<double> *xe0 = nullptr;
    std::vector<double> *ye0 = nullptr;
    std::vector<double> *ze0 = nullptr;

    tree3->SetBranchAddress("event", &event);
    tree3->SetBranchAddress("xe0", &xe0);
    tree3->SetBranchAddress("ye0", &ye0);
    tree3->SetBranchAddress("ze0", &ze0);

    // ====== 准备2D分布图目录（仅前100个事件） ======
    TDirectory *dir_elec2d_xz = nullptr;
    TDirectory *dir_elec2d_yz = nullptr;
    if (dirOut) {
        dir_elec2d_xz = dirOut->mkdir("electron_2d_xz");
        dir_elec2d_yz = dirOut->mkdir("electron_2d_yz");
    }

    // ====== 创建输出 TTree ======
    TTree *tout = new TTree("tree_electron_layers", "Electron counts and mean positions per layer per event");

    std::vector<int> layer_electronCount(n_layer, 0);
    std::vector<double> layer_electronMeanX(n_layer, 0.0);
    std::vector<double> layer_electronMeanY(n_layer, 0.0);
    std::vector<double> layer_electronMeanZ(n_layer, 0.0);

    tout->Branch("event", &event, "event/I");
    tout->Branch("layer_electronCount", &layer_electronCount);
    tout->Branch("layer_electronMeanX", &layer_electronMeanX);
    tout->Branch("layer_electronMeanY", &layer_electronMeanY);
    tout->Branch("layer_electronMeanZ", &layer_electronMeanZ);

    const double dz = thickness / n_layer;
    const int nDrawEvents = 100;  // 只画前100个事件

    // ====== 遍历每个 event ======
    for (int entry = 0; entry < tree3->GetEntries(); entry++) {
        tree3->GetEntry(entry);

        if (entry % 100 == 0) {
            std::cout << "AnalyzeElectronLayers: Processing event " << event << " (entry " << entry << ")" << std::endl;
        }

        // 每层累积量：电子数, sum_x, sum_y, sum_z
        std::vector<int> count(n_layer, 0);
        std::vector<double> sumX(n_layer, 0.0);
        std::vector<double> sumY(n_layer, 0.0);
        std::vector<double> sumZ(n_layer, 0.0);

        // 按 z 坐标将电子分配到各层
        for (size_t i = 0; i < ze0->size(); i++) {
            double z = (*ze0)[i]*10;
            int layer = (int)(z / dz);
            if (layer < 0) layer = 0;
            if (layer >= n_layer) layer = n_layer - 1;

            count[layer]++;
            sumX[layer] += (*xe0)[i]*10;
            sumY[layer] += (*ye0)[i]*10;
            sumZ[layer] += z;
        }

        // 计算每层平均值
        for (int l = 0; l < n_layer; l++) {
            layer_electronCount[l] = count[l];
            if (count[l] > 0) {
                layer_electronMeanX[l] = sumX[l] / count[l];
                layer_electronMeanY[l] = sumY[l] / count[l];
                layer_electronMeanZ[l] = sumZ[l] / count[l];
            } else {
                layer_electronMeanX[l] = std::numeric_limits<double>::infinity();
                layer_electronMeanY[l] = std::numeric_limits<double>::infinity();
                layer_electronMeanZ[l] = std::numeric_limits<double>::infinity();
            }
        }

        // ====== 前100个事件画2D电子分布图 (x vs z 和 y vs z) ======
        if (entry < nDrawEvents && dir_elec2d_xz && dir_elec2d_yz && ze0->size() > 0) {
            TH2D *h_elec2d_xz = new TH2D(
                Form("h_elec2d_xz_evt%d", event),
                Form("Electron 2D Distribution (x vs z) Event %d;x [mm];z [mm]", event), 20, 0, 8,
                n_tier, 0, thickness);
            TH2D *h_elec2d_yz = new TH2D(
                Form("h_elec2d_yz_evt%d", event),
                Form("Electron 2D Distribution (y vs z) Event %d;y [mm];z [mm]", event), 25, 0, 10,
                n_tier, 0, thickness);
            for (size_t i = 0; i < ze0->size(); i++) {
                h_elec2d_xz->Fill((*xe0)[i]*10, (*ze0)[i]*10);
                h_elec2d_yz->Fill((*ye0)[i]*10, (*ze0)[i]*10);
            }
            dir_elec2d_xz->WriteTObject(h_elec2d_xz);
            dir_elec2d_yz->WriteTObject(h_elec2d_yz);
            delete h_elec2d_xz;
            delete h_elec2d_yz;
        }

        tout->Fill();
    }

    std::cout << "AnalyzeElectronLayers: Completed, " << tout->GetEntries() << " events processed." << std::endl;
    return tout;
}

int main(){//分析模拟数据
    const int n_tier = 5;
    const int n_sig_samples = 15;  // must match DeconvBaseMatrix n_sig_track_sample
    int voltage = 500;

    // ====== 模板选择开关 ======
    // true:  使用 raw/matrix.root 的 5 层独立信号模板（COMSOL仿真，每层信号形状不同）
    // false: 使用 gain_and_mean_sample.root 的单模板 + Toeplitz 移位（传统方法，假设每层信号相同）
    const bool useMatrixRoot = true;
    double scale = 4200; //这个是增益的缩放因子，主要用于将模拟信号的幅值调整到实际实验的范围内，确保反卷积和后续分析的准确性。不过目前只是反卷积矩阵有这个
    DeconvCuts m_cuts;
    m_cuts.n_tier = n_tier;
    m_cuts.thickness = thickness;
    m_cuts.z_ref_position = thickness / 2.0;
    m_cuts.clusterMinEnergy[0] = 0;
    m_cuts.clusterMinEnergy[1] = 0;
    m_cuts.apvsamples = n_sig_samples;
    m_cuts.clusterMinMaxAmp = 0;
    m_cuts.clusterChi2ReLimit = 10000;
    m_cuts.clusterMinTierCharge = 0;
    double stripchargecut = 0;
    double rnormcut = 5;

    // ====== Phase 1: 构建反卷积矩阵 ======
    Eigen::MatrixXd deconvMatrix[2];

    // if (useMatrixRoot) {
        // ----- 方案 A: 使用 matrix.root 的 5 层独立信号模板 -----
        // 每层信号由 COMSOL 仿真给出，形状可能不同，比单模板+移位更精确
        std::cout << "[Phase 1] Using matrix.root (5-layer independent templates)" << std::endl;
        TFile *fmatrix = TFile::Open("../../simulation/analysis/matrix.root", "READ");
        if (!fmatrix || fmatrix->IsZombie()) {
            std::cerr << "Error: Cannot open matrix.root" << std::endl;
            return 0;
        }
        TTree *tmatrix = (TTree*)fmatrix->Get("matrix_tree");
        if (!tmatrix) {
            std::cerr << "Error: Cannot find matrix_tree in matrix.root" << std::endl;
            return 0;
        }
        std::vector<double> *sigx_conv = nullptr;
        std::vector<double> *sigy_conv = nullptr;
        tmatrix->SetBranchAddress("sigx_conv", &sigx_conv);
        tmatrix->SetBranchAddress("sigy_conv", &sigy_conv);

        std::vector<std::vector<double>> sigzero_samples[2];  // [0]=X方向5层, [1]=Y方向5层
        if (tmatrix->GetEntries() < n_tier) {
            std::cerr << "Error: matrix_tree has fewer than " << n_tier << " templates" << std::endl;
            return 0;
        }
        for (int entry = 0; entry < n_tier; entry++) {
            tmatrix->GetEntry(entry);
            if (!sigx_conv || !sigy_conv ||
                sigx_conv->size() < static_cast<size_t>(n_sig_samples * 5) ||
                sigy_conv->size() < static_cast<size_t>(n_sig_samples * 5)) {
                std::cerr << "Error: invalid template size at matrix entry " << entry << std::endl;
                return 0;
            }
            std::vector<double> sig_temp[2];
            for(int i=0;i<sigx_conv->size();i++){
              if(i%5==0 && i < n_sig_samples*5){  // 每5个采样点对应一个层，截断到 n_sig_samples
                sig_temp[0].push_back((*sigx_conv)[i]*scale);
                sig_temp[1].push_back(-(*sigy_conv)[i]*scale);
              }
            }
            sigzero_samples[0].push_back(sig_temp[0]);  // X方向第 entry 层的信号模板
            sigzero_samples[1].push_back(sig_temp[1]);  // Y方向第 entry 层的信号模板
        }
        fmatrix->Close();
        deconvMatrix[0] = AnalysisUtils::DeconvBaseMatrix(sigzero_samples[0], n_sig_samples);
        deconvMatrix[1] = AnalysisUtils::DeconvBaseMatrix(sigzero_samples[1], n_sig_samples);
    // } else {
    //     // ----- 方案 B: 传统单模板 + Toeplitz 移位矩阵 -----
    //     std::cout << "[Phase 1] Using gain_and_mean_sample.root (single template + shift)" << std::endl;
    //     TFile *fmean_sample;
    //     if (voltage == 250) {
    //         fmean_sample = TFile::Open("../raw/gain_and_mean_sample250.root", "READ");
    //     } else {
    //         fmean_sample = TFile::Open("../raw/gain_and_mean_sample.root", "READ");
    //     }

    //     if (!fmean_sample || fmean_sample->IsZombie()) {
    //         std::cerr << "Error: Cannot open mean sample file" << std::endl;
    //         return 0;
    //     }

    //     TTree *tree_mean_sample = (TTree *)fmean_sample->Get("tree_mean_sample");
    //     const char *treeName[2] = {"tree_x", "tree_y"};
    //     TTree *tree_sample[2] = {nullptr, nullptr};
    //     std::vector<double> *sigzero[2] = {nullptr, nullptr};
    //     std::vector<double> *sigtrackmean[2] = {nullptr, nullptr};
    //     std::vector<double> *sigtrackmean_conv[2] = {nullptr, nullptr};
    //     TGraph *gr_trackmean_tran[2] = {nullptr, nullptr};
    //     gr_trackmean_tran[0] = (TGraph *)fmean_sample->Get("gr_trackmean_tran_d0");
    //     gr_trackmean_tran[1] = (TGraph *)fmean_sample->Get("gr_trackmean_tran_d1");
    //     TGraph *gr_trackmean_conv_tran[2] = {nullptr, nullptr};
    //     gr_trackmean_conv_tran[0] = (TGraph *)fmean_sample->Get("gr_trackmean_conv_tran_d0");
    //     gr_trackmean_conv_tran[1] = (TGraph *)fmean_sample->Get("gr_trackmean_conv_tran_d1");
    //     if (!gr_trackmean_tran[0] || !gr_trackmean_tran[1] || !gr_trackmean_conv_tran[0] || !gr_trackmean_conv_tran[1]) {
    //         std::cerr << "Error: Cannot find gr_trackmean_tran or gr_trackmean_conv_tran in mean sample file" << std::endl;
    //         return 0;
    //     }
    //     for (int d = 0; d < 2; ++d) {
    //         tree_sample[d] = (TTree *)fmean_sample->Get(treeName[d]);
    //         tree_sample[d]->SetBranchAddress("sigzero", &sigzero[d]);
    //         tree_sample[d]->SetBranchAddress("sigtrackmean", &sigtrackmean[d]);
    //         tree_sample[d]->SetBranchAddress("sigtrackmean_conv", &sigtrackmean_conv[d]);
    //         tree_sample[d]->GetEntry(0);
    //     }

    //     std::vector<double> sigzero_sample[2];
    //     sigzero_sample[0] = std::vector<double>(sigzero[0]->begin(), sigzero[0]->end());
    //     sigzero_sample[1] = std::vector<double>(sigzero[1]->begin(), sigzero[1]->end());
    //     deconvMatrix[0] = AnalysisUtils::DeconvBaseMatrix(sigzero_sample[0], n_sig_samples, n_tier);
    //     deconvMatrix[1] = AnalysisUtils::DeconvBaseMatrix(sigzero_sample[1], n_sig_samples, n_tier);

    //     fmean_sample->Close();
    // }

    if (deconvMatrix[0].size() == 0 || deconvMatrix[1].size() == 0) {
        std::cerr << "Error: Failed to initialize deconvolution base matrix!" << std::endl;
        return 0;
    }


    // ====== Phase 2: 分析模拟的数据 ======
    TFile *fin = TFile::Open("../../simulation/result/nodelta/20test.root", "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: Cannot open test1.root" << std::endl;
        return 0;
    }
    TTree *tree = (TTree *)fin->Get("tree");
    if (!tree) {
        std::cerr << "Error: Cannot find tree in 0test1.root" << std::endl;
        return 0;
    }
    int event;
    std::vector<double> *t = nullptr;
    std::vector<int> *types = nullptr;
    std::vector<std::vector<double>> *waveforms = nullptr;
    std::vector<int> *stripIDs = nullptr;
    double k[2], b[2];
    tree->SetBranchAddress("event", &event);
    tree->SetBranchAddress("t", &t);
    tree->SetBranchAddress("types", &types);
    tree->SetBranchAddress("waveforms_conv", &waveforms);
    tree->SetBranchAddress("stripIDs", &stripIDs);
    tree->SetBranchAddress("kx", &k[0]);
    tree->SetBranchAddress("ky", &k[1]);
    tree->SetBranchAddress("bx", &b[0]);
    tree->SetBranchAddress("by", &b[1]);

    std::vector<std::vector<double>> vec_ApvValue_temp[2];        // 保存每个条的每个apv点的adc值，用于后续画分布图
    std::vector<std::vector<double>> vec_tracksignal_temp[2];      // 保存每个条的每个apv点的拟合信号值，用于后续画分布图
    std::vector<StripHitDeconv> vec_StripHitDeconv_temp[2]; // 按照类型保存，后续聚类成clusterdeconv时用
    std::vector<int> vec_StripHitDeconvid_temp[2];          // 同上，但保存的是条号和类型信息，用于后续聚类成clusterdeconv

    std::vector<StripHitDeconv> StripHitDeconvs; // 保存每个事件的条信息，按类型保存
    std::vector<ClusterDeconv> ClusterDeconvs; // 保存每个事件的cluster信息，按类型保存
    int clusternum[2] = {0, 0}; // 保存每个事件的cluster数量，按类型保存
    TFile *fout = TFile::Open("../result/20test_nodelta_midestbig.root", "RECREATE");
    TTree *tout = new TTree("tout", "tout");
    tout->Branch("event", &event);
    tout->Branch("ClusterDeconvs", &ClusterDeconvs);
    tout->Branch("StripHitDeconvs", &StripHitDeconvs);
    tout->Branch("clusternumx", &clusternum[0]);
    tout->Branch("clusternumy", &clusternum[1]);

    // ====== Phase 3: 用tree3统计每层电子的数量和平均位置 ======
    std::cout << "===== Phase 3: Analyzing electron layers from tree3 =====" << std::endl;
    TTree *tree_electron = AnalyzeElectronLayers(fin, fout, n_tier, 5.0);
    if (tree_electron) {
        tree_electron->SetDirectory(fout);  // 关联到输出文件，后续 Close 时自动写入
    }

    TDirectory *dir_becut = fout->mkdir("waveforms_noovertherhold");
    TDirectory *dir = fout->mkdir("waveforms_overtherhold");
    TH1D *h_res[2][n_tier];
    TH1D *h_res_microTPC[2];
    TH1D *h_res_angle[2];
    TH2D *h_correlation[2];
    TH2D *h_cov[2];
    for(int type = 0; type < 2; ++type) {
        h_res_microTPC[type] = new TH1D(Form("h_res_type_%d", type), Form("h_res_type_%d", type), 200, -1, 1);
        h_res_angle[type] = new TH1D(Form("h_res_angle_type_%d", type), Form("h_res_angle_type_%d", type), 200, -1, 1);
        h_correlation[type] = new TH2D(Form("h_correlation_type_%d", type), Form("h_correlation_type_%d", type), n_tier, 0, n_tier, n_tier, 0, n_tier);
        h_cov[type] = new TH2D(Form("h_cov_type_%d", type), Form("h_cov_type_%d", type), n_tier, 0, n_tier, n_tier, 0, n_tier);
        for(int tier = 0; tier < n_tier; ++tier) {
            h_res[type][tier] = new TH1D(Form("h_res_type_%d_tier_%d", type, tier), Form("h_res_type_%d_tier_%d", type, tier), 200, -1, 1);
        }
    }
    TDirectory *dir_charge2d = fout->mkdir("charge2d");
    TDirectory *dir_waveform2d = fout->mkdir("waveform2d");
    TDirectory *dir_charge_one[2] = {fout->mkdir("charge_one_x"), fout->mkdir("charge_one_y")};
    TDirectory *dir_waveform_one[2] = {fout->mkdir("waveform_one_x"), fout->mkdir("waveform_one_y")};
    
    TDirectory *dir_correlation = fout->mkdir("correlation");
    TDirectory *dir_correlation_one[2] = {dir_correlation->mkdir("correlation_one_x"), dir_correlation->mkdir("correlation_one_y")};
    TDirectory *dir_residual_diagnostics = dir_correlation->mkdir("residual_diagnostics");

    // X 分层残差诊断：区分仅要求 X 有单 cluster 与同时要求 X/Y 有单 cluster，
    // 并在通常的 X/Y 单 cluster 样本中按本层电荷的高低各取一半。
    TH1D *h_res_xonly[n_tier];
    TH1D *h_res_charge_low[n_tier];
    TH1D *h_res_charge_high[n_tier];
    std::vector<std::pair<double, double>> res_charge_x[n_tier]; // (residual [mm], tier charge)
    for (int tier = 0; tier < n_tier; ++tier) {
      h_res_xonly[tier] = new TH1D(
          Form("h_res_xonly_tier_%d", tier),
          Form("X residual, X-only single-cluster selection, tier %d;residual [mm];entries", tier),
          200, -1, 1);
      h_res_charge_low[tier] = new TH1D(
          Form("h_res_x_lowcharge_tier_%d", tier),
          Form("X residual, lower half of tier charge, tier %d;residual [mm];entries", tier),
          200, -1, 1);
      h_res_charge_high[tier] = new TH1D(
          Form("h_res_x_highcharge_tier_%d", tier),
          Form("X residual, upper half of tier charge, tier %d;residual [mm];entries", tier),
          200, -1, 1);
    }


    std::vector<double> res_arr[2][n_tier];

    for(int entry = 0; entry < tree->GetEntries(); entry++){
      if(entry % 100 == 0) std::cout << "Processing entry: " << entry << std::endl;
        // 每轮事件清空临时向量，避免跨事件累积
        StripHitDeconvs.clear();
        ClusterDeconvs.clear();
        for (int t = 0; t < 2; ++t) {
            vec_ApvValue_temp[t].clear();
            vec_tracksignal_temp[t].clear();
            vec_StripHitDeconv_temp[t].clear();
            vec_StripHitDeconvid_temp[t].clear();
            clusternum[t] = 0;
        }
        tree->GetEntry(entry);
        if (!waveforms || !types || !stripIDs) {
            std::cerr << "Error: Null branch pointer after GetEntry — check branch names in test.root" << std::endl;
            continue;
        }
        if (waveforms->size() != types->size() || waveforms->size() != stripIDs->size()) {
            std::cerr << "Error: inconsistent strip branch sizes in event " << event << std::endl;
            continue;
        }
        for(int i = 0; i < types->size(); i++){
            std::vector<double> waveform1 = waveforms->at(i);
            std::vector<double> waveform;
            int stripID = stripIDs->at(i);
            int type = types->at(i);
            if(type==1){
              for(int j = 0; j < waveform1.size(); j++) {
                waveform1[j] = waveform1[j] * -1; //y方向的信号是反过来的
              }
            }
            if(waveform1.size() < n_sig_samples * 5) {
                std::cerr << "Warning: Waveform size is smaller than expected for strip ID " << stripIDs->at(i) << std::endl;
                continue;
            }
            for(int j = 0; j < waveform1.size(); j++) {
              if(j%5==0 && j < n_sig_samples*5){
                waveform.push_back(waveform1[j]*1000);
              }
            }
            if (waveform.size() != static_cast<size_t>(n_sig_samples)) {
                std::cerr << "Warning: sampled waveform has " << waveform.size()
                          << " points for strip ID " << stripID << std::endl;
                continue;
            }

            StripHitDeconv sh;
            sh.isBad = false;
            sh.ID = stripID;
            sh.type = type;            
            if (!AnalyzeWaveformFeatures(waveform, 0, 0, &sh)) continue;
            // Resize waveform to match deconvolution matrix row count (n_sig_samples)
            auto [x_vec, rnorm] = AnalysisUtils::SolveNNLSLasso(waveform, deconvMatrix[sh.type], 0.00);
            sh.chargetier = x_vec;
      
            // std::cout << x_vec[0] << "," << x_vec[1] << "," << x_vec[2] << "," << x_vec[3] << "," << x_vec[4] << std::endl;
            const double charge = std::accumulate(x_vec.begin(), x_vec.end(), 0.0);            
            if(charge == 0) {
                // std::cerr << "Warning: Charge is zero for " << "type " << type << ",strip ID "<< stripID << ",event " << event << std::endl;
                TGraph *gr = new TGraph(waveform.size());
                for (size_t k = 0; k < waveform.size(); ++k) {
                    gr->SetPoint(k, k, waveform[k]);
                }
                gr->SetName(Form("gr_waveform_stripID_%d_type_%d(event_%d)", stripID, type, event));
                gr->SetTitle(Form("Waveform for Strip ID %d (Type %d,Event %d)", stripID, type, event));
                dir_becut->WriteTObject(gr);
                delete gr;
                continue;
            }
            else{
                TCanvas *c = new TCanvas(Form("c_waveform_stripID_%d_type_%d(event_%d)", stripID, type, event), Form("Waveform for Strip ID %d (Type %d,Event %d)", stripID, type, event), 800, 600);
                std::vector<double> fit_waveform = FitWaveform(deconvMatrix[sh.type], sh);
                TGraph *gr = new TGraph(waveform.size());
                TGraph *gr_fit = new TGraph(fit_waveform.size());
                for (size_t k = 0; k < waveform.size(); ++k) {
                    gr->SetPoint(k, k, waveform[k]);
                    gr_fit->SetPoint(k, k, fit_waveform[k]);
                }
                gr->SetName(Form("gr_waveform_stripID_%d_type_%d(event_%d)", stripID, type, event));
                gr->SetTitle(Form("Waveform for Strip ID %d (Type %d,Event %d)", stripID, type, event));
                gr->SetLineColor(kBlue);
                gr_fit->SetLineColor(kRed);
                TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
                legend->AddEntry(gr, "Original Waveform", "l");
                legend->AddEntry(gr_fit, "Fitted Waveform", "l");
                gr->Draw("AL");
                gr_fit->Draw("L SAME");
                legend->Draw();
                dir->WriteTObject(c);
                delete c;
                delete legend;
                delete gr_fit;
                delete gr;              
            }
            for(int j = 0; j < n_tier; ++j) {
                if (x_vec[j] < stripchargecut) {
                     x_vec[j] = 0; // Set to zero if below threshold
                }
            }
            sh.chargetier = x_vec;
            sh.rnorm = rnorm / charge;
            if(sh.rnorm > rnormcut) {
              continue; // Skip this strip if rnorm exceeds the threshold
            }
            sh.rnorm_no = rnorm;   
            vec_StripHitDeconv_temp[sh.type].push_back(sh);
            vec_StripHitDeconvid_temp[sh.type].push_back(sh.ID);
            vec_ApvValue_temp[sh.type].push_back(waveform);
            vec_tracksignal_temp[sh.type].push_back(waveform);
            StripHitDeconvs.push_back(sh);
        }
    

      // ======phase 6: 聚类，把相邻的条聚成cluster ======
      int hasCluster[2] = {0, 0};      
      std::vector<ClusterDeconv> vec_ClusterDeconv_temp; // cluster的临时变量
      for (int type = 0; type < 2; type++) {
        if (vec_StripHitDeconvid_temp[type].size() == 0)
          continue;
        for (int m = 0; m < vec_StripHitDeconvid_temp[type].size() - 1; m++) {
          for (int n = m + 1; n < vec_StripHitDeconvid_temp[type].size(); n++) {
            if (vec_StripHitDeconvid_temp[type][m] > vec_StripHitDeconvid_temp[type][n]) {
              std::swap(vec_StripHitDeconvid_temp[type][m], vec_StripHitDeconvid_temp[type][n]);
              std::swap(vec_StripHitDeconv_temp[type][m], vec_StripHitDeconv_temp[type][n]);
              std::swap(vec_tracksignal_temp[type][m], vec_tracksignal_temp[type][n]);
              std::swap(vec_ApvValue_temp[type][m], vec_ApvValue_temp[type][n]);
            }
          }
        }
      }
      // 聚类（调用独立函数）
      for (int type = 0; type < 2; type++) {
        ClusterStripHits(vec_StripHitDeconv_temp[type], nullptr, m_cuts, event, type,
                         hasCluster[type], vec_ClusterDeconv_temp);
      }
      ClusterDeconvs.clear();
      std::vector<ClusterDeconv> ClusterDeconvs_one[2];
      for(const auto& cluster : vec_ClusterDeconv_temp) {
        if (!cluster.isValid) continue;
        ClusterDeconvs.push_back(cluster);
        clusternum[cluster.type]++;
        ClusterDeconvs_one[cluster.type].push_back(cluster);
      }

      // 对照 1：仅要求 X 方向有一个有效 cluster，不对 Y 重建施加选择。
      if (clusternum[0] == 1) {
        const ClusterDeconv &clusterX = ClusterDeconvs_one[0][0];
        for (int tier = 0; tier < n_tier; ++tier) {
          if (!std::isfinite(clusterX.pos[tier])) continue;
          const double res = clusterX.pos[tier] * 0.4 + 0.2
              - b[0] * 10.0
              - k[0] * thickness / n_tier * (tier + 0.5);
          h_res_xonly[tier]->Fill(res);
        }
      }

      //analysis and draw the events with only one cluster in each direction

      if(clusternum[0] == 1&& clusternum[1] == 1) {
        for(int type = 0;type < 2;type++){
          TH2D *h_waveform_one = new TH2D(Form("h_waveform_one_type_%d", type), Form("h_waveform_one_type_%d", type), 20, 0, 20, n_sig_samples, 0, n_sig_samples);
          TH2D *h_charge_one = new TH2D(Form("h_charge_one_type_%d", type), Form("h_charge_one_type_%d", type), 20, 0, 20, n_tier, 0, thickness);
          TGraph *gr_hit = new TGraph();
          gr_hit->SetName(Form("gr_track_type_%d", type));
          gr_hit->SetTitle(Form("gr_track_type_%d", type));
          TGraph *gr_pred = new TGraph();
          gr_pred->SetName(Form("gr_track_pred_type_%d", type));
          gr_pred->SetTitle(Form("gr_track_pred_type_%d", type));
          TGraphErrors *gr_position = new TGraphErrors();
          gr_position->SetName(Form("gr_position_type_%d", type));
          gr_position->SetTitle(Form("gr_position_type_%d", type));

          for(int i = 0;i<vec_StripHitDeconv_temp[type].size();i++){
            StripHitDeconv sh = vec_StripHitDeconv_temp[type][i];
            std::vector<double> waveform = vec_ApvValue_temp[type][i];
            for(int j = 0;j<waveform.size();j++){
              h_waveform_one->SetBinContent(sh.ID+1,j+1,waveform[j]);
            }
            for(int j = 0;j<n_tier;j++){
              h_charge_one->SetBinContent(sh.ID+1,j+1,sh.chargetier[j]);
            }
          }

          for(int i = 0;i<ClusterDeconvs_one[type].size();i++){
            ClusterDeconv cluster = ClusterDeconvs_one[type][i];

            // 计算残差并填充直方图
            h_res_microTPC[cluster.type]->Fill(cluster.microTPCposition*0.4+0.2-b[cluster.type]*10-k[cluster.type]*thickness/2);//0.2是第0个strip的bias
            h_res_angle[cluster.type]->Fill(cluster.k*0.4-k[cluster.type]);

            for(int tier = 0; tier < n_tier; ++tier) {
              if(cluster.pos[tier] == std::numeric_limits<double>::infinity()) {
                res_arr[cluster.type][tier].push_back(std::numeric_limits<double>::infinity());
                continue;
              }
              double res = cluster.pos[tier]*0.4+0.2-b[cluster.type]*10-k[cluster.type]*thickness/n_tier*(tier+0.5);//0.2是第0个strip的bias
              gr_position->SetPoint(gr_position->GetN(), cluster.pos[tier]+0.5, thickness/n_tier*(tier+0.5));
              gr_position->SetPointError(gr_position->GetN(), sqrt(thickness/n_tier*(tier+0.5)/cluster.chargetier[tier]*1000), 0);
              res_arr[cluster.type][tier].push_back(res);
              h_res[cluster.type][tier]->Fill(res);
              if (cluster.type == 0) {
                res_charge_x[tier].emplace_back(res, cluster.chargetier[tier]);
              }
            }

            double predk = cluster.k;
            double predb = cluster.b+0.2*2.5;//0.2是第0个strip的bias
            for(int j = 0;j<n_tier;j++){
              double z = thickness/n_tier*(j+0.5);
              double pred = predk*z+predb;
              double hit = k[type]*z*2.5+b[type]*25;
              gr_pred->SetPoint(j,pred,z);
              gr_hit->SetPoint(j, hit, z);
            }
          }
          TCanvas *c_compare = new TCanvas(Form("c_compare_evt%d_type%d", event, type), Form("Event %d %s Comparison", event, type == 0 ? "X" : "Y"), 800, 600);
          c_compare->cd();
          h_charge_one->Draw("COLZ");
          gr_pred->SetLineColor(kRed);
          gr_pred->SetLineWidth(2);
          gr_pred->Draw("L SAME");
          gr_hit->SetLineColor(kBlue);
          gr_hit->SetLineWidth(2);
          gr_hit->Draw("L SAME");
          gr_position->SetMarkerStyle(20);
          gr_position->SetMarkerColor(kRed);
          gr_position->Draw("P SAME");

          // TLegend: 各 tier 残差 + 角度差
          {
            ClusterDeconv &c = ClusterDeconvs_one[type][0];  // 只有一个 cluster
            TLegend *leg = new TLegend(0.12, 0.65, 0.48, 0.90);
            leg->SetName(Form("leg_evt%d_type%d", event, type));
            leg->SetHeader(Form("Event %d  #color[4]{#bullet} pred  #color[2]{line} true", event), "C");
            for(int tier = 0; tier < n_tier; ++tier) {
              if (c.pos[tier] != std::numeric_limits<double>::infinity()) {
                double res = c.pos[tier]*0.4+0.2-b[c.type]*10-k[c.type]*thickness/n_tier*(tier+0.5);
                leg->AddEntry((TObject*)0, Form("Tier%d res = %.3f mm", tier, res), "");
              } else {
                leg->AddEntry((TObject*)0, Form("Tier%d: no data", tier), "");
              }
            }
            double ang_diff = c.k*0.4-k[c.type];
            leg->AddEntry((TObject*)0, Form("#Deltak = %.4f (pred %.4f vs true %.4f)", ang_diff, c.k*0.4, k[c.type]), "");
            double mtpc_res = c.microTPCposition*0.4+0.2-b[c.type]*10-k[c.type]*thickness/2;
            leg->AddEntry((TObject*)0, Form("microTPC res = %.3f mm", mtpc_res), "");
            leg->Draw();
          }
          dir_charge_one[type]->WriteTObject(c_compare);
          dir_waveform_one[type]->WriteTObject(h_waveform_one);
          delete h_charge_one;
          delete h_waveform_one;
          delete gr_pred;
          delete gr_hit;
          delete c_compare;
          delete gr_position;
        }
      }
      tout->Fill();
    }
    std::cout << "Simulation analysis completed." << std::endl;

    // 每个 tier 对只使用两者均有效的同一批 event 计算协方差。
    double sum_i[2][n_tier][n_tier] = {{{0}}};
    double sum_j[2][n_tier][n_tier] = {{{0}}};
    double sum_i2[2][n_tier][n_tier] = {{{0}}};
    double sum_j2[2][n_tier][n_tier] = {{{0}}};
    double sum_ij[2][n_tier][n_tier] = {{{0}}};
    int count[2][n_tier][n_tier] = {{{0}}};

    for(int type = 0; type < 2; ++type) {
      for(int i = 0; i < res_arr[type][0].size(); ++i) {
        for(int tier = 0; tier < n_tier; ++tier) {
          for(int tier2 = tier; tier2 < n_tier; ++tier2){
            const double res1 = res_arr[type][tier][i];
            const double res2 = res_arr[type][tier2][i];
            if(!std::isfinite(res1) || !std::isfinite(res2)) continue;
            sum_i[type][tier][tier2] += res1;
            sum_j[type][tier][tier2] += res2;
            sum_i2[type][tier][tier2] += res1 * res1;
            sum_j2[type][tier][tier2] += res2 * res2;
            sum_ij[type][tier][tier2] += res1 * res2;
            count[type][tier][tier2]++;
          }
        }
      }

      for(int tier = 0; tier < n_tier; ++tier) {     
        for(int tier2 = tier; tier2 < n_tier; ++tier2) {
          if(count[type][tier][tier2] > 0) {
            const double n = count[type][tier][tier2];
            const double mean_i = sum_i[type][tier][tier2] / n;
            const double mean_j = sum_j[type][tier][tier2] / n;
            const double cov = sum_ij[type][tier][tier2] / n - mean_i * mean_j;
            const double var_i = sum_i2[type][tier][tier2] / n - mean_i * mean_i;
            const double var_j = sum_j2[type][tier][tier2] / n - mean_j * mean_j;
            const double corr = (var_i > 0.0 && var_j > 0.0)
                ? cov / std::sqrt(var_i * var_j) : 0.0;
            h_cov[type]->SetBinContent(tier+1, tier2+1, cov);
            h_cov[type]->SetBinContent(tier2+1, tier+1, cov);  // 对称矩阵
            h_correlation[type]->SetBinContent(tier+1, tier2+1, corr);
            h_correlation[type]->SetBinContent(tier2+1, tier+1, corr);  // 对称矩阵
          }
        }
      }
    }

    // 对照 2：按每层电荷排序，低/高各保留一半，避免硬编码电荷阈值。
    for (int tier = 0; tier < n_tier; ++tier) {
      auto ordered = res_charge_x[tier];
      std::sort(ordered.begin(), ordered.end(),
                [](const auto &lhs, const auto &rhs) { return lhs.second < rhs.second; });
      const size_t split = ordered.size() / 2;
      for (size_t i = 0; i < ordered.size(); ++i) {
        if (i < split) h_res_charge_low[tier]->Fill(ordered[i].first);
        else h_res_charge_high[tier]->Fill(ordered[i].first);
      }
      if (!ordered.empty()) {
        const double boundary = ordered[split].second;
        std::cout << "X tier " << tier << " charge split: " << boundary
                  << " (low/high entries " << split << "/"
                  << ordered.size() - split << ")" << std::endl;
      }
      dir_residual_diagnostics->WriteTObject(h_res_xonly[tier]);
      dir_residual_diagnostics->WriteTObject(h_res_charge_low[tier]);
      dir_residual_diagnostics->WriteTObject(h_res_charge_high[tier]);
      delete h_res_xonly[tier];
      delete h_res_charge_low[tier];
      delete h_res_charge_high[tier];
    }


    for(int type = 0; type < 2; ++type) {
        dir_correlation_one[type]->WriteTObject(h_correlation[type]);
        dir_correlation_one[type]->WriteTObject(h_cov[type]);
        for(int tier = 0; tier < n_tier; ++tier) {
            dir_correlation_one[type]->WriteTObject(h_res[type][tier]);
            delete h_res[type][tier];
        }
        delete h_correlation[type];
        delete h_cov[type];
    }
    tout->Write();
    fout->Write();
    fout->Close();
    delete fout;
    delete fin;
    return 0;
}
