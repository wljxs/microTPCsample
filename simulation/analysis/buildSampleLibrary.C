// ============================================================
// sampleconstruction.C — 抽样库构建
//
//   buildSampleLibrary() : 从 avalancheunit.root 构建 samplelibrary.root
// ============================================================

#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <TGraph.h>

// ---------- 全局参数 ----------
const int    kNbins        = 10;       // 每维 10 个 bin，共 100 区域
const double kEz1Target    = 0.01;     // 筛选 ez1 = 0.01
const double kEz1Tolerance = 0.000001;    // 浮点容差

// ---- 信号通道名称表（8 个，与 microTPC.C 一致）----
static const int n_sig = 8;
static const char* sig_name[n_sig] = {
        "sigxtestminus1","sigxtest1", 
    "sigxtestminus2",  "sigzero","sigxtest2",   
     "sigytestminus1", "sigy","sigytest1" 
};

// ============================================================
// buildSampleLibrary — 构建抽样库
//   从 avalancheunit.root 读取 tree3/tree，
//   按 ez1=0.01 筛选，用 ex1×ey2 分 100 区域，输出到 samplelibrarysmall.root
// ============================================================
void buildSampleLibrary(
    const char *inputFile  = "merged1.root",
    const char *outputFile = "samplelibrarysmall.root")
{
    // ----- 打开输入文件 -----
    TFile *fin = TFile::Open(inputFile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: Cannot open " << inputFile << std::endl;
        return;
    }

    TTree *tree3 = (TTree *)fin->Get("tree3");
    TTree *tree1 = (TTree *)fin->Get("tree");               // 原始信号

    if (!tree3 || !tree1) {
        std::cerr << "Error: Missing trees in input file." << std::endl;
        fin->Close();
        return;
    }

    // ----- 绑定 tree3 分支（用于分类）-----
    // 与 microTPC.C 中 tree3 的分支一致：ex1, ey1, ez1, et1, ee1, it2
    // 另读 ex2, ey2, ez2 等（如文件中存在的话）
    double x0, y0, z0, t0, e0;

    // 先尝试标准 microTPC.C 的分支
    tree3->SetBranchAddress("x0", &x0);
    tree3->SetBranchAddress("y0", &y0);
    tree3->SetBranchAddress("z0", &z0);
    tree3->SetBranchAddress("t0", &t0);
    tree3->SetBranchAddress("e0", &e0);


    // ----- 绑定 tree (tree1) 的分支 — 与 microTPC.C 完全一致 -----
    int event;
    std::vector<double> *t = nullptr;
    std::vector<double> *sig_in[n_sig] = {};

    tree1->SetBranchAddress("event", &event);
    tree1->SetBranchAddress("t",     &t);
    for (int is = 0; is < n_sig; ++is)
        tree1->SetBranchAddress(sig_name[is], &sig_in[is]);

    // ----- 第一遍：扫描 tree2 收集满足 ez1=0.01 的条目 -----
    Long64_t nTotal = tree3->GetEntries();
    std::vector<double>   vec_x0, vec_y0;
    std::vector<Long64_t> vec_entry;  // 同时记录 entry 索引，便于后续读取 tree1/tree2 的信号数据

    std::cout << "[buildSampleLibrary] Scanning " << nTotal
              << " entries in tree3 for z0=" << kEz1Target << " ..." << std::endl;
    int entry = 0;


    double x_min = 0.056;
    double x_max = 0.064;
    double y_min = 0.04;
    double y_max = 0.08;
    double dx = (x_max - x_min) / kNbins;
    double dy = (y_max - y_min) / kNbins;


    // ----- 按区域分类 entry 索引 -----
    std::vector<Long64_t> region_entries[kNbins][kNbins]; 
    for (Long64_t i = 0; i < nTotal; ++i) {
        tree3->GetEntry(i);
        vec_x0.push_back(x0);
        vec_y0.push_back(y0);
        vec_entry.push_back(entry);

        int ix = (int)((x0 - x_min) / dx);
        int iy = (int)((y0 - y_min) / dy);
        if (ix < 0 || ix >= kNbins || iy < 0 || iy >= kNbins) {
            std::cout << "Warning: Entry " << entry << " with (x0, y0) = (" << x0 << ", " << y0 << ") is out of bounds and will be ignored." << std::endl;
            continue;
        }
        region_entries[ix][iy].push_back(entry);
        entry++;
    }

    Long64_t nValid = vec_entry.size();
    std::cout << "Found " << nValid << " valid entries." << std::endl;
    if (nValid == 0) { fin->Close(); return; }

    std::cout << "ex1 range: [" << x_min << ", " << x_max << "]" << std::endl;
    std::cout << "ey2 range: [" << y_min << ", " << y_max << "]" << std::endl;
    std::cout << "Grid: " << kNbins << "x" << kNbins << "  bin width: ex1=" << dx << ", ey2=" << dy << std::endl;

    // ----- 创建输出文件 + 抽样库 tree -----
    TFile *fout = new TFile(outputFile, "RECREATE");

    // 保存分箱元数据 + 共享时间轴
    TTree *meta = new TTree("meta", "Binning & time metadata");
    double meta_xmin, meta_xmax, meta_ymin, meta_ymax;
    int    meta_nbins;
    std::vector<double> *meta_t = nullptr;
    meta->Branch("xmin",  &meta_xmin,  "xmin/D");
    meta->Branch("xmax", &meta_xmax, "xmax/D");
    meta->Branch("ymin",  &meta_ymin,  "ymin/D");
    meta->Branch("ymax", &meta_ymax, "ymax/D");
    meta->Branch("nbins",    &meta_nbins,    "nbins/I");
    meta->Branch("t",        &meta_t);

    meta_xmin  = x_min;
    meta_xmax = x_max;
    meta_ymin  = y_min;
    meta_ymax = y_max;
    meta_nbins    = kNbins;

    // 读取第一条有效 entry，获取共享时间轴
    tree3->GetEntry(vec_entry[0]);
    tree1->GetEntry(vec_entry[0]);
    meta_t = t;
    meta->Fill();
    meta->Write();

    // ----- 逐区域创建独立树，写入该区域的信号 -----
    Long64_t totalInLib = 0;
    int nRegionsWithData = 0;


    TTree *regiontree[10][10];
    for (int i = 0; i < kNbins; ++i){
        for (int j = 0; j < kNbins; ++j){
            regiontree[i][j] = new TTree(Form("region_%d_%d", i, j), Form("Signals in region (%d, %d)", i, j));
            regiontree[i][j]->Branch("x0", &x0);
            regiontree[i][j]->Branch("y0", &y0);
            for (int is = 0; is < n_sig; ++is) {
                regiontree[i][j]->Branch(sig_name[is], &sig_in[is]);}
        }
    }

    for (int ix = 0; ix < kNbins; ++ix) {
        for (int iy = 0; iy < kNbins; ++iy) {
            const auto &entries = region_entries[ix][iy];
            if (entries.empty()) continue;

            std::cout << "Region (" << ix << ", " << iy << ") has " << entries.size() << " entries." << std::endl;
            for (int j=0; j < entries.size(); ++j) {
                Long64_t entryIdx = entries[j];
                std::cout << "  Entry " << entryIdx << ": x0=" << vec_x0[entryIdx]
                          << ", y0=" << vec_y0[entryIdx] << std::endl;
                tree3->GetEntry(entryIdx);
                tree1->GetEntry(entryIdx);
                regiontree[ix][iy]->Fill();
                ++totalInLib;
            }
            regiontree[ix][iy]->Write();
            ++nRegionsWithData;
        }
    }
    fout->Close();
    fin->Close();

    std::cout << "\n[buildSampleLibrary] Done!  " << totalInLib
              << " entries in " << nRegionsWithData << " region trees -> " << outputFile << std::endl;
    std::cout << "  Meta tree:    \"meta\" (boundaries + shared time axis)" << std::endl;
    std::cout << "  Region trees: \"region_0\" .. \"region_" << (kNbins * kNbins - 1)
              << "\" (" << n_sig << " signal channels each, nSamples = tree->GetEntries())" << std::endl;
}