// matrix.C
// 功能：
//   1. 统计当前目录下 mean*.root 文件的个数
//   2. 将每个文件中卷积好的两个方向信号 (sigzero, sigy) 汇集到 std::vector<std::vector<double>>
//   3. 将文件个数和信号矩阵存入一个 TTree

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <string>

void matrix() {
    // ============================================================
    // Step 1: 统计当前目录下 mean*.root 文件的个数
    // ============================================================
    // 方法：从 0 开始递增尝试打开 mean%d.root，直到文件不存在
    std::vector<int> fileIndices;
    int idx = 0;
    while (true) {
        TString fname = TString::Format("mean%d.root", idx);
        if (gSystem->AccessPathName(fname)) {
            // AccessPathName 返回非0表示文件不存在
            break;
        }
        fileIndices.push_back(idx);
        idx++;
    }

    int nFiles = (int)fileIndices.size();
    std::cout << "========================================" << std::endl;
    std::cout << "  找到 " << nFiles << " 个 mean root 文件: " << std::endl;
    for (int i = 0; i < nFiles; i++) {
        std::cout << "    mean" << fileIndices[i] << ".root" << std::endl;
    }
    std::cout << "========================================" << std::endl;

    if (nFiles == 0) {
        std::cerr << "错误: 当前目录下没有找到 mean*.root 文件!" << std::endl;
        return;
    }

    // ============================================================
    // Step 2: 读取每个文件的卷积信号，汇集到 vector<vector<double>>
    // ============================================================
    // sigx_all[i] = 第 i 个文件的 sigzero (x方向卷积信号)
    // sigy_all[i] = 第 i 个文件的 sigy    (y方向卷积信号)
    std::vector<std::vector<double>> sigx_all;
    std::vector<std::vector<double>> sigy_all;

    for (int i = 0; i < nFiles; i++) {
        TString fname = TString::Format("mean%d.root", fileIndices[i]);
        TFile *fin = TFile::Open(fname, "READ");
        if (!fin || fin->IsZombie()) {
            std::cerr << "警告: 无法打开 " << fname << ", 跳过." << std::endl;
            continue;
        }

        TTree *tree = (TTree *)fin->Get("tree_mean");
        if (!tree) {
            std::cerr << "警告: " << fname << " 中没有 tree_mean, 跳过." << std::endl;
            fin->Close();
            continue;
        }

        // 绑定分支 (注意: mean.C 中卷积信号分支名为 sigzero 和 sigy)
        std::vector<double> *sigx_conv = nullptr;
        std::vector<double> *sigy_conv = nullptr;
        tree->SetBranchAddress("sigzero", &sigx_conv);
        tree->SetBranchAddress("sigy",    &sigy_conv);

        tree->GetEntry(0);  // tree_mean 只有 1 个 entry

        // 将信号向量存入二维矩阵
        sigx_all.push_back(*sigx_conv);
        sigy_all.push_back(*sigy_conv);

        std::cout << "  读取 " << fname << ": sigzero 长度 = " << sigx_conv->size()
                  << ", sigy 长度 = " << sigy_conv->size() << std::endl;

        fin->Close();
    }

    int nActualFiles = (int)sigx_all.size();
    std::cout << "  成功读取 " << nActualFiles << " 个文件的卷积信号." << std::endl;

    // ============================================================
    // Step 3: 将信息存入 TTree
    // ============================================================
    // TTree 设计:
    //   - 每个 entry 对应一个 mean 文件
    //   - 分支: nFiles      (Int_t)        文件总数 (每个 entry 值相同)
    //   - 分支: fileIndex   (Int_t)        原始文件编号
    //   - 分支: sigx_conv   (vector<double>)  x 方向卷积信号
    //   - 分支: sigy_conv   (vector<double>)  y 方向卷积信号

    TFile *fout = new TFile("matrix.root", "RECREATE");
    TTree *tout = new TTree("matrix_tree", "Matrix of convolved signals from mean files");

    // 需要在 Fill 前设置分支，且每个 entry 更新这些变量
    Int_t   nFiles_branch = nActualFiles;
    Int_t   fileIndex;
    std::vector<double> sigx_vec;
    std::vector<double> sigy_vec;

    tout->Branch("nFiles",    &nFiles_branch, "nFiles/I");
    tout->Branch("fileIndex", &fileIndex,     "fileIndex/I");
    tout->Branch("sigx_conv", &sigx_vec);
    tout->Branch("sigy_conv", &sigy_vec);

    for (int i = 0; i < nActualFiles; i++) {
        fileIndex = fileIndices[i];
        sigx_vec  = sigx_all[i];
        sigy_vec  = sigy_all[i];
        tout->Fill();
    }

    fout->Write();
    fout->Close();

    std::cout << "========================================" << std::endl;
    std::cout << "  输出文件: matrix.root" << std::endl;
    std::cout << "  TTree: matrix_tree, 共 " << nActualFiles << " entries" << std::endl;
    std::cout << "  分支: nFiles, fileIndex, sigx_conv, sigy_conv" << std::endl;
    std::cout << "========================================" << std::endl;
}