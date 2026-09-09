#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>

void merge(int maxFiles = 4000) {
    const TString inputDir = "../build/run/result";

    // avalanche 文件必须从 avalanche0.root 开始连续编号。
    int nFiles = 0;
    while (nFiles < maxFiles && !gSystem->AccessPathName(
        TString::Format("%s/avalanche%d.root", inputDir.Data(), nFiles))) {
        ++nFiles;
    }

    if (nFiles == 0) {
        std::cerr << "没有找到 avalanche*.root 文件。" << std::endl;
        return;
    }

    TFile output("merged1.root", "RECREATE");

    for (const char* treeName : {"tree", "tree_after_conv", "tree2", "tree3"}) {
        TChain chain(treeName);

        for (int i = 0; i < nFiles; ++i) {
            chain.Add(TString::Format(
                "%s/avalanche%d.root", inputDir.Data(), i));
        }

        std::cout << "合并 " << treeName << "："
                  << chain.GetEntries() << " entries" << std::endl;
        output.cd();
        chain.Merge(&output, 0, "keep");
    }

    std::cout << "已生成 merged1.root，共使用 "
              << nFiles << " 个文件。" << std::endl;
}
