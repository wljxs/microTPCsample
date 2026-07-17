// merge.C - 合并 ../build/run/result1/ 下所有 root 文件中的 TTree
// 只合并树，忽略直方图

#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <dirent.h>

#include <regex>

// 从文件名中提取数字，用于自然排序（如 avalanche1.root, avalanche2.root, ..., avalanche10.root）
int extractNumber(const std::string &filename) {
  std::regex re("(\\d+)");
  std::smatch match;
  if (std::regex_search(filename, match, re)) {
    return std::stoi(match.str());
  }
  return 0;
}

void merge(int maxFiles = 4000) {
  // 输入目录和输出文件
  const std::string inputDir = "../build/run/result/";
  const std::string outputFile = "merged1.root";

  // 要合并的树名列表（来自 microTPC.C）
  std::vector<std::string> treeNames = {"tree", "tree_after_conv", "tree2", "tree3"};

  // 收集并排序所有 root 文件
  std::vector<std::string> rootFiles;
  DIR *dir = opendir(inputDir.c_str());
  if (dir) {
    struct dirent *entry;
    while ((entry = readdir(dir)) != nullptr) {
      std::string name = entry->d_name;
      // 确保是 .root 文件且不是 merged.root
      if (name.size() > 5 && name.substr(name.size() - 5) == ".root" 
          && name != "merged.root") {
        rootFiles.push_back(inputDir + name);
      }
    }
    closedir(dir);
  }

  // 按文件名中的数字自然排序（避免 avalanche10.root 排在 avalanche2.root 前面）
  std::sort(rootFiles.begin(), rootFiles.end(),
            [](const std::string &a, const std::string &b) {
              return extractNumber(a) < extractNumber(b);
            });

  int nAvailable = (int)rootFiles.size();
  int nUse = std::min(nAvailable, maxFiles);
  std::cout << "Found " << nAvailable << " root files, using first " 
            << nUse << " files." << std::endl;

  // 创建输出文件
  TFile *fout = new TFile(outputFile.c_str(), "RECREATE");

  for (const auto &treeName : treeNames) {
    std::cout << "Merging tree: " << treeName << std::endl;

    // 创建 TChain
    TChain chain(treeName.c_str());

    // 只添加前 nUse 个文件
    int nAdded = 0;
    for (int i = 0; i < nUse; i++) {
      chain.Add(rootFiles[i].c_str());
      nAdded++;
    }

    if (nAdded == 0) {
      std::cerr << "Warning: no files added for tree " << treeName << std::endl;
      continue;
    }

    std::cout << "  Added " << nAdded << " files, "
              << chain.GetEntries() << " total entries." << std::endl;

    // 将 chain 合并写入输出文件（只写树，不写直方图）
    fout->cd();
    chain.Merge(fout, 0, "keep");
  }

  fout->Close();
  delete fout;

  std::cout << "\nMerge complete! Output: " << outputFile << std::endl;
}