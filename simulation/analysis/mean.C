#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <ios>
#include <iostream>
#include <TGraph.h>
#include <vector>
void mean(int num = 0){
    TString inputFile = TString::Format("../build/tiersimulation/tier%d.root", num);
    TString outputFile = TString::Format("mean%d.root", num);
    TFile *fin = TFile::Open(inputFile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Failed to open file" << std::endl;
        return;
    }
    TTree *tree = (TTree *)fin->Get("tree");
    if (!tree) {
        std::cerr << "Failed to get tree" << std::endl;
        return;
    }
    std::vector<double> *t = nullptr;
    std::vector<double> *sigx2_conv = nullptr;
    std::vector<double> *sigy2_conv = nullptr;
    std::vector<double> *sigx2 = nullptr;
    std::vector<double> *sigy2 = nullptr;
    tree->SetBranchAddress("t", &t);           // 传递指针的地址
    tree->SetBranchAddress("sigx2_conv", &sigx2_conv);   // 传递指针的地址
    tree->SetBranchAddress("sigy2_conv", &sigy2_conv);   // 传递指针的地址
    tree->SetBranchAddress("sigx2", &sigx2);   // 传递指针的地址
    tree->SetBranchAddress("sigy2", &sigy2);   // 传递指针的地址
    int nEntries = tree->GetEntries();
    if (nEntries == 0) {
        std::cerr << "Tree has no entries" << std::endl;
        fin->Close();
        return;
    }

    // 先读取第一条数据以获取 vector 的大小
    tree->GetEntry(0);
    size_t nPoints = t->size();
    std::vector<double> mean_sigx;
    std::vector<double> mean_sigy;
    std::vector<double> mean_sigx_conv;
    std::vector<double> mean_sigy_conv;
    TFile *fout = new TFile(outputFile, "RECREATE");
    TTree *tout = new TTree("tree_mean", "Mean signal");
    tout->Branch("t", &t);
    tout->Branch("sigzero_ori", &mean_sigx);
    tout->Branch("sigy_ori", &mean_sigy);
    tout->Branch("sigzero", &mean_sigx_conv);
    tout->Branch("sigy", &mean_sigy_conv);

    mean_sigx.resize(nPoints, 0);
    mean_sigy.resize(nPoints, 0);
    mean_sigx_conv.resize(nPoints, 0);
    mean_sigy_conv.resize(nPoints, 0);
    for (int i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        for (size_t j = 0; j < nPoints; ++j) {
            mean_sigx[j] += sigx2->at(j);
            mean_sigy[j] += sigy2->at(j);
            mean_sigx_conv[j] += sigx2_conv->at(j);
            mean_sigy_conv[j] += sigy2_conv->at(j);
            
        }
    }
    for (size_t j = 0; j < nPoints; ++j) {
        mean_sigx[j] /= nEntries;
        mean_sigy[j] /= nEntries;
        mean_sigx_conv[j] /= nEntries;
        mean_sigy_conv[j] /= nEntries;
        // std::cout << "t: " << t->at(j) << ", mean_sigx4: " << mean_sigx4[j] << ", mean_sigy4: " << mean_sigy4[j] << std::endl;
    }

    TGraph *grx = new TGraph(nPoints);
    TGraph *gry = new TGraph(nPoints);
    TGraph *grx_conv = new TGraph(nPoints);
    TGraph *gry_conv = new TGraph(nPoints);
    grx->SetTitle("Mean sigx;Time (ns);Mean Signal");
    gry->SetTitle("Mean sigy;Time (ns);Mean Signal");
    grx_conv->SetTitle("Mean sigx_conv;Time (ns);Mean Signal");
    gry_conv->SetTitle("Mean sigy_conv;Time (ns);Mean Signal");
    for (size_t j = 0; j < nPoints; ++j) {
        grx->SetPoint(j, t->at(j), mean_sigx[j]);
        gry->SetPoint(j, t->at(j), mean_sigy[j]);
        grx_conv->SetPoint(j, t->at(j), mean_sigx_conv[j]);
        gry_conv->SetPoint(j, t->at(j), mean_sigy_conv[j]);
    }
    fout->WriteTObject(grx, "mean_sigx");
    fout->WriteTObject(gry, "mean_sigy");
    fout->WriteTObject(grx_conv, "mean_sigx_conv");
    fout->WriteTObject(gry_conv, "mean_sigy_conv");

    tout->Fill();
    fin->Close();
    fout->Write();
    fout->Close();
}