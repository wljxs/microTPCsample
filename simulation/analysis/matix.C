#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TParameter.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

void matrix(TString inputDir = ".", TString outputFile = "matrix.root") {
    // mean 文件必须从 mean0.root 开始连续编号。
    int nFiles = 0;
    while (!gSystem->AccessPathName(
        TString::Format("%s/mean%d.root", inputDir.Data(), nFiles))) {
        ++nFiles;
    }

    if (nFiles == 0) {
        std::cerr << "没有找到 mean*.root 文件。" << std::endl;
        return;
    }

    TFile output(outputFile, "RECREATE");
    TTree matrixTree("matrix_tree", "Matrix of convolved signals");

    int fileIndex = 0;
    double template_z_mm = 0.;
    std::vector<double> sigx_conv;
    std::vector<double> sigy_conv;

    matrixTree.Branch("nFiles", &nFiles);
    matrixTree.Branch("fileIndex", &fileIndex);
    matrixTree.Branch("template_z_mm", &template_z_mm);
    matrixTree.Branch("sigx_conv", &sigx_conv);
    matrixTree.Branch("sigy_conv", &sigy_conv);

    TMultiGraph meanX;
    TMultiGraph meanY;
    meanX.SetTitle("X direction;Time (ns);Mean signal");
    meanY.SetTitle("Y direction;Time (ns);Mean signal");

    TLegend legendX(0.62, 0.65, 0.89, 0.89);
    TLegend legendY(0.62, 0.65, 0.89, 0.89);
    legendX.SetBorderSize(0);
    legendY.SetBorderSize(0);

    std::vector<TGraph*> graphs;
    gStyle->SetPalette(kRainBow);

    for (fileIndex = 0; fileIndex < nFiles; ++fileIndex) {
        const TString filename = TString::Format(
            "%s/mean%d.root", inputDir.Data(), fileIndex);
        TFile input(filename, "READ");
        auto* meanTree = static_cast<TTree*>(input.Get("tree_mean"));

        if (input.IsZombie() || !meanTree) {
            std::cerr << "无法读取 " << filename << std::endl;
            return;
        }

        std::vector<double>* t = nullptr;
        std::vector<double>* sigzero = nullptr;
        std::vector<double>* sigy = nullptr;
        meanTree->SetBranchAddress("t", &t);
        meanTree->SetBranchAddress("sigzero", &sigzero);
        meanTree->SetBranchAddress("sigy", &sigy);
        meanTree->GetEntry(0);

        if (!t || !sigzero || !sigy) {
            std::cerr << filename << " 中的平均波形为空。" << std::endl;
            return;
        }

        sigx_conv = *sigzero;
        sigy_conv = *sigy;

        auto* sourceZ = static_cast<TParameter<double>*>(
            input.Get("templateZ"));
        template_z_mm = sourceZ
            ? sourceZ->GetVal() * 10.
            : std::numeric_limits<double>::quiet_NaN();

        const size_t nPoints = std::min({
            t->size(), sigzero->size(), sigy->size()});
        if (nPoints == 0) {
            std::cerr << filename << " 中没有可绘制的点。" << std::endl;
            return;
        }

        const int colorIndex = nFiles == 1
            ? 0
            : fileIndex * (TColor::GetNumberOfColors() - 1) / (nFiles - 1);
        const int color = TColor::GetColorPalette(colorIndex);

        auto* graphX = new TGraph(
            static_cast<int>(nPoints), t->data(), sigzero->data());
        auto* graphY = new TGraph(
            static_cast<int>(nPoints), t->data(), sigy->data());
        graphX->SetLineColor(color);
        graphY->SetLineColor(color);
        graphX->SetLineWidth(2);
        graphY->SetLineWidth(2);
        meanX.Add(graphX, "L");
        meanY.Add(graphY, "L");
        graphs.push_back(graphX);
        graphs.push_back(graphY);

        const TString label = std::isfinite(template_z_mm)
            ? TString::Format("tier %d, z = %.3f mm", fileIndex, template_z_mm)
            : TString::Format("tier %d", fileIndex);
        legendX.AddEntry(graphX, label, "l");
        legendY.AddEntry(graphY, label, "l");

        matrixTree.Fill();
        std::cout << "读取 " << filename << std::endl;
    }

    output.cd();
    matrixTree.Write();

    TCanvas canvas("meanWaveforms", "Mean waveforms", 1200, 600);
    canvas.Divide(2, 1);
    canvas.cd(1);
    gPad->SetGrid();
    meanX.Draw("A");
    legendX.Draw();
    canvas.cd(2);
    gPad->SetGrid();
    meanY.Draw("A");
    legendY.Draw();
    canvas.Write();

    for (auto* graph : graphs) delete graph;
    std::cout << "已生成 " << outputFile
              << "，共 " << nFiles
              << " 层，平均波形画布已写入 meanWaveforms。"
              << std::endl;
}
