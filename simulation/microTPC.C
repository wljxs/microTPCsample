#include <Rtypes.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2.h>
#include <TSystem.h>
#include <TTree.h>
#include <TRandom3.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Garfield/AvalancheGrid.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
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


using namespace Garfield;

int main(int argc, char *argv[]) {
  int runid;
  if(argc > 1) {
    runid = std::atoi(argv[1]);
  } else {
    runid = 0;//0表示测试，1表示正式运行
  }
  // TApplication app("app", &argc, argv);
  std::string detector = "microTPC5mmsmall";
  TRandom3 rnd(0);
  // MediumMagboltz gas("ar", 95., "ic4h10", 5.);
  // gas.SetTemperature(293.15);
  // gas.SetPressure(760.);
  MediumMagboltz gas;
  gas.LoadGasFile("../" + detector + "/ar_45_cf4_40_co2_15100-100000.gas");

  constexpr double rPenning = 0.20;
  constexpr double lambdaPenning = 0.;
  gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");

  constexpr double pitch = 0.04;
  constexpr unsigned int nEvents = 100;
  gas.SetMaxElectronEnergy(200.);
  std::string U = std::to_string(500);
  // Load the field map.

  ComponentComsol fm;
  //fm.SetImportRange(0, 0.6*pitch, 0, 1.1*pitch, 0, 1.01);
  fm.Initialise(("../" + detector + "/mesh.mphtxt").c_str(),
                ("../" + detector + "/mplist.txt").c_str(),
                ("../" + detector + "/field.txt").c_str(),
                "mm");

                
  fm.EnablePeriodicityX();
  fm.EnablePeriodicityY();
  fm.PrintRange();
  std::string xlabel = "xlabel", ylabel = "ylabel";
  std::string xlabel1 = "xlabeltest1",xlabel2 = "xlabeltest2";
  std::string xlabelminus1 = "xlabeltest-1";
  std::string xlabelminus2 = "xlabeltest-2";  
  std::string ylabelminus1 = "ylabeltest-1";
  std::string ylabel1 = "ylabeltest1";
  fm.SetWeightingPotential("../" + detector + "/xlabel.txt", xlabel);
  fm.SetWeightingPotential("../" + detector + "/xlabel1.txt", xlabel1);
  fm.SetWeightingPotential("../" + detector + "/xlabel2.txt", xlabel2);
  fm.SetWeightingPotential("../" + detector + "/xlabel-2.txt", xlabelminus2);
  fm.SetWeightingPotential("../" + detector + "/xlabel-1.txt", xlabelminus1);

  // fm.CopyWeightingPotential("xlabel-1", xlabel, -1./2*pitch, 0, 0, 0, 0, 0);
  // fm.CopyWeightingPotential("xlabel-2", xlabel, -1* pitch, 0, 0, 0, 0, 0);
  // fm.CopyWeightingPotential("xlabel-3", xlabel, -1.5*pitch, 0, 0, 0, 0, 0);
  // fm.CopyWeightingPotential("xlabel-4",xlabel,-2*pitch,0,0,0,0,0);
  // fm.CopyWeightingPotential("xlabel+1", xlabel, 1./2*pitch, 0, 0, 0, 0, 0);
  // fm.CopyWeightingPotential("xlabel+2", xlabel, 1 * pitch, 0, 0, 0, 0, 0);
  // fm.CopyWeightingPotential("xlabel+3", xlabel, 1.5*pitch, 0, 0, 0, 0, 0);
  // fm.CopyWeightingPotential("xlabel+4",xlabel,2*pitch,0,0,0,0,0);
  fm.SetWeightingPotential("../" + detector + "/ylabel.txt", ylabel);
  fm.SetWeightingPotential("../" + detector + "/ylabel1.txt", ylabel1);
  fm.SetWeightingPotential("../" + detector + "/ylabel-1.txt", ylabelminus1);
  // fm.CopyWeightingPotential("ylabel-1", ylabel, 0, -1.*pitch, 0, 0, 0, 0);
  // fm.CopyWeightingPotential("ylabel-2", ylabel, 0, -2.*pitch, 0, 0, 0, 0);
  // fm.CopyWeightingPotential("ylabel+1", ylabel, 0, 1.*pitch, 0, 0, 0, 0);
  // fm.CopyWeightingPotential("ylabel+2", ylabel, 0, 2.*pitch, 0, 0, 0, 0);

  // Dimensions of the microTPC [cm]

  fm.SetGas(&gas);
  fm.PrintMaterials();

  ViewField fieldView(&fm);
  ViewFEMesh meshView(&fm);
  constexpr bool plotField = false;



  // Create the sensor.
  Sensor sensor(&fm);
  sensor.SetArea(0 * pitch, 0 * pitch, -0.003, 3. * pitch, 3 * pitch, 0.501);
  std::vector<std::string> buffer(nEvents);
  sensor.AddElectrode(&fm, xlabel);
  sensor.AddElectrode(&fm, ylabel);
  sensor.AddElectrode(&fm, xlabel1);
  sensor.AddElectrode(&fm, xlabel2);
  sensor.AddElectrode(&fm, xlabelminus1);
  sensor.AddElectrode(&fm, xlabelminus2);
  sensor.AddElectrode(&fm, ylabel1);
  sensor.AddElectrode(&fm, ylabelminus1);

  const unsigned int nTimeBins = 1400;
  const double tmin = 0.;
  const double tmax = 700;
  const double tstep = (tmax - tmin) / nTimeBins;
  sensor.SetTimeWindow(tmin, tstep, nTimeBins);
  Shaper shaper(1, 50., 1., "unipolar");
  sensor.SetTransferFunction(shaper);
  constexpr bool plotSignal = true;
  ViewSignal *signalView = nullptr;
  TCanvas *cSignal = nullptr;
  TCanvas *cSignal1 = nullptr;
  if (plotSignal) {
    cSignal = new TCanvas("cSignal", "", 600, 600);
    cSignal1 = new TCanvas("cSignal1", "", 600, 600);
    signalView = new ViewSignal(&sensor);
    signalView->SetCanvas(cSignal);
  }

  AvalancheMicroscopic aval(&sensor);

  ViewDrift driftView;
  constexpr bool plotDrift = false;
  if (plotDrift) {
    aval.EnablePlotting(&driftView, 10);
  }
  
  AvalancheMC drift(&sensor);
  //drift.SetTimeSteps(50);
  drift.SetTimeWindow(0, 1400);
  drift.SetDistanceSteps(1.e-3);
  //drift.EnableRKFSteps(true);
  TFile *myfile = new TFile(("run/result/avalanche"+std::to_string(runid)+".root").c_str(), "RECREATE");
  TTree *tree1 = new TTree("tree", "Avalanche data");
  int event, ne, ni;
  std::vector<double> t, sigzero, sigy,sigtest1,sigtest2,sigtestminu1,sigtestminus2,sigytest1,sigytestminus1;
  tree1->Branch("event", &event, "event/I");
  tree1->Branch("t", &t);
  tree1->Branch("sigzero", &sigzero);
  tree1->Branch("sigy", &sigy);
  tree1->Branch("sigxtest1", &sigtest1);
  tree1->Branch("sigxtest2", &sigtest2);
  tree1->Branch("sigxtestminus1", &sigtestminu1);
  tree1->Branch("sigxtestminus2", &sigtestminus2);
  tree1->Branch("sigytest1", &sigytest1);
  tree1->Branch("sigytestminus1", &sigytestminus1);


  TTree *tree_after_conv = new TTree("tree_after_conv", "The tree_after_conv Title");
  std::vector<double> sigzero_after_conv, sigy_after_conv,sigytest1_after_conv,sigytestminus1_after_conv,sigtest1_after_conv,sigtest2_after_conv,sigtestminu1_after_conv,sigtestminus2_after_conv;
  tree_after_conv->Branch("event", &event, "event/I");
  tree_after_conv->Branch("ne", &ne, "ne/I");
  tree_after_conv->Branch("ni", &ni, "ni/I");
  tree_after_conv->Branch("t", &t);
  
  tree_after_conv->Branch("sigzero", &sigzero_after_conv);
  tree_after_conv->Branch("sigxtest1_after_conv", &sigtest1_after_conv);
  tree_after_conv->Branch("sigxtest2_after_conv", &sigtest2_after_conv);
  tree_after_conv->Branch("sigxtestminus1_after_conv", &sigtestminu1_after_conv);
  tree_after_conv->Branch("sigxtestminus2_after_conv", &sigtestminus2_after_conv);

  tree_after_conv->Branch("sigy", &sigy_after_conv);  
  tree_after_conv->Branch("sigytest1_after_conv", &sigytest1_after_conv);
  tree_after_conv->Branch("sigytestminus1_after_conv", &sigytestminus1_after_conv);


  TTree *tree2 = new TTree("tree2", "The Tree2 Title");
  Float_t ex1, ey1, ez1, et1, ee1, ex2, ey2, ez2, et2, ee2, ix1, iy1, iz1, it1,
      ix2, iy2, iz2, it2;
  std::vector<float> xyzte;
  tree2->Branch("event", &event, "event/I");
  tree2->Branch("ex1", &ex1, "ex1/F");
  tree2->Branch("ey1", &ey1, "ey1/F");
  tree2->Branch("ez1", &ez1, "ez1/F");
  tree2->Branch("et1", &et1, "et1/F");
  tree2->Branch("ee1", &ee1, "ee1/F");
  tree2->Branch("it2", &it2, "it2/F");
  int autosaveEvery = 10;

  double x0, y0, z0, t0, e0;
  TTree *tree3 = new TTree("tree3", "The Tree3 Title");
  tree3->Branch("event", &event, "event/I");
  tree3->Branch("x0", &x0, "x0/D");
  tree3->Branch("y0", &y0, "y0/D");
  tree3->Branch("z0", &z0, "z0/D");
  tree3->Branch("t0", &t0, "t0/D");
  tree3->Branch("e0", &e0, "e0/D");
  for (unsigned int i = 0; i < nEvents; ++i) {
    const auto event_start_time = std::chrono::high_resolution_clock::now();
    
    t.clear();
    sigzero.clear();
    sigy.clear();
    sigzero_after_conv.clear();
    sigy_after_conv.clear();
    sigtest1.clear();
    sigtestminus2.clear();
    sigtestminu1.clear();
    sigtest1_after_conv.clear();
    sigtestminus2_after_conv.clear();
    sigtestminu1_after_conv.clear();
    sigytest1.clear();
    sigytestminus1.clear();
    sigytest1_after_conv.clear();
    sigytestminus1_after_conv.clear();
    ne = 0;
    ni = 0;

    x0 = 0.008*rnd.Uniform(-0.5, 0.5)+0.06;
    y0 = pitch*rnd.Uniform(-0.5, 0.5)+0.06;
  
    z0 = 0.01;
    t0 = 0.;
    e0 = 0.1;
    sensor.ClearSignal();
    if (plotDrift) driftView.Clear();  // 清空上一事件的漂移线
    std::cout << "x0 = " << x0 << ", y0 = " << y0 << ", z0 = " << z0 << std::endl;

    aval.AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    auto np = aval.GetNumberOfElectronEndpoints();
    std::cout << "np = " << np << std::endl;
    event = i+runid*nEvents;
    for (int j = 0; j < np; j++) {
      double xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, xi1, yi1, zi1, ti1,
          xi2, yi2, zi2, ti2;
      int status;
      aval.GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2,
                               e2, status);
      ex1 = xe1;
      ey1 = ye1;
      ez1 = ze1, et1 = te1;
      ee1 = e1;
      ex2 = xe2;
      ey2 = ye2;
      ez2 = ze2, et2 = te2;
      ee2 = e2;
      drift.DriftIon(xe1, ye1, ze1, te1);
      drift.GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, status);
      tree2->Fill();
    }
    
    aval.GetAvalancheSize(ne, ni);
    //driftView.Plot();
    std::cout << ne << "," << ni << std::endl;

    for (int j = 0; j < nTimeBins; j++) {
      if(j%10!=0) continue;
      if (j % 10000 == 0) std::cout << "Time bin: " << j << "/" << nTimeBins << std::endl;
      t.push_back(j * tstep + tmin);
      sigzero.push_back(sensor.GetSignal("xlabel", j));
      sigy.push_back(sensor.GetSignal(ylabel, j));
      sigtest1.push_back(sensor.GetSignal(xlabel1, j));
      sigtest2.push_back(sensor.GetSignal(xlabel2, j));
      sigtestminu1.push_back(sensor.GetSignal(xlabelminus1, j));
      sigtestminus2.push_back(sensor.GetSignal(xlabelminus2, j));

      sigytest1.push_back(sensor.GetSignal(ylabel1, j));
      sigytestminus1.push_back(sensor.GetSignal(ylabelminus1, j));

    }
    sensor.ConvoluteSignals();

    for (int j = 0; j < nTimeBins; j++) {
      if(j%10!=0) continue;
      sigzero_after_conv.push_back(sensor.GetSignal("xlabel", j));
      sigtest1_after_conv.push_back(sensor.GetSignal(xlabel1, j));
      sigtest2_after_conv.push_back(sensor.GetSignal(xlabel2, j));
      sigtestminu1_after_conv.push_back(sensor.GetSignal(xlabelminus1, j));
      sigtestminus2_after_conv.push_back(sensor.GetSignal(xlabelminus2, j));


      sigy_after_conv.push_back(sensor.GetSignal(ylabel, j));      
      sigytest1_after_conv.push_back(sensor.GetSignal(ylabel1, j));
      sigytestminus1_after_conv.push_back(sensor.GetSignal(ylabelminus1, j));
    }

    tree1->Fill();
    tree_after_conv->Fill();    
    tree3->Fill();

    if (plotDrift) {
      TCanvas *cd = new TCanvas(Form("drift_event%d", i), "", 800, 600);
      constexpr bool plotMesh = true;
      if (plotMesh) {
        constexpr bool twod = true;
        meshView.SetPlane(0, -1, 0, 0, 0, 0);  // 法向沿 y -> 显示 xz 平面
        meshView.SetArea(0. * pitch, -0.05, 3. * pitch, 0.51);
        meshView.SetFillMesh(true);
        meshView.SetColor(0, kGray);
        // Set the color of the kapton.
        meshView.SetColor(2, kYellow + 3);
        meshView.EnableAxes();
        meshView.SetViewDrift(&driftView);
        meshView.SetCanvas(cd);
        const bool outline = false;
        meshView.Plot(twod, outline);
      }
      cd->Write();
      delete cd;
    }
    const auto event_end_time = std::chrono::high_resolution_clock::now();
    const auto event_duration = std::chrono::duration_cast<std::chrono::milliseconds>(event_end_time - event_start_time);
    std::cout << "Event " << i << " took " << event_duration.count() << " ms" << std::endl;
  }

  if (plotField) {
    // Set the normal vector of the viewing plane (xz plane).
    fieldView.SetPlane(0, -1, 0, 0, 0, 0);
    fieldView.SetNumberOfContours(50);
    // Set the plot limits in the current viewing plane.
    fieldView.SetArea(0. * pitch, -0.05, 3. * pitch, 0.51);
    fieldView.SetVoltageRange(-300,500.);
    TCanvas *cf = new TCanvas("cf", "", 600, 600);
    cf->SetLeftMargin(0.16);
    fieldView.SetCanvas(cf);
    fieldView.PlotContour("v");
    TCanvas *cf1 = new TCanvas("cf1", "", 600, 600);
    cf1->SetLeftMargin(0.16);
    fieldView.SetCanvas(cf1);
    fieldView.PlotContourWeightingField(xlabel, "v");

    TCanvas *cf2 = new TCanvas("cf2", "", 600, 600);
    cf2->SetLeftMargin(0.16);
    fieldView.SetCanvas(cf2);
    fieldView.PlotContourWeightingField(ylabel, "v");

    // TCanvas *cf3 = new TCanvas("cf3", "", 600, 600);
    // cf3->SetLeftMargin(0.16);
    // fieldView.SetCanvas(cf3);
    // fieldView.PlotContourWeightingField("xlabel-1", "v");

    TCanvas *cf3 = new TCanvas("cf3", "", 600, 600);
    cf3->SetLeftMargin(0.16);
    fieldView.SetCanvas(cf3);
    fieldView.PlotContourWeightingField("xlabel-1", "v");
    cf3->SaveAs("xlabel-1_weighting_field.png");
    TCanvas *cf4 = new TCanvas("cf4", "", 600, 600);
    cf4->SetLeftMargin(0.16);
    fieldView.SetCanvas(cf4);
    fieldView.PlotContourWeightingField("xlabel+1", "v");
    cf4->SaveAs("xlabel+1_weighting_field.png");

    meshView.SetPlane(0, -1, 0, 0, 0, 0);
    meshView.SetArea(0. * pitch, -0.05, 3. * pitch, 0.51);
    meshView.SetCanvas(cf);
    meshView.SetFillMesh(true);
    meshView.SetColor(0, kWhite);
    meshView.SetColor(1,kRed);
    meshView.SetColor(2, kGray);
    meshView.SetColor(3,kYellow);
    myfile->WriteTObject(cf);
    myfile->WriteTObject(cf1);
    myfile->WriteTObject(cf2);
    myfile->WriteTObject(cf3);
    myfile->WriteTObject(cf4);
    meshView.Plot(true);

    // ===== 绘制两个权场的差 =====
    // 设定采样范围和分辨率
    const int nBinsX = 200;
    const int nBinsZ = 200;
    const double xMin = 0. * pitch, xMax = 3. * pitch;
    const double zMin = -0.05, zMax = 0.501;
    const double yFix = 1.; // 固定在 y=0 平面

    // 选择要比较的两个权场名称
    std::string wf1 = "xlabel-1";      // 第一个权场
    std::string wf2 = "xlabeltest-1";     // 第二个权场

    TH2F *hDiff = new TH2F("hDiff",
                           Form("Weighting Field Difference: %s - %s ;x [cm];z [cm];#DeltaV",
                                wf1.c_str(), wf2.c_str()),
                           nBinsX, xMin, xMax, nBinsZ, zMin, zMax);

    for (int ix = 0; ix <= nBinsX; ++ix) {
      double x = xMin + ix * (xMax - xMin) / nBinsX;
      for (int iz = 0; iz <= nBinsZ; ++iz) {
        double z = zMin + iz * (zMax - zMin) / nBinsZ;

        double wx1, wy1, wz1, wx2, wy2, wz2;
        std::string active1, active2;
        // 采样第一个权场
        fm.WeightingField( x, yFix, z, wx1, wy1, wz1,wf1);
        // 采样第二个权场
        fm.WeightingField( x, yFix, z, wx2, wy2, wz2, wf2);

        // 用权场的电势值 (标量) 做差；也可改用 wx 分量等
        double diff = wz1 - wz2;  // 或改成 wx, wy, wz 任一分量
        hDiff->SetBinContent(ix + 1, iz + 1, diff);
      }
    }

    TCanvas *cDiff = new TCanvas("cDiff", "Weighting Field Difference", 700, 600);
    cDiff->SetLeftMargin(0.14);
    cDiff->SetRightMargin(0.14);
    hDiff->SetStats(0);
    hDiff->Draw("colz");
    cDiff->SaveAs(Form("diff_%s_%s.png", wf1.c_str(), wf2.c_str()));
    myfile->WriteTObject(hDiff);
    myfile->WriteTObject(cDiff);
  }

  myfile->WriteTObject(tree1);
  myfile->WriteTObject(tree2);
  myfile->WriteTObject(tree3);
  myfile->WriteTObject(tree_after_conv);
  myfile->Close();
  std::cout << "Finished!" << std::endl;
  // app.Run(false);
}