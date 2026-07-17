#include "Script/DUTAnalysisScript.h"
#include "Detector/DetectorFactory.h"
#include "Event/DataModel.h"
#include "Event/DetectorFrame.h"
#include "RtypesCore.h"
#include "Script/Base/ScriptFactory.h"
#include "Script/Base/ScriptManager.h"
#include "Script/Base/RawDataParser.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TTree.h>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>

using namespace std;
bool IsPredOnBadStrip(const Detector& det, const TVector3& predL) {
    const auto* mask = det.GetBadStripMask();
    if (!mask) return false;

    const auto& cfg = det.getConfig();

    for (int type : cfg.readoutPlaneType) {
        auto itPitch = cfg.readoutPlanePitch.find(type);
        auto itNum = cfg.readoutPlaneStripNumber.find(type);
        if (itPitch == cfg.readoutPlanePitch.end() || itNum == cfg.readoutPlaneStripNumber.end()) continue;

        const double pitch = itPitch->second;
        const int nStrips = itNum->second;
        if (pitch <= 0 || nStrips <= 0) continue;

        double coord = 0.0;
        if (type == DUTAnalysisConfig::kTypeX) {
            coord = predL.X();
        } else if (type == DUTAnalysisConfig::kTypeY) {
            coord = predL.Y();
        } else {
            continue;
        }

        // Strip IDs in this framework are effectively 1..N (consistent with bad_strips.json)
        const int stripID = static_cast<int>(std::lround(coord / pitch));
        if (stripID < 1 || stripID > nStrips) continue;

        if (mask->IsBad(type, stripID)) return true;
    }
    return false;
}
void RunDUTAlign(const std::vector<Event>& events, std::shared_ptr<Detector> detector, int detID, bool firstAlign);
LocalHit CalcuDutResidual(std::shared_ptr<Detector> detector, const std::vector<Cluster>& clusters, const TVector3& predL, double& residualX, double& residualY,bool isAlign=false);
double DUTChi2Objective(const double* par, const std::vector<Event>& events, std::shared_ptr<Detector> detector, int detID);
double DUTChi2Objective(const double* par, const std::vector<Event>& events, std::shared_ptr<Detector> detector, int detID, bool isFirstAlign);

void DUTAnalysisScript::LoadConfig(const json& config) {
    m_runAlignment = config.value("runAlignment", true);
    m_saveNoiseData = config.value("saveNoiseData", true);
    m_progressInterval = config.value("progressInterval", 1000);
    m_maxEvents = config.value("maxEvents", -1);
}

void DUTAnalysisScript::Print() const {
    cout << "DUTAnalysisScript Configuration:" << endl;
    cout << "  Run Alignment: " << (m_runAlignment ? "Yes" : "No") << endl;
    cout << "  Save Noise Data: " << (m_saveNoiseData ? "Yes" : "No") << endl;
    cout << "  Progress Interval: " << m_progressInterval << endl;
    cout << "  Max Events: " << (m_maxEvents > 0 ? to_string(m_maxEvents) : "All") << endl;
}

bool DUTAnalysisScript::Execute() {
    auto t0 = chrono::high_resolution_clock::now();
   // m_runAlignment = false; // 临时关闭对齐功能，后续完善后再打开
    cout << "\n========================================" << endl;
    cout << " DUT Analysis Script" << endl;
    cout << "========================================" << endl;

    auto parser = GetParser();
    if (!parser) {
        cerr << "Error: Parser not set!" << endl;
        return false;
    }

    auto& factory = DetectorFactory::GetInstance();

    // 加载track信息
    string trackFile = GetOutputDir() + "TrackInfo.root";
    cout << "\nLoading track info..." << endl;
    cout << "File: " << trackFile << endl;

    TFile* f = TFile::Open(trackFile.c_str(), "READ");
    if (!f || f->IsZombie()) {
        cerr << "Error: Cannot open " << trackFile << endl;
        cerr << "Please run Track Analysis first!" << endl;
        return false;
    }

    TTree* trackTree = (TTree*)f->Get("Tracks");
    if (!trackTree) {
        cerr << "Error: No Tracks tree!" << endl;
        f->Close();
        return false;
    }

    Int_t eventID;
    Track* track = nullptr;
    double sigTime;
    trackTree->SetBranchAddress("eventID", &eventID);
    trackTree->SetBranchAddress("track", &track);
    trackTree->SetBranchAddress("t0", &sigTime);

    cout << "\nProcessing DUT data..." << endl;
    std::vector<Event> events;  // Script本地数据
    std::vector<Event> single;  // 用于发现单个cluster（其中cluster_size要大于2）的事例，

    int processed = 0;
    Long64_t nEntries = trackTree->GetEntries();
    if (m_maxEvents > 0 && nEntries > m_maxEvents) {
        nEntries = m_maxEvents;
    }

    const auto& duts = factory.GetDetectorsByRole(Detector::Role::DUT);

    // 只读一次 pedestals.apv_pedstd，并用 shared_ptr 供每个 DetectorFrame 引用
    // 结构：detID -> shared_ptr(type -> [strip1..N] pedstd)
    std::map<int, std::shared_ptr<const std::map<int, std::vector<float>>>> pedStdByDet;
    for (const auto& det : duts) {
        const int detID = det->GetID();
        auto pedByType = std::make_shared<std::map<int, std::vector<float>>>();
        for (int type : det->getConfig().readoutPlaneType) {
            (*pedByType)[type] = parser->GetPedstd(detID, type);
        }
        pedStdByDet[detID] = pedByType;
    }

    string dutFile = GetOutputDir() + "DUTInfo.root";
    TFile* fDut = new TFile(dutFile.c_str(), "RECREATE");
    TTree* tDut = new TTree("DUTTree", "DUT data");

    TTree* noiseTree = nullptr;
    int noise_eventID, noise_stripID, noise_stripType;
    double noise_value;


    struct NoiseStat {
        long long n = 0;
        double sum = 0.0;
        double sum2 = 0.0;
    };
    map<pair<int, int>, NoiseStat> noiseStats;

    TTree *tDutori = new TTree("DUTTreeOri", "Original DUT data before alignment");
    Int_t dutID;
    vector<Cluster> clusters;
     tDutori->Branch("eventID", &eventID);
     tDutori->Branch("dutID", &dutID);
     tDutori->Branch("clusters", &clusters);

    if (m_saveNoiseData) {
        noiseTree = new TTree("Noise", "Per-strip noise samples");
        noiseTree->Branch("eventID", &noise_eventID);
        noiseTree->Branch("stripID", &noise_stripID);
        noiseTree->Branch("stripType", &noise_stripType);
        noiseTree->Branch("noise", &noise_value);
    }

    // 加载事件并处理DUT数据
    for (Long64_t i = 0; i < nEntries; ++i) {
        trackTree->GetEntry(i);

        if (processed % m_progressInterval == 0) {
            cout << "\r  Processing: " << processed << "/" << nEntries << flush;
        }

        auto rawHits = parser->LoadEvent(eventID);
        if (rawHits.empty()) continue;

        Event evt{.eventID = int(eventID), .track = *track, .t0 = sigTime};
        int clusterCountx = 0;
        int clusterCounty = 0;
        bool issingleClusterEvent=false;

        for (auto& det : duts) {//这里目前默认是一个DUT,后续如果有多个DUT的话需要改成循环处理每个DUT的事件数据
            clusterCountx=0;
            issingleClusterEvent=false;
            const int id = det->GetID();

            auto detEvt = make_shared<DetectorFrame>(*det);
            detEvt->SetRawData(rawHits[det->GetID()]);

            // 注入每条 strip 的 pedstd，WaveformProcessor 会用它作为阈值（apv_pedstd）
            auto itPed = pedStdByDet.find(det->GetID());
            if (itPed != pedStdByDet.end()) {
                detEvt->SetStripPedestalStd(itPed->second);
            }
            detEvt->Process(evt.t0);
            // 保存对齐前的 DUT 数据（主要是 cluster 信息）
            dutID = id;
            clusters = detEvt->Clusters();
            tDutori->Fill();

            for(auto cluster : detEvt->Clusters()){
                if(cluster.maxAmp>200 && cluster.type==DUTAnalysisConfig::kTypeX){
                     clusterCountx++;
                }
                if(cluster.maxAmp>200 && cluster.type==DUTAnalysisConfig::kTypeY){
                     clusterCounty++;
                }
            }
            if(clusterCountx==1 && clusterCounty==1){
                issingleClusterEvent=true;
            }

            evt.detectorFramesMap[id] = std::move(detEvt);

            if (m_saveNoiseData) {
                noise_eventID = eventID;
                for (const auto& raw : evt.detectorFramesMap[id]->Raw()) {
                    noise_stripID = raw.stripID;
                    noise_stripType = raw.type;
                    noise_value = raw.adc[0];
                    noiseTree->Fill();

                    auto key = make_pair(noise_stripType, noise_stripID);
                    auto& stat = noiseStats[key];//就是对每一个条的噪声进行统计
                    stat.n++;
                    stat.sum += noise_value;
                    stat.sum2 += noise_value * noise_value;
                }
            }
        }
        if(issingleClusterEvent){
        single.push_back(evt);  // 添加到单事件调试容器
        }        
        events.push_back(std::move(evt));

        processed++;
    }
    cout << endl;
    cout << "  Processed " << events.size() << " DUT events," << single.size() << " single-cluster events" << endl;

    // 运行DUT对齐
    std::vector<Event> alignedEvents,alignedsingleEvents;
    int badEventCount = 0;
    if (m_runAlignment) {
        for (auto& det : duts) {
            RunDUTAlign(single, det, det->GetID(), true);
          //  RunDUTAlign(events, det, det->GetID(), true);

            for(auto det : duts){
                int id = det->GetID();
                for(auto evt : events){
                eventID = evt.eventID;

                TVector3 predG = det->CalcHitFromTrack(evt.track);
                TVector3 predL = det->GlobalToLocal(predG);

            // 跳过落在坏道上的预测位置
            if (IsPredOnBadStrip(*det, predL)) {
                // badEventCount++;
                continue;
            }
            alignedEvents.push_back(evt);
            }
                for(auto evt : single){
                    eventID = evt.eventID;
    
                    TVector3 predG = det->CalcHitFromTrack(evt.track);
                    TVector3 predL = det->GlobalToLocal(predG);
            if (IsPredOnBadStrip(*det, predL)) {
                //std::cout << "Event " << eventID << ", DUT " << id << ": Prediction on bad strip, skipping." << std::endl;
                continue;
            }
                alignedsingleEvents.push_back(evt);
            }
            }
            std::cout << "after alignment, aligned events: " << alignedEvents.size() << ", aligned single-cluster events: " << alignedsingleEvents.size() << std::endl;
            RunDUTAlign(alignedsingleEvents, det, det->GetID(), true);
            RunDUTAlign(alignedsingleEvents, det, det->GetID(), false);
           // RunDUTAlign(events, det, det->GetID(), false);
        }
    }

    // 初始化配置和统计
    DUTAnalysisConfig analysisConfig;
    DUTStatistics statistics(analysisConfig);

    cout << "\nSaving DUT results..." << endl;

    // 保存噪声统计
    if (m_saveNoiseData) {
        TTree* noiseStatTree = new TTree("NoiseStat", "Per-strip noise variance");
        int stat_stripID, stat_stripType;
        long long stat_n;
        double stat_mean, stat_variance, stat_sigma;

        noiseStatTree->Branch("stripID", &stat_stripID);
        noiseStatTree->Branch("stripType", &stat_stripType);
        noiseStatTree->Branch("nSamples", &stat_n);
        noiseStatTree->Branch("mean", &stat_mean);
        noiseStatTree->Branch("variance", &stat_variance);
        noiseStatTree->Branch("sigma", &stat_sigma);

        for (const auto& [key, stat] : noiseStats) {
            stat_stripType = key.first;
            stat_stripID = key.second;
            stat_n = stat.n;

            if (stat.n < 2) continue;

            stat_mean = stat.sum / stat.n;
            stat_variance = (stat.sum2 - stat.n * stat_mean * stat_mean) / (stat.n - 1);
            stat_sigma = sqrt(stat_variance);

            noiseStatTree->Fill();
        }

        noiseTree->Write();
        noiseStatTree->Write();
    }

    // DUT数据分支
    Double_t resX, resY, predX, predY;
    Double_t hitX, hitY;
    vector<Int_t> clusterIndices;
    vector<StripHit> stripHits;
    
    Int_t hitFlag;
    Cluster clusterX, clusterY;
    vector<StripHit> stripHitsX, stripHitsY;
    TH2D *h2hitpred = new TH2D("h2hitpred", "Hit vs Prediction;PredX (mm);PredY (mm)", 200, 0, 100, 200, 0, 100);

    tDut->Branch("eventID", &eventID);
    tDut->Branch("dutID", &dutID);
    tDut->Branch("predX", &predX);
    tDut->Branch("predY", &predY);
    tDut->Branch("resX", &resX);
    tDut->Branch("resY", &resY);
    tDut->Branch("hitX", &hitX);
    tDut->Branch("hitY", &hitY);
    tDut->Branch("clusterIndices", &clusterIndices);
    tDut->Branch("stripHits", &stripHits);
    tDut->Branch("clusters", &clusters);
    tDut->Branch("hitFlag", &hitFlag);
    tDut->Branch("clusterX", &clusterX);
    tDut->Branch("clusterY", &clusterY);
    tDut->Branch("stripHitsX", &stripHitsX);
    tDut->Branch("stripHitsY", &stripHitsY);
    TH1D *hResX = new TH1D("hResX", "Residual X;resX (mm);Entries", 5000, -10, 10);
    TH1D *hResY = new TH1D("hResY", "Residual Y;resY (mm);Entries", 5000, -10, 10);
    TH2D *hRes2D = new TH2D("hRes2D", "Residual 2D;resX (mm);resY (mm)", 1000, -2, 2, 1000, -2, 2);
    TH2D *hhit2D = new TH2D("hhit2D", "Hit position 2D;hitX (mm);hitY (mm)", 200,0,100,200,0,100);
    TF1 *fResX = new TF1("fResX", "gaus", -5, 5);
    TF1 *fResY = new TF1("fResY", "gaus", -5, 5);
    int effectiveHitCount = 0;
    // 填充DUT数据
    for (auto& det : duts) {
        int id = det->GetID();
        for (auto& evt : events) {
            eventID = evt.eventID;
            dutID = id;

            TVector3 predG = det->CalcHitFromTrack(evt.track);
            TVector3 predL = det->GlobalToLocal(predG);
            predX = predL.X();
            predY = predL.Y();
            h2hitpred->Fill(predX, predY);

            // 跳过落在坏道上的预测位置
            if (IsPredOnBadStrip(*det, predL)&& predY<45) {
                badEventCount++;
                //std::cout << "Event " << eventID << ", DUT " << id << ": Prediction on bad strip, skipping." << std::endl;
                continue;
            }

            auto frameIt = evt.detectorFramesMap.find(id);
            if (frameIt != evt.detectorFramesMap.end()) {
                const auto& detFrame = frameIt->second;

                stripHits = detFrame->StripHits();
                clusters = detFrame->Clusters();

                if (!clusters.empty()) {
                    LocalHit localHit = CalcuDutResidual(det, clusters, predL, resX, resY);
                    const auto& clusterIdx = localHit.clusterIndices;

                    int idxX = clusterIdx[0];//对应X方向的cluster索引,-1表示没有
                    int idxY = clusterIdx[1];//对应Y方向的cluster索引,-1表示没有

                    clusterX = (idxX >= 0) ? clusters[idxX] : CreateInvalidCluster(DUTAnalysisConfig::kTypeX);
                    clusterY = (idxY >= 0) ? clusters[idxY] : CreateInvalidCluster(DUTAnalysisConfig::kTypeY);

                    hitX = localHit.localPos.X();
                    hitY = localHit.localPos.Y();

                    bool hasX = (idxX >= 0);
                    bool hasY = (idxY >= 0);
                    if (hasX && hasY)
                        hitFlag = 3;
                    else if (hasX)
                        hitFlag = 1;
                    else if (hasY)
                        hitFlag = 2;
                    else
                        hitFlag = 0;

                    clusterIndices = clusterIdx;

                    auto [binX, binY] = statistics.GetBinIndices(predX, predY);
                    bool hasValidHit = (idxX >= 0 && idxY >= 0);
                    statistics.AddBinData(id, binX, binY, hasValidHit, resX, resY);
                    if(predY<45){
                        effectiveHitCount++;  
                    }
                                      
                    if (hasValidHit) {
                    if( predY<45){
                        hResX->Fill(resX);
                        hResY->Fill(resY);
                        hRes2D->Fill(resX, resY);
                        hhit2D->Fill(hitX, hitY);
                    }
                    }
                } else {
                    clusterX = CreateInvalidCluster(DUTAnalysisConfig::kTypeX);
                    clusterY = CreateInvalidCluster(DUTAnalysisConfig::kTypeY);
                    hitX = DUTAnalysisConfig::kInvalidValue;
                    hitY = DUTAnalysisConfig::kInvalidValue;
                    resX = DUTAnalysisConfig::kInvalidValue;
                    resY = DUTAnalysisConfig::kInvalidValue;
                    hitFlag = 0;
                    clusterIndices.clear();

                    auto [binX, binY] = statistics.GetBinIndices(predX, predY);
                    statistics.AddBinData(id, binX, binY, false, 0, 0);
                }
            }
        //   tDut->Fill();
            // if(hitFlag==3){
            //     tDut->Fill();}
        //    if(hitFlag==3&&predX>20&&predX<60&&predY>15&&predY<40){
        //     tDut->Fill();}
            if(predY<45){
                tDut->Fill();
            }
        }
    }
    hResX->Fit(fResX, "QR");
    hResY->Fit(fResY, "QR");

    //统计效率,认为在10倍sigma范围内的残差是有效hit
    int BedetectCount = 0;
    int BedetectCountx = 0;
    int BedetectCounty = 0;
    for(int entry=0;entry<tDut->GetEntries();entry++){
        tDut->GetEntry(entry);
        if(hitFlag==3 && resX<fResX->GetParameter(2)*10 && resX>-fResX->GetParameter(2)*10 && resY<fResY->GetParameter(2)*10 && resY>-fResY->GetParameter(2)*10){
            if(predY<45){
            BedetectCount++;
            }
        }
        if(resX<fResX->GetParameter(2)*10 && resX>-fResX->GetParameter(2)*10){
            if(predY<45){
                BedetectCountx++;
            }
        }
        if(resY<fResY->GetParameter(2)*10 && resY>-fResY->GetParameter(2)*10){
            if(predY<45){
                BedetectCounty++;
            }
        }
    }
    std::cout << hRes2D->GetCovariance(1,2) << std::endl;

    hResX->Write();
    hResY->Write();
    hRes2D->Write();
    hhit2D->Write();
    h2hitpred->Write();
    tDutori->Write();
    tDut->Write();
    fDut->Close();
    delete fDut;

    f->Close();
    delete f;

    cout << "DUT data saved to: " << dutFile << endl;
    cout << "Bad event count: " << badEventCount << endl;
    cout << "Effective hit count: " << effectiveHitCount << endl;
    cout << "BeDetect hit count: " << BedetectCount << endl;
    cout << "BeDetect hit count X: " << BedetectCountx << endl;
    cout << "BeDetect hit count Y: " << BedetectCounty << endl;

    auto t1 = chrono::high_resolution_clock::now();
    double sec = chrono::duration<double>(t1 - t0).count();

    cout << "\n========================================" << endl;
    cout << "DUT Analysis Complete" << endl;
    cout << "Time: " << fixed << setprecision(2) << sec << " seconds" << endl;
    cout << "========================================" << endl;

    return true;
}

void RunDUTAlign(const std::vector<Event>& events, std::shared_ptr<Detector> detector, int detID, bool firstAlign) {

    if (events.empty()) {
        std::cerr << "[DUT Alignment] No events to analyze for DUT " << detID << "!" << std::endl;
        return;
    }

    std::cout << "[DUT " << detID << "] Aligning (6-parameter)..." << std::endl;
    auto minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetTolerance(0.005);
    minimizer->SetPrintLevel(0);
    UInt_t nPar = firstAlign ? 3 : 6;//第一次对齐只调整x，y和rotz

    // 使用lambda捕获this和参数
    auto chi2Func = [&events, &detector, detID,firstAlign](const double* par) -> double {
         {if(firstAlign){return DUTChi2Objective(par, events, detector, detID,true);}  
          else{return DUTChi2Objective(par, events, detector, detID);}
         }};

    ROOT::Math::Functor f(chi2Func, nPar);
    minimizer->SetFunction(f);

    // 使用当前探测器对齐作为初始猜测，避免第二次对齐从(0,0,0)开始覆盖粗对齐结果
        TVector3 curAlignPos = detector->GetAlignPos();
        TVector3 curAlignRot = detector->GetAlignRot();

    if (firstAlign) {
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
        //minimizer->SetVariableLimits(0,curAlignPos.X()-2,curAlignPos.X()+2); // 限制dx在合理范围，避免过度调整
        //minimizer->SetVariableLimits(1,curAlignPos.Y()-2,curAlignPos.Y()+2); // 限制dy在合理范围，避免过度调整
        minimizer->SetVariableLimits(2,curAlignPos.Z()-10,curAlignPos.Z()+10); // 限制dz在合理范围，避免过度调整
        minimizer->SetVariableLimits(3, curAlignRot.X()-0.1, curAlignRot.X()+0.1); // 限制rotX在合理范围，避免过度调整
        minimizer->SetVariableLimits(4, curAlignRot.Y()-0.1, curAlignRot.Y()+0.1); // 限制rotY在合理范围，避免过度调整
        minimizer->SetVariableLimits(5, curAlignRot.Z()-0.1, curAlignRot.Z()+0.1); // 限制rotZ在合理范围，避免过度调整
    }
    minimizer->Minimize();

    // 应用结果
    if(firstAlign) {
        const double* result = minimizer->X();
        double dx = result[0];
        double dy = result[1];
        double rotZ = result[2];

        detector->SetAlignment(dx, dy, 0, 0, 0, rotZ);

        TVector3 pos = detector->GetPos();
        TVector3 rot = detector->GetRot();

        std::cout << "DUT " << detID << " coarse alignment: "
                  << std::fixed << std::setprecision(5)
                  << "\"position\": [" << pos.X() << "," << pos.Y() << "," << pos.Z() << "],"
                  << "\"rotation\": [" << rot.X() << "," << rot.Y() << "," << rot.Z() << "]"
                  << std::endl;
        std::cout << "DUT " << detID << " coarse alignment corrections: "
                  << "dx=" << dx << ", dy=" << dy << ", rotZ=" << rotZ
                  << std::endl;

        std::cout << "[DUT " << detID << "] Coarse alignment chi2: " << minimizer->MinValue() << std::endl;
    }
    else{
    const double* result = minimizer->X();
    double dx = result[0];
    double dy = result[1];
    double dz = result[2];
    double rotX = result[3];
    double rotY = result[4];
    double rotZ = result[5];

    detector->SetAlignment(dx, dy, dz, rotX, rotY, rotZ);

    TVector3 pos = detector->GetPos();
    TVector3 rot = detector->GetRot();

    std::cout << "DUT " << detID << " alignment: "
              << std::fixed << std::setprecision(5)
              << "\"position\": [" << pos.X() << "," << pos.Y() << "," << pos.Z() << "],"
              << "\"rotation\": [" << rot.X() << "," << rot.Y() << "," << rot.Z() << "]"
              << std::endl;
    std::cout << "DUT " << detID << " alignment corrections: "
              << "dx=" << dx << ", dy=" << dy << ", dz=" << dz << ", "
              << "rotX=" << rotX << ", rotY=" << rotY << ", rotZ=" << rotZ
              << std::endl;

    std::cout << "[DUT " << detID << "] Final chi2: " << minimizer->MinValue() << std::endl;

    delete minimizer;}   
}

// ========== DUT对齐私有方法 ==========

LocalHit CalcuDutResidual(std::shared_ptr<Detector> detector, const std::vector<Cluster>& clusters, const TVector3& predL, double& residualX, double& residualY,bool isAlign) {

    double predX = predL.X();
    double predY = predL.Y();

    // 从探测器配置读取参数，消除魔法数字
    const auto& config = detector->getConfig();
    const int typeX = DUTAnalysisConfig::kTypeX;
    const int typeY = DUTAnalysisConfig::kTypeY;

    // X方向处理：找到最优cluster
    int bestClusterXIndex = -1;
    double minResX = std::numeric_limits<double>::infinity();
    double bestPosX = DUTAnalysisConfig::kInvalidValue;
    if (config.readoutPlanePitch.find(typeX) != config.readoutPlanePitch.end()) {
        double pitchX = config.readoutPlanePitch.at(typeX);
        // bool isBadclusterX = false;
        for (size_t i = 0; i < clusters.size(); ++i) {
            if (clusters[i].type == typeX) {
                if(clusters[i].isBad && isAlign) continue; // 跳过包含坏道的cluster
                if(clusters[i].maxAmp<150) continue; // 分析的时候就考虑较大的cluster.
                if(isAlign && clusters[i].maxAmp<200) continue; // 对齐时只考虑振幅较大的cluster，减少噪声影响                
                double currentResX = std::abs(clusters[i].pos * pitchX - predX);
                if (currentResX < minResX) {
                    minResX = currentResX;
                    bestClusterXIndex = static_cast<int>(i);
                    bestPosX = clusters[i].pos * pitchX;
                    // if(clusters[i].isBad){
                    //     isBadclusterX = true;
                    // }
                }
            }
        }
        // if(isBadclusterX){
        //     bestClusterXIndex = -1; // 如果最优cluster包含坏道，则认为没有有效hit
        //     bestPosX = DUTAnalysisConfig::kInvalidValue;
        // }
    }

    int bestClusterYIndex = -1;
    double minResY = std::numeric_limits<double>::infinity();
    double bestPosY = DUTAnalysisConfig::kInvalidValue;
    if (config.readoutPlanePitch.find(typeY) != config.readoutPlanePitch.end()) {
        double pitchY = config.readoutPlanePitch.at(typeY);
        // bool isBadclusterY = false;
        for (size_t i = 0; i < clusters.size(); ++i) {
            if (clusters[i].type == typeY) {
                if(clusters[i].isBad && isAlign) continue; // 跳过包含坏道的cluster
                if(clusters[i].maxAmp<130) continue; // 这个阈值是根据单事件调试中cluster振幅分布观察到的，后续可以优化
                if(isAlign && clusters[i].maxAmp<200) continue; // 对齐时只考虑振幅较大的cluster，减少噪声影响
                double currentResY = std::abs(clusters[i].pos * pitchY - predY);
                if (currentResY < minResY) {
                    minResY = currentResY;
                    bestClusterYIndex = static_cast<int>(i);
                    bestPosY = clusters[i].pos * pitchY;
                    // if(clusters[i].isBad){
                    //     isBadclusterY = true;
                    // }
                }
            }
        }
        // if(isBadclusterY){
        //     bestClusterYIndex = -1; // 如果最优cluster包含坏道，则认为没有有效hit
        //     bestPosY = DUTAnalysisConfig::kInvalidValue;
        // }
    }

    // 构建LocalHit
    LocalHit localHit;

    if (bestClusterXIndex != -1) {
        residualX = bestPosX - predX;
    } else {
        bestPosX = DUTAnalysisConfig::kInvalidValue;
        residualX = 100;
    }

    if (bestClusterYIndex != -1) {
        residualY = bestPosY - predY;
    } else {
        bestPosY = DUTAnalysisConfig::kInvalidValue;
        residualY = 100;
    }

    localHit.localPos.SetXYZ(bestPosX, bestPosY, 0);
    localHit.clusterIndices = {bestClusterXIndex, bestClusterYIndex};
    return localHit;
}

double DUTChi2Objective(const double* par, const std::vector<Event>& events, std::shared_ptr<Detector> detector, int detID) {

    const double dx = par[0];
    const double dy = par[1];
    const double dz = par[2];
    const double rotX = par[3];
    const double rotY = par[4];
    const double rotZ = par[5];

    detector->SetAlignment(dx, dy, dz, rotX, rotY, rotZ);

    // 对所有事件求平均 χ²
    double chi2 = 0.0;
    int nEvents = 0;

    for (const auto& evt : events) {
        double residualX = 0.0;
        double residualY = 0.0;

        // 使用detectorFramesMap获取Clusters
        auto frameIt = evt.detectorFramesMap.find(detID);
        if (frameIt == evt.detectorFramesMap.end() || frameIt->second->Clusters().empty()) {
            continue;
        }

        const auto& clusters = frameIt->second->Clusters();
        TVector3 predG = detector->CalcHitFromTrack(evt.track);
        TVector3 predL = detector->GlobalToLocal(predG);

        double predX = predL.X();
        double predY = predL.Y();
        // if (predY < 55 || predY > 80 || predX < -30 || predX > -15) continue;
        // 调用新的CalcuDutResidual方法
        LocalHit localHit = CalcuDutResidual(detector, clusters, predL, residualX, residualY,true);

        double res = residualX * residualX + residualY * residualY;
        if (res > 3) continue;
        chi2 += res;
        nEvents++;
    }

    // 返回平均 χ²
    return (nEvents > 0) ? chi2 / nEvents : 1e9;
}


double DUTChi2Objective(const double* par, const std::vector<Event>& events, std::shared_ptr<Detector> detector, int detID, bool isFirstAlign) {

    const double dx = par[0];
    const double dy = par[1];
    const double rotZ = par[2];

    detector->SetAlignment(dx, dy, 0, 0, 0, rotZ);

    // 对所有事件求平均 χ²
    double chi2 = 0.0;
    int nEvents = 0;

    for (const auto& evt : events) {
        double residualX = 0.0;
        double residualY = 0.0;

        // 使用detectorFramesMap获取Clusters
        auto frameIt = evt.detectorFramesMap.find(detID);
        if (frameIt == evt.detectorFramesMap.end() || frameIt->second->Clusters().empty()) {
            continue;
        }

        const auto& clusters = frameIt->second->Clusters();
        TVector3 predG = detector->CalcHitFromTrack(evt.track);
        TVector3 predL = detector->GlobalToLocal(predG);

        double predX = predL.X();
        double predY = predL.Y();
        // if (predY < 55 || predY > 80 || predX < -30 || predX > -15) continue;
        // 调用新的CalcuDutResidual方法
        LocalHit localHit = CalcuDutResidual(detector, clusters, predL, residualX, residualY,true);

        double res = residualX * residualX + residualY * residualY;
        if (res > 3) continue;
        chi2 += res;
        nEvents++;
    }

    // 返回平均 χ²
    return (nEvents > 0) ? chi2 / nEvents : 1e9;
}

// ========== DUTStatistics方法实现 ==========

std::pair<int, int> DUTStatistics::GetBinIndices(double predX, double predY) const {
    const auto& binning = m_config.binning;

    if (predX < binning.predX_min || predX > binning.predX_max ||
        predY < binning.predY_min || predY > binning.predY_max) {
        return {-1, -1};
    }

    int binX = static_cast<int>((predX - binning.predX_min) /
                                (binning.predX_max - binning.predX_min) * binning.nBinsX);
    int binY = static_cast<int>((predY - binning.predY_min) /
                                (binning.predY_max - binning.predY_min) * binning.nBinsY);

    // 处理边界情况
    if (binX >= binning.nBinsX) binX = binning.nBinsX - 1;
    if (binY >= binning.nBinsY) binY = binning.nBinsY - 1;

    return {binX, binY};
}

void DUTStatistics::AddBinData(int dutID, int binX, int binY,
                               bool hasValidHit, double resX, double resY) {
    if (binX < 0 || binY < 0) return;

    auto& binData = m_binDataMap[dutID][{binX, binY}];
    binData.totalEvents++;

    if (hasValidHit) {
        binData.hitEvents++;
        binData.resX_values.push_back(resX);
        binData.resY_values.push_back(resY);
    }
}

// ========== DUTAnalysisScript静态方法实现 ==========

Cluster DUTAnalysisScript::CreateInvalidCluster(int type) {
    Cluster invalidCluster;
    invalidCluster.type = type;
    invalidCluster.size = DUTAnalysisConfig::kInvalidSize;
    invalidCluster.range = DUTAnalysisConfig::kInvalidSize;
    invalidCluster.charge = DUTAnalysisConfig::kInvalidValue;
    invalidCluster.maxAmp = DUTAnalysisConfig::kInvalidValue;
    invalidCluster.time = DUTAnalysisConfig::kInvalidValue;
    invalidCluster.centroid = DUTAnalysisConfig::kInvalidValue;
    invalidCluster.pos = DUTAnalysisConfig::kInvalidValue;
    invalidCluster.stripHitIndices.clear();
    return invalidCluster;
}

REGISTER_SCRIPT("DUTAnalysis", DUTAnalysisScript);
