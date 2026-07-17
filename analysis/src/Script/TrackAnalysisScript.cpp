#include "Script/TrackAnalysisScript.h"
#include "Algorithm/AnalysisUtils.h"
#include "Detector/DetectorFactory.h"
#include "Event/DataModel.h"
#include "Event/DetectorFrame.h"
#include "Script/Base/RawDataParser.h"
#include "Script/Base/ScriptFactory.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TDirectory.h"
#include "TVector3.h"
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <memory>

using namespace std;
using AnalysisUtils::FitTrack;
using AnalysisUtils::GetRange;
int clustermaxampthreshold=100;
int tracker_num=0;
//自己定义一个函数来做线性拟合，考虑误差权重
std::array<double,3> TrackFit(std::vector<double> z_vec1,std::vector<double> x_vec1,std::vector<double> x_error_vec1){//目前这个没考虑x方向的误差
  std::vector<double> weights(z_vec1.size());
  double weight_sum = 0.0;
  for (size_t i = 0; i < z_vec1.size(); ++i) {
    weights[i] = 1.0 / (x_error_vec1[i] * x_error_vec1[i]);
    weight_sum += weights[i];
  }
  double chi2=0.0;
  double x_mean=0.0;
  double z_mean=0.0;
  int n=x_vec1.size();
  for(int i=0;i<n;i++){
    x_mean+=x_vec1[i]*weights[i]/weight_sum;
    z_mean+=z_vec1[i]*weights[i]/weight_sum;
  }
  double lxx=0.0, lzz=0.0, lxz=0.0;
  for(int i=0;i<n;i++){
    lxx+=weights[i]*(x_vec1[i]-x_mean)*(x_vec1[i]-x_mean);
    lzz+=weights[i]*(z_vec1[i]-z_mean)*(z_vec1[i]-z_mean);
    lxz+=weights[i]*(x_vec1[i]-x_mean)*(z_vec1[i]-z_mean);
  }
  double slope = lxz / lzz;
  double intercept = x_mean - slope * z_mean;
  chi2 = 0.0;
  for (int i = 0; i < n; i++) {
    double x_fit = slope * z_vec1[i] + intercept;
    chi2 += weights[i] * (x_vec1[i] - x_fit) * (x_vec1[i] - x_fit);
  }
  return std::array<double,3>{chi2, slope, intercept};
}


void TrackAnalysisScript::LoadConfig(const json& config) {
    m_saveValidationData = config.value("saveValidationData", true);
    m_progressInterval = config.value("progressInterval", 10000);
}

void TrackAnalysisScript::Print() const {
    cout << "TrackAnalysisScript Configuration:" << endl;
    cout << "  Save Validation Data: " << (m_saveValidationData ? "Yes" : "No") << endl;
    cout << "  Progress Interval: " << m_progressInterval << endl;
}



bool TrackAnalysisScript::Execute() {
    using Clock = chrono::high_resolution_clock;
    auto t1 = Clock::now();

    cout << "\n========================================" << endl;
    cout << " Track Analysis Script" << endl;
    cout << "========================================" << endl;

    auto parser = GetParser();
    if (!parser) {
        cerr << "Error: Parser not set!" << endl;
        return false;
    }

    auto& factory = DetectorFactory::GetInstance();
    const auto trackers = factory.GetDetectorsByRole(Detector::Role::Tracker);

    // 初始化m_trackerIDs
    m_trackerIDs.clear();
    for (const auto& tracker : trackers) {
        m_trackerIDs.push_back(tracker->GetID());
        std::cout << "Found tracker: ID=" << tracker->GetID() << ", Name=" << tracker->GetName() << std::endl;
        tracker_num++;
    }
    sort(m_trackerIDs.begin(), m_trackerIDs.end());
    cout << "[Track Analysis] Initialized with " << m_trackerIDs.size() << " trackers" << endl;

    // 这个是用来读取每个通道的std的，注入到DetectorFrame里用的，可以用来卡阈值
    std::map<int, std::shared_ptr<const std::map<int, std::vector<float>>>> pedStdByDet;
    for (const auto& det : trackers) {
        const int detID = det->GetID();
        auto pedByType = std::make_shared<std::map<int, std::vector<float>>>();
        for (int type : det->getConfig().readoutPlaneType) {
            (*pedByType)[type] = parser->GetPedstd(detID, type);
        }
        pedStdByDet[detID] = pedByType;
    }

    const Long64_t total = parser->GetTotalEvents();//读取事件数
    std::vector<Event> events;
    events.reserve(total);

    // Event Loop
    cout << "\nProcessing events..." << endl;
    for (Long64_t i = 0; i < total; ++i) {
        if ((i % m_progressInterval) == 0 || i == total - 1) {
            cout << "\r[Track Analysis] Processed " << (i + 1) << "/" << total << flush;
        }

        auto rawHits = parser->LoadEvent(i);
        if (rawHits.empty())
            continue;

        Event evt{.eventID = int(i)};
        bool validEvent = true;

        for (auto& det : trackers) {//将原始数据变成strips，cluster，hit的形式保存在DetectorFrame里,而且要3个都有的事例
            auto detEvt = make_shared<DetectorFrame>(*det);
            detEvt->SetRawData(rawHits[det->GetID()]);

            // 把每条 strip 的 pedstd 注入到 DetectorFrame（不拷贝大向量）
            auto itPed = pedStdByDet.find(det->GetID());
            if (itPed != pedStdByDet.end()) {
                detEvt->SetStripPedestalStd(itPed->second);
            }

            if (!detEvt->Process()) {
                // std::cout << "\n[Track Analysis] Warning: Failed to process event " << i << " for detector " << det->GetID() << std::endl;
                validEvent = false;
                break;
            }
            const int id = det->GetID();
            // std::cout << "\n[Track Analysis] Event " << i << ", Detector " << id << ": " << detEvt->LocalHits().size() << " local hits" << std::endl;
            evt.detectorFramesMap[id] = std::move(detEvt);
        }

        if (validEvent) events.push_back(std::move(evt));
    }

    cout << "\n[Track Analysis] Valid events: " << events.size() << endl;

    if (events.empty()) {
        cerr << "Error: No valid events!" << endl;
        return false;
    }

    string trackFile = GetOutputDir() + "TrackInfo.root";
    TFile* f = new TFile(trackFile.c_str(), "RECREATE");

    // 直接调用私有方法
    RunTrackerAlign(events, f);

    f->cd();

    // Track Tree
    TTree* tTrack = new TTree("Tracks", "Track info");
    Int_t eventID;
    Track track;
    double t0;
    tTrack->Branch("eventID", &eventID);
    tTrack->Branch("track", &track);
    tTrack->Branch("t0", &t0);

    // Validation Tree
    TTree* tVal = nullptr;
    Int_t detID;
    Double_t resX, resY;
    Double_t hitX, hitY;
    vector<Int_t> clusterIndices;
    vector<StripHit> stripHits;
    vector<Cluster> clusters;

    if (m_saveValidationData) {
        tVal = new TTree("TrackerValidation", "Tracker QA");
        tVal->Branch("eventID", &eventID);
        tVal->Branch("detID", &detID);
        tVal->Branch("resX", &resX);
        tVal->Branch("resY", &resY);
        tVal->Branch("hitX", &hitX);
        tVal->Branch("hitY", &hitY);
        tVal->Branch("clusterIndices", &clusterIndices);
        tVal->Branch("stripHits", &stripHits);
        tVal->Branch("clusters", &clusters);
    }

    int saved = 0;
    int totalEvents = events.size();
    TH2D *hRes_det[tracker_num];
    TH1D *hRes_detx[tracker_num];
    TH1D *hRes_dety[tracker_num];
    for(int i=0;i<tracker_num;i++){
        hRes_det[i] = new TH2D(Form("hRes_det%d", i+1), Form("Residuals for Detector %d;ResX [mm];ResY [mm]", i+1), 50, -2, 2, 50, -2, 2);
        hRes_detx[i] = new TH1D(Form("hRes_det%dx", i+1), Form("Residual X for Detector %d;ResX [mm];Entries", i+1), 50, -2, 2);
        hRes_dety[i] = new TH1D(Form("hRes_det%dy", i+1), Form("Residual Y for Detector %d;ResY [mm];Entries", i+1), 50, -2, 2);
    }
    std::vector<LocalHit> localHits(trackers.size());
    cout << "\nSaving track data..." << endl;
    for (Event& evt : events) {
        vector<double> t0Vec;
        auto [bestTrack, hitIndices, success] = FindBestTrack(evt);

        if (!success)
            continue;
        for(int i=0;i<tracker_num;i++){
            auto det = trackers[i];
            int tid = det->GetID();
            const auto& detFrame = evt.detectorFramesMap[tid];
            localHits[i] = detFrame->LocalHits()[hitIndices[tid]];
        }
        std::vector<TVector3> predGlobalHits;
        for(int i=0;i<tracker_num;i++){//用其他的tracker来预言这个
            int tid = trackers[i]->GetID();
            predGlobalHits.clear();
            for(int j=0;j<tracker_num;j++){
                int tid2 = trackers[j]->GetID();
                if(tid==tid2) continue;
                auto detector = factory.GetDetector(tid2);
                predGlobalHits.push_back(detector->LocalToGlobal(localHits[j].localPos));
            }
            auto track1 = FitTrack(predGlobalHits);
            auto detector = trackers[i];
            resX = localHits[i].localPos.X() - detector->GlobalToLocal(detector->CalcHitFromTrack(track1)).X();
            resY = localHits[i].localPos.Y() - detector->GlobalToLocal(detector->CalcHitFromTrack(track1)).Y();
            hRes_det[i]->Fill(resX, resY);
            hRes_detx[i]->Fill(resX);
            hRes_dety[i]->Fill(resY);
        }


        if (m_saveValidationData) {
            for (const auto& det : trackers) {
                int tid = det->GetID();

                int hitIdx = hitIndices[tid];
                const auto& detFrame = evt.detectorFramesMap[tid];

                stripHits = detFrame->StripHits();
                clusters = detFrame->Clusters();
                const LocalHit& localHit = detFrame->LocalHits()[hitIdx];

                for (auto cluster : clusters)
                    t0Vec.push_back(cluster.time);

                auto detector = factory.GetDetector(tid);
                TVector3 predG = detector->CalcHitFromTrack(bestTrack);
                TVector3 predL = detector->GlobalToLocal(predG);

                hitX = localHit.localPos.X();
                hitY = localHit.localPos.Y();
                clusterIndices.clear();
                clusterIndices = localHit.clusterIndices;

                resX = hitX - predL.X();
                resY = hitY - predL.Y();

                tVal->Fill();
            }
        } else {
            for (const auto& det : trackers) {
                const auto& detFrame = evt.detectorFramesMap[det->GetID()];
                for (auto cluster : detFrame->Clusters()) {
                    t0Vec.push_back(cluster.time);
                }
            }
        }

        eventID = evt.eventID;
        track = bestTrack;
        // t0 = accumulate(t0Vec.begin(), t0Vec.end(), 0.0) / t0Vec.size();
        t0 = t0Vec[0];
        tTrack->Fill();
        saved++;
    }
    TF1 *fResX[tracker_num];
    TF1 *fResY[tracker_num];
    for(int i=0;i<tracker_num;i++){
        fResX[i] = new TF1(Form("fResX_det%d", i+1), "gaus", -2, 2);
        fResY[i] = new TF1(Form("fResY_det%d", i+1), "gaus", -2, 2);
        hRes_detx[i]->Fit(fResX[i], "Q");
        hRes_dety[i]->Fit(fResY[i], "Q");
    }
    TDirectory *dir = f->mkdir("Residuals");
    for(int i=0;i<tracker_num;i++){
        dir->WriteTObject(hRes_det[i]);
        dir->WriteTObject(hRes_detx[i]);
        dir->WriteTObject(hRes_dety[i]);
    }
    tTrack->Write();
    if (tVal) tVal->Write();
    f->Close();
    delete f;
    cout << "Dis sigmaX: ";
    for(int i=0;i<tracker_num;i++){
        cout << "trackid: " << m_trackerIDs[i] << " " << fResX[i]->GetParameter(2) << " ";
    }
    cout << "\n Dis sigmaY: ";
    for(int i=0;i<tracker_num;i++){
        cout << "trackid: " << m_trackerIDs[i] << " " << fResY[i]->GetParameter(2) << " ";
    }


    double z_vec[tracker_num];
    for(int i=0;i<tracker_num;i++){
        auto det = trackers[i];
        z_vec[i] = det->GetPos().Z();
    }
    std::vector<double> sigmachange(tracker_num);//这个是消除径迹误差的，假设每个探测器的本征分辨一样，看拟合径迹的时候，各位置sigma的变化率来估计本征分辨率,而且认为各sigma之间没有相关性
    std::vector<double> error;//这个是拟合径迹的时候用的误差
    for(int i=0;i<tracker_num;i++){
        error.push_back(1);
    }
    double error_matrix[tracker_num][tracker_num];
    for(int i=0;i<tracker_num;i++){
        error_matrix[i][i] = 1;
        double z = z_vec[i];
        std::vector<double> z_other(tracker_num-1);
        std::vector<double> weight_other(tracker_num-1);
        double weight_sum = 0;
        for(size_t j=0;j<tracker_num-1;j++){
            z_other[j] = z_vec[i>j ? j : j+1];
        }

        std::vector<double> x_other(tracker_num-1);
        std::vector<double> error_other(tracker_num-1);
        for(int j=0;j<tracker_num-1;j++){
            x_other[j] = 0;
            error_other[j] = 1;
        }
        auto [chi2, slope, intercept] = TrackFit(z_other, x_other, error_other);
        double baseline_track = intercept + slope * z;
        for(int j=0;j<tracker_num-1;j++){
            for(int k=0;k<tracker_num-1;k++){
                x_other[k] = 0;
            }
            x_other[j] = 1;
            auto [chi2, slope, intercept] = TrackFit(z_other, x_other, error_other);
            error_matrix[i][i>j?j:j+1] = slope*z+intercept-baseline_track;
            std::cout << error_matrix[i][i>j?j:j+1] << "i:j;"<<std::endl;
        }
        for(int j=0;j<tracker_num-1;j++){
            sigmachange[i] += error_matrix[i][j]*error_matrix[i][j];
        }
        sigmachange[i] = sqrt(sigmachange[i]);
    }
    cout << "\nError matrix: " << endl;
    for(int i=0;i<tracker_num;i++){
        for(int j=0;j<tracker_num;j++){
            cout << error_matrix[i][j] << " ";
        }        cout << endl;
    }
    cout << "\nSigma change factor: ";
    for(int i=0;i<tracker_num;i++){
        cout << "trackid: " << m_trackerIDs[i] << " " << sigmachange[i] << " ";
    }
    cout << "\nEstimated intrinsic resolutionx (mm): ";
    for(int i=0;i<tracker_num;i++){
        cout << "trackid: " << m_trackerIDs[i] << " " << fResX[i]->GetParameter(2)/sigmachange[i] << " ";
    }   
    cout << "\nEstimated intrinsic resolutiony (mm): ";
    for(int i=0;i<tracker_num;i++){
        cout << "trackid: " << m_trackerIDs[i] << " " << fResY[i]->GetParameter(2)/sigmachange[i] << " ";
    }



    cout << "[Track Analysis] Saved " << saved << " / " << totalEvents << " events" << endl;
    cout << "Output: " << trackFile << endl;

    auto t2 = Clock::now();
    double sec = chrono::duration<double>(t2 - t1).count();

    cout << "\n========================================" << endl;
    cout << "Track Analysis Complete" << endl;
    cout << "Time: " << fixed << setprecision(2) << sec << " seconds" << endl;
    cout << "========================================" << endl;

    return true;
}

// ========== 从TrackAnalysis迁移的方法实现 ==========

map<int, pair<double, double>> TrackAnalysisScript::ComputeTrackError(const vector<Event>& events, TFile* file) {

    auto& factory = DetectorFactory::GetInstance();
    map<int, vector<double>> residX, residY;

    // 计算残差,这个是用cluster中能量最高的hit来拟合轨迹的残差
    for (auto& e : events) {
        vector<GlobalHit> hits;
        for (int tid : m_trackerIDs) {
            auto detector = factory.GetDetector(tid);
            const LocalHit& hit = e.detectorFramesMap.at(tid)->LocalHits().at(0);
            hits.push_back(detector->LocalToGlobal(hit.localPos));
        }

        Track t = FitTrack(hits);

        for (size_t i = 0; i < m_trackerIDs.size(); ++i) {
            int tid = m_trackerIDs[i];
            auto detector = factory.GetDetector(tid);

            GlobalHit pred = detector->CalcHitFromTrack(t);
            GlobalHit meas = hits[i];

            residX[tid].push_back(meas.X() - pred.X());
            residY[tid].push_back(meas.Y() - pred.Y());
        }
    }

    // create directory
    int index = 0;
    std::string dirname = "TrackError";
    while (file->GetDirectory(dirname.c_str()) != nullptr) {
        ++index;
        dirname = "TrackError_" + std::to_string(index);
    }

    file->mkdir(dirname.c_str());
    file->cd(dirname.c_str());

    map<int, pair<double, double>> sigmas;

    for (int tid : m_trackerIDs) {
        auto& vx = residX[tid];
        auto& vy = residY[tid];

        auto [xmin, xmax] = GetRange(vx);
        TH1D* hx = new TH1D(Form("hFinalResX_%d", tid),
                            Form("Tracker %d Final Residual X;#DeltaX [mm];Events", tid),
                            200, xmin, xmax);
        for (double r : vx) hx->Fill(r);
        hx->Fit("gaus", "Q");
        hx->Write();
        sigmas[tid].first = hx->GetFunction("gaus")->GetParameter(2);

        auto [ymin, ymax] = GetRange(vy);
        TH1D* hy = new TH1D(Form("hFinalResY_%d", tid),
                            Form("Tracker %d Final Residual Y;#DeltaY [mm];Events", tid),
                            200, ymin, ymax);
        for (double r : vy) hy->Fill(r);
        hy->Fit("gaus", "Q");
        hy->Write();
        sigmas[tid].second = hy->GetFunction("gaus")->GetParameter(2);
    }

    return sigmas;
}

void TrackAnalysisScript::AlignTrackers(const vector<Event>& events) {
    auto& factory = DetectorFactory::GetInstance();

    const int nParPerDet = 3;  // dx, dy, dRotZ
    const UInt_t nPar = (m_trackerIDs.size() - 1) * nParPerDet;

    auto minim = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minim->SetTolerance(0.00005);
    minim->SetPrintLevel(0);

    // χ² 定义
    auto chi2 = [&](const double* par) {
        for (size_t i = 1; i < m_trackerIDs.size(); ++i) {
            int tid = m_trackerIDs[i];
            double dx = par[(i - 1) * nParPerDet];
            double dy = par[(i - 1) * nParPerDet + 1];
            double dRotZ = par[(i - 1) * nParPerDet + 2];
            auto detector = factory.GetDetector(tid);
            detector->SetAlignment(dx, dy, 0, 0, 0, dRotZ);
        }

        double chi2Sum = 0;
        int nevt = 0;

        for (auto& e : events) {
            vector<GlobalHit> hits;
            for (int tid : m_trackerIDs) {
                auto detector = factory.GetDetector(tid);
                hits.push_back(detector->LocalToGlobal(e.detectorFramesMap.at(tid)->LocalHits().at(0).localPos));

            }
            Track t = FitTrack(hits);
            chi2Sum += t.chi2;
            nevt++;
        }
        return chi2Sum > 0 ? chi2Sum / nevt : 1e9;
    };

    ROOT::Math::Functor f(chi2, nPar);
    minim->SetFunction(f);


    for(size_t i=1;i<m_trackerIDs.size();++i) {
        int tid = m_trackerIDs[i];
        minim->SetVariable((i - 1) * nParPerDet, Form("dx_%d", tid), factory.GetDetector(tid)->GetAlignPos().X(), 0.001);
        minim->SetVariable((i - 1) * nParPerDet + 1, Form("dy_%d", tid), factory.GetDetector(tid)->GetAlignPos().Y(), 0.001);
        minim->SetVariable((i - 1) * nParPerDet + 2, Form("dRotZ_%d", tid), factory.GetDetector(tid)->GetAlignRot().Z(), 0.001);
    }


    // for (UInt_t i = 0; i < nPar; i++)
    //     minim->SetVariable(i, Form("p%d", i), 0, 0.001);
    minim->Minimize();

    const double* par = minim->X();
    for (size_t i = 1; i < m_trackerIDs.size(); ++i) {
        int tid = m_trackerIDs[i];
        double dx = par[(i - 1) * 3];
        double dy = par[(i - 1) * 3 + 1];
        double dRotZ = par[(i - 1) * 3 + 2];
        auto detector = factory.GetDetector(tid);
        detector->SetAlignment(dx, dy, 0, 0, 0, dRotZ);
    }
    delete minim;
}


void TrackAnalysisScript::AlignTrackersfurther(const vector<Event>& events) {
    auto& factory = DetectorFactory::GetInstance();

    const int nParPerDet = 6;  // dx, dy, dz, dRotX, dRotY, dRotZ
    const UInt_t nPar = (m_trackerIDs.size() - 1) * nParPerDet;

    auto minim = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minim->SetTolerance(500);
    minim->SetPrintLevel(0);

    // χ² 定义
    auto chi2 = [&](const double* par) {
        for (size_t i = 1; i < m_trackerIDs.size(); ++i) {
            int tid = m_trackerIDs[i];
            double dx = par[(i - 1) * nParPerDet];
            double dy = par[(i - 1) * nParPerDet + 1];
            double dz = par[(i - 1) * nParPerDet + 2];
            double dRotX = par[(i - 1) * nParPerDet + 3];
            double dRotY = par[(i - 1) * nParPerDet + 4];
            double dRotZ = par[(i - 1) * nParPerDet + 5];
            auto detector = factory.GetDetector(tid);
            detector->SetAlignment(dx, dy, dz, dRotX, dRotY, dRotZ);
        }

        double chi2Sum = 0;
        int nevt = 0;

        for (auto& e : events) {
            vector<GlobalHit> hits;
            for (int tid : m_trackerIDs) {
                auto detector = factory.GetDetector(tid);
                hits.push_back(detector->LocalToGlobal(e.detectorFramesMap.at(tid)->LocalHits().at(0).localPos));

            }
            Track t = FitTrack(hits);
            chi2Sum += t.chi2;
            nevt++;
        }
        return chi2Sum > 0 ? chi2Sum / nevt : 1e9;
    };

    ROOT::Math::Functor f(chi2, nPar);
    minim->SetFunction(f);


    for(size_t i=1;i<m_trackerIDs.size();++i) {
        int tid = m_trackerIDs[i];
        minim->SetVariable((i - 1) * nParPerDet, Form("dx_%d", tid), factory.GetDetector(tid)->GetAlignPos().X(), 0.001);
        minim->SetVariable((i - 1) * nParPerDet + 1, Form("dy_%d", tid), factory.GetDetector(tid)->GetAlignPos().Y(), 0.001);
        minim->SetVariable((i - 1) * nParPerDet + 2, Form("dz_%d", tid), factory.GetDetector(tid)->GetAlignPos().Z(), 0.001);
        minim->SetVariable((i - 1) * nParPerDet + 3, Form("dRotX_%d", tid), factory.GetDetector(tid)->GetAlignRot().X(), 0.001);
        minim->SetVariable((i - 1) * nParPerDet + 4, Form("dRotY_%d", tid), factory.GetDetector(tid)->GetAlignRot().Y(), 0.0001);
        minim->SetVariable((i - 1) * nParPerDet + 5, Form("dRotZ_%d", tid), factory.GetDetector(tid)->GetAlignRot().Z(), 0.0001);
        minim->SetVariableLimits((i - 1) * nParPerDet, -20, 20);       // dx
        minim->SetVariableLimits((i - 1) * nParPerDet + 1, -20, 20);   // dy
        minim->SetVariableLimits((i - 1) * nParPerDet + 2, -10, 10);   // dz
        minim->SetVariableLimits((i - 1) * nParPerDet + 3, -0.15, 0.15); // dRotX
        minim->SetVariableLimits((i - 1) * nParPerDet + 4, -0.15, 0.15); // dRotY
        minim->SetVariableLimits((i - 1) * nParPerDet + 5, -0.15, 0.15); // dRotZ
    }


    // for (UInt_t i = 0; i < nPar; i++)
    //     minim->SetVariable(i, Form("p%d", i), 0, 0.001);
    minim->Minimize();

    const double* par = minim->X();
    for (size_t i = 1; i < m_trackerIDs.size(); ++i) {
        int tid = m_trackerIDs[i];
        double dx = par[(i - 1) * 6];
        double dy = par[(i - 1) * 6 + 1];
        double dz = par[(i - 1) * 6 + 2];
        double dRotX = par[(i - 1) * 6 + 3];
        double dRotY = par[(i - 1) * 6 + 4];
        double dRotZ = par[(i - 1) * 6 + 5];
        auto detector = factory.GetDetector(tid);
        detector->SetAlignment(dx, dy, dz, dRotX, dRotY, dRotZ);
    }
    delete minim;
}

pair<double, double> TrackAnalysisScript::ComputePredictionError(int targetDetID) {
    auto& factory = DetectorFactory::GetInstance();

    // Seed tracker IDs
    const int seed1 = m_seedTrackerIDs[0];
    const int seed2 = m_seedTrackerIDs[1];

    // Their intrinsic resolutions (assumed independent)
    const auto [s1x, s1y] = m_sigmaMap[seed1];
    const auto [s2x, s2y] = m_sigmaMap[seed2];

    // Target intrinsic resolution
    const auto [stx, sty] = m_sigmaMap[targetDetID];

    // Z positions
    auto det1 = factory.GetDetector(seed1);
    auto det2 = factory.GetDetector(seed2);
    auto detT = factory.GetDetector(targetDetID);

    const double z1 = det1->GetPos().Z();
    const double z2 = det2->GetPos().Z();
    const double zt = detT->GetPos().Z();

    // Linear interpolation/extrapolation factor α
    const double L = (z2 - z1);
    if (std::abs(L) < 1e-12) {
        std::cerr << "[ComputePredictionError] Seed trackers have identical z!" << std::endl;
        return {0, 0};
    }
    const double alpha = (zt - z1) / L;

    // === Exact track prediction variance (no approximation) ===
    const double varTrackX = (1 - alpha) * (1 - alpha) * s1x * s1x + alpha * alpha * s2x * s2x;
    const double varTrackY = (1 - alpha) * (1 - alpha) * s1y * s1y + alpha * alpha * s2y * s2y;

    // === Combine with target detector intrinsic resolution ===
    const double sigmaPredX = std::sqrt(varTrackX + stx * stx);
    const double sigmaPredY = std::sqrt(varTrackY + sty * sty);

    return {sigmaPredX, sigmaPredY};
}

tuple<Track, map<int, int>, bool> TrackAnalysisScript::FindBestTrack(const Event& event) {
    auto& factory = DetectorFactory::GetInstance();

    // 选择seed tracker, hit数量最小
    vector<pair<int, int>> hitCounts;  // (hitCount, trackerIndex)
    for (size_t i = 0; i < m_trackerIDs.size(); ++i) {
        int tid = m_trackerIDs[i];
        int count = event.detectorFramesMap.at(tid)->LocalHits().size();
        hitCounts.push_back({count, i});
    }

    sort(hitCounts.begin(), hitCounts.end());

    // 选择击中数最少的两个tracker
    int seedTrackerId1 = m_trackerIDs[hitCounts[0].second];
    int seedTrackerId2 = m_trackerIDs[hitCounts[1].second];

    m_seedTrackerIDs = {seedTrackerId1, seedTrackerId2};

    const auto& hits_seed1 = event.detectorFramesMap.at(seedTrackerId1)->LocalHits();
    const auto& hits_seed2 = event.detectorFramesMap.at(seedTrackerId2)->LocalHits();

    auto det1 = factory.GetDetector(seedTrackerId1);
    auto det2 = factory.GetDetector(seedTrackerId2);
    std::vector<Track> candidateTracks;
    std::vector<map<int, int>> candidateHitIndices;

    // 贪婪搜索：遍历seed组合
    for (size_t i1 = 0; i1 < hits_seed1.size(); ++i1) {
        for (size_t i2 = 0; i2 < hits_seed2.size(); ++i2) {
            if(hits_seed1[i1].isValid == false || hits_seed2[i2].isValid == false) continue; // 如果seed hit无效，跳过这个组合
            GlobalHit globalHit1 = det1->LocalToGlobal(hits_seed1[i1].localPos);
            GlobalHit globalHit2 = det2->LocalToGlobal(hits_seed2[i2].localPos);

            Track currentTrack = FitTrack({globalHit1, globalHit2});

            // 初始化击中索引数组
            std::map<int, int> hitIndices;
            hitIndices[seedTrackerId1] = i1;
            hitIndices[seedTrackerId2] = i2;

            bool allValid = true;

            for (size_t i = 0; i < m_trackerIDs.size(); ++i) {

                int tid = m_trackerIDs[i];
                if (tid == seedTrackerId1 || tid == seedTrackerId2) continue;

                auto detector = factory.GetDetector(tid);

                // 预测击中位置
                GlobalHit predictedGlobal = detector->CalcHitFromTrack(currentTrack);

                // 计算预测误差
                auto [sigma_pred_X, sigma_pred_Y] = ComputePredictionError(tid);

                // 寻找最近击中
                const auto& hits_i = event.detectorFramesMap.at(tid)->LocalHits();
                double minDist = numeric_limits<double>::infinity();
                int bestIdx = -1;

                for (size_t j = 0; j < hits_i.size(); ++j) {
                    if(hits_i[j].isValid == false) continue; // 如果hit无效，跳过这个hit
                    GlobalHit globalHit = detector->LocalToGlobal(hits_i[j].localPos);
                    double resX = globalHit.X() - predictedGlobal.X();
                    double resY = globalHit.Y() - predictedGlobal.Y();

                    double normDist = sqrt((resX / sigma_pred_X) * (resX / sigma_pred_X) + (resY / sigma_pred_Y) * (resY / sigma_pred_Y));

                    if (normDist < minDist) {
                        minDist = normDist;
                        bestIdx = j;
                    }
                }

                if (minDist > 10.0) {
                    allValid = false;
                    break;
                }

                hitIndices[tid] = bestIdx;
            }

            if (allValid) {

                vector<GlobalHit> allGlobalHits;
                for (size_t i = 0; i < m_trackerIDs.size(); ++i) {
                    int tid = m_trackerIDs[i];
                    int hitIdx = hitIndices[tid];
                    LocalHit localHit = event.detectorFramesMap.at(tid)->LocalHits().at(hitIdx);
                    auto detector = factory.GetDetector(tid);
                    GlobalHit globalHit = detector->LocalToGlobal(localHit.localPos);
                    allGlobalHits.push_back(globalHit);
                }

                Track finalTrack = FitTrack(allGlobalHits);
                candidateTracks.push_back(finalTrack);
                candidateHitIndices.push_back(hitIndices);
                //return {finalTrack, hitIndices, true};
            }
        }
    }
    if(candidateTracks.empty()) {
        return {Track{}, std::map<int, int>(), false};
    }
    else{
        for(int i=0;i<candidateTracks.size();++i){
            candidateTracks[i].chi2 /= (m_trackerIDs.size() - 2);
        }
    }
    auto minIt = min_element(candidateTracks.begin(), candidateTracks.end(), [](const Track& a, const Track& b) {
        return a.chi2 < b.chi2;
    });
    int minIndex = minIt - candidateTracks.begin();
    return {candidateTracks[minIndex], candidateHitIndices[minIndex], true};
}




void TrackAnalysisScript::RunTrackerAlign(const vector<Event>& events, TFile* file) {
    auto& factory = DetectorFactory::GetInstance();
    const auto& trackers = factory.GetDetectorsByRole(Detector::Role::Tracker);

    cout << "[TrackerAlign] Perform tracker alignment? (y/n): ";

    char choice;
    cin >> choice;
    cin.ignore();

    if (choice != 'y' && choice != 'Y') return;

    // 1) Collect single-hit events
    vector<Event> single;
    single.reserve(events.size());

    for (const auto& evt : events) {//筛选所有的tracker在两个方向上都只有一个cluster且maxAmp>75的事件，且这个cluster不能包含坏strip
        bool ok = true;
        for (auto& det : trackers) {//所有的tracker
            int id = det->GetID();
            for (int type : det->getConfig().readoutPlaneType) {//两个方向
                int numberofoverthereshold = 0;
                bool isbadstrip=false;
                for (auto cluster : evt.detectorFramesMap.at(id)->Clusters(type)) {
                    if (cluster.maxAmp > 75) {
                        numberofoverthereshold++;
                        if(cluster.isBad){
                            isbadstrip=true;
                            break;
                        }
                    }
                }
                if(numberofoverthereshold != 1 || isbadstrip){
                    ok = false;
                    break;
                }
            }
        }
        if (ok) single.push_back(evt);
    }
    cout << "[TrackerAlign] Found " << single.size() << " single-hit events." << endl;
    if (single.size() < 20) {
        cerr << "[Alignment] WARNING: Too few single-hit events: " << single.size() << endl;
    }

    // 2) Coarse alignment (residual mean shifts)
    int refIndex = 0;
    int refID = m_trackerIDs[refIndex];
    TH2D *hitmap_align[tracker_num];
    for(int i=0;i<tracker_num;i++){
        hitmap_align[i] = new TH2D(Form("hitmap_align_det%d", i+1), Form("Hitmap for Alignment - Detector %d;X [mm];Y [mm]", i+1), 100, 0, 100, 100, 0,100);
    }
    map<int, vector<double>> resX, resY;

    for (const auto& evt : single) {
        auto refDet = factory.GetDetector(refID);
        auto ref = refDet->LocalToGlobal(evt.detectorFramesMap.at(refID)->LocalHits().at(0).localPos);
        hitmap_align[refIndex]->Fill(ref.X(), ref.Y());
        for (size_t i = 0; i < m_trackerIDs.size(); ++i) {
            int tid = m_trackerIDs[i];
            if(tid == refID) continue;
            auto detector = factory.GetDetector(tid);
            auto hit = detector->LocalToGlobal(evt.detectorFramesMap.at(tid)->LocalHits().at(0).localPos);
            resX[tid].push_back(hit.X() - ref.X());
            resY[tid].push_back(hit.Y() - ref.Y());
            // std::cout << hit.X() << " " << hit.Y() << std::endl;
            hitmap_align[i]->Fill(hit.X(), hit.Y());
        }
    }
    TDirectory *dirhitmapalign = file->mkdir("HitmapsForAlignment");
    for(int i=0;i<tracker_num;i++){
        dirhitmapalign->WriteTObject(hitmap_align[i]);
    }
    TDirectory *dircoarse = file->mkdir("CoarseAlignment");

    for (size_t i = 0; i < m_trackerIDs.size(); ++i) {
        int tid = m_trackerIDs[i];
        if(tid == refID) continue;
        // X
        auto [xmin, xmax] = GetRange(resX[tid]);
        TH1D* hX = new TH1D(Form("hResX_%d", tid),
                            Form("Tracker %d Residual X;#DeltaX [mm];Events", tid),
                            200, xmin, xmax);
        for (double r : resX[tid]) hX->Fill(r);
        hX->Fit("gaus", "Q");
        dircoarse->WriteTObject(hX);

        double dx = hX->GetFunction("gaus")->GetParameter(1);

        // Y
        auto [ymin, ymax] = GetRange(resY[tid]);
        TH1D* hY = new TH1D(Form("hResY_%d", tid),
                            Form("Tracker %d Residual Y;#DeltaY [mm];Events", tid),
                            200, ymin, ymax);
        for (double r : resY[tid]) hY->Fill(r);
        hY->Fit("gaus", "Q");
        dircoarse->WriteTObject(hY);

        double dy = hY->GetFunction("gaus")->GetParameter(1);

        // apply coarse shifts
        auto detector = factory.GetDetector(tid);
        detector->SetAlignment(-dx, -dy, 0, 0, 0, 0);//注意这里是反向补偿
        std::cout << "[TrackerAlign] Coarse alignment applied to Tracker " << tid
                  << ": dx = " << dx << " mm, dy = " << dy << " mm" << std::endl;
    }
    
    AlignTrackers(single);
    AlignTrackersfurther(single);
    // 3) first σ measurement
    auto sigmaMap = ComputeTrackError(single, file);

    // // 4) filter events by 3σ
    vector<Event> filtered = single;
    filtered.reserve(single.size());

    for (auto& e : single) {
        bool good = true;

        vector<GlobalHit> hits;
        for (int tid : m_trackerIDs) {
            auto detector = factory.GetDetector(tid);
            hits.push_back(detector->LocalToGlobal(e.detectorFramesMap.at(tid)->LocalHits().at(0).localPos));
        }

        Track t = FitTrack(hits);

        for (size_t i = 0; i < m_trackerIDs.size(); ++i) {
            int tid = m_trackerIDs[i];
            auto detector = factory.GetDetector(tid);
            auto pred = detector->CalcHitFromTrack(t);
            auto meas = hits[i];

            double dx = meas.X() - pred.X();
            double dy = meas.Y() - pred.Y();

            if (fabs(dx) > 3 * sigmaMap[tid].first || fabs(dy) > 3 * sigmaMap[tid].second) {
                good = false;
                break;
            }
        }
        if (good) filtered.push_back(e);
    }

    // 5) second fine alignment
    AlignTrackers(filtered);
    
    // 6) final σ measurement
    m_sigmaMap = ComputeTrackError(filtered, file);
    std::cout << "\n[TrackerAlign] Final intrinsic resolutions after alignment:" << std::endl;
    for (int tid : m_trackerIDs) {
        std::cout << "  Tracker " << tid << ": σx = " << m_sigmaMap[tid].first << " mm, σy = " << m_sigmaMap[tid].second << " mm" << std::endl;
    }

    cout << "\n[Tracker Alignment] Final detector positions:" << endl;
    for (int tid : m_trackerIDs) {
        auto det = factory.GetDetector(tid);
        TVector3 pos = det->GetPos();
        TVector3 rot = det->GetRot();
        cout << "  " << det->GetName() << " (ID=" << tid << "):" << endl;
        cout << "    \"position\": [" << fixed << setprecision(5)
             << pos.X() << ", " << pos.Y() << ", " << pos.Z() << "]," << endl;
        cout << "    \"rotation\": ["
             << rot.X() << ", " << rot.Y() << ", " << rot.Z() << "]," << endl;
    }
    cout << "\n[Tracker Alignment] done." << endl;
}

REGISTER_SCRIPT("TrackAnalysis", TrackAnalysisScript);
