// ============================================================
// microTPCsimulation.C — 基于抽样库的可配置角度 Track 信号模拟
//
//   流程：
//     1. TrackHeed 生成 20° 入射角的 muon track
//     2. 对 track 中每个 cluster 的每个初级电子，漂移到读出平面
//     3. 按漂移终点 (ex1, ey1) 在抽样库中查表获取信号
//     4. 聚合所有电子信号，输出 stripid + 波形 + track 斜率/截距
//
//   依赖：需先运行 analysis/sampleconstruction.C 中的 buildSampleLibrary()
// ============================================================

#include <Rtypes.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TTree.h>
#include <TRandom3.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

#include "Garfield/ComponentComsol.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/RandomEngine.hh"
#include "Garfield/RandomEngineRoot.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/Shaper.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/TrackHeed.hh"
#include "Config.h"

using namespace Garfield;

// ---------- 抽样库元数据（从 samplelibrary.root 的 meta tree 读取） ----------
static double gLib_xmin, gLib_xmax, gLib_ymin, gLib_ymax;
static int    gLib_nbins;
static std::vector<double> gLib_t;  // 共享时间轴（从 meta 读取一次）
static const int kNbins = 10;       // 与 buildSampleLibrary.C 中 kNbins 一致

// ---- 信号通道名称表（8 个，与 microTPC.C 一致）----
static const int N_SIG = 8;
static const char* SIG_NAME[N_SIG] = {
        "sigxtestminus1","sigxtest1", 
    "sigxtestminus2",  "sigzero","sigxtest2",   
     "sigytestminus1", "sigy","sigytest1"
};//前2个信号是x方向落到两个条中间时相邻2个条产生的，中间3个信号是x方向落到条上时相邻3个条产生的，后3个信号是y方向上相邻3个条产生的

// ============================================================
// main — 主模拟程序
// ============================================================
int main(int argc, char *argv[])
{
    // ---------- 参数 ----------
    Config config;

    if (argc >= 2) {
        config = ReadConfig(argv[1]);
    } else {
        config = ReadConfig("../config/default.jsonc");
    }

    const std::string outputFile = config.outputDir + "/events.root";

    constexpr double pitch = 0.04;
    const double trackAngleRad = config.simulation.theta * M_PI / 180.0;
    const double trackPhiRad = config.simulation.phi * M_PI / 180.0;
    constexpr int N_STRIPS = 20;                 // x/y 方向各 20 条
    constexpr double zBottom = 0.01;             // 漂移区下边界，单位 cm
    const double z0 = config.simulation.gap;     // 漂移区上边界，单位 cm
    const double zMid = 0.5 * z0;

    if (config.simulation.nEvents == 0 ||
        config.simulation.sampleInterval <= 0. ||
        config.electronics.sampleInterval <= 0. || config.totalTime <= 0. ||
        config.electronics.n == 0 || config.electronics.tau <= 0. ||
        z0 <= zBottom) {
        std::cerr << "Error: invalid simulation/electronics configuration." << std::endl;
        return 1;
    }

    //设置随机数
    TRandom3 rnd(config.simulation.seed);
    RandomEngineRoot garfieldRandomEngine;
    garfieldRandomEngine.SetSeed(config.simulation.seed);
    Random::SetEngine(garfieldRandomEngine);

    std::cout << "Simulation configuration: theta=" << config.simulation.theta
              << " deg, phi=" << config.simulation.phi
              << " deg, events=" << config.simulation.nEvents
              << ", seed=" << config.simulation.seed
              << ", position-mode=" << config.simulation.mode << std::endl;

    // ---------- 加载抽样库 ----------
    const char *libFile = "../analysis/samplelibrarysmall.root";
    TFile *fLib = TFile::Open(libFile, "READ");
    if (!fLib || fLib->IsZombie()) {
        std::cerr << "Error: Cannot open " << libFile << std::endl;
        return 1;
    }

    // —— 读取分箱元数据 + 共享时间轴 ——
    TTree *meta = (TTree *)fLib->Get("meta");
    if (!meta) {
        std::cerr << "Error: No meta tree found in " << libFile << std::endl;
        fLib->Close();
        return 1;
    }
    {
        double xmin, xmax, ymin, ymax;
        int    nbins;
        std::vector<double> *meta_t = nullptr;
        meta->SetBranchAddress("xmin",  &xmin);
        meta->SetBranchAddress("xmax",  &xmax);
        meta->SetBranchAddress("ymin",  &ymin);
        meta->SetBranchAddress("ymax",  &ymax);
        meta->SetBranchAddress("nbins", &nbins);
        meta->SetBranchAddress("t",     &meta_t);
        meta->GetEntry(0);
        gLib_xmin  = xmin;
        gLib_xmax  = xmax;
        gLib_ymin  = ymin;
        gLib_ymax  = ymax;
        gLib_nbins = nbins;
        gLib_t     = *meta_t;

        if (gLib_nbins != kNbins) {
            std::cerr << "Error: sample library nbins=" << gLib_nbins
                      << " does not match compiled kNbins=" << kNbins << std::endl;
            fLib->Close();
            return 1;
        }

        std::cout << "[loadSampleLibrary] Binning: nbins=" << nbins
                  << "  x=[" << xmin << ", " << xmax << "]"
                  << "  y=[" << ymin << ", " << ymax << "]"
                  << "  t.size=" << gLib_t.size() << std::endl;
    }

    // —— 逐一打开各区域的独立树，预绑定分支 ——
    TTree *regionTrees[kNbins][kNbins] = {};
    std::vector<double> *region_sig[kNbins][kNbins][N_SIG] = {};
    Long64_t region_n[kNbins][kNbins] = {};
    {
        int nNonEmpty = 0;
        for (int ix = 0; ix < gLib_nbins; ++ix) {
            for (int iy = 0; iy < gLib_nbins; ++iy) {
                std::string tname = "region_" + std::to_string(ix) + "_" + std::to_string(iy);
                TTree *rt = (TTree *)fLib->Get(tname.c_str());
                if (!rt) continue;

                regionTrees[ix][iy] = rt;
                region_n[ix][iy]   = rt->GetEntries();
                for (int is = 0; is < N_SIG; ++is)
                    rt->SetBranchAddress(SIG_NAME[is], &region_sig[ix][iy][is]);

                ++nNonEmpty;
                if (nNonEmpty <= 10)
                    std::cout << "  region (" << ix << "," << iy << "): "
                              << region_n[ix][iy] << " samples" << std::endl;
            }
        }
        std::cout << "  Total non-empty regions: " << nNonEmpty
                  << " / " << gLib_nbins * gLib_nbins << std::endl;
    }
    // NOTE: fLib 不关闭，各 region tree 需要它保持打开

    // ---------- 气体 ----------
    MediumMagboltz gas;
    gas.LoadGasFile(
        config.simulation.detectorDir + "/" + config.simulation.gasFile);
    constexpr double rPenning = 0.20;
    constexpr double lambdaPenning = 0.;
    gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");
    gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");
    gas.SetMaxElectronEnergy(200.);

    // ---------- 电场 ----------
    ComponentComsol fm;
    fm.Initialise(
        (config.simulation.detectorDir + "/mesh.mphtxt").c_str(),
        (config.simulation.detectorDir + "/mplist.txt").c_str(),
        (config.simulation.detectorDir + "/" +
         config.simulation.fieldFile).c_str(),
        "mm");
    fm.EnablePeriodicityX();
    fm.EnablePeriodicityY();
    fm.SetGas(&gas);

    // 在漂移区中心取局部电场，用 Garfield 计算z方向电子漂移速度。
    double ex = 0., ey = 0., ez = 0.;
    double vx = 0., vy = 0., vz = 0.;
    int fieldStatus = 0;
    Medium* medium = nullptr;
    fm.ElectricField(0.5 * N_STRIPS * pitch, 0.5 * N_STRIPS * pitch, zMid,
                     ex, ey, ez, medium, fieldStatus);
    if (fieldStatus != 0 || medium != &gas ||
        !gas.ElectronVelocity(ex, ey, ez, 0., 0., 0., vx, vy, vz)) {
        std::cerr << "Error: cannot calculate the drift velocity at the "
                  << "center of the drift region." << std::endl;
        return 1;
    }

    const double dz =
        std::abs(vz) * config.electronics.sampleInterval;
    if (dz <= 0.) {
        std::cerr << "Error: the z drift velocity is zero." << std::endl;
        return 1;
    }

    // 先从 z=0 开始按 dz 分层，再用 zBottom 和 gap 截断。
    // 因此最下层可以是不完整的，模板位置取该层自身的中点。
    std::vector<double> tierCenter;
    for (unsigned int tier = 0; tier * dz < z0; ++tier) {
        const double lower = std::max(zBottom, tier * dz);
        const double upper = std::min((tier + 1) * dz, z0);
        if (upper > lower) {
            tierCenter.push_back(0.5 * (lower + upper));
        }
    }

    std::cout << "Drift velocity z: " << std::abs(vz) << " cm/ns\n"
              << "Tier thickness: " << dz << " cm\n"
              << "Number of tiers: " << tierCenter.size() << std::endl;

    // ---------- Sensor ----------
    Sensor sensor(&fm);
    sensor.SetArea(0, 0, 0.01, N_STRIPS * pitch, N_STRIPS * pitch, 0.51);

    const unsigned int nTimeBins = static_cast<unsigned int>(
        std::ceil(config.totalTime / config.simulation.sampleInterval));
    const double tmin = 0.;
    const double tstep = config.simulation.sampleInterval;
    sensor.SetTimeWindow(tmin, tstep, nTimeBins);
    Shaper shaper(config.electronics.n, config.electronics.tau, 1., "unipolar");
    // 抽样库保存的是卷积前信号；所有电子叠加后在此处统一成形。
    // 直接调用 Garfield 的 Shaper，避免手写核与 Shaper 参数发生漂移。
    std::vector<double> kernel_vec(nTimeBins);
    for (int i = 0; i < nTimeBins; ++i) {
        double t = tmin + i * tstep;
        kernel_vec[i] = shaper.Shape(t);
    }
 
    // ---------- TrackHeed ----------
    TrackHeed track(&sensor);
    track.SetParticle("muon");
    track.SetEnergy(150.e9);  // 150 GeV muon
    if (config.simulation.deltaElectron) {
        track.EnableDeltaElectronTransport();
    } else {
        track.DisableDeltaElectronTransport();
    }

    // ---------- Drift ----------
    AvalancheMC drift(&sensor);
    drift.SetTimeWindow(0, 1400);
    drift.SetDistanceSteps(1.e-3);

    // ---------- 输出 Tree ----------
    std::filesystem::create_directories(config.outputDir);
    TFile *fout = new TFile(outputFile.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << outputFile << std::endl;
        return 1;
    }

    std::vector<double> vec_t;
    for(int i = 0; i < nTimeBins; ++i) vec_t.push_back(tmin + i * tstep);

    int event;
    // 内部信号累加器（固定 N_STRIPS 条 x + N_STRIPS 条 y）
    std::vector<double> vec_sigx[N_STRIPS];
    std::vector<double> vec_sigy[N_STRIPS];
    std::vector<double> vec_sigx_conv[N_STRIPS];
    std::vector<double> vec_sigy_conv[N_STRIPS];
    
    // 每 event 的输出：strip ID 向量 + 波形向量 + track 斜率/截距
    std::vector<int>                  stripIDs;
    std::vector<std::vector<double>>  waveforms;
    std::vector<std::vector<double>>  waveforms_conv;
    std::vector<int>                  types;  // 0: x strip, 1: y strip
    double kx, bx;
    double ky, by;
    double trackX0 = 0.0, trackY0 = 0.0;
    double trackMidX = 0.0, trackMidY = 0.0;

    TTree *tout = new TTree("tree", "Track Sampled data");
    tout->Branch("event",         &event,          "event/I");
    tout->Branch("t",             &vec_t);
    tout->Branch("types",         &types);
    tout->Branch("stripIDs",     &stripIDs);
    tout->Branch("waveforms",     &waveforms);
    tout->Branch("waveforms_conv",&waveforms_conv);
    tout->Branch("kx",      &kx,       "kx/D");
    tout->Branch("bx",  &bx,   "bx/D");
    tout->Branch("ky",      &ky,       "ky/D");
    tout->Branch("by",  &by,   "by/D");
    tout->Branch("track_x0", &trackX0, "track_x0/D");
    tout->Branch("track_y0", &trackY0, "track_y0/D");
    tout->Branch("track_mid_x", &trackMidX, "track_mid_x/D");
    tout->Branch("track_mid_y", &trackMidY, "track_mid_y/D");

    // 每个 event 的电子终点位置（tree2，与 tree 一一对应）
    TTree *tpos = new TTree("tree2", "Electron endpoint positions per event");
    std::vector<double> vec_xe1, vec_ye1, vec_ze1, vec_te1;
    std::vector<int>    vec_xe_stripid, vec_ye_stripid;
    tpos->Branch("event",        &event,           "event/I");
    tpos->Branch("xe1",          &vec_xe1);
    tpos->Branch("ye1",          &vec_ye1);
    tpos->Branch("ze1",          &vec_ze1);
    tpos->Branch("te1",          &vec_te1);
    tpos->Branch("xe_stripid",   &vec_xe_stripid);
    tpos->Branch("ye_stripid",   &vec_ye_stripid);

    // 每个 event 的电子初始电离位置（tree3，与 tree 一一对应）
    TTree *tinit = new TTree("tree3", "Electron initial ionization positions per event");
    std::vector<double> vec_xe0, vec_ye0, vec_ze0, vec_te0, vec_ee0;
    std::vector<int>    vec_xe0_stripid, vec_ye0_stripid;
    tinit->Branch("event",        &event,           "event/I");
    tinit->Branch("xe0",          &vec_xe0);
    tinit->Branch("ye0",          &vec_ye0);
    tinit->Branch("ze0",          &vec_ze0);
    tinit->Branch("te0",          &vec_te0);
    tinit->Branch("ee0",          &vec_ee0);
    tinit->Branch("xe0_stripid",  &vec_xe0_stripid);
    tinit->Branch("ye0_stripid",  &vec_ye0_stripid);

    // 统计漂移终点到抽样库的接受率，便于检查 (ix, iy) 是否越界。
    Long64_t nPrimary = 0;
    Long64_t nDriftSuccess = 0;
    Long64_t nReadoutInRange = 0;
    Long64_t nLibraryOutOfRange = 0;
    Long64_t nEmptyRegion = 0;
    Long64_t nSampledElectrons = 0;
    // 统计总的Garfield 漂移、ROOT 模板读取/叠加和电子学卷积的时间
    double drifttime = 0.;
    double readtime = 0.;
    double convolutiontime = 0.;
    const auto eventLoopStart = std::chrono::steady_clock::now();

    // ==================== 事件循环 ====================

    for (unsigned int i = 0; i < config.simulation.nEvents; ++i) {

        event = i + config.simulation.runid * config.simulation.nEvents;
        // —— 信号清零 ——
        for (int indexx = 0; indexx < N_STRIPS; ++indexx) {
            vec_sigx[indexx].assign(nTimeBins, 0);
            vec_sigy[indexx].assign(nTimeBins, 0);
            vec_sigx_conv[indexx].assign(nTimeBins, 0);
            vec_sigy_conv[indexx].assign(nTimeBins, 0);
        }

        if(i % 100 == 0) {
            std::cout << "Event: " << i << "/"
                      << config.simulation.nEvents << std::endl;
        }
        sensor.ClearSignal();

        // —— 清空本 event 的电子信息向量 ——
        vec_xe1.clear(); vec_ye1.clear(); vec_ze1.clear(); vec_te1.clear();
        vec_xe_stripid.clear(); vec_ye_stripid.clear();
        vec_xe0.clear(); vec_ye0.clear(); vec_ze0.clear(); vec_te0.clear(); vec_ee0.clear();
        vec_xe0_stripid.clear(); vec_ye0_stripid.clear();

        // 固定轨迹在漂移区中点附近穿过读出中心，避免角度扫描受到边界距离影响。
        // — 计算轨迹中点位置 ——
        trackMidX = 0.5 * N_STRIPS * pitch;
        trackMidY = 0.5 * N_STRIPS * pitch;
        if (config.simulation.mode == "uniform") {
            trackMidX += rnd.Uniform(-0.5 * pitch, 0.5 * pitch);
            trackMidY += rnd.Uniform(-0.5 * pitch, 0.5 * pitch);
        }

        const double dx = std::sin(trackAngleRad) * std::cos(trackPhiRad);
        const double dy = std::sin(trackAngleRad) * std::sin(trackPhiRad);
        double dz = -std::cos(trackAngleRad);  // 向下

        // —— 计算斜率 ——
        kx = dx / dz;
        ky = dy / dz;

        // —— 轨迹起点 ——
        double t0 = 0.;
        trackX0 = trackMidX + kx * (z0 - zMid);
        trackY0 = trackMidY + ky * (z0 - zMid);

        // — 轨迹终点（z=0.01 cm） ——
        const double trackBottomX = trackMidX + kx * (zBottom - zMid);
        const double trackBottomY = trackMidY + ky * (zBottom - zMid);

        // 检查轨迹是否越界
        const auto outsideReadout = [=](double x, double y) {
          return x < 0.0 || x >= N_STRIPS * pitch ||
                 y < 0.0 || y >= N_STRIPS * pitch;
        };
        if (outsideReadout(trackX0, trackY0) ||
            outsideReadout(trackBottomX, trackBottomY)) {
            std::cerr << "Error: configured track crosses outside readout at event "
                      << i << ": top=(" << trackX0 << ", " << trackY0
                      << "), bottom=(" << trackBottomX << ", "
                      << trackBottomY << ") cm"
                      << std::endl;
            return 1;
        }

        // —— 计算截距 ——
        bx = trackX0 - kx * z0;
        by = trackY0 - ky * z0;

        // —— 生成 track ——
        track.NewTrack(trackX0, trackY0, z0, t0, dx, dy, dz);
        
        // —— 获取 clusters ——
        double xc, yc, zc, tc, ec, extra;
        int ne;
        
        // 遍历每个 cluster
        while (track.GetCluster(xc, yc, zc, tc, ne, ec, extra)) {
            // 遍历 cluster 中的每个初级电子
            for (int ie = 0; ie < ne; ++ie) {
                double xe0, ye0, ze0, te0, ee0,dxe, dye, dze;
                track.GetElectron(ie, xe0, ye0, ze0, te0, ee0, dxe, dye, dze);
                ++nPrimary;
                
                // —— 记录初始电离位置（tree3） ——
                vec_xe0.push_back(xe0);
                vec_ye0.push_back(ye0);
                vec_ze0.push_back(ze0);
                vec_te0.push_back(te0);
                vec_ee0.push_back(ee0);
                vec_xe0_stripid.push_back((int)(xe0 / pitch));
                vec_ye0_stripid.push_back((int)(ye0 / pitch));
                
                // 漂移模拟
                const auto driftStart = std::chrono::steady_clock::now();
                const bool driftOk = drift.DriftElectron(xe0, ye0, ze0, te0, 1);
                drifttime += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - driftStart).count();

                if (!driftOk || drift.GetNumberOfElectronEndpoints() == 0) continue;
                
                int status;
                double xe1, ye1, ze1, te1;
                drift.GetElectronEndpoint(0, xe0, ye0, ze0, te0,
                                          xe1, ye1, ze1, te1, status);
                
                if(status != -1) {
                    // 漂移失败，跳过
                    continue;
                }
                constexpr double readoutZTolerance = 1.e-4;
                if (std::abs(ze1 - zBottom) > readoutZTolerance) continue;
                ++nDriftSuccess;
                
                // 计算 strip 编号
                int xe_stripid = (int)(xe1 / pitch);
                int ye_stripid = (int)(ye1 / pitch);
                
                // 边界检查
                if(xe_stripid < 0 || xe_stripid >= N_STRIPS ||
                   ye_stripid < 0 || ye_stripid >= N_STRIPS) continue;
                ++nReadoutInRange;
                
                // 记录每个电子的终点（tree2）
                vec_xe1.push_back(xe1);
                vec_ye1.push_back(ye1);
                vec_ze1.push_back(ze1);
                vec_te1.push_back(te1);
                vec_xe_stripid.push_back(xe_stripid);
                vec_ye_stripid.push_back(ye_stripid);
                
                // —— 区域映射（与 microTPCsample.C 一致） ——
                double deltax = gLib_xmax - gLib_xmin;
                double binWx = (gLib_xmax - gLib_xmin) / gLib_nbins;//x抽样bin宽
                double binWy = (gLib_ymax - gLib_ymin) / gLib_nbins;//y抽样bin宽
                int ix, iy;
                bool iscenter = false;
                int isleft = 0;
                int stripidx = (int)(xe1 / pitch);//第几个通道
                int signaljudege = (int)(xe1 / (pitch/4)) % 4;//因为是两个条合成一个通道，所以有4种情况，0表示落在左边条的左边，1表示落在左边条的右边，2表示落在右边条的左边，3表示落在右边条的右边
                int stripidy = (int)(ye1 / pitch);//第几个Y通道，Y方向是一个条一个通道
                //这里加上deltax/2是因为前面取余后去掉周期后，是参考原点的，而我们希望参考是-deltax/2。
                if(signaljudege == 0){
                    ix = (int)((std::fmod(xe1, pitch) + deltax/2) / binWx);
                    isleft = 1;
                }
                else if(signaljudege == 1){
                    ix = (int)((std::fmod(xe1, pitch) - 0.5*pitch + deltax/2) / binWx);
                    iscenter = true;
                }
                else if(signaljudege == 2){
                    ix = (int)((std::fmod(xe1, pitch) - 0.5*pitch + deltax/2) / binWx);
                    iscenter = true;
                }
                else{
                    ix = (int)((std::fmod(xe1, pitch) - pitch + deltax/2) / binWx);
                }
                iy = (int)((std::fmod(ye1, pitch)) / binWy);
                
                // 边界检查
                if(ix < 0 || ix >= gLib_nbins || iy < 0 || iy >= gLib_nbins) {
                    ++nLibraryOutOfRange;
                    continue;
                }
                
                // —— 从抽样库获取信号 ——
                if (!regionTrees[ix][iy] || region_n[ix][iy] == 0) {
                    ++nEmptyRegion;
                    continue;
                }
                
                const auto libraryStart = std::chrono::steady_clock::now();
                Long64_t entry_choose = rnd.Integer(region_n[ix][iy]);
                regionTrees[ix][iy]->GetEntry(entry_choose);
                
                // 检查信号向量非空
                bool sigEmpty = false;
                if (iscenter) {
                    for (int is = 2; is <= 4; ++is) {
                        if (!region_sig[ix][iy][is] || region_sig[ix][iy][is]->empty()) {
                            sigEmpty = true; break;
                        }
                    }
                } else {
                    for (int is = 0; is <= 1; ++is) {
                        if (!region_sig[ix][iy][is] || region_sig[ix][iy][is]->empty()) {
                            sigEmpty = true; break;
                        }
                    }
                }
                for (int is = 5; is <= 7; ++is) {
                    if (!region_sig[ix][iy][is] || region_sig[ix][iy][is]->empty()) {
                        sigEmpty = true; break;
                    }
                }
                if (sigEmpty) {
                    std::cerr << "Warning: empty signal vector in region (" << ix << "," << iy
                              << ") at event " << i << ", entry " << entry_choose
                              << std::endl;
                    continue;}
                
                // 浮点索引 + 线性插值，解决 te1 不是 tstep 整数倍时的截断误差
                double fracBegin = te1 / tstep;
                int    beginindex = (int)std::floor(fracBegin);//取整
                double alpha      = fracBegin - beginindex;//小数部分
                double w0 = 1.0 - alpha, w1 = alpha;

                if (beginindex < 0 || beginindex >= static_cast<int>(nTimeBins)) continue;

                // 所有实际使用的通道必须具有同一段可插值的有效长度，目前就是用信号的长度
                int samplesize = std::numeric_limits<int>::max();
                const int xFirst = iscenter ? 2 : 0;
                const int xLast = iscenter ? 4 : 1;
                for (int is = xFirst; is <= xLast; ++is) {
                    samplesize = std::min(samplesize, static_cast<int>(region_sig[ix][iy][is]->size()));
                }
                for (int is = 5; is <= 7; ++is) {
                    samplesize = std::min(samplesize, static_cast<int>(region_sig[ix][iy][is]->size()));
                }
                if (samplesize < 2) continue;
                ++nSampledElectrons;

                for (int k = beginindex; k < nTimeBins; ++k) {
                    // 当 j == 0 时，模板位置为 -alpha，仍在信号开始之前，因此置零。
                    // 当 j >= 1 时，j-alpha 位于模板的第 j-1 和第 j 个采样点之间：
                    //     signal(j-alpha)
                    //       = alpha * signal[j-1] + (1-alpha) * signal[j]
                    //       = w1    * signal[j-1] + w0          * signal[j]
                    int j = k - beginindex;
                    if (j >= samplesize) break;

                    const auto delayedSample = [&](const int signalidx) {
                        double value = 0.0;
                        if(j==0) {
                            return value;
                        }

                        value += w1 * region_sig[ix][iy][signalidx]->at(j - 1) + w0 * region_sig[ix][iy][signalidx]->at(j);
                        return value;

                    };

                    // "sigxtestminus1","sigxtest1","sigxtestminus2","sigzero","sigxtest2", "sigytestminus1", "sigy","sigytest1"
                    if(iscenter){//中间的话，就取2-4号信号
                        int signalidx = 2;
                        for(int indexx = stripidx-1; indexx <= stripidx+1; ++indexx) {
                            const double val = delayedSample(signalidx);
                            if (indexx >= 0 && indexx < N_STRIPS) vec_sigx[indexx][k] += val;
                            signalidx++;
                        }
                    }
                    else{//左边或者右边的话，就取0-1号信号
                        int signalidx = 0;
                        for(int indexx = stripidx-isleft; indexx <= stripidx-isleft+1; ++indexx) {
                            const double val = delayedSample(signalidx);
                            if (indexx >= 0 && indexx < N_STRIPS) vec_sigx[indexx][k] += val;
                            signalidx++;
                        }
                    }
                    // y方向的信号
                    int signalidy = 5;
                    for(int indexy = stripidy-1; indexy <= stripidy+1; ++indexy) {
                        const double val = delayedSample(signalidy);
                        if (indexy >= 0 && indexy < N_STRIPS) vec_sigy[indexy][k] += val;
                        signalidy++;
                    }
                }
                readtime += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - libraryStart).count();
            } // end electron loop
        } // end cluster loop
        
        // —— 卷积 ——
        const auto convolutionStart = std::chrono::steady_clock::now();
        for(int indexx = 0; indexx < N_STRIPS; ++indexx) {
            for (int k = 0; k < nTimeBins; ++k) {
                double convx = 0;
                double convy = 0;
                for (int j = 0; j <= k; ++j) {
                    int idx = k - j;
                    convx += vec_sigx[indexx][idx] * kernel_vec[j];
                    convy += vec_sigy[indexx][idx] * kernel_vec[j];
                }
                vec_sigx_conv[indexx][k] = convx * tstep;
                vec_sigy_conv[indexx][k] = convy * tstep;
            }
        }
        convolutiontime += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - convolutionStart).count();
        
        // —— 转换为输出格式：收集非零信号的 strip ——
        types.clear();
        stripIDs.clear();
        waveforms.clear();
        waveforms_conv.clear();
        for(int idx = 0; idx < N_STRIPS; ++idx) {
            // 检查 x 方向是否有信号
            bool hasSignalX = false;
            for(int k = 0; k < nTimeBins; ++k) {
                if(vec_sigx_conv[idx][k] != 0) { hasSignalX = true; break; }
            }
            if(hasSignalX) {
                stripIDs.push_back(idx);           // x strip: 0..N_STRIPS-1
                types.push_back(0);  // 标记为 x strip
                waveforms.push_back(vec_sigx[idx]);
                waveforms_conv.push_back(vec_sigx_conv[idx]);
            }
        }
        for(int idy = 0; idy < N_STRIPS; ++idy) {
            // 检查 y 方向是否有信号
            bool hasSignalY = false;
            for(int k = 0; k < nTimeBins; ++k) {
                if(vec_sigy_conv[idy][k] != 0) { hasSignalY = true; break; }
            }
            if(hasSignalY) {
                stripIDs.push_back(idy);      // y strip: N_STRIPS..2*N_STRIPS-1
                types.push_back(1);  // 标记为 y strip
                waveforms.push_back(vec_sigy[idy]);
                waveforms_conv.push_back(vec_sigy_conv[idy]);
            }
        }
        
        tpos->Fill();
        tinit->Fill();
        tout->Fill();
    }

    TTree metadata("simulation_metadata", "Angle-scan simulation configuration");
    int metadataNStrips = N_STRIPS;
    double metadataPitchCm = pitch;
    std::string metadataLibrary = "../analysis/samplelibrarysmall.root";
    metadata.Branch("theta", &config.simulation.theta);
    metadata.Branch("phi", &config.simulation.phi);
    metadata.Branch("seed", &config.simulation.seed);
    metadata.Branch("nEvents", &config.simulation.nEvents);
    metadata.Branch("runid", &config.simulation.runid);
    metadata.Branch("mode", &config.simulation.mode);
    metadata.Branch("deltaElectron", &config.simulation.deltaElectron);
    metadata.Branch("detectorDir", &config.simulation.detectorDir);
    metadata.Branch("gasFile", &config.simulation.gasFile);
    metadata.Branch("fieldFile", &config.simulation.fieldFile);
    metadata.Branch("simulationSampleInterval",
                    &config.simulation.sampleInterval);
    metadata.Branch("n", &config.electronics.n);
    metadata.Branch("tau", &config.electronics.tau);
    metadata.Branch("electronicsSampleInterval",
                    &config.electronics.sampleInterval);
    metadata.Branch("totalTime", &config.totalTime);
    metadata.Branch("outputDir", &config.outputDir);
    metadata.Branch("gap", &config.simulation.gap);
    metadata.Branch("tierCenter", &tierCenter);
    metadata.Branch("nStrips", &metadataNStrips);
    metadata.Branch("pitch", &metadataPitchCm);
    metadata.Branch("sampleLibrary", &metadataLibrary);
    metadata.Branch("primaryElectrons", &nPrimary);
    metadata.Branch("driftSuccess", &nDriftSuccess);
    metadata.Branch("readoutInRange", &nReadoutInRange);
    metadata.Branch("libraryOutOfRange", &nLibraryOutOfRange);
    metadata.Branch("emptyLibraryRegion", &nEmptyRegion);
    metadata.Branch("sampledElectrons", &nSampledElectrons);
    metadata.Fill();

    // 写入输出
    tout->Write();
    tpos->Write();
    tinit->Write();
    metadata.Write();
    TNamed signalCoordinates("signal_coordinates", "physical_strip");
    signalCoordinates.Write();
    fout->Close();
    fLib->Close();
    const double eventLoopSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - eventLoopStart).count();

    std::cout << "\nSample-library acceptance:" << std::endl;
    std::cout << "  primary electrons: " << nPrimary << std::endl;
    std::cout << "  drift successful: " << nDriftSuccess << std::endl;
    std::cout << "  readout in range: " << nReadoutInRange << std::endl;
    std::cout << "  library (ix, iy) out of range: " << nLibraryOutOfRange << std::endl;
    std::cout << "  empty library region: " << nEmptyRegion << std::endl;
    std::cout << "  sampled electrons: " << nSampledElectrons << std::endl;
    std::cout << "\nTiming (event loop): " << eventLoopSeconds << " s" << std::endl;
    std::cout << "  DriftElectron: " << drifttime << " s"
              << " (" << 100. * drifttime / eventLoopSeconds << "%)" << std::endl;
    std::cout << "  template read + accumulation: " << readtime << " s"
              << " (" << 100. * readtime / eventLoopSeconds << "%)" << std::endl;
    std::cout << "  waveform convolution: " << convolutiontime << " s"
              << " (" << 100. * convolutiontime / eventLoopSeconds << "%)" << std::endl;
    std::cout << "Done! Output written to " << outputFile << std::endl;
    return 0;
}
