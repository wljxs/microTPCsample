// ============================================================
// microTPCsimulationxray.C — 基于抽样库的可配置角度 X-ray 信号模拟
//
//   流程：
//     1. TrackHeed 输运 8 keV X-ray，获取电离电子
//     2. 对光子产生的每个电离电子，漂移到读出平面
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

using namespace Garfield;

// ---------- 抽样库元数据（从 samplelibrary.root 的 meta tree 读取） ----------
static double gLib_xmin, gLib_xmax, gLib_ymin, gLib_ymax;
static int    gLib_nbins;
static std::vector<double> gLib_t;  // 共享时间轴（从 meta 读取一次）
static const int kNbins = 10;       // 与 buildSampleLibrary.C 中 kNbins 一致

// ---- 信号通道名称表（14 个，与 microTPC.C 一致）----
static const int N_SIG = 8;
static const char* SIG_NAME[N_SIG] = {
        "sigxtestminus1","sigxtest1",
    "sigxtestminus2",  "sigzero","sigxtest2",
     "sigytestminus1", "sigy","sigytest1"
};//前4个信号是x方向落到两个条中间时相邻4个条产生的，中间5个信号是x方向落到条上时相邻5个条产生的，后5个信号是y方向上相邻5个条产生的

namespace {

struct SimulationConfig {
    int runid = 0;
    std::string outputFile;
    unsigned int nEvents = 50;
    double thetaDeg = 0;
    double phiDeg = 0;
    unsigned int seed = 42;
    std::string positionMode = "center";
};

void PrintUsage(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n"
        << "  --run-id N\n"
        << "  --output PATH\n"
        << "  --events N\n"
        << "  --theta-deg V       polar angle from downward z (default: 0)\n"
        << "  --phi-deg V         azimuth from +x toward +y (default: 0)\n"
        << "  --seed N            non-zero deterministic seed (default: 42)\n"
        << "  --position-mode center|uniform\n"
        << "\nLegacy positional form is still accepted:\n"
        << "  " << program << " [runid] [output-root-file] [n-events]\n";
}

SimulationConfig ParseConfig(int argc, char** argv) {
    SimulationConfig config;
    if (argc > 1 && argv[1][0] != '-') {
        config.runid = std::stoi(argv[1]);
        if (argc > 2) config.outputFile = argv[2];
        if (argc > 3) config.nEvents = static_cast<unsigned int>(std::stoul(argv[3]));
        config.seed = config.runid == 0
            ? 42u : static_cast<unsigned int>(config.runid) * 12345u;
    } else {
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--help") {
                PrintUsage(argv[0]);
                std::exit(0);
            }
            if (index + 1 >= argc) {
                throw std::invalid_argument("missing value for " + argument);
            }
            const std::string value = argv[++index];
            if (argument == "--run-id") config.runid = std::stoi(value);
            else if (argument == "--output") config.outputFile = value;
            else if (argument == "--events") {
                config.nEvents = static_cast<unsigned int>(std::stoul(value));
            } else if (argument == "--theta-deg") config.thetaDeg = std::stod(value);
            else if (argument == "--phi-deg") config.phiDeg = std::stod(value);
            else if (argument == "--seed") {
                config.seed = static_cast<unsigned int>(std::stoul(value));
            } else if (argument == "--position-mode") config.positionMode = value;
            else throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (config.nEvents == 0 || config.seed == 0 ||
        !std::isfinite(config.thetaDeg) || config.thetaDeg < 0.0 ||
        config.thetaDeg >= 90.0 || !std::isfinite(config.phiDeg) ||
        (config.positionMode != "center" && config.positionMode != "uniform")) {
        throw std::invalid_argument("invalid events, seed, angle, or position mode");
    }
    if (config.outputFile.empty()) {
        std::ostringstream name;
        const int phiLabel = static_cast<int>(std::lround(
            std::fmod(std::fmod(config.phiDeg, 360.0) + 360.0, 360.0)));
        name << "../result/angle_scan/theta" << std::setw(3) << std::setfill('0')
             << static_cast<int>(std::lround(config.thetaDeg))
             << "_phi" << std::setw(3) << std::setfill('0') << phiLabel
             << "_seed" << config.seed << ".root";
        config.outputFile = name.str();
    }
    return config;
}

}  // namespace

// ============================================================
// main — 主模拟程序
// ============================================================
int main(int argc, char *argv[])
{
    // ---------- 参数 ----------
    SimulationConfig config;
    try {
        config = ParseConfig(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << std::endl;
        PrintUsage(argv[0]);
        return 1;
    }
    const int runid = config.runid;
    const std::string outputFile = config.outputFile;

    std::string detector = "microTPC5mmone";
    constexpr double pitch = 0.04;
    const unsigned int nEvents = config.nEvents;
    constexpr int N_STRIPS = 20;                 // x/y 方向各 20 条
    const double trackAngleRad = config.thetaDeg * M_PI / 180.0;
    const double trackPhiRad = config.phiDeg * M_PI / 180.0;
    TRandom3 rnd(config.seed);
    RandomEngineRoot garfieldRandomEngine(config.seed);
    Random::SetEngine(garfieldRandomEngine);
    std::cout << "Simulation configuration: theta=" << config.thetaDeg
              << " deg, phi=" << config.phiDeg << " deg, events=" << nEvents
              << ", seed=" << config.seed
              << ", position-mode=" << config.positionMode << std::endl;

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
    gas.LoadGasFile("../" + detector + "/ar_45_cf4_40_co2_15100-100000.gas");
    constexpr double rPenning = 0.20;
    constexpr double lambdaPenning = 0.;
    gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");
    gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");
    gas.SetMaxElectronEnergy(200.);

    // ---------- 电场 ----------
    ComponentComsol fm;
    fm.Initialise(("../" + detector + "/mesh.mphtxt").c_str(),
                  ("../" + detector + "/mplist.txt").c_str(),
                  ("../" + detector + "/field.txt").c_str(), "mm");
    fm.EnablePeriodicityX();
    fm.EnablePeriodicityY();
    fm.SetGas(&gas);

    // ---------- Sensor ----------
    Sensor sensor(&fm);
    sensor.SetArea(0, 0, 0.01, N_STRIPS * pitch, N_STRIPS * pitch, 0.51);


    const unsigned int nTimeBins = 140;
    const double tmin = 0.;
    const double tmax = 700;
    const double tstep = (tmax - tmin) / nTimeBins;
    sensor.SetTimeWindow(tmin, tstep, nTimeBins);
    Shaper shaper(1, 50., 1., "unipolar");
    // 抽样库保存的是卷积前信号；所有电子叠加后在此处统一成形。
    // 直接调用 Garfield 的 Shaper，避免手写核与 Shaper 参数发生漂移。
    std::vector<double> kernel_vec(nTimeBins);
    for (int i = 0; i < nTimeBins; ++i) {
        double t = tmin + i * tstep;
        kernel_vec[i] = shaper.Shape(t);
    }

    // ---------- TrackHeed（8 keV X-ray） ----------
    TrackHeed track(&sensor);

    // ---------- Drift ----------
    AvalancheMC drift(&sensor);
    drift.SetTimeWindow(0, 1400);
    drift.SetDistanceSteps(1.e-3);

    // ---------- 输出 Tree ----------
    const std::filesystem::path outputPath(outputFile);
    if (!outputPath.parent_path().empty()) {
        std::filesystem::create_directories(outputPath.parent_path());
    }
    TFile *fout = new TFile(outputFile.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << outputFile << std::endl;
        return 1;
    }

    int event;
    // 时间轴
    std::vector<double> vec_t;
    for(int i = 0; i < nTimeBins; ++i) vec_t.push_back(tmin + i * tstep);

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
    // 性能统计：区分 Garfield 漂移、ROOT 模板读取/叠加和电子学卷积。
    double driftSeconds = 0.;
    double librarySeconds = 0.;
    double convolutionSeconds = 0.;
    const auto eventLoopStart = std::chrono::steady_clock::now();
    // ==================== 事件循环 ====================
    // track.DisableDeltaElectronTransport();
    for (unsigned int i = 0; i < nEvents; ++i) {

        // —— 重置信号累加器 ——
        for (int indexx = 0; indexx < N_STRIPS; ++indexx) {
            vec_sigx[indexx].assign(nTimeBins, 0);
            vec_sigy[indexx].assign(nTimeBins, 0);
            vec_sigx_conv[indexx].assign(nTimeBins, 0);
            vec_sigy_conv[indexx].assign(nTimeBins, 0);
        }

        if(i % 100 == 0) std::cout << "Event: " << i << "/" << nEvents << std::endl;
        sensor.ClearSignal();

        // —— 清空本 event 的电子位置向量 ——
        vec_xe1.clear(); vec_ye1.clear(); vec_ze1.clear(); vec_te1.clear();
        vec_xe_stripid.clear(); vec_ye_stripid.clear();
        vec_xe0.clear(); vec_ye0.clear(); vec_ze0.clear(); vec_te0.clear(); vec_ee0.clear();
        vec_xe0_stripid.clear(); vec_ye0_stripid.clear();

        // 固定轨迹在漂移区中点附近穿过读出中心，避免角度扫描受到边界距离影响。
        double z0 = 0.49;   // 从探测器顶部附近开始
        constexpr double zMid = 0.25;
        double t0 = 0.;

        trackMidX = 0.5 * N_STRIPS * pitch;
        trackMidY = 0.5 * N_STRIPS * pitch;
        if (config.positionMode == "uniform") {
            trackMidX += rnd.Uniform(-0.5 * pitch, 0.5 * pitch);
            trackMidY += rnd.Uniform(-0.5 * pitch, 0.5 * pitch);
        }
        const double dx = std::sin(trackAngleRad) * std::cos(trackPhiRad);
        const double dy = std::sin(trackAngleRad) * std::sin(trackPhiRad);
        const double dz = -std::cos(trackAngleRad);  // 向下，默认垂直入射
        constexpr double eGamma = 8.0e3; // 8 keV X-ray
        // —— 计算斜率与截距（z 为自变量） ——
        kx = dx / dz;
        ky = dy / dz;
        trackX0 = trackMidX + kx * (z0 - zMid);
        trackY0 = trackMidY + ky * (z0 - zMid);
        constexpr double zBottom = 0.01;
        const double trackBottomX = trackMidX + kx * (zBottom - zMid);
        const double trackBottomY = trackMidY + ky * (zBottom - zMid);
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
        bx = trackX0 - kx * z0;
        by = trackY0 - ky * z0;

        event = i + runid * nEvents;

        // 返回值包含本次光子输运产生的电子，不使用 NewTrack 的 cluster 迭代接口。
        const auto cluster = track.TransportPhoton(
            trackX0, trackY0, z0, t0, eGamma, dx, dy, dz);
        std::cout << "  Photon event " << event << ": "
                  << cluster.electrons.size() << " ionization electrons" << std::endl;
        // 未发生电离的光子保留为空事件。
        for (const auto& electron : cluster.electrons) {
            double xe0 = electron.x, ye0 = electron.y, ze0 = electron.z;
            double te0 = electron.t;
            const double ee0 = electron.e;
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
            driftSeconds += std::chrono::duration<double>(
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
            double binWx = (gLib_xmax - gLib_xmin) / gLib_nbins;
            double binWy = (gLib_ymax - gLib_ymin) / gLib_nbins;
            int ix, iy;
            bool iscenter = false;
            int isleft = 0;
            int stripidx = (int)(xe1 / pitch);
            int signaljudege = (int)(xe1 / (pitch/4)) % 4;
            int stripidy = (int)(ye1 / pitch);

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
            if (sigEmpty) continue;

            // 浮点索引 + 线性插值，解决 te1 不是 tstep 整数倍时的截断误差
            double fracBegin = te1 / tstep;
            int    beginindex = (int)std::floor(fracBegin);
            double alpha      = fracBegin - beginindex;
            double w0 = 1.0 - alpha, w1 = alpha;

            if (beginindex < 0 || beginindex >= static_cast<int>(nTimeBins)) continue;

            // 所有实际使用的通道必须具有同一段可插值的有效长度。
            int libSize = std::numeric_limits<int>::max();
            const int xFirst = iscenter ? 2 : 0;
            const int xLast = iscenter ? 4 : 1;
            for (int is = xFirst; is <= xLast; ++is) {
                libSize = std::min(libSize, static_cast<int>(region_sig[ix][iy][is]->size()));
            }
            for (int is = 5; is <= 7; ++is) {
                libSize = std::min(libSize, static_cast<int>(region_sig[ix][iy][is]->size()));
            }
            if (libSize < 2) continue;
            ++nSampledElectrons;

            for (int k = beginindex; k < nTimeBins; ++k) {
                int j = k - beginindex;               // 库内索引
                if (j >= libSize) break;

                // 将模板延迟 alpha 个采样点。模板起点之前视为零；使用
                // sample[j - 1] 而不是 sample[j + 1]，避免把波形向前推进。
                const auto delayedSample = [&](const int signalidx) {
                    double value = w0 * region_sig[ix][iy][signalidx]->at(j);
                    if (j > 0) {
                        value += w1 * region_sig[ix][iy][signalidx]->at(j - 1);
                    }
                    return value;
                };

                if(iscenter){
                    int signalidx = 2;
                    for(int indexx = stripidx-1; indexx <= stripidx+1; ++indexx) {
                        const double val = delayedSample(signalidx);
                        if (indexx >= 0 && indexx < N_STRIPS) vec_sigx[indexx][k] += val;
                        signalidx++;
                    }
                }
                else{
                    int signalidx = 0;
                    for(int indexx = stripidx-isleft; indexx <= stripidx-isleft+1; ++indexx) {
                        const double val = delayedSample(signalidx);
                        if (indexx >= 0 && indexx < N_STRIPS) vec_sigx[indexx][k] += val;
                        signalidx++;
                    }
                }
                int signalidy = 5;
                for(int indexy = stripidy-1; indexy <= stripidy+1; ++indexy) {
                    const double val = delayedSample(signalidy);
                    if (indexy >= 0 && indexy < N_STRIPS) vec_sigy[indexy][k] += val;
                    signalidy++;
                }
            }
            librarySeconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - libraryStart).count();
        } // end electron loop

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
        convolutionSeconds += std::chrono::duration<double>(
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
    double metadataThetaDeg = config.thetaDeg;
    double metadataPhiDeg = config.phiDeg;
    unsigned int metadataSeed = config.seed;
    unsigned int metadataEvents = config.nEvents;
    int metadataRunid = config.runid;
    int metadataNStrips = N_STRIPS;
    double metadataPitchCm = pitch;
    double metadataZMidCm = 0.25;
    std::string metadataPositionMode = config.positionMode;
    std::string metadataLibrary = "../analysis/samplelibrarysmall.root";
    metadata.Branch("theta_deg", &metadataThetaDeg);
    metadata.Branch("phi_deg", &metadataPhiDeg);
    metadata.Branch("seed", &metadataSeed);
    metadata.Branch("events", &metadataEvents);
    metadata.Branch("run_id", &metadataRunid);
    metadata.Branch("n_strips", &metadataNStrips);
    metadata.Branch("pitch_cm", &metadataPitchCm);
    metadata.Branch("z_mid_cm", &metadataZMidCm);
    metadata.Branch("position_mode", &metadataPositionMode);
    metadata.Branch("sample_library", &metadataLibrary);
    metadata.Branch("primary_electrons", &nPrimary);
    metadata.Branch("drift_success", &nDriftSuccess);
    metadata.Branch("readout_in_range", &nReadoutInRange);
    metadata.Branch("library_out_of_range", &nLibraryOutOfRange);
    metadata.Branch("empty_library_region", &nEmptyRegion);
    metadata.Branch("sampled_electrons", &nSampledElectrons);
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
    std::cout << "  DriftElectron: " << driftSeconds << " s"
              << " (" << 100. * driftSeconds / eventLoopSeconds << "%)" << std::endl;
    std::cout << "  template read + accumulation: " << librarySeconds << " s"
              << " (" << 100. * librarySeconds / eventLoopSeconds << "%)" << std::endl;
    std::cout << "  waveform convolution: " << convolutionSeconds << " s"
              << " (" << 100. * convolutionSeconds / eventLoopSeconds << "%)" << std::endl;
    std::cout << "Done! Output written to " << outputFile << std::endl;
    return 0;
}
