// ============================================================
// microTPCsample.C — 基于抽样库的快速信号模拟
//
//   流程：
//     1. 源电子在指定深度层中心产生并漂移至读出面
//     2. 对每个雪崩电子，按其(ex1, ey2)在抽样库中查表
//     3. 将抽样库中的信号作为该电子漂移到 z=0.01 之后的部分
//     4. 聚合所有电子信号写入新 ROOT 文件
//
//   依赖：需先运行 analysis/sampleconstruction.C 中的 buildSampleLibrary()
// ============================================================

#include <Rtypes.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TTree.h>
#include <TRandom3.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <random>
#include <string>

#include "Garfield/AvalancheMicroscopic.hh"
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
#include "Config.h"

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

// ============================================================
// main — 主模拟程序
// ============================================================
int main(int argc, char *argv[])
{
    // 与 microTPCsimulation 使用同一套配置和命名。
    const std::string configFile = argc >= 2 ? argv[1] : "../config/default.jsonc";
    Config config = ReadConfig(configFile);

    const std::string simulationFile = config.outputDir + "/events.root";
    TFile simulationOutput(simulationFile.c_str(), "READ");
    auto* simulationMetadata = static_cast<TTree*>(
        simulationOutput.Get("simulation_metadata"));
    if (simulationOutput.IsZombie() || !simulationMetadata) {
        std::cerr << "Error: cannot read simulation metadata from "
                  << simulationFile << std::endl;
        return 1;
    }

    std::vector<double>* tierCenterFromMetadata = nullptr;
    simulationMetadata->SetBranchAddress("tierCenter", &tierCenterFromMetadata);
    simulationMetadata->GetEntry(0);

    if (config.simulation.nEvents == 0 ||
        config.simulation.sampleInterval <= 0. ||
        config.electronics.n == 0 || config.electronics.tau <= 0. ||
        config.totalTime <= 0. ||
        !tierCenterFromMetadata || tierCenterFromMetadata->empty()) {
        std::cerr << "Error: invalid simulation metadata in "
                  << simulationFile << std::endl;
        return 1;
    }
    const std::vector<double> tierCenter = *tierCenterFromMetadata;
    simulationOutput.Close();

    constexpr double pitch_cm = 0.04;
    RandomEngineRoot garfieldRandomEngine;
    garfieldRandomEngine.SetSeed(config.simulation.seed);
    Random::SetEngine(garfieldRandomEngine);

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

        if (gLib_nbins <= 0 || gLib_nbins > kNbins || gLib_t.empty()) {
            std::cerr << "Error: unsupported sample-library metadata: nbins="
                      << gLib_nbins << ", t.size=" << gLib_t.size() << std::endl;
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

    // ---------- Sensor ----------
    Sensor sensor(&fm);
    sensor.SetArea(0, 0, 0.01, 10. * pitch_cm, 10 * pitch_cm, 0.51);

    const unsigned int nTimeBins = static_cast<unsigned int>(
        std::ceil(config.totalTime / config.simulation.sampleInterval));
    const double tmin = 0.;
    const double tstep = config.simulation.sampleInterval;
    sensor.SetTimeWindow(tmin, tstep, nTimeBins);
    Shaper shaper(config.electronics.n, config.electronics.tau, 1., "unipolar");
    std::vector<double> kernel_vec(nTimeBins);
    for (int i = 0; i < nTimeBins; ++i) {
        double t = tmin + i * tstep;
        kernel_vec[i] = shaper.Shape(t);
    }

    // ---------- Avalanche ----------
    AvalancheMicroscopic aval(&sensor);
    // ---------- Drift ----------

    AvalancheMC drift(&sensor);
    //drift.SetTimeSteps(50);
    drift.SetTimeWindow(0, 1400);
    drift.SetDistanceSteps(1.e-3);
    // metadata 中有几个中心位置，就一次生成几个 tier 文件。
    const std::string tierOutputDir = config.outputDir + "/tiers";
    std::filesystem::create_directories(tierOutputDir);
    for (size_t tier = 0; tier < tierCenter.size(); ++tier) {
    const double template_z_cm = tierCenter[tier];
    TRandom3 rnd(config.simulation.seed +
                 1009u * static_cast<unsigned int>(tier));

    // ---------- 输出 Tree ----------
    const std::string outputFile = tierOutputDir + "/tier" +
                                   std::to_string(tier) + ".root";
    TFile *fout = TFile::Open(outputFile.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: cannot create " << outputFile << std::endl;
        return 1;
    }

    int event;
    // 信号累加器（10 条 x + 10 条 y，以及时间轴）。输出坐标以最终落点
    // 所在 strip 为中心：sigx2/sigy2 是落点中心响应，相邻编号保存邻条响应。
    std::vector<double> t_total;
    std::vector<double> vec_t;
    std::vector<double> vec_sigx[10];
    std::vector<double> vec_sigy[10];
    std::vector<double> vec_sigx_conv[10];
    std::vector<double> vec_sigy_conv[10];
    TTree *tout = new TTree("tree", "Sampled Avalanche data");
    tout->Branch("event", &event, "event/I");
    tout->Branch("t",     &vec_t);
    for(int indexx = 0; indexx < 10; ++indexx) {
        tout->Branch(Form("sigx%d", indexx), &vec_sigx[indexx]);
        tout->Branch(Form("sigy%d", indexx), &vec_sigy[indexx]);
        tout->Branch(Form("sigx%d_conv", indexx), &vec_sigx_conv[indexx]);
        tout->Branch(Form("sigy%d_conv", indexx), &vec_sigy_conv[indexx]);
    }
    for(int i=0;i<nTimeBins;++i) vec_t.push_back(tmin + i*tstep);

    // 记录电子初始/终止位置的 tree
    TTree *tpos = new TTree("tree2", "Electron positions");
    double xe1, ye1, ze1, te1;
    tpos->Branch("event", &event, "event/I");
    tpos->Branch("xe1",   &xe1,   "xe1/D");
    tpos->Branch("ye1",   &ye1,   "ye1/D");
    tpos->Branch("ze1",   &ze1,   "ze1/D");
    tpos->Branch("te1",   &te1,   "te1/D");
    // ==================== 事件循环 ====================
    for (unsigned int i = 0; i < config.simulation.nEvents; ++i) {
        for (int indexx = 0; indexx < 10; ++indexx) {
            vec_sigx[indexx].clear();
            vec_sigy[indexx].clear();
            vec_sigx_conv[indexx].clear();
            vec_sigy_conv[indexx].clear();
        }
        for (int indexx = 0; indexx < 10; ++indexx) {
            vec_sigx[indexx].resize(nTimeBins, 0);
            vec_sigy[indexx].resize(nTimeBins, 0);
            vec_sigx_conv[indexx].resize(nTimeBins, 0);
            vec_sigy_conv[indexx].resize(nTimeBins, 0);
        }
        if(i % 100 == 0) {
            std::cout << "Event: " << i << "/"
                      << config.simulation.nEvents << std::endl;
        }

        // 清除累加器
        t_total.clear();


        sensor.ClearSignal();

        // 随机电子起始位置
        double x0 = 0.1 + pitch_cm * rnd.Uniform(-0.5, 0.5);
        double y0 = 0.1 + pitch_cm * rnd.Uniform(-0.5, 0.5);
        double z0 = template_z_cm;
        double t0 = 0.;
        // 漂移模拟

        const bool driftOk = drift.DriftElectron(x0, y0, z0, t0, 1);
        if (!driftOk || drift.GetNumberOfElectronEndpoints() == 0) continue;

        int status;
        drift.GetElectronEndpoint(0, x0, y0, z0, t0, xe1, ye1, ze1, te1,status);
        if (status != -1) continue;
        constexpr double readoutZ = 0.01;
        constexpr double readoutZTolerance = 1.e-4;
        if (std::abs(ze1 - readoutZ) > readoutZTolerance) continue;

        event = i;

        // 源电子位置 → 区域索引
        double deltax = gLib_xmax - gLib_xmin;
        double deltay = gLib_ymax - gLib_ymin;
        double dx = (gLib_xmax - gLib_xmin) / gLib_nbins;
        double dy = (gLib_ymax - gLib_ymin) / gLib_nbins;
        int ix,iy;
        bool iscenter = false;
        int isleft = 0;
        int stripidx = (int)(xe1/pitch_cm);//x方向上一个pitch有两个条，先判断在哪个条，再判断在条的哪一侧
        int signaljudege = (int)(xe1/(pitch_cm/4))%4;//判断在电极的哪一侧（左边，中间，右边）
        int stripidy = (int)(ye1/pitch_cm);//y方向上一个pitch一个条

        if(signaljudege == 0){
            ix=  (int)((std::fmod(xe1, pitch_cm)+deltax/2) /dx);
            isleft = 1;
        }
        else if(signaljudege == 1){
            ix=  (int)((std::fmod(xe1,pitch_cm)-0.5*pitch_cm+deltax/2)/dx);
            iscenter = true;
        }
        else if(signaljudege == 2){
            ix=  (int)((std::fmod(xe1,pitch_cm)-0.5*pitch_cm+deltax/2)/dx);
            iscenter = true;
        }
        else{
            ix=  (int)((std::fmod(xe1,pitch_cm)-pitch_cm+deltax/2)/dx);
        }

        iy = int((std::fmod(ye1,pitch_cm))/dy); //因为y方向上一个pitch一个条
        // std::cout << "  Avalanche electron at (" << xe1 << ", " << ye1 << ", " << ze1 << " , " << te1 << ")" << std::endl;
        // std::cout << stripidx << " " << signaljudege << " " << stripidy << std::endl;
        // std::cout << std::fmod(xe1, pitch_cm) << " " << std::fmod(ye1, pitch_cm) << " " << dx << " " << dy << std::endl;
        // std::cout << "  Mapped to region (" << ix << ", " << iy << ")" << std::endl;

        // 从抽样库获取信号
        if (ix < 0 || ix >= gLib_nbins || iy < 0 || iy >= gLib_nbins ||
            stripidx < 1 || stripidx >= 9 || stripidy < 1 || stripidy >= 9) {
            std::cerr << "Warning: endpoint outside template/strip acceptance: ("
                      << xe1 << ", " << ye1 << ") -> region (" << ix << ", " << iy << ")" << std::endl;
            continue;
        }
        if (!regionTrees[ix][iy] || region_n[ix][iy] == 0) {
            std::cerr << "  Warning: empty region (" << ix << "," << iy
                      << "), skipping event" << std::endl;
            continue;
        }
        Long64_t entry_choose = rnd.Integer(region_n[ix][iy]);
        regionTrees[ix][iy]->GetEntry(entry_choose);

        // 检查抽样库中的信号向量是否非空
        bool sigEmpty = false;
        if (iscenter) {
            for (int is = 2; is <= 4; ++is) {
                if (!region_sig[ix][iy][is] || region_sig[ix][iy][is]->empty()) {
                    sigEmpty = true; break;
                }
            }
        } else {
            for (int is = isleft ? 0 : 0; is <= (isleft ? 1 : 1); ++is) {
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
            std::cerr << "  Warning: empty signal vector in region (" << ix << "," << iy
                      << "), skipping event" << std::endl;
            continue;
        }

        // 浮点索引 + 线性插值，解决 te1 不是 tstep 整数倍时的截断误差
        double fracBegin = te1 / tstep;
        int    beginindex = (int)fracBegin;
        double alpha      = fracBegin - beginindex;
        double w0 = 1.0 - alpha, w1 = alpha;

        // 所有实际使用的通道必须具有同一段有效长度。
        int libSize = std::numeric_limits<int>::max();
        const int xFirst = iscenter ? 2 : 0;
        const int xLast = iscenter ? 4 : 1;
        for (int is = xFirst; is <= xLast; ++is) {
            libSize = std::min(libSize, static_cast<int>(region_sig[ix][iy][is]->size()));
        }
        for (int is = 5; is <= 7; ++is) {
            libSize = std::min(libSize, static_cast<int>(region_sig[ix][iy][is]->size()));
        }
        if (beginindex < 0 || beginindex >= static_cast<int>(nTimeBins) || libSize < 1) {
            continue;
        }

        for (int k = beginindex; k < nTimeBins; ++k) {
            int j = k - beginindex;               // 库内索引
            if (j >= libSize) break;

            // 将模板延迟 alpha 个采样点；模板起点之前视为零。
            const auto delayedSample = [&](const int signalidx) {
                double value = w0 * region_sig[ix][iy][signalidx]->at(j);
                if (j > 0) {
                    value += w1 * region_sig[ix][iy][signalidx]->at(j - 1);
                }
                return value;
            };

            if(iscenter){
                int signalidx = 2;
                for(int outputStrip = 1; outputStrip <= 3; ++outputStrip) {
                    vec_sigx[outputStrip][k] += delayedSample(signalidx);
                    signalidx++;
                }
            }
            else{
                int signalidx = 0;
                const int outputStart = isleft ? 1 : 2;
                for(int outputStrip = outputStart; outputStrip <= outputStart + 1; ++outputStrip) {
                    vec_sigx[outputStrip][k] += delayedSample(signalidx);
                    signalidx++;
                }
            }
            int signalidy = 5;
            for(int outputStrip = 1; outputStrip <= 3; ++outputStrip) {
                vec_sigy[outputStrip][k] += delayedSample(signalidy);
                signalidy++;
            }
        }

        for(int indexx = 0; indexx < 10; ++indexx) {
            for (int k = 0; k < nTimeBins; ++k) {
                double convx = 0;
                double convy = 0;
                for (int j = 0; j <= k; ++j) {
                    int idx = k - j;
                    convx += vec_sigx[indexx][idx] * kernel_vec[j];
                    convy += vec_sigy[indexx][idx] * kernel_vec[j];

                }
                vec_sigx_conv[indexx][k] = convx*tstep;
                vec_sigy_conv[indexx][k] = convy*tstep;
            }
        }

        tout->Fill();
        // // 遍历每个雪崩电子（填 tpos tree）
        // // 时间轴（直接复用抽样库的共享时间向量）
        // t_total = gLib_t;

        tpos->Fill();
    }

    // 写入输出
    tout->Write();
    tpos->Write();
    TNamed signalCoordinates("signalCoordinates", "endpoint_centered");
    TNamed schemaVersion("tierSchemaVersion", "2");
    TParameter<double> templateZ("templateZ", template_z_cm);
    signalCoordinates.Write();
    schemaVersion.Write();
    templateZ.Write();
    fout->Close();

    // 抽样库的 TFile 由 regionTrees 引用，不能关闭
    std::cout << "\nDone! Direct template for tier " << tier << "/"
              << tierCenter.size()
              << " at z=" << template_z_cm * 10. << " mm written to "
              << outputFile << std::endl;
    }

    fLib->Close();
    std::cout << "Generated " << tierCenter.size()
              << " tier files in " << tierOutputDir << std::endl;
    return 0;
}
