#ifndef MICROTPC_CONFIG_H
#define MICROTPC_CONFIG_H

#include <fstream>
#include <iostream>
#include <string>

#include <nlohmann/json.hpp>

struct Simulation {
    int runid = 0;
    unsigned int nEvents = 5000;
    std::string detectorDir = "../microTPC5mmone"; // 探测器网格和电场文件所在目录
    std::string gasFile = "ar_45_cf4_40_co2_15100-100000.gas"; // 气体文件，相对于 detectorDir
    std::string fieldFile = "field.txt"; // 电场文件，相对于 detectorDir
    double theta = 20.;            // 极角，单位 degree
    double phi = 90.;              // 方位角，单位 degree
    unsigned int seed = 42;
    std::string mode = "center";   // 入射位置模式："center" 或 "uniform"
    bool deltaElectron = false;    // 是否启用 delta 电子输运
    double sampleInterval = 5.;    // 模拟波形采样间隔，单位 ns
    double gap = 0.5;              // 漂移区上边界z坐标，单位 cm
};

struct Electronics {
    unsigned int n = 1;            // 成形阶数
    double tau = 50.;              // 成形时间常数，单位 ns
    double sampleInterval = 25.;   // 电子学采样间隔，单位 ns
};

struct Config {
    std::string outputDir = "../result/default";
    double totalTime = 700.;       // 波形总时间，单位 ns
    Simulation simulation;
    Electronics electronics;
};

inline Config ReadConfig(const std::string& filename) {
    Config config;
    std::ifstream file(filename);
    if (!file) {
        std::cout << "Config file " << filename
                  << " not found; using default values." << std::endl;
        return config;
    }

    // 最后一个 true 表示允许配置文件中使用注释。
    const nlohmann::json j = nlohmann::json::parse(
        file, nullptr, true, true);
    config.outputDir = j.value("outputDir", config.outputDir);
    config.totalTime = j.value("totalTime", config.totalTime);

    if (j.contains("simulation")) {
        const auto& simulation = j["simulation"];
        config.simulation.runid = simulation.value("runid", config.simulation.runid);
        config.simulation.nEvents = simulation.value("nEvents", config.simulation.nEvents);
        config.simulation.detectorDir = simulation.value(
            "detectorDir", config.simulation.detectorDir);
        config.simulation.gasFile = simulation.value(
            "gasFile", config.simulation.gasFile);
        config.simulation.fieldFile = simulation.value(
            "fieldFile", config.simulation.fieldFile);
        config.simulation.theta = simulation.value("theta", config.simulation.theta);
        config.simulation.phi = simulation.value("phi", config.simulation.phi);
        config.simulation.seed = simulation.value("seed", config.simulation.seed);
        config.simulation.mode = simulation.value("mode", config.simulation.mode);
        config.simulation.deltaElectron = simulation.value(
            "deltaElectron", config.simulation.deltaElectron);
        config.simulation.sampleInterval = simulation.value(
            "sampleInterval", config.simulation.sampleInterval);
        config.simulation.gap = simulation.value("gap", config.simulation.gap);
    }

    if (j.contains("electronics")) {
        const auto& electronics = j["electronics"];
        config.electronics.n = electronics.value("n", config.electronics.n);
        config.electronics.tau = electronics.value("tau", config.electronics.tau);
        config.electronics.sampleInterval = electronics.value(
            "sampleInterval", config.electronics.sampleInterval);
    }

    return config;
}

#endif
