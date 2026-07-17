#pragma once

#include <nlohmann/json.hpp>
#include <iostream>
#include <algorithm>

using json = nlohmann::json;

// ========== 反卷积分析参数集（可通过 config JSON 的 "cuts" 字段配置）==========
// 新增参数：在 struct 中加一行成员 + 在 LoadFromJson 中加一行读取 + 在 Print 中加一行输出
struct DeconvCuts {
// === 可配置参数（通过 JSON 的 "cuts" 字段覆盖默认值） ===
    int    n_tier           = 5;    // 电荷层数
    int    apvsamples       = 27;   // APV25 采样点数

    // === 探测器物理参数 ===
    double pitchwidth[2]   = {0.4, 0.4};   // [X, Y] strip pitch (mm)
    int    voltage          = 250;          // 工作电压 (V)
    double thickness        = 5.0;          // 增益层厚度 (mm)
    double vdrift           = 0.04;         // 电子漂移速度 (mm/ns)
    double drift_ratio       = 1.0;          // 漂移速度/厚度比
    double z_ref_position   = 2.5;          // 参考 z 位置 (mm)，通常取厚度中间

    // === 时序参数 ===
    double delayindex2d[2]  = {97/25.0, 88/25.0}; // [X, Y] 粗延迟 (采样点)
    int    delayindex        = 0;            // 粗对齐补偿
    int    fit_point         = 22;           // 成形函数拟合点数

    // === StripHit 级 cut ===
    double rnormlim[2]      = {0.5, 0.5};   // [X, Y] 反卷积拟合 rnorm 上限

    // === Cluster 级 cut ===
    int    clusterMinEnergy[2]  = {300, 410}; // [X, Y] cluster 最小能量
    int    clusterMinMaxAmp     = 40;         // cluster 最小最大幅度
    double clusterChi2ReLimit   = 0.18;       // cluster chi2_re 上限
    int    clusterMinTierCharge = 4;          // cluster 最少有效电荷层数

    // ========== 从 JSON 加载 ==========
    void LoadFromJson(const json& j) {
        if (!j.contains("cuts") || !j["cuts"].is_object()) return;
        const auto& c = j["cuts"];

        n_tier           = c.value("n_tier",           n_tier);
        apvsamples       = c.value("apvsamples",       apvsamples);

        // 探测器物理参数
        if (c.contains("pitchwidth") && c["pitchwidth"].is_array() && c["pitchwidth"].size() >= 2) {
            pitchwidth[0] = c["pitchwidth"][0].get<double>();
            pitchwidth[1] = c["pitchwidth"][1].get<double>();
        }
        voltage          = c.value("voltage",          voltage);
        thickness        = c.value("thickness",        thickness);
        vdrift           = c.value("vdrift",           vdrift);
        drift_ratio       = c.value("drift_ratio",      drift_ratio);
        z_ref_position   = c.value("z_ref_position",   z_ref_position);

        // 时序参数
        if (c.contains("delayindex2d") && c["delayindex2d"].is_array() && c["delayindex2d"].size() >= 2) {
            delayindex2d[0] = c["delayindex2d"][0].get<double>();
            delayindex2d[1] = c["delayindex2d"][1].get<double>();
        }
        delayindex = c.value("delayindex", delayindex);
        fit_point  = c.value("fit_point",  fit_point);

        // StripHit 级 cut
        if (c.contains("rnormlim") && c["rnormlim"].is_array() && c["rnormlim"].size() >= 2) {
            rnormlim[0] = c["rnormlim"][0].get<double>();
            rnormlim[1] = c["rnormlim"][1].get<double>();
        }

        // Cluster 级 cut
        if (c.contains("clusterMinEnergy") && c["clusterMinEnergy"].is_array() && c["clusterMinEnergy"].size() >= 2) {
            clusterMinEnergy[0] = c["clusterMinEnergy"][0].get<int>();
            clusterMinEnergy[1] = c["clusterMinEnergy"][1].get<int>();
        }
        clusterMinMaxAmp     = c.value("clusterMinMaxAmp",     clusterMinMaxAmp);
        clusterChi2ReLimit   = c.value("clusterChi2ReLimit",   clusterChi2ReLimit);
        clusterMinTierCharge = c.value("clusterMinTierCharge", clusterMinTierCharge);
    }

    // ========== 打印当前配置 ==========
    void Print() const {
        std::cout << "  === DeconvCuts ===" << std::endl;
        std::cout << "  n_tier: " << n_tier << ", apvsamples: " << apvsamples << std::endl;
        std::cout << "  pitchwidth: [" << pitchwidth[0] << ", " << pitchwidth[1] << "]" << std::endl;
        std::cout << "  voltage: " << voltage << ", thickness: " << thickness
                  << ", vdrift: " << vdrift << ", drift_ratio: " << drift_ratio << std::endl;
        std::cout << "  z_ref_position: " << z_ref_position << std::endl;
        std::cout << "  delayindex2d: [" << delayindex2d[0] << ", " << delayindex2d[1] << "]"
                  << ", delayindex: " << delayindex << ", fit_point: " << fit_point << std::endl;
        std::cout << "  --- cuts ---" << std::endl;
        std::cout << "  rnormlim: [" << rnormlim[0] << ", " << rnormlim[1] << "]" << std::endl;
        std::cout << "  clusterMinEnergy: [" << clusterMinEnergy[0] << ", " << clusterMinEnergy[1] << "]" << std::endl;
        std::cout << "  clusterMinMaxAmp: " << clusterMinMaxAmp
                  << ", clusterChi2ReLimit: " << clusterChi2ReLimit
                  << ", clusterMinTierCharge: " << clusterMinTierCharge << std::endl;
    }
};