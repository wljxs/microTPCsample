#pragma once

#ifndef RAW_WAVEFORM_DUMP_SCRIPT_H
#define RAW_WAVEFORM_DUMP_SCRIPT_H

#include "Script/Base/IScript.h"
#include <string>
#include <vector>

/**
 * @brief 直接对 raw 波形( RawData.adc )做分析与转储的脚本
 *
 * 目标：绕开 StripHit/Cluster 等重建流程，直接把 raw 波形拿出来：
 * - 可选计算：baseline/peak/FWHM
 * - 输出：ROOT 文件 + TTree，包含 adc 向量与关键指标
 */
class RawWaveformDumpScript : public IScript {
public:
    RawWaveformDumpScript() = default;
    ~RawWaveformDumpScript() override = default;

    std::string GetName() const override { return "RawWaveformDumpScript"; }

    std::string GetDescription() const override {
        return "Dump and analyze raw waveforms (adc samples) into a ROOT TTree";
    }

    void LoadConfig(const json& config) override;
    void Print() const override;
    bool Execute() override;

private:
    // 配置
    int m_startEvent = 0;
    int m_maxEvents = 1000;          // <=0 means all
    int m_progressInterval = 1000;

    // 选择过滤（为空则不过滤）
    std::vector<int> m_detIDs;
    std::vector<int> m_stripTypes;
    std::vector<int> m_stripIDs;

    // 波形分析参数
    int m_baselineNSamples = 1;      // baseline = mean(first N)
    double m_noiseThreshold = 0.0;   // optional: peak-baseline < noiseTh => mark fail

    // 是否只输出“坏波形”
    bool m_dumpOnlyFailed = false;

    // 输出
    std::string m_outputFile = "RawWaveforms.root";

private:
    bool passFilter(int detID, int stripType, int stripID) const;
};

#endif // RAW_WAVEFORM_DUMP_SCRIPT_H
