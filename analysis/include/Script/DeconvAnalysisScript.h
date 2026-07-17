#pragma once

#ifndef DECONV_ANALYSIS_SCRIPT_H
#define DECONV_ANALYSIS_SCRIPT_H

#include "Script/Base/IScript.h"
#include "Algorithm/DeconvCuts.h"
#include "Event/DataModel.h"
#include <string>


/**
 * @brief 示例脚本
 * 
 * 演示如何创建自定义分析脚本
 * 展示基本的数据访问和输出功能
 */
class DeconvAnalysisScript : public IScript {
public:
    DeconvAnalysisScript() = default;
    ~DeconvAnalysisScript() override = default;

    std::string GetName() const override { return "DeconvAnalysisScript"; }
    
    std::string GetDescription() const override { 
        return "DeconvAnalysis script demonstrating the script framework"; 
    }

    void LoadConfig(const json& config) override;
    void Print() const override;
    bool Execute() override;

private:
    bool RunDeconvolution();  // Phase 1: 反卷积+聚类+画图+平均波形+坏道检测
    bool RunAlignment();      // Phase 2: 从DeConvInfo.root读取数据 → 对齐 → 位置分辨

    std::string m_message;
    int m_eventLimit;

    // 运行模式: "deconv" / "align" / "full" (默认)
    std::string m_mode = "full";
    // 是否使用 delay.root 的逐通道 t0 缓存
    bool m_usePerChannelT0 = true;
    // 是否在反卷积阶段检测坏道
    bool m_detectBadStrips = true;

    // 可配置的 cut 参数
    DeconvCuts m_cuts;
};

#endif // DECONV_ANALYSIS_SCRIPT_H
