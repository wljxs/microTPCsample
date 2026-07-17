#pragma once

#ifndef BAD_STRIP_SCRIPT_H
#define BAD_STRIP_SCRIPT_H

#include "Script/Base/IScript.h"
#include <string>
#include <map>
#include <vector>

/**
 * @brief 异常通道噪声分析脚本,改动启轩的ChannelNoiseScript
 * 
 * 分析每个探测器通道的噪声特性，包括：
 * - 异常通道识别
 * 
 * 输出ROOT文件包含：
 * - 探测器的hitmap图
 * - 每个探测器的异常通道分布图
 */
class BadStripScript : public IScript {
public:
    BadStripScript() = default;
    ~BadStripScript() override = default;

    std::string GetName() const override { return "BadStripScript"; }
    
    std::string GetVersion() const override { return "1.0.0"; }
    
    std::string GetDescription() const override { 
        return "Analyze detector bad channel and output to ROOT file"; 
    }

    void LoadConfig(const json& config) override;
    void Print() const override;
    bool Execute() override;

private:
    /**
     * @brief 异常通道噪声统计数据结构
     */
    struct ChannelNoise {
        int detectorID;
        int stripID;
        int type;  // 0=X, 1=Y
        
        double meanBaseline;     // 平均基线
        double rmsBaseline;      // 基线RMS
        double meanNoise;        // 平均噪声（采样点标准差）
        double maxNoise;         // 最大噪声
        int nSamples;            // 采样事例数

        bool isBadChannel;       // 是否为异常通道

        
        std::vector<double> baselineList;  // 保存原始基线数据用于绘制直方图
    };

    /**
     * @brief 分析所有事例的通道噪声
     */
    void AnalyzeChannelNoise();

    /**
     * @brief 保存结果到ROOT文件
     */
    bool SaveResults();

    // 配置参数
    int m_maxEvents;           // 分析的最大事例数（0表示全部）
    int m_progressInterval;    // 进度显示间隔
    int m_detnum;              // 探测器数量
    int m_stripnum;            // 每个探测器的strip数量
    int m_typenum;            // strip类型数量
    double m_badStripThreshold; // 坏道定义：hit占比超过这个值（例如0.1表示10%）
    double m_Threshold;   // 异常通道阈值（ADC）

    // 分析结果
    std::map<int, std::map<int, std::map<int, ChannelNoise>>> m_channelNoiseMap;  
    // [detectorID][type][stripID] -> ChannelNoise
};

#endif // BAD_STRIP_SCRIPT_H