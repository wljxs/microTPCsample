#ifndef WAVEFORM_PROCESSOR_H
#define WAVEFORM_PROCESSOR_H

#include "AlgorithmFactory.h"
#include "IAlgorithm.h"

#include "Config.h"
#include "DataModel.h"

#include <Eigen/Dense>

// 前向声明，避免循环依赖
class DetectorFrame;

/**
 * @brief 波形处理算法 - 负责RawData到StripHit的转换
 *
 * 功能：
 * - 处理原始波形数据（ADC采样点）
 * - 提取关键物理量（峰值、电荷、时间等）
 * - 判断信号有效性
 */
class WaveformProcessor : public IAlgorithm {
   public:
    WaveformProcessor() = default;
    virtual ~WaveformProcessor() = default;

    // 实现IAlgorithm接口
    std::string GetName() const override { return "WaveformProcessor"; }
    std::string GetVersion() const override { return "1.0.0"; }

    void LoadConfig(const json& config) override {
        m_config.loadFrom(config);
    }

    void Print() const override {
        std::cout << "[" << GetName() << " v" << GetVersion() << "]" << std::endl;
        m_config.print();
    }

    // 统一接口: 处理DetectorFrame中的所有RawData
    bool Process(DetectorFrame& frame) override;

    StripHit ProcessWaveform(const RawData& rawData);

    // 设置反卷积矩阵
    void SetDeconvMatrix(const Eigen::MatrixXd& A) { m_deconvMatrix = A; }

   private:
    WaveformConfig m_config;

    // 私有处理方法

    // 无副作用：对每条 strip 使用独立阈值
    StripHit ProcessWaveform(const RawData& rawData, double noiseTh);

    StripHit processWaveformDefault(const RawData& rawData, double noiseTh);
    /*另一种成形函数的得到*/    
    StripHit processWaveformLeadingEdgeFit(const RawData& rawData, double noiseTh);
    /*用garfield的成形函数拟合来得到波形数据的*/
    StripHit processWaveformMode1(const RawData& rawData, double noiseTh);
    /*用deconvolution的方法来处理波形数据，得到更好的时间分辨率*/
    StripHitDeconv processWaveformDeconv(const RawData& rawData, double noiseTh);

    Eigen::MatrixXd m_deconvMatrix;  // 反卷积矩阵，由外部设置
};

#endif  // WAVEFORM_PROCESSOR_H
