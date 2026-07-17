#pragma once

#include "DeconvCuts.h"
#include "Event/DataModel.h"
#include <Eigen/Dense>
#include <array>
#include <vector>

// 前向声明
class BadStripMask;

// ========== 全局物理参数（定义在 DeconvAnalysisScript.cpp）==========
extern int n_tier;
extern int voltage;
extern double thickness;
extern double z_ref_position;
extern double pitchwidth[2];

// ========== 加权直线拟合（用于 microTPC 位置重建）==========
std::array<double, 3> TrackFitself(
    std::vector<double> z_vec,
    std::vector<double> x_vec,
    std::vector<double> x_error_vec);

// ========== APV25 单条波形特征提取 ==========
// 提取: 半高宽(FWHM)、峰值(amp)、上升时间(riseTime)、过阈时间(overthresholdTime)、达峰时间(peakTime)
// 返回: true 如果波形有效（过阈且特征成功提取）, false 否则
bool AnalyzeWaveformFeatures(
    const std::vector<double>& waveform,
    double noiseTh,
    double t0_ns,
    StripHitDeconv* sh,
    double* peakTime_ns = nullptr);

// vector<short> 重载：内部转为 double 后调用主函数
inline bool AnalyzeWaveformFeatures(
    const std::vector<short>& waveform,
    double noiseTh,
    double t0_ns,
    StripHitDeconv* sh,
    double* peakTime_ns = nullptr)
{
    std::vector<double> wf(waveform.begin(), waveform.end());
    return AnalyzeWaveformFeatures(wf, noiseTh, t0_ns, sh, peakTime_ns);
}

// ========== 反卷积条聚类（将相邻 StripHitDeconv 合并为 ClusterDeconv）==========
// 输入 stripHits 必须已按 ID 升序排列
// hasCluster 累加有效 cluster 数量，outClusters 追加聚类结果
void ClusterStripHits(
    const std::vector<StripHitDeconv>& stripHits,
    const BadStripMask* badStripMask,
    const DeconvCuts& cuts,
    int eventID,
    int type,
    int& hasCluster,
    std::vector<ClusterDeconv>& outClusters);

std::vector<double> FitWaveform(Eigen::MatrixXd& A,StripHitDeconv& sh);

// ========== 单个 Cluster 的残差结果 ==========
struct SingleClusterResidual {
    std::vector<double> tier_residuals;  // 各 tier 位置残差 (mm), 无电荷为 inf
    double microtpc_residual = std::numeric_limits<double>::infinity();
    int n_valid_tiers = 0;               // 有效 tier 数
};

// ========== 计算单个 Cluster 相对于真值的各 tier 残差 ==========
SingleClusterResidual ComputeTierResiduals(
    const ClusterDeconv& cluster,
    double true_k[2], double true_b[2]);