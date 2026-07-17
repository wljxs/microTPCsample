#include "WaveformAnalysis.h"
#include "BadStripMask.h"
#include "AnalysisUtils.h"
#include "DataModel.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>

// ========== 加权直线拟合（用于 microTPC 位置重建）==========
std::array<double, 3> TrackFitself(
    std::vector<double> z_vec1,
    std::vector<double> x_vec1,
    std::vector<double> x_error_vec1)
{
  const double nanv = std::numeric_limits<double>::quiet_NaN();
  if (z_vec1.empty() || x_vec1.empty() || x_error_vec1.empty() ||
      z_vec1.size() != x_vec1.size() || x_vec1.size() != x_error_vec1.size()) {
    return std::array<double, 3>{nanv, nanv, nanv};
  }
  std::vector<double> weights(z_vec1.size());
  double weight_sum = 0.0;
  for (size_t i = 0; i < z_vec1.size(); ++i) {
    if (!std::isfinite(x_error_vec1[i]) || x_error_vec1[i] <= 0.0) {
      std::cout << "Warning: Invalid error value for point " << i
                << ", setting weight to 0. Error: " << x_error_vec1[i] << std::endl;
      return std::array<double, 3>{nanv, nanv, nanv};
    }
    weights[i] = 1.0 / (x_error_vec1[i] * x_error_vec1[i]);
    weight_sum += weights[i];
  }
  if (weight_sum <= 0.0 || !std::isfinite(weight_sum)) {
    std::cout << "Warning: Non-positive or non-finite total weight, cannot perform fit. "
              << "Weight sum: " << weight_sum << std::endl;
    return std::array<double, 3>{nanv, nanv, nanv};
  }
  double x_mean = 0.0;
  double z_mean = 0.0;
  int n = static_cast<int>(x_vec1.size());
  for (int i = 0; i < n; i++) {
    x_mean += x_vec1[i] * weights[i] / weight_sum;
    z_mean += z_vec1[i] * weights[i] / weight_sum;
  }
  double lxx = 0.0, lzz = 0.0, lxz = 0.0;
  for (int i = 0; i < n; i++) {
    lxx += weights[i] * (x_vec1[i] - x_mean) * (x_vec1[i] - x_mean);
    lzz += weights[i] * (z_vec1[i] - z_mean) * (z_vec1[i] - z_mean);
    lxz += weights[i] * (x_vec1[i] - x_mean) * (z_vec1[i] - z_mean);
  }
  if (!std::isfinite(lzz) || lzz == 0.0) {
    std::cout << "Warning: Invalid variance in z, cannot perform fit. Lzz: " << lzz << std::endl;
    return std::array<double, 3>{nanv, nanv, nanv};
  }
  double slope = lxz / lzz;
  double intercept = x_mean - slope * z_mean;
  double chi2 = 0.0;
  for (int i = 0; i < n; i++) {
    double x_fit = slope * z_vec1[i] + intercept;
    chi2 += weights[i] * (x_vec1[i] - x_fit) * (x_vec1[i] - x_fit);
  }
  return std::array<double, 3>{chi2, slope, intercept};
}

// ========== APV25 单条波形特征提取 ==========
bool AnalyzeWaveformFeatures(
    const std::vector<double>& waveform,
    double noiseTh,
    double t0_ns,
    StripHitDeconv* sh,
    double* peakTime_ns)
{
  const size_t nSamples = waveform.size();
  if (nSamples == 0) return false;

  // 1. 复制波形
  std::vector<double> smooth(nSamples, 0.0);
  for (size_t j = 0; j < nSamples; ++j) {
    smooth[j] = waveform[j];
  }

  // 2. 找峰值
  double peakAmpD = -1.0;
  int peakIdx = -1;
  for (size_t j = 0; j < nSamples; ++j) {
    if (smooth[j] > peakAmpD) {
      peakAmpD = smooth[j];
      peakIdx = static_cast<int>(j);
    }
  }
  // 3. 阈值与边界检查
  int peakAmp = static_cast<int>(std::round(peakAmpD));
  if (peakAmp < noiseTh) return false;
  if (peakIdx <= 0 || peakIdx >= static_cast<int>(nSamples) - 1 || peakIdx > 14) return false;

  // 4. 达峰时间（可选输出）
  if (peakTime_ns) *peakTime_ns = peakIdx * 25.0;

  // 5. 半高宽 (FWHM)
  const double baseline = 0.0;
  const double halfHeight = (peakAmpD - baseline) * 0.5;
  int leftIdx = 30, rightIdx = -1;
  for (int j = peakIdx; j >= 0; --j) {
    if ((smooth[j] - baseline) < halfHeight) { leftIdx = j; break; }
  }
  for (size_t j = peakIdx; j < nSamples; ++j) {
    if ((smooth[j] - baseline) < halfHeight) { rightIdx = static_cast<int>(j); break; }
  }
  if (leftIdx >= rightIdx || leftIdx < 0 || rightIdx < 0 || (rightIdx-leftIdx)>11) return false;

  // 6. 上升时间 (10% → 90%)
  double risetimeindex1 = 30, risetimeindex2 = 0;
  for (int j = peakIdx; j >= 0; --j) {
    if ((smooth[j] - baseline) < (peakAmpD - baseline) * 0.1) {
      risetimeindex1 = j; break;
    }
  }
  for (int j = peakIdx; j >= 0; --j) {
    if ((smooth[j] - baseline) < (peakAmpD - baseline) * 0.9) {
      risetimeindex2 = j; break;
    }
  }
  if (risetimeindex1 > risetimeindex2 || (risetimeindex2-risetimeindex1)>7) return false;

  // 7. 过阈时间 (50% 前沿穿越 + 线性插值)
  const double targetY = baseline + (peakAmpD - baseline) * 0.5;
  int crossingIdx = -1;
  for (int j = 1; j <= peakIdx; ++j) {
    if (smooth[j - 1] < targetY && smooth[j] >= targetY) {
      crossingIdx = j; break;
    }
  }
  if (crossingIdx == -1) return false;

  const double frac = (targetY - smooth[crossingIdx - 1])
          / (smooth[crossingIdx] - smooth[crossingIdx - 1])
          + (crossingIdx - 1);

  // 8. 积分电荷
  double integralCharge = 0.0;
  for (size_t j = 0; j < nSamples; ++j) {
    integralCharge += smooth[j] - baseline;
  }
  // 9. 填充结果（时间单位：ns，采样间隔 25ns）
  sh->amp               = peakAmpD;
  sh->peakTime          = peakIdx;
  sh->Integral          = integralCharge;
  sh->FwHM              = (rightIdx - leftIdx) * 25.0;
  sh->riseTime          = (risetimeindex2 - risetimeindex1) * 25.0;
  sh->overthresholdTime = frac * 25.0;
  sh->t0                = t0_ns;

  return true;
}

// ========== 反卷积条聚类（将相邻 StripHitDeconv 合并为 ClusterDeconv）==========
void ClusterStripHits(
    const std::vector<StripHitDeconv>& stripHits,
    const BadStripMask* badStripMask,
    const DeconvCuts& cuts,
    int eventID,
    int type,
    int& hasCluster,
    std::vector<ClusterDeconv>& outClusters)
{
  if (stripHits.empty()) return;

  // ---- 辅助 lambda：计算 cluster 的 microTPC 位置和质量 ----
  auto FinalizeCluster = [&](ClusterDeconv& cluster) -> bool {
    std::vector<double> x_vec;
    std::vector<double> x_error_vec;
    std::vector<double> z_vec;
    double x_mean = 0.0;
    double z_mean = 0.0;
    double weight_sum = 0.0;
    cluster.tier_hascharge = 0;

    for (int k = 0; k < cuts.n_tier; k++) {
      if (cluster.chargetier[k] == 0) {
        cluster.pos[k] = std::numeric_limits<double>::infinity();
        continue;
      }
      cluster.pos[k] /= cluster.chargetier[k];
      if (!std::isfinite(cluster.pos[k])) {
        std::cout << "liu:" << cluster.pos[k] << "," << cluster.chargetier[k] << std::endl;
      }
      x_vec.push_back(cluster.pos[k]);
      z_vec.push_back(cuts.thickness / cuts.n_tier * (k + 0.5));
      x_error_vec.push_back(std::sqrt(( k + 0.5) / cluster.chargetier[k]));  // 0.4mm 是 strip pitch
      x_mean += cluster.pos[k] * cluster.chargetier[k];
      z_mean += (cuts.thickness / cuts.n_tier * (k + 0.5)) * cluster.chargetier[k];
      weight_sum += cluster.chargetier[k];
      cluster.tier_hascharge++;
    }

    if (cluster.tier_hascharge == 0) {
      std::cout << "Warning: Cluster has no charge tiers in event " << eventID << std::endl;
      return false;
    }
    if (cluster.tier_hascharge > 1) {
      auto result = TrackFitself(z_vec, x_vec, x_error_vec);
      cluster.chi2 = result[0];
      cluster.k = result[1];
      cluster.b = result[2];
      cluster.microTPCposition = result[2] + result[1] * cuts.z_ref_position;
      if (cluster.tier_hascharge > 2) {
        cluster.chi2_re = cluster.chi2 / (cluster.tier_hascharge - 2);
      }
    } else {
      cluster.microTPCposition = x_vec[0];
    }

    if (cluster.microTPCposition < 0 || cluster.microTPCposition > 256 ||
        !std::isfinite(cluster.microTPCposition)) {
      std::cout << "Warning: Cluster microTPC position out of bounds in event "
                << eventID << std::endl;
      for (size_t k = 0; k < x_vec.size(); k++) {
        std::cout << "Tier " << k << ": pos = " << x_vec[k]
                  << ", error = " << x_error_vec[k] << ", z = " << z_vec[k] << std::endl;
      }
      std::cout << "Cluster microTPC position: " << cluster.microTPCposition << std::endl;
      return false;
    }

    if (cluster.energy == 0) {
      std::cout << "Warning: Cluster energy is zero in event " << eventID << std::endl;
      cluster.chargeposition = 0.0;
    } else {
      cluster.chargeposition /= cluster.energy;
    }

    if (cluster.energy > cuts.clusterMinEnergy[cluster.type] &&
        cluster.maxAmp > cuts.clusterMinMaxAmp && !cluster.isBad &&
        cluster.chi2_re < cuts.clusterChi2ReLimit &&
        cluster.tier_hascharge >= cuts.clusterMinTierCharge) {
      cluster.isValid = true;
      hasCluster++;
    }
    outClusters.push_back(cluster);
    return true;
  };

  // ---- 初始化第一个 cluster ----
  ClusterDeconv cluster;
  cluster.isValid = false;
  cluster.isBad = stripHits[0].isBad;
  cluster.chargeposition = 0.0;
  cluster.chi2 = std::numeric_limits<double>::infinity();
  cluster.chi2_re = std::numeric_limits<double>::infinity();
  cluster.k = std::numeric_limits<double>::infinity();
  cluster.b = std::numeric_limits<double>::infinity();
  cluster.chargetier.resize(cuts.n_tier);
  cluster.pos.resize(cuts.n_tier);
  cluster.type = type;
  for (int k = 0; k < cuts.n_tier; ++k) {
    cluster.pos[k] = stripHits[0].ID * stripHits[0].chargetier[k];
    cluster.chargetier[k] = stripHits[0].chargetier[k];
  }
  cluster.energy = stripHits[0].amp;
  cluster.maxAmp = stripHits[0].amp;
  cluster.chargeposition = stripHits[0].ID * stripHits[0].amp;
  cluster.stripHitIndices.push_back(0);
  cluster.stripIDs.push_back(stripHits[0].ID);
  cluster.size = 1;

  // ---- 遍历剩余条 ----
  for (size_t m = 1; m < stripHits.size(); m++) {
    if (stripHits[m].isBad) {
      cluster.isBad = true;
    }

    bool shouldMerge = (stripHits[m].ID <= stripHits[m - 1].ID + 2);
    if (!shouldMerge) {
      // 间隙不连续，检查是否因为中间全是坏道导致
      shouldMerge = true;
      for (int stripid = stripHits[m - 1].ID + 1; stripid < stripHits[m].ID; stripid++) {
        if (!(badStripMask && badStripMask->IsBad(cluster.type, stripid))) {
          shouldMerge = false;  // 间隙中有好道，不能合并
          break;
        }
      }
      if (shouldMerge) {
        cluster.isBad = true;  // 间隙全是坏道，跨过合并，标记为坏cluster
      }
    }

    if (shouldMerge) {
      // 合并到当前 cluster
      cluster.stripHitIndices.push_back(static_cast<int>(m));
      cluster.stripIDs.push_back(stripHits[m].ID);
      for (int k = 0; k < cuts.n_tier; ++k) {
        cluster.pos[k] += stripHits[m].ID * stripHits[m].chargetier[k];
        cluster.chargetier[k] += stripHits[m].chargetier[k];
      }
      cluster.energy += stripHits[m].amp;
      cluster.chargeposition += stripHits[m].ID * stripHits[m].amp;
      if (stripHits[m].amp > cluster.maxAmp)
        cluster.maxAmp = stripHits[m].amp;
      cluster.size++;
    } else {
      // 当前 cluster 结束，处理并保存
      FinalizeCluster(cluster);

      // 开始新的 cluster
      cluster.stripHitIndices.clear();
      cluster.stripIDs.clear();
      cluster.isValid = false;
      cluster.isBad = stripHits[m].isBad;
      cluster.chargeposition = 0.0;
      cluster.chi2 = std::numeric_limits<double>::infinity();
      cluster.chi2_re = std::numeric_limits<double>::infinity();
      cluster.k = std::numeric_limits<double>::infinity();
      cluster.b = std::numeric_limits<double>::infinity();
      cluster.type = type;
      for (int k = 0; k < cuts.n_tier; ++k) {
        cluster.pos[k] = stripHits[m].ID * stripHits[m].chargetier[k];
        cluster.chargetier[k] = stripHits[m].chargetier[k];
      }
      cluster.chargeposition = stripHits[m].ID * stripHits[m].amp;
      cluster.energy = stripHits[m].amp;
      cluster.maxAmp = stripHits[m].amp;
      cluster.stripHitIndices.push_back(static_cast<int>(m));
      cluster.stripIDs.push_back(stripHits[m].ID);
      cluster.size = 1;
    }
  }

  // ---- 最后一个 cluster ----
  FinalizeCluster(cluster);
}

std::vector<double> FitWaveform(Eigen::MatrixXd& A,StripHitDeconv& sh){
  Eigen::VectorXd b(A.cols());
  for(int i = 0;i<A.cols();i++){
    b(i) = sh.chargetier[i];
  }
  Eigen::VectorXd x = A * b;
  std::vector<double> result(x.data(), x.data() + x.size());
  return result;
}

// ========== 计算单个 Cluster 相对于真值的各 tier 残差 ==========
SingleClusterResidual ComputeTierResiduals(
    const ClusterDeconv& cluster,
    double true_k[2], double true_b[2])
{
  SingleClusterResidual result;
  const int type = cluster.type;
  result.tier_residuals.resize(n_tier);

  for (int tier = 0; tier < n_tier; ++tier) {
    if (cluster.pos[tier] == std::numeric_limits<double>::infinity()) {
      result.tier_residuals[tier] = std::numeric_limits<double>::infinity();
      continue;
    }
    // 重建位置 (mm) = strip编号 × 0.4mm + 0.2mm偏移
    double reco = cluster.pos[tier] * 0.4 + 0.2;
    // 真值位置 (mm)
    double z = thickness / n_tier * (tier + 0.5);
    double truth = true_b[type] * 10.0 + true_k[type] * z;
    result.tier_residuals[tier] = reco - truth;
    result.n_valid_tiers++;
  }

  // microTPC 位置残差
  if (std::isfinite(cluster.microTPCposition)) {
    double reco = cluster.microTPCposition * 0.4 + 0.2;
    double truth = true_b[type] * 10.0 + true_k[type] * thickness / 2.0;
    result.microtpc_residual = reco - truth;
  }

  return result;
}