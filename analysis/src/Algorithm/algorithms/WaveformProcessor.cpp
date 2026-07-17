#include "algorithms/WaveformProcessor.h"
#include "AlgorithmFactory.h"
#include "DataModel.h"
#include "DetectorFrame.h"
#include "TVirtualFFT.h"
#include <AnalysisUtils.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>

REGISTER_ALGORITHM("WaveformProcessor", WaveformProcessor)

bool WaveformProcessor::Process(DetectorFrame &frame) {
  const auto &rawData = frame.Raw();
  if (rawData.empty())
    return false;

  const auto *det = &frame.det();
  const auto *badMask = det ? det->GetBadStripMask() : nullptr;

  auto &stripHitsDeconv = frame.GetMutableStripHitsDeconv();
  stripHitsDeconv.clear();
  stripHitsDeconv.reserve(rawData.size());

  auto &stripHits = frame.GetMutableStripHits();
  stripHits.clear();
  stripHits.reserve(rawData.size());

  const auto *pedStdByType = frame.StripPedestalStd();

  // 处理每个RawData，生成StripHit
  for (size_t i = 0; i < rawData.size(); ++i) {
    if (badMask && badMask->IsBad(rawData[i].type, rawData[i].stripID)) {
      //      continue;  // 屏蔽坏道
    }

    // 默认阈值：配置里的 noiseThreshold
    double noiseTh = m_config.noiseThreshold;

    // 如有 pedestals.apv_pedstd，则使用该 strip 的值作为阈值
    if (pedStdByType) {
      auto it = pedStdByType->find(rawData[i].type);
      if (it != pedStdByType->end()) {
        const auto &v = it->second;
        const int idx = rawData[i].stripID - 1; // stripID 从 1 开始
        if (idx >= 0 && idx < static_cast<int>(v.size()) && v[idx] > 0.0f) {
          noiseTh = static_cast<double>(v[idx]) * 3;
        }
      }
    }
    if (m_config.mode == "DeConv") {
      // std::cout << 2 << std::endl;
      StripHitDeconv sh = processWaveformDeconv(rawData[i], noiseTh);
      sh.rawIndices = static_cast<int>(i);
      stripHitsDeconv.push_back(sh);
      continue;
    } else {
      StripHit sh = processWaveformDefault(rawData[i], noiseTh);
      sh.rawIndices = static_cast<int>(i);
      if(sh.isValid) stripHits.push_back(sh);
      continue;
    }
  }

  // 排序：首先按type升序，相同type内按stripID升序
  std::sort(stripHits.begin(), stripHits.end(),
            [](const StripHit &a, const StripHit &b) {
              if (a.type != b.type)
                return a.type < b.type;
              return a.ID < b.ID;
            });

  return true;
}


StripHit WaveformProcessor::processWaveformDefault(const RawData &rawData,
                                                   double noiseTh) {
  const auto &waveform = rawData.adc;
  const size_t nSamples = waveform.size();

  StripHit stripData;
  stripData.isValid = true;
  stripData.ID = rawData.stripID;
  stripData.type = rawData.type;

  if (nSamples == 0) {
    stripData.isValid = false;
    return stripData;
  }

  static const double weights[5] = {0.25, 0.5, 1.0, 0.5, 0.25};
  std::vector<double> smooth(nSamples, 0.0);

  // 是否平滑滤波,目前认为不开，因为信号的点比较靠前，平滑影响基线
  if (false) {
    for (size_t i = 0; i < nSamples; ++i) {
      double sumw = 0;
      double s = 0;
      for (int m = -2; m <= 2; ++m) {
        int idx = static_cast<int>(i) + m;
        if (idx >= 0 && idx < static_cast<int>(nSamples)) {
          double w = weights[m + 2];
          s += w * static_cast<double>(waveform[idx]);
          sumw += w;
        }
      }
      smooth[i] = (sumw > 0) ? (s / sumw) : 0.0;
    }
  } else {
    for (size_t i = 0; i < nSamples; ++i) {
      smooth[i] = waveform[i];
    }
  }

  // 找峰值
  double peakAmpD = -1.0;
  int peakIdx = -1;
  for (size_t i = 0; i < nSamples; ++i) {
    if (smooth[i] > peakAmpD) {
      peakAmpD = smooth[i];
      peakIdx = static_cast<int>(i);
    }
  }

  int peakAmp = static_cast<int>(std::round(peakAmpD));

  if (peakAmp < noiseTh) {
    stripData.isValid = false;
  }
  if (peakIdx <= 0 || peakIdx >= static_cast<int>(nSamples) - 1) {
    stripData.isValid = false;
  }

  // 计算基线!!!!目前认为扣光了
  double baseline = 0.0;
  // size_t baselineSamples = std::min<size_t>(1, nSamples);
  // double bsum = 0;
  // for (size_t i = 0; i < baselineSamples; ++i)
  //   bsum += waveform[i];
  // baseline = bsum / baselineSamples;
  // if(baseline>60 || baseline<-60) {
  //   stripData.isValid = false;
  // }
  // //看峰值是不是在左右两边
  // double halfHeight = (peakAmpD-baseline) *0.9;
  // int leftIdx = 15, rightIdx = -1;
  // for (int i = peakIdx; i >= 0; --i) {
  //   if ((smooth[i] - baseline) < halfHeight) {
  //     leftIdx = i;
  //     break;
  //   }
  // }
  // for (size_t i = peakIdx; i < nSamples; ++i) {
  //   if ((smooth[i] - baseline) < halfHeight) {
  //     rightIdx = static_cast<int>(i);
  //     break;
  //   }
  // }
  // if(leftIdx >= rightIdx || leftIdx < 0 || rightIdx < 0) {
  //   stripData.isValid = false;
  // }

  // 计算电荷//积分
  double inducedCharge = 0.0;
  for (size_t i = 0; i < nSamples; ++i) {
    if (waveform[i] > baseline) {
      inducedCharge += (waveform[i] - baseline);
    }
  }

  // CFD时间提取
  double targetY = baseline + m_config.cfdFraction * (peakAmpD - baseline);
  double minSmooth = *std::min_element(smooth.begin(), smooth.end());
  double maxSmooth = *std::max_element(smooth.begin(), smooth.end());

  if (targetY < minSmooth || targetY > maxSmooth) {
    stripData.isValid = false;
  }

  // 找穿越点
  int crossingIdx = -1;
  for (int i = 1; i <= peakIdx; ++i) {
    if (smooth[i] >= targetY && smooth[i - 1] < targetY) {
      crossingIdx = i;
      break;
    }
  }

  if (crossingIdx == -1) {
    stripData.isValid = false;
  }

  // 线性插值
  double fitTime = 0.0;
  if (stripData.isValid) {
    double y1 = smooth[crossingIdx - 1];
    double y2 = smooth[crossingIdx];
    double frac = (targetY - y1) / (y2 - y1);
    fitTime = (crossingIdx - 1) + frac;
  }

  // 填充结果
  stripData.amp = peakAmp - baseline;
  stripData.charge = inducedCharge;
  stripData.peakTime = peakIdx * m_config.timePitch;
  stripData.time = fitTime * m_config.timePitch;
  stripData.baseline = baseline;
  stripData.riseTime = 0;
  stripData.timeError = 0;
  stripData.isSaturated = (peakAmp > m_config.saturationLevel);

  return stripData;
}

StripHitDeconv WaveformProcessor::processWaveformDeconv(const RawData &rawData,
                                                        double noiseTh) {

  const auto &waveform = rawData.adc;
  const size_t nSamples = waveform.size();
  // std::cout << nSamples << std::endl;
  StripHitDeconv stripData;
  bool isValid = true;
  stripData.ID = rawData.stripID;
  stripData.type = rawData.type;

  if (nSamples == 0) {
    isValid = false;
    return stripData;
  }

  static const double weights[5] = {0.25, 0.5, 1.0, 0.5, 0.25};
  std::vector<double> smooth(nSamples, 0.0);

  for (int i = 0; i < nSamples; ++i) {
    smooth[i] = waveform[i];
  }

  // 找峰值
  double peakAmpD = -1.0;
  int peakIdx = -1;
  for (size_t i = 0; i < nSamples; ++i) {
    if (smooth[i] > peakAmpD) {
      peakAmpD = smooth[i];
      peakIdx = static_cast<int>(i);
    }
  }

  int peakAmp = static_cast<int>(std::round(peakAmpD));
  if (peakAmp < noiseTh) {
    isValid = false;
  }
  if (peakIdx <= 0 || peakIdx >= static_cast<int>(nSamples) - 1) {
    isValid = false;
  }

  // 计算基线!!!!目前认为扣光了
  double baseline = 0.0;
  // //看峰值是不是在左右两边
  double halfHeight = (peakAmpD-baseline) *0.9;
  int leftIdx = 15, rightIdx = -1;
  for (int i = peakIdx; i >= 0; --i) {
    if ((smooth[i] - baseline) < halfHeight) {
      leftIdx = i;
      break;
    }
  }
  for (size_t i = peakIdx; i < nSamples; ++i) {
    if ((smooth[i] - baseline) < halfHeight) {
      rightIdx = static_cast<int>(i);
      break;
    }
  }
  if(leftIdx >= rightIdx || leftIdx < 0 || rightIdx < 0) {
    isValid = false;
  }

  // 计算电荷//积分
  double inducedCharge = 0.0;
  for (size_t i = 0; i < nSamples; ++i) {
    if (waveform[i] > baseline) {
      inducedCharge += (waveform[i] - baseline);
    }
  }

  // CFD时间提取
  double targetY = baseline + m_config.cfdFraction * (peakAmpD - baseline);
  double minSmooth = *std::min_element(smooth.begin(), smooth.end());
  double maxSmooth = *std::max_element(smooth.begin(), smooth.end());

  if (targetY < minSmooth || targetY > maxSmooth) {
    isValid = false;
  }

  // 找穿越点
  int crossingIdx = -1;
  for (int i = 1; i <= peakIdx; ++i) {
    if (smooth[i] >= targetY && smooth[i - 1] < targetY) {
      crossingIdx = i;
      break;
    }
  }

  if (crossingIdx == -1) {
    isValid = false;
  }

  // 线性插值
  double fitTime = 0.0;
  if (isValid) {
    double y1 = smooth[crossingIdx - 1];
    double y2 = smooth[crossingIdx];
    double frac = (targetY - y1) / (y2 - y1);
    fitTime = (crossingIdx - 1) + frac;
  }

  if (isValid) {
    std::vector<double> b_vec(waveform.begin(), waveform.end());
    auto [x_vec, rnorm] = AnalysisUtils::SolveNNLSLawHanson(b_vec, m_deconvMatrix);
    stripData.chargetier = x_vec;
    stripData.rnorm = rnorm;
  }

  return stripData;
}
