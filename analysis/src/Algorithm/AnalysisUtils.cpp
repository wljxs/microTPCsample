#include "Algorithm/AnalysisUtils.h"
#include "Event/DetectorFrame.h"
#include "TTree.h"

#include <TFile.h>
#include <TH1D.h>
#include <TVirtualFFT.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <ostream>

using namespace std;


namespace AnalysisUtils {

Track FitTrack(const vector<TVector3> &hits) {
  Track t{};
  size_t n = hits.size();
  if (n < 2)
    return t;

  double sum_x = 0, sum_y = 0, sum_z = 0;
  for (auto &p : hits) {
    sum_x += p.X();
    sum_y += p.Y();
    sum_z += p.Z();
  }
  double mean_x = sum_x / n, mean_y = sum_y / n, mean_z = sum_z / n;

  double szz = 0, szx = 0, szy = 0;
  for (auto &p : hits) {
    double dz = p.Z() - mean_z;
    szz += dz * dz;
    szx += dz * (p.X() - mean_x);
    szy += dz * (p.Y() - mean_y);
  }

  if (szz == 0.0) {
    std::cerr << "[AnalysisUtils] FitTrack: all hits have the same z coordinate, cannot fit." << std::endl;
    return t;
  }
  t.kx = szx / szz;
  t.ky = szy / szz;
  t.bx = mean_x - t.kx * mean_z;
  t.by = mean_y - t.ky * mean_z;

  double chi2 = 0;
  for (auto &p : hits) {
    double dx = p.X() - (t.kx * p.Z() + t.bx);
    double dy = p.Y() - (t.ky * p.Z() + t.by);
    chi2 += dx * dx + dy * dy;
  }
  t.chi2 = chi2 / (2 * n - 4);
  return t;
}

std::pair<double, double> GetRange(const std::vector<double> &v) {
  if (v.size() < 3)
    return std::make_pair(0.0, 1.0);
  const double k = 5;

  // -----------------------------
  // 1st pass: raw mean / sigma
  // -----------------------------
  double sum1 = 0, sq1 = 0;
  for (double x : v)
    sum1 += x;
  double mean1 = sum1 / v.size();

  for (double x : v)
    sq1 += (x - mean1) * (x - mean1);
  double sigma1 = std::sqrt(sq1 / v.size());

  double low1 = mean1 - k * sigma1;
  double high1 = mean1 + k * sigma1;

  double sum2 = 0, sq2 = 0;
  int n2 = 0;

  for (double x : v) {
    if (x >= low1 && x <= high1) {
      sum2 += x;
      n2++;
    }
  }

  double mean2 = (n2 > 0) ? sum2 / n2 : mean1;
  if (n2 == 0) return std::make_pair(mean1, mean1);

  for (double x : v) {
    if (x >= low1 && x <= high1)
      sq2 += (x - mean2) * (x - mean2);
  }

  double sigma2 = std::sqrt(sq2 / n2);

  return std::make_pair(mean2 - k * sigma2, mean2 + k * sigma2);
}

void FFTAnalyzer(Cluster &cluster, Event &evt, int det_id) {

  // ========= 静态累积对象 =========
  static TH1D *h_spectrum = nullptr;
  static int n_accumulated = 0;
  constexpr int n_samples = 26;

  // ---------- 1. 初始化直方图（只做一次） ----------
  if (!h_spectrum) {
    const double dt = 25;                      // ns
    const double df = 1.0 / (n_samples * dt);  // GHz
    const int n_freq = n_samples / 2 + 1;

    h_spectrum =
        new TH1D("h_spectrum", "Accumulated Power Spectrum;Frequency;Power",
                 n_freq, 0, n_freq * df);
  }

  for (int strip_index : cluster.stripHitIndices) {
    StripHit hit = evt.detectorFramesMap[det_id]->GetStripHit(strip_index);
    auto raw_data = evt.detectorFramesMap[det_id]->GetRawFromStrip(hit);
    const auto &waveform = raw_data->adc;

    if (waveform.size() < static_cast<size_t>(n_samples)) continue;

    std::vector<double> time_signal(n_samples);
    for (int i = 0; i < n_samples; ++i)
      time_signal[i] = static_cast<double>(waveform[i]);

    int n_fft = n_samples;
    TVirtualFFT *fft = TVirtualFFT::FFT(1, &n_fft, "R2C M K");
    fft->SetPoints(time_signal.data());
    fft->Transform();

    const int n_freq = n_samples / 2 + 1;

    // ---------- 3. 累积功率谱 ----------
    double re = 0.0, im = 0.0;
    for (int k = 0; k < n_freq; ++k) {
      fft->GetPointComplex(k, re, im);

      double power = re * re + im * im;

      // 单边谱修正
      if (k != 0 && k != n_samples / 2)
        power *= 2.0;

      power /= (n_samples * n_samples);

      h_spectrum->AddBinContent(k + 1, power);
    }

    delete fft;
    ++n_accumulated;
  }

  // ---------- 2. FFT ----------

  if (n_accumulated % 10000 == 0) {
    h_spectrum->Scale(1.0 / n_accumulated);

    TFile f_out("spectrum.root", "RECREATE");
    h_spectrum->Write();
    f_out.Close();

    std::cout << "[Mode1] Power spectrum written, N = " << n_accumulated
              << std::endl;
  }
}

double CalculateMean(const std::vector<double> &values) {
  if (values.empty())
    return 0.0;

  double sum = 0.0;
  for (double val : values) {
    sum += val;
  }
  return sum / values.size();
}

double CalculateRMS(const std::vector<double> &values) {
  if (values.empty())
    return 0.0;

  double mean = CalculateMean(values);
  double sum_sq = 0.0;
  for (double val : values) {
    sum_sq += (val - mean) * (val - mean);
  }
  return std::sqrt(sum_sq / values.size());
}

Eigen::MatrixXd DeconvBaseMatrix(std::vector<double> sigzero_sample,int n_sig_track_sample, int n_tier) {
    if (sigzero_sample.empty()) {
      std::cerr << "[AnalysisUtils] sigzero_sample is empty." << std::endl;
      return {};
    }
    if (n_tier <= 0) {
      std::cerr << "[AnalysisUtils] invalid n_tier for deconvolution base matrix."
                << std::endl;
      return {};
    }

  int n_sample_signal = sigzero_sample.size();
  Eigen::MatrixXd A(n_sig_track_sample, n_tier);
  A.setZero();

  for (int i = 0; i < n_tier; i++) {
    for (int k = i; k < n_sig_track_sample; ++k) {
      // Use zero-based template index to avoid shifting the response.
      const int idx = k - i;
      if (idx >= n_sample_signal) {
        break;
      }
      A(k, i) = sigzero_sample.at(idx);
    }
  }
  return A;
}

Eigen::MatrixXd DeconvBaseMatrix(std::vector<std::vector<double>> sigzero_samples,int n_sig_track_sample){
  if(sigzero_samples.empty()){
      std::cerr << "[AnalysisUtils] sigzero_samples is empty." << std::endl;
      return {};
  }
  int n_tier = sigzero_samples.size();
  std::cout << "DeconvBaseMatrix has " << n_sig_track_sample << " raw " << n_tier << " column" <<std::endl;
  Eigen::MatrixXd A(n_sig_track_sample,n_tier);
  A.setZero();

  for(int i=0;i<n_sig_track_sample;i++){
    for(int j=0;j<n_tier;j++){
      A(i,j) = sigzero_samples[j][i];
    }
  }
  return A;
}

std::pair<std::vector<double>, double> SolveNNLSLawHanson(
    std::vector<double> b_vec, const Eigen::MatrixXd &A)
{
  Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(b_vec.data(), b_vec.size());
  Eigen::VectorXd x;
  double rnorm;
  const int m = static_cast<int>(A.rows());
  const int n = static_cast<int>(A.cols());
  if (n == 0) {
    x.resize(0);
    rnorm = b.norm();
    return {{}, rnorm};
  }

  const double tol = 1e-12;
  x = Eigen::VectorXd::Zero(n);

  std::vector<int> passive;
  std::vector<int> active;
  active.reserve(n);
  for (int i = 0; i < n; ++i) {
    active.push_back(i);
  }

  auto remove_index = [](std::vector<int> &v, int idx) {
    v.erase(std::remove(v.begin(), v.end(), idx), v.end());
  };

  Eigen::VectorXd w = A.transpose() * (b - A * x);
  while (true) {
    int t = -1;
    double wmax = tol;
    for (int idx : active) {
      if (w(idx) > wmax) {
        wmax = w(idx);
        t = idx;
      }
    }
    if (t < 0) {
      break;
    }

    remove_index(active, t);
    passive.push_back(t);

    while (true) {
      if (passive.empty()) {
        break;
      }

      Eigen::MatrixXd A_p(m, static_cast<int>(passive.size()));
      for (size_t k = 0; k < passive.size(); ++k) {
        A_p.col(static_cast<int>(k)) = A.col(passive[k]);
      }

      Eigen::VectorXd z = A_p.colPivHouseholderQr().solve(b);
      Eigen::VectorXd x_new = Eigen::VectorXd::Zero(n);
      for (size_t k = 0; k < passive.size(); ++k) {
        x_new(passive[k]) = z(static_cast<int>(k));
      }

      bool all_positive = true;
      for (size_t k = 0; k < passive.size(); ++k) {
        if (x_new(passive[k]) <= 0.0) {
          all_positive = false;
          break;
        }
      }

      if (all_positive) {
        x = x_new;
        break;
      }

      double alpha = 1.0;
      for (size_t k = 0; k < passive.size(); ++k) {
        const int idx = passive[k];
        if (x_new(idx) <= 0.0) {
          const double denom = x(idx) - x_new(idx);
          if (denom > 0.0) {
            alpha = std::min(alpha, x(idx) / denom);
          }
        }
      }

      x = x + alpha * (x_new - x);

      std::vector<int> to_move;
      for (int idx : passive) {
        if (x(idx) <= tol) {
          x(idx) = 0.0;
          to_move.push_back(idx);
        }
      }
      for (int idx : to_move) {
        remove_index(passive, idx);
        active.push_back(idx);
      }
    }

    w = A.transpose() * (b - A * x);
  }

  rnorm = (b - A * x).norm();

  std::vector<double> x_vec(x.data(), x.data() + x.size()); 
  return {x_vec, rnorm};
}

std::pair<std::vector<double>, double> SolveNNLSTikhonov(
    std::vector<double> b_vec, const Eigen::MatrixXd &A, double lambda)
{
  const int m = static_cast<int>(A.rows());  // 采样点数
  const int n = static_cast<int>(A.cols());  // tier 数
  if (n <= 1 || lambda <= 0.0) {
    // 退化: tier 数不足或正则化为零 → 直接调用原始 NNLS
    return SolveNNLSLawHanson(std::move(b_vec), A);
  }

  // ---- 构造一阶差分矩阵 L : (n-1) × n ----
  // L = [ -1  1  0 ...  0 ]
  //     [  0 -1  1 ...  0 ]
  //     [  .............. ]
  //     [  0 ...  0 -1  1 ]
  Eigen::MatrixXd L(n - 1, n);
  L.setZero();
  for (int i = 0; i < n - 1; ++i) {
    L(i, i)     = -1.0;
    L(i, i + 1) =  1.0;
  }

  const double sqrt_lambda = std::sqrt(lambda);

  // ---- 构造增广矩阵 A_aug = [A; sqrt(lambda) * L] ----
  Eigen::MatrixXd A_aug(m + n - 1, n);
  A_aug.topRows(m) = A;
  A_aug.bottomRows(n - 1) = sqrt_lambda * L;

  // ---- 构造增广向量 b_aug = [b; 0] ----
  std::vector<double> b_aug(m + n - 1, 0.0);
  std::copy(b_vec.begin(), b_vec.end(), b_aug.begin());

  return SolveNNLSLawHanson(std::move(b_aug), A_aug);
}

std::pair<std::vector<double>, double> SolveNNLSLasso(
    std::vector<double> b_vec, const Eigen::MatrixXd &A, double lambda)
{
  const int m = static_cast<int>(A.rows());
  const int n = static_cast<int>(A.cols());
  if (n == 0) {
    double rnorm = Eigen::Map<Eigen::VectorXd>(b_vec.data(), b_vec.size()).norm();
    return {{}, rnorm};
  }
  if (lambda <= 0.0) {
    return SolveNNLSLawHanson(std::move(b_vec), A);
  }

  Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(b_vec.data(), b_vec.size());
  const double tol = 1e-12;
  const double lambda_half = lambda / 2.0;  // KKT 条件中 λ/2

  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);

  std::vector<int> passive;
  std::vector<int> active;
  active.reserve(n);
  for (int i = 0; i < n; ++i) {
    active.push_back(i);
  }

  auto remove_index = [](std::vector<int> &v, int idx) {
    v.erase(std::remove(v.begin(), v.end(), idx), v.end());
  };

  // 对偶变量 w = A^T(b - Ax)，但判断时用 w_j - λ/2 > tol
  Eigen::VectorXd w = A.transpose() * (b - A * x);

  while (true) {
    int t = -1;
    double wmax = tol;
    for (int idx : active) {
      if (w(idx) - lambda_half > wmax) {
        wmax = w(idx) - lambda_half;
        t = idx;
      }
    }
    if (t < 0) {
      break;  // 最优解
    }

    remove_index(active, t);
    passive.push_back(t);

    while (true) {
      if (passive.empty()) {
        break;
      }

      const int np = static_cast<int>(passive.size());
      Eigen::MatrixXd A_p(m, np);
      for (int k = 0; k < np; ++k) {
        A_p.col(k) = A.col(passive[k]);
      }

      // 求解 A_P^T A_P z = A_P^T b - (λ/2)·1_P
      Eigen::VectorXd rhs = A_p.transpose() * b;
      rhs.array() -= lambda_half;
      Eigen::MatrixXd AtA = A_p.transpose() * A_p;
      Eigen::VectorXd z = AtA.ldlt().solve(rhs);

      Eigen::VectorXd x_new = Eigen::VectorXd::Zero(n);
      for (int k = 0; k < np; ++k) {
        x_new(passive[k]) = z(k);
      }

      bool all_positive = true;
      for (int k = 0; k < np; ++k) {
        if (x_new(passive[k]) <= 0.0) {
          all_positive = false;
          break;
        }
      }

      if (all_positive) {
        x = x_new;
        break;
      }

      // 计算 α，使 x + α*(x_new - x) 刚好碰到非负边界
      double alpha = 1.0;
      for (int k = 0; k < np; ++k) {
        const int idx = passive[k];
        if (x_new(idx) <= 0.0) {
          const double denom = x(idx) - x_new(idx);
          if (denom > 0.0) {
            alpha = std::min(alpha, x(idx) / denom);
          }
        }
      }

      x = x + alpha * (x_new - x);

      std::vector<int> to_move;
      for (int idx : passive) {
        if (x(idx) <= tol) {
          x(idx) = 0.0;
          to_move.push_back(idx);
        }
      }
      for (int idx : to_move) {
        remove_index(passive, idx);
        active.push_back(idx);
      }
    }

    w = A.transpose() * (b - A * x);
  }

  double rnorm = (b - A * x).norm();

  std::vector<double> x_vec(x.data(), x.data() + x.size());
  return {x_vec, rnorm};
}


template <typename T> std::vector<double> DelaySignal(const std::vector<T> &signal, double delay_index) {
  std::vector<double> delay_signal(signal.size(), 0.0);

  for (size_t i = 0; i < signal.size(); i++) {
    double idx = static_cast<double>(i) + delay_index;
    double val = 0.0;
      
    InterpLinear(signal, idx, val);

    delay_signal[i] = val;
  }

  return delay_signal;
}


template <typename T> bool InterpLinear(const std::vector<T> &v, double index, double &out) {
  if (v.empty())
    return false;

  const double eps = 1e-12;
  const double max_idx = static_cast<double>(v.size() - 1);

  if (index < 0.0 && index > -eps)
    index = 0.0;

  if (index > max_idx && index < max_idx + eps)
    index = max_idx;

  if (index < 0.0 || index > max_idx)
    return false;

  const size_t i0 = static_cast<size_t>(std::floor(index));

  const size_t i1 = (i0 + 1 >= v.size()) ? i0 : i0 + 1;

  const double frac = index - static_cast<double>(i0);

  out = static_cast<double>(v[i0]) +
        (static_cast<double>(v[i1]) - static_cast<double>(v[i0])) * frac;

  return true;
}

static DelayFitResult FitDelayLeastSquaresImpl(const std::vector<double> &signal,
                                               const std::vector<double> &reference,
                                               double dt,
                                               double search_min,
                                               double search_max,
                                               double step,
                                               bool fit_offset) {
  DelayFitResult best{};
  best.chi2 = std::numeric_limits<double>::infinity();

  if (signal.size() < 3 || reference.size() < 3 || dt <= 0.0) {
    return best;
  }
  if (step <= 0.0) {
    step = dt;
  }
  if (search_min > search_max) {
    std::swap(search_min, search_max);
  }

  const double ref_max_t = (reference.size() - 1) * dt;

  for (double delay = search_min; delay <= search_max; delay += step) {
    double Sx = 0.0, Sy = 0.0, Sxx = 0.0, Sxy = 0.0;
    int n = 0;

    for (size_t i = 0; i < signal.size(); ++i) {
      const double t = static_cast<double>(i) * dt;
      const double tref = t - delay;
      // if (tref < 0.0 || tref > ref_max_t) {
      //   continue;
      // }

      double ref_val = 0.0;
      const double ref_idx = tref / dt;
      if (!InterpLinear(reference, ref_idx, ref_val)) {
        ref_val = 0;
      }

      const double y = signal[i];
      const double x = ref_val;
      Sx += x;
      Sy += y;
      Sxx += x * x;
      Sxy += x * y;
      ++n;
    }

    if (n < 3) {
      continue;
    }

    double a = 0.0;
    double b = 0.0;
    if (fit_offset) {
      const double denom = n * Sxx - Sx * Sx;
      if (std::fabs(denom) > 0.0) {
        a = (n * Sxy - Sx * Sy) / denom;
        b = (Sy - a * Sx) / n;
      } else {
        b = Sy / n;
      }
    } else if (std::fabs(Sxx) > 0.0) {
      a = Sxy / Sxx;
    }
    double chi2 = 0.0;
    double y_mean = Sy / n;
    double sst = 0.0;

    for (size_t i = 0; i < signal.size(); ++i) {
      const double t = static_cast<double>(i) * dt;
      const double tref = t - delay;
      // if (tref < 0.0 || tref > ref_max_t) {
      //   continue;
      // }

      double ref_val = 0.0;
      const double ref_idx = tref / dt;
      if (!InterpLinear(reference, ref_idx, ref_val)) {
        ref_val = 0;
      }

      const double y = signal[i];
      const double y_fit = a * ref_val + b;
      const double resid = y - y_fit;
      chi2 += resid * resid;

      const double dy = y - y_mean;
      sst += dy * dy;
    }
    // std::cout << delay << "," << chi2 << ","<< a << std::endl;
    if (chi2 < best.chi2) {
      best.delay = delay;
      best.scale = a;
      best.offset = b;
      best.chi2 = chi2;
      best.ndf = static_cast<double>(n - (fit_offset ? 2 : 1));
      best.r2 = (sst > 0.0) ? (1.0 - chi2 / sst) : 0.0;
    }
  }
  if (best.scale < 0) {
    std::cerr << "[AnalysisUtils] FitDelayLeastSquares: negative best scale=" << best.scale
              << " delay=" << best.delay << " chi2=" << best.chi2
              << " r2=" << best.r2 << std::endl;
  }
  return best;
}

DelayFitResult FitDelayLeastSquares(const std::vector<double> &signal,
                                    const std::vector<double> &reference,
                                    double dt,
                                    double search_min,
                                    double search_max,
                                    double step) {
  return FitDelayLeastSquaresImpl(signal, reference, dt, search_min, search_max,
                                  step, false);
}

DelayFitResult FitDelayLeastSquaresWithOffset(const std::vector<double> &signal,
                                              const std::vector<double> &reference,
                                              double dt,
                                              double search_min,
                                              double search_max,
                                              double step) {
  return FitDelayLeastSquaresImpl(signal, reference, dt, search_min, search_max,
                                  step, true);
}

DelayFitResult FitDelayLeastSquaresImpl_TGraph(const std::vector<double> &signal,
                                               const TGraph *gRef,
                                               double dt,
                                               double search_min,
                                               double search_max,
                                               double step,
                                               bool fit_offset,
                                               int fit_start,int fit_end) {
  DelayFitResult best{};
  best.chi2 = std::numeric_limits<double>::infinity();

  for (double delay = search_min; delay <= search_max; delay += step) {
    double Sx = 0.0, Sy = 0.0, Sxx = 0.0, Sxy = 0.0;
    int n = 0;

    for (size_t i = fit_start; i < fit_end; ++i) {
      const double t = static_cast<double>(i) * dt;
      const double tref = t - delay;

      double ref_val = 0.0;

      double xmin = gRef->GetX()[0];
      double xmax = gRef->GetX()[gRef->GetN() - 1];
      if (tref < xmin || tref > xmax) {
        ref_val = 0;
      } else {
        ref_val = gRef->Eval(tref);
      }

      const double y = signal[i];
      const double x = ref_val;
      Sx += x;
      Sy += y;
      Sxx += x * x;
      Sxy += x * y;
      ++n;
    }

    if (n < 3) {
      continue;
    }

    double a = 0.0;
    double b = 0.0;
    if (fit_offset) {
      const double denom = n * Sxx - Sx * Sx;
      if (std::fabs(denom) > 0.0) {
        a = (n * Sxy - Sx * Sy) / denom;
        b = (Sy - a * Sx) / n;
      } else {
        b = Sy / n;
      }
    } else if (std::fabs(Sxx) > 0.0) {
      a = Sxy / Sxx;
    }
    double chi2 = 0.0;
    double y_mean = Sy / n;
    double sst = 0.0;

    for (size_t i = fit_start; i < fit_end; ++i) {
      const double t = static_cast<double>(i) * dt;
      const double tref = t - delay;
      // if (tref < 0.0 || tref > ref_max_t) {
      //   continue;
      // }

      double ref_val = 0.0;

      double xmin = gRef->GetX()[0];
      double xmax = gRef->GetX()[gRef->GetN() - 1];
      if (tref < xmin || tref > xmax) {
        ref_val = 0;
      } else {
        ref_val = gRef->Eval(tref);
      }

      const double y = signal[i];
      const double y_fit = a * ref_val + b;
      const double resid = y - y_fit;
      chi2 += resid * resid;

      const double dy = y - y_mean;
      sst += dy * dy;
    }
    // std::cout << delay << "," << chi2 << ","<< a << std::endl;
    if (chi2 < best.chi2) {
      best.delay = delay;
      best.scale = a;
      best.offset = b;
      best.chi2 = chi2;
      best.ndf = static_cast<double>(n - (fit_offset ? 2 : 1));
      best.r2 = (sst > 0.0) ? (1.0 - chi2 / sst) : 0.0;
    }
  }
  if (best.scale < 0) {
    std::cerr << "[AnalysisUtils] FitDelayLeastSquares: negative best scale=" << best.scale
              << " delay=" << best.delay << " chi2=" << best.chi2
              << " r2=" << best.r2 << std::endl;
  }
  return best;
}


double FindT0(const std::vector<double> &adc, TGraph *gTemplate) {
  double bestT0 = 0;
  double bestCorr = -1e30;

  for (double t0 = -100; t0 < 300; t0 += 0.1) {
    double corr = Correlation(adc, gTemplate, t0);

    if (corr > bestCorr) {
      bestCorr = corr;
      bestT0 = t0;
    }
  }

  return bestT0;
}

double Correlation(const std::vector<double> &adc, TGraph *gTemplate, double t0) {
  double sumDT = 0.0;
  double sumD2 = 0.0;
  double sumT2 = 0.0;

  for (size_t i = 0; i < adc.size(); i++) {
    double tData = i * 25.0; // ns

    double D = adc[i];

    double T = gTemplate->Eval(tData - t0);

    sumDT += D * T;
    sumD2 += D * D;
    sumT2 += T * T;
  }

  if (sumD2 <= 0 || sumT2 <= 0)
    return 0;

  return sumDT / std::sqrt(sumD2 * sumT2);
}

// 显式模板实例化，供其他翻译单元（如 DeconvAnalysisScript.cpp）链接
template std::vector<double> DelaySignal<short>(const std::vector<short>&, double);
template std::vector<double> DelaySignal<double>(const std::vector<double>&, double);
template bool InterpLinear<short>(const std::vector<short>&, double, double&);
template bool InterpLinear<double>(const std::vector<double>&, double, double&);

}


