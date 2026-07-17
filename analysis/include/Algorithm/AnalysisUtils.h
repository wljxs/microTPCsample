#pragma once

#include "Event/DataModel.h"
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <TFile.h>
#include <TGraph.h>

// 前向声明
class Event;
class Cluster;

/**
 * @brief 分析工具函数命名空间
 *
 * 提供通用的分析算法和辅助函数，供脚本和分析模块使用
 */
namespace AnalysisUtils {

Track FitTrack(const std::vector<TVector3>& hits);

std::pair<double, double> GetRange(const std::vector<double>& v);

double CalculateMean(const std::vector<double>& values);

double CalculateRMS(const std::vector<double>& values);

struct DelayFitResult {
	double delay = 0.0;
	double scale = 0.0;
	double offset = 0.0;
	double chi2 = 0.0;
	double ndf = 0.0;
	double r2 = 0.0;
};

template <typename T> std::vector<double> DelaySignal(const std::vector<T> &signal, double delay_index);

DelayFitResult FitDelayLeastSquares(const std::vector<double>& signal,
									const std::vector<double>& reference,
									double dt,
									double search_min,
									double search_max,
									double step);// 给定一个信号和一个参考波形，搜索最佳的时间延迟，使得将参考波形按该延迟平移后与信号的拟合度（R^2）最高，返回最佳延迟、缩放因子、拟合优度等信息

DelayFitResult FitDelayLeastSquaresWithOffset(const std::vector<double>& signal,
										const std::vector<double>& reference,
										double dt,
										double search_min,
										double search_max,
										double step);// 和上一个函数一样，但多了一个offset参数，也就是拟合函数是y=a*ref+b，而不是y=a*ref，这样可以更好地拟合有基线的信号

DelayFitResult FitDelayLeastSquaresImpl_TGraph(const std::vector<double> &signal,
                                               const TGraph *gRef,
                                               double dt,
                                               double search_min,
                                               double search_max,
                                               double step,
                                               bool fit_offset,
                                               int fit_start,int fit_end);// 给定一个信号和一个参考波形，搜索最佳的时间延迟，使得将参考波形按该延迟平移后与信号的拟合度（R^2）最高，返回最佳延迟、缩放因子、拟合优度等信息。这个函数是上面两个函数的实现，区别是参考波形是一个TGraph，而不是一个vector<double>，并且可以指定拟合的起始和结束点
void FFTAnalyzer(Cluster& cluster, Event& evt, int det_id);

template <typename T> bool InterpLinear(const std::vector<T> &v, double index, double &out);
std::pair<std::vector<double>,double> SolveNNLSLawHanson(std::vector<double> vec_b, const Eigen::MatrixXd& A);//基于NNLS给定一个b和矩阵A，求解Ax=b中x的非负解

/**
 * @brief 带一阶差分Tikhonov正则化的NNLS反卷积
 * @param b_vec  观测信号（ADC波形片段）
 * @param A      卷积核矩阵（由DeconvBaseMatrix构造）
 * @param lambda 正则化强度，越大解越平滑；0等价于原始NNLS
 * @return pair<解向量x（非负电荷分布）, 增广系统的残差范数>
 *
 * 求解:  min ||A*x - b||^2 + lambda * ||L*x||^2,  s.t. x >= 0
 * 其中 L 是一阶差分矩阵，惩罚相邻 tier 之间的电荷跳变。
 * 通过构造增广矩阵 [A; sqrt(lambda)*L] 和 [b; 0] 后调用 NNLS 实现。
 */
std::pair<std::vector<double>, double> SolveNNLSTikhonov(
    std::vector<double> b_vec, const Eigen::MatrixXd &A, double lambda);

/**
 * @brief 带L1正则化（LASSO）的NNLS反卷积
 * @param b_vec  观测信号（ADC波形片段）
 * @param A      卷积核矩阵（由DeconvBaseMatrix构造）
 * @param lambda 正则化强度，越大解越稀疏（更多 component 被压到 0）；0等价于原始NNLS
 * @return pair<解向量x（非负稀疏电荷分布）, 残差范数>
 *
 * 求解:  min ||A*x - b||^2 + lambda * ||x||_1,  s.t. x >= 0
 * 因非负约束，||x||_1 = sum(x_i)，通过修改Lawson-Hanson的KKT条件实现：
 * 将标准NNLS的对偶变量阈值从0上移至lambda/2，解方程时从右端项减去lambda/2。
 */
std::pair<std::vector<double>, double> SolveNNLSLasso(
    std::vector<double> b_vec, const Eigen::MatrixXd &A, double lambda);

Eigen::MatrixXd DeconvBaseMatrix(std::vector<double> vec_sigzero, int n_sig_track, int n_tier);//从文件中读到要取的波形，按照层数生成一个矩阵
Eigen::MatrixXd DeconvBaseMatrix(std::vector<std::vector<double>> vec_sigzeros,int n_sig_track);//用多个波形构建矩阵
double FindT0(const std::vector<double>& adc, TGraph* gTemplate);//给定一个adc波形和一个模板波形，找到adc波形的t0，也就是adc波形的上升沿与模板波形的上升沿对齐的时间点
double Correlation(const std::vector<double>& adc, TGraph* gTemplate, double t0);//给定一个adc波形和一个模板波形，计算adc波形与模板波形在t0时间点的相关系数
}  // namespace AnalysisUtils
