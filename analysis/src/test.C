#include "AnalysisUtils.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TRandom3.h"

#include <Eigen/Dense>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

int main() {
    const int n_tier = 5; 
    const int n_sig_samples = 10;
    double scale = 4200; //这个是增益的缩放因子，主要用于将模拟信号的幅值调整到实际实验的范围内，确保反卷积和后续分析的准确性。不过目前只是反卷积矩阵有这个
    // ====== Phase 1: 构建反卷积矩阵 ======
    std::cout << "[Phase 1] Building deconvolution matrix from raw/matrix.root" << std::endl;

    TFile *fmatrix = TFile::Open("../../simulation/analysis/matrix.root", "READ");
    if (!fmatrix || fmatrix->IsZombie()) {
        std::cerr << "Error: Cannot open matrix.root" << std::endl;
        return 1;
    }
    TTree *tmatrix = (TTree *)fmatrix->Get("matrix_tree");
    if (!tmatrix) {
        std::cerr << "Error: Cannot find matrix_tree" << std::endl;
        return 1;
    }
    std::vector<double> *sigx_conv = nullptr;
    std::vector<double> *sigy_conv = nullptr;
    tmatrix->SetBranchAddress("sigx_conv", &sigx_conv);
    tmatrix->SetBranchAddress("sigy_conv", &sigy_conv);

    std::vector<std::vector<double>> sigzero_samples[2];
    for (int entry = 0; entry < n_tier; entry++) {
        tmatrix->GetEntry(entry);
        std::vector<double> sig_temp[2];
        for (size_t i = 0; i < sigx_conv->size(); i++) {
            if (i % 5 == 0 && i < n_sig_samples * 5 - 1) {
                sig_temp[0].push_back((*sigx_conv)[i] * scale);
                sig_temp[1].push_back(-(*sigy_conv)[i] * scale);
            }
        }
        sigzero_samples[0].push_back(sig_temp[0]);
        sigzero_samples[1].push_back(sig_temp[1]);
    }
    fmatrix->Close();

    Eigen::MatrixXd deconvMatrix[2];
    deconvMatrix[0] = AnalysisUtils::DeconvBaseMatrix(sigzero_samples[0], n_sig_samples);
    deconvMatrix[1] = AnalysisUtils::DeconvBaseMatrix(sigzero_samples[1], n_sig_samples);

    std::cout << "Matrix A (X) size: " << deconvMatrix[0].rows() << "x" << deconvMatrix[0].cols() << std::endl;

    // ====== Phase 2: 精度测试（10000轮） ======
    std::cout << "\n========== Phase 2: Deconvolution Accuracy Test (10000 rounds) ==========\n" << std::endl;

    const int n_tests = 10000;
    const double noise_level = 1;
    const double lambda_lasso = 20;
    const double lambda_tikhonov = 20;

    TRandom3 rng(42);

    // 每层残差直方图 (残差 = x_pred - x_true)
    TH1D *h_res_nnls[n_tier];
    TH1D *h_res_lasso[n_tier];
    TH1D *h_res_tikhonov[n_tier];
    for (int tier = 0; tier < n_tier; tier++) {
        h_res_nnls[tier] = new TH1D(
            Form("h_res_nnls_tier%d", tier),
            Form("NNLS Residual Tier %d;x_pred - x_true;Entries", tier),
            1000, -5, 5);
        h_res_lasso[tier] = new TH1D(
            Form("h_res_lasso_tier%d", tier),
            Form("Lasso Residual Tier %d;x_pred - x_true;Entries", tier),
            1000, -5, 5);
        h_res_tikhonov[tier] = new TH1D(
            Form("h_res_tikhonov_tier%d", tier),
            Form("Tikhonov Residual Tier %d;x_pred - x_true;Entries", tier),
            1000, -5, 5);
    }

    // 每层相对误差直方图
    TH1D *h_relerr_nnls[n_tier];
    TH1D *h_relerr_lasso[n_tier];
    TH1D *h_relerr_tikhonov[n_tier];
    for (int tier = 0; tier < n_tier; tier++) {
        h_relerr_nnls[tier] = new TH1D(
            Form("h_relerr_nnls_tier%d", tier),
            Form("NNLS Relative Error Tier %d;|pred-true|/true;Entries", tier),
            1000, -5, 5);
        h_relerr_lasso[tier] = new TH1D(
            Form("h_relerr_lasso_tier%d", tier),
            Form("Lasso Relative Error Tier %d;|pred-true|/true;Entries", tier),
            100, -5, 5);
        h_relerr_tikhonov[tier] = new TH1D(
            Form("h_relerr_tikhonov_tier%d", tier),
            Form("Tikhonov Relative Error Tier %d;|pred-true|/true;Entries", tier),
            100, -5, 5);
    }

    // 每层反卷积电荷直方图（x_pred 和 x_true）
    TH1D *h_charge_true[n_tier];
    TH1D *h_charge_nnls[n_tier];
    TH1D *h_charge_lasso[n_tier];
    TH1D *h_charge_tikhonov[n_tier];
    for (int tier = 0; tier < n_tier; tier++) {
        h_charge_true[tier] = new TH1D(
            Form("h_charge_true_tier%d", tier),
            Form("True Charge Tier %d;Charge;Entries", tier),
            60, 0, 30);
        h_charge_nnls[tier] = new TH1D(
            Form("h_charge_nnls_tier%d", tier),
            Form("NNLS Charge Tier %d;Charge;Entries", tier),
            60, 0, 30);
        h_charge_lasso[tier] = new TH1D(
            Form("h_charge_lasso_tier%d", tier),
            Form("Lasso Charge Tier %d;Charge;Entries", tier),
            60, 0, 30);
        h_charge_tikhonov[tier] = new TH1D(
            Form("h_charge_tikhonov_tier%d", tier),
            Form("Tikhonov Charge Tier %d;Charge;Entries", tier),
            60, 0, 30);
    }

    // 收集所有轮次的残差向量（用于协方差矩阵）
    std::vector<Eigen::VectorXd> all_residuals_nnls;
    std::vector<Eigen::VectorXd> all_residuals_lasso;
    std::vector<Eigen::VectorXd> all_residuals_tikhonov;
    all_residuals_nnls.reserve(n_tests);
    all_residuals_lasso.reserve(n_tests);
    all_residuals_tikhonov.reserve(n_tests);

    TGraph *gr_b = new TGraph();    
    for (int test = 0; test < n_tests; test++) {
        if (test % 1000 == 0)
            std::cout << "Progress: " << test << "/" << n_tests << std::endl;

        // 生成随机 x_true
        Eigen::VectorXd x_true_eig(n_tier);
        for (int i = 0; i < n_tier; i++)
            x_true_eig(i) = rng.PoissonD(6);  // 随机生成每层的真实信号强度

        // 正向 b = A * x_true
        Eigen::VectorXd b = deconvMatrix[0] * x_true_eig;
        std::vector<double> b_vec(b.size());
        for (int i = 0; i < b.size(); i++)
            b_vec[i] = b(i);

        if(test==0){
            for(int i=0;i<b_vec.size();i++){
                gr_b->SetPoint(i,i,b_vec[i]);
            }
        }
        // 加噪声
        if (noise_level > 0) {
            double signal_rms = 12;
            for (int i = 0; i < (int)b_vec.size(); i++)
                b_vec[i] += rng.Gaus(0.0, noise_level * signal_rms);
        }

        // NNLS 反卷积
        auto [x_nnls, rnorm_nnls] = AnalysisUtils::SolveNNLSLawHanson(b_vec, deconvMatrix[0]);
        // Lasso 反卷积
        auto [x_lasso, rnorm_lasso] = AnalysisUtils::SolveNNLSLasso(b_vec, deconvMatrix[0], lambda_lasso);
        auto [x_tikhonov, rnorm_tikhonov] = AnalysisUtils::SolveNNLSTikhonov(b_vec, deconvMatrix[0], lambda_tikhonov);

        // 收集残差向量
        Eigen::VectorXd res_vec_nnls(n_tier), res_vec_lasso(n_tier), res_vec_tikhonov(n_tier);
        for (int i = 0; i < n_tier; i++) {
            res_vec_nnls(i) = x_nnls[i] - x_true_eig(i);
            res_vec_lasso(i) = x_lasso[i] - x_true_eig(i);
            res_vec_tikhonov(i) = x_tikhonov[i] - x_true_eig(i);
        }
        all_residuals_nnls.push_back(res_vec_nnls);
        all_residuals_lasso.push_back(res_vec_lasso);
        all_residuals_tikhonov.push_back(res_vec_tikhonov);

        // 填充直方图
        for (int i = 0; i < n_tier; i++) {
            double res_nnls = x_nnls[i] - x_true_eig(i);
            double res_lasso = x_lasso[i] - x_true_eig(i);
            double res_tikhonov = x_tikhonov[i] - x_true_eig(i);
            double rel_nnls = (x_true_eig(i) > 1e-9) ? fabs(res_nnls) / x_true_eig(i) : 0.0;
            double rel_lasso = (x_true_eig(i) > 1e-9) ? fabs(res_lasso) / x_true_eig(i) : 0.0;
            double rel_tikhonov = (x_true_eig(i) > 1e-9) ? fabs(res_tikhonov) / x_true_eig(i) : 0.0;

            h_res_nnls[i]->Fill(res_nnls);
            h_res_lasso[i]->Fill(res_lasso);
            h_res_tikhonov[i]->Fill(res_tikhonov);
            h_relerr_nnls[i]->Fill(rel_nnls);
            h_relerr_lasso[i]->Fill(rel_lasso);
            h_relerr_tikhonov[i]->Fill(rel_tikhonov);

            h_charge_true[i]->Fill(x_true_eig(i));
            h_charge_nnls[i]->Fill(x_nnls[i]);
            h_charge_lasso[i]->Fill(x_lasso[i]);
            h_charge_tikhonov[i]->Fill(x_tikhonov[i]);
        }
    }

    std::cout << "Done. " << n_tests << " tests completed.\n" << std::endl;

    // ====== Phase 2.5: 计算协方差矩阵和关联矩阵 ======
    auto compute_cov_corr = [&](const std::vector<Eigen::VectorXd>& residuals,
                                 const std::string& label) {
        // 计算每层残差的均值
        Eigen::VectorXd mean_res = Eigen::VectorXd::Zero(n_tier);
        for (const auto& r : residuals)
            mean_res += r;
        mean_res /= residuals.size();

        // 计算协方差矩阵 cov[i][j] = E[(r_i - mu_i)(r_j - mu_j)]
        Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(n_tier, n_tier);
        for (const auto& r : residuals) {
            Eigen::VectorXd centered = r - mean_res;
            cov += centered * centered.transpose();
        }
        cov /= (residuals.size() - 1);

        // 计算关联矩阵 corr[i][j] = cov[i][j] / (sigma_i * sigma_j)
        Eigen::MatrixXd corr = Eigen::MatrixXd::Zero(n_tier, n_tier);
        for (int i = 0; i < n_tier; i++) {
            for (int j = 0; j < n_tier; j++) {
                double denom = std::sqrt(cov(i, i) * cov(j, j));
                corr(i, j) = (denom > 0) ? cov(i, j) / denom : 0.0;
            }
        }

        // 创建 TH2D 直方图（用 TString::Format 避免 Form 静态缓冲区冲突）
        TString name_cov = TString::Format("h_cov_%s", label.c_str());
        TString name_corr = TString::Format("h_corr_%s", label.c_str());
        TString title_cov = TString::Format("%s Residual Covariance Matrix;Tier i;Tier j", label.c_str());
        TString title_corr = TString::Format("%s Residual Correlation Matrix;Tier i;Tier j", label.c_str());

        TH2D *h_cov = new TH2D(name_cov, title_cov,
                                n_tier, 0, n_tier, n_tier, 0, n_tier);
        TH2D *h_corr = new TH2D(name_corr, title_corr,
                                 n_tier, 0, n_tier, n_tier, 0, n_tier);

        for (int i = 0; i < n_tier; i++) {
            for (int j = 0; j < n_tier; j++) {
                h_cov->SetBinContent(i + 1, j + 1, cov(i, j));
                h_corr->SetBinContent(i + 1, j + 1, corr(i, j));
            }
        }

        // 打印矩阵
        std::cout << "\n--- " << label << " Covariance Matrix ---" << std::endl;
        std::cout << std::fixed << std::setprecision(4);
        for (int i = 0; i < n_tier; i++) {
            for (int j = 0; j < n_tier; j++)
                std::cout << std::setw(10) << cov(i, j);
            std::cout << std::endl;
        }
        std::cout << "\n--- " << label << " Correlation Matrix ---" << std::endl;
        for (int i = 0; i < n_tier; i++) {
            for (int j = 0; j < n_tier; j++)
                std::cout << std::setw(10) << corr(i, j);
            std::cout << std::endl;
        }

        return std::make_tuple(h_cov, h_corr);
    };

    auto [h_cov_nnls, h_corr_nnls] = compute_cov_corr(all_residuals_nnls, "NNLS");
    auto [h_cov_lasso, h_corr_lasso] = compute_cov_corr(all_residuals_lasso, "Lasso");
    auto [h_cov_tikhonov, h_corr_tikhonov] = compute_cov_corr(all_residuals_tikhonov, "Tikhonov");

    // ====== Phase 3: 输出直方图统计 ======
    std::cout << "========== Residual Summary ==========" << std::endl;
    for (int tier = 0; tier < n_tier; tier++) {
        std::cout << "Tier " << tier << ": "
                  << "NNLS mean=" << h_res_nnls[tier]->GetMean()
                  << " RMS=" << h_res_nnls[tier]->GetRMS()
                  << " | Lasso mean=" << h_res_lasso[tier]->GetMean()
                  << " RMS=" << h_res_lasso[tier]->GetRMS() 
                  << " | Tikhonov mean=" << h_res_tikhonov[tier]->GetMean()
                  << " RMS=" << h_res_tikhonov[tier]->GetRMS()
                  << std::endl;
    }
    std::cout << "\n========== Relative Error Summary ==========" << std::endl;
    for (int tier = 0; tier < n_tier; tier++) {
        std::cout << "Tier " << tier << ": "
                  << "NNLS mean=" << h_relerr_nnls[tier]->GetMean() * 100 << "%"
                  << " | Lasso mean=" << h_relerr_lasso[tier]->GetMean() * 100 << "%"
                  << " | Tikhonov mean=" << h_relerr_tikhonov[tier]->GetMean() * 100 << "%"
                  << std::endl;
    }

    std::cout << "\n========== Charge Distribution Summary ==========" << std::endl;
    for (int tier = 0; tier < n_tier; tier++) {
        std::cout << "Tier " << tier << ": "
                  << "True mean=" << h_charge_true[tier]->GetMean()
                  << " RMS=" << h_charge_true[tier]->GetRMS()
                  << " | NNLS mean=" << h_charge_nnls[tier]->GetMean()
                  << " RMS=" << h_charge_nnls[tier]->GetRMS()
                  << " | Lasso mean=" << h_charge_lasso[tier]->GetMean()
                  << " RMS=" << h_charge_lasso[tier]->GetRMS()
                  << " | Tikhonov mean=" << h_charge_tikhonov[tier]->GetMean()
                  << " RMS=" << h_charge_tikhonov[tier]->GetRMS()
                  << std::endl;
    }

    // ====== Phase 4: 保存到 ROOT 文件（按类型分目录） ======
    TFile *fout = TFile::Open("test_deconv_accuracy.root", "RECREATE");

    fout->WriteTObject(gr_b, "gr_b");
    delete gr_b;

    TDirectory *dir_residuals = fout->mkdir("residuals");
    TDirectory *dir_relerr    = fout->mkdir("relerr");
    TDirectory *dir_charge    = fout->mkdir("charge");
    TDirectory *dir_matrix    = fout->mkdir("matrix");

    // 保存协方差和关联矩阵
    dir_matrix->cd();
    h_cov_nnls->Write();
    h_corr_nnls->Write();
    h_cov_lasso->Write();
    h_corr_lasso->Write();
    h_cov_tikhonov->Write();
    h_corr_tikhonov->Write();
    delete h_cov_nnls; delete h_corr_nnls;
    delete h_cov_lasso; delete h_corr_lasso;
    delete h_cov_tikhonov; delete h_corr_tikhonov;

    for (int tier = 0; tier < n_tier; tier++) {
        dir_residuals->cd();
        h_res_nnls[tier]->Write();
        h_res_lasso[tier]->Write();
        h_res_tikhonov[tier]->Write();

        dir_relerr->cd();
        h_relerr_nnls[tier]->Write();
        h_relerr_lasso[tier]->Write();
        h_relerr_tikhonov[tier]->Write();

        dir_charge->cd();
        h_charge_true[tier]->Write();
        h_charge_nnls[tier]->Write();
        h_charge_lasso[tier]->Write();
        h_charge_tikhonov[tier]->Write();

        delete h_res_nnls[tier];
        delete h_res_lasso[tier];
        delete h_res_tikhonov[tier];
        delete h_relerr_nnls[tier];
        delete h_relerr_lasso[tier];
        delete h_relerr_tikhonov[tier];
        delete h_charge_true[tier];
        delete h_charge_nnls[tier];
        delete h_charge_lasso[tier];
        delete h_charge_tikhonov[tier];
    }
    fout->Close();
    std::cout << "\nHistograms saved to test_deconv_accuracy.root" << std::endl;

    return 0;
}