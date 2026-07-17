#include "AnalysisUtils.h"

#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TDirectory.h"

#include <Eigen/Dense>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

int main() {
    const int n_tier = 5;
    const int n_sig_samples = 20;
    double scale = 4200;

    // ====== Phase 1: build deconvolution matrix ======
    std::cout << "[Phase 1] Building deconvolution matrix" << std::endl;
    TFile *fmatrix = TFile::Open("../../simulation/analysis/matrix.root", "READ");
    if (!fmatrix || fmatrix->IsZombie()) {
        std::cerr << "Error: Cannot open matrix.root" << std::endl;
        return 1;
    }
    TTree *tmatrix = (TTree *)fmatrix->Get("matrix_tree");
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

    Eigen::MatrixXd A = AnalysisUtils::DeconvBaseMatrix(sigzero_samples[0], n_sig_samples);
    std::cout << "Matrix A size: " << A.rows() << "x" << A.cols() << std::endl;

    // ====== Phase 2: Lambda scan ======
    const int n_tests = 5000;
    const double noise_level = 1.0;
    const double signal_rms = 12.0;

    // lambda list: 0 (NNLS) + log-spaced
    std::vector<double> lambdas;
    lambdas.push_back(0.0);
    for (double lam = 0.01; lam <= 200.0; lam *= 2.0)
        lambdas.push_back(lam);

    struct Stats {
        std::vector<double> sum_x, sum2_x, sum_err, sum2_err;
        int count = 0;
        Stats() : sum_x(n_tier, 0), sum2_x(n_tier, 0),
                  sum_err(n_tier, 0), sum2_err(n_tier, 0) {}
    };

    std::vector<Stats> stats_tik(lambdas.size());
    std::vector<Stats> stats_lasso(lambdas.size());
    Stats stats_nnls;

    TRandom3 rng(42);

    for (int test = 0; test < n_tests; test++) {
        if (test % 1000 == 0)
            std::cout << "Progress: " << test << "/" << n_tests << std::endl;

        Eigen::VectorXd x_true(n_tier);
        for (int i = 0; i < n_tier; i++)
            if(i==2)   x_true(i) = rng.PoissonD(6);
            else x_true(i) = 0;

        Eigen::VectorXd b = A * x_true;
        std::vector<double> b_vec(b.size());
        for (int i = 0; i < b.size(); i++)
            b_vec[i] = b(i);
        for (int i = 0; i < (int)b_vec.size(); i++)
            b_vec[i] += rng.Gaus(0.0, noise_level * signal_rms);

        // NNLS
        {
            auto [x_vec, rn] = AnalysisUtils::SolveNNLSLawHanson(b_vec, A);
            for (int i = 0; i < n_tier; i++) {
                double xx = x_vec[i], err = xx - x_true(i);
                stats_nnls.sum_x[i] += xx; stats_nnls.sum2_x[i] += xx * xx;
                stats_nnls.sum_err[i] += err; stats_nnls.sum2_err[i] += err * err;
            }
            stats_nnls.count++;
        }

        for (size_t ilam = 1; ilam < lambdas.size(); ilam++) {
            double lam = lambdas[ilam];

            auto [xt, rt] = AnalysisUtils::SolveNNLSTikhonov(b_vec, A, lam);
            for (int i = 0; i < n_tier; i++) {
                double xx = xt[i], err = xx - x_true(i);
                stats_tik[ilam].sum_x[i] += xx; stats_tik[ilam].sum2_x[i] += xx * xx;
                stats_tik[ilam].sum_err[i] += err; stats_tik[ilam].sum2_err[i] += err * err;
            }
            stats_tik[ilam].count++;

            auto [xl, rl] = AnalysisUtils::SolveNNLSLasso(b_vec, A, lam);
            for (int i = 0; i < n_tier; i++) {
                double xx = xl[i], err = xx - x_true(i);
                stats_lasso[ilam].sum_x[i] += xx; stats_lasso[ilam].sum2_x[i] += xx * xx;
                stats_lasso[ilam].sum_err[i] += err; stats_lasso[ilam].sum2_err[i] += err * err;
            }
            stats_lasso[ilam].count++;
        }
    }

    std::cout << "Done. " << n_tests << " tests completed.\n" << std::endl;

    // ====== Phase 3: compute bias^2, variance, MSE ======
    auto decomp = [&](const Stats &s) {
        std::vector<double> bias2(n_tier), variance(n_tier), mse(n_tier);
        for (int i = 0; i < n_tier; i++) {
            double n = s.count;
            double mp = s.sum_x[i] / n;
            double mt = 6.0;  // Poisson(6) true mean
            bias2[i] = (mp - mt) * (mp - mt);
            variance[i] = (s.sum2_x[i] / n) - (mp * mp);
            if (variance[i] < 0) variance[i] = 0;
            mse[i] = s.sum2_err[i] / n;
        }
        return std::make_tuple(bias2, variance, mse);
    };

    auto [b2_nnls, var_nnls, mse_nnls] = decomp(stats_nnls);

    // ====== Phase 4: output & plot ======
    TFile *fout = TFile::Open("lambda_scan.root", "RECREATE");

    std::cout << "\n========== NNLS Baseline ==========" << std::endl;
    std::cout << "Tier |  bias^2   variance    MSE" << std::endl;
    for (int i = 0; i < n_tier; i++)
        std::cout << "  " << i << "  | " << std::fixed << std::setprecision(4)
                  << std::setw(8) << b2_nnls[i] << "  " << var_nnls[i] << "  " << mse_nnls[i] << std::endl;

    std::cout << "\n========== Tikhonov vs Lambda ==========" << std::endl;
    std::cout << "lambda  | avg_bias^2  avg_variance  avg_MSE" << std::endl;
    for (size_t ilam = 1; ilam < lambdas.size(); ilam++) {
        auto [b2, v, m] = decomp(stats_tik[ilam]);
        double ab = 0, av = 0, am = 0;
        for (int i = 0; i < n_tier; i++) { ab += b2[i]; av += v[i]; am += m[i]; }
        std::cout << std::setw(6) << lambdas[ilam] << " | "
                  << std::setw(10) << std::setprecision(4) << ab / n_tier << "  "
                  << std::setw(12) << av / n_tier << "  "
                  << std::setw(9) << am / n_tier << std::endl;
    }

    std::cout << "\n========== Lasso vs Lambda ==========" << std::endl;
    std::cout << "lambda  | avg_bias^2  avg_variance  avg_MSE" << std::endl;
    for (size_t ilam = 1; ilam < lambdas.size(); ilam++) {
        auto [b2, v, m] = decomp(stats_lasso[ilam]);
        double ab = 0, av = 0, am = 0;
        for (int i = 0; i < n_tier; i++) { ab += b2[i]; av += v[i]; am += m[i]; }
        std::cout << std::setw(6) << lambdas[ilam] << " | "
                  << std::setw(10) << std::setprecision(4) << ab / n_tier << "  "
                  << std::setw(12) << av / n_tier << "  "
                  << std::setw(9) << am / n_tier << std::endl;
    }

    // ---- Per-tier Tikhonov: bias^2, variance, MSE vs lambda ----
    TDirectory *dir_tik = fout->mkdir("Tikhonov");
    for (int tier = 0; tier < n_tier; tier++) {
        TGraph *g_b = new TGraph(), *g_v = new TGraph(), *g_m = new TGraph();
        for (size_t ilam = 1; ilam < lambdas.size(); ilam++) {
            auto [b2, v, m] = decomp(stats_tik[ilam]);
            int n = g_b->GetN();
            g_b->SetPoint(n, lambdas[ilam], b2[tier]);
            g_v->SetPoint(n, lambdas[ilam], v[tier]);
            g_m->SetPoint(n, lambdas[ilam], m[tier]);
        }
        TMultiGraph *mg = new TMultiGraph();
        g_b->SetLineColor(kRed);   g_b->SetMarkerColor(kRed);   g_b->SetMarkerStyle(21);
        g_v->SetLineColor(kBlue);  g_v->SetMarkerColor(kBlue);  g_v->SetMarkerStyle(22);
        g_m->SetLineColor(kBlack); g_m->SetMarkerColor(kBlack); g_m->SetMarkerStyle(20);
        mg->Add(g_b); mg->Add(g_v); mg->Add(g_m);
        mg->SetTitle(TString::Format("Tikhonov Tier %d: Bias-Variance;#lambda;Bias^{2}/Variance/MSE", tier));

        TCanvas *c = new TCanvas(TString::Format("c_tik_t%d", tier),
                                 TString::Format("Tikhonov Tier %d", tier), 800, 600);
        mg->Draw("ALP"); c->SetLogx();
        TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(g_b, "Bias^{2}", "lp");
        leg->AddEntry(g_v, "Variance", "lp");
        leg->AddEntry(g_m, "MSE", "lp");
        leg->Draw();
        dir_tik->WriteTObject(c);
    }

    // ---- Per-tier Lasso ----
    TDirectory *dir_lasso = fout->mkdir("Lasso");
    for (int tier = 0; tier < n_tier; tier++) {
        TGraph *g_b = new TGraph(), *g_v = new TGraph(), *g_m = new TGraph();
        for (size_t ilam = 1; ilam < lambdas.size(); ilam++) {
            auto [b2, v, m] = decomp(stats_lasso[ilam]);
            int n = g_b->GetN();
            g_b->SetPoint(n, lambdas[ilam], b2[tier]);
            g_v->SetPoint(n, lambdas[ilam], v[tier]);
            g_m->SetPoint(n, lambdas[ilam], m[tier]);
        }
        TMultiGraph *mg = new TMultiGraph();
        g_b->SetLineColor(kRed);   g_b->SetMarkerColor(kRed);   g_b->SetMarkerStyle(21);
        g_v->SetLineColor(kBlue);  g_v->SetMarkerColor(kBlue);  g_v->SetMarkerStyle(22);
        g_m->SetLineColor(kBlack); g_m->SetMarkerColor(kBlack); g_m->SetMarkerStyle(20);
        mg->Add(g_b); mg->Add(g_v); mg->Add(g_m);
        mg->SetTitle(TString::Format("Lasso Tier %d: Bias-Variance;#lambda;Bias^{2}/Variance/MSE", tier));

        TCanvas *c = new TCanvas(TString::Format("c_lasso_t%d", tier),
                                 TString::Format("Lasso Tier %d", tier), 800, 600);
        mg->Draw("ALP"); c->SetLogx();
        TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(g_b, "Bias^{2}", "lp");
        leg->AddEntry(g_v, "Variance", "lp");
        leg->AddEntry(g_m, "MSE", "lp");
        leg->Draw();
        dir_lasso->WriteTObject(c);
    }

    // ---- Summary: average MSE vs lambda (Tikhonov vs Lasso vs NNLS baseline) ----
    TCanvas *c_sum = new TCanvas("c_summary", "Average MSE vs Lambda", 900, 600);
    TGraph *g_tm = new TGraph(), *g_lm = new TGraph();
    double nnls_avg = 0;
    for (int i = 0; i < n_tier; i++) nnls_avg += mse_nnls[i];
    nnls_avg /= n_tier;
    TGraph *g_nline = new TGraph();
    g_nline->SetPoint(0, lambdas[1] * 0.5, nnls_avg);
    g_nline->SetPoint(1, lambdas.back() * 2, nnls_avg);
    g_nline->SetLineColor(kGreen + 2); g_nline->SetLineStyle(2);

    for (size_t ilam = 1; ilam < lambdas.size(); ilam++) {
        auto [bt, vt, mt] = decomp(stats_tik[ilam]);
        auto [bl, vl, ml] = decomp(stats_lasso[ilam]);
        double amt = 0, aml = 0;
        for (int i = 0; i < n_tier; i++) { amt += mt[i]; aml += ml[i]; }
        g_tm->SetPoint(g_tm->GetN(), lambdas[ilam], amt / n_tier);
        g_lm->SetPoint(g_lm->GetN(), lambdas[ilam], aml / n_tier);
    }
    g_tm->SetLineColor(kRed);  g_tm->SetMarkerColor(kRed);  g_tm->SetMarkerStyle(20);
    g_lm->SetLineColor(kBlue); g_lm->SetMarkerColor(kBlue); g_lm->SetMarkerStyle(21);

    TMultiGraph *mg_sum = new TMultiGraph();
    mg_sum->SetTitle("Average MSE over 5 Tiers vs #lambda;#lambda;Average MSE");
    mg_sum->Add(g_tm); mg_sum->Add(g_lm);
    mg_sum->Draw("ALP");
    g_nline->Draw("L SAME");
    c_sum->SetLogx();
    TLegend *leg_sum = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg_sum->AddEntry(g_tm, "Tikhonov", "lp");
    leg_sum->AddEntry(g_lm, "Lasso", "lp");
    leg_sum->AddEntry(g_nline, "NNLS baseline", "l");
    leg_sum->Draw();
    fout->WriteTObject(c_sum);

    fout->Write();
    fout->Close();
    std::cout << "\nResults saved to lambda_scan.root" << std::endl;
    return 0;
}