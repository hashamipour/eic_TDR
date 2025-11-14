//g++ -o DDIS_Plot_Combined $(root-config --cflags --glibs) DDIS_Plot_Combined.cpp Plotting.cpp

#include "Plotting.hpp"
#include <TFile.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include "TError.h"
#include <TCanvas.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TProfile2D.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TString.h>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root>" << std::endl;
        return 1;
    }
    
    gErrorIgnoreLevel = kWarning; // Suppress ROOT infos

    TString inputFileName = argv[1];

    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(120);
    gStyle->SetPalette(kBlueRedYellow);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleSize(0.048, "XYZ");
    gStyle->SetLabelSize(0.038, "XYZ");
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
    gStyle->SetOptFit(111);

    TFile* inputFile = TFile::Open(inputFileName);

    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << inputFileName << std::endl;
        return 1;
    }

    std::vector<PlotOptions*> plots;
    PlotOptions1D* plot_ptr = nullptr;

    //=================================================================
    // Q2/xy ANALYSIS PLOTS
    //=================================================================

    plot_ptr = new PlotOptions1D(
        {"h_Q2_truth", "h_Q2_EM", "h_Q2_DA","h_Q2_Sigma"},    // hist names
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Sigma"}, // legend entries
        {"hist", "pe", "pe","pe"},                             // draw options
        "Q^{2} Reconstruction Methods",                        // canvas title
        "Q^{2}",                                               // x label
        "# of events",                                         // y label
        "figs/Q2_hist.png",                                    // save name
        true,                                                  // isLogX (default is false)
        true                                                   // isLogY (default is false)
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);
    
    plots.push_back(new PlotOptionsResponseMatrix (
        "Corr_Q2_EM",
        "Q^{2} (true) [GeV]",
        "Q^{2} (EM) [GeV]",
        "figs/response_matrix_Q2_EM.png",
        true,  // isLogX
        true,  // isLogY
        std::make_pair(1.0, 300),
        std::make_pair(1.0, 300)
    ));

    plots.push_back(new PlotOptionsResponseMatrix (
        "Corr_Q2_DA",
        "Q^{2} (true) [GeV]",
        "Q^{2} (DA) [GeV]",
        "figs/response_matrix_Q2_DA.png",
        true,  // isLogX
        true,  // isLogY
        std::make_pair(1.0, 300),
        std::make_pair(1.0, 300)
    ));

    plots.push_back(new PlotOptionsResponseMatrix (
        "Corr_Q2_Sigma",
        "Q^{2} (true) [GeV]",
        "Q^{2} (Sigma) [GeV]",
        "figs/response_matrix_Q2_Esigma.png",
        true,  // isLogX
        true,  // isLogY
        std::make_pair(1.0, 300),
        std::make_pair(1.0, 300)
    ));

    // Response Matrices for x_{Bj}
    plots.push_back(new PlotOptionsResponseMatrix(
        "x_Corr_EM",
        "x_{Bj} (true)",
        "x_{Bj} (EM)",
        "figs/response_matrix_x_EM.png",
        true,  // isLogX
        true,   // isLogY
        {1e-3,0.3},
        {1e-3,0.3}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "x_Corr_DA",
        "x_{Bj} (true)",
        "x_{Bj} (DA)",
        "figs/response_matrix_x_DA.png",
        true,  // isLogX
        true,   // isLogY
        {1e-3,0.3},
        {1e-3,0.3}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "x_Corr_Sigma",
        "x_{Bj} (true)",
        "x_{Bj} (Sigma)",
        "figs/response_matrix_x_Sigma.png",
        true,  // isLogX
        true,   // isLogY
        {1e-3,0.3},
        {1e-3,0.5}
    ));

    // E-pz distribution plots
    plot_ptr = new PlotOptions1D(
        {"h_EPz_truth", "h_EPz"},
        {"MC Truth"   , "Reconstruction"},
        {"hist"       , "pe"},
        "Hadronic Final State E-p_{z}",
        "#Sigma(E-p_{z}) [GeV]",
        "Counts",
        "figs/EPz_distribution.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.2, 0.7, 0.4, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"h_EPz_truth", "h_EPz"},
        {"MC Truth"   , "Reconstruction"},
        {"hist"       , "pe"},
        "Hadronic Final State E-p_{z}",
        "#Sigma(E-p_{z}) [GeV]",
        "Counts",
        "figs/EPz_distribution_logY.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.2, 0.7, 0.4, 0.9);
    plots.push_back(plot_ptr);

    // eta_max distribution plots
    plot_ptr = new PlotOptions1D(
        {"h_eta_max_truth", "h_eta_max"},
        {"MC Truth", "Reconstruction"},
        {"hist", "pe"},
        "Maximum Pseudorapidity per Event",
        "#eta_{max}",
        "Counts",
        "figs/eta_max_distribution.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.2, 0.7, 0.4, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"h_eta_max_truth", "h_eta_max"},
        {"MC Truth", "Reconstruction"},
        {"hist", "pe"},
        "Maximum Pseudorapidity per Event",
        "#eta_{max}",
        "Counts",
        "figs/eta_max_distribution_logY.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.1, 0.7, 0.3, 0.9);
    plots.push_back(plot_ptr);

    // M_X^2 distribution plots
    plot_ptr = new PlotOptions1D(
        {"h_MX2_truth", "h_MX2"},
        {"MC Truth", "Reconstruction"},
        {"hist", "pe"},
        "Hadronic Invariant Mass Squared",
        "M_{X}^{2} [GeV^{2}]",
        "Counts",
        "figs/MX2_distribution.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.2, 0.7, 0.4, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"h_MX2_truth", "h_MX2"},
        {"MC Truth", "Reconstruction"},
        {"hist", "pe"},
        "Hadronic Invariant Mass Squared",
        "M_{X}^{2} [GeV^{2}]",
        "Counts",
        "figs/MX2_distribution_logY.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.2, 0.7, 0.4, 0.9);
    plots.push_back(plot_ptr);

    //=================================================================
    // MANDELSTAM t ANALYSIS PLOTS
    //=================================================================

    // Response matrices for t
    plots.push_back(new PlotOptionsResponseMatrix(
        "t_corr_B0",
        "Truth |t| [GeV^{2}]",
        "B0 Reco |t| [GeV^{2}]",
        "figs/response_matrix_t_B0.png",
        false,
        false,
        {0.0, 2.0},
        {0.0, 2.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "t_corr_RP",
        "Truth |t| [GeV^{2}]",
        "RP Reco |t| [GeV^{2}]",
        "figs/response_matrix_t_RP.png",
        true,
        true,
        {0.0, 0.5},
        {0.0, 0.5}
    ));

    // x_L response matrices
    plots.push_back(new PlotOptionsResponseMatrix(
        "xL_corr_B0",
        "Truth x_L",
        "B0 Reco x_L",
        "figs/response_matrix_xL_B0.png",
        false,
        false,
        {0.0, 2.0},
        {0.0, 2.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xL_corr_RP",
        "Truth x_L",
        "RP Reco x_L",
        "figs/response_matrix_xL_RP.png",
        false,
        false,
        {0.0, 0.5},
        {0.0, 0.5}
    ));

    // x_pom correlation matrices (from x_L method)
    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_corr_B0",
        "Truth x_{pom} (1-x_L)",
        "B0 Reco x_{pom} (1-x_L)",
        "figs/response_matrix_xpom_B0.png",
        true,   // logX
        true,   // logY
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_corr_RP",
        "Truth x_{pom} (1-x_L)",
        "RP Reco x_{pom} (1-x_L)",
        "figs/response_matrix_xpom_RP.png",
        true,   // logX
        true,   // logY
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    //=================================================================
    // X_POM COMPARISON PLOTS (Definition vs from x_L)
    //=================================================================

    // x_pom comparison: 1D overlays
    plot_ptr = new PlotOptions1D(
        {"xpom_MC", "xpom_def_MC"},
        {"x_{pom} = 1 - x_{L}", "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})"},
        {"hist", "hist"},
        "MC Truth x_{pom} Comparison",
        "x_{pom}",
        "Counts",
        "figs/xpom_comparison_MC_logxy.png",
        true,   // isLogX
        true    // isLogY
    );
    plot_ptr->SetLegendPosition(0.1, 0.7, 0.5, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xpom_B0", "xpom_def_B0"},
        {"x_{pom} = 1 - x_{L}", "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})"},
        {"hist", "hist"},
        "B0 Reco x_{pom} Comparison",
        "x_{pom}",
        "Counts",
        "figs/xpom_comparison_B0_logxy.png",
        true,   // isLogX
        true    // isLogY
    );
    plot_ptr->SetLegendPosition(0.1, 0.7, 0.5, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xpom_RP", "xpom_def_RP"},
        {"x_{pom} = 1 - x_{L}", "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})"},
        {"hist", "hist"},
        "RP Reco x_{pom} Comparison",
        "x_{pom}",
        "Counts",
        "figs/xpom_comparison_RP_logxy.png",
        true,   // isLogX
        true    // isLogY
    );
    plot_ptr->SetLegendPosition(0.1, 0.7, 0.5, 0.9);
    plots.push_back(plot_ptr);


        plot_ptr = new PlotOptions1D(
        {"xpom_def_MC", "xpom_def_B0", "xpom_def_RP"},
        {"x_{pom} MC", "x_{pom} B0 Reco", "x_{pom} RP Reco"},
        {"hist", "hist", "hist"},
        "RP Reco x_{pom} Comparison",
        "x_{pom}",
        "Counts",
        "figs/xpom_comparison_all_logxy.png",
        true,   // isLogX
        true    // isLogY
    );
    plot_ptr->SetLegendPosition(0.1, 0.7, 0.5, 0.9);
    plots.push_back(plot_ptr);

    // x_pom comparison: 2D correlation plots
    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_comp_MC",
        "x_{pom} = 1 - x_{L}",
        "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})",
        "figs/xpom_2D_comparison_MC.png",
        true,   // logX
        true,   // logY
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_comp_B0",
        "x_{pom} = 1 - x_{L}",
        "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})",
        "figs/xpom_2D_comparison_B0.png",
        true,   // logX
        true,   // logY
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_comp_RP",
        "x_{pom} = 1 - x_{L}",
        "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})",
        "figs/xpom_2D_comparison_RP.png",
        true,   // logX
        true,   // logY
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    //=================================================================
    // PROTON THETA ANGLE COMPARISON (All TS protons vs B0 accepted)
    //=================================================================

    // Theta angle comparison: all truth-seeded protons vs B0 accepted
    plot_ptr = new PlotOptions1D(
        {"theta_all_TS", "theta_B0"},
        {"All Truth-Seeded Protons", "B0 Accepted Protons"},
        {"hist", "hist"},
        "Proton Scattering Angle Distribution",
        "#theta [mrad]",
        "Counts",
        "figs/theta_comparison_B0_acceptance.png",
        false,   // isLogX (logarithmic x-axis for log-binned histogram)
        false   // isLogY
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"theta_all_TS", "theta_B0"},
        {"All Truth-Seeded Protons", "B0 Accepted Protons"},
        {"hist", "hist"},
        "Proton Scattering Angle Distribution",
        "#theta [mrad]",
        "Counts",
        "figs/theta_comparison_B0_acceptance_logxy.png",
        false,   // isLogX (logarithmic x-axis for log-binned histogram)
        true    // isLogY
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    //=================================================================
    // BETA ANALYSIS PLOTS (beta = x_Bj / x_pom from definition)
    //=================================================================

    // Beta distributions with log y-axis for better visibility
    PlotOptions1D* plot_beta_logy = new PlotOptions1D(
        {"beta_MC", "beta_B0", "beta_RP"},
        {"MC Truth", "B0 Reco", "RP Reco"},
        {"hist", "hist", "hist"},
        "#beta = x_{Bj} / x_{pom} Distributions (from x_{pom} definition)",
        "#beta",
        "Counts",
        "figs/beta_distributions_logy.png",
        false,  // isLogX
        true    // isLogY
    );
    plot_beta_logy->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plot_beta_logy->SetRangeY(1, 4e4); // Set Y axis range

    plots.push_back(plot_beta_logy);

    // Beta resolution plots
    plots.push_back(new PlotOptions1D(
        {"beta_res_B0"},
        {"B0 Reco"},
        {"hist"},
        "B0 #beta Resolution",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        "Counts",
        "figs/beta_resolution_B0.png",
        false,
        false
    ));

    plots.push_back(new PlotOptions1D(
        {"beta_res_RP"},
        {"RP Reco"},
        {"hist"},
        "RP #beta Resolution",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        "Counts",
        "figs/beta_resolution_RP.png",
        false,
        false
    ));

    // Beta response matrices
    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_corr_B0",
        "Truth #beta",
        "B0 Reco #beta",
        "figs/response_matrix_beta_B0.png",
        false,  // isLogX
        false,  // isLogY
        {0.0, 1.0},
        {0.0, 1.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_corr_RP",
        "Truth #beta",
        "RP Reco #beta",
        "figs/response_matrix_beta_RP.png",
        false,  // isLogX
        false,  // isLogY
        {0.0, 1.0},
        {0.0, 1.0}
    ));

    //=================================================================
    // DIFFERENTIAL CROSS SECTION d(sigma)/dt PLOTS
    //=================================================================

    // Create sum of B0 and RP reconstructed dsigma/dt
    TH1D* h_dsigma_dt_B0_temp = (TH1D*)inputFile->Get("dsigma_dt_B0");
    TH1D* h_dsigma_dt_RP_temp = (TH1D*)inputFile->Get("dsigma_dt_RP");

    if (h_dsigma_dt_B0_temp && h_dsigma_dt_RP_temp) {
        TH1D* h_dsigma_dt_Sum = (TH1D*)h_dsigma_dt_B0_temp->Clone("dsigma_dt_Sum");
        h_dsigma_dt_Sum->SetTitle("B0+RP Sum d#sigma/dt;|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}]");
        h_dsigma_dt_Sum->Add(h_dsigma_dt_RP_temp);
        // Keep in gDirectory so plotting framework can find it
    }

    // dsigma/dt plot with linear y-axis
    plot_ptr = new PlotOptions1D(
        {"dsigma_dt_MC", "dsigma_dt_B0", "dsigma_dt_RP", "dsigma_dt_Sum"},
        {"MC Truth", "B0 Reco", "RP Reco", "B0+RP Sum"},
        {"hist", "pe", "pe", "pe"},
        "Differential Cross Section d#sigma/dt",
        "|t| [GeV^{2}]",
        "d#sigma/dt [nb/GeV^{2}]",
        "figs/dsigma_dt.png",
        true,   // isLogX
        false   // isLogY
    );
    plot_ptr->SetLegendPosition(0.6, 0.65, 0.85, 0.9);
    plots.push_back(plot_ptr);

    // dsigma/dt plot with log y-axis for better visibility
    plot_ptr = new PlotOptions1D(
        {"dsigma_dt_MC", "dsigma_dt_B0", "dsigma_dt_RP", "dsigma_dt_Sum"},
        {"MC Truth", "B0 Reco", "RP Reco", "B0+RP Sum"},
        {"hist", "pe", "pe", "pe"},
        "Differential Cross Section d#sigma/dt",
        "|t| [GeV^{2}]",
        "d#sigma/dt [nb/GeV^{2}]",
        "figs/dsigma_dt_logy.png",
        true,   // isLogX
        true    // isLogY
    );
    plot_ptr->SetLegendPosition(0.1, 0.65, 0.3, 0.9);
    plots.push_back(plot_ptr);

    //=================================================================
    // CREATE Q2 RESOLUTION VISUALIZATIONS WITH HOLLOW CIRCLES
    //=================================================================
    std::cout << "Creating Q2 resolution visualizations with hollow circles..." << std::endl;

    // Function to create hollow circle plots for Q2 resolution
    // Profile is binned in (x_Bj, y) but displayed in (x_Bj, Q2) space
    auto createResolutionCirclePlot = [](TProfile2D* profile, const char* outputName, const char* title, double s, bool showEntries = false) {
        if (!profile) {
            std::cerr << "Error: Profile histogram not found!" << std::endl;
            return;
        }

        TCanvas* c = new TCanvas("c_temp", title, 1600, 1200);
        gStyle->SetOptStat(0);

        // Create empty 2D histogram for axes (x_Bj on x-axis, Q2 on y-axis)
        TH2D* h_axes = new TH2D("h_axes_temp", title,
                                10, 0.0001, 1.0, 10, 1.0, 150.0);
        h_axes->GetXaxis()->SetTitle("x_{Bj}");
        h_axes->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
        h_axes->Draw("AXIS");

        // Define plot range
        double x_min = 0.0001;
        double x_max = 1.0;
        double q2_min = 1.0;
        double q2_max = 150.0;

        // Draw hollow circles for each bin
        // Profile is binned in (x_Bj, y) but we display at (x_Bj, Q2 = s*x*y)
        for(int i = 1; i <= profile->GetNbinsX(); ++i) {
            for(int j = 1; j <= profile->GetNbinsY(); ++j) {
                double x = profile->GetXaxis()->GetBinCenter(i);  // x_Bj
                double y = profile->GetYaxis()->GetBinCenter(j);  // inelasticity

                // Transform to Q2 space: Q2 = s * x * y
                double q2 = s * x * y;

                // Skip circles outside the plot range
                if(x < x_min || x > x_max || q2 < q2_min || q2 > q2_max) continue;

                // Get RMS (spread) of relative resolution in this bin
                double resolution_rms = profile->GetBinError(i, j);
                int entries = profile->GetBinEntries(profile->GetBin(i, j));

                // Skip empty bins
                if(resolution_rms <= 0 || entries < 5) continue;

                // Scale circle size by resolution (larger resolution = larger circle)
                // Typical Q2 resolution RMS is 0.01-0.10
                double size = 0.1 + 15.0 * resolution_rms;
                if(size > 4.0) size = 4.0; // Cap maximum size

                // Choose marker style: filled circle for < 100 entries, hollow for >= 100
                int marker_style = (entries < 100) ? 20 : 24;  // 20=filled, 24=hollow

                TMarker* marker = new TMarker(x, q2, marker_style);
                marker->SetMarkerColor(kGreen+2);  // Forest green for all
                marker->SetMarkerSize(size);
                marker->Draw("SAME");

                // Optionally draw entry count in the center of each bin
                if(showEntries) {
                    TLatex* entry_label = new TLatex();
                    entry_label->SetTextSize(0.02);
                    entry_label->SetTextAlign(22); // center-aligned
                    entry_label->SetTextColor(kBlack);
                    entry_label->DrawLatex(x, q2, Form("%d", entries));
                }
            }
        }

        c->SetLogx();
        c->SetLogy();

        // Draw combined legend for marker styles and sizes
        double legendX = 0.15;
        double legendY = 0.15;
        TLegend* legend = new TLegend(legendX, legendY, legendX + 0.18, legendY + 0.20);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.03);

        TMarker* m_hollow = new TMarker(0, 0, 24);
        m_hollow->SetMarkerColor(kGreen+2);
        m_hollow->SetMarkerSize(2.0);
        legend->AddEntry(m_hollow, "Entries #geq 100", "p");

        TMarker* m_filled = new TMarker(0, 0, 20);
        m_filled->SetMarkerColor(kGreen+2);
        m_filled->SetMarkerSize(2.0);
        legend->AddEntry(m_filled, "Entries < 100", "p");

        // Add size scale examples
        std::vector<double> example_res = {0.02, 0.05, 0.10};
        for(size_t i = 0; i < example_res.size(); ++i) {
            double res_val = example_res[i];
            double size = 0.3 + 25.0 * res_val;
            if(size > 4.0) size = 4.0;

            TMarker* m_size = new TMarker(0, 0, 24);
            m_size->SetMarkerColor(kGreen+2);
            m_size->SetMarkerSize(size);
            legend->AddEntry(m_size, Form("%.2f", res_val), "p");
        }

        legend->Draw();

        // Add ePIC simulation labels
        TLatex latex;
        latex.SetTextSize(0.04);
        latex.SetNDC();
        latex.SetTextColor(kBlack);
        latex.DrawLatex(0.2, 0.92, "#bf{ePIC} Simulation (100k events)");
        latex.DrawLatex(0.65, 0.92, "#bf{Diff. DIS} 10x100 GeV");

        c->SaveAs(outputName);
        delete c;
    };

    // Function to create best method comparison plot
    // Shows the best (smallest) resolution per bin, color-coded by method
    auto createBestMethodPlot = [](TProfile2D* prof_EM, TProfile2D* prof_DA, TProfile2D* prof_Sigma,
                                    const char* outputName, const char* title, double s, bool showEntries = false) {
        if (!prof_EM || !prof_DA || !prof_Sigma) {
            std::cerr << "Error: One or more profile histograms not found!" << std::endl;
            return;
        }

        TCanvas* c = new TCanvas("c_best", title, 1600, 1200);
        gStyle->SetOptStat(0);

        // Create empty 2D histogram for axes
        TH2D* h_axes = new TH2D("h_axes_best", title,
                                10, 0.0001, 1.0, 10, 1.0, 150.0);
        h_axes->GetXaxis()->SetTitle("x_{Bj}");
        h_axes->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
        h_axes->Draw("AXIS");

        // Define plot range
        double x_min = 0.0001;
        double x_max = 1.0;
        double q2_min = 1.0;
        double q2_max = 150.0;

        // Draw circles for best method in each bin
        for(int i = 1; i <= prof_EM->GetNbinsX(); ++i) {
            for(int j = 1; j <= prof_EM->GetNbinsY(); ++j) {
                double x = prof_EM->GetXaxis()->GetBinCenter(i);
                double y = prof_EM->GetYaxis()->GetBinCenter(j);
                double q2 = s * x * y;

                // Skip circles outside the plot range
                if(x < x_min || x > x_max || q2 < q2_min || q2 > q2_max) continue;

                // Get resolution RMS for each method
                double res_EM = prof_EM->GetBinError(i, j);
                double res_DA = prof_DA->GetBinError(i, j);
                double res_Sigma = prof_Sigma->GetBinError(i, j);

                // Get entries for each method
                int entries_EM = prof_EM->GetBinEntries(prof_EM->GetBin(i, j));
                int entries_DA = prof_DA->GetBinEntries(prof_DA->GetBin(i, j));
                int entries_Sigma = prof_Sigma->GetBinEntries(prof_Sigma->GetBin(i, j));

                // Find the method with smallest (best) resolution
                double best_res = 1e9;
                int best_method = -1;  // 0=EM, 1=DA, 2=Sigma
                int best_entries = 0;

                if(res_EM > 0 && entries_EM >= 5 && res_EM < best_res) {
                    best_res = res_EM;
                    best_method = 0;
                    best_entries = entries_EM;
                }
                if(res_DA > 0 && entries_DA >= 5 && res_DA < best_res) {
                    best_res = res_DA;
                    best_method = 1;
                    best_entries = entries_DA;
                }
                if(res_Sigma > 0 && entries_Sigma >= 5 && res_Sigma < best_res) {
                    best_res = res_Sigma;
                    best_method = 2;
                    best_entries = entries_Sigma;
                }

                // Skip if no valid method found
                if(best_method == -1) continue;

                // Scale circle size by resolution
                double size = 0.3 + 25.0 * best_res;
                if(size > 4.0) size = 4.0;

                // Choose marker style: filled circle for < 100 entries, hollow for >= 100
                int marker_style = (best_entries < 100) ? 20 : 24;  // 20=filled, 24=hollow

                // Create marker with color based on best method
                TMarker* marker = new TMarker(x, q2, marker_style);
                if(best_method == 0) {
                    marker->SetMarkerColor(kRed);      // Electron method
                } else if(best_method == 1) {
                    marker->SetMarkerColor(kBlue);     // DA method
                } else {
                    marker->SetMarkerColor(kGreen+2);  // Sigma method
                }
                marker->SetMarkerSize(size);
                marker->Draw("SAME");

                // Optionally draw entry count
                if(showEntries) {
                    TLatex* entry_label = new TLatex();
                    entry_label->SetTextSize(0.02);
                    entry_label->SetTextAlign(22);
                    entry_label->SetTextColor(kBlack);
                    entry_label->DrawLatex(x, q2, Form("%d", best_entries));
                }
            }
        }

        c->SetLogx();
        c->SetLogy();

        // Draw legend for methods
        double legendX = 0.7;
        double legendY = 0.15;
        TLegend* legend = new TLegend(legendX, legendY, legendX + 0.18, legendY + 0.12);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.03);

        TMarker* m_EM = new TMarker(0, 0, 24);
        m_EM->SetMarkerColor(kRed);
        m_EM->SetMarkerSize(2.0);
        legend->AddEntry(m_EM, "Electron Method", "p");

        TMarker* m_DA = new TMarker(0, 0, 24);
        m_DA->SetMarkerColor(kBlue);
        m_DA->SetMarkerSize(2.0);
        legend->AddEntry(m_DA, "DA Method", "p");

        TMarker* m_Sigma = new TMarker(0, 0, 24);
        m_Sigma->SetMarkerColor(kGreen+2);
        m_Sigma->SetMarkerSize(2.0);
        legend->AddEntry(m_Sigma, "Sigma Method", "p");

        legend->Draw();

        // Draw combined legend for marker styles and sizes
        double legendX2 = 0.15;
        double legendY2 = 0.15;
        TLegend* legend2 = new TLegend(legendX2, legendY2, legendX2 + 0.18, legendY2 + 0.20);
        legend2->SetFillStyle(0);
        legend2->SetBorderSize(0);
        legend2->SetTextSize(0.03);

        TMarker* m_hollow = new TMarker(0, 0, 24);
        m_hollow->SetMarkerColor(kBlack);
        m_hollow->SetMarkerSize(2.0);
        legend2->AddEntry(m_hollow, "Entries #geq 100", "p");

        TMarker* m_filled = new TMarker(0, 0, 20);
        m_filled->SetMarkerColor(kBlack);
        m_filled->SetMarkerSize(2.0);
        legend2->AddEntry(m_filled, "Entries < 100", "p");

        // Add size scale examples
        std::vector<double> example_res = {0.02, 0.05, 0.10};
        for(size_t i = 0; i < example_res.size(); ++i) {
            double res_val = example_res[i];
            double size = 0.3 + 25.0 * res_val;
            if(size > 4.0) size = 4.0;

            TMarker* m_size = new TMarker(0, 0, 24);
            m_size->SetMarkerColor(kBlack);
            m_size->SetMarkerSize(size);
            legend2->AddEntry(m_size, Form("%.2f", res_val), "p");
        }

        legend2->Draw();

        // Add ePIC simulation labels
        TLatex latex;
        latex.SetTextSize(0.04);
        latex.SetNDC();
        latex.SetTextColor(kBlack);
        latex.DrawLatex(0.2, 0.92, "#bf{ePIC} Simulation (100k events)");
        latex.DrawLatex(0.65, 0.92, "#bf{Diff. DIS} 10x100 GeV");

        c->SaveAs(outputName);
        delete c;
    };

    // Create Q2 resolution circle plots
    // Histograms are binned in (x_Bj, y) but displayed in (x_Bj, Q2) space
    // s = 4*E_e*E_p (center-of-mass energy squared)
    // For 18x275 GeV: s = 19800 GeV²
    // For 10x100 GeV: s = 4000 GeV²
    // For 5x41 GeV: s = 820 GeV²
    double s = 4000.0; // Default

    TProfile2D* prof_EM = (TProfile2D*)inputFile->Get("Q2_RelRes_vs_xy_EM");
    TProfile2D* prof_DA = (TProfile2D*)inputFile->Get("Q2_RelRes_vs_xy_DA");
    TProfile2D* prof_Sigma = (TProfile2D*)inputFile->Get("Q2_RelRes_vs_xy_Sigma");

    if (prof_EM) {
        createResolutionCirclePlot(prof_EM, "figs/Q2_RelRes_Q2x_EM.png",
                                  "Q^{2} Relative Resolution vs x_{Bj} and Q^{2} (Electron Method)", s);
    }
    if (prof_DA) {
        createResolutionCirclePlot(prof_DA, "figs/Q2_RelRes_Q2x_DA.png",
                                  "Q^{2} Relative Resolution vs x_{Bj} and Q^{2} (Double Angle Method)", s);
    }
    if (prof_Sigma) {
        createResolutionCirclePlot(prof_Sigma, "figs/Q2_RelRes_Q2x_Sigma.png",
                                  "Q^{2} Relative Resolution vs x_{Bj} and Q^{2} (Sigma Method)", s);
    }

    // Create best method comparison plot
    if (prof_EM && prof_DA && prof_Sigma) {
        std::cout << "Creating best method comparison plot..." << std::endl;
        createBestMethodPlot(prof_EM, prof_DA, prof_Sigma,
                            "figs/Q2_RelRes_Q2x_BestMethod.png",
                            "Q^{2} Relative Resolution - Best Method per Bin", s);
    }

    //=================================================================
    // TRIPLE DIFFERENTIAL CROSS SECTION d³σ/(dQ² dβ dx_pom) PLOTS
    //=================================================================
    std::cout << "Creating triple differential cross section 1D slices..." << std::endl;

    // Load the 3D histograms
    TH3D* h_d3sigma_MC = (TH3D*)inputFile->Get("d3sigma_dQ2dbeta_dxpom_MC");
    TH3D* h_d3sigma_B0 = (TH3D*)inputFile->Get("d3sigma_dQ2dbeta_dxpom_B0");
    TH3D* h_d3sigma_RP = (TH3D*)inputFile->Get("d3sigma_dQ2dbeta_dxpom_RP");

    if (h_d3sigma_MC && h_d3sigma_B0 && h_d3sigma_RP) {
        int n_Q2_bins = h_d3sigma_MC->GetNbinsX();
        int n_beta_bins = h_d3sigma_MC->GetNbinsY();
        int n_xpom_bins = h_d3sigma_MC->GetNbinsZ();

        //=================================================================
        // PLOT 1: d³σ vs β for each (Q², x_pom) bin
        // Grid: columns = Q² bins, rows = x_pom bins
        // Compare MC (black hist), B0 (red points), RP (blue points), Sum (orange points)
        //=================================================================
        std::cout << "Creating d3sigma vs beta comparison plots..." << std::endl;

        // Create canvas with grid layout: n_xpom_bins rows x n_Q2_bins columns
        TCanvas* c_beta = new TCanvas("c_d3sigma_vs_beta", "d3sigma vs beta", 400*n_Q2_bins, 300*n_xpom_bins);
        c_beta->Divide(n_Q2_bins, n_xpom_bins);

        for (int xpom_bin = 1; xpom_bin <= n_xpom_bins; xpom_bin++) {
            for (int Q2_bin = 1; Q2_bin <= n_Q2_bins; Q2_bin++) {
                int pad_num = (xpom_bin - 1) * n_Q2_bins + Q2_bin;

                double Q2_low = h_d3sigma_MC->GetXaxis()->GetBinLowEdge(Q2_bin);
                double Q2_high = h_d3sigma_MC->GetXaxis()->GetBinUpEdge(Q2_bin);
                double xpom_low = h_d3sigma_MC->GetZaxis()->GetBinLowEdge(xpom_bin);
                double xpom_high = h_d3sigma_MC->GetZaxis()->GetBinUpEdge(xpom_bin);

                // Create 1D histograms for this (Q², x_pom) slice
                TH1D* h_slice_MC = new TH1D(TString::Format("slice_beta_MC_Q%d_xp%d", Q2_bin, xpom_bin),
                    TString::Format("Q^{2}=[%.1f,%.1f], x_{pom}=[%.3f,%.3f]", Q2_low, Q2_high, xpom_low, xpom_high),
                    n_beta_bins, h_d3sigma_MC->GetYaxis()->GetXmin(), h_d3sigma_MC->GetYaxis()->GetXmax());

                TH1D* h_slice_B0 = new TH1D(TString::Format("slice_beta_B0_Q%d_xp%d", Q2_bin, xpom_bin),
                    TString::Format("Q^{2}=[%.1f,%.1f], x_{pom}=[%.3f,%.3f]", Q2_low, Q2_high, xpom_low, xpom_high),
                    n_beta_bins, h_d3sigma_B0->GetYaxis()->GetXmin(), h_d3sigma_B0->GetYaxis()->GetXmax());

                TH1D* h_slice_RP = new TH1D(TString::Format("slice_beta_RP_Q%d_xp%d", Q2_bin, xpom_bin),
                    TString::Format("Q^{2}=[%.1f,%.1f], x_{pom}=[%.3f,%.3f]", Q2_low, Q2_high, xpom_low, xpom_high),
                    n_beta_bins, h_d3sigma_RP->GetYaxis()->GetXmin(), h_d3sigma_RP->GetYaxis()->GetXmax());

                TH1D* h_slice_Sum = new TH1D(TString::Format("slice_beta_Sum_Q%d_xp%d", Q2_bin, xpom_bin),
                    TString::Format("Q^{2}=[%.1f,%.1f], x_{pom}=[%.3f,%.3f]", Q2_low, Q2_high, xpom_low, xpom_high),
                    n_beta_bins, h_d3sigma_MC->GetYaxis()->GetXmin(), h_d3sigma_MC->GetYaxis()->GetXmax());

                // Fill the 1D histograms
                for (int beta_bin = 1; beta_bin <= n_beta_bins; beta_bin++) {
                    double value_MC = h_d3sigma_MC->GetBinContent(Q2_bin, beta_bin, xpom_bin);
                    double error_MC = h_d3sigma_MC->GetBinError(Q2_bin, beta_bin, xpom_bin);
                    h_slice_MC->SetBinContent(beta_bin, value_MC);
                    h_slice_MC->SetBinError(beta_bin, error_MC);

                    double value_B0 = h_d3sigma_B0->GetBinContent(Q2_bin, beta_bin, xpom_bin);
                    double error_B0 = h_d3sigma_B0->GetBinError(Q2_bin, beta_bin, xpom_bin);
                    h_slice_B0->SetBinContent(beta_bin, value_B0);
                    h_slice_B0->SetBinError(beta_bin, error_B0);

                    double value_RP = h_d3sigma_RP->GetBinContent(Q2_bin, beta_bin, xpom_bin);
                    double error_RP = h_d3sigma_RP->GetBinError(Q2_bin, beta_bin, xpom_bin);
                    h_slice_RP->SetBinContent(beta_bin, value_RP);
                    h_slice_RP->SetBinError(beta_bin, error_RP);

                    // Create sum of B0 and RP
                    double value_Sum = value_B0 + value_RP;
                    double error_Sum = TMath::Sqrt(error_B0*error_B0 + error_RP*error_RP);
                    h_slice_Sum->SetBinContent(beta_bin, value_Sum);
                    h_slice_Sum->SetBinError(beta_bin, error_Sum);
                }

                // Set axis titles
                h_slice_MC->GetXaxis()->SetTitle("#beta");
                h_slice_MC->GetYaxis()->SetTitle("d^{3}#sigma/(dQ^{2}d#betadx_{pom}) [nb/GeV^{2}]");

                // Set styles: MC = black histogram
                h_slice_MC->SetLineColor(kBlack);
                h_slice_MC->SetLineWidth(2);

                // B0 = red points
                h_slice_B0->SetMarkerStyle(20);
                h_slice_B0->SetMarkerColor(kRed);
                h_slice_B0->SetLineColor(kRed);
                h_slice_B0->SetMarkerSize(0.8);

                // RP = blue points
                h_slice_RP->SetMarkerStyle(20);
                h_slice_RP->SetMarkerColor(kBlue);
                h_slice_RP->SetLineColor(kBlue);
                h_slice_RP->SetMarkerSize(0.8);

                // Sum = orange points
                h_slice_Sum->SetMarkerStyle(20);
                h_slice_Sum->SetMarkerColor(kOrange+7);
                h_slice_Sum->SetLineColor(kOrange+7);
                h_slice_Sum->SetMarkerSize(0.8);

                // Draw all on the same pad
                c_beta->cd(pad_num);
                gPad->SetLogy();
                h_slice_MC->Draw("HIST");
                h_slice_B0->Draw("PE SAME");
                h_slice_RP->Draw("PE SAME");
                h_slice_Sum->Draw("PE SAME");

                // Add legend only in the first pad
                if (pad_num == 1) {
                    TLegend* leg = new TLegend(0.5, 0.65, 0.88, 0.88);
                    leg->AddEntry(h_slice_MC, "MC Truth", "l");
                    leg->AddEntry(h_slice_B0, "B0 Reco", "pe");
                    leg->AddEntry(h_slice_RP, "RP Reco", "pe");
                    leg->AddEntry(h_slice_Sum, "B0+RP Sum", "pe");
                    leg->SetTextSize(0.03);
                    leg->Draw();
                }
            }
        }

        c_beta->SaveAs("figs/d3sigma_vs_beta.png");
        delete c_beta;

        //=================================================================
        // PLOT 2: d³σ vs x_pom for each (Q², β) bin
        // Grid: columns = Q² bins, rows = β bins
        // Compare MC (black hist), B0 (red points), RP (blue points), Sum (orange points)
        //=================================================================
        std::cout << "Creating d3sigma vs x_pom comparison plots..." << std::endl;

        TCanvas* c_xpom = new TCanvas("c_d3sigma_vs_xpom", "d3sigma vs xpom", 400*n_Q2_bins, 300*n_beta_bins);
        c_xpom->Divide(n_Q2_bins, n_beta_bins);

        for (int beta_bin = 1; beta_bin <= n_beta_bins; beta_bin++) {
            for (int Q2_bin = 1; Q2_bin <= n_Q2_bins; Q2_bin++) {
                int pad_num = (beta_bin - 1) * n_Q2_bins + Q2_bin;

                double Q2_low = h_d3sigma_MC->GetXaxis()->GetBinLowEdge(Q2_bin);
                double Q2_high = h_d3sigma_MC->GetXaxis()->GetBinUpEdge(Q2_bin);
                double beta_low = h_d3sigma_MC->GetYaxis()->GetBinLowEdge(beta_bin);
                double beta_high = h_d3sigma_MC->GetYaxis()->GetBinUpEdge(beta_bin);

                // Create 1D histograms for this (Q², β) slice - with log binning
                std::vector<double> xpom_edges;
                for (int i = 0; i <= n_xpom_bins; i++) {
                    xpom_edges.push_back(h_d3sigma_MC->GetZaxis()->GetBinLowEdge(i+1));
                }

                TH1D* h_slice_MC = new TH1D(TString::Format("slice_xpom_MC_Q%d_b%d", Q2_bin, beta_bin),
                    TString::Format("Q^{2}=[%.1f,%.1f], #beta=[%.2f,%.2f]", Q2_low, Q2_high, beta_low, beta_high),
                    n_xpom_bins, xpom_edges.data());

                TH1D* h_slice_B0 = new TH1D(TString::Format("slice_xpom_B0_Q%d_b%d", Q2_bin, beta_bin),
                    TString::Format("Q^{2}=[%.1f,%.1f], #beta=[%.2f,%.2f]", Q2_low, Q2_high, beta_low, beta_high),
                    n_xpom_bins, xpom_edges.data());

                TH1D* h_slice_RP = new TH1D(TString::Format("slice_xpom_RP_Q%d_b%d", Q2_bin, beta_bin),
                    TString::Format("Q^{2}=[%.1f,%.1f], #beta=[%.2f,%.2f]", Q2_low, Q2_high, beta_low, beta_high),
                    n_xpom_bins, xpom_edges.data());

                TH1D* h_slice_Sum = new TH1D(TString::Format("slice_xpom_Sum_Q%d_b%d", Q2_bin, beta_bin),
                    TString::Format("Q^{2}=[%.1f,%.1f], #beta=[%.2f,%.2f]", Q2_low, Q2_high, beta_low, beta_high),
                    n_xpom_bins, xpom_edges.data());

                // Fill the 1D histograms
                for (int xpom_bin = 1; xpom_bin <= n_xpom_bins; xpom_bin++) {
                    double value_MC = h_d3sigma_MC->GetBinContent(Q2_bin, beta_bin, xpom_bin);
                    double error_MC = h_d3sigma_MC->GetBinError(Q2_bin, beta_bin, xpom_bin);
                    h_slice_MC->SetBinContent(xpom_bin, value_MC);
                    h_slice_MC->SetBinError(xpom_bin, error_MC);

                    double value_B0 = h_d3sigma_B0->GetBinContent(Q2_bin, beta_bin, xpom_bin);
                    double error_B0 = h_d3sigma_B0->GetBinError(Q2_bin, beta_bin, xpom_bin);
                    h_slice_B0->SetBinContent(xpom_bin, value_B0);
                    h_slice_B0->SetBinError(xpom_bin, error_B0);

                    double value_RP = h_d3sigma_RP->GetBinContent(Q2_bin, beta_bin, xpom_bin);
                    double error_RP = h_d3sigma_RP->GetBinError(Q2_bin, beta_bin, xpom_bin);
                    h_slice_RP->SetBinContent(xpom_bin, value_RP);
                    h_slice_RP->SetBinError(xpom_bin, error_RP);

                    // Create sum of B0 and RP
                    double value_Sum = value_B0 + value_RP;
                    double error_Sum = TMath::Sqrt(error_B0*error_B0 + error_RP*error_RP);
                    h_slice_Sum->SetBinContent(xpom_bin, value_Sum);
                    h_slice_Sum->SetBinError(xpom_bin, error_Sum);
                }

                // Set axis titles
                h_slice_MC->GetXaxis()->SetTitle("x_{pom}");
                h_slice_MC->GetYaxis()->SetTitle("d^{3}#sigma/(dQ^{2}d#betadx_{pom}) [nb/GeV^{2}]");

                // Set styles: MC = black histogram
                h_slice_MC->SetLineColor(kBlack);
                h_slice_MC->SetLineWidth(2);

                // B0 = red points
                h_slice_B0->SetMarkerStyle(20);
                h_slice_B0->SetMarkerColor(kRed);
                h_slice_B0->SetLineColor(kRed);
                h_slice_B0->SetMarkerSize(0.8);

                // RP = blue points
                h_slice_RP->SetMarkerStyle(20);
                h_slice_RP->SetMarkerColor(kBlue);
                h_slice_RP->SetLineColor(kBlue);
                h_slice_RP->SetMarkerSize(0.8);

                // Sum = orange points
                h_slice_Sum->SetMarkerStyle(20);
                h_slice_Sum->SetMarkerColor(kOrange+7);
                h_slice_Sum->SetLineColor(kOrange+7);
                h_slice_Sum->SetMarkerSize(0.8);

                // Draw all on the same pad
                c_xpom->cd(pad_num);
                gPad->SetLogx();
                gPad->SetLogy();
                h_slice_MC->Draw("HIST");
                h_slice_B0->Draw("PE SAME");
                h_slice_RP->Draw("PE SAME");
                h_slice_Sum->Draw("PE SAME");

                // Add legend only in the first pad
                if (pad_num == 1) {
                    TLegend* leg = new TLegend(0.5, 0.65, 0.88, 0.88);
                    leg->AddEntry(h_slice_MC, "MC Truth", "l");
                    leg->AddEntry(h_slice_B0, "B0 Reco", "pe");
                    leg->AddEntry(h_slice_RP, "RP Reco", "pe");
                    leg->AddEntry(h_slice_Sum, "B0+RP Sum", "pe");
                    leg->SetTextSize(0.03);
                    leg->Draw();
                }
            }
        }

        c_xpom->SaveAs("figs/d3sigma_vs_xpom.png");
        delete c_xpom;

        std::cout << "Triple differential cross section 1D slice plots created!" << std::endl;
    } else {
        std::cout << "Warning: Could not find d3sigma histograms" << std::endl;
    }

    //=================================================================
    // PLOT ALL CONFIGURATIONS
    //=================================================================

    gSystem->mkdir("figs", kTRUE);

    for (const auto& plot : plots) {
        plot->Plot(inputFile);
    }
    
    for (const auto& plot : plots) {
        delete plot;
    }

    inputFile->Close();
    delete inputFile;

    std::cout << "\nPlotting complete! All plots saved to figs/ directory" << std::endl;
    std::cout << "\nKey plots created:" << std::endl;
    std::cout << "  Q2/xy analysis: response matrices, E-pz, eta_max, M_X^2" << std::endl;
    std::cout << "  Q2 resolution: hollow circle plots binned in (x_Bj, y), displayed in (x_Bj, Q2) plane" << std::endl;
    std::cout << "                 Creates parallelogram pattern due to Q2 = s*x*y transformation" << std::endl;
    std::cout << "                 Best method comparison plot: color-coded by best resolution method" << std::endl;
    std::cout << "                 (Red=Electron, Blue=DA, Green=Sigma)" << std::endl;
    std::cout << "  Mandelstam t analysis: t, x_L, x_pom response matrices" << std::endl;
    std::cout << "  x_pom comparison: Definition vs (1-x_L) for MC, B0, RP" << std::endl;
    std::cout << "  Beta analysis: distributions, resolutions, response matrices (MC, B0, RP)" << std::endl;
    std::cout << "  Differential cross section: d(sigma)/dt for MC, B0, RP, and B0+RP sum (linear and log-y)" << std::endl;
    std::cout << "                              Reconstructed data shown as points with error bars" << std::endl;
    std::cout << "  Triple differential cross section: d^3(sigma)/(dQ^2 dbeta dx_pom) 1D slices" << std::endl;
    std::cout << "                                      - vs beta (grid: Q2 columns x x_pom rows)" << std::endl;
    std::cout << "                                      - vs x_pom (grid: Q2 columns x beta rows)" << std::endl;
    std::cout << "                                      MC (black hist), B0 (red pts), RP (blue pts), Sum (orange pts)" << std::endl;

    return 0;
}
