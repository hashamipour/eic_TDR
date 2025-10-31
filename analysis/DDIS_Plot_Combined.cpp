//g++ -o DDIS_Plot_Combined $(root-config --cflags --glibs) DDIS_Plot_Combined.cpp Plotting.cpp

#include "Plotting.hpp"
#include <TFile.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include "TError.h"

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

    // Response matrices for Q2
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
        "Corr_Q2_ESigma",
        "Q^{2} (true) [GeV]",
        "Q^{2} (ESigma) [GeV]",
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
        "x_Corr_ESigma",
        "x_{Bj} (true)",
        "x_{Bj} (ESigma)",
        "figs/response_matrix_x_ESigma.png",
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
    plot_ptr->SetLegendPosition(0.15, 0.7, 0.55, 0.9);
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
    plot_ptr->SetLegendPosition(0.15, 0.7, 0.55, 0.9);
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
    plot_ptr->SetLegendPosition(0.15, 0.7, 0.55, 0.9);
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
    std::cout << "  Mandelstam t analysis: t, x_L, x_pom response matrices" << std::endl;
    std::cout << "  x_pom comparison: Definition vs (1-x_L) for MC, B0, RP" << std::endl;

    return 0;
}
