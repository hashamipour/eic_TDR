//g++ -o DDIS_Plot_t $(root-config --cflags --glibs) DDIS_Plot_t.cpp Plotting.cpp Utility.cpp

#include "Plotting.hpp"
#include "Utility.hpp"
#include <TFile.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root>" << std::endl;
        return 1;
    }
    
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

    // 1D histogram plot configuration for Mandelstam t
    plots.push_back(new PlotOptions1D(
        {"t_MC", "t_B0", "t_RP_histo"},
        {"MC Truth", "Reco B0 (5.5-20 mrad)", "Reco RP (<5 mrad)"},
        {"hist", "pe", "pe"},
        "Mandelstam |t| Distributions",
        "|t| [GeV^{2}]",
        "Counts",
        "figs/t_distributions.png",
        false,
        false
    ));

    // Normalized PDF comparison
    plots.push_back(new PlotOptions1D(
        {"t_MC", "t_B0", "t_RP_histo"},
        {"MC Truth", "Reco B0 (5.5-20 mrad)", "Reco RP (<5 mrad)"},
        {"hist", "pe", "pe"},
        "Mandelstam |t| PDF Comparison",
        "|t| [GeV^{2}]",
        "PDF",
        "figs/t_pdf_comparison.png",
        true,
        true,
        true
    ));

    // Angular distributions
    plots.push_back(new PlotOptions1D(
        {"theta_MC", "theta_B0", "theta_RP"},
        {"MC Truth", "Reco B0", "Reco RP"},
        {"hist", "pe", "pe"},
        "Proton Scattering Angles",
        "#theta [mrad]",
        "Counts",
        "figs/theta_distributions.png",
        false,
        true
    ));

    // Resolution plots
    plots.push_back(new PlotOptionsRelRes(
        "t_res_B0",
        "(|t|_{reco} - |t|_{truth})/|t|_{truth}",
        "Counts",
        -999., -999.,
        "figs/t_resolution_B0.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "t_res_RP",
        "(|t|_{reco} - |t|_{truth})/|t|_{truth}",
        "Counts",
        -999., -999.,
        "figs/t_resolution_RP.png"
    ));

    // Combined correlation plot using TGraphs (no more blocky appearance!)
    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"t_corr_B0_graph", "t_corr_RP_graph"},  // TGraph names
        {"B0 Protons (5.5-20 mrad)", "Roman Pot (<5 mrad)"},
        {kBlue, kCyan+1},
        {20, 21},
        "Combined Truth vs Reco |t| Correlation",
        "Truth |t| [GeV^{2}]",
        "Reco |t| [GeV^{2}]",
        "figs/t_correlation_combined.png",
        {0.0, 2.0},
        {0.0, 2.0}
    ));

    // Response matrices (keep TH2D for these)
    plots.push_back(new PlotOptionsResponseMatrix(
        "t_corr_B0",
        "Truth |t| [GeV^{2}]",
        "B0 Reco |t| [GeV^{2}]",
        "figs/response_matrix_B0.png",
        false,
        false,
        {0.0, 2.0},
        {0.0, 2.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "t_corr_RP",
        "Truth |t| [GeV^{2}]",
        "RP Reco |t| [GeV^{2}]",
        "figs/response_matrix_RP.png",
        false,
        false,
        {0.0, 0.5},
        {0.0, 0.5}
    ));

    // After the RP resolution plot (~line 70), add:

    // eX method plots
    plots.push_back(new PlotOptions1D(
        {"t_MC", "t_eX"},
        {"MC Truth", "eX Method"},
        {"hist", "pe"},
        "eX Method |t| Distribution",
        "|t| [GeV^{2}]",
        "Counts",
        "figs/t_distribution_eX.png",
        true,
        true
    ));

    plots.push_back(new PlotOptionsRelRes(
        "t_res_eX",
        "(|t|_{reco} - |t|_{truth})/|t|_{truth}",
        "Counts",
        -999., -999.,
        "figs/t_resolution_eX.png"
    ));

    // correlation eX:
    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"t_corr_eX_graph"},  // Add eX
        {"eX Method"},
        {kGreen+2},  // Add green
        {22},
        "Combined Truth vs Reco |t| Correlation",
        "Truth |t| [GeV^{2}]",
        "Reco |t| [GeV^{2}]",
        "figs/t_correlation_combined_eX.png",
        {0.0, 2.0},
        {0.0, 2.0}
    ));


    /// debug
    plots.push_back(new PlotOptions1D(
        {"Q2_EICRecon", "Q2_calc"},
        {"EICRecon", "Our Calculation"},
        {"hist", "hist"},
        "Q^{2} Comparison",
        "Q^{2} [GeV^{2}]",
        "Counts",
        "figs/Q2_comparison.png",
        false, true
    ));

    gSystem->mkdir("figs", kTRUE);

    for (const auto& plot : plots) {
        plot->Plot(inputFile);
    }
    
    for (const auto& plot : plots) {
        delete plot;
    }

    inputFile->Close();
    delete inputFile;

    std::cout << "\nPlotting complete! All plots saved to figs/ directory:" << std::endl;
    std::cout << "- t_distributions.png" << std::endl;
    std::cout << "- t_pdf_comparison.png" << std::endl;
    std::cout << "- theta_distributions.png" << std::endl;
    std::cout << "- t_resolution_B0.png" << std::endl;
    std::cout << "- t_resolution_RP.png" << std::endl;
    std::cout << "- t_correlation_combined.png (using TGraphs - smooth!)" << std::endl;
    std::cout << "- response_matrix_B0.png" << std::endl;
    std::cout << "- response_matrix_RP.png" << std::endl;

    return 0;
}