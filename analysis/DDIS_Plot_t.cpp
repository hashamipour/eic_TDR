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
    // plots.push_back(new PlotOptions1D(
    //     {"t_MC", "t_B0", "t_RP_histo"},
    //     {"MC Truth", "Reco B0 (5.5-20 mrad)", "Reco RP (<5 mrad)"},
    //     {"hist", "pe", "pe"},
    //     "Mandelstam |t| Distributions",
    //     "|t| [GeV^{2}]",
    //     "Counts",
    //     "figs/t_distributions.png",
    //     false,
    //     false
    // ));

    // PlotOptions1D* plot_ptr = new PlotOptions1D(
    //     {"t_MC", "t_B0", "t_RP_histo"},
    //     {"MC Truth", "Reco B0 (5.5-20 mrad)", "Reco RP (<5 mrad)"},
    //     {"hist", "pe", "pe"},
    //     "Mandelstam |t| Distributions",
    //     "|t| [GeV^{2}]",
    //     "Counts",
    //     "figs/t_distributions_logy.png",
    //     false,
    //     true // logY
    // );
    // plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    // plots.push_back(plot_ptr);

//    // Normalized PDF comparison
    // plots.push_back(new PlotOptions1D(
    //     {"t_MC", "t_B0", "t_RP_histo"},
    //     {"MC Truth", "Reco B0 (5.5-20 mrad)", "Reco RP (<5 mrad)"},
    //     {"hist", "pe", "pe"},
    //     "Mandelstam |t| PDF Comparison",
    //     "|t| [GeV^{2}]",
    //     "PDF",
    //     "figs/t_pdf_comparison.png",
    //     true,
    //     true,
    //     true
    // ));

  //  // Angular distributions
    // plots.push_back(new PlotOptions1D(
    //     {"theta_MC", "theta_B0", "theta_RP"},
    //     {"MC Truth", "Reco B0", "Reco RP"},
    //     {"hist", "pe", "pe"},
    //     "Proton Scattering Angles",
    //     "#theta [mrad]",
    //     "Counts",
    //     "figs/theta_distributions.png",
    //     false,
    //     true
    // ));

    // // Resolution plots
    // plots.push_back(new PlotOptionsRelRes(
    //     "t_res_B0",
    //     "(|t|_{reco} - |t|_{truth})/|t|_{truth}",
    //     "Counts",
    //     -999., -999.,
    //     "figs/t_resolution_B0.png"
    // ));

    // plots.push_back(new PlotOptionsRelRes(
    //     "t_res_RP",
    //     "(|t|_{reco} - |t|_{truth})/|t|_{truth}",
    //     "Counts",
    //     -999., -999.,
    //     "figs/t_resolution_RP.png"
    // ));

    // // x_L resolution plots
    // plots.push_back(new PlotOptionsRelRes(
    //     "xL_res_B0",
    //     "(x_{L,reco} - x_{L,truth})/x_{L,truth}",
    //     "Counts",
    //     -999., -999.,
    //     "figs/xL_resolution_B0.png"
    // ));

    // plots.push_back(new PlotOptionsRelRes(
    //     "xL_res_RP",
    //     "(x_{L,reco} - x_{L,truth})/x_{L,truth}",
    //     "Counts",
    //     -999., -999.,
    //     "figs/xL_resolution_RP.png"
    // ));

    // // x_pom resolution plots
    // plots.push_back(new PlotOptionsRelRes(
    //     "xpom_res_B0",
    //     "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
    //     "Counts",
    //     -999., -999.,
    //     "figs/xpom_resolution_B0.png"
    // ));

    // plots.push_back(new PlotOptionsRelRes(
    //     "xpom_res_RP",
    //     "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
    //     "Counts",
    //     -999., -999.,
    //     "figs/xpom_resolution_RP.png"
    // ));

    // Combined correlation plot using TGraphs (no more blocky appearance!)
    // plots.push_back(new PlotOptionsCombinedCorrelation(
    //     {"t_corr_B0_graph", "t_corr_RP_graph"},  // TGraph names
    //     {"B0 Protons (5.5-20 mrad)", "Roman Pot (<5 mrad)"},
    //     {kBlue, kCyan+1},
    //     {20, 21},
    //     "Combined Truth vs Reco |t| Correlation",
    //     "Truth |t| [GeV^{2}]",
    //     "Reco |t| [GeV^{2}]",
    //     "figs/t_correlation_combined.png",
    //     {0.0, 2.0},
    //     {0.0, 2.0}
    // ));

    // Response matrices (keep TH2D for these)
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

    // // debug
    // plots.push_back(new PlotOptions1D(
    //     {"Q2_EICRecon", "Q2_calc"},
    //     {"EICRecon", "Our Calculation"},
    //     {"hist", "hist"},
    //     "Q^{2} Comparison",
    //     "Q^{2} [GeV^{2}]",
    //     "Counts",
    //     "figs/Q2_comparison.png",
    //     false, true
    // ));

    // plot_ptr = new PlotOptions1D(
    //     {"xL_MC", "xL_B0", "xL_RP"},//, "xL_eX"},
    //     {"Truth", "B0", "Roman Pot", "eX"},
    //     {"hist", "hist", "hist", "hist"},
    //     "x_{L}",
    //     "x_{L}",
    //     "Counts",
    //     "figs/x_L_comparison_logy.png",
    //     false,  // isLogX
    //     true,   // isLogY
    //     false   // normalize
    // );
    // plot_ptr->SetLegendPosition(0.15, 0.7, 0.35, 0.9);
    // plots.push_back(plot_ptr);


    // plot_ptr = new PlotOptions1D(
    //     {"xL_MC", "xL_B0", "xL_RP"},//, "xL_eX"},
    //     {"Truth", "B0", "Roman Pot", "eX"},
    //     {"hist", "hist", "hist", "hist"},
    //     "x_{L}",
    //     "x_{L}",
    //     "Counts",
    //     "figs/x_L_comparison.png",
    //     false,  // isLogX
    //     false,   // isLogY
    //     false   // normalize
    // );
    // plot_ptr->SetLegendPosition(0.15, 0.7, 0.35, 0.9);
    // plots.push_back(plot_ptr);

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

    // x_pom plots with log x and log y
    // PlotOptions1D* plot_xpom_logy = new PlotOptions1D(
    //     {"xpom_MC", "xpom_B0", "xpom_RP"},
    //     {"Truth", "B0", "Roman Pot"},
    //     {"hist", "hist", "hist"},
    //     "x_{pom}",
    //     "x_{pom}",
    //     "Counts",
    //     "figs/xpom_comparison_logxy.png",
    //     true,   // isLogX
    //     true,   // isLogY
    //     false   // normalize
    // );
    // plot_xpom_logy->SetLegendPosition(0.15, 0.7, 0.35, 0.9);
    // plots.push_back(plot_xpom_logy);

    // // x_pom plots with log x only
    // PlotOptions1D* plot_xpom = new PlotOptions1D(
    //     {"xpom_MC", "xpom_B0", "xpom_RP"},
    //     {"Truth", "B0", "Roman Pot"},
    //     {"hist", "hist", "hist"},
    //     "x_{pom}",
    //     "x_{pom}",
    //     "Counts",
    //     "figs/xpom_comparison_logx.png",
    //     true,   // isLogX
    //     false,  // isLogY
    //     false   // normalize
    // );
    // plot_xpom->SetLegendPosition(0.15, 0.7, 0.35, 0.9);
    // plots.push_back(plot_xpom);

    // x_pom correlation matrices
    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_corr_B0",
        "Truth x_{pom}",
        "B0 Reco x_{pom}",
        "figs/response_matrix_xpom_B0.png",
        true,   // logX
        true,   // logY
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_corr_RP",
        "Truth x_{pom}",
        "RP Reco x_{pom}",
        "figs/response_matrix_xpom_RP.png",
        true,   // logX
        true,   // logY
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    //////////////////////////////////////////////////////////////////////
    // B0 Comparison: Original vs First Bin Cut
    //////////////////////////////////////////////////////////////////////

    // |t| distribution: B0 vs B0_cutFirstBin
    // plot_ptr = new PlotOptions1D(
    //     {"t_B0", "t_B0_cutFirstBin"},
    //     {"B0 (all bins)", "B0 (first bin cut)"},
    //     {"hist", "hist"},
    //     "B0 |t| Distribution: Effect of First Bin Cut",
    //     "|t| [GeV^{2}]",
    //     "Counts",
    //     "figs/t_B0_comparison_firstBinCut.png",
    //     false,  // isLogX
    //     true    // isLogY
    // );
    // plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    // plots.push_back(plot_ptr);

    // x_L distribution: B0 vs B0_cutFirstBin
    // plot_ptr = new PlotOptions1D(
    //     {"xL_B0", "xL_B0_cutFirstBin"},
    //     {"B0 (all bins)", "B0 (first bin cut)"},
    //     {"hist", "hist"},
    //     "B0 x_{L} Distribution: Effect of First Bin Cut",
    //     "x_{L}",
    //     "Counts",
    //     "figs/xL_B0_comparison_firstBinCut.png",
    //     false,  // isLogX
    //     true    // isLogY
    // );
    // plot_ptr->SetLegendPosition(0.15, 0.7, 0.4, 0.9);
    // plots.push_back(plot_ptr);

    // x_pom distribution: B0 vs B0_cutFirstBin
    // plot_ptr = new PlotOptions1D(
    //     {"xpom_B0", "xpom_B0_cutFirstBin"},
    //     {"B0 (all bins)", "B0 (first bin cut)"},
    //     {"hist", "hist"},
    //     "B0 x_{pom} Distribution: Effect of First Bin Cut",
    //     "x_{pom}",
    //     "Counts",
    //     "figs/xpom_B0_comparison_firstBinCut.png",
    //     true,   // isLogX
    //     true    // isLogY
    // );
    // plot_ptr->SetLegendPosition(0.15, 0.7, 0.4, 0.9);
    // plots.push_back(plot_ptr);

    // Response matrices for B0_cutFirstBin
    // plots.push_back(new PlotOptionsResponseMatrix(
    //     "t_corr_B0_cutFirstBin",
    //     "Truth |t| [GeV^{2}]",
    //     "B0 Reco |t| [GeV^{2}] (first bin cut)",
    //     "figs/response_matrix_B0_cutFirstBin.png",
    //     false,
    //     false,
    //     {0.0, 2.0},
    //     {0.0, 2.0}
    // ));

    // plots.push_back(new PlotOptionsResponseMatrix(
    //     "xL_corr_B0_cutFirstBin",
    //     "Truth x_L",
    //     "B0 Reco x_L (first bin cut)",
    //     "figs/response_matrix_xL_B0_cutFirstBin.png",
    //     false,
    //     false,
    //     {0.0, 2.0},
    //     {0.0, 2.0}
    // ));

    // plots.push_back(new PlotOptionsResponseMatrix(
    //     "xpom_corr_B0_cutFirstBin",
    //     "Truth x_{pom}",
    //     "B0 Reco x_{pom} (first bin cut)",
    //     "figs/response_matrix_xpom_B0_cutFirstBin.png",
    //     true,   // logX
    //     true,   // logY
    //     {1e-4, 0.4},
    //     {1e-4, 0.4}
    // ));

//     plot_ptr = new PlotOptions1D(
//     {"MX2_MC","MX2_eX"},//, "xL_eX"},
//     {"Truth", "B0", "Roman Pot", "eX"},
//     {"hist", "hist", "hist", "hist"},
//     "M_{X}^{2}",
//     "M_{X}^{2}",
//     "Counts",
//     "figs/MX2_comparison.png",
//     false,  // isLogX
//     true,   // isLogY
//     false   // normalize
// );
//     plot_ptr->SetLegendPosition(0.75, 0.7, 0.95, 0.9);
//     plots.push_back(plot_ptr);



/////////////////////////////////////////////////////////////////////////////////
    gSystem->mkdir("figs", kTRUE);

    for (const auto& plot : plots) {
        plot->Plot(inputFile);
    }
    
    for (const auto& plot : plots) {
        delete plot;
    }

    // delete plot_ptr;

    inputFile->Close();
    delete inputFile;

    std::cout << "\nPlotting complete! All plots saved to figs/ directory:" << std::endl;
    // std::cout << "- t_distributions.png" << std::endl;
    // std::cout << "- t_pdf_comparison.png" << std::endl;
    // std::cout << "- theta_distributions.png" << std::endl;
    // std::cout << "- t_resolution_B0.png" << std::endl;
    // std::cout << "- t_resolution_RP.png" << std::endl;
    // std::cout << "- t_correlation_combined.png (using TGraphs - smooth!)" << std::endl;
    // std::cout << "- response_matrix_B0.png" << std::endl;
    // std::cout << "- response_matrix_RP.png" << std::endl;

    return 0;
}