//g++ -o DDIS_Plots_Q2 $(root-config --cflags --glibs) DDIS_Plots_Q2.cpp Plotting.cpp

#include "Plotting.hpp"
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

    // Q2 1D histogram plot configuration
    // plots.push_back(new PlotOptions1D(
    //     {"h_Q2_truth", "h_Q2_EM", "h_Q2_DA","h_Q2_ESigma"},    // hist names
    //     {"MC: truth", "Reco. EM", "Reco. DA", "Reco. ESigma"}, // legend entries
    //     {"hist", "pe", "pe","pe"},                             // draw options
    //     "Q^{2} Reconstruction Methods",                        // canvas title
    //     "Q^{2}",                                               // x label
    //     "# of events",                                         // y label
    //     "figs/Q2_hist.png",                                    // save name
    //     true,                                                  // isLogX (default is false)
    //     true                                                   // isLogY (default is false)
    // ));


    // // Q2 normalized to PDF plot 
    // plots.push_back(new PlotOptions1D( 
    //     {"h_Q2_truth", "h_Q2_EM", "h_Q2_DA","h_Q2_ESigma"},
    //     {"MC: truth", "Reco. EM", "Reco. DA", "Reco. ESigma"},
    //     {"hist", "pe", "pe","pe"},
    //     "Q^{2} PDF Comparison",
    //     "Q^{2}",
    //     "PDF",
    //     "figs/Q2_pdf.png",
    //     true,  // logX
    //     true,  // logY
    //     true   // normalize to PDF (default is false)
    // ));

    
    // // ---- x_Bj: 1D distributions (counts) ----
    // plots.push_back(new PlotOptions1D(
    //     {"x_truth", "x_EM", "x_DA", "x_ESigma"},
    //     {"MC: truth", "Reco. EM", "Reco. DA", "Reco. ESigma"},
    //     {"hist", "pe", "pe", "pe"},
    //     "x_{Bj} Reconstruction Methods",
    //     "x_{Bj}",
    //     "# of events",
    //     "figs/x_hist.png",
    //     true,   // isLogX
    //     true    // isLogY
    // ));

    // // ---- x_Bj: normalized to PDF ----
    // plots.push_back(new PlotOptions1D(
    //     {"x_truth", "x_EM", "x_DA", "x_ESigma"},
    //     {"MC: truth", "Reco. EM", "Reco. DA", "Reco. ESigma"},
    //     {"hist", "pe", "pe", "pe"},
    //     "x_{Bj} PDF Comparison",
    //     "x_{Bj}",
    //     "PDF",
    //     "figs/x_pdf.png",
    //     true,  // logX
    //     true,  // logY
    //     true   // normalize to PDF
    // ));

    // // ---- y (inelasticity): 1D distributions (counts) ----
    // plots.push_back(new PlotOptions1D(
    //     {"y_truth", "y_EM", "y_DA", "y_ESigma"},
    //     {"MC: truth", "Reco. EM", "Reco. DA", "Reco. ESigma"},
    //     {"hist", "pe", "pe", "pe"},
    //     "y (inelasticity) Reconstruction Methods",
    //     "y",
    //     "# of events",
    //     "figs/y_hist.png",
    //     false,  // isLogX -> linear in [0,1]
    //     true    // isLogY (optional; keep for dynamic range)
    // ));

    // // ---- y (inelasticity): normalized to PDF ----
    // plots.push_back(new PlotOptions1D(
    //     {"y_truth", "y_EM", "y_DA", "y_ESigma"},
    //     {"MC: truth", "Reco. EM", "Reco. DA", "Reco. ESigma"},
    //     {"hist", "pe", "pe", "pe"},
    //     "y (inelasticity) PDF Comparison",
    //     "y",
    //     "PDF",
    //     "figs/y_pdf.png",
    //     false, // logX
    //     false, // logY
    //     true   // normalize to PDF
    // ));

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Relative resolution plot configurations
    // plots.push_back(new PlotOptionsRelRes(
    //     "Q2_RelRes_EM",// histogram name
    //     "#frac{Q^{2}_{EM} - Q^{2}_{MC}}{ Q^{2}_{MC}}",// x label
    //     "Counts",// y label
    //     -0.01, 0.01,// fit range
    //     "figs/DDIS_Q2RelRes_EM.png"// save name
    // ));
    // plots.push_back(new PlotOptionsRelRes(
    //     "Q2_RelRes_DA",
    //     "Q^{2}_{DA} - Q^{2}_{MC} / Q^{2}_{MC}",
    //     "Counts",
    //     -0.005, 0.02,
    //     "figs/DDIS_Q2RelRes_DA.png"
    // ));
    // plots.push_back(new PlotOptionsRelRes(
    //     "Q2_RelRes_ESigma",
    //     "Q^{2}_{reco} - Q^{2}_{MC} / Q^{2}_{MC}",
    //     "Counts",
    //     -0.01, 0.01,// fit range
    //     "figs/DDIS_Q2RelRes_ESigma.png"
    // ));

    // // Relative resolution plot configurations for x_Bj
    // plots.push_back(new PlotOptionsRelRes(
    //     "x_RelRes_EM",// histogram name
    //     "#frac{x_{EM} - x_{MC}}{ x_{MC}}",// x label
    //     "Counts",// y label
    //     -0.025, 0.02,// fit range
    //     "figs/DDIS_RelRes_xBj_EM.png"// save name
    // ));
    // plots.push_back(new PlotOptionsRelRes(
    //     "x_RelRes_DA",
    //     "#frac{x_{DA} - x_{MC}}{X_{MC}}",
    //     "Counts",
    //     0., 0.,// Skip fitting and just save the histogram
    //     "figs/DDIS_RelRes_xBj_DA.png"
    // ));
    // plots.push_back(new PlotOptionsRelRes(
    //     "x_RelRes_ESigma",
    //     "#frac{x_{E#Sigma} - x_{MC}}{X_{MC}}",
    //     "Counts",
    //     0., 0.,// Skip fitting and just save the histogram
    //     "figs/DDIS_RelRes_x_ESigma.png"
    // ));

    // // Relative resolution plot configurations for 'y'
    // plots.push_back(new PlotOptionsRelRes(
    //     "y_RelRes_EM",// histogram name
    //     "#frac{y_{EM} - y_{MC}}{ y_{MC}}",// x label
    //     "Counts",// y label
    //     -0.009, 0.009,// fit range
    //     // -999, -999.,// auto fit range
    //     "figs/DDIS_RelRes_y_EM.png"// save name
    // ));
    // plots.push_back(new PlotOptionsRelRes(
    //     "y_RelRes_DA",
    //     "#frac{y_{DA} - y_{MC}}{y_{MC}}",
    //     "Counts",
    //     0., 0.,// Skip fitting and just save the histogram
    //     "figs/DDIS_RelRes_y_DA.png"
    // ));
    // plots.push_back(new PlotOptionsRelRes(
    //     "y_RelRes_ESigma",
    //     "#frac{y_{E#Sigma} - y_{MC}}{y_{MC}}",
    //     "Counts",
    //     0., 0.,// Skip fitting and just save the histogram
    //     "figs/DDIS_RelRes_y_ESigma.png"
    // ));

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

    // Binned relative resolution plot configurations
    // plots.push_back(new PlotOptionsBinnedRelRes(
    //     "Q2_RelRes_binned_EM",
    //     "Relative bin by bin resolution (EM);Q^{2}_{MC};#frac{Q^{2}_{EM} - Q^{2}_{MC}}{Q^{2}_{MC}}",// title, x label, y label
    //     "Q^{2}_{EM}",
    //     "",
    //     {
    //      {-0.0, 0.0},/*2*/ {-0.022, 0.02},{-0.02, 0.02},{-0.02, 0.02},{-0.02, 0.02},
    //      {-0.015, 0.015},/*7*/{-0.015, 0.015},/*8*/{-0.014, 0.015},{-0.025, 0.025},{-0.01, 0.012},
    //      {-0.027, 0.028},{-0.018, 0.02},{-0.022, 0.02},{-0.02, 0.015},{-0.018, 0.02},
    //      {-0.02, 0.015},{-0.02, 0.017},{-0.017, 0.02},{-0.02, 0.02},{-0.04, 0.04},
    //      {-0.025, 0.03},{-0.015, 0.025},{-0.05, 0.06}
    //     },
    //     "figs/DDIS_Q2RelRes_binned_EM.png",
    //     "DDIS_Q2RelRes_binned_EM"
    // ));
    // plots.push_back(new PlotOptionsBinnedRelRes(
    //     "Q2_RelRes_binned_DA",
    //     "Relative bin by bin resolution (DA);Q^{2}_{MC};#frac{Q^{2}_{DA} - Q^{2}_{MC}}{Q^{2}_{MC}}", // title, x label, y label
    //     "Q^{2}_{DA}",
    //     "",
    //     {
    //       {-0.0, 0.0}, {-0.005, 0.025},{-0.01, 0.025},/*4*/{-0.006, 0.02},{-0.01, 0.02},
    //       {-0.009, 0.02},/*7*/{-0.008, 0.02},{-0.01, 0.025},{-0.04, 0.04},{-0.03, 0.05},
    //       {-0.015, 0.03},{-0.009, 0.02},{-0.01, 0.02},/*14*/{-0.01, 0.02},{-0.01, 0.02},
    //       {-0.01, 0.02},/*17*/{-0.004, 0.02},{-0.017, 0.027},{-0.025, 0.03},{-0.08, 0.08},
    //       {-0.05, 0.06},{-0.05, 0.065},{-0.05, 0.06}
    //     },
    //     "figs/DDIS_Q2RelRes_binned_DA.png",
    //     "DDIS_Q2RelRes_binned_DA"
    // ));
    // plots.push_back(new PlotOptionsBinnedRelRes(
    //     "Q2_RelRes_binned_ESigma",
    //     ";Q^{2}_{MC};#frac{Q^{2}_{E#Sigma} - Q^{2}_{MC}}{Q^{2}_{MC}}",
    //     "Q^{2}_{E#Sigma}",
    //     "",
    //     {
    //      {-0.0, 0.0},/*2*/ {-0.022, 0.02},{-0.02, 0.02},{-0.02, 0.02},{-0.02, 0.02},
    //      {-0.015, 0.015},/*7*/{-0.015, 0.015},/*8*/{-0.014, 0.015},{-0.025, 0.025},{-0.01, 0.012},
    //      {-0.027, 0.028},{-0.018, 0.02},{-0.022, 0.02},{-0.02, 0.015},{-0.018, 0.02},
    //      {-0.02, 0.015},{-0.02, 0.017},{-0.017, 0.02},{-0.02, 0.02},{-0.04, 0.04},
    //      {-0.025, 0.03},{-0.015, 0.025},{-0.05, 0.06}
    //     },
    //     "figs/DDIS_Q2RelRes_binned_ESigma.png",
    //     "DDIS_Q2RelRes_binned_ESigma"
    // ));

    plots.push_back(new PlotOptionsBinnedRelRes(
        "x_RelRes_binned_EM",
        ";x_{MC};#frac{x_{EM} - x_{MC}}{x_{MC}}",
        "x_{EM}",
        "",
        {
         {-0.0, 0.0},/*2*/ {-0.022, 0.02},{-0.02, 0.02},{-0.02, 0.02},{-0.02, 0.02},
         {-0.015, 0.015},/*7*/{-0.015, 0.015},/*8*/{-0.014, 0.015},{-0.025, 0.025},{-0.01, 0.012},
         {-0.027, 0.028},{-0.018, 0.02},{-0.022, 0.02},{-0.02, 0.015},{-0.018, 0.02},
         {-0.02, 0.015},{-0.02, 0.017},{-0.017, 0.02},{-0.02, 0.02},{-0.04, 0.04},
         {-0.025, 0.03},{-0.015, 0.025},{-0.05, 0.06}
        },
        "figs/DDIS_RelRes_binned_x_EM.png",
        "DDIS_RelRes_binned_x_EM"
    ));
    // plots.push_back(new PlotOptionsBinnedRelRes(
    //     "x_RelRes_binned_DA",
    //     ";x_{MC};#frac{x_{DA} - x_{MC}}{x_{MC}}",
    //     "x_{DA}",
    //     "",
    //     {},
    //     "figs/DDIS_RelRes_binned_x_DA.png",
    //     "DDIS_RelRes_binned_x_DA"
    // ));
    // plots.push_back(new PlotOptionsBinnedRelRes(
    //     "x_RelRes_binned_ESigma",
    //     ";x_{MC};#frac{x_{E#Sigma} - x_{MC}}{x_{MC}}",
    //     "x_{ESigma}",
    //     "",
    //     {},
    //     "figs/DDIS_RelRes_binned_x_ESigma.png",
    //     "DDIS_RelRes_binned_x_ESigma"
    // ));

    // plots.push_back(new PlotOptionsBinnedRelRes(
    //     "y_RelRes_binned_EM",
    //     ";y_{MC};#frac{y_{EM} - y_{MC}}{y_{MC}}",
    //     "y_{EM}",
    //     "",
    //     {},
    //     "figs/DDIS_RelRes_binned_y_EM.png",
    //     "DDIS_RelRes_binned_y_EM"
    // ));
    // plots.push_back(new PlotOptionsBinnedRelRes(
    //     "y_RelRes_binned_DA",
    //     ";y_{MC};#frac{y_{DA} - y_{MC}}{y_{MC}}",
    //     "y_{DA}",
    //     "",
    //     {},
    //     "figs/DDIS_RelRes_binned_y_DA.png",
    //     "DDIS_RelRes_binned_y_DA"
    // ));
    // plots.push_back(new PlotOptionsBinnedRelRes(
    //     "y_RelRes_binned_ESigma",
    //     ";y_{MC};#frac{y_{E#Sigma} - y_{MC}}{y_{MC}}",
    //     "y_{ESigma}",
    //     "",
    //     {},
    //     "figs/DDIS_RelRes_binned_y_ESigma.png",
    //     "DDIS_RelRes_binned_y_ESigma"
    // ));
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///////// Response Matrix plot configurations

    // Create an instance of the new plotting class
    // Parameters: histName, xLabel, yLabel, saveName


    // plots.push_back(new PlotOptionsResponseMatrix (
    //     "Corr_Q2_EM",       // name of 2D histogram
    //     "Q^{2} (true) [GeV]",   // X-axis label
    //     "Q^{2} (EM) [GeV]",   // Y-axis label
    //     "figs/response_matrix_EM.png" // Output file name
    // ));

    // plots.push_back(new PlotOptionsResponseMatrix (
    //     "Corr_Q2_DA",       // name of 2D histogram
    //     "Q^{2} (true) [GeV]",   // X-axis label
    //     "Q^{2} (DA) [GeV]",   // Y-axis label
    //     "figs/response_matrix_DA.png" // Output file name
    // ));

    // plots.push_back(new PlotOptionsResponseMatrix (
    //     "Corr_Q2_ESigma",       // name of 2D histogram
    //     "Q^{2} (true) [GeV]",   // X-axis label
    //     "Q^{2} (ESigma) [GeV]",   // Y-axis label
    //     "figs/response_matrix_Esigma.png" // Output file name
    // ));

    
    // // ---- Response Matrices for x_{Bj} (log-log) ----
    // plots.push_back(new PlotOptionsResponseMatrix(
    //     "x_Corr_EM",
    //     "x_{Bj} (true)",
    //     "x_{Bj} (EM)",
    //     "figs/response_matrix_x_EM.png",
    //     true,  // isLogX
    //     true,   // isLogY
    //     {1e-4,0.5},
    //     {1e-4,0.5}
    // ));
    // plots.push_back(new PlotOptionsResponseMatrix(
    //     "x_Corr_DA",
    //     "x_{Bj} (true)",
    //     "x_{Bj} (DA)",
    //     "figs/response_matrix_x_DA.png",
    //     true,  // isLogX
    //     true,   // isLogY
    //     {1e-4,0.5},
    //     {1e-4,0.5}
    // ));
    // plots.push_back(new PlotOptionsResponseMatrix(
    //     "x_Corr_ESigma",
    //     "x_{Bj} (true)",
    //     "x_{Bj} (ESigma)",
    //     "figs/response_matrix_x_ESigma.png",
    //     true,  // isLogX
    //     true,   // isLogY
    //     {1e-4,0.5},
    //     {1e-4,0.5}
    // ));

    // // // ---- Response Matrices for y (linear-linear) ----
    // plots.push_back(new PlotOptionsResponseMatrix(
    //     "y_Corr_EM",
    //     "y (true)",
    //     "y (EM)",
    //     "figs/response_matrix_y_EM.png",
    //     false, // isLogX
    //     false,  // isLogY
    //     {0.,1.}, // x range {0.,1.}
    //     {0.,1.}  // y range {0.,1.}
    // ));

    // plots.push_back(new PlotOptionsResponseMatrix(
    //     "y_Corr_DA",
    //     "y (true)",
    //     "y (DA)",
    //     "figs/response_matrix_y_DA.png",
    //     false, // isLogX
    //     false,  // isLogY
    //     {0.,1.}, // x range {0.,1.}
    //     {0.,1.} // y range {0.,1.}
    // ));

    // plots.push_back(new PlotOptionsResponseMatrix(
    //     "y_Corr_ESigma",
    //     "y (true)",
    //     "y (ESigma)",
    //     "figs/response_matrix_y_ESigma.png",
    //     false, // isLogX
    //     false,  // isLogY
    //     {0.,1.}, // x range {0.,1.}
    //     {0.,1.} // y range {0.,1.}
    // ));



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/////////////////// plot all the configurations //////////////////////
// Loop through the vector and plot each configuration
    for (const auto& plot : plots) {
        plot->Plot(inputFile);
    }
    
    // Clean up
    for (const auto& plot : plots) {
        delete plot;
    }

    inputFile->Close();
    delete inputFile;

    return 0;
}