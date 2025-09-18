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

    // 1D histogram plot configuration
    plots.push_back(new PlotOptions1D(
        {"h_Q2_truth", "h_Q2_EM", "h_Q2_DA","h_Q2_ESigma"},                // hist names
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. ESigma"}, // legend entries
        {"hist", "pe", "pe","pe"},                                // draw options
        "Q^{2} Reconstruction Methods",                      // canvas title
        "Q^{2}",                                             // x label
        "# of events",                                       // y label
        "figs/Q2_hist.png",                                  // save name
        true,                                                // isLogX (default is false)
        true                                                 // isLogY (default is false)
    ));


    // normalized to PDF plot 
    plots.push_back(new PlotOptions1D( 
        {"h_Q2_truth", "h_Q2_EM", "h_Q2_DA","h_Q2_ESigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. ESigma"},
        {"hist", "pe", "pe","pe"},
        "Q^{2} PDF Comparison",
        "Q^{2}",
        "PDF",
        "figs/Q2_pdf.png",
        true,  // logX
        true, // logY
        true   // normalize to PDF (default is false)
    ));

    
    // ---- x_Bj: 1D distributions (counts) ----
    plots.push_back(new PlotOptions1D(
        {"x_truth", "x_EM", "x_DA", "x_ESigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. ESigma"},
        {"hist", "pe", "pe", "pe"},
        "x_{Bj} Reconstruction Methods",
        "x_{Bj}",
        "# of events",
        "figs/x_hist.png",
        true,   // isLogX
        true    // isLogY
    ));

    // ---- x_Bj: normalized to PDF ----
    plots.push_back(new PlotOptions1D(
        {"x_truth", "x_EM", "x_DA", "x_ESigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. ESigma"},
        {"hist", "pe", "pe", "pe"},
        "x_{Bj} PDF Comparison",
        "x_{Bj}",
        "PDF",
        "figs/x_pdf.png",
        true,  // logX
        true,  // logY
        true   // normalize to PDF
    ));

    // ---- y (inelasticity): 1D distributions (counts) ----
    plots.push_back(new PlotOptions1D(
        {"y_truth", "y_EM", "y_DA", "y_ESigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. ESigma"},
        {"hist", "pe", "pe", "pe"},
        "y (inelasticity) Reconstruction Methods",
        "y",
        "# of events",
        "figs/y_hist.png",
        false,  // isLogX -> linear in [0,1]
        true    // isLogY (optional; keep for dynamic range)
    ));

    // ---- y (inelasticity): normalized to PDF ----
    plots.push_back(new PlotOptions1D(
        {"y_truth", "y_EM", "y_DA", "y_ESigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. ESigma"},
        {"hist", "pe", "pe", "pe"},
        "y (inelasticity) PDF Comparison",
        "y",
        "PDF",
        "figs/y_pdf.png",
        false, // logX
        false, // logY
        true   // normalize to PDF
    ));
// Relative resolution plot configurations
    // plots.push_back(new PlotOptionsRelRes(
    //     "Q2_RelRes_EM",// histogram name
    //     "#frac{Q^{2}_{EM} - Q^{2}_{MC}}{ Q^{2}_{MC}}",// x label
    //     "Counts",// y label
    //     -0.025, 0.02,// fit range
    //     "figs/DDIS_Q2RelRes_EM.png"// save name
    // ));
    // plots.push_back(new PlotOptionsRelRes(
    //     "Q2_RelRes_DA",
    //     "Q^{2}_{DA} - Q^{2}_{MC} / Q^{2}_{MC}",
    //     "Counts",
    //     -0.02, 0.03,
    //     "figs/DDIS_Q2RelRes_DA.png"
    // ));
    // plots.push_back(new PlotOptionsRelRes(
    //     "Q2_RelRes_ESigma",
    //     "Q^{2}_{reco} - Q^{2}_{MC} / Q^{2}_{MC}",
    //     "Counts",
    //     -0.02, 0.03,
    //     "figs/DDIS_Q2RelRes_ESigma.png"
    // ));

    // Binned relative resolution plot configurations
    plots.push_back(new PlotOptionsBinnedRelRes(
        "Q2_RelRes_binned_EM",
        "Relative bin by bin resolution (EM);Q^{2}_{MC};#frac{Q^{2}_{EM} - Q^{2}_{MC}}{Q^{2}_{MC}}",// title, x label, y label
        "Q^{2}_{EM}",
        "",
        {
         {-0.0, 0.0}, {-0.0, 0.0},{-0.0, 0.0},{-0.0, 0.0},{-0.0, 0.0},
         {-0.0, 0.0},{-0.0, 0.0},{-0.02, 0.02},{-0.02, 0.025},{-0.022, 0.025},
         {-0.027, 0.028},{-0.018, 0.02},{-0.022, 0.02},{-0.02, 0.015},{-0.018, 0.02},
         {-0.02, 0.015},{-0.02, 0.017},{-0.017, 0.02},{-0.02, 0.02},{-0.04, 0.04},
         {-0.025, 0.03},{-0.015, 0.025},{-0.05, 0.06}
        },
        "figs/DDIS_Q2RelRes_binned_EM.png",
        "DDIS_Q2RelRes_binned_EM"
    ));
    // plots.push_back(new PlotOptionsBinnedRelRes(
    //     "Q2_RelRes_binned_DA",
    //     "Relative bin by bin resolution (DA);Q^{2}_{MC};#frac{Q^{2}_{DA} - Q^{2}_{MC}}{Q^{2}_{MC}}", // title, x label, y label
    //     "Q^{2}_{DA}",
    //     "",
    //     {{-0.0, 0.0}, {-0.0, 0.0},{-0.0, 0.0},{-0.0, 0.0},{-0.0, 0.0},
    //         {-0.0, 0.0},{-0.0, 0.0},{-0.01, 0.03},{-0.01, 0.025},{-0.01, 0.025},
    //         {-0.015, 0.035},{-0.01, 0.025},{-0.01, 0.025},{-0.019, 0.029},{-0.01, 0.02},
    //         {-0.015, 0.03},{-0.01, 0.0275},{-0.017, 0.027},{-0.025, 0.03},{-0.08, 0.08},
    //         {-0.05, 0.06},{-0.05, 0.065},{-0.05, 0.06}},
    //     "figs/DDIS_Q2RelRes_binned_DA.png",
    //     "DDIS_Q2RelRes_binned_DA"
    // ));
    // plots.push_back(new PlotOptionsBinnedRelRes(
    //     "Q2_RelRes_binned_ESigma",
    //     ";Q^{2}_{MC};#frac{Q^{2}_{ESigma} - Q^{2}_{MC}}{Q^{2}_{MC}}",
    //     "Q^{2}_{ESigma}",
    //     "",
    //     {},
    //     "figs/DDIS_Q2RelRes_binned_ESigma.png",
    //     "DDIS_Q2RelRes_binned_ESigma"
    // ));

    // Create an instance of the new plotting class
    // Parameters: histName, xLabel, yLabel, saveName
    plots.push_back(new PlotOptionsResponseMatrix (
        "Corr_Q2_EM",       // name of 2D histogram
        "Q^{2} (true) [GeV]",   // X-axis label
        "Q^{2} (EM) [GeV]",   // Y-axis label
        "figs/response_matrix_EM.png" // Output file name
    ));

        plots.push_back(new PlotOptionsResponseMatrix (
        "Corr_Q2_DA",       // name of 2D histogram
        "Q^{2} (true) [GeV]",   // X-axis label
        "Q^{2} (DA) [GeV]",   // Y-axis label
        "figs/response_matrix_DA.png" // Output file name
    ));

        plots.push_back(new PlotOptionsResponseMatrix (
        "Corr_Q2_ESigma",       // name of 2D histogram
        "Q^{2} (true) [GeV]",   // X-axis label
        "Q^{2} (ESigma) [GeV]",   // Y-axis label
        "figs/response_matrix_Esigma.png" // Output file name
    ));

    
    // ---- Response Matrices for x_{Bj} (log-log) ----
    plots.push_back(new PlotOptionsResponseMatrix(
        "x_Corr_EM",
        "x_{Bj} (true)",
        "x_{Bj} (EM)",
        "figs/response_matrix_x_EM.png",
        true,  // isLogX
        true   // isLogY
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "x_Corr_DA",
        "x_{Bj} (true)",
        "x_{Bj} (DA)",
        "figs/response_matrix_x_DA.png",
        true,  // isLogX
        true   // isLogY
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "x_Corr_ESigma",
        "x_{Bj} (true)",
        "x_{Bj} (ESigma)",
        "figs/response_matrix_x_ESigma.png",
        true,  // isLogX
        true   // isLogY
    ));

    // ---- Response Matrices for y (linear-linear) ----
    plots.push_back(new PlotOptionsResponseMatrix(
        "y_Corr_EM",
        "y (true)",
        "y (EM)",
        "figs/response_matrix_y_EM.png",
        false, // isLogX
        false,  // isLogY
        {0.,1.}, // x range {0.,1.}
        {0.,1.} // y range {0.,1.}

    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "y_Corr_DA",
        "y (true)",
        "y (DA)",
        "figs/response_matrix_y_DA.png",
        false, // isLogX
        false,  // isLogY
        {0.,1.}, // x range {0.,1.}
        {0.,1.} // y range {0.,1.}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "y_Corr_ESigma",
        "y (true)",
        "y (ESigma)",
        "figs/response_matrix_y_ESigma.png",
        false, // isLogX
        false,  // isLogY
        {0.,1.}, // x range {0.,1.}
        {0.,1.} // y range {0.,1.}
    ));




    // === Relative resolution (Δ/true) 1D for Q2, x_Bj, y ===
    plots.push_back(new PlotOptionsRelRes(
        "Q2_RelRes_EM", "#Delta Q^{2}/Q^{2}_{true}", "Counts", -0.2, 0.2, "figs/Q2RelRes_EM.png"));
    plots.push_back(new PlotOptionsRelRes(
        "Q2_RelRes_DA", "#Delta Q^{2}/Q^{2}_{true}", "Counts", -0.2, 0.2, "figs/Q2RelRes_DA.png"));
    plots.push_back(new PlotOptionsRelRes(
        "Q2_RelRes_ESigma", "#Delta Q^{2}/Q^{2}_{true}", "Counts", -0.2, 0.2, "figs/Q2RelRes_ESigma.png"));

    plots.push_back(new PlotOptionsRelRes(
        "x_RelRes_EM", "#Delta x / x_{true}", "Counts", -0.5, 0.5, "figs/xRelRes_EM.png"));
    plots.push_back(new PlotOptionsRelRes(
        "x_RelRes_DA", "#Delta x / x_{true}", "Counts", -0.5, 0.5, "figs/xRelRes_DA.png"));
    plots.push_back(new PlotOptionsRelRes(
        "x_RelRes_ESigma", "#Delta x / x_{true}", "Counts", -0.5, 0.5, "figs/xRelRes_ESigma.png"));

    plots.push_back(new PlotOptionsRelRes(
        "y_RelRes_EM", "#Delta y / y_{true}", "Counts", -0.5, 0.5, "figs/yRelRes_EM.png"));
    plots.push_back(new PlotOptionsRelRes(
        "y_RelRes_DA", "#Delta y / y_{true}", "Counts", -0.5, 0.5, "figs/yRelRes_DA.png"));
    plots.push_back(new PlotOptionsRelRes(
        "y_RelRes_ESigma", "#Delta y / y_{true}", "Counts", -0.5, 0.5, "figs/yRelRes_ESigma.png"));

    // === Binned relative resolution vs truth (auto-fit ranges) ===
    // Q2: log-x
    plots.push_back(new PlotOptionsBinnedRelRes(
        "Q2_RelRes_binned_EM",
        "Relative bin-by-bin resolution (EM);Q^{2}_{true};#Delta Q^{2}/Q^{2}_{true}",
        "Q^{2}_{EM}", "", { /* auto */ }, "figs/Q2RelRes_binned_EM.png", "Q2RelRes_binned_EM"));
    plots.push_back(new PlotOptionsBinnedRelRes(
        "Q2_RelRes_binned_DA",
        "Relative bin-by-bin resolution (DA);Q^{2}_{true};#Delta Q^{2}/Q^{2}_{true}",
        "Q^{2}_{DA}", "", { /* auto */ }, "figs/Q2RelRes_binned_DA.png", "Q2RelRes_binned_DA"));
    plots.push_back(new PlotOptionsBinnedRelRes(
        "Q2_RelRes_binned_ESigma",
        "Relative bin-by-bin resolution (ESigma);Q^{2}_{true};#Delta Q^{2}/Q^{2}_{true}",
        "Q^{2}_{ESigma}", "", { /* auto */ }, "figs/Q2RelRes_binned_ESigma.png", "Q2RelRes_binned_ESigma"));

    // x_Bj: log-x
    plots.push_back(new PlotOptionsBinnedRelRes(
        "x_RelRes_binned_EM",
        "Relative bin-by-bin resolution (EM);x_{Bj}^{true};#Delta x / x_{true}",
        "x_{Bj}^{EM}", "", { /* auto */ }, "figs/xRelRes_binned_EM.png", "xRelRes_binned_EM"));
    plots.push_back(new PlotOptionsBinnedRelRes(
        "x_RelRes_binned_DA",
        "Relative bin-by-bin resolution (DA);x_{Bj}^{true};#Delta x / x_{true}",
        "x_{Bj}^{DA}", "", { /* auto */ }, "figs/xRelRes_binned_DA.png", "xRelRes_binned_DA"));
    plots.push_back(new PlotOptionsBinnedRelRes(
        "x_RelRes_binned_ESigma",
        "Relative bin-by-bin resolution (ESigma);x_{Bj}^{true};#Delta x / x_{true}",
        "x_{Bj}^{ESigma}", "", { /* auto */ }, "figs/xRelRes_binned_ESigma.png", "xRelRes_binned_ESigma"));

    // y: linear x
    plots.push_back(new PlotOptionsBinnedRelRes(
        "y_RelRes_binned_EM",
        "Relative bin-by-bin resolution (EM);y^{true};#Delta y / y_{true}",
        "y^{EM}", "", { /* auto */ }, "figs/yRelRes_binned_EM.png", "yRelRes_binned_EM"));
    plots.push_back(new PlotOptionsBinnedRelRes(
        "y_RelRes_binned_DA",
        "Relative bin-by-bin resolution (DA);y^{true};#Delta y / y_{true}",
        "y^{DA}", "", { /* auto */ }, "figs/yRelRes_binned_DA.png", "yRelRes_binned_DA"));
    plots.push_back(new PlotOptionsBinnedRelRes(
        "y_RelRes_binned_ESigma",
        "Relative bin-by-bin resolution (ESigma);y^{true};#Delta y / y_{true}",
        "y^{ESigma}", "", { /* auto */ }, "figs/yRelRes_binned_ESigma.png", "yRelRes_binned_ESigma"));
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