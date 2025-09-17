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
        {"MC: truth", "Reco.: electron method", "Reco. DA", "Reco. ESigma"}, // legend entries
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
        {"MC: truth", "Reco.: electron method", "Reco. DA", "Reco. ESigma"},
        {"hist", "pe", "pe","pe"},
        "Q^{2} PDF Comparison",
        "Q^{2}",
        "PDF",
        "figs/Q2_pdf.png",
        true,  // logX
        true, // logY
        true   // normalize to PDF (default is false)
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