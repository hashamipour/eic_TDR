// g++ -o DDIS_Plot_AllKinematics $(root-config --cflags --glibs) DDIS_Plot_AllKinematics.cpp Plotting.cpp Utility.cpp
// Usage: ./DDIS_Plot_AllKinematics DDIS_AllKinematics.root

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
    const char* inputFileName = argv[1];

    // ROOT style similar to your existing plotters
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(120);
    gStyle->SetPalette(kBlueRedYellow);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
    gStyle->SetOptFit(0);

    TFile* inputFile = TFile::Open(inputFileName);
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << inputFileName << std::endl;
        return 1;
    }

    std::vector<PlotOptions*> plots;

    // -----------------------------
    // Proton-side variables (MC, B0, RP)
    // -----------------------------
    {
        // |t| overlays
        PlotOptions1D* p = new PlotOptions1D(
            {"proton/t_MC", "proton/t_B0", "proton/t_RP"},
            {"MC Truth", "Reco B0 (5.5-20 mrad)", "Reco RP (<5 mrad)"},
            {"hist", "pe", "pe"},
            "Mandelstam |t| Distributions",
            "|t| [GeV^{2}]",
            "Counts",
            "figs/t_distributions_logy.png",
            false, // logX
            true   // logY
        );
        p->SetLegendPosition(0.60, 0.70, 0.82, 0.90);
        plots.push_back(p);

        // |t| PDF (normalized)
        p = new PlotOptions1D(
            {"proton/t_MC", "proton/t_B0", "proton/t_RP"},
            {"MC Truth", "Reco B0 (5.5-20 mrad)", "Reco RP (<5 mrad)"},
            {"hist", "pe", "pe"},
            "Mandelstam |t| — PDF Comparison",
            "|t| [GeV^{2}]",
            "PDF",
            "figs/t_pdf_comparison.png",
            false, // logX
            true,  // logY
            true   // normalize to PDF
        );
        p->SetLegendPosition(0.60, 0.70, 0.82, 0.90);
        plots.push_back(p);

        // θ overlays (mrad)
        p = new PlotOptions1D(
            {"proton/theta_MC", "proton/theta_B0", "proton/theta_RP"},
            {"MC Truth", "Reco B0 (5.5-20 mrad)", "Reco RP (<5 mrad)"},
            {"hist", "pe", "pe"},
            "Scattered Proton Polar Angle",
            "#theta [mrad]",
            "Counts",
            "figs/theta_distributions.png",
            false, // logX
            false  // logY
        );
        p->SetLegendPosition(0.60, 0.70, 0.82, 0.90);
        plots.push_back(p);

        // xL overlays
        p = new PlotOptions1D(
            {"proton/xL_MC", "proton/xL_B0", "proton/xL_RP"},
            {"MC Truth", "Reco B0 (5.5-20 mrad)", "Reco RP (<5 mrad)"},
            {"hist", "pe", "pe"},
            "x_{L} Distributions",
            "x_{L}",
            "Counts",
            "figs/xL_distributions.png",
            false, // logX
            false  // logY
        );
        p->SetLegendPosition(0.15, 0.70, 0.35, 0.90);
        plots.push_back(p);

        // MX2 overlays
        p = new PlotOptions1D(
            {"proton/MX2_MC", "proton/MX2_B0", "proton/MX2_RP"},
            {"MC Truth", "Reco B0 (5.5-20 mrad)", "Reco RP (<5 mrad)"},
            {"hist", "pe", "pe"},
            "M_{X}^{2} Distributions",
            "M_{X}^{2} [GeV^{2}]",
            "Counts",
            "figs/MX2_distributions.png",
            false, // logX
            true   // logY
        );
        p->SetLegendPosition(0.60, 0.70, 0.82, 0.90);
        plots.push_back(p);

        // xP overlays
        p = new PlotOptions1D(
            {"proton/xP_MC", "proton/xP_B0", "proton/xP_RP"},
            {"MC Truth", "Reco B0 (5.5-20 mrad)", "Reco RP (<5 mrad)"},
            {"hist", "pe", "pe"},
            "x_{P} Distributions",
            "x_{P}",
            "Counts",
            "figs/xP_distributions.png",
            false, // logX
            true   // logY
        );
        p->SetLegendPosition(0.60, 0.70, 0.82, 0.90);
        plots.push_back(p);

        // Response matrices / correlations (truth vs reco)
        plots.push_back(new PlotOptionsResponseMatrix(
            "proton/t_corr_B0",
            "Truth |t| [GeV^{2}]",
            "Reco |t| [GeV^{2}]",
            "figs/response_matrix_t_B0.png",
            false, false,
            {0.0, 2.0}, {0.0, 2.0}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "proton/t_corr_RP",
            "Truth |t| [GeV^{2}]",
            "Reco |t| [GeV^{2}]",
            "figs/response_matrix_t_RP.png",
            false, false,
            {0.0, 0.5}, {0.0, 0.5}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "proton/xL_corr_B0",
            "Truth x_{L}",
            "Reco x_{L}",
            "figs/response_matrix_xL_B0.png",
            false, false,
            {0.90, 1.02}, {0.90, 1.02}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "proton/xL_corr_RP",
            "Truth x_{L}",
            "Reco x_{L}",
            "figs/response_matrix_xL_RP.png",
            false, false,
            {0.90, 1.02}, {0.90, 1.02}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "proton/MX2_corr_B0",
            "Truth M_{X}^{2} [GeV^{2}]",
            "Reco M_{X}^{2} [GeV^{2}]",
            "figs/response_matrix_MX2_B0.png",
            false, true,
            {0.0, 100.0}, {0.0, 100.0}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "proton/MX2_corr_RP",
            "Truth M_{X}^{2} [GeV^{2}]",
            "Reco M_{X}^{2} [GeV^{2}]",
            "figs/response_matrix_MX2_RP.png",
            false, true,
            {0.0, 100.0}, {0.0, 100.0}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "proton/xP_corr_B0",
            "Truth x_{P}",
            "Reco x_{P}",
            "figs/response_matrix_xP_B0.png",
            false, true,
            {0.0, 0.1}, {0.0, 0.1}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "proton/xP_corr_RP",
            "Truth x_{P}",
            "Reco x_{P}",
            "figs/response_matrix_xP_RP.png",
            false, true,
            {0.0, 0.1}, {0.0, 0.1}
        ));

        // t resolution (1D) — separate for B0 and RP
        plots.push_back(new PlotOptionsRelRes(
            "proton/t_res_B0",
            "(|t|_{reco} - |t|_{truth})/|t|_{truth}",
            "Counts",
            -0.5, 0.5,
            "figs/t_resolution_B0.png"
        ));
        plots.push_back(new PlotOptionsRelRes(
            "proton/t_res_RP",
            "(|t|_{reco} - |t|_{truth})/|t|_{truth}",
            "Counts",
            -0.5, 0.5,
            "figs/t_resolution_RP.png"
        ));
    }

    // -----------------------------
    // Global checks: Sum(E-pz) and Mx (truth vs reco)
    // -----------------------------
    {
        PlotOptions1D* p = new PlotOptions1D(
            {"proton/h_SumEPz_truth", "proton/h_SumEPz_reco"},
            {"Truth", "Reco"},
            {"hist", "pe"},
            "#Sigma(E - p_{z})",
            "#Sigma(E - p_{z}) [GeV]",
            "Counts",
            "figs/SumEPz_comparison.png",
            false, // logX
            false  // logY
        );
        p->SetLegendPosition(0.65, 0.75, 0.85, 0.90);
        plots.push_back(p);

        p = new PlotOptions1D(
            {"proton/h_SumEPz_truth", "proton/h_SumEPz_reco"},
            {"Truth", "Reco"},
            {"hist", "pe"},
            "#Sigma(E - p_{z})",
            "#Sigma(E - p_{z}) [GeV]",
            "Counts",
            "figs/SumEPz_comparison_logy.png",
            false, // logX
            true   // logY
        );
        p->SetLegendPosition(0.65, 0.75, 0.85, 0.90);
        plots.push_back(p);

        p = new PlotOptions1D(
            {"proton/h_Mx_truth", "proton/h_Mx_reco"},
            {"Truth (4-vector)", "Reco (4-vector)"},
            {"hist", "pe"},
            "M_{X} from 4-vector sum",
            "M_{X} [GeV]",
            "Counts",
            "figs/Mx_comparison.png",
            false, // logX
            false  // logY
        );
        p->SetLegendPosition(0.65, 0.75, 0.85, 0.90);
        plots.push_back(p);

        p = new PlotOptions1D(
            {"proton/h_Mx_truth", "proton/h_Mx_reco"},
            {"Truth (4-vector)", "Reco (4-vector)"},
            {"hist", "pe"},
            "M_{X} from 4-vector sum",
            "M_{X} [GeV]",
            "Counts",
            "figs/Mx_comparison_logy.png",
            false, // logX
            true   // logY
        );
        p->SetLegendPosition(0.65, 0.75, 0.85, 0.90);
        plots.push_back(p);
    }

    // -----------------------------
    // Inclusive DIS: Q2, x, y (Truth vs EM/DA/ESigma/Sigma)
    // -----------------------------
    {
        PlotOptions1D* p = nullptr;
        // Q2 — counts
        p = new PlotOptions1D(
            {"inclusive/h_Q2_truth", "inclusive/h_Q2_EM", "inclusive/h_Q2_DA", "inclusive/h_Q2_ESigma", "inclusive/h_Q2_Sigma"},
            {"MC: truth", "Reco: EM", "Reco: DA", "Reco: ESigma", "Reco: Sigma"},
            {"hist", "pe", "pe", "pe", "pe"},
            "Q^{2} Distributions",
            "Q^{2} [GeV^{2}]",
            "Counts",
            "figs/Q2_counts.png",
            true,  // logX
            true   // logY
        );
        p->SetLegendPosition(0.70, 0.70, 0.90, 0.90);
        plots.push_back(p);

        // Q2 — PDFs (normalized)
        p = new PlotOptions1D(
            {"inclusive/h_Q2_truth", "inclusive/h_Q2_EM", "inclusive/h_Q2_DA", "inclusive/h_Q2_ESigma", "inclusive/h_Q2_Sigma"},
            {"MC: truth", "Reco: EM", "Reco: DA", "Reco: ESigma", "Reco: Sigma"},
            {"hist", "pe", "pe", "pe", "pe"},
            "Q^{2} PDF Comparison",
            "Q^{2} [GeV^{2}]",
            "PDF",
            "figs/Q2_pdf.png",
            true,  // logX
            true,  // logY
            true   // normalize to PDF
        );
        p->SetLegendPosition(0.70, 0.70, 0.90, 0.90);
        plots.push_back(p);

        // x — counts
        p = new PlotOptions1D(
            {"inclusive/h_x_truth", "inclusive/h_x_EM", "inclusive/h_x_DA", "inclusive/h_x_ESigma", "inclusive/h_x_Sigma"},
            {"MC: truth", "Reco: EM", "Reco: DA", "Reco: ESigma", "Reco: Sigma"},
            {"hist", "pe", "pe", "pe", "pe"},
            "x Distributions",
            "x",
            "Counts",
            "figs/x_counts.png",
            true,  // logX
            true   // logY
        );
        p->SetLegendPosition(0.70, 0.70, 0.90, 0.90);
        plots.push_back(p);

        // x — PDFs (normalized)
        p = new PlotOptions1D(
            {"inclusive/h_x_truth", "inclusive/h_x_EM", "inclusive/h_x_DA", "inclusive/h_x_ESigma", "inclusive/h_x_Sigma"},
            {"MC: truth", "Reco: EM", "Reco: DA", "Reco: ESigma", "Reco: Sigma"},
            {"hist", "pe", "pe", "pe", "pe"},
            "x PDF Comparison",
            "x",
            "PDF",
            "figs/x_pdf.png",
            true,  // logX
            true,  // logY
            true   // normalize to PDF
        );
        p->SetLegendPosition(0.70, 0.70, 0.90, 0.90);
        plots.push_back(p);

        // y — counts
        p = new PlotOptions1D(
            {"inclusive/h_y_truth", "inclusive/h_y_EM", "inclusive/h_y_DA", "inclusive/h_y_ESigma", "inclusive/h_y_Sigma"},
            {"MC: truth", "Reco: EM", "Reco: DA", "Reco: ESigma", "Reco: Sigma"},
            {"hist", "pe", "pe", "pe", "pe"},
            "y Distributions",
            "y",
            "Counts",
            "figs/y_counts.png",
            false, // logX
            true   // logY
        );
        p->SetLegendPosition(0.70, 0.70, 0.90, 0.90);
        plots.push_back(p);

        // y — PDFs (normalized)
        p = new PlotOptions1D(
            {"inclusive/h_y_truth", "inclusive/h_y_EM", "inclusive/h_y_DA", "inclusive/h_y_ESigma", "inclusive/h_y_Sigma"},
            {"MC: truth", "Reco: EM", "Reco: DA", "Reco: ESigma", "Reco: Sigma"},
            {"hist", "pe", "pe", "pe", "pe"},
            "y PDF Comparison",
            "y",
            "PDF",
            "figs/y_pdf.png",
            false, // logX
            false, // logY
            true   // normalize to PDF
        );
        p->SetLegendPosition(0.70, 0.70, 0.90, 0.90);
        plots.push_back(p);

        // Relative Q2 resolutions (each separately)
        plots.push_back(new PlotOptionsRelRes(
            "inclusive/Q2_RelRes_EM",
            "(Q^{2}_{EM} - Q^{2}_{truth})/Q^{2}_{truth}",
            "Counts",
            -0.15, 0.15,
            "figs/Q2_relres_EM.png"
        ));
        plots.push_back(new PlotOptionsRelRes(
            "inclusive/Q2_RelRes_DA",
            "(Q^{2}_{DA} - Q^{2}_{truth})/Q^{2}_{truth}",
            "Counts",
            -0.15, 0.15,
            "figs/Q2_relres_DA.png"
        ));
        plots.push_back(new PlotOptionsRelRes(
            "inclusive/Q2_RelRes_ESigma",
            "(Q^{2}_{ESigma} - Q^{2}_{truth})/Q^{2}_{truth}",
            "Counts",
            -0.15, 0.15,
            "figs/Q2_relres_ESigma.png"
        ));
        plots.push_back(new PlotOptionsRelRes(
            "inclusive/Q2_RelRes_Sigma",
            "(Q^{2}_{Sigma} - Q^{2}_{truth})/Q^{2}_{truth}",
            "Counts",
            -0.15, 0.15,
            "figs/Q2_relres_Sigma.png"
        ));

        // Truth vs Reco correlations for Q2/x/y
        plots.push_back(new PlotOptionsResponseMatrix(
            "inclusive/Corr_Q2_EM",
            "Q^{2}_{truth} [GeV^{2}]",
            "Q^{2}_{EM} [GeV^{2}]",
            "figs/Corr_Q2_EM.png",
            true, true,
            {3.4, 150.0}, {3.4, 150.0}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "inclusive/Corr_Q2_DA",
            "Q^{2}_{truth} [GeV^{2}]",
            "Q^{2}_{DA} [GeV^{2}]",
            "figs/Corr_Q2_DA.png",
            true, true,
            {3.4, 150.0}, {3.4, 150.0}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "inclusive/Corr_Q2_ESigma",
            "Q^{2}_{truth} [GeV^{2}]",
            "Q^{2}_{ESigma} [GeV^{2}]",
            "figs/Corr_Q2_ESigma.png",
            true, true,
            {3.4, 150.0}, {3.4, 150.0}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "inclusive/Corr_Q2_Sigma",
            "Q^{2}_{truth} [GeV^{2}]",
            "Q^{2}_{Sigma} [GeV^{2}]",
            "figs/Corr_Q2_Sigma.png",
            true, true,
            {3.4, 150.0}, {3.4, 150.0}
        ));

        plots.push_back(new PlotOptionsResponseMatrix(
            "inclusive/Corr_x_EM",
            "x_{truth}",
            "x_{EM}",
            "figs/Corr_x_EM.png",
            true, true,
            {1e-4, 1.0}, {1e-4, 1.0}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "inclusive/Corr_x_DA",
            "x_{truth}",
            "x_{DA}",
            "figs/Corr_x_DA.png",
            true, true,
            {1e-4, 1.0}, {1e-4, 1.0}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "inclusive/Corr_x_ESigma",
            "x_{truth}",
            "x_{ESigma}",
            "figs/Corr_x_ESigma.png",
            true, true,
            {1e-4, 1.0}, {1e-4, 1.0}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "inclusive/Corr_x_Sigma",
            "x_{truth}",
            "x_{Sigma}",
            "figs/Corr_x_Sigma.png",
            true, true,
            {1e-4, 1.0}, {1e-4, 1.0}
        ));

        plots.push_back(new PlotOptionsResponseMatrix(
            "inclusive/Corr_y_EM",
            "y_{truth}",
            "y_{EM}",
            "figs/Corr_y_EM.png",
            false, false,
            {0.0, 1.0}, {0.0, 1.0}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "inclusive/Corr_y_DA",
            "y_{truth}",
            "y_{DA}",
            "figs/Corr_y_DA.png",
            false, false,
            {0.0, 1.0}, {0.0, 1.0}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "inclusive/Corr_y_ESigma",
            "y_{truth}",
            "y_{ESigma}",
            "figs/Corr_y_ESigma.png",
            false, false,
            {0.0, 1.0}, {0.0, 1.0}
        ));
        plots.push_back(new PlotOptionsResponseMatrix(
            "inclusive/Corr_y_Sigma",
            "y_{truth}",
            "y_{Sigma}",
            "figs/Corr_y_Sigma.png",
            false, false,
            {0.0, 1.0}, {0.0, 1.0}
        ));
    }

    // -----------------------------
    // Execute
    // -----------------------------
    gSystem->mkdir("figs", kTRUE);
    for (auto* plot : plots) {
        plot->Plot(inputFile);
    }
    for (auto* plot : plots) {
        delete plot;
    }

    inputFile->Close();
    delete inputFile;

    std::cout << "\nPlotting complete! All plots saved to figs/ directory." << std::endl;
    return 0;
}
