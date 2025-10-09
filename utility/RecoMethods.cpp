#include "RecoMethods.hpp"
#include "TMath.h"
#include <iostream>


MethodHistograms::MethodHistograms(const std::string& method_name, const std::vector<Double_t>& t_bins) 
    : name(method_name) {
    
    TString n(name.c_str());  // Convert to TString once
    TString hist_suffix = (name == "B0") ? "" : "_histo";
    
    h_t_reco = new TH1D(
        ("t_" + n + hist_suffix).Data(),
        (n + " Reco Mandelstam t;|t| [GeV^{2}];Counts").Data(),
        t_bins.size()-1, t_bins.data()
    );
    
    h_t_corr = new TH2D(
        ("t_corr_" + n).Data(),
        (n + ": Truth vs Reco;Truth |t| [GeV^{2}];Reco |t| [GeV^{2}]").Data(),
        t_bins.size()-1, t_bins.data(), t_bins.size()-1, t_bins.data()
    );
    
    h_t_res = new TH1D(
        ("t_res_" + n).Data(),
        (n + " |t| Resolution;(|t|_{reco} - |t|_{truth})/|t|_{truth};Counts").Data(),
        101, -0.5, 0.5
    );
    
    h_theta = new TH1D(
        ("theta_" + n).Data(),
        (n + " Scattering Angle;#theta [mrad];Counts").Data(),
        100, 0.0, 25.0
    );

    h_xL = new TH1D(Form("xL_%s", name.data()), 
                Form("%s x_{L};x_{L};Counts", name.data()), 100, 0.75, 1.05);
    
    h_xL_corr = new TH2D(Form("xL_corr_%s", name.data()),
                     Form("%s x_{L} Correlation;Truth x_{L};Reco x_{L}", name.data()),
                     10, 0.75, 1.05, 10, 0.75, 1.05);

    h_MX2 = new TH1D(Form("MX2_%s", name.data()), 
                 Form("%s M_{X}^{2};M_{X}^{2} [GeV^{2}];Counts", name.data()), 
                 100, 0, 200);

    h_MX2_corr = new TH2D(Form("MX2_corr_%s", name.data()),
                      Form("%s M_{X}^{2} Correlation;Truth M_{X}^{2} [GeV^{2}];LPS M_{X}^{2} [GeV^{2}]", name.data()),
                      100, 0, 100, 100, 0, 100);
}

TGraph* MethodHistograms::MakeCorrelationGraph() {
    if(t_truth_vec.empty()) return nullptr;
    
    TString n(name.c_str());
    TGraph* g = new TGraph(t_truth_vec.size(), t_truth_vec.data(), t_reco_vec.data());
    g->SetName(("t_corr_" + n + "_graph").Data());
    g->SetTitle((n + ": Truth vs Reco |t|;Truth |t| [GeV^{2}];Reco |t| [GeV^{2}]").Data());
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.5);
    return g;
}

void MethodHistograms::FillCorrelation(double t_truth, double t_reco) {
    double abs_t_truth = TMath::Abs(t_truth);
    double abs_t_reco = TMath::Abs(t_reco);
    
    h_t_reco->Fill(abs_t_reco);
    h_t_corr->Fill(abs_t_truth, abs_t_reco);
    t_truth_vec.push_back(abs_t_truth);
    t_reco_vec.push_back(abs_t_reco);
    
    if(abs_t_truth > 1e-6) {
        h_t_res->Fill((abs_t_reco - abs_t_truth) / abs_t_truth);
    }
}

void MethodHistograms::Write() {
    h_t_reco->Write();
    h_t_corr->Write();
    h_t_res->Write();
    
    h_theta->Write();

    h_xL->Write();
    h_xL_corr->Write();
    h_MX2->Write();
    h_MX2_corr->Write();
    
    TGraph* g = MakeCorrelationGraph();
    if(g) {
        g->Write();
        delete g;
    }
}



ScatteredElectronInfo GetScatteredElectron(
    TTreeReaderArray<int>& electron_index_array,
    TTreeReaderArray<float>& rpf_px, TTreeReaderArray<float>& rpf_py,
    TTreeReaderArray<float>& rpf_pz, double electron_mass
) {
    ScatteredElectronInfo info;
    info.found = false;
    info.index = 0;
    
    if(electron_index_array.GetSize() == 0) return info;
    
    unsigned int e_idx = electron_index_array[0];
    
    if(e_idx >= (unsigned)rpf_px.GetSize()) return info;
    
    info.p4.SetCoordinates(rpf_px[e_idx], rpf_py[e_idx], rpf_pz[e_idx], electron_mass);
    info.found = true;
    info.index = e_idx;
    
    return info;
}