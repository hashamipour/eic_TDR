#ifndef RECO_METHODS_HPP
#define RECO_METHODS_HPP

#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "Math/Vector4D.h"
#include <vector>
#include <string>
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

using P3MVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzMVector>;

// Structure to hold histograms for a reconstruction method
struct MethodHistograms {
    std::string name;
    TH1D* h_t_reco;
    TH2D* h_t_corr;
    TH1D* h_t_res;
    TH1D* h_theta;
    TH1D* h_xL;
    TH2D* h_xL_corr;
    TH1D* h_MX2;
    TH2D* h_MX2_corr;
    std::vector<double> t_truth_vec;
    std::vector<double> t_reco_vec;
    
    MethodHistograms(const std::string& method_name, const std::vector<Double_t>& t_bins);
    TGraph* MakeCorrelationGraph();
    void FillCorrelation(double t_truth, double t_reco);
    void Write();
};

// Structure to hold beam information
struct BeamInfo {
    P3MVector e_beam;
    P3MVector p_beam;
    double fMass_electron = 0.000511;
    double fMass_proton   = 0.938272;
};

// Calculate |t| from momentum transfer: BABE method
inline Double_t CalcT(const P3MVector& p_initial, const P3MVector& p_final) {
    return (p_final - p_initial).M2();
}

// Calculate |t| for eX method: t = (q - PX)^2
inline Double_t CalcT_eX(const P3MVector& q_gamma, const P3MVector& X_system) {
    return (q_gamma - X_system).M2();
}

// Calculate x_L = Pz_proton / E_proton
inline double CalcXL(const P3MVector& proton){
    return proton.Pz() / proton.E();
}

inline double CalcMX2(const P3MVector& X_system){
    return X_system.M2();
}

inline double CalcMX2_LPS(double xL, double xBj, double W){
    return (1.0 - xL*(1.0 + xBj)) * W * W;
}

////////////////////////////////////////////////////////////////////////
// Find and return scattered electron using EICRecon's identification
struct ScatteredElectronInfo {
    P3MVector p4;
    bool found;
    unsigned int index;
};

ScatteredElectronInfo GetScatteredElectron(
    TTreeReaderArray<int>& electron_index_array,
    TTreeReaderArray<float>& rpf_px, TTreeReaderArray<float>& rpf_py,
    TTreeReaderArray<float>& rpf_pz, double electron_mass
);

#endif // RECO_METHODS_HPP