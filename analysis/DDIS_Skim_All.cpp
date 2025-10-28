// Combined Analysis: Proton Mandelstam t + Q2/x/y Kinematic Variables
// Compile: g++ DDIS_Skim_Full.cpp -o DDIS_Skim_Full $(root-config --cflags --glibs)
// Usage: ./DDIS_Skim_Full filelist.txt

#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <filesystem>

#include "TFile.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/GenVector/Boost.h"
#include "Math/VectorUtil.h"

#include "Utility.hpp"
#include "DDIS_Util.hpp"
#include "RecoMethods.hpp"

using std::cout; using std::endl;
using ROOT::Math::VectorUtil::boost;
using ROOT::Math::VectorUtil::Angle;
using ROOT::Math::RotationX;
using ROOT::Math::RotationY;
using P3EVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;
using P3MVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzMVector>;
using MomVector = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>, ROOT::Math::DefaultCoordinateSystemTag>;

// Global afterburner correction parameters
Float_t fXAngle{-0.025};
RotationX rotAboutX;
RotationY rotAboutY;
MomVector vBoostToCoM;
MomVector vBoostToHoF;

const Float_t fMass_proton{0.938272};
const Float_t fMass_electron{0.000511};

void undoAfterburnAndCalc(P3MVector& p, P3MVector& k){
    P3MVector p_beam(fXAngle*p.E(), 0., p.E(), p.M());
    P3MVector e_beam(0., 0., -k.E(), k.M());
    
    P3MVector CoM_boost = p_beam + e_beam;
    vBoostToCoM.SetXYZ(-CoM_boost.X()/CoM_boost.E(), -CoM_boost.Y()/CoM_boost.E(), -CoM_boost.Z()/CoM_boost.E());
    
    p_beam = boost(p_beam, vBoostToCoM);
    e_beam = boost(e_beam, vBoostToCoM);
    
    Float_t fRotY = -1.0*TMath::ATan2(p_beam.X(), p_beam.Z());
    Float_t fRotX = 1.0*TMath::ATan2(p_beam.Y(), p_beam.Z());
    
    rotAboutY = RotationY(fRotY);
    rotAboutX = RotationX(fRotX);
    
    p_beam = rotAboutY(p_beam);
    p_beam = rotAboutX(p_beam);
    e_beam = rotAboutY(e_beam);
    e_beam = rotAboutX(e_beam);
    
    P3EVector HoF_boost(0., 0., CoM_boost.Z(), CoM_boost.E());
    vBoostToHoF.SetXYZ(HoF_boost.X()/HoF_boost.E(), HoF_boost.Y()/HoF_boost.E(), HoF_boost.Z()/HoF_boost.E());
    
    p_beam = boost(p_beam, vBoostToHoF);
    e_beam = boost(e_beam, vBoostToHoF);
    
    p.SetPxPyPzE(p_beam.X(), p_beam.Y(), p_beam.Z(), p_beam.E());
    k.SetPxPyPzE(e_beam.X(), e_beam.Y(), e_beam.Z(), e_beam.E());
}

void undoAfterburn(P3MVector& a){
    a = boost(a, vBoostToCoM);
    a = rotAboutY(a);
    a = rotAboutX(a);
    a = boost(a, vBoostToHoF);
}

// Process truth protons from MC
std::vector<P3MVector> ProcessTruthProtons(
    TTreeReaderArray<double>& mc_px, TTreeReaderArray<double>& mc_py, TTreeReaderArray<double>& mc_pz,
    TTreeReaderArray<double>& mc_mass, TTreeReaderArray<int>& mc_gen_status, TTreeReaderArray<int>& mc_pdg,
    const BeamInfo& beams, TH1D* h_theta, TH1D* h_t, TH1D* h_xL
) {
    std::vector<P3MVector> truth_protons;
    
    for(int i = 0; i < mc_px.GetSize(); i++){
        if(mc_gen_status[i] == 1 && mc_pdg[i] == 2212){
            P3MVector p(mc_px[i], mc_py[i], mc_pz[i], mc_mass[i]);
            undoAfterburn(p);
            
            truth_protons.push_back(p);
            h_theta->Fill(p.Theta() * 1000.0);
            h_t->Fill(TMath::Abs(CalcT(beams.p_beam, p)));
            
            // Calculate x_L = |P'|/|P|
            double xL = p.P() / beams.p_beam.P();
            h_xL->Fill(xL);
        }
    }
    
    return truth_protons;
}

// Calculate truth M_X^2 from all final state particles except scattered e- and p
P3MVector CalculateTruthXSystem(
    TTreeReaderArray<double>& mc_px, TTreeReaderArray<double>& mc_py, TTreeReaderArray<double>& mc_pz,
    TTreeReaderArray<double>& mc_mass, TTreeReaderArray<int>& mc_gen_status, TTreeReaderArray<int>& mc_pdg
) {
    P3MVector X_system(0,0,0,0);
    
    for(int i = 0; i < mc_px.GetSize(); i++){
        if(mc_gen_status[i] != 1) continue;  // Only final state
        if(mc_pdg[i] == 11 || mc_pdg[i] == 2212) continue;  // Exclude scattered e- and p
        
        P3MVector p(mc_px[i], mc_py[i], mc_pz[i], mc_mass[i]);
        undoAfterburn(p);
        X_system += p;
    }
    
    return X_system;
}

// Calculate sum(E - pz) for hadronic final state (MC Truth)
double CalculateSumEPz_Truth(
    TTreeReaderArray<double>& mc_px,
    TTreeReaderArray<double>& mc_py,
    TTreeReaderArray<double>& mc_pz,
    TTreeReaderArray<double>& mc_mass,
    TTreeReaderArray<int>& mc_genStatus,
    TTreeReaderArray<int>& mc_pdg
) {
    double sumEPz = 0.0;
    
    for(int i = 0; i < mc_px.GetSize(); i++){
        if(mc_genStatus[i] != 1) continue;
        if(mc_pdg[i] == 11 || mc_pdg[i] == 2212) continue;
        
        double px = mc_px[i];
        double py = mc_py[i];
        double pz = mc_pz[i];
        double m = mc_mass[i];
        double E = TMath::Sqrt(px*px + py*py + pz*pz + m*m);
        
        sumEPz += (E - pz);
    }
    
    return sumEPz;
}

// Calculate sum(E - pz) for hadronic final state (Reconstructed)
double CalculateSumEPz_Reco(
    TTreeReaderArray<float>& re_px,
    TTreeReaderArray<float>& re_py,
    TTreeReaderArray<float>& re_pz,
    TTreeReaderArray<float>& re_energy,
    TTreeReaderArray<int>& electron_scat_index
) {
    double sumEPz = 0.0;
    
    int scat_e_idx = -1;
    if(electron_scat_index.GetSize() > 0) {
        scat_e_idx = electron_scat_index[0];
    }
    
    for(int i = 0; i < re_px.GetSize(); i++){
        if(i == scat_e_idx) continue;
        
        double px = re_px[i];
        double py = re_py[i];
        double pz = re_pz[i];
        double E = re_energy[i];
        
        sumEPz += (E - pz);
    }
    
    return sumEPz;
}


void runFullAnalysis(TString fileList){
    cout << "====================================" << endl;
    cout << " Combined DDIS Analysis" << endl;
    cout << " - Proton Mandelstam t" << endl;
    cout << " - Q2/x/y Kinematics" << endl;
    cout << "====================================" << endl;
    cout << "\nInput filelist: " << fileList << endl;

    BeamInfo beams;
    
    // Setup input chain
    TChain* events = new TChain("events");
    Int_t nFiles = 0;
    
    std::ifstream fileListStream(fileList.Data());
    std::string fileName;
    while(getline(fileListStream, fileName)){
        if (!std::filesystem::exists(fileName)) {
            std::cerr << "Warning: File does not exist: " << fileName << std::endl;
            continue;
        }
        events->Add((TString)fileName);
        nFiles++;
    }
    
    cout << "\nNo. of files: " << nFiles << "; no. of events: " << events->GetEntries() << endl;
    
    //=================================================================
    // SETUP HISTOGRAMS - PROTON MANDELSTAM T ANALYSIS
    //=================================================================
    
    // Truth histograms - combined log/linear binning
    std::vector<Double_t> t_bins_log = GetLogBins(1e-4, 0.5, 30);
    std::vector<Double_t> t_bins_lin = GetLinBins(0.5, 1.6, 11);
    std::vector<Double_t> t_bins_MC;
    t_bins_MC.insert(t_bins_MC.end(), t_bins_log.begin(), t_bins_log.end());
    t_bins_MC.insert(t_bins_MC.end(), t_bins_lin.begin()+1, t_bins_lin.end());

    TH1D* h_t_MC = new TH1D("t_MC", "Truth Mandelstam t;|t| [GeV^{2}];Counts", 
                            t_bins_MC.size()-1, t_bins_MC.data());
    TH1D* h_theta_MC = new TH1D("theta_MC", "MC Proton Scattering Angle;#theta [mrad];Counts", 
                                100, 0.0, 25.0);
    TH1D* h_xL_MC = new TH1D("xL_MC", "Truth x_{L};x_{L};Counts", 200, 0.92, 1.02);
    TH1D* h_MX2_MC = new TH1D("MX2_MC", "Truth M_{X}^{2};M_{X}^{2} [GeV^{2}];Counts", 100, 0, 100);
    TH1D* h_xP_MC = new TH1D("xP_MC", "Truth x_{P};x_{P};Counts", 100, 0, 0.1);
    
    // B0 histograms
    TH1D* h_theta_B0 = new TH1D("theta_B0", "B0 Proton Scattering Angle;#theta [mrad];Counts", 
                                100, 0.0, 25.0);
    std::vector<Double_t> t_bins_B0 = t_bins_lin;
    TH1D* h_t_B0 = new TH1D("t_B0", "B0 Mandelstam t;|t| [GeV^{2}];Counts", 
                            t_bins_B0.size()-1, t_bins_B0.data());
    TH2D* h_t_corr_B0 = new TH2D("t_corr_B0", "B0 t Correlation;Truth |t| [GeV^{2}];Reco |t| [GeV^{2}]",
                                 t_bins_B0.size()-1, t_bins_B0.data(),t_bins_B0.size()-1, t_bins_B0.data());
    TH1D* h_xL_B0 = new TH1D("xL_B0", "B0 x_{L};x_{L};Counts", 200, 0.92, 1.02);
    TH2D* h_xL_corr_B0 = new TH2D("xL_corr_B0", "B0 x_{L} Correlation;Truth x_{L};Reco x_{L}",
                                  200, 0.92, 1.02, 200, 0.92, 1.02);
    TH1D* h_MX2_B0 = new TH1D("MX2_B0", "B0 M_{X}^{2};M_{X}^{2} [GeV^{2}];Counts", 100, 0, 100);
    TH2D* h_MX2_corr_B0 = new TH2D("MX2_corr_B0", "B0 M_{X}^{2} Correlation;Truth M_{X}^{2} [GeV^{2}];Reco M_{X}^{2} [GeV^{2}]",
                                   100, 0, 100, 100, 0, 100);
    TH1D* h_xP_B0 = new TH1D("xP_B0", "B0 x_{P};x_{P};Counts", 100, 0, 0.1);
    TH2D* h_xP_corr_B0 = new TH2D("xP_corr_B0", "B0 x_{P} Correlation;Truth x_{P};Reco x_{P}",
                                  100, 0, 0.1, 100, 0, 0.1);
    TH1D* h_t_res_B0 = new TH1D("t_res_B0", "B0 t Resolution;(|t|_{reco} - |t|_{truth})/|t|_{truth};Counts",
                                100, -.5, .5);
    
    // RP histograms
    TH1D* h_theta_RP = new TH1D("theta_RP", "RP Proton Scattering Angle;#theta [mrad];Counts", 
                                100, 0.0, 25.0);
    std::vector<Double_t> t_bins_RP = GetLogBins(1e-4, 0.5, 30);
    TH1D* h_t_RP = new TH1D("t_RP", "RP Mandelstam t;|t| [GeV^{2}];Counts", 
                            t_bins_RP.size()-1, t_bins_RP.data());
    TH2D* h_t_corr_RP = new TH2D("t_corr_RP", "RP t Correlation;Truth |t| [GeV^{2}];Reco |t| [GeV^{2}]",
                                 t_bins_RP.size()-1, t_bins_RP.data(),t_bins_RP.size()-1, t_bins_RP.data());
    TH1D* h_xL_RP = new TH1D("xL_RP", "RP x_{L};x_{L};Counts", 200, 0.92, 1.02);
    TH2D* h_xL_corr_RP = new TH2D("xL_corr_RP", "RP x_{L} Correlation;Truth x_{L};Reco x_{L}",
                                  200, 0.92, 1.02, 200, 0.92, 1.02);
    TH1D* h_MX2_RP = new TH1D("MX2_RP", "RP M_{X}^{2};M_{X}^{2} [GeV^{2}];Counts", 100, 0, 100);
    TH2D* h_MX2_corr_RP = new TH2D("MX2_corr_RP", "RP M_{X}^{2} Correlation;Truth M_{X}^{2} [GeV^{2}];Reco M_{X}^{2} [GeV^{2}]",
                                   100, 0, 100, 100, 0, 100);
    TH1D* h_xP_RP = new TH1D("xP_RP", "RP x_{P};x_{P};Counts", 100, 0, 0.1);
    TH2D* h_xP_corr_RP = new TH2D("xP_corr_RP", "RP x_{P} Correlation;Truth x_{P};Reco x_{P}",
                                  100, 0, 0.1, 100, 0, 0.1);
    TH1D* h_t_res_RP = new TH1D("t_res_RP", "RP t Resolution;(|t|_{reco} - |t|_{truth})/|t|_{truth};Counts",
                                100, -.5, .5);
    
    // Sum(E-pz) and M_X histograms (from Combined version)
    TH1D* h_SumEPz_truth = new TH1D("SumEPz_truth", "Sum(E-pz) Truth;Sum(E-p_{z}) [GeV];Counts", 100, 0, 50);
    TH1D* h_SumEPz_reco = new TH1D("SumEPz_reco", "Sum(E-pz) Reco;Sum(E-p_{z}) [GeV];Counts", 100, 0, 50);
    TH1D* h_Mx_truth = new TH1D("Mx_truth", "M_{X} (4-vector) Truth;M_{X} [GeV];Counts", 100, 0, 20);
    TH1D* h_Mx_reco = new TH1D("Mx_reco", "M_{X} (4-vector) Reco;M_{X} [GeV];Counts", 100, 0, 20);
    
    //=================================================================
    // SETUP HISTOGRAMS - Q2/x/y ANALYSIS
    //=================================================================
    
    int n_bins_Q2 = 10;
    std::vector<Double_t> bin_edges_Q2 = GetRoundedLogBins(3.4, 150.0, n_bins_Q2);
    n_bins_Q2 = bin_edges_Q2.size()-1;
    
    // x histograms
    std::vector<Double_t> x_bins = GetLogBins(1.0e-4, 1.0, 25);
    TH1D* h_x_EM     = new TH1D("x_EM",     "electron method;x_{Bj}", x_bins.size()-1, x_bins.data());
    TH1D* h_x_DA     = new TH1D("x_DA",     "DA method;x_{Bj}",       x_bins.size()-1, x_bins.data());
    TH1D* h_x_ESigma = new TH1D("x_ESigma", "ESigma method;x_{Bj}",   x_bins.size()-1, x_bins.data());
    TH1D* h_x_Sigma  = new TH1D("x_Sigma",  "Sigma method;x_{Bj}",    x_bins.size()-1, x_bins.data());
    TH1D* h_x_truth  = new TH1D("x_truth",  "truth;x_{Bj}",           x_bins.size()-1, x_bins.data());

    // y histograms
    int n_y_bins = 20;
    TH1D* h_y_EM     = new TH1D("y_EM",     "electron method;y", n_y_bins, 0.0, 1.0);
    TH1D* h_y_DA     = new TH1D("y_DA",     "DA method;y",       n_y_bins, 0.0, 1.0);
    TH1D* h_y_ESigma = new TH1D("y_ESigma", "ESigma method;y",   n_y_bins, 0.0, 1.0);
    TH1D* h_y_Sigma  = new TH1D("y_Sigma",  "Sigma method;y",    n_y_bins, 0.0, 1.0);
    TH1D* h_y_truth  = new TH1D("y_truth",  "truth;y",           n_y_bins, 0.0, 1.0);

    // Relative resolution for x and y
    TH1D* h_RelRes_x_EM     = new TH1D("x_RelRes_EM",     "electron method;#frac{x(Reco)-x(MC)}{x(MC)}", 101, -0.15, 0.15);
    TH1D* h_RelRes_x_DA     = new TH1D("x_RelRes_DA",     "DA method;#frac{x(Reco)-x(MC)}{x(MC)}",       101, -0.15, 0.15);
    TH1D* h_RelRes_x_ESigma = new TH1D("x_RelRes_ESigma", "ESigma method;#frac{x(Reco)-x(MC)}{x(MC)}",   101, -0.15, 0.15);
    TH1D* h_RelRes_x_Sigma  = new TH1D("x_RelRes_Sigma",  "Sigma method;#frac{x(Reco)-x(MC)}{x(MC)}",    101, -0.15, 0.15);

    TH1D* h_RelRes_y_EM     = new TH1D("y_RelRes_EM",     "electron method;#frac{y(Reco)-y(MC)}{y(MC)}", 101, -0.15, 0.15);
    TH1D* h_RelRes_y_DA     = new TH1D("y_RelRes_DA",     "DA method;#frac{y(Reco)-y(MC)}{y(MC)}",       101, -0.15, 0.15);
    TH1D* h_RelRes_y_ESigma = new TH1D("y_RelRes_ESigma", "ESigma method;#frac{y(Reco)-y(MC)}{y(MC)}",   101, -0.15, 0.15);
    TH1D* h_RelRes_y_Sigma  = new TH1D("y_RelRes_Sigma",  "Sigma method;#frac{y(Reco)-y(MC)}{y(MC)}",    101, -0.15, 0.15);

    // 2D binned relres vs truth
    int n_binned = 51;
    TH2D* h_RelRes_x_binned_EM     = new TH2D("x_RelRes_binned_EM",     "x: truth vs rel. res (EM);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",       
        x_bins.size()-1, x_bins.data(), n_binned, -0.15, 0.15);
    TH2D* h_RelRes_x_binned_DA     = new TH2D("x_RelRes_binned_DA",     "x: truth vs rel. res (DA);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",       
        x_bins.size()-1, x_bins.data(), n_binned, -0.15, 0.15);
    TH2D* h_RelRes_x_binned_ESigma = new TH2D("x_RelRes_binned_ESigma", "x: truth vs rel. res (ESigma);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",   
        x_bins.size()-1, x_bins.data(), n_binned, -0.15, 0.15);
    TH2D* h_RelRes_x_binned_Sigma  = new TH2D("x_RelRes_binned_Sigma",  "x: truth vs rel. res (Sigma);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",    
        x_bins.size()-1, x_bins.data(), n_binned, -0.15, 0.15);

    TH2D* h_RelRes_y_binned_EM     = new TH2D("y_RelRes_binned_EM",     "y: truth vs rel. res (EM);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",       n_y_bins, 0.0, 1.0, n_binned, -0.15, 0.15);
    TH2D* h_RelRes_y_binned_DA     = new TH2D("y_RelRes_binned_DA",     "y: truth vs rel. res (DA);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",       n_y_bins, 0.0, 1.0, n_binned, -0.15, 0.15);
    TH2D* h_RelRes_y_binned_ESigma = new TH2D("y_RelRes_binned_ESigma", "y: truth vs rel. res (ESigma);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",   n_y_bins, 0.0, 1.0, n_binned, -0.15, 0.15);
    TH2D* h_RelRes_y_binned_Sigma  = new TH2D("y_RelRes_binned_Sigma",  "y: truth vs rel. res (Sigma);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",    n_y_bins, 0.0, 1.0, n_binned, -0.15, 0.15);

    // Correlation plots (truth vs reco)
    TH2D* h_Corr_x_EM     = new TH2D("x_Corr_EM",     "x correlation (EM);x_{truth};x_{EM}",           x_bins.size()-1, x_bins.data(), x_bins.size()-1, x_bins.data());
    TH2D* h_Corr_x_DA     = new TH2D("x_Corr_DA",     "x correlation (DA);x_{truth};x_{DA}",           x_bins.size()-1, x_bins.data(), x_bins.size()-1, x_bins.data());
    TH2D* h_Corr_x_ESigma = new TH2D("x_Corr_ESigma", "x correlation (ESigma);x_{truth};x_{ESigma}",   x_bins.size()-1, x_bins.data(), x_bins.size()-1, x_bins.data());
    TH2D* h_Corr_x_Sigma  = new TH2D("x_Corr_Sigma",  "x correlation (Sigma);x_{truth};x_{Sigma}",     x_bins.size()-1, x_bins.data(), x_bins.size()-1, x_bins.data());

    TH2D* h_Corr_y_EM     = new TH2D("y_Corr_EM",     "y correlation (EM);y_{truth};y_{EM}",           n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);
    TH2D* h_Corr_y_DA     = new TH2D("y_Corr_DA",     "y correlation (DA);y_{truth};y_{DA}",           n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);
    TH2D* h_Corr_y_ESigma = new TH2D("y_Corr_ESigma", "y correlation (ESigma);y_{truth};y_{ESigma}",   n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);
    TH2D* h_Corr_y_Sigma  = new TH2D("y_Corr_Sigma",  "y correlation (Sigma);y_{truth};y_{Sigma}",     n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);

    // Q2 histograms
    TH1D* h_RelRes_Q2_EM     = new TH1D("Q2_RelRes_EM",     "electron method;#frac{Q^{2}(Reco)-Q^{2}(MC)}{Q^{2}(MC)}",  101,-0.15,0.15);
    TH1D* h_RelRes_Q2_DA     = new TH1D("Q2_RelRes_DA",     "DA method;#frac{Q^{2}(DA)-Q^{2}(MC)}{Q^{2}(MC)}",          101,-0.15,0.15);
    TH1D* h_RelRes_Q2_ESigma = new TH1D("Q2_RelRes_ESigma", "ESigma method;#frac{Q^{2}(ESigma)-Q^{2}(MC)}{Q^{2}(MC)}", 101,-0.15,0.15);
    TH1D* h_RelRes_Q2_Sigma  = new TH1D("Q2_RelRes_Sigma",  "Sigma method;#frac{Q^{2}(Sigma)-Q^{2}(MC)}{Q^{2}(MC)}",    101,-0.15,0.15);
    
    TH2D* h_RelRes_Q2_binned_EM     = new TH2D("Q2_RelRes_binned_EM",     ";Q^{2} [GeV^{2}];#frac{Q^{2}(EM)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins_Q2, bin_edges_Q2.data(), n_binned,-0.15,0.15);
    TH2D* h_RelRes_Q2_binned_DA     = new TH2D("Q2_RelRes_binned_DA",     ";Q^{2} [GeV^{2}];#frac{Q^{2}(DA)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins_Q2, bin_edges_Q2.data(), n_binned,-0.15,0.15);
    TH2D* h_RelRes_Q2_binned_ESigma = new TH2D("Q2_RelRes_binned_ESigma", ";Q^{2} [GeV^{2}];#frac{Q^{2}(ESigma)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins_Q2, bin_edges_Q2.data(), n_binned,-0.15,0.15);
    TH2D* h_RelRes_Q2_binned_Sigma  = new TH2D("Q2_RelRes_binned_Sigma",  ";Q^{2} [GeV^{2}];#frac{Q^{2}(Sigma)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins_Q2, bin_edges_Q2.data(), n_binned,-0.15,0.15);

    TH2D* h_Corr_Q2_EM     = new TH2D("Corr_Q2_EM",     ";Q^{2}_{MC};Q^{2}_{EM}",
                                  n_bins_Q2, bin_edges_Q2.data(), n_bins_Q2, bin_edges_Q2.data());
    TH2D* h_Corr_Q2_DA     = new TH2D("Corr_Q2_DA",     ";Q^{2}_{MC};Q^{2}_{DA}",
                                    n_bins_Q2, bin_edges_Q2.data(), n_bins_Q2, bin_edges_Q2.data());
    TH2D* h_Corr_Q2_ESigma = new TH2D("Corr_Q2_ESigma", ";Q^{2}_{MC};Q^{2}_{ESigma}",
                                    n_bins_Q2, bin_edges_Q2.data(), n_bins_Q2, bin_edges_Q2.data());
    TH2D* h_Corr_Q2_Sigma  = new TH2D("Corr_Q2_Sigma",  ";Q^{2}_{MC};Q^{2}_{Sigma}",
                                    n_bins_Q2, bin_edges_Q2.data(), n_bins_Q2, bin_edges_Q2.data());

    TH1D* h_Q2_truth  = new TH1D("h_Q2_truth",  "Q^2;# of events", n_bins_Q2, bin_edges_Q2.data());
    TH1D* h_Q2_EM     = new TH1D("h_Q2_EM",     ";Q^{2}",          n_bins_Q2, bin_edges_Q2.data());
    TH1D* h_Q2_DA     = new TH1D("h_Q2_DA",     ";Q^{2}",          n_bins_Q2, bin_edges_Q2.data());
    TH1D* h_Q2_ESigma = new TH1D("h_Q2_ESigma", ";Q^{2}",          n_bins_Q2, bin_edges_Q2.data());
    TH1D* h_Q2_Sigma  = new TH1D("h_Q2_Sigma",  ";Q^{2}",          n_bins_Q2, bin_edges_Q2.data());

    
    //=================================================================
    // PRINT RESULTS AND SAVE TO FILE
    //=================================================================
    
    cout << "\n=== PROTON ANALYSIS RESULTS ===" << endl;
    cout << "Truth protons found: " << n_truth_protons << endl;
    cout << "B0 matched protons: " << n_b0_matches << endl;
    cout << "RP matched protons: " << n_rp_matches << endl;
    if(n_truth_protons > 0){
        cout << "B0 matching efficiency: " << (double)n_b0_matches/n_truth_protons*100 << "%" << endl;
        cout << "RP matching efficiency: " << (double)n_rp_matches/n_truth_protons*100 << "%" << endl;
    }
    
    // Save to file
    gStyle->SetOptStat(0);
    
    TFile* outfile = new TFile("DDIS_Skim_Full_output.root", "RECREATE");
    
    cout << "\nWriting histograms to output file..." << endl;
    
    // Write proton t analysis histograms
    h_t_MC->Write();
    h_theta_MC->Write();
    h_xL_MC->Write();
    h_MX2_MC->Write();
    h_xP_MC->Write();
    
    h_theta_B0->Write();
    h_t_B0->Write();
    h_t_corr_B0->Write();
    h_t_res_B0->Write();
    h_xL_B0->Write();
    h_xL_corr_B0->Write();
    h_MX2_B0->Write();
    h_MX2_corr_B0->Write();
    h_xP_B0->Write();
    h_xP_corr_B0->Write();
    
    h_theta_RP->Write();
    h_t_RP->Write();
    h_t_corr_RP->Write();
    h_t_res_RP->Write();
    h_xL_RP->Write();
    h_xL_corr_RP->Write();
    h_MX2_RP->Write();
    h_MX2_corr_RP->Write();
    h_xP_RP->Write();
    h_xP_corr_RP->Write();
    
    // Write Sum(E-pz) and M_X histograms
    h_SumEPz_truth->Write();
    h_SumEPz_reco->Write();
    h_Mx_truth->Write();
    h_Mx_reco->Write();
    
    // Write Q2/x/y analysis histograms
    h_x_truth->Write();
    h_x_EM->Write();
    h_x_DA->Write();
    h_x_ESigma->Write();
    h_x_Sigma->Write();
    
    h_y_truth->Write();
    h_y_EM->Write();
    h_y_DA->Write();
    h_y_ESigma->Write();
    h_y_Sigma->Write();
    
    h_RelRes_x_EM->Write();
    h_RelRes_x_DA->Write();
    h_RelRes_x_ESigma->Write();
    h_RelRes_x_Sigma->Write();
    
    h_RelRes_y_EM->Write();
    h_RelRes_y_DA->Write();
    h_RelRes_y_ESigma->Write();
    h_RelRes_y_Sigma->Write();
    
    h_RelRes_x_binned_EM->Write();
    h_RelRes_x_binned_DA->Write();
    h_RelRes_x_binned_ESigma->Write();
    h_RelRes_x_binned_Sigma->Write();
    
    h_RelRes_y_binned_EM->Write();
    h_RelRes_y_binned_DA->Write();
    h_RelRes_y_binned_ESigma->Write();
    h_RelRes_y_binned_Sigma->Write();
    
    h_Corr_x_EM->Write();
    h_Corr_x_DA->Write();
    h_Corr_x_ESigma->Write();
    h_Corr_x_Sigma->Write();
    
    h_Corr_y_EM->Write();
    h_Corr_y_DA->Write();
    h_Corr_y_ESigma->Write();
    h_Corr_y_Sigma->Write();
    
    h_Q2_truth->Write();
    h_Q2_EM->Write();
    h_Q2_DA->Write();
    h_Q2_ESigma->Write();
    h_Q2_Sigma->Write();
    
    h_RelRes_Q2_EM->Write();
    h_RelRes_Q2_DA->Write();
    h_RelRes_Q2_ESigma->Write();
    h_RelRes_Q2_Sigma->Write();
    
    h_RelRes_Q2_binned_EM->Write();
    h_RelRes_Q2_binned_DA->Write();
    h_RelRes_Q2_binned_ESigma->Write();
    h_RelRes_Q2_binned_Sigma->Write();
    
    h_Corr_Q2_EM->Write();
    h_Corr_Q2_DA->Write();
    h_Corr_Q2_ESigma->Write();
    h_Corr_Q2_Sigma->Write();
    
    outfile->Close();
    
    cout << "\n=== ANALYSIS COMPLETE ===" << endl;
    cout << "Output saved to: DDIS_Skim_Full_output.root" << endl;
    cout << "\nHistograms included:" << endl;
    cout << "  - Proton Mandelstam t analysis (MC truth, B0, RP)" << endl;
    cout << "  - M_X^2 calculated using: M_X^2 = s*y*(1 - x_L - x_Bj)" << endl;
    cout << "  - x_P calculated using: x_P = (M_X^2 + Q^2 - t) / (W^2 + Q^2 - M_p^2)" << endl;
    cout << "  - Sum(E-pz) and M_X (4-vector) histograms" << endl;
    cout << "  - Q2/x/y kinematics (EM, DA, ESigma, Sigma methods)" << endl;
    cout << "  - Resolutions and correlations for all variables" << endl;
}

int main(int argc, char** argv){
    if(argc < 2){
        std::cerr << "Usage: " << argv[0] << " <filelist.txt>\n";
        return 1;
    }
    
    runFullAnalysis(argv[1]);
    return 0;
}
hAddress("InclusiveKinematicsElectron.y", &electron_y_EM);
    events->SetBranchAddress("InclusiveKinematicsDA.y",       &electron_y_DA);
    events->SetBranchAddress("InclusiveKinematicsESigma.y",   &electron_y_ESigma);
    events->SetBranchAddress("InclusiveKinematicsSigma.y",    &electron_y_Sigma);
    events->SetBranchAddress("InclusiveKinematicsTruth.y",    &electron_y_truth);
    
    //=================================================================
    // EVENT LOOP
    //=================================================================
    
    const double s = beams.s;
    const double Mp2 = fMass_proton * fMass_proton;
    
    int n_truth_protons = 0;
    int n_b0_matches = 0;
    int n_rp_matches = 0;
    
    Long64_t nentries = events->GetEntries();
    cout << "\nProcessing " << nentries << " events..." << endl;
    
    for(Long64_t i = 0; i < nentries; i++){
        if(i % 1000 == 0){
            printf("\rProcessing event %lld of %lld; %.2f%% done.", i, nentries, 100.0*i/nentries);
            fflush(stdout);
        }
        
        tree_reader.Next();
        events->GetEntry(i);
        
        // Process truth protons
        std::vector<P3MVector> truth_protons = ProcessTruthProtons(
            mc_px_array, mc_py_array, mc_pz_array, mc_mass_array,
            mc_genStatus_array, mc_pdg_array, beams,
            h_theta_MC, h_t_MC, h_xL_MC
        );
        
        n_truth_protons += truth_protons.size();
        
        // Calculate truth M_X^2
        P3MVector X_truth = CalculateTruthXSystem(
            mc_px_array, mc_py_array, mc_pz_array, mc_mass_array,
            mc_genStatus_array, mc_pdg_array
        );
        double MX2_truth = CalcMX2(X_truth);
        h_MX2_MC->Fill(MX2_truth);
        
        // Calculate Sum(E-pz) - truth
        double sumEPz_truth = CalculateSumEPz_Truth(
            mc_px_array, mc_py_array, mc_pz_array, mc_mass_array,
            mc_genStatus_array, mc_pdg_array
        );
        h_SumEPz_truth->Fill(sumEPz_truth);
        
        // Calculate M_X (4-vector) - truth
        h_Mx_truth->Fill(X_truth.M());
        
        // Calculate Sum(E-pz) - reco
        double sumEPz_reco = CalculateSumEPz_Reco(
            re_px_array, re_py_array, re_pz_array, re_energy_array,
            electron_scat_index
        );
        h_SumEPz_reco->Fill(sumEPz_reco);
        
        // Fill Q2/x/y histograms
        h_x_EM->Fill(electron_x_EM);
        h_x_DA->Fill(electron_x_DA);
        h_x_ESigma->Fill(electron_x_ESigma);
        h_x_Sigma->Fill(electron_x_Sigma);
        h_x_truth->Fill(electron_x_truth);
        
        h_y_EM->Fill(electron_y_EM);
        h_y_DA->Fill(electron_y_DA);
        h_y_ESigma->Fill(electron_y_ESigma);
        h_y_Sigma->Fill(electron_y_Sigma);
        h_y_truth->Fill(electron_y_truth);
        
        // Relative resolutions
        if(electron_x_truth != 0.0f){
            h_RelRes_x_EM->Fill((electron_x_EM - electron_x_truth) / electron_x_truth);
            h_RelRes_x_DA->Fill((electron_x_DA - electron_x_truth) / electron_x_truth);
            h_RelRes_x_ESigma->Fill((electron_x_ESigma - electron_x_truth) / electron_x_truth);
            h_RelRes_x_Sigma->Fill((electron_x_Sigma - electron_x_truth) / electron_x_truth);
            h_RelRes_x_binned_EM->Fill(electron_x_truth, (electron_x_EM - electron_x_truth)/electron_x_truth);
            h_RelRes_x_binned_DA->Fill(electron_x_truth, (electron_x_DA - electron_x_truth)/electron_x_truth);
            h_RelRes_x_binned_ESigma->Fill(electron_x_truth, (electron_x_ESigma - electron_x_truth)/electron_x_truth);
            h_RelRes_x_binned_Sigma->Fill(electron_x_truth, (electron_x_Sigma - electron_x_truth)/electron_x_truth);
        }
        
        if(electron_y_truth != 0.0f){
            h_RelRes_y_EM->Fill((electron_y_EM - electron_y_truth) / electron_y_truth);
            h_RelRes_y_DA->Fill((electron_y_DA - electron_y_truth) / electron_y_truth);
            h_RelRes_y_ESigma->Fill((electron_y_ESigma - electron_y_truth) / electron_y_truth);
            h_RelRes_y_Sigma->Fill((electron_y_Sigma - electron_y_truth) / electron_y_truth);
            h_RelRes_y_binned_EM->Fill(electron_y_truth, (electron_y_EM - electron_y_truth)/electron_y_truth);
            h_RelRes_y_binned_DA->Fill(electron_y_truth, (electron_y_DA - electron_y_truth)/electron_y_truth);
            h_RelRes_y_binned_ESigma->Fill(electron_y_truth, (electron_y_ESigma - electron_y_truth)/electron_y_truth);
            h_RelRes_y_binned_Sigma->Fill(electron_y_truth, (electron_y_Sigma - electron_y_truth)/electron_y_truth);
        }
        
        // Correlations
        h_Corr_x_EM->Fill(electron_x_truth, electron_x_EM);
        h_Corr_x_DA->Fill(electron_x_truth, electron_x_DA);
        h_Corr_x_ESigma->Fill(electron_x_truth, electron_x_ESigma);
        h_Corr_x_Sigma->Fill(electron_x_truth, electron_x_Sigma);
        
        h_Corr_y_EM->Fill(electron_y_truth, electron_y_EM);
        h_Corr_y_DA->Fill(electron_y_truth, electron_y_DA);
        h_Corr_y_ESigma->Fill(electron_y_truth, electron_y_ESigma);
        h_Corr_y_Sigma->Fill(electron_y_truth, electron_y_Sigma);
        
        // Q2 histograms
        h_Q2_EM->Fill(electron_Q2_EM);
        h_Q2_DA->Fill(electron_Q2_DA);
        h_Q2_ESigma->Fill(electron_Q2_ESigma);
        h_Q2_Sigma->Fill(electron_Q2_Sigma);
        h_Q2_truth->Fill(electron_Q2_truth);
        
        h_RelRes_Q2_EM->Fill((electron_Q2_EM - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_DA->Fill((electron_Q2_DA - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_ESigma->Fill((electron_Q2_ESigma - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_Sigma->Fill((electron_Q2_Sigma - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_binned_EM->Fill(electron_Q2_truth, (electron_Q2_EM - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_binned_DA->Fill(electron_Q2_truth, (electron_Q2_DA - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_binned_ESigma->Fill(electron_Q2_truth, (electron_Q2_ESigma - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_binned_Sigma->Fill(electron_Q2_truth, (electron_Q2_Sigma - electron_Q2_truth)/electron_Q2_truth);
        
        h_Corr_Q2_EM->Fill(electron_Q2_truth, electron_Q2_EM);
        h_Corr_Q2_DA->Fill(electron_Q2_truth, electron_Q2_DA);
        h_Corr_Q2_ESigma->Fill(electron_Q2_truth, electron_Q2_ESigma);
        h_Corr_Q2_Sigma->Fill(electron_Q2_truth, electron_Q2_Sigma);
        
        // Process B0 and RP protons (if kinematics available)
        if(y_electron_array.GetSize() > 0 && x_electron_array.GetSize() > 0 && 
           q2_electron_array.GetSize() > 0 && w_electron_array.GetSize() > 0 &&
           electron_scat_index.GetSize() > 0 && truth_protons.size() > 0){
            
            double y = y_electron_array[0];
            double xBj = x_electron_array[0];
            double Q2 = q2_electron_array[0];
            double W = w_electron_array[0];
            double W2 = W * W;
            
            double xL_truth = truth_protons[0].P() / beams.p_beam.P();
            double t_truth = CalcT(beams.p_beam, truth_protons[0]);
            
            // Calculate x_P for truth
            double xP_truth = (MX2_truth + Q2 - TMath::Abs(t_truth)) / (W2 + Q2 - Mp2);
            h_xP_MC->Fill(xP_truth);
            
            // Process B0 protons
            for(unsigned int j = 0; j < tsassoc_rec_id.GetSize(); j++){
                auto mc_idx = tsassoc_sim_id[j];
                
                if(mc_idx >= (unsigned)mc_pdg_array.GetSize() || mc_genStatus_array[mc_idx] != 1 || 
                   mc_pdg_array[mc_idx] != 2212)
                    continue;
                    
                P3MVector p_reco(tsre_px_array[j], tsre_py_array[j], tsre_pz_array[j], mc_mass_array[mc_idx]);
                undoAfterburn(p_reco);
                
                // B0 angular acceptance
                if(p_reco.Theta() <= 0.0055 || p_reco.Theta() >= 0.02) continue;
                
                h_theta_B0->Fill(p_reco.Theta() * 1000.0);
                
                // Get truth for correlation
                P3MVector p_truth(mc_px_array[mc_idx], mc_py_array[mc_idx], mc_pz_array[mc_idx], 
                                 mc_mass_array[mc_idx]);
                undoAfterburn(p_truth);
                
                double t_reco = CalcT(beams.p_beam, p_reco);
                double t_truth_B0 = CalcT(beams.p_beam, p_truth);
                h_t_B0->Fill(TMath::Abs(t_reco));
                h_t_corr_B0->Fill(TMath::Abs(t_truth_B0), TMath::Abs(t_reco));
                
                if(TMath::Abs(t_truth_B0) > 1e-6){
                    double t_rel_res = (TMath::Abs(t_reco) - TMath::Abs(t_truth_B0)) / TMath::Abs(t_truth_B0);
                    h_t_res_B0->Fill(t_rel_res);
                }
                
                double xL_reco = p_reco.P() / beams.p_beam.P();
                double xL_truth_B0 = p_truth.P() / beams.p_beam.P();
                h_xL_B0->Fill(xL_reco);
                h_xL_corr_B0->Fill(xL_truth_B0, xL_reco);
                
                double MX2_reco = s * y * (1.0 - xL_reco - xBj);
                h_MX2_B0->Fill(MX2_reco);
                h_MX2_corr_B0->Fill(MX2_truth, MX2_reco);
                
                double xP_reco = (MX2_reco + Q2 - TMath::Abs(t_reco)) / (W2 + Q2 - Mp2);
                h_xP_B0->Fill(xP_reco);
                h_xP_corr_B0->Fill(xP_truth, xP_reco);
                
                n_b0_matches++;
            }
            
            // Process RP protons
            for(int j = 0; j < rp_px_array.GetSize(); j++){
                if(rp_pdg_array[j] != 2212) continue;
                
                P3MVector p_rp(rp_px_array[j], rp_py_array[j], rp_pz_array[j], rp_mass_array[j]);
                
                // RP angular acceptance
                if(p_rp.Theta() >= 0.005) continue;
                
                h_theta_RP->Fill(p_rp.Theta() * 1000.0);
                
                // Match to truth by angular proximity
                double min_dtheta = 1e9;
                int best_match = -1;
                
                for(size_t k = 0; k < truth_protons.size(); k++){
                    if(truth_protons[k].Theta() >= 0.006) continue;
                    
                    double dtheta = TMath::Abs(p_rp.Theta() - truth_protons[k].Theta());
                    if(dtheta < min_dtheta){
                        min_dtheta = dtheta;
                        best_match = k;
                    }
                }
                
                if(best_match >= 0 && min_dtheta < 0.001){
                    double t_reco = CalcT(beams.p_beam, p_rp);
                    double t_truth_RP = CalcT(beams.p_beam, truth_protons[best_match]);
                    h_t_RP->Fill(TMath::Abs(t_reco));
                    h_t_corr_RP->Fill(TMath::Abs(t_truth_RP), TMath::Abs(t_reco));
                    
                    if(TMath::Abs(t_truth_RP) > 1e-6){
                        double t_rel_res = (TMath::Abs(t_reco) - TMath::Abs(t_truth_RP)) / TMath::Abs(t_truth_RP);
                        h_t_res_RP->Fill(t_rel_res);
                    }
                    
                    double xL_reco = p_rp.P() / beams.p_beam.P();
                    double xL_truth_RP = truth_protons[best_match].P() / beams.p_beam.P();
                    h_xL_RP->Fill(xL_reco);
                    h_xL_corr_RP->Fill(xL_truth_RP, xL_reco);
                    
                    double MX2_reco = s * y * (1.0 - xL_reco - xBj);
                    h_MX2_RP->Fill(MX2_reco);
                    h_MX2_corr_RP->Fill(MX2_truth, MX2_reco);
                    
                    double xP_reco = (MX2_reco + Q2 - TMath::Abs(t_reco)) / (W2 + Q2 - Mp2);
                    h_xP_RP->Fill(xP_reco);
                    h_xP_corr_RP->Fill(xP_truth, xP_reco);
                    
                    n_rp_matches++;
                }
            }
        }
    }
    
    cout << "\n\nProcessing complete!" << endl;