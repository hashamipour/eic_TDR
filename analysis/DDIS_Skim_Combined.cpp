#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "TFile.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "Math/Vector4D.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/GenVector/Boost.h"
#include "Math/VectorUtil.h"

#include "Utility.hpp"
#include "RecoMethods.hpp"

using std::cout; using std::endl;
using ROOT::Math::VectorUtil::boost;
using ROOT::Math::RotationX;
using ROOT::Math::RotationY;
using P3EVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;
using MomVector = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>, ROOT::Math::DefaultCoordinateSystemTag>;

// Global afterburner correction parameters
Float_t fXAngle{-0.025};
RotationX rotAboutX;
RotationY rotAboutY;
MomVector vBoostToCoM;
MomVector vBoostToHoF;

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

/////////////////////////////////////////////////////////////////////////////////////////////
// NEW FUNCTION: Calculate sum(E - pz) for hadronic final state (MC Truth)
/////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Calculate sum(E - pz) for all MC final state particles excluding scattered electron and proton
 * 
 * This variable is used to suppress Initial State Radiation (ISR) background.
 * For a standard DIS event, sum(E - pz) should be close to 2*E_electron.
 * 
 * NOTE: Calculated in LAB FRAME (no afterburner corrections applied) to preserve the 
 * expected peak at 2*E_electron. For 10x100 GeV beams, expect peak around 20 GeV.
 * 
 * @param mc_px MC particle px array
 * @param mc_py MC particle py array
 * @param mc_pz MC particle pz array
 * @param mc_mass MC particle mass array
 * @param mc_genStatus MC particle generator status array
 * @param mc_pdg MC particle PDG code array
 * @return Sum of (E - pz) for hadronic final state in lab frame
 */
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
        // Only final state particles (generatorStatus == 1)
        if(mc_genStatus[i] != 1) continue;
        
        // Exclude scattered electron (PDG == 11) and scattered proton (PDG == 2212)
        if(mc_pdg[i] == 11 || mc_pdg[i] == 2212) continue;
        
        // Create 4-vector - NO afterburner correction for E-pz calculation
        // We want this in the lab frame to get the expected peak at 2*E_electron
        double px = mc_px[i];
        double py = mc_py[i];
        double pz = mc_pz[i];
        double m = mc_mass[i];
        double E = TMath::Sqrt(px*px + py*py + pz*pz + m*m);
        
        sumEPz += (E - pz);
    }
    
    return sumEPz;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// NEW FUNCTION: Calculate sum(E - pz) for hadronic final state (Reconstructed)
/////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Calculate sum(E - pz) for reconstructed particles excluding scattered electron
 * 
 * NOTE: Calculated in LAB FRAME (no afterburner corrections applied) to preserve the 
 * expected peak at 2*E_electron. For 10x100 GeV beams, expect peak around 20 GeV.
 * 
 * @param re_px Reconstructed particle px array
 * @param re_py Reconstructed particle py array
 * @param re_pz Reconstructed particle pz array
 * @param re_energy Reconstructed particle energy array
 * @param electron_scat_index Index of the scattered electron in ReconstructedParticles
 * @return Sum of (E - pz) for hadronic final state in lab frame
 */
double CalculateSumEPz_Reco(
    TTreeReaderArray<float>& re_px,
    TTreeReaderArray<float>& re_py,
    TTreeReaderArray<float>& re_pz,
    TTreeReaderArray<float>& re_energy,
    TTreeReaderArray<int>& electron_scat_index
) {
    double sumEPz = 0.0;
    
    // Get the index of the scattered electron (if available)
    int scat_e_idx = -1;
    if(electron_scat_index.GetSize() > 0) {
        scat_e_idx = electron_scat_index[0];
    }
    
    for(int i = 0; i < re_px.GetSize(); i++){
        // Skip the scattered electron
        if(i == scat_e_idx) continue;
        
        // Use 4-vector components directly - NO afterburner correction
        // We want this in the lab frame to get the expected peak at 2*E_electron
        double px = re_px[i];
        double py = re_py[i];
        double pz = re_pz[i];
        double E = re_energy[i];
        
        sumEPz += (E - pz);
    }
    
    return sumEPz;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// NEW FUNCTION: Calculate M_X using 4-vector sum definition (MC Truth)
/////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Calculate M_X using the 4-vector sum: M_X^2 = (sum E)^2 - (sum p)^2
 * Returns M_X (not M_X^2)
 */
double CalculateMX_Formula_Truth(
    TTreeReaderArray<double>& mc_px,
    TTreeReaderArray<double>& mc_py,
    TTreeReaderArray<double>& mc_pz,
    TTreeReaderArray<double>& mc_mass,
    TTreeReaderArray<int>& mc_genStatus,
    TTreeReaderArray<int>& mc_pdg
) {
    P3MVector X_system = CalculateTruthXSystem(mc_px, mc_py, mc_pz, mc_mass, mc_genStatus, mc_pdg);
    double MX2 = X_system.M2();
    return (MX2 > 0) ? TMath::Sqrt(MX2) : 0.0;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// NEW FUNCTION: Calculate M_X using 4-vector sum (Reconstructed)
/////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Calculate M_X using the 4-vector sum for reconstructed particles
 * Returns M_X (not M_X^2)
 */
double CalculateMX_Formula_Reco(
    TTreeReaderArray<float>& re_px,
    TTreeReaderArray<float>& re_py,
    TTreeReaderArray<float>& re_pz,
    TTreeReaderArray<float>& re_energy,
    TTreeReaderArray<int>& electron_scat_index
) {
    // Get the index of the scattered electron
    int scat_e_idx = -1;
    if(electron_scat_index.GetSize() > 0) {
        scat_e_idx = electron_scat_index[0];
    }
    
    P3MVector X_system(0,0,0,0);
    
    for(int i = 0; i < re_px.GetSize(); i++){
        if(i == scat_e_idx) continue;
        
        double px = re_px[i];
        double py = re_py[i];
        double pz = re_pz[i];
        double E = re_energy[i];
        
        P3MVector p(px, py, pz, TMath::Sqrt(E*E - px*px - py*py - pz*pz));
        undoAfterburn(p);
        X_system += p;
    }
    
    double MX2 = X_system.M2();
    return (MX2 > 0) ? TMath::Sqrt(MX2) : 0.0;
}

void analyzeProtonsMandelstamT(TString fileList){
    cout << "----------------------------" << endl;
    cout << "  Proton Mandelstam t Analysis" << endl;
    cout << "  with E-pz and M_X" << endl;
    cout << "----------------------------" << endl;
    cout << "\nInput filelist: " << fileList << endl;

    BeamInfo beams;
    
    // Setup input chain
    TChain* events = new TChain("events");
    Int_t nFiles = 0;
    
    std::ifstream fileListStream(fileList.Data());
    std::string fileName;
    while(getline(fileListStream, fileName)){
        events->Add((TString)fileName);
        nFiles++;
    }
    
    cout << "\nNo. of files: " << nFiles << "; no. of events: " << events->GetEntries() << endl;
    
    // Setup binning
    std::vector<Double_t> t_bins = GetLinBins(0.5, 1.5, 10);
    
    // Truth histograms
    // MC Truth: Combined (log up to 0.5, then linear to 1.6)
    std::vector<Double_t> t_bins_log = GetLogBins(1e-4, 0.5, 30);
    std::vector<Double_t> t_bins_lin = GetLinBins(0.5, 1.6, 11);
    std::vector<Double_t> t_bins_MC;
    t_bins_MC.insert(t_bins_MC.end(), t_bins_log.begin(), t_bins_log.end());
    // Skip first element of linear bins to avoid duplicate at 0.5
    t_bins_MC.insert(t_bins_MC.end(), t_bins_lin.begin()+1, t_bins_lin.end());

    cout << "bin edges for t MC: ";
    for(const auto& edge : t_bins_MC) cout << edge << " ";
    cout << endl;

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
    
    //---------------------------------------------------------
    // NEW HISTOGRAMS: Sum(E - pz) and M_X
    //---------------------------------------------------------
    
    // Sum(E - pz) histograms
    TH1D* h_SumEPz_truth = new TH1D("h_SumEPz_truth", 
        "Truth #Sigma(E-p_{z}) Hadronic;#Sigma(E-p_{z}) [GeV];Counts", 
        100, 0.0, 100.0);
    TH1D* h_SumEPz_reco = new TH1D("h_SumEPz_reco", 
        "Reco #Sigma(E-p_{z}) Hadronic;#Sigma(E-p_{z}) [GeV];Counts", 
        100, 0.0, 100.0);
    
    // M_X histograms (using 4-vector formula: M_X^2 = (sum E)^2 - (sum p)^2)
    TH1D* h_Mx_truth = new TH1D("h_Mx_truth",
        "Truth M_{X} (4-vector);M_{X} [GeV];Counts",
        50, 1.0, 50.0);
    TH1D* h_Mx_reco = new TH1D("h_Mx_reco",
        "Reco M_{X} (4-vector);M_{X} [GeV];Counts",
        50, 1.0, 50.0);
    
    cout << "\n=== NEW HISTOGRAMS CREATED ===" << endl;
    cout << "h_SumEPz_truth: Sum(E-pz) for MC truth hadronic final state" << endl;
    cout << "h_SumEPz_reco:  Sum(E-pz) for reconstructed hadronic final state" << endl;
    cout << "h_Mx_truth:     M_X using 4-vector sum (truth)" << endl;
    cout << "h_Mx_reco:      M_X using 4-vector sum (reco)" << endl;
    
    // Setup tree reader
    TTreeReader tree_reader(events);
    
    TTreeReaderArray<double> mc_px_array(tree_reader, "MCParticles.momentum.x");
    TTreeReaderArray<double> mc_py_array(tree_reader, "MCParticles.momentum.y");
    TTreeReaderArray<double> mc_pz_array(tree_reader, "MCParticles.momentum.z");
    TTreeReaderArray<double> mc_mass_array(tree_reader, "MCParticles.mass");
    TTreeReaderArray<int> mc_genStatus_array(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<int> mc_pdg_array(tree_reader, "MCParticles.PDG");
    
    TTreeReaderArray<unsigned int> tsassoc_rec_id(tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.recID");
    TTreeReaderArray<unsigned int> tsassoc_sim_id(tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.simID");
    
    TTreeReaderArray<float> tsre_px_array(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x");
    TTreeReaderArray<float> tsre_py_array(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y");
    TTreeReaderArray<float> tsre_pz_array(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z");
    
    TTreeReaderArray<float> rp_px_array(tree_reader, "ForwardRomanPotRecParticles.momentum.x");
    TTreeReaderArray<float> rp_py_array(tree_reader, "ForwardRomanPotRecParticles.momentum.y");
    TTreeReaderArray<float> rp_pz_array(tree_reader, "ForwardRomanPotRecParticles.momentum.z");
    TTreeReaderArray<float> rp_mass_array(tree_reader, "ForwardRomanPotRecParticles.mass");
    TTreeReaderArray<int> rp_pdg_array(tree_reader, "ForwardRomanPotRecParticles.PDG");
    
    // Kinematic variables from EICRecon
    TTreeReaderArray<float> y_electron_array(tree_reader, "InclusiveKinematicsElectron.y");
    TTreeReaderArray<float> x_electron_array(tree_reader, "InclusiveKinematicsElectron.x");
    TTreeReaderArray<float> q2_electron_array(tree_reader, "InclusiveKinematicsElectron.Q2");
    TTreeReaderArray<float> w_electron_array(tree_reader, "InclusiveKinematicsElectron.W");
    TTreeReaderArray<int> electron_scat_index(tree_reader, "_InclusiveKinematicsElectron_scat.index");
    
    // NEW: ReconstructedParticles arrays for Sum(E-pz) and M_X calculations
    TTreeReaderArray<float> re_px_array(tree_reader, "ReconstructedParticles.momentum.x");
    TTreeReaderArray<float> re_py_array(tree_reader, "ReconstructedParticles.momentum.y");
    TTreeReaderArray<float> re_pz_array(tree_reader, "ReconstructedParticles.momentum.z");
    TTreeReaderArray<float> re_energy_array(tree_reader, "ReconstructedParticles.energy");
    
    // Find beam particles
    cout << "Finding beam particles..." << endl;
    
    P3MVector beame4_acc(0,0,0,0), beamp4_acc(0,0,0,0);
    
    while(tree_reader.Next()){
        for(int i = 0; i < mc_px_array.GetSize(); i++){
            if(mc_genStatus_array[i] != 4) continue;
            
            if(mc_pdg_array[i] == 2212){
                P3MVector p(mc_px_array[i], mc_py_array[i], mc_pz_array[i], beams.fMass_proton);
                beamp4_acc += p;
            }
            else if(mc_pdg_array[i] == 11){
                P3MVector p(mc_px_array[i], mc_py_array[i], mc_pz_array[i], beams.fMass_electron);
                beame4_acc += p;
            }
        }
    }
    
    auto nEntries = std::max<Long64_t>(1, events->GetEntries());
    beams.e_beam.SetCoordinates(
        beame4_acc.X()/nEntries, 
        beame4_acc.Y()/nEntries, 
        beame4_acc.Z()/nEntries, 
        beams.fMass_electron
    );
    beams.p_beam.SetCoordinates(
        beamp4_acc.X()/nEntries, 
        beamp4_acc.Y()/nEntries, 
        beamp4_acc.Z()/nEntries, 
        beams.fMass_proton
    );
    
    cout << "Found beam energies " << beams.e_beam.E() << "x" << beams.p_beam.E() << " GeV" << endl;
    
    undoAfterburnAndCalc(beams.p_beam, beams.e_beam);
    
    tree_reader.Restart();
    
    // Calculate s (center-of-mass energy squared)
    P3MVector total_beam = beams.e_beam + beams.p_beam;
    double s = total_beam.M2();
    double Mp2 = beams.fMass_proton * beams.fMass_proton;
    cout << "Center-of-mass energy squared s = " << s << " GeV^2" << endl;
    cout << "Proton mass squared Mp^2 = " << Mp2 << " GeV^2" << endl;
    
    // Event loop
    int n_truth_protons = 0, n_b0_matches = 0, n_rp_matches = 0;
    
    while(tree_reader.Next()){
        // Process truth protons
        auto truth_protons = ProcessTruthProtons(
            mc_px_array, mc_py_array, mc_pz_array, mc_mass_array,
            mc_genStatus_array, mc_pdg_array, beams, h_theta_MC, h_t_MC, h_xL_MC
        );
        n_truth_protons += truth_protons.size();
        
        //---------------------------------------------------------
        // CALCULATE NEW VARIABLES: Sum(E-pz) and M_X
        //---------------------------------------------------------
        
        // Calculate sum(E - pz) for truth
        double sumEPz_truth = CalculateSumEPz_Truth(
            mc_px_array, mc_py_array, mc_pz_array, mc_mass_array,
            mc_genStatus_array, mc_pdg_array
        );
        h_SumEPz_truth->Fill(sumEPz_truth);
        
        // Calculate M_X using 4-vector formula (truth)
        double Mx_truth = CalculateMX_Formula_Truth(
            mc_px_array, mc_py_array, mc_pz_array, mc_mass_array,
            mc_genStatus_array, mc_pdg_array
        );
        h_Mx_truth->Fill(Mx_truth);
        
        // Calculate sum(E - pz) and M_X for reco (if ReconstructedParticles available)
        if(re_px_array.GetSize() > 0 && re_energy_array.GetSize() > 0) {
            double sumEPz_reco = CalculateSumEPz_Reco(
                re_px_array, re_py_array, re_pz_array, re_energy_array, electron_scat_index
            );
            h_SumEPz_reco->Fill(sumEPz_reco);
            
            double Mx_reco = CalculateMX_Formula_Reco(
                re_px_array, re_py_array, re_pz_array, re_energy_array, electron_scat_index
            );
            h_Mx_reco->Fill(Mx_reco);
        }
        
        //---------------------------------------------------------
        // EXISTING CODE: Calculate truth M_X^2 (kinematic formula)
        //---------------------------------------------------------
        
        // Calculate truth M_X^2
        P3MVector X_truth = CalculateTruthXSystem(
            mc_px_array, mc_py_array, mc_pz_array, mc_mass_array,
            mc_genStatus_array, mc_pdg_array
        );
        double MX2_truth = CalcMX2(X_truth);
        h_MX2_MC->Fill(MX2_truth);
        
        // Get kinematic variables from EICRecon
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
            for(unsigned int i = 0; i < tsassoc_rec_id.GetSize(); i++){
                auto mc_idx = tsassoc_sim_id[i];
                
                if(mc_idx >= (unsigned)mc_pdg_array.GetSize() || mc_genStatus_array[mc_idx] != 1 || 
                   mc_pdg_array[mc_idx] != 2212)
                    continue;
                    
                P3MVector p_reco(tsre_px_array[i], tsre_py_array[i], tsre_pz_array[i], mc_mass_array[mc_idx]);
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
                
                // Fill resolution histogram
                if(TMath::Abs(t_truth_B0) > 1e-6) {
                    double t_rel_res = (TMath::Abs(t_reco) - TMath::Abs(t_truth_B0)) / TMath::Abs(t_truth_B0);
                    h_t_res_B0->Fill(t_rel_res);
                }
                
                // Calculate x_L and M_X^2
                double xL_reco = p_reco.P() / beams.p_beam.P();
                double xL_truth_B0 = p_truth.P() / beams.p_beam.P();
                h_xL_B0->Fill(xL_reco);
                h_xL_corr_B0->Fill(xL_truth_B0, xL_reco);
                
                double MX2_reco = s * y * (1.0 - xL_reco - xBj);
                h_MX2_B0->Fill(MX2_reco);
                h_MX2_corr_B0->Fill(MX2_truth, MX2_reco);
                
                // Calculate x_P
                double xP_reco = (MX2_reco + Q2 - TMath::Abs(t_reco)) / (W2 + Q2 - Mp2);
                h_xP_B0->Fill(xP_reco);
                h_xP_corr_B0->Fill(xP_truth, xP_reco);
                
                n_b0_matches++;
            }
            
            // Process RP protons
            for(int i = 0; i < rp_px_array.GetSize(); i++){
                if(rp_pdg_array[i] != 2212) continue;
                
                P3MVector p_rp(rp_px_array[i], rp_py_array[i], rp_pz_array[i], rp_mass_array[i]);
                
                // RP angular acceptance
                if(p_rp.Theta() >= 0.005) continue;
                
                h_theta_RP->Fill(p_rp.Theta() * 1000.0);
                
                // Match to truth by angular proximity
                double min_dtheta = 1e9;
                int best_match = -1;
                
                for(size_t j = 0; j < truth_protons.size(); j++){
                    if(truth_protons[j].Theta() >= 0.006) continue;
                    
                    double dtheta = TMath::Abs(p_rp.Theta() - truth_protons[j].Theta());
                    if(dtheta < min_dtheta){
                        min_dtheta = dtheta;
                        best_match = j;
                    }
                }
                
                if(best_match >= 0 && min_dtheta < 0.001){
                    double t_reco = CalcT(beams.p_beam, p_rp);
                    double t_truth_RP = CalcT(beams.p_beam, truth_protons[best_match]);
                    h_t_RP->Fill(TMath::Abs(t_reco));
                    h_t_corr_RP->Fill(TMath::Abs(t_truth_RP), TMath::Abs(t_reco));
                    
                    // Fill resolution histogram
                    if(TMath::Abs(t_truth_RP) > 1e-6) {
                        double t_rel_res = (TMath::Abs(t_reco) - TMath::Abs(t_truth_RP)) / TMath::Abs(t_truth_RP);
                        h_t_res_RP->Fill(t_rel_res);
                    }
                    
                    // Calculate x_L and M_X^2
                    double xL_reco = p_rp.P() / beams.p_beam.P();
                    double xL_truth_RP = truth_protons[best_match].P() / beams.p_beam.P();
                    h_xL_RP->Fill(xL_reco);
                    h_xL_corr_RP->Fill(xL_truth_RP, xL_reco);
                    
                    double MX2_reco = s * y * (1.0 - xL_reco - xBj);
                    h_MX2_RP->Fill(MX2_reco);
                    h_MX2_corr_RP->Fill(MX2_truth, MX2_reco);
                    
                    // Calculate x_P
                    double xP_reco = (MX2_reco + Q2 - TMath::Abs(t_reco)) / (W2 + Q2 - Mp2);
                    h_xP_RP->Fill(xP_reco);
                    h_xP_corr_RP->Fill(xP_truth, xP_reco);
                    
                    n_rp_matches++;
                }
            }
        }
    }
    
    // Print results
    cout << "\n=== ANALYSIS RESULTS ===" << endl;
    cout << "Truth protons found: " << n_truth_protons << endl;
    cout << "B0 matched protons: " << n_b0_matches << endl;
    cout << "RP matched protons: " << n_rp_matches << endl;
    if(n_truth_protons > 0){
        cout << "B0 matching efficiency: " << (double)n_b0_matches/n_truth_protons*100 << "%" << endl;
        cout << "RP matching efficiency: " << (double)n_rp_matches/n_truth_protons*100 << "%" << endl;
    }
    
    // Save to file
    gStyle->SetOptStat(0);
    
    TFile* outfile = new TFile("proton_mandelstam_analysis.root", "RECREATE");
    
    // Write truth histograms
    h_t_MC->Write();
    h_theta_MC->Write();
    h_xL_MC->Write();
    h_MX2_MC->Write();
    h_xP_MC->Write();
    
    // Write B0 histograms
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
    
    // Write RP histograms
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
    
    //---------------------------------------------------------
    // WRITE NEW HISTOGRAMS
    //---------------------------------------------------------
    h_SumEPz_truth->Write();
    h_SumEPz_reco->Write();
    h_Mx_truth->Write();
    h_Mx_reco->Write();
    
    cout << "\n=== NEW HISTOGRAMS WRITTEN ===" << endl;
    
    outfile->Close();
    
    cout << "\nAnalysis complete! Output saved to proton_mandelstam_analysis.root" << endl;
    cout << "M_X^2 calculated using formula: M_X^2 = s*y*(1 - x_L - x_Bj)" << endl;
    cout << "x_P calculated using formula: x_P = (M_X^2 + Q^2 - t) / (W^2 + Q^2 - M_p^2)" << endl;
    cout << "\nNEW: Sum(E-pz) and M_X (4-vector) histograms added" << endl;
}

int main(int argc, char** argv){
    if(argc < 2){
        std::cerr << "Usage: skim_t <filelist.txt>\n";
        return 1;
    }
    
    analyzeProtonsMandelstamT(argv[1]);
    return 0;
}