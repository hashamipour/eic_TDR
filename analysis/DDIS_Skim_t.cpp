#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "TVector3.h"
#include "Math/Vector4D.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/GenVector/Boost.h"
#include "Math/VectorUtil.h"
#include "TLine.h"

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
            h_xL->Fill(CalcXL(p));
        }
    }
    
    return truth_protons;
}

// Process B0 (truth-seeded) protons
void ProcessB0Protons(
    TTreeReaderArray<unsigned int>& tsassoc_rec_id, TTreeReaderArray<unsigned int>& tsassoc_sim_id,
    TTreeReaderArray<float>& tsre_px, TTreeReaderArray<float>& tsre_py, TTreeReaderArray<float>& tsre_pz,
    TTreeReaderArray<double>& mc_px, TTreeReaderArray<double>& mc_py, TTreeReaderArray<double>& mc_pz,
    TTreeReaderArray<double>& mc_mass, TTreeReaderArray<int>& mc_gen_status, TTreeReaderArray<int>& mc_pdg,
    const BeamInfo& beams, MethodHistograms& hist, int& n_matches
) {
    for(unsigned int i = 0; i < tsassoc_rec_id.GetSize(); i++){
        auto mc_idx = tsassoc_sim_id[i];
        
        if(mc_idx >= (unsigned)mc_pdg.GetSize() || mc_gen_status[mc_idx] != 1 || mc_pdg[mc_idx] != 2212)
            continue;
            
        P3MVector p_reco(tsre_px[i], tsre_py[i], tsre_pz[i], mc_mass[mc_idx]);
        undoAfterburn(p_reco);
        
        // B0 angular acceptance
        if(p_reco.Theta() <= 0.0055 || p_reco.Theta() >= 0.02) continue;
        
        hist.h_theta->Fill(p_reco.Theta() * 1000.0);
        
        // Get truth for correlation
        P3MVector p_truth(mc_px[mc_idx], mc_py[mc_idx], mc_pz[mc_idx], mc_mass[mc_idx]);
        undoAfterburn(p_truth);
        
        double t_reco = CalcT(beams.p_beam, p_reco);
        double t_truth = CalcT(beams.p_beam, p_truth);
        
        hist.FillCorrelation(t_truth, t_reco);
        
        double xL_reco = CalcXL(p_reco);
        double xL_truth = CalcXL(p_truth);
        hist.h_xL->Fill(xL_reco);
        hist.h_xL_corr->Fill(xL_truth, xL_reco);
        
        n_matches++;
    }
}

// Process Roman Pot protons
void ProcessRPProtons(
    TTreeReaderArray<float>& rp_px, TTreeReaderArray<float>& rp_py, TTreeReaderArray<float>& rp_pz,
    TTreeReaderArray<float>& rp_mass, TTreeReaderArray<int>& rp_pdg,
    const std::vector<P3MVector>& truth_protons, const BeamInfo& beams,
    MethodHistograms& hist, int& n_matches
) {
    for(int i = 0; i < rp_px.GetSize(); i++){
        if(rp_pdg[i] != 2212) continue;
        
        P3MVector p_rp(rp_px[i], rp_py[i], rp_pz[i], rp_mass[i]);
        
        // RP angular acceptance
        if(p_rp.Theta() >= 0.005) continue;
        
        hist.h_theta->Fill(p_rp.Theta() * 1000.0);
        
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
            double t_truth = CalcT(beams.p_beam, truth_protons[best_match]);
            
            hist.FillCorrelation(t_truth, t_reco);
            
            double xL_reco = CalcXL(p_rp);
            double xL_truth = CalcXL(truth_protons[best_match]);
            hist.h_xL->Fill(xL_reco);
            hist.h_xL_corr->Fill(xL_truth, xL_reco);
            
            n_matches++;
        }
    }
}

void analyzeProtonsMandelstamT(TString fileList){
    cout << "----------------------------" << endl;
    cout << "  Proton Mandelstam t Analysis" << endl;
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
    std::vector<Double_t> t_bins = GetLinBins(0.0, 1.6, 10);
    cout << "bin edges for t: ";
    for(const auto& edge : t_bins) cout << edge << " ";
    cout << endl;
    
    // Create histograms for each method
    TH1D* h_t_MC = new TH1D("t_MC", "Truth Mandelstam t;|t| [GeV^{2}];Counts", 
                            t_bins.size()-1, t_bins.data());
    TH1D* h_theta_MC = new TH1D("theta_MC", "MC Proton Scattering Angle;#theta [mrad];Counts", 
                                100, 0.0, 25.0);
    TH1D* h_xL_MC = new TH1D("xL_MC", "Truth x_{L};x_{L};Counts", 20, 0.92, 1.02);   
    
    // Q² comparison histograms
    TH1D* h_Q2_EICRecon = new TH1D("Q2_EICRecon", "EICRecon Q^{2};Q^{2} [GeV^{2}];Counts", 100, 0, 20);
    TH1D* h_Q2_calc = new TH1D("Q2_calc", "Calculated Q^{2} from e^{-};Q^{2} [GeV^{2}];Counts", 100, 0, 20);
    TH2D* h_Q2_corr = new TH2D("Q2_corr", "Q^{2} Correlation;EICRecon Q^{2} [GeV^{2}];Calc Q^{2} [GeV^{2}]", 
                               100, 0, 20, 100, 0, 20);
    
    MethodHistograms hist_B0("B0", t_bins);
    MethodHistograms hist_RP("RP", t_bins);
    MethodHistograms hist_eX("eX", t_bins);
    
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

    // Particle Flow collection (includes charged + neutral)
    TTreeReaderArray<float> rpf_px_array(tree_reader, "ReconstructedParticles.momentum.x");
    TTreeReaderArray<float> rpf_py_array(tree_reader, "ReconstructedParticles.momentum.y");
    TTreeReaderArray<float> rpf_pz_array(tree_reader, "ReconstructedParticles.momentum.z");
    TTreeReaderArray<float> rpf_e_array(tree_reader, "ReconstructedParticles.energy");
    TTreeReaderArray<int>   rpf_pdg_array(tree_reader, "ReconstructedParticles.PDG");
    
    // Q² from EICRecon
    TTreeReaderArray<float> q2_electron_array(tree_reader, "InclusiveKinematicsElectron.Q2");

    TTreeReaderArray<int> electron_scat_index(tree_reader, "_InclusiveKinematicsElectron_scat.index");

    
    // Find beam particles
    cout << "Finding beam particles." << endl;
    
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
    
    cout << "[DEBUG] Found beam energies " << beams.e_beam.E() << "x" << beams.p_beam.E() << " GeV" << endl;
    
    undoAfterburnAndCalc(beams.p_beam, beams.e_beam);
    
    tree_reader.Restart();
    
    // Event loop
    int n_truth_protons = 0, n_b0_matches = 0, n_rp_matches = 0;
    
    while(tree_reader.Next()){
        // Process truth protons
        auto truth_protons = ProcessTruthProtons(
            mc_px_array, mc_py_array, mc_pz_array, mc_mass_array,
            mc_genStatus_array, mc_pdg_array, beams, h_theta_MC, h_t_MC, h_xL_MC
        );
        n_truth_protons += truth_protons.size();
        
        // Process B0 protons
        ProcessB0Protons(
            tsassoc_rec_id, tsassoc_sim_id,
            tsre_px_array, tsre_py_array, tsre_pz_array,
            mc_px_array, mc_py_array, mc_pz_array, mc_mass_array,
            mc_genStatus_array, mc_pdg_array,
            beams, hist_B0, n_b0_matches
        );
        
        // Process RP protons
        ProcessRPProtons(
            rp_px_array, rp_py_array, rp_pz_array, rp_mass_array, rp_pdg_array,
            truth_protons, beams, hist_RP, n_rp_matches
        );

        // Process eX method with Q² diagnostic
        auto elec_info = GetScatteredElectron(electron_scat_index, rpf_px_array, 
                                      rpf_py_array, rpf_pz_array, beams.fMass_electron);

        if(elec_info.found && truth_protons.size() > 0){
            P3MVector e_scattered = elec_info.p4;
            undoAfterburn(e_scattered);
            
            P3MVector q_gamma = beams.e_beam - e_scattered;
            double Q2_calc = -q_gamma.M2();
            
            // Fill Q² histograms...
            
            // Build X excluding electron and protons
            P3MVector X_system(0,0,0,0);
            for(int i = 0; i < rpf_px_array.GetSize(); i++){
                if(i == (int)elec_info.index) continue;  // Use stored index
                if(rpf_pdg_array[i] == 2212) continue;
                
                P3MVector p(rpf_px_array[i], rpf_py_array[i], rpf_pz_array[i], rpf_e_array[i]);
                undoAfterburn(p);
                X_system += p;
            }
            
            double t_eX = CalcT_eX(q_gamma, X_system);
            double t_truth = CalcT(beams.p_beam, truth_protons[0]);
            hist_eX.FillCorrelation(t_truth, t_eX);
            
            double xL_truth = CalcXL(truth_protons[0]);
            double xL_eX = CalcXL(beams.p_beam - q_gamma - X_system);
            hist_eX.h_xL->Fill(xL_eX);
            hist_eX.h_xL_corr->Fill(xL_truth, xL_eX);
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
    h_t_MC->Write();
    h_theta_MC->Write();
    h_xL_MC->Write();
    h_Q2_EICRecon->Write();
    h_Q2_calc->Write();
    h_Q2_corr->Write();
    hist_B0.Write();
    hist_RP.Write();
    hist_eX.Write();
    outfile->Close();
    
    cout << "\nAnalysis complete! Output saved to proton_mandelstam_analysis.root" << endl;
}

int main(int argc, char** argv){
    if(argc < 2){
        std::cerr << "Usage: skim_t <filelist.txt>\n";
        return 1;
    }
    
    analyzeProtonsMandelstamT(argv[1]);
    return 0;
}