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
#include "TVector2.h"

#include "Utility.hpp"
#include "RecoMethods.hpp"

using std::cout; using std::endl;
using ROOT::Math::VectorUtil::boost;
using ROOT::Math::RotationX;
using ROOT::Math::RotationY;
using P3EVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;
using MomVector = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>, ROOT::Math::DefaultCoordinateSystemTag>;

using std::cout; using std::endl;








using std::cout; using std::endl;
using P3MVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzMVector>;

// Dedicated function to calculate Delta R
double CalculateDeltaR(const P3MVector& p1, const P3MVector& p2) {
    double eta1 = p1.Eta();
    double phi1 = p1.Phi();
    double eta2 = p2.Eta();
    double phi2 = p2.Phi();
    
    double deta = eta1 - eta2;
    double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
    
    return TMath::Sqrt(deta*deta + dphi*dphi);
}


// Copy the undoAfterburn functions from original code
Float_t fXAngle{-0.025};
ROOT::Math::RotationX rotAboutX;
ROOT::Math::RotationY rotAboutY;
ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>, ROOT::Math::DefaultCoordinateSystemTag> vBoostToCoM;
ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>, ROOT::Math::DefaultCoordinateSystemTag> vBoostToHoF;

void undoAfterburnAndCalc(P3MVector& p, P3MVector& k){
    P3MVector p_beam(fXAngle*p.E(), 0., p.E(), p.M());
    P3MVector e_beam(0., 0., -k.E(), k.M());
    
    P3MVector CoM_boost = p_beam + e_beam;
    vBoostToCoM.SetXYZ(-CoM_boost.X()/CoM_boost.E(), -CoM_boost.Y()/CoM_boost.E(), -CoM_boost.Z()/CoM_boost.E());
    
    p_beam = ROOT::Math::VectorUtil::boost(p_beam, vBoostToCoM);
    e_beam = ROOT::Math::VectorUtil::boost(e_beam, vBoostToCoM);
    
    Float_t fRotY = -1.0*TMath::ATan2(p_beam.X(), p_beam.Z());
    Float_t fRotX = 1.0*TMath::ATan2(p_beam.Y(), p_beam.Z());
    
    rotAboutY = ROOT::Math::RotationY(fRotY);
    rotAboutX = ROOT::Math::RotationX(fRotX);
    
    p_beam = rotAboutY(p_beam);
    p_beam = rotAboutX(p_beam);
    e_beam = rotAboutY(e_beam);
    e_beam = rotAboutX(e_beam);
    
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector> HoF_boost(0., 0., CoM_boost.Z(), CoM_boost.E());
    vBoostToHoF.SetXYZ(HoF_boost.X()/HoF_boost.E(), HoF_boost.Y()/HoF_boost.E(), HoF_boost.Z()/HoF_boost.E());
    
    p_beam = ROOT::Math::VectorUtil::boost(p_beam, vBoostToHoF);
    e_beam = ROOT::Math::VectorUtil::boost(e_beam, vBoostToHoF);
    
    p.SetPxPyPzE(p_beam.X(), p_beam.Y(), p_beam.Z(), p_beam.E());
    k.SetPxPyPzE(e_beam.X(), e_beam.Y(), e_beam.Z(), e_beam.E());
}

void undoAfterburn(P3MVector& a){
    a = ROOT::Math::VectorUtil::boost(a, vBoostToCoM);
    a = rotAboutY(a);
    a = rotAboutX(a);
    a = ROOT::Math::VectorUtil::boost(a, vBoostToHoF);
}

void analyzeRPProtons(const std::string& fileList) {
    cout << "=== Enhanced RP Proton Analysis ===" << endl;
    
    // Setup input chain
    TChain* events = new TChain("events");
    
    std::ifstream fileListStream(fileList);
    std::string fileName;
    int nFiles = 0;
    while(getline(fileListStream, fileName)){
        events->Add(fileName.c_str());
        nFiles++;
    }
    
    cout << "Files: " << nFiles << ", Events: " << events->GetEntries() << endl;
    
    // Create histograms
    TH1D* h_deltaR_all = new TH1D("deltaR_all", "All RP-Truth #DeltaR matches;#DeltaR;Entries", 
                                  100, 0, 0.2);
    TH1D* h_deltaR_best = new TH1D("deltaR_best", "Best RP-Truth #DeltaR matches;#DeltaR;Entries", 
                                   100, 0, 0.05);
    TH1D* h_theta_mc_rp = new TH1D("theta_mc_rp", "MC Protons in RP acceptance;#theta [mrad];Entries", 
                                   100, 0, 6);
    TH1D* h_theta_reco_rp = new TH1D("theta_reco_rp", "Reconstructed RP Protons;#theta [mrad];Entries", 
                                     100, 0, 6);
    TH1D* h_eta_mc_rp = new TH1D("eta_mc_rp", "MC Protons in RP acceptance;#eta;Entries", 
                                 100, 4, 12);
    TH1D* h_eta_reco_rp = new TH1D("eta_reco_rp", "Reconstructed RP Protons;#eta;Entries", 
                                   100, 4, 12);
    TH2D* h_deltaR_vs_theta = new TH2D("deltaR_vs_theta", "#DeltaR vs #theta;MC #theta [mrad];#DeltaR", 
                                       50, 0, 6, 50, 0, 0.05);
    TH1D* h_rp_multiplicity = new TH1D("rp_multiplicity", "RP Proton Multiplicity;N_{RP protons};Events", 
                                       10, 0, 10);
    TH1D* h_mc_multiplicity = new TH1D("mc_multiplicity", "MC RP Proton Multiplicity;N_{MC protons in RP};Events", 
                                       10, 0, 10);
    
    // Setup tree reader
    TTreeReader tree_reader(events);
    TTreeReaderArray<double> mc_px_array(tree_reader, "MCParticles.momentum.x");
    TTreeReaderArray<double> mc_py_array(tree_reader, "MCParticles.momentum.y");
    TTreeReaderArray<double> mc_pz_array(tree_reader, "MCParticles.momentum.z");
    TTreeReaderArray<double> mc_mass_array(tree_reader, "MCParticles.mass");
    TTreeReaderArray<int> mc_genStatus_array(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<int> mc_pdg_array(tree_reader, "MCParticles.PDG");
    
    TTreeReaderArray<float> rp_px_array(tree_reader, "ForwardRomanPotRecParticles.momentum.x");
    TTreeReaderArray<float> rp_py_array(tree_reader, "ForwardRomanPotRecParticles.momentum.y");
    TTreeReaderArray<float> rp_pz_array(tree_reader, "ForwardRomanPotRecParticles.momentum.z");
    TTreeReaderArray<float> rp_mass_array(tree_reader, "ForwardRomanPotRecParticles.mass");
    TTreeReaderArray<int> rp_pdg_array(tree_reader, "ForwardRomanPotRecParticles.PDG");
    
    // Find beam particles first (same logic as original)
    P3MVector beame4_acc(0,0,0,0), beamp4_acc(0,0,0,0);
    const double fMass_proton = 0.938272;
    const double fMass_electron = 0.000511;
    
    while(tree_reader.Next()){
        for(int i = 0; i < mc_px_array.GetSize(); i++){
            if(mc_genStatus_array[i] != 4) continue;
            
            if(mc_pdg_array[i] == 2212){
                P3MVector p(mc_px_array[i], mc_py_array[i], mc_pz_array[i], fMass_proton);
                beamp4_acc += p;
            }
            else if(mc_pdg_array[i] == 11){
                P3MVector p(mc_px_array[i], mc_py_array[i], mc_pz_array[i], fMass_electron);
                beame4_acc += p;
            }
        }
    }
    
    auto nEntries = std::max<Long64_t>(1, events->GetEntries());
    P3MVector p_beam(beamp4_acc.X()/nEntries, beamp4_acc.Y()/nEntries, 
                     beamp4_acc.Z()/nEntries, fMass_proton);
    P3MVector e_beam(beame4_acc.X()/nEntries, beame4_acc.Y()/nEntries, 
                     beame4_acc.Z()/nEntries, fMass_electron);
    
    undoAfterburnAndCalc(p_beam, e_beam);
    
    // Reset reader for analysis
    tree_reader.Restart();
    
    // Analysis counters
    int total_events = 0;
    int events_with_mc_rp_protons = 0;
    int events_with_reco_rp_protons = 0;
    int events_with_multiple_mc_rp = 0;
    int events_with_multiple_reco_rp = 0;
    int total_mc_rp_protons = 0;
    int total_reco_rp_protons = 0;
    int total_matches = 0;
    
    const double RP_THETA_CUT = 0.005; // 5 mrad
    const double DELTA_R_CUT = 0.05;   // generous cut for all matches
    
    while(tree_reader.Next()){
        total_events++;
        
        // Collect MC protons in RP acceptance
        std::vector<P3MVector> mc_rp_protons;
        for(int i = 0; i < mc_px_array.GetSize(); i++){
            if(mc_genStatus_array[i] == 1 && mc_pdg_array[i] == 2212){
                P3MVector p(mc_px_array[i], mc_py_array[i], mc_pz_array[i], mc_mass_array[i]);
                undoAfterburn(p);
                
                if(p.Theta() < RP_THETA_CUT){
                    mc_rp_protons.push_back(p);
                    h_theta_mc_rp->Fill(p.Theta() * 1000.0);
                    h_eta_mc_rp->Fill(p.Eta());
                }
            }
        }
        
        // Collect reconstructed RP protons
        std::vector<P3MVector> reco_rp_protons;
        for(int i = 0; i < rp_px_array.GetSize(); i++){
            if(rp_pdg_array[i] == 2212){
                P3MVector p(rp_px_array[i], rp_py_array[i], rp_pz_array[i], rp_mass_array[i]);
                
                if(p.Theta() < RP_THETA_CUT){
                    reco_rp_protons.push_back(p);
                    h_theta_reco_rp->Fill(p.Theta() * 1000.0);
                    h_eta_reco_rp->Fill(p.Eta());
                }
            }
        }
        
        // Fill multiplicity histograms
        h_mc_multiplicity->Fill(mc_rp_protons.size());
        h_rp_multiplicity->Fill(reco_rp_protons.size());
        
        // Update counters
        total_mc_rp_protons += mc_rp_protons.size();
        total_reco_rp_protons += reco_rp_protons.size();
        
        if(mc_rp_protons.size() > 0) events_with_mc_rp_protons++;
        if(reco_rp_protons.size() > 0) events_with_reco_rp_protons++;
        if(mc_rp_protons.size() > 1) events_with_multiple_mc_rp++;
        if(reco_rp_protons.size() > 1) events_with_multiple_reco_rp++;
        
        // Calculate ΔR for all MC-Reco combinations
        for(size_t i = 0; i < reco_rp_protons.size(); i++){
            double best_deltaR = 1e9;
            int best_match = -1;
            
            for(size_t j = 0; j < mc_rp_protons.size(); j++){
                double deltaR = CalculateDeltaR(reco_rp_protons[i], mc_rp_protons[j]);
                h_deltaR_all->Fill(deltaR);
                h_deltaR_vs_theta->Fill(mc_rp_protons[j].Theta() * 1000.0, deltaR);
                
                if(deltaR < best_deltaR){
                    best_deltaR = deltaR;
                    best_match = j;
                }
            }
            
            if(best_match >= 0){
                h_deltaR_best->Fill(best_deltaR);
                if(best_deltaR < DELTA_R_CUT){
                    total_matches++;
                }
            }
        }
        
        // Progress indicator
        if(total_events % 10000 == 0){
            cout << "Processed " << total_events << " events..." << endl;
        }
    }
    
    // Print comprehensive results
    cout << "\n=== COMPREHENSIVE RP PROTON ANALYSIS RESULTS ===" << endl;
    cout << "Total events analyzed: " << total_events << endl;
    cout << "\n--- MC Truth Information ---" << endl;
    cout << "Events with MC RP protons: " << events_with_mc_rp_protons 
         << " (" << (double)events_with_mc_rp_protons/total_events*100 << "%)" << endl;
    cout << "Events with multiple MC RP protons: " << events_with_multiple_mc_rp 
         << " (" << (double)events_with_multiple_mc_rp/total_events*100 << "%)" << endl;
    cout << "Total MC RP protons: " << total_mc_rp_protons << endl;
    if(events_with_mc_rp_protons > 0){
        cout << "Average MC RP protons per event (given any): " 
             << (double)total_mc_rp_protons/events_with_mc_rp_protons << endl;
    }
    
    cout << "\n--- Reconstructed Information ---" << endl;
    cout << "Events with reconstructed RP protons: " << events_with_reco_rp_protons 
         << " (" << (double)events_with_reco_rp_protons/total_events*100 << "%)" << endl;
    cout << "Events with multiple reconstructed RP protons: " << events_with_multiple_reco_rp 
         << " (" << (double)events_with_multiple_reco_rp/total_events*100 << "%)" << endl;
    cout << "Total reconstructed RP protons: " << total_reco_rp_protons << endl;
    if(events_with_reco_rp_protons > 0){
        cout << "Average reconstructed RP protons per event (given any): " 
             << (double)total_reco_rp_protons/events_with_reco_rp_protons << endl;
    }
    
    cout << "\n--- Matching Information ---" << endl;
    cout << "Total good matches (#DeltaR < " << DELTA_R_CUT << "): " << total_matches << endl;
    if(total_reco_rp_protons > 0){
        cout << "Matching efficiency: " 
             << (double)total_matches/total_reco_rp_protons*100 << "%" << endl;
    }
    if(total_mc_rp_protons > 0){
        cout << "Reconstruction efficiency: " 
             << (double)total_matches/total_mc_rp_protons*100 << "%" << endl;
    }
    
    cout << "\n--- Analysis Settings ---" << endl;
    cout << "RP acceptance criterion: theta < " << RP_THETA_CUT << " rad (" 
         << RP_THETA_CUT*1000 << " mrad)" << endl;
    cout << "Delta R matching cut: " << DELTA_R_CUT << endl;
    
    // Save histograms
    TFile* outfile = new TFile("rp_proton_analysis.root", "RECREATE");
    h_deltaR_all->Write();
    h_deltaR_best->Write();
    h_theta_mc_rp->Write();
    h_theta_reco_rp->Write();
    h_eta_mc_rp->Write();
    h_eta_reco_rp->Write();
    h_deltaR_vs_theta->Write();
    h_rp_multiplicity->Write();
    h_mc_multiplicity->Write();
    outfile->Close();
    
    cout << "\nHistograms saved to rp_proton_analysis.root" << endl;
    cout << "Key plots:" << endl;
    cout << "- deltaR_all: All RP-MC DeltaR combinations" << endl;
    cout << "- deltaR_best: Best DeltaR match per RP proton" << endl;
    cout << "- theta_mc_rp vs theta_reco_rp: Angular distributions" << endl;
    cout << "- rp_multiplicity vs mc_multiplicity: Event multiplicities" << endl;
    
    delete events;
}

int main(int argc, char** argv){
    if(argc < 2){
        std::cerr << "Usage: analyze_rp_protons <filelist.txt>\n";
        return 1;
    }
    
    analyzeRPProtons(argv[1]);
    return 0;
}