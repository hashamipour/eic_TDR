#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>
#include <string>

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
#include "Math/Vector3D.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/GenVector/Boost.h"
#include "Math/VectorUtil.h"
#include <TLine.h>
#include <TLatex.h>

using std::cout; using std::endl;
using ROOT::Math::VectorUtil::boost;
using ROOT::Math::RotationX;
using ROOT::Math::RotationY;
using P3MVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzMVector>;
using P3EVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;
using MomVector = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>, ROOT::Math::DefaultCoordinateSystemTag>;
using ROOT::Math::PxPyPzEVector;

// FIX 1: Add missing coordinate transformation functions from DDIS_TDR.cpp
// Objects for undoing afterburn boost
Float_t fXAngle{-0.025}; // Crossing angle in radians
Float_t fRotX{};
RotationX rotAboutX;
Float_t fRotY{};
RotationY rotAboutY;
MomVector vBoostToCoM;
MomVector vBoostToHoF;

// Undo AB and calculate boost vectors - DO THIS FIRST FOR EACH EVENT
void undoAfterburnAndCalc(P3MVector& p, P3MVector& k){
  // Holding vectors for beam - undoing crossing angle ONLY
  P3MVector p_beam(fXAngle*p.E(), 0., p.E(), p.M());
  P3MVector e_beam(0., 0., -k.E(), k.M());
  
  // Define boost vector to CoM frame
  P3MVector CoM_boost = p_beam+e_beam;
  vBoostToCoM.SetXYZ(-CoM_boost.X()/CoM_boost.E(), -CoM_boost.Y()/CoM_boost.E(), -CoM_boost.Z()/CoM_boost.E());
  
  // Apply boost to beam vectors
  p_beam = boost(p_beam, vBoostToCoM);
  e_beam = boost(e_beam, vBoostToCoM);
  
  // Calculate rotation angles and create rotation objects
  fRotY = -1.0*TMath::ATan2(p_beam.X(), p_beam.Z());
  fRotX = 1.0*TMath::ATan2(p_beam.Y(), p_beam.Z());

  rotAboutY = RotationY(fRotY);
  rotAboutX = RotationX(fRotX);

  // Apply rotation to beam vectors
  p_beam = rotAboutY(p_beam);
  p_beam = rotAboutX(p_beam);
  e_beam = rotAboutY(e_beam);
  e_beam = rotAboutX(e_beam);

  // Define boost vector back to head-on frame
  P3EVector HoF_boost(0., 0., CoM_boost.Z(), CoM_boost.E());
  vBoostToHoF.SetXYZ(HoF_boost.X()/HoF_boost.E(), HoF_boost.Y()/HoF_boost.E(), HoF_boost.Z()/HoF_boost.E());

  // Apply boost back to head on frame to beam vectors
  p_beam = boost(p_beam, vBoostToHoF);
  e_beam = boost(e_beam, vBoostToHoF);

  // Make changes to input vectors
  p.SetPxPyPzE(p_beam.X(), p_beam.Y(), p_beam.Z(), p_beam.E());
  k.SetPxPyPzE(e_beam.X(), e_beam.Y(), e_beam.Z(), e_beam.E());
}

// Undo afterburn procedure only
void undoAfterburn(P3MVector& a){
  // Undo AB procedure for single vector, a^{mu}
  a = boost(a, vBoostToCoM); // BOOST TO COM FRAME
  a = rotAboutY(a);          // ROTATE TO Z-AXIS
  a = rotAboutX(a);          // ROTATE TO Z-AXIS
  a = boost(a, vBoostToHoF); // BOOST BACK TO HEAD ON FRAME
}

void analyzeProtonsMandelstamT(TString fileList){
    cout<<"----------------------------"<<endl;
    cout<<"                            "<<endl;
    cout<<"   Proton Mandelstam t Analysis   "<<endl;
    cout<<"                            "<<endl;
    cout<<"----------------------------"<<endl;
    cout<<"\nInput filelist: "<<fileList<<endl;

    // Constants
    const Double_t fMass_proton   = 0.938272;
    const Double_t fMass_electron = 0.000511;

    //---------------------------------------------------------
    // CREATE TCHAIN FROM INPUT ROOT FILES
    //---------------------------------------------------------
    TChain* events = new TChain("events");
    Int_t nFiles = 0;

    std::ifstream fileListStream;
    fileListStream.open(fileList);
    std::string fileName;

    while(getline(fileListStream, fileName)){
        events->Add((TString)fileName);
        nFiles++;
    }
    cout<<"\nNo. of files: "<<nFiles<<"; no. of events: "<<events->GetEntries()<<endl;

    //---------------------------------------------------------
    // DECLARE OUTPUT HISTOGRAMS
    //---------------------------------------------------------
    // 2D correlation plots (Truth vs Reco Mandelstam t)
    TH2D* h_t_corr_B0 = new TH2D("t_corr_B0", "B0 Protons: Truth vs Reco Mandelstam t;Truth |t| [GeV^{2}];Reco |t| [GeV^{2}]", 
                                  100, 0.0, 2.0, 100, 0.0, 2.0);
    TH2D* h_t_corr_RP = new TH2D("t_corr_RP", "Roman Pot Protons: Truth vs Reco Mandelstam t;Truth |t| [GeV^{2}];Reco |t| [GeV^{2}]", 
                                  100, 0.0, 0.5, 100, 0.0, 0.5);

    // 1D histograms for superimposed comparison
    TH1D* h_t_MC      = new TH1D("t_MC",      "Truth Mandelstam t;|t| [GeV^{2}];Counts", 100, 0.0, 2.0);
    TH1D* h_t_B0      = new TH1D("t_B0",      "B0 Reco Mandelstam t;|t| [GeV^{2}];Counts", 100, 0.0, 2.0);
    TH1D* h_t_RP_histo= new TH1D("t_RP_histo","Roman Pot Reco Mandelstam t;|t| [GeV^{2}];Counts", 100, 0.0, 2.0);

    // Resolution plots
    TH1D* h_t_res_B0 = new TH1D("t_res_B0", "B0 Mandelstam t Resolution;(|t|_{reco} - |t|_{truth})/|t|_{truth}", 100, -2.0, 2.0);
    TH1D* h_t_res_RP = new TH1D("t_res_RP", "RP Mandelstam t Resolution;(|t|_{reco} - |t|_{truth})/|t|_{truth}", 100, -2.0, 2.0);

    // Angular distributions for verification
    TH1D* h_theta_MC = new TH1D("theta_MC", "MC Proton Scattering Angle;#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_B0 = new TH1D("theta_B0", "B0 Proton Scattering Angle;#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_RP = new TH1D("theta_RP", "RP Proton Scattering Angle;#theta [mrad];Counts", 100, 0.0, 25.0);

    //---------------------------------------------------------
    // DECLARE TTREEREADER AND BRANCHES TO USE
    //---------------------------------------------------------
    TTreeReader tree_reader(events);

    // --- MC particles (edm4hep/podio split branches) ---
    // Use Double_t because the file reports Double_t leaves for momentum.*
    TTreeReaderArray<double> mc_px_array        = {tree_reader, "MCParticles.momentum.x"};
    TTreeReaderArray<double> mc_py_array        = {tree_reader, "MCParticles.momentum.y"};
    TTreeReaderArray<double> mc_pz_array        = {tree_reader, "MCParticles.momentum.z"};
    TTreeReaderArray<double> mc_mass_array      = {tree_reader, "MCParticles.mass"};
    TTreeReaderArray<int>    mc_genStatus_array = {tree_reader, "MCParticles.generatorStatus"};
    TTreeReaderArray<int>    mc_pdg_array       = {tree_reader, "MCParticles.PDG"};

    // --- Truth-seeded charged RECO (B0 acceptance window later by angle) ---
    TTreeReaderArray<unsigned int> tsassoc_rec_id = {tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.recID"};
    TTreeReaderArray<unsigned int> tsassoc_sim_id = {tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.simID"};

    TTreeReaderArray<float> tsre_px_array     = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x"};
    TTreeReaderArray<float> tsre_py_array     = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y"};
    TTreeReaderArray<float> tsre_pz_array     = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z"};
    TTreeReaderArray<float> tsre_e_array      = {tree_reader, "ReconstructedTruthSeededChargedParticles.energy"};

    // --- Roman Pot RECO (very small angles) ---
    TTreeReaderArray<float> rp_px_array     = {tree_reader, "ForwardRomanPotRecParticles.momentum.x"};
    TTreeReaderArray<float> rp_py_array     = {tree_reader, "ForwardRomanPotRecParticles.momentum.y"};
    TTreeReaderArray<float> rp_pz_array     = {tree_reader, "ForwardRomanPotRecParticles.momentum.z"};
    TTreeReaderArray<float> rp_mass_array   = {tree_reader, "ForwardRomanPotRecParticles.mass"};
    TTreeReaderArray<int>    rp_pdg_array    = {tree_reader, "ForwardRomanPotRecParticles.PDG"};

    //---------------------------------------------------------
    // FIX 2: PROPER TWO-PASS APPROACH - FIRST PASS: FIND BEAM PARTICLES
    //---------------------------------------------------------
    cout << "Finding beam particles." << endl;

    P3MVector beame4(0,0,0,-1), beamp4(0,0,0,-1);
    P3MVector beame4_acc(0,0,0,-1), beamp4_acc(0,0,0,-1);

    while (tree_reader.Next()){
        P3MVector beame4_evt(0,0,0,-1), beamp4_evt(0,0,0,-1);
        TVector3 mctrk;

        for(int imc=0; imc<mc_px_array.GetSize(); imc++){
            mctrk.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);
            if(mc_genStatus_array[imc] == 4 && mc_pdg_array[imc] == 2212) {
                beamp4_evt.SetCoordinates(mctrk.X(), mctrk.Y(), mctrk.Z(), fMass_proton);
            }
            if(mc_genStatus_array[imc] == 4 && mc_pdg_array[imc] == 11) {
                beame4_evt.SetCoordinates(mctrk.X(), mctrk.Y(), mctrk.Z(), fMass_electron);
            }
        }
        beame4_acc += beame4_evt;
        beamp4_acc += beamp4_evt;
    }

    // Average beam momenta over entries
    auto nEntries = std::max<Long64_t>(1, events->GetEntries());
    beame4.SetCoordinates(beame4_acc.X()/nEntries, beame4_acc.Y()/nEntries, beame4_acc.Z()/nEntries, beame4_acc.M()/nEntries);
    beamp4.SetCoordinates(beamp4_acc.X()/nEntries, beamp4_acc.Y()/nEntries, beamp4_acc.Z()/nEntries, beamp4_acc.M()/nEntries);

    cout << "[DEBUG] Found beam energies " << beame4.E() << "x" << beamp4.E() << " GeV" << endl;

    // FIX 3: Apply coordinate transformation setup
    undoAfterburnAndCalc(beamp4, beame4);

    //---------------------------------------------------------
    // SECOND PASS: CALCULATE |t| FOR TRUTH, B0, AND RP
    //---------------------------------------------------------
    tree_reader.Restart();

    int n_truth_protons = 0;
    int n_b0_matches = 0;
    int n_rp_matches = 0;

    while (tree_reader.Next()){
        std::vector<P3MVector> scatp4_gen; // truth
        std::vector<P3MVector> scatp4_b0;  // reco B0
        std::vector<P3MVector> scatp4_rp;  // reco RP

        TVector3 mctrk, recotrk;

        // FIX 4: Process truth particles FIRST to fill scatp4_gen
        for(int imc=0; imc<mc_px_array.GetSize(); imc++){
            if(mc_genStatus_array[imc] == 1 && mc_pdg_array[imc] == 2212){
                mctrk.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);
                P3MVector q_truth(mctrk.X(), mctrk.Y(), mctrk.Z(), mc_mass_array[imc]);
                
                // FIX 5: Apply coordinate transformation
                undoAfterburn(q_truth);
                
                scatp4_gen.push_back(q_truth);
                h_theta_MC->Fill(q_truth.Theta()*1000.0);
                auto t_truth = (q_truth - beamp4).M2();
                h_t_MC->Fill(TMath::Abs(t_truth));
            }
        }
        n_truth_protons += (int)scatp4_gen.size();

        // 2) reco B0 via truth-seeded association (5.5-20 mrad)
        for(unsigned int i=0; i<tsassoc_rec_id.GetSize(); i++){
            auto mc_idx = tsassoc_sim_id[i];
            if(mc_idx < (unsigned)mc_pdg_array.GetSize() && mc_genStatus_array[mc_idx]==1 && mc_pdg_array[mc_idx]==2212){
                recotrk.SetXYZ(tsre_px_array[i], tsre_py_array[i], tsre_pz_array[i]);
                P3MVector q_reco(recotrk.X(), recotrk.Y(), recotrk.Z(), mc_mass_array[mc_idx]);
                
                // FIX 6: Apply coordinate transformation
                undoAfterburn(q_reco);
                
                if(q_reco.Theta()>0.0055 && q_reco.Theta()<0.02){
                    scatp4_b0.push_back(q_reco);
                    h_theta_B0->Fill(q_reco.Theta()*1000.0);
                    auto t_reco = (q_reco - beamp4).M2();
                    h_t_B0->Fill(TMath::Abs(t_reco));
                    
                    // truth partner
                    mctrk.SetXYZ(mc_px_array[mc_idx], mc_py_array[mc_idx], mc_pz_array[mc_idx]);
                    P3MVector q_truth(mctrk.X(), mctrk.Y(), mctrk.Z(), mc_mass_array[mc_idx]);
                    
                    // FIX 7: Apply coordinate transformation to truth partner
                    undoAfterburn(q_truth);
                    
                    auto t_truth = (q_truth - beamp4).M2();
                    h_t_corr_B0->Fill(TMath::Abs(t_truth), TMath::Abs(t_reco));
                    if(TMath::Abs(t_truth)>1e-6) h_t_res_B0->Fill( (TMath::Abs(t_reco)-TMath::Abs(t_truth))/TMath::Abs(t_truth) );
                    n_b0_matches++;
                }
            }
        }

        // FIX 8: 3) reco RP - NO COORDINATE TRANSFORMATION (as per DDIS_TDR.cpp comment)
        for(int ir=0; ir<rp_px_array.GetSize(); ir++){
            if(rp_pdg_array[ir]==2212){
                recotrk.SetXYZ(rp_px_array[ir], rp_py_array[ir], rp_pz_array[ir]);
                P3MVector q_rp(recotrk.X(), recotrk.Y(), recotrk.Z(), rp_mass_array[ir]);
                // Note: NO undoAfterburn for RP as per DDIS_TDR.cpp: "NO NEED TO UNDO AFTERBURNER FOR FF DETECTORS"
                
                if(q_rp.Theta()<0.005){
                    scatp4_rp.push_back(q_rp);
                    h_theta_RP->Fill(q_rp.Theta()*1000.0);
                    auto t_reco = (q_rp - beamp4).M2();
                    h_t_RP_histo->Fill(TMath::Abs(t_reco));
                    
                    // FIX 9: match to truth by angle - now scatp4_gen is properly filled
                    double min_dtheta=1e9; int best=-1;
                    for(size_t it=0; it<scatp4_gen.size(); ++it){
                        double dth = TMath::Abs(q_rp.Theta()-scatp4_gen[it].Theta());
                        if(dth < min_dtheta && scatp4_gen[it].Theta()<0.006){
                            min_dtheta=dth; best=(int)it;
                        }
                    }
                    if(best>=0 && min_dtheta<0.001){
                        auto t_truth = (scatp4_gen[best] - beamp4).M2();
                        h_t_corr_RP->Fill(TMath::Abs(t_truth), TMath::Abs(t_reco));
                        if(TMath::Abs(t_truth)>1e-6) h_t_res_RP->Fill( (TMath::Abs(t_reco)-TMath::Abs(t_truth))/TMath::Abs(t_truth) );
                        n_rp_matches++;
                    }
                }
            }
        }
    }

    //---------------------------------------------------------
    // PRINT STATS
    //---------------------------------------------------------
    cout << "\n=== ANALYSIS RESULTS ===" << endl;
    cout << "Truth protons found: " << n_truth_protons << endl;
    cout << "B0 matched protons: " << n_b0_matches << endl;
    cout << "RP matched protons: " << n_rp_matches << endl;
    if(n_truth_protons>0){
        cout << "B0 matching efficiency: " << (double)n_b0_matches/n_truth_protons*100 << "%" << endl;
        cout << "RP matching efficiency: " << (double)n_rp_matches/n_truth_protons*100 << "%" << endl;
    }

    //---------------------------------------------------------
    // DRAW
    //---------------------------------------------------------
    gStyle->SetOptStat(0);

    // Correlations
    TCanvas* c1 = new TCanvas("c1","Mandelstam t Correlations",1400,600);
    c1->Divide(2,1);
    c1->cd(1);
    h_t_corr_B0->Draw("COLZ");
    h_t_corr_B0->SetTitle("B0 Protons: Truth vs Reco |t|;Truth |t| [GeV^{2}];Reco |t| [GeV^{2}]");
    { TLine* l=new TLine(0,0,2,2); l->SetLineColor(kRed); l->SetLineWidth(2); l->SetLineStyle(2); l->Draw("same"); }
    c1->cd(2);
    h_t_corr_RP->Draw("COLZ");
    h_t_corr_RP->SetTitle("Roman Pot Protons: Truth vs Reco |t|;Truth |t| [GeV^{2}];Reco |t| [GeV^{2}]");
    { TLine* l=new TLine(0,0,0.5,0.5); l->SetLineColor(kRed); l->SetLineWidth(2); l->SetLineStyle(2); l->Draw("same"); }

    // Distributions
    TCanvas* c2 = new TCanvas("c2","Mandelstam t Distributions",1200,800);
    gPad->SetLogy(1);
    h_t_MC->SetMinimum(1);
    h_t_MC->SetLineColor(kBlack); h_t_MC->SetLineWidth(2); h_t_MC->GetXaxis()->SetTitle("|t| [GeV^{2}]"); h_t_MC->GetYaxis()->SetTitle("Counts");
    h_t_MC->Draw("HIST");
    h_t_B0->SetLineColor(kBlue); h_t_B0->SetMarkerColor(kBlue); h_t_B0->SetMarkerStyle(20); h_t_B0->Draw("PE SAME");
    h_t_RP_histo->SetLineColor(kCyan+1); h_t_RP_histo->SetMarkerColor(kCyan+1); h_t_RP_histo->SetMarkerStyle(20); h_t_RP_histo->Draw("PE SAME");
    { TLegend* leg = new TLegend(0.6,0.7,0.89,0.89); leg->SetLineColorAlpha(kWhite,0); leg->SetFillColorAlpha(kWhite,0);
      leg->AddEntry(h_t_MC, "Truth MC", "l");
      leg->AddEntry(h_t_B0, "Reco B0 (5.5-20 mrad)", "lp");
      leg->AddEntry(h_t_RP_histo, "Reco RP (<5 mrad)", "lp");
      leg->Draw(); }

    // Resolution
    TCanvas* c4 = new TCanvas("c4","Mandelstam t Resolution",1200,600);
    c4->Divide(2,1);
    c4->cd(1); h_t_res_B0->SetLineColor(kBlue); h_t_res_B0->SetLineWidth(2); h_t_res_B0->Draw("HIST"); h_t_res_B0->SetTitle("B0 |t| Resolution;(|t|_{reco} - |t|_{truth})/|t|_{truth};Counts");
    c4->cd(2); h_t_res_RP->SetLineColor(kCyan+1); h_t_res_RP->SetLineWidth(2); h_t_res_RP->Draw("HIST"); h_t_res_RP->SetTitle("RP |t| Resolution;(|t|_{reco} - |t|_{truth})/|t|_{truth};Counts");

    //---------------------------------------------------------
    // SAVE RESULTS
    //---------------------------------------------------------
    TFile* outfile = new TFile("proton_mandelstam_analysis.root", "RECREATE");
    h_t_corr_B0->Write(); h_t_corr_RP->Write(); h_t_MC->Write(); h_t_B0->Write(); h_t_RP_histo->Write();
    h_theta_MC->Write(); h_theta_B0->Write(); h_theta_RP->Write(); h_t_res_B0->Write(); h_t_res_RP->Write();
    outfile->Close();
    c1->SaveAs("mandelstam_t_correlations.png");
    c2->SaveAs("mandelstam_t_distributions.png");
    c4->SaveAs("mandelstam_t_resolution.png");

    cout << "\nAnalysis complete! Files saved:" << endl;
    cout << "- proton_mandelstam_analysis.root (histograms)" << endl;
    cout << "- mandelstam_t_correlations.png" << endl;
    cout << "- mandelstam_t_distributions.png" << endl;
    cout << "- mandelstam_t_resolution.png" << endl;

    return;
}

// -------------------------
// Entry point
// -------------------------
int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: skim_t <filelist.txt>\n";
    return 1;
  }

  cout << "AFTERBURNER_FIX_APPLIED" << endl;
  TString fileList = argv[1];
  analyzeProtonsMandelstamT(fileList);

  return 0;
}