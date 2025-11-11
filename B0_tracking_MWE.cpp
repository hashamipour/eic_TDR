// Minimal Working Example: B0 Track Extraction


#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TMath.h>
#include <TH1D.h>
#include <iostream>

// ROOT Math vectors
#include "Math/Vector4D.h"
using P3MVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzMVector>;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root>" << std::endl;
        return 1;
    }

    TFile* inputFile = TFile::Open(argv[1]);
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << argv[1] << std::endl;
        return 1;
    }

    TTree* events = (TTree*)inputFile->Get("events");
    if (!events) {
        std::cerr << "Error: Could not find 'events' tree" << std::endl;
        return 1;
    }

    // Set up TTreeReader for the branches we need
    TTreeReader tree_reader(events);

    // MC particle information
    TTreeReaderArray<double> mc_px_array(tree_reader, "MCParticles.momentum.x");
    TTreeReaderArray<double> mc_py_array(tree_reader, "MCParticles.momentum.y");
    TTreeReaderArray<double> mc_pz_array(tree_reader, "MCParticles.momentum.z");
    TTreeReaderArray<double> mc_mass_array(tree_reader, "MCParticles.mass");
    TTreeReaderArray<int> mc_genStatus_array(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<int> mc_pdg_array(tree_reader, "MCParticles.PDG");

    // Truth-seeded reconstructed particles (used for B0)
    TTreeReaderArray<unsigned int> tsassoc_rec_id(tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.recID");
    TTreeReaderArray<unsigned int> tsassoc_sim_id(tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.simID");
    TTreeReaderArray<float> tsre_px_array(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x");
    TTreeReaderArray<float> tsre_py_array(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y");
    TTreeReaderArray<float> tsre_pz_array(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z");

    // Histograms
    TH1D* h_theta_all_TS = new TH1D("theta_all_TS", "All Truth-Seeded Protons;#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_B0 = new TH1D("theta_B0", "B0 Accepted Protons;#theta [mrad];Counts", 100, 0.0, 25.0);

    // Counters
    int n_total_protons = 0;
    int n_b0_accepted = 0;

    std::cout << "Processing events..." << std::endl;

    // Loop over events
    Long64_t nentries = events->GetEntries();
    while(tree_reader.Next()) {

        // Loop over truth-seeded reconstructed particles
        for(unsigned int j = 0; j < tsassoc_rec_id.GetSize(); j++) {
            auto mc_idx = tsassoc_sim_id[j];

            // Check if this is a final-state proton
            if(mc_idx >= (unsigned)mc_pdg_array.GetSize() ||
               mc_genStatus_array[mc_idx] != 1 ||
               mc_pdg_array[mc_idx] != 2212)
                continue;

            // Get reconstructed momentum
            P3MVector p_reco(tsre_px_array[j],
                            tsre_py_array[j],
                            tsre_pz_array[j],
                            mc_mass_array[mc_idx]);

            // Calculate theta angle in radians, convert to mrad
            double theta_rad = p_reco.Theta();
            double theta_mrad = theta_rad * 1000.0;

            n_total_protons++;
            h_theta_all_TS->Fill(theta_mrad);

            // B0 angular acceptance cut
            // B0 acceptance: 5.5 mrad < theta < 20 mrad
            if(theta_rad > 0.0055 && theta_rad < 0.02) {
                n_b0_accepted++;
                h_theta_B0->Fill(theta_mrad);
            }
        }
    }

    std::cout << "\n=== B0 Tracking Statistics ===" << std::endl;
    std::cout << "Total events processed: " << nentries << std::endl;
    std::cout << "Total truth-seeded protons found: " << n_total_protons << std::endl;
    std::cout << "Protons in B0 acceptance (5.5-20 mrad): " << n_b0_accepted << std::endl;
    std::cout << "B0 acceptance efficiency: " << (n_total_protons > 0 ? 100.0*n_b0_accepted/n_total_protons : 0.0) << "%" << std::endl;

    std::cout << "\nBranch used for B0 tracks: ReconstructedTruthSeededChargedParticles" << std::endl;
    std::cout << "Association branch: ReconstructedTruthSeededChargedParticleAssociations" << std::endl;

    // Save histograms
    TFile* outputFile = new TFile("B0_tracking_check.root", "RECREATE");
    h_theta_all_TS->Write();
    h_theta_B0->Write();
    outputFile->Close();

    std::cout << "\nOutput saved to: B0_tracking_check.root" << std::endl;

    inputFile->Close();
    return 0;
}
