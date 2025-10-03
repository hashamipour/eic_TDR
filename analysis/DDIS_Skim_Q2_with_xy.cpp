// g++ DDIS_Skim_Q2.cpp -o DDIS_Skim_Q2 $(root-config --cflags --glibs)
// ./DDIS_Skim_Q2 filelist.txt

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TMath.h>
#include <TObject.h>
#include <TString.h>
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include "Utility.hpp"

// These are the crucial headers for the ROOT::Math objects you are using.
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"

using ROOT::Math::VectorUtil::boost;
using ROOT::Math::VectorUtil::Angle;
// using ROOT::Math::RotationX;
// using ROOT::Math::RotationY;
using P3MVector=ROOT::Math::LorentzVector<ROOT::Math::PxPyPzMVector>;
using P3EVector=ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;
using MomVector=ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>,ROOT::Math::DefaultCoordinateSystemTag>;

// Functions: kinematic quantities
Double_t calcT(P3MVector p, P3MVector pprime);
Double_t calcQ2(P3MVector k, P3MVector kprime);
Double_t calcBjorkenX(P3MVector k, P3MVector kprime, P3MVector p);
Double_t calcM2Miss_3Body(P3MVector a, P3MVector b, P3MVector c, P3MVector d, P3MVector f);

// Functions: undo afterburner
void undoAfterburnAndCalc(P3MVector& p, P3MVector& k);
void undoAfterburn(P3MVector& p);
double rapidity(const P3MVector& p4);
double calculatePseudorapidity(const P3MVector& p4);

const Float_t fMass_proton{0.938272};
const Float_t fMass_electron{0.000511};
double MASS_PROTON   = fMass_proton;
double MASS_ELECTRON = fMass_electron;

static inline void SetNiceStyle(){
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(120);
  gStyle->SetPalette(kBlueRedYellow);
  gStyle->SetTitleFont(42, "XYZ"); gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleSize(0.048, "XYZ"); gStyle->SetLabelSize(0.038, "XYZ");
  gStyle->SetPadGridX(1); gStyle->SetPadGridY(1);
}

// Function to generate the logarithmic-reflected bin edges
std::vector<Double_t> GetLogReflectedBins(int n_log_bins, double min_range, double max_range) {
    
    // Create the positive side bin edges
    std::vector<Double_t> positive_edges;
    positive_edges.push_back(0.0); // Start at zero

    const double min_pos = min_range;
    const double max_pos = max_range;

    for (int i = 0; i <= n_log_bins; ++i) {
        // Logarithmic scale formula
        double log_edge = TMath::Log10(min_pos) + (TMath::Log10(max_pos) - TMath::Log10(min_pos)) * i / n_log_bins;
        positive_edges.push_back(TMath::Power(10, log_edge));
    }

    // Create the negative side by reflecting the positive edges
    std::vector<Double_t> negative_edges;
    for (size_t i = 1; i < positive_edges.size(); ++i) {
        negative_edges.push_back(-positive_edges[i]);
    }
    std::sort(negative_edges.begin(), negative_edges.end());

    // Combine the negative and positive edges
    std::vector<Double_t> all_edges;
    all_edges.insert(all_edges.end(), negative_edges.begin(), negative_edges.end());
    all_edges.insert(all_edges.end(), positive_edges.begin(), positive_edges.end());

    return all_edges;
}



/**
 * @brief Generates hybrid bin edges with two linear regions,
 * ensuring a bin is centered at zero. The central region is denser than the outer region.
 *
 * @param min_val The overall minimum value for the histogram range (e.g., -0.5).
 * @param max_val The overall maximum value for the histogram range (e.g., 0.5).
 * @param transition_point The absolute value at which to switch from the dense to the coarse linear binning (e.g., 0.05).
 * @param n_dense_bins_per_side The number of linear bins for each central, dense region. Must be an even number.
 * @param n_coarse_bins_per_side The number of linear bins for the outer, coarse region.
 * @return A vector of doubles containing the correctly ordered bin edges.
 */
std::vector<Double_t> GetHybridBins(
    double min_val, double max_val, double transition_point,
    int n_dense_bins_per_side, int n_coarse_bins_per_side) {
    
    // Ensure the number of dense bins for the central region is even to center a bin at zero
    if (n_dense_bins_per_side % 2 != 0) {
        n_dense_bins_per_side++;
    }

    // 1. Create all positive bin edges first
    std::vector<Double_t> positive_edges;
    
    // Dense linear part: from 0.0 to transition_point
    double dense_linear_step = transition_point / n_dense_bins_per_side;
    for (int i = 0; i <= n_dense_bins_per_side; i++) {
        positive_edges.push_back(i * dense_linear_step);
    }
    
    // Coarse linear part: from transition_point to max_val
    double coarse_linear_step = (max_val - transition_point) / n_coarse_bins_per_side;
    for (int i = 1; i <= n_coarse_bins_per_side; i++) {
        positive_edges.push_back(transition_point + i * coarse_linear_step);
    }
    
    // Sort the positive edges to ensure they are strictly increasing
    std::sort(positive_edges.begin(), positive_edges.end());
    // Remove any duplicates or near-duplicates at the transition point
    positive_edges.erase(std::unique(positive_edges.begin(), positive_edges.end()), positive_edges.end());
    
    // 2. Reflect the positive edges to create the negative ones
    std::vector<Double_t> negative_edges;
    for (size_t i = 0; i < positive_edges.size(); ++i) {
        negative_edges.push_back(-positive_edges[i]);
    }
    // Sort the negative edges to be in increasing order
    std::sort(negative_edges.begin(), negative_edges.end());

    // 3. Combine the two sets of edges with the central zero point
    std::vector<Double_t> all_edges;
    all_edges.insert(all_edges.end(), negative_edges.begin(), negative_edges.end());
    all_edges.insert(all_edges.end(), positive_edges.begin(), positive_edges.end());
    
    // CRITICAL: A final sort and unique to guarantee the entire vector is ordered.
    // This resolves any floating-point imprecision issues at the transitions.
    std::sort(all_edges.begin(), all_edges.end());
    all_edges.erase(std::unique(all_edges.begin(), all_edges.end()), all_edges.end());

    // std::cout << "Generated " << all_edges.size() - 1 << " bins with edges:\n";
    // for (const auto& edge : all_edges) {
        // std::cout << edge << " ";
    // }
    // std::cout << "\n";

    return all_edges;
}


/////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Converts a std::vector of doubles to a dynamically allocated C-style array.
 * * This function allocates new memory for the C-style array and copies the elements
 * from the input vector. It is the caller's responsibility to free this memory
 * using `delete[]` when it is no longer needed to prevent a memory leak.
 * * @param vec The input vector to be converted. It is passed by const reference
 * to avoid an unnecessary copy.
 * @param size A reference to an integer that will be updated with the size
 * of the created array. This is the number of elements in the array.
 * @return A pointer to the newly allocated C-style array of doubles.
 */
Double_t* vectorToCArray(const std::vector<Double_t>& vec, Int_t& size) {
    // Determine the size of the vector
    size = static_cast<Int_t>(vec.size());

    // If the vector is empty, return a null pointer and size 0.
    if (size == 0) {
        return nullptr;
    }

    // Dynamically allocate a new C-style array to hold the elements
    Double_t* cArray = new Double_t[size];

    // Copy the elements from the vector to the array
    for (Int_t i = 0; i < size; ++i) {
        cArray[i] = vec[i];
    }
    
    // Return the pointer to the new array
    return cArray;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <fileList.txt>" << std::endl;
        return 1;
    }
    TString fileList = argv[1];

    std::cout<< " __ __ __ __ __ __ __ __ __ __" <<std::endl;
    std::cout<< "|                             |"<<std::endl;
    std::cout<< "|     ePIC DDIS analysis      |"<<std::endl;
    std::cout<< "|     Inclusive Skimmer       |"<<std::endl;
    std::cout<< "|__ __ __ __ __ __ __ __ __ __|"<<std::endl;
    std::cout<< "\nInput filelist: " << fileList <<std::endl;

    std::ifstream fileListStream;
    fileListStream.open(fileList);
    std::string fileName;

    //---------------------------------------------------------
    // CREATE TCHAIN AND OUTPUT ROOT FILE
    //---------------------------------------------------------
    TChain* events = new TChain("events");
    Int_t nFiles{0};
    while(getline(fileListStream, fileName)){
        TString tmp = fileName;
        if (!std::filesystem::exists(tmp.Data())) {
        std::cerr << "\nError: File does not exist: " << fileName << std::endl;
        continue;
        }
        auto inputRootFile = TFile::Open(tmp);
        events->Add((TString)fileName);
        inputRootFile->Close();
        nFiles++;
    }
    std::cout<<"\nNo. of files: "<<nFiles<<"; no. of events: "<<events->GetEntries()<<std::endl;
    
    // Create output file
    TFile* outputFile = new TFile("DDIS_Skim_Q2_output.root", "RECREATE");

    // Create TTree for event-level data
    TTree* tree = new TTree("Q2_tree", "Q2 Kinematics Data");
    // Add x and y branches to output tree
    float out_x_truth, out_x_EM, out_x_DA, out_x_ESigma;
    float out_y_truth, out_y_EM, out_y_DA, out_y_ESigma;
    tree->Branch("x_truth",    &out_x_truth,   "x_truth/F");
    tree->Branch("x_EM",       &out_x_EM,      "x_EM/F");
    tree->Branch("x_DA",       &out_x_DA,      "x_DA/F");
    tree->Branch("x_ESigma",   &out_x_ESigma,  "x_ESigma/F");
    tree->Branch("y_truth",    &out_y_truth,   "y_truth/F");
    tree->Branch("y_EM",       &out_y_EM,      "y_EM/F");
    tree->Branch("y_DA",       &out_y_DA,      "y_DA/F");
    tree->Branch("y_ESigma",   &out_y_ESigma,  "y_ESigma/F");


    //---------------------------------------------------------
    // DECLARE OUTPUT HISTOGRAMS
    //---------------------------------------------------------

    int n_bins = 10;
    // std::vector<Double_t> bin_edges_Q2 = GetLogBins(1.0, 200.0, n_bins); // 20 bins from 0.5 to 200 GeV^2
    std::vector<Double_t> bin_edges_Q2 = GetRoundedLogBins(3.4, 150.0, n_bins); // 20 bins from 0.5 to 200 GeV^2
    // std::vector<Double_t> bin_edges_Q2 = GetManualQ2Bins(); // manually defined Q2 bins
    n_bins = bin_edges_Q2.size()-1; // update n_bins in case of rounding adjustments
    std::cout << "Number of Q2 bins: " << n_bins << std::endl;
    std::cout << "Q2 bin edges: ";
    for (const auto& edge : bin_edges_Q2) {
        std::cout << edge << " ";
    }
    std::cout << std::endl;
    
    // TH1D* h_Res_Q2_EM    = new TH1D("Q2_Res_EM",";Q^{2}(Reco)-Q^{2}(MC) [GeV^{2}]",n_bins, bin_edges_vec.data()); // .data() gives pointer to underlying array, bc TH1D constructor needs Double_t*
    TH1D* h_RelRes_Q2_EM = new TH1D("Q2_RelRes_EM","electron method;#frac{Q^{2}(Reco)-Q^{2}(MC)}{Q^{2}(MC)}",101,-0.15,0.15);

    // === New histograms for Bjorken-x and inelasticity y ===
    // x histograms (0..1)
    std::vector<Double_t> x_bins = GetLogBins(1.0e-4, 1.0, 25);
    std::cout << "bin edges for x: ";
    for (const auto& edge : x_bins) {
        std::cout << edge << " ";
    }
    std::cout << std::endl;
    TH1D* h_x_EM     = new TH1D("x_EM",     "electron method;x_{Bj}", x_bins.size()-1, x_bins.data());
    TH1D* h_x_DA     = new TH1D("x_DA",     "DA method;x_{Bj}",       x_bins.size()-1, x_bins.data());
    TH1D* h_x_ESigma = new TH1D("x_ESigma", "ESigma method;x_{Bj}",   x_bins.size()-1, x_bins.data());
    TH1D* h_x_truth  = new TH1D("x_truth",  "truth;x_{Bj}",           x_bins.size()-1, x_bins.data());

    // y histograms (0..1)
    int n_y_bins = 20;
    TH1D* h_y_EM     = new TH1D("y_EM",     "electron method;y", n_y_bins, 0.0, 1.0);
    TH1D* h_y_DA     = new TH1D("y_DA",     "DA method;y",       n_y_bins, 0.0, 1.0);
    TH1D* h_y_ESigma = new TH1D("y_ESigma", "ESigma method;y",   n_y_bins, 0.0, 1.0);
    TH1D* h_y_truth  = new TH1D("y_truth",  "truth;y",           n_y_bins, 0.0, 1.0);

    // Relative resolution for x and y
    TH1D* h_RelRes_x_EM     = new TH1D("x_RelRes_EM",     "electron method;#frac{x(Reco)-x(MC)}{x(MC)}", 101, -0.15, 0.15);
    TH1D* h_RelRes_x_DA     = new TH1D("x_RelRes_DA",     "DA method;#frac{x(Reco)-x(MC)}{x(MC)}",       101, -0.15, 0.15);
    TH1D* h_RelRes_x_ESigma = new TH1D("x_RelRes_ESigma", "ESigma method;#frac{x(Reco)-x(MC)}{x(MC)}",   101, -0.15, 0.15);

    TH1D* h_RelRes_y_EM     = new TH1D("y_RelRes_EM",     "electron method;#frac{y(Reco)-y(MC)}{y(MC)}", 101, -0.15, 0.15);
    TH1D* h_RelRes_y_DA     = new TH1D("y_RelRes_DA",     "DA method;#frac{y(Reco)-y(MC)}{y(MC)}",       101, -0.15, 0.15);
    TH1D* h_RelRes_y_ESigma = new TH1D("y_RelRes_ESigma", "ESigma method;#frac{y(Reco)-y(MC)}{y(MC)}",   101, -0.15, 0.15);

    // 2D binned relres vs truth
    // Use linear bins for x and y truth (0..1)
    int n_binned = 51;
    TH2D* h_RelRes_x_binned_EM     = new TH2D("x_RelRes_binned_EM",     "x: truth vs rel. res (EM);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",       
        x_bins.size()-1, x_bins.data(), n_binned, -0.15, 0.15);
    TH2D* h_RelRes_x_binned_DA     = new TH2D("x_RelRes_binned_DA",     "x: truth vs rel. res (DA);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",       
        x_bins.size()-1, x_bins.data(), n_binned, -0.15, 0.15);
    TH2D* h_RelRes_x_binned_ESigma = new TH2D("x_RelRes_binned_ESigma", "x: truth vs rel. res (ESigma);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",   
        x_bins.size()-1, x_bins.data(), n_binned, -0.15, 0.15);

    TH2D* h_RelRes_y_binned_EM     = new TH2D("y_RelRes_binned_EM",     "y: truth vs rel. res (EM);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",       n_y_bins, 0.0, 1.0, n_binned, -0.15, 0.15);
    TH2D* h_RelRes_y_binned_DA     = new TH2D("y_RelRes_binned_DA",     "y: truth vs rel. res (DA);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",       n_y_bins, 0.0, 1.0, n_binned, -0.15, 0.15);
    TH2D* h_RelRes_y_binned_ESigma = new TH2D("y_RelRes_binned_ESigma", "y: truth vs rel. res (ESigma);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",   n_y_bins, 0.0, 1.0, n_binned, -0.15, 0.15);

    // Correlation plots (truth vs reco)
    TH2D* h_Corr_x_EM     = new TH2D("x_Corr_EM",     "x correlation (EM);x_{truth};x_{EM}",           x_bins.size()-1, x_bins.data(), x_bins.size()-1, x_bins.data());
    TH2D* h_Corr_x_DA     = new TH2D("x_Corr_DA",     "x correlation (DA);x_{truth};x_{DA}",           x_bins.size()-1, x_bins.data(), x_bins.size()-1, x_bins.data());
    TH2D* h_Corr_x_ESigma = new TH2D("x_Corr_ESigma", "x correlation (ESigma);x_{truth};x_{ESigma}",   x_bins.size()-1, x_bins.data(), x_bins.size()-1, x_bins.data());

    TH2D* h_Corr_y_EM     = new TH2D("y_Corr_EM",     "y correlation (EM);y_{truth};y_{EM}",           n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);
    TH2D* h_Corr_y_DA     = new TH2D("y_Corr_DA",     "y correlation (DA);y_{truth};y_{DA}",           n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);
    TH2D* h_Corr_y_ESigma = new TH2D("y_Corr_ESigma", "y correlation (ESigma);y_{truth};y_{ESigma}",   n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);
    TH1D* h_RelRes_Q2_DA = new TH1D("Q2_RelRes_DA","DA method;#frac{Q^{2}(DA)-Q^{2}(MC)}{Q^{2}(MC)}",101,-0.15,0.15);
    TH1D* h_RelRes_Q2_ESigma = new TH1D("Q2_RelRes_ESigma","ESigma method;#frac{Q^{2}(ESigma)-Q^{2}(MC)}{Q^{2}(MC)}",101,-0.15,0.15);
    
    TH2D* h_RelRes_Q2_binned_EM = new TH2D("Q2_RelRes_binned_EM",";Q^{2} [GeV^{2}];#frac{Q^{2}(EM)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins, bin_edges_Q2.data(), n_binned,-0.15,0.15);
    TH2D* h_RelRes_Q2_binned_DA = new TH2D("Q2_RelRes_binned_DA",";Q^{2} [GeV^{2}];#frac{Q^{2}(DA)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins, bin_edges_Q2.data(), n_binned,-0.15,0.15);
    TH2D* h_RelRes_Q2_binned_ESigma = new TH2D("Q2_RelRes_binned_ESigma",";Q^{2} [GeV^{2}];#frac{Q^{2}(ESigma)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins, bin_edges_Q2.data(), n_binned,-0.15,0.15);

    
    TH2D* h_Corr_Q2_EM = new TH2D("Corr_Q2_EM", ";Q^{2}_{MC};Q^{2}_{EM}",
                                  n_bins, bin_edges_Q2.data(),
                                  n_bins, bin_edges_Q2.data());
    TH2D* h_Corr_Q2_DA = new TH2D("Corr_Q2_DA", ";Q^{2}_{MC};Q^{2}_{DA}",
                                    n_bins, bin_edges_Q2.data(),
                                    n_bins, bin_edges_Q2.data());
    TH2D* h_Corr_Q2_ESigma = new TH2D("Corr_Q2_ESigma", ";Q^{2}_{MC};Q^{2}_{ESigma}",
                                    n_bins, bin_edges_Q2.data(),
                                    n_bins, bin_edges_Q2.data());

    TH1D* h_Q2_truth    = new TH1D("h_Q2_truth","Q^2;# of events",n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_EM       = new TH1D("h_Q2_EM",";Q^{2}",n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_DA       = new TH1D("h_Q2_DA",";Q^{2}",n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_ESigma   = new TH1D("h_Q2_ESigma",";Q^{2}",n_bins, bin_edges_Q2.data());

    //---------------------------------------------------------
    // DECLARE TTREEREADER AND BRANCHES TO USE
    //---------------------------------------------------------
    TTreeReader tree_reader(events);
    TTreeReaderArray<double>  mc_px_array         = {tree_reader, "MCParticles.momentum.x"};
    TTreeReaderArray<double>  mc_py_array         = {tree_reader, "MCParticles.momentum.y"};
    TTreeReaderArray<double>  mc_pz_array         = {tree_reader, "MCParticles.momentum.z"};
    TTreeReaderArray<double> mc_mass_array        = {tree_reader, "MCParticles.mass"};
    TTreeReaderArray<int>    mc_genStatus_array   = {tree_reader, "MCParticles.generatorStatus"};
    TTreeReaderArray<int>    mc_pdg_array         = {tree_reader, "MCParticles.PDG"};
    TTreeReaderArray<unsigned int> assoc_rec_id   = {tree_reader, "ReconstructedParticleAssociations.recID"};
    TTreeReaderArray<unsigned int> assoc_sim_id   = {tree_reader, "ReconstructedParticleAssociations.simID"};
    TTreeReaderArray<float>  re_px_array          = {tree_reader, "ReconstructedParticles.momentum.x"};
    TTreeReaderArray<float>  re_py_array          = {tree_reader, "ReconstructedParticles.momentum.y"};
    TTreeReaderArray<float>  re_pz_array          = {tree_reader, "ReconstructedParticles.momentum.z"};
    TTreeReaderArray<int>    re_pdg_array         = {tree_reader, "ReconstructedParticles.PDG"};
    TTreeReaderArray<unsigned int> tsassoc_rec_id = {tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.recID"};
    TTreeReaderArray<unsigned int> tsassoc_sim_id = {tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.simID"};
    TTreeReaderArray<float>  tsre_px_array        = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x"};
    TTreeReaderArray<float>  tsre_py_array        = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y"};
    TTreeReaderArray<float>  tsre_pz_array        = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z"};
    TTreeReaderArray<float> rp_px_array           = {tree_reader, "ForwardRomanPotRecParticles.momentum.x"};
    TTreeReaderArray<float> rp_py_array           = {tree_reader, "ForwardRomanPotRecParticles.momentum.y"};
    TTreeReaderArray<float> rp_pz_array           = {tree_reader, "ForwardRomanPotRecParticles.momentum.z"};
    TTreeReaderArray<float> rp_mass_array         = {tree_reader, "ForwardRomanPotRecParticles.mass"};
    TTreeReaderArray<int>   rp_pdg_array          = {tree_reader, "ForwardRomanPotRecParticles.PDG"};

    // Declare variables to store branch data
    Float_t electron_Q2_EM;
    Float_t electron_Q2_DA;
    Float_t electron_Q2_ESigma;
    Float_t electron_Q2_truth;

    // New: Bjorken-x and inelasticity y (EM, DA, ESigma, Truth)
    Float_t electron_x_EM;
    Float_t electron_x_DA;
    Float_t electron_x_ESigma;
    Float_t electron_x_truth;

    Float_t electron_y_EM;
    Float_t electron_y_DA;
    Float_t electron_y_ESigma;
    Float_t electron_y_truth;

    // Set all branch addresses once
    events->SetBranchAddress("InclusiveKinematicsElectron.Q2", &electron_Q2_EM);
    events->SetBranchAddress("InclusiveKinematicsDA.Q2"      , &electron_Q2_DA);
    events->SetBranchAddress("InclusiveKinematicsESigma.Q2"  , &electron_Q2_ESigma);
    events->SetBranchAddress("InclusiveKinematicsTruth.Q2"   , &electron_Q2_truth);

    // New: x and y from the same InclusiveKinematics blocks
    events->SetBranchAddress("InclusiveKinematicsElectron.x", &electron_x_EM);
    events->SetBranchAddress("InclusiveKinematicsDA.x"      , &electron_x_DA);
    events->SetBranchAddress("InclusiveKinematicsESigma.x"  , &electron_x_ESigma);
    events->SetBranchAddress("InclusiveKinematicsTruth.x"   , &electron_x_truth);

    events->SetBranchAddress("InclusiveKinematicsElectron.y", &electron_y_EM);
    events->SetBranchAddress("InclusiveKinematicsDA.y"      , &electron_y_DA);
    events->SetBranchAddress("InclusiveKinematicsESigma.y"  , &electron_y_ESigma);
    events->SetBranchAddress("InclusiveKinematicsTruth.y"   , &electron_y_truth);

    // getting the maximum of true Q2 for binning the resolution vs true Q2
    // events->Draw("InclusiveKinematicsTruth.Q2 >> htemp", "", "goff");
    // std::cout << "Max Q2:" << ((TH1F*)gDirectory->Get("htemp"))->GetMaximum() << std::endl;

    //---------------------------------------------------------
    // RUN OVER TTREEREADER
    //---------------------------------------------------------
    Long64_t nentries = events->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        // Update the counter every 1000 events
        if (i % 1000 == 0) {
        printf("\rProcessing event %lld of %lld; %.2f percent done.", i, nentries, 100.0*i/nentries);
        fflush(stdout);
        }
        tree_reader.Next();

        // Loop over events **once** and fill all histograms
        events->GetEntry(i);

        // Fill histograms w/ different methods
        // x, y distributions
        h_x_EM->Fill(electron_x_EM);
        h_x_DA->Fill(electron_x_DA);
        h_x_ESigma->Fill(electron_x_ESigma);
        h_x_truth->Fill(electron_x_truth);

        h_y_EM->Fill(electron_y_EM);
        h_y_DA->Fill(electron_y_DA);
        h_y_ESigma->Fill(electron_y_ESigma);
        h_y_truth->Fill(electron_y_truth);

        // Relative resolutions vs truth (guard against zero)
        if (electron_x_truth != 0.0f) {
            h_RelRes_x_EM->Fill(    (electron_x_EM     - electron_x_truth)    / electron_x_truth);
            h_RelRes_x_DA->Fill(    (electron_x_DA     - electron_x_truth)    / electron_x_truth);
            h_RelRes_x_ESigma->Fill((electron_x_ESigma - electron_x_truth)    / electron_x_truth);
            h_RelRes_x_binned_EM->Fill(    electron_x_truth, (electron_x_EM     - electron_x_truth)/electron_x_truth);
            h_RelRes_x_binned_DA->Fill(    electron_x_truth, (electron_x_DA     - electron_x_truth)/electron_x_truth);
            h_RelRes_x_binned_ESigma->Fill(electron_x_truth, (electron_x_ESigma - electron_x_truth)/electron_x_truth);
        }
        if (electron_y_truth != 0.0f) {
            h_RelRes_y_EM->Fill(    (electron_y_EM     - electron_y_truth)    / electron_y_truth);
            h_RelRes_y_DA->Fill(    (electron_y_DA     - electron_y_truth)    / electron_y_truth);
            h_RelRes_y_ESigma->Fill((electron_y_ESigma - electron_y_truth)    / electron_y_truth);
            h_RelRes_y_binned_EM->Fill(    electron_y_truth, (electron_y_EM     - electron_y_truth)/electron_y_truth);
            h_RelRes_y_binned_DA->Fill(    electron_y_truth, (electron_y_DA     - electron_y_truth)/electron_y_truth);
            h_RelRes_y_binned_ESigma->Fill(electron_y_truth, (electron_y_ESigma - electron_y_truth)/electron_y_truth);
        }

        // Correlations
        h_Corr_x_EM->Fill(electron_x_truth, electron_x_EM);
        h_Corr_x_DA->Fill(electron_x_truth, electron_x_DA);
        h_Corr_x_ESigma->Fill(electron_x_truth, electron_x_ESigma);

        h_Corr_y_EM->Fill(electron_y_truth, electron_y_EM);
        h_Corr_y_DA->Fill(electron_y_truth, electron_y_DA);
        h_Corr_y_ESigma->Fill(electron_y_truth, electron_y_ESigma);

        h_Q2_EM->Fill(electron_Q2_EM);
        h_Q2_DA->Fill(electron_Q2_DA);
        h_Q2_ESigma->Fill(electron_Q2_ESigma);


        // Also fill the output tree for x and y
        out_x_truth  = electron_x_truth;
        out_x_EM     = electron_x_EM;
        out_x_DA     = electron_x_DA;
        out_x_ESigma = electron_x_ESigma;
        out_y_truth  = electron_y_truth;
        out_y_EM     = electron_y_EM;
        out_y_DA     = electron_y_DA;
        out_y_ESigma = electron_y_ESigma;
        tree->Fill();
        // Fill histograms for truth
        h_Q2_truth->Fill(electron_Q2_truth);

        // h_Res_Q2_EM->Fill(electron_Q2_EM - electron_Q2_truth);
        h_RelRes_Q2_EM->Fill((electron_Q2_EM - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_DA->Fill((electron_Q2_DA - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_ESigma->Fill((electron_Q2_ESigma - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_binned_EM->Fill(electron_Q2_truth, (electron_Q2_EM - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_binned_DA->Fill(electron_Q2_truth, (electron_Q2_DA - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_binned_ESigma->Fill(electron_Q2_truth, (electron_Q2_ESigma - electron_Q2_truth)/electron_Q2_truth);

        // Fill correlation histograms
        h_Corr_Q2_EM->Fill(electron_Q2_truth, electron_Q2_EM);
        h_Corr_Q2_DA->Fill(electron_Q2_truth, electron_Q2_DA);
        h_Corr_Q2_ESigma->Fill(electron_Q2_truth, electron_Q2_ESigma);
    }
    std::cout<<"\nDone looping over events.\n"<<std::endl;

    // Write all histograms and TTree to the output file
    outputFile->cd();
    tree->Write();
    // h_Res_Q2_EM->Write();

    h_RelRes_Q2_EM->Write();
    h_RelRes_Q2_DA->Write();
    h_RelRes_Q2_ESigma->Write();

    h_RelRes_Q2_binned_EM->Write();
    h_RelRes_Q2_binned_DA->Write();
    h_RelRes_Q2_binned_ESigma->Write();

    h_Corr_Q2_EM->Write();
    h_Corr_Q2_DA->Write();
    h_Corr_Q2_ESigma->Write();

    h_Q2_truth->Write();
    h_Q2_EM->Write();
    h_Q2_DA->Write();
    h_Q2_ESigma->Write();
    // Write x/y histograms
    h_x_truth->Write();  h_x_EM->Write();  h_x_DA->Write();  h_x_ESigma->Write();
    h_y_truth->Write();  h_y_EM->Write();  h_y_DA->Write();  h_y_ESigma->Write();

    h_RelRes_x_EM->Write();
    h_RelRes_x_DA->Write();
    h_RelRes_x_ESigma->Write();

    h_RelRes_y_EM->Write();
    h_RelRes_y_DA->Write();
    h_RelRes_y_ESigma->Write();

    h_RelRes_x_binned_EM->Write();
    h_RelRes_x_binned_DA->Write();
    h_RelRes_x_binned_ESigma->Write();

    h_RelRes_y_binned_EM->Write();
    h_RelRes_y_binned_DA->Write();
    h_RelRes_y_binned_ESigma->Write();

    h_Corr_x_EM->Write();
    h_Corr_x_DA->Write();
    h_Corr_x_ESigma->Write();

    h_Corr_y_EM->Write();
    h_Corr_y_DA->Write();
    h_Corr_y_ESigma->Write();


    outputFile->Close();
    delete events;
    delete outputFile;
    return 0;
}