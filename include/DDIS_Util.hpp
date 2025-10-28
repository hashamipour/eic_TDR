// ===============================
// File: DDIS_Util.hpp
// ===============================
#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>

#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"

// --- Short type aliases ---
using ROOT::Math::VectorUtil::boost;
using ROOT::Math::RotationX;
using ROOT::Math::RotationY;
using P3MVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzMVector>; // (px,py,pz; M)
using P3EVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>; // (px,py,pz; E)
using MomVector = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>, ROOT::Math::DefaultCoordinateSystemTag>;

static constexpr double MASS_PROTON   = 0.9382720813;  // GeV
static constexpr double MASS_ELECTRON = 0.00051099895; // GeV

// Global afterburner correction parameters (kept identical to your skimmers)
static float g_XAngle = -0.025f;            // rad
static RotationX g_rotX;                    // rotations are event-constant once beams set
static RotationY g_rotY;
static MomVector g_boostToCOM;              // boost to the Collider COM
static MomVector g_boostToHoF;              // boost to Head-on Frame (post-rotations)

// Convenience container for beam 4-vectors after undo-afterburner alignment
struct Beams {
  P3MVector p_beam; // incoming proton
  P3MVector e_beam; // incoming electron
  double    s  = 0.0; // (p+e)^2
  double    Mp2= MASS_PROTON*MASS_PROTON;
};

// -------------------------------------------------------------
// Numeric helpers
// -------------------------------------------------------------
static inline std::vector<Double_t> GetLogBins(double lo, double hi, int nbins){
  std::vector<Double_t> edges; edges.reserve(nbins+1);
  const double loglo = std::log10(lo), loghi = std::log10(hi);
  const double step  = (loghi-loglo)/nbins;
  for(int i=0;i<=nbins;++i) edges.push_back(std::pow(10.0, loglo + i*step));
  return edges;
}
static inline std::vector<Double_t> GetLinBins(double lo, double hi, int nbins){
  std::vector<Double_t> edges; edges.reserve(nbins+1);
  const double step = (hi-lo)/nbins;
  for(int i=0;i<=nbins;++i) edges.push_back(lo + i*step);
  return edges;
}
static inline std::vector<Double_t> GetRoundedLogBins(double lo, double hi, int nbins){
  // simple wrapper: compute log bins then round to one decimal for nicer labels
  auto v = GetLogBins(lo, hi, nbins);
  for(auto &x: v){ x = std::round(x*10.0)/10.0; }
  // ensure strictly increasing (repair accidental duplicates from rounding)
  for(size_t i=1;i<v.size();++i){ if(v[i]<=v[i-1]) v[i] = std::nextafter(v[i-1], std::numeric_limits<double>::infinity()); }
  return v;
}

// -------------------------------------------------------------
// Beam finding & afterburner handling
// -------------------------------------------------------------
static inline Beams FindAndAlignBeams(
  TTreeReaderArray<double>& mc_px,
  TTreeReaderArray<double>& mc_py,
  TTreeReaderArray<double>& mc_pz,
  TTreeReaderArray<double>& mc_mass,
  TTreeReaderArray<int>&    mc_genStatus,
  TTreeReaderArray<int>&    mc_pdg
){
  // Accumulate beam 4-vectors from status==4, PDG=11/2212
  P3MVector acc_p(0,0,0,0), acc_e(0,0,0,0);
  int np=0, ne=0;
  for(int i=0;i<mc_px.GetSize();++i){
    if(mc_genStatus[i]!=4) continue;
    if(mc_pdg[i]==2212){ acc_p += P3MVector(mc_px[i],mc_py[i],mc_pz[i],mc_mass[i]); ++np; }
    if(mc_pdg[i]==11  ){ acc_e += P3MVector(mc_px[i],mc_py[i],mc_pz[i],mc_mass[i]); ++ne; }
  }
  if(np>0) acc_p = P3MVector(acc_p.Px()/np, acc_p.Py()/np, acc_p.Pz()/np, MASS_PROTON);
  if(ne>0) acc_e = P3MVector(acc_e.Px()/ne, acc_e.Py()/ne, acc_e.Pz()/ne, MASS_ELECTRON);

  // Build idealized beam 4-vectors in lab (approx like your code)
  P3MVector p_beam(g_XAngle*acc_p.E(), 0.0,  acc_p.E(), MASS_PROTON);
  P3MVector e_beam(0.0,              0.0, -acc_e.E(), MASS_ELECTRON);

  // 1) boost to collider COM
  P3MVector com = p_beam + e_beam;
  g_boostToCOM.SetXYZ(-com.X()/com.E(), -com.Y()/com.E(), -com.Z()/com.E());
  p_beam = boost(p_beam, g_boostToCOM);
  e_beam = boost(e_beam, g_boostToCOM);

  // 2) align z with proton direction, then x with scattering plane (as in your skimmers)
  const double rotY = -1.0*TMath::ATan2(p_beam.X(), p_beam.Z());
  const double rotX =  1.0*TMath::ATan2(p_beam.Y(), p_beam.Z());
  g_rotY = RotationY(rotY);
  g_rotX = RotationX(rotX);
  p_beam = g_rotY(p_beam); p_beam = g_rotX(p_beam);
  e_beam = g_rotY(e_beam); e_beam = g_rotX(e_beam);

  // 3) final tiny boost to the head-on frame along +z of proton
  P3EVector hof(0,0,p_beam.Z(), p_beam.E());
  g_boostToHoF.SetXYZ(hof.X()/hof.E(), hof.Y()/hof.E(), hof.Z()/hof.E());
  p_beam = boost(p_beam, g_boostToHoF);
  e_beam = boost(e_beam, g_boostToHoF);

  Beams out; out.p_beam=p_beam; out.e_beam=e_beam; out.s=(p_beam+e_beam).M2();
  return out;
}

// static inline void undoAfterburn(P3MVector &p){
//   // rotate & boost the particle to head-on-like frame (same as beams)
//   p = g_rotY(p); p = g_rotX(p);
//   p = boost(p, g_boostToHoF);
// }

// -------------------------------------------------------------
// Event-level building blocks
// -------------------------------------------------------------
// static inline P3MVector BuildTruthXSystem(
//   TTreeReaderArray<double>& mc_px,
//   TTreeReaderArray<double>& mc_py,
//   TTreeReaderArray<double>& mc_pz,
//   TTreeReaderArray<double>& mc_mass,
//   TTreeReaderArray<int>&    mc_genStatus,
//   TTreeReaderArray<int>&    mc_pdg
// ){
//   P3MVector X(0,0,0,0);
//   for(int i=0;i<mc_px.GetSize();++i){
//     if(mc_genStatus[i]!=1) continue;                      // final state only
//     if(mc_pdg[i]==11 || mc_pdg[i]==2212) continue;       // exclude scattered e- and proton
//     P3MVector v(mc_px[i],mc_py[i],mc_pz[i],mc_mass[i]);
//     undoAfterburn(v);
//     X += v;
//   }
//   return X;
// }

// static inline double SumEPz_Truth_Lab(
//   TTreeReaderArray<double>& mc_px,
//   TTreeReaderArray<double>& mc_py,
//   TTreeReaderArray<double>& mc_pz,
//   TTreeReaderArray<double>& mc_mass,
//   TTreeReaderArray<int>&    mc_genStatus,
//   TTreeReaderArray<int>&    mc_pdg
// ){
//   double sum=0.0; // LAB FRAME (no corrections) to keep the expected ~2*E_e peak
//   for(int i=0;i<mc_px.GetSize();++i){
//     if(mc_genStatus[i]!=1) continue;
//     if(mc_pdg[i]==11 || mc_pdg[i]==2212) continue;
//     const double px=mc_px[i], py=mc_py[i], pz=mc_pz[i], m=mc_mass[i];
//     const double E = std::sqrt(px*px+py*py+pz*pz + m*m);
//     sum += (E - pz);
//   }
//   return sum;
// }

// static inline double SumEPz_Reco_Lab(
//   TTreeReaderArray<float>& re_px,
//   TTreeReaderArray<float>& re_py,
//   TTreeReaderArray<float>& re_pz,
//   TTreeReaderArray<float>& re_E,
//   TTreeReaderArray<int>*   re_pdg // optional (may be nullptr)
// ){
//   double sum=0.0;
//   const bool havePDG = (re_pdg!=nullptr);
//   for(int i=0;i<re_px.GetSize();++i){
//     if(havePDG){
//       const int pdg = (*re_pdg)[i];
//       if(pdg==11 || pdg==2212) continue; // exclude scattered e- and leading p if PDG is available
//     }
//     sum += (re_E[i] - re_pz[i]);
//   }
//   return sum;
// }

// // -------------------------------------------------------------
// // Kinematic calculators
// // -------------------------------------------------------------
// static inline double CalcT(const P3MVector &p_in, const P3MVector &p_out){
//   // Mandelstam t = (p_out - p_in)^2 with our (+,-,-,-) metric => negative; return |t|
//   const P3MVector q = p_out - p_in;
//   return -q.M2();
// }

// static inline double Calc_xL(const P3MVector &p_out, const P3MVector &p_beam){
//   return p_out.P() / p_beam.P();
// }

// static inline double Calc_MX2_from_formula(double s, double y, double xBj, double xL){
//   return s * y * (1.0 - xL - xBj);
// }

// static inline double Calc_xP(double MX2, double Q2, double t_abs, double W2){
//   // x_P = (M_X^2 + Q^2 - t) / (W^2 + Q^2 - M_p^2)
//   return (MX2 + Q2 - t_abs) / (W2 + Q2 - MASS_PROTON*MASS_PROTON);
// }

// // -------------------------------------------------------------
// // Matching helpers
// // -------------------------------------------------------------
// static inline int FindNearestTruthByTheta(const std::vector<P3MVector>& truths, const P3MVector &cand, double &dtheta){
//   int best=-1; dtheta = 1e9;
//   for(size_t i=0;i<truths.size();++i){
//     const double d = std::abs(truths[i].Theta() - cand.Theta());
//     if(d < dtheta){ dtheta = d; best = static_cast<int>(i); }
//   }
//   return best;
// }

// // Gather all final-state truth protons (after undoAfterburn)
// static inline std::vector<P3MVector> CollectTruthProtons(
//   TTreeReaderArray<double>& mc_px,
//   TTreeReaderArray<double>& mc_py,
//   TTreeReaderArray<double>& mc_pz,
//   TTreeReaderArray<double>& mc_mass,
//   TTreeReaderArray<int>&    mc_genStatus,
//   TTreeReaderArray<int>&    mc_pdg
// ){
//   std::vector<P3MVector> out; out.reserve(2);
//   for(int i=0;i<mc_px.GetSize();++i){
//     if(mc_genStatus[i]!=1 || mc_pdg[i]!=2212) continue;
//     P3MVector p(mc_px[i],mc_py[i],mc_pz[i],mc_mass[i]);
//     undoAfterburn(p);
//     out.push_back(p);
//   }
//   return out;
// }
