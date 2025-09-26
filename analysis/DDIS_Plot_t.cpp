// Plot_t_from_skim.cpp
// Read the skimmed file (tree "skim") and make 1D t/|t| plots and 2D truth↔reco correlations.
// Usage:
//   g++ -O2 -std=c++17 Plot_t_from_skim.cpp -o plot_t `root-config --cflags --libs`
//   ./plot_t DDIS_Skim_t_xp_output.root plots_t.pdf histos_t.root

#include <iostream>
#include <cmath>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>

static inline bool isfinitef(float x){ return std::isfinite(x); }
static inline double absd(double x){ return x<0? -x : x; }

int main(int argc, char** argv){
  if(argc < 2){
    std::cerr << "Usage: " << argv[0] << " skim.root [out.pdf] [out.root]\n";
    return 1;
  }
  const std::string inroot  = argv[1];
  const std::string outpdf  = (argc>2 ? argv[2] : "plots_t.pdf");
  const std::string outroot = (argc>3 ? argv[3] : "histos_t.root");

  // open
  TFile fin(inroot.c_str(), "READ");
  if(fin.IsZombie()){
    std::cerr << "ERROR: cannot open input: " << inroot << "\n";
    return 2;
  }
  TTree* skim = dynamic_cast<TTree*>(fin.Get("skim"));
  if(!skim){
    std::cerr << "ERROR: tree 'skim' not found in " << inroot << "\n";
    return 3;
  }

  // reader
  TTreeReader r(skim);
  // truth
  TTreeReaderValue<float> t_true_after   (r, "t_true_after");
  TTreeReaderValue<float> t_abs_true_after(r, "t_abs_true_after");
  TTreeReaderValue<float> t_true_post    (r, "t_true_post");
  TTreeReaderValue<float> t_abs_true_post(r, "t_abs_true_post");
  // reco RP
  TTreeReaderValue<float> t_rp_after     (r, "t_rp_after");
  TTreeReaderValue<float> t_abs_rp_after (r, "t_abs_rp_after");
  TTreeReaderValue<float> t_rp_post      (r, "t_rp_post");
  TTreeReaderValue<float> t_abs_rp_post  (r, "t_abs_rp_post");
  // reco B0
  TTreeReaderValue<float> t_b0_after     (r, "t_b0_after");
  TTreeReaderValue<float> t_abs_b0_after (r, "t_abs_b0_after");
  TTreeReaderValue<float> t_b0_post      (r, "t_b0_post");
  TTreeReaderValue<float> t_abs_b0_post  (r, "t_abs_b0_post");
  // flags
  TTreeReaderValue<Bool_t> has_rp_reco   (r, "has_rp_reco");
  TTreeReaderValue<Bool_t> has_b0_reco   (r, "has_b0_reco");

  // style
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(255);

  // binning (t in GeV^2)
  const int    nb = 60;
  const double tmin = -1.0;   // signed t (can be negative)
  const double tmax =  0.2;
  const int    nbabs = 60;
  const double abmin = 0.0;   // |t|
  const double abmax = 1.5;

  // 1D histos
  auto h_t_true_after   = new TH1D("h_t_true_after",   "t (truth afterburner);t [GeV^{2}];Events", nb, tmin, tmax);
  auto h_t_true_post    = new TH1D("h_t_true_post",    "t (truth postburn);t [GeV^{2}];Events",    nb, tmin, tmax);
  auto h_t_rp_after     = new TH1D("h_t_rp_after",     "t (RP reco, afterburner);t [GeV^{2}];Events", nb, tmin, tmax);
  auto h_t_rp_post      = new TH1D("h_t_rp_post",      "t (RP reco, postburn);t [GeV^{2}];Events",    nb, tmin, tmax);
  auto h_t_b0_after     = new TH1D("h_t_b0_after",     "t (B0 reco, afterburner);t [GeV^{2}];Events", nb, tmin, tmax);
  auto h_t_b0_post      = new TH1D("h_t_b0_post",      "t (B0 reco, postburn);t [GeV^{2}];Events",    nb, tmin, tmax);

  auto h_tab_true_after = new TH1D("h_tab_true_after", "|t| (truth afterburner);|t| [GeV^{2}];Events", nbabs, abmin, abmax);
  auto h_tab_true_post  = new TH1D("h_tab_true_post",  "|t| (truth postburn);|t| [GeV^{2}];Events",    nbabs, abmin, abmax);
  auto h_tab_rp_after   = new TH1D("h_tab_rp_after",   "|t| (RP reco, afterburner);|t| [GeV^{2}];Events", nbabs, abmin, abmax);
  auto h_tab_rp_post    = new TH1D("h_tab_rp_post",    "|t| (RP reco, postburn);|t| [GeV^{2}];Events",    nbabs, abmin, abmax);
  auto h_tab_b0_after   = new TH1D("h_tab_b0_after",   "|t| (B0 reco, afterburner);|t| [GeV^{2}];Events", nbabs, abmin, abmax);
  auto h_tab_b0_post    = new TH1D("h_tab_b0_post",    "|t| (B0 reco, postburn);|t| [GeV^{2}];Events",    nbabs, abmin, abmax);

  // 2D correlations (truth vs reco)
  const int nb2x = 60, nb2y = 60;
  auto h2_rp_after = new TH2D("h2_rp_after", "RP: |t|_{truth(after)} vs |t|_{reco(after)};|t|^{truth} [GeV^{2}];|t|^{reco} [GeV^{2}]",
                              nb2x, abmin, abmax, nb2y, abmin, abmax);
  auto h2_rp_post  = new TH2D("h2_rp_post",  "RP: |t|_{truth(post)} vs |t|_{reco(post)};|t|^{truth} [GeV^{2}];|t|^{reco} [GeV^{2}]",
                              nb2x, abmin, abmax, nb2y, abmin, abmax);
  auto h2_b0_after = new TH2D("h2_b0_after", "B0: |t|_{truth(after)} vs |t|_{reco(after)};|t|^{truth} [GeV^{2}];|t|^{reco} [GeV^{2}]",
                              nb2x, abmin, abmax, nb2y, abmin, abmax);
  auto h2_b0_post  = new TH2D("h2_b0_post",  "B0: |t|_{truth(post)} vs |t|_{reco(post)};|t|^{truth} [GeV^{2}];|t|^{reco} [GeV^{2}]",
                              nb2x, abmin, abmax, nb2y, abmin, abmax);

  // loop
  Long64_t n = 0;
  while(r.Next()){
    ++n;
    // truth
    if(isfinitef(*t_true_after))    h_t_true_after->Fill(*t_true_after);
    if(isfinitef(*t_true_post))     h_t_true_post->Fill(*t_true_post);
    if(isfinitef(*t_abs_true_after))h_tab_true_after->Fill(*t_abs_true_after);
    if(isfinitef(*t_abs_true_post)) h_tab_true_post->Fill(*t_abs_true_post);

    // RP reco
    if(*has_rp_reco){
      if(isfinitef(*t_rp_after))      h_t_rp_after->Fill(*t_rp_after);
      if(isfinitef(*t_rp_post))       h_t_rp_post->Fill(*t_rp_post);
      if(isfinitef(*t_abs_rp_after))  h_tab_rp_after->Fill(*t_abs_rp_after);
      if(isfinitef(*t_abs_rp_post))   h_tab_rp_post->Fill(*t_abs_rp_post);

      if(isfinitef(*t_abs_true_after) && isfinitef(*t_abs_rp_after))
        h2_rp_after->Fill(*t_abs_true_after, *t_abs_rp_after);
      if(isfinitef(*t_abs_true_post) && isfinitef(*t_abs_rp_post))
        h2_rp_post->Fill(*t_abs_true_post, *t_abs_rp_post);
    }

    // B0 reco
    if(*has_b0_reco){
      if(isfinitef(*t_b0_after))      h_t_b0_after->Fill(*t_b0_after);
      if(isfinitef(*t_b0_post))       h_t_b0_post->Fill(*t_b0_post);
      if(isfinitef(*t_abs_b0_after))  h_tab_b0_after->Fill(*t_abs_b0_after);
      if(isfinitef(*t_abs_b0_post))   h_tab_b0_post->Fill(*t_abs_b0_post);

      if(isfinitef(*t_abs_true_after) && isfinitef(*t_abs_b0_after))
        h2_b0_after->Fill(*t_abs_true_after, *t_abs_b0_after);
      if(isfinitef(*t_abs_true_post) && isfinitef(*t_abs_b0_post))
        h2_b0_post->Fill(*t_abs_true_post, *t_abs_b0_post);
    }
  }

  // output ROOT file
  TFile fout(outroot.c_str(), "RECREATE");
  h_t_true_after->Write();   h_t_true_post->Write();
  h_t_rp_after->Write();     h_t_rp_post->Write();
  h_t_b0_after->Write();     h_t_b0_post->Write();
  h_tab_true_after->Write(); h_tab_true_post->Write();
  h_tab_rp_after->Write();   h_tab_rp_post->Write();
  h_tab_b0_after->Write();   h_tab_b0_post->Write();
  h2_rp_after->Write();      h2_rp_post->Write();
  h2_b0_after->Write();      h2_b0_post->Write();
  fout.Close();

  // PDF (multi-page)
  TCanvas c("c","c",900,700);
  c.Print((outpdf+"[").c_str());

  // 1D signed t (truth)
  c.Clear(); h_t_true_after->SetLineWidth(2); h_t_true_after->Draw("HIST");
  h_t_true_post->SetLineWidth(2); h_t_true_post->Draw("HIST SAME");
  {
    auto leg=new TLegend(0.58,0.74,0.88,0.88);
    leg->AddEntry(h_t_true_after,"truth (after)","l");
    leg->AddEntry(h_t_true_post, "truth (post)","l");
    leg->Draw();
  }
  c.Print(outpdf.c_str());

  // 1D signed t (RP/B0)
  c.Clear(); h_t_rp_after->SetLineWidth(2); h_t_rp_after->Draw("HIST");
  h_t_rp_post->SetLineWidth(2);  h_t_rp_post->Draw("HIST SAME");
  {
    auto leg=new TLegend(0.58,0.74,0.88,0.88);
    leg->AddEntry(h_t_rp_after,"RP reco (after)","l");
    leg->AddEntry(h_t_rp_post, "RP reco (post)","l");
    leg->Draw();
  }
  c.Print(outpdf.c_str());

  c.Clear(); h_t_b0_after->SetLineWidth(2); h_t_b0_after->Draw("HIST");
  h_t_b0_post->SetLineWidth(2);  h_t_b0_post->Draw("HIST SAME");
  {
    auto leg=new TLegend(0.58,0.74,0.88,0.88);
    leg->AddEntry(h_t_b0_after,"B0 reco (after)","l");
    leg->AddEntry(h_t_b0_post, "B0 reco (post)","l");
    leg->Draw();
  }
  c.Print(outpdf.c_str());

  // 1D |t| overlays (truth vs reco)
  c.Clear(); h_tab_true_after->SetLineWidth(2); h_tab_true_after->Draw("HIST");
  h_tab_rp_after->SetLineWidth(2);   h_tab_rp_after->Draw("HIST SAME");
  h_tab_b0_after->SetLineWidth(2);   h_tab_b0_after->Draw("HIST SAME");
  {
    auto leg=new TLegend(0.58,0.70,0.88,0.88);
    leg->AddEntry(h_tab_true_after,"|t| truth (after)","l");
    leg->AddEntry(h_tab_rp_after,  "|t| RP reco (after)","l");
    leg->AddEntry(h_tab_b0_after,  "|t| B0 reco (after)","l");
    leg->Draw();
  }
  c.Print(outpdf.c_str());

  c.Clear(); h_tab_true_post->SetLineWidth(2); h_tab_true_post->Draw("HIST");
  h_tab_rp_post->SetLineWidth(2);   h_tab_rp_post->Draw("HIST SAME");
  h_tab_b0_post->SetLineWidth(2);   h_tab_b0_post->Draw("HIST SAME");
  {
    auto leg=new TLegend(0.58,0.70,0.88,0.88);
    leg->AddEntry(h_tab_true_post,"|t| truth (post)","l");
    leg->AddEntry(h_tab_rp_post,  "|t| RP reco (post)","l");
    leg->AddEntry(h_tab_b0_post,  "|t| B0 reco (post)","l");
    leg->Draw();
  }
  c.Print(outpdf.c_str());

  // 2D: truth vs reco (RP)
  c.Clear(); h2_rp_after->Draw("COLZ"); c.Print(outpdf.c_str());
  c.Clear(); h2_rp_post->Draw("COLZ");  c.Print(outpdf.c_str());

  // 2D: truth vs reco (B0)
  c.Clear(); h2_b0_after->Draw("COLZ"); c.Print(outpdf.c_str());
  c.Clear(); h2_b0_post->Draw("COLZ");  c.Print(outpdf.c_str());

  c.Print((outpdf+"]").c_str());

  std::cout << "Wrote: " << outpdf << " and " << outroot << "\n";
  return 0;
}
