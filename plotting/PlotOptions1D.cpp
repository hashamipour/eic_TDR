#include "Plotting.hpp"
#include <TF1.h>
#include <TLatex.h>
#include <iostream>

PlotOptions1D::PlotOptions1D(const std::vector<TString>& histNames,
                             const std::vector<const char*>& legendEntries,
                             const std::vector<const char*>& drawOptions,
                             const char* canvasTitle,
                             const char* xLabel,
                             const char* yLabel,
                             const char* saveName,
                             const bool  isLogX,
                             const bool  isLogY,
                             const bool  normalizeToPDF)
    : m_histNames(histNames),
      m_legendEntries(legendEntries),
      m_drawOptions(drawOptions),
      m_canvasTitle(canvasTitle),
      m_xLabel(xLabel),
      m_yLabel(yLabel),
      m_saveName(saveName),
      m_isLogX(isLogX),
      m_isLogY(isLogY),
      m_normalizeToPDF(normalizeToPDF) {}

void PlotOptions1D::Plot(TFile* inputFile) {
    TCanvas* c = new TCanvas("my_canvas", m_canvasTitle, 1200, 900);
    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    TH1* firstHist = nullptr;

    for (size_t i = 0; i < m_histNames.size(); ++i) {
        TH1* hist = (TH1*)inputFile->Get(m_histNames[i]);
        if (!hist) {
            std::cerr << "Error: Histogram " << m_histNames[i] << " not found." << std::endl;
            continue;
        }

        if (m_normalizeToPDF) {
            hist->Sumw2();
            double integral_w = hist->Integral("width");
            if (integral_w > 0) {
                hist->Scale(1.0 / integral_w);
            }
        }

        hist->SetLineWidth((i == 0) ? 2 : 1);
        if (m_histNames[i].Contains("truth")) {
            hist->SetLineColor(kBlack);
        } else if (m_histNames[i].Contains("EM")) {
            hist->SetLineColor(kRed);
            hist->SetMarkerColor(kRed);
            hist->SetMarkerStyle(20);
        } else if (m_histNames[i].Contains("DA")) {
            hist->SetLineColor(kBlue);
            hist->SetMarkerColor(kBlue);
            hist->SetMarkerStyle(20);
        } else if (m_histNames[i].Contains("ESigma")) {
            hist->SetLineColor(kGreen + 2);
            hist->SetMarkerColor(kGreen + 2);
            hist->SetMarkerStyle(20);
        }

        if (i == 0) {
            TString drawOption = m_drawOptions[i];
            hist->Draw(drawOption);
            firstHist = hist;
        } else {
            TString drawOption = m_drawOptions[i];
            hist->Draw(drawOption + "SAME");
        }
        legend->AddEntry(hist, m_legendEntries[i], m_drawOptions[i]);
    }

    if (firstHist) {
        firstHist->GetXaxis()->SetTitle(m_xLabel);
        firstHist->GetYaxis()->SetTitle(m_yLabel);
        firstHist->GetXaxis()->SetTitleSize(0.04);
        firstHist->GetYaxis()->SetTitleSize(0.04);
        firstHist->GetXaxis()->SetTitleOffset(1.0);
        firstHist->GetYaxis()->SetTitleOffset(1.0);
    }

    if (m_isLogX) c->SetLogx();
    if (m_isLogY) c->SetLogy();

    legend->Draw();
    c->Update();
    c->SaveAs(m_saveName);
    delete c;
}
