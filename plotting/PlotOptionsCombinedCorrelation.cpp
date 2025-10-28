#include "Plotting.hpp"
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <iostream>

PlotOptionsCombinedCorrelation::PlotOptionsCombinedCorrelation(const std::vector<TString>& histNames,
                                                               const std::vector<const char*>& legendEntries,
                                                               const std::vector<Color_t>& colors,
                                                               const std::vector<Style_t>& markerStyles,
                                                               const char* canvasTitle,
                                                               const char* xLabel,
                                                               const char* yLabel,
                                                               const char* saveName,
                                                               const std::pair<double, double>& xRange,
                                                               const std::pair<double, double>& yRange)
    : m_histNames(histNames),
      m_legendEntries(legendEntries),
      m_colors(colors),
      m_markerStyles(markerStyles),
      m_canvasTitle(canvasTitle),
      m_xLabel(xLabel),
      m_yLabel(yLabel),
      m_saveName(saveName),
      m_xRange(xRange),
      m_yRange(yRange) {}

void PlotOptionsCombinedCorrelation::Plot(TFile* inputFile) {
    TCanvas* c = new TCanvas("c_combined_corr", m_canvasTitle, 1200, 1000);
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);
    c->SetTopMargin(0.1);
    c->SetBottomMargin(0.1);
    
    gStyle->SetOptStat(0);
    
    bool firstPlot = true;
    
    for (size_t i = 0; i < m_histNames.size(); ++i) {
        TGraph* graph = (TGraph*)inputFile->Get(m_histNames[i]);
        
        if (!graph) {
            std::cerr << "Error: TGraph " << m_histNames[i] << " not found." << std::endl;
            continue;
        }
        
        graph->SetMarkerColor(m_colors[i]);
        graph->SetMarkerStyle(m_markerStyles[i]);
        graph->SetMarkerSize(0.3);
        
        if (firstPlot) {
            graph->SetTitle("");
            graph->GetXaxis()->SetTitle(m_xLabel);
            graph->GetYaxis()->SetTitle(m_yLabel);
            graph->GetXaxis()->SetTitleOffset(1.3);
            graph->GetYaxis()->SetTitleOffset(1.3);
            
            if (m_xRange.first != -999. && m_xRange.second != -999.) {
                graph->GetXaxis()->SetLimits(m_xRange.first, m_xRange.second);
            }
            if (m_yRange.first != -999. && m_yRange.second != -999.) {
                graph->GetHistogram()->SetMinimum(m_yRange.first);
                graph->GetHistogram()->SetMaximum(m_yRange.second);
            }
            
            graph->Draw("AP");
            firstPlot = false;
        } else {
            graph->Draw("P SAME");
        }
    }
    
    double xmin = (m_xRange.first != -999.) ? m_xRange.first : 0.0;
    double xmax = (m_xRange.second != -999.) ? m_xRange.second : 2.0;
    double ymin = (m_yRange.first != -999.) ? m_yRange.first : 0.0;
    double ymax = (m_yRange.second != -999.) ? m_yRange.second : 2.0;
    
    double line_min = TMath::Max(xmin, ymin);
    double line_max = TMath::Min(xmax, ymax);
    
    TLine* diagLine = new TLine(line_min, line_min, line_max, line_max);
    diagLine->SetLineColor(kRed);
    diagLine->SetLineWidth(2);
    diagLine->SetLineStyle(2);
    diagLine->Draw("same");
    
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.17, 0.92, "#bf{ePIC} Simulation (100k events)");
    latex.DrawLatex(0.65, 0.92, "#bf{Diff. DIS} 10x100 GeV");
    
    c->Update();
    c->SaveAs(m_saveName);
    
    delete diagLine;
    delete c;
}