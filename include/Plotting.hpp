#ifndef PLOTTING_HPP
#define PLOTTING_HPP

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <vector>
#include <TStyle.h>
#include <TColor.h>
#include "Utility.hpp"  // Use your existing utility functions

// Base class for plot options
class PlotOptions {
public:
    virtual ~PlotOptions();
    virtual void Plot(TFile* inputFile) = 0;
};

// 1D histogram plot options
class PlotOptions1D : public PlotOptions {
private:
    std::vector<TString> m_histNames;
    std::vector<const char*> m_legendEntries;
    std::vector<const char*> m_drawOptions;
    const char* m_canvasTitle;
    const char* m_xLabel;
    const char* m_yLabel;
    const char* m_saveName;
    bool m_isLogX;
    bool m_isLogY;
    bool m_normalizeToPDF;

public:
    PlotOptions1D(const std::vector<TString>& histNames,
                  const std::vector<const char*>& legendEntries,
                  const std::vector<const char*>& drawOptions,
                  const char* canvasTitle,
                  const char* xLabel,
                  const char* yLabel,
                  const char* saveName,
                  const bool  isLogX = false,
                  const bool  isLogY = false,
                  const bool  normalizeToPDF = false);
    
    void Plot(TFile* inputFile) override;
};

// Relative resolution plot options
class PlotOptionsRelRes : public PlotOptions {
private:
    TString m_histName;
    const char* m_xLabel;
    const char* m_yLabel;
    double m_xMinFit;
    double m_xMaxFit;
    const char* m_saveName;
    // optimization variables for fit
    double m_bestMean;
    double m_bestSigma;
    double m_bestAmplitude;

    void SetFitRangeByBins(TH1D* hist);

public:
    PlotOptionsRelRes(const TString& histName,
                      const char* xLabel,
                      const char* yLabel,
                      double xMinFit,
                      double xMaxFit,
                      const char* saveName);
    
    void Plot(TFile* inputFile) override;
};

// Binned relative resolution plot options
class PlotOptionsBinnedRelRes : public PlotOptions {
private:
    TString m_histName;
    const char* m_title;
    const char* m_xLabel;
    const char* m_yLabel;
    std::vector<std::pair<double, double>> m_fitRanges;
    double m_xMinFit;
    double m_xMaxFit;
    const char* m_saveName;
    const char* m_binSavePrefix;

    void SetFitRangeByBins(TH1D* hist);

public:
    PlotOptionsBinnedRelRes(const TString& histName,
                            const char* title,
                            const char* xLabel,
                            const char* yLabel,
                            const std::vector<std::pair<double, double>>& fitRanges,
                            const char* saveName,
                            const char* binSavePrefix);
    
    void Plot(TFile* inputFile) override;
};

// Response matrix plot options
class PlotOptionsResponseMatrix : public PlotOptions {
private:
    TString m_histName;
    const char* m_xLabel;
    const char* m_yLabel;
    const char* m_saveName;
    bool m_isLogX;
    bool m_isLogY;
    std::pair<double, double> m_xRange;
    std::pair<double, double> m_yRange;

public:
    PlotOptionsResponseMatrix(const TString& histName,
                              const char* xLabel,
                              const char* yLabel,
                              const char* saveName,
                              const bool isLogX = false,
                              const bool isLogY = false,
                              const std::pair<double, double>& xRange = {-999., -999.},
                              const std::pair<double, double>& yRange = {-999., -999.});
    
    void Plot(TFile* inputFile) override;
};

// NEW CLASS: Combined correlation plot options
class PlotOptionsCombinedCorrelation : public PlotOptions {
private:
    std::vector<TString> m_histNames;
    std::vector<const char*> m_legendEntries;
    std::vector<Color_t> m_colors;
    std::vector<Style_t> m_markerStyles;
    const char* m_canvasTitle;
    const char* m_xLabel;
    const char* m_yLabel;
    const char* m_saveName;
    std::pair<double, double> m_xRange;
    std::pair<double, double> m_yRange;

public:
    PlotOptionsCombinedCorrelation(const std::vector<TString>& histNames,
                                   const std::vector<const char*>& legendEntries,
                                   const std::vector<Color_t>& colors,
                                   const std::vector<Style_t>& markerStyles,
                                   const char* canvasTitle,
                                   const char* xLabel,
                                   const char* yLabel,
                                   const char* saveName,
                                   const std::pair<double, double>& xRange = {-999., -999.},
                                   const std::pair<double, double>& yRange = {-999., -999.});
    
    void Plot(TFile* inputFile) override;
};

#endif // PLOTTING_HPP