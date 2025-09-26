#ifndef PLOTTING_HPP
#define PLOTTING_HPP

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <vector>
#include <utility>

// Forward declaration of TFile
class TFile;

// Base class for plotting options
class PlotOptions {
public:
    virtual ~PlotOptions();
    virtual void Plot(TFile* inputFile) = 0;
protected:
    PlotOptions() = default;
    bool m_isLogX = false;
    bool m_isLogY = false;
    std::pair<double, double> m_xRange = {-999., -999.};
    std::pair<double, double> m_yRange = {-999., -999.};
};

// Derived class for 1D histogram plotting
class PlotOptions1D : public PlotOptions {
public:
        PlotOptions1D(const std::vector<TString>& histNames,
                    const std::vector<const char*>& legendEntries,
                    const std::vector<const char*>& drawOptions,
                    const char* canvasTitle,
                    const char* xLabel,
                    const char* yLabel,
                    const char* saveName,
                    const bool isLogX = false,
                    const bool isLogY = false,
                    const bool normalizeToPDF = false
                );
                  
    void Plot(TFile* inputFile) override;

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
    // Member variables for axis ranges
    std::pair<double, double> m_xRange;
    std::pair<double, double> m_yRange;
};

// Derived class for 1D relative resolution plots with Gaussian fit
class PlotOptionsRelRes : public PlotOptions {
public:
    PlotOptionsRelRes(const TString& histName,
                      const char* xLabel,
                      const char* yLabel,
                      double xMinFit=-999.,
                      double xMaxFit=-999.,
                      const char* saveName="relres.png");
    void SetFitRangeByBins(TH1D* hist);
    void Plot(TFile* inputFile) override;
    double GetBestSymmetryValue() const { return m_bestSymmetryValue; };


private:
    TString m_histName;
    const char* m_xLabel;
    const char* m_yLabel;
    double m_xMinFit;
    double m_xMaxFit;
    const char* m_saveName;
    double m_bestMean;
    double m_bestSigma;
    double m_bestAmplitude;
    double m_bestSymmetryValue;
};

// Derived class for binned relative resolution plots
class PlotOptionsBinnedRelRes : public PlotOptions {
public:
    PlotOptionsBinnedRelRes(const TString& histName,
                            const char* title,
                            const char* xLabel,
                            const char* yLabel,
                            const std::vector<std::pair<double, double>>& fitRanges,
                            const char* saveName,
                            const char* binSavePrefix);
    void SetFitRangeByRMS(TH1D* hist);
    void SetFitRangeByBins(TH1D* hist);
    void Plot(TFile* inputFile) override;

private:
    TString m_histName;
    const char* m_title;
    const char* m_xLabel;
    const char* m_yLabel;
    std::vector<std::pair<double, double>> m_fitRanges;
    double m_xMinFit;  // Current bin's fit range min
    double m_xMaxFit;  // Current bin's fit range max
    const char* m_saveName;
    const char* m_binSavePrefix;
};

// Derived class for response matrix plotting
class PlotOptionsResponseMatrix : public PlotOptions {
public:
    PlotOptionsResponseMatrix(const TString& histName,
                              const char* xLabel,
                              const char* yLabel,
                              const char* saveName,
                              const bool  isLogX=true,
                              const bool  isLogY=true,
                              const std::pair<double, double>& xRange = {-999., -999.},
                              const std::pair<double, double>& yRange = {-999., -999.}
                             );
    
    void Plot(TFile* inputFile) override;

private:
    TString m_histName;
    const char* m_xLabel;
    const char* m_yLabel;
    const char* m_saveName;
    // bool m_isLogX;
    // bool m_isLogY;
};

#endif // PLOTTING_HPP