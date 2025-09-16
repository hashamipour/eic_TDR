#include "Plotting.hpp"
#include <TF1.h>
#include <TLatex.h>
#include <TPave.h>
#include <TLine.h>
#include <iostream>
#include <TStyle.h>  // for gStyle
#include <TPaveText.h>


// Function to define and set a custom color palette
void SetCustomPalette() {
    // Define the number of primary colors and the total number of color levels in the palette
    const Int_t NRGBs = 5; // We have 5 distinct colors
    const Int_t NCont = 255; // Number of smoothly interpolated colors in the palette

    // Define the hexadecimal color codes from your image
    // You can also add more colors for finer gradient control if needed
    // These are the "anchor" points for the gradient.
    // https://digitalsynopsis.com/design/beautiful-color-gradient-palettes/
    // NUMBER 11
    // we can call it "SolarSplash"
    // UInt_t colors[NRGBs] = {
    //     0x35BBCA, // color 1
    //     0x0191B4, // color 2
    //     0xF8D90F, // color 3
    //     0xD3DD18, // color 4
    //     0xFE7A15  // color 5
    //     // 0xD46C4E  // color 5, a bit softer than FE7A15
    // };


    // another one I found online; https://www.schemecolor.com/best-gradient.php
    // we can call it "CandyEclipse"
    //     UInt_t colors[NRGBs] = {
    //     0xFEE27A, // color 1
    //     0xFEA959, // color 2
    //     0xFF605D, // color 3
    //     0xF13484, // color 4
    //     0x8D379E  // color 5
    // };


    // yet another one, from https://www.schemecolor.com/yellow-green-gradient-color-scheme.php
    // we can call it "LemonLush"
    //     UInt_t colors[NRGBs] = {
    //     0xFEFE69, // color 1
    //     0xDDF969, // color 2
    //     0xA9F36A, // color 3
    //     0x78EC6C, // color 4
    //     0x57E86B  // color 5
    // };

    // another choice, from https://www.schemecolor.com/green-yellow-gradient-color-scheme.php
    // with modifications to make it more green
    // we can call it "CitrusFade"
        // UInt_t colors[NRGBs] = {
        //     0xFFFF00, // Fluorescent Yellow
        //     0xC6F94D, // Citrus Lime
        //     0x7ED957, // Leafy Green
        //     0x4CCB5B, // Spring Green
        //     0x249E57  // Lush Mid-Green
        // };


    // a three color gradient inspired by the former
    // we can call it "PetalFlare"
UInt_t colors[NRGBs] = {
    0xE85A78, // Rosy Crimson (with pink undertone)
    0xF57C3A, // Early Orange (10% mark)
    0xFFE94F, // Bright Yellow
    0xA4D65E, // Zesty Green
    0x2FAE5E  // Slightly Darker Green
};

    // Convert hex colors to RGB values (ROOT expects RGB from 0.0 to 1.0)
    // TColor::SetRGB takes unsigned char (0-255) for R, G, B
    Double_t stops[NRGBs]; // Position of each color in the gradient (0.0 to 1.0)
    Double_t red[NRGBs], green[NRGBs], blue[NRGBs];

    for (int i = 0; i < NRGBs; ++i) {
        // Calculate stop position (evenly distributed for now)
        stops[i] = (Double_t)i / (NRGBs - 1.0);

        // Extract RGB components from hex
        UChar_t r = (colors[i] >> 16) & 0xFF;
        UChar_t g = (colors[i] >> 8) & 0xFF;
        UChar_t b = colors[i] & 0xFF;

        // Convert to Double_t for TColor::CreateGradientColorTable
        red[i] = (Double_t)r / 255.0;
        green[i] = (Double_t)g / 255.0;
        blue[i] = (Double_t)b / 255.0;
    }

    // Create the gradient color table
    // This defines the smooth transition between the specified colors
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

    // Set the number of contours (color levels) for the palette
    gStyle->SetNumberContours(NCont);

    // You might want to call gStyle->SetPalette(NCont) as well,
    // but CreateGradientColorTable often implicitly sets the palette.
    // If not, explicitly set it:
    // gStyle->SetPalette(NCont, 0); // 0 means use the current gradient table
}


// PlotOptions destructor
PlotOptions::~PlotOptions() {}

// PlotOptions1D implementation
PlotOptions1D::PlotOptions1D(const std::vector<TString>& histNames,
                             const std::vector<const char*>& legendEntries,
                             const std::vector<const char*>& drawOptions,
                             const char* canvasTitle,
                             const char* xLabel,
                             const char* yLabel,
                             const char* saveName,
                             const bool isLogX,
                             const bool isLogY,
                             const bool normalizeToPDF)
    : m_histNames(histNames),
      m_legendEntries(legendEntries),
      m_drawOptions(drawOptions),
      m_canvasTitle(canvasTitle),
      m_xLabel(xLabel),
      m_yLabel(yLabel),
      m_saveName(saveName)  ,
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
        hist->Sumw2(); // ensure errors are tracked
        // compute integral taking bin width into account
        double integral_w = hist->Integral("width"); // sum(c_i * width_i)
        if (integral_w > 0) {
            hist->Scale(1.0 / integral_w); // now histogram integrates to 1 (PDF)
        }
    }

        hist->SetLineWidth( (i == 0) ? 2 : 1);
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
        hist->Draw(drawOption);       // draw first with its option (no "SAME")
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


/////////////////////////////////////////////////////////////////
// PlotOptionsRelRes implementation
PlotOptionsRelRes::PlotOptionsRelRes(const TString& histName,
                                     const char* xLabel,
                                     const char* yLabel,
                                     double xMinFit,
                                     double xMaxFit,
                                     const char* saveName)
    : m_histName(histName),
      m_xLabel(xLabel),
      m_yLabel(yLabel),
      m_xMinFit(xMinFit),
      m_xMaxFit(xMaxFit),
      m_saveName(saveName),
      m_bestMean(0),
      m_bestSigma(0),
      m_bestAmplitude(0) {}

void PlotOptionsRelRes::SetFitRangeByBins(TH1D* hist) {
    // Find the peak bin and its center
    int peakBin = hist->GetMaximumBin();
    double peakCenter = hist->GetBinCenter(peakBin);

    // Structure to store fit results
    struct FitResult {
        int nBinsLeft;
        int nBinsRight;
        int totalBins;
        double chi2_ndf;
        double yValueDiff;
        double yValueDiffPercent;
        double y_left;
        double y_right;
        double amplitude;
        double mean;
        double sigma;
        double xMinFit;
        double xMaxFit;
    };

    std::vector<FitResult> validFits;

    // Print header for the list of fits
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "   nBins_L | nBins_R | Total | Chi2/NDF | Y_L | Y_R | |f(x_L) - f(x_R)|% " << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    
    // Loop through different numbers of bins
    for (int nBinsLeft = 3; nBinsLeft <= 15; ++nBinsLeft) {
        for (int nBinsRight = 3; nBinsRight <= 15; ++nBinsRight) {
            
            int minBinFit = peakBin - nBinsLeft;
            int maxBinFit = peakBin + nBinsRight;
            
            if (minBinFit < 1 || maxBinFit > hist->GetNbinsX()) {
                continue;
            }
            
            double x_min_fit = hist->GetBinLowEdge(minBinFit);
            double x_max_fit = hist->GetBinLowEdge(maxBinFit + 1);
            
            // Get bin centers for y-value evaluation
            double x_left_bin_center = hist->GetBinCenter(minBinFit);
            double x_right_bin_center = hist->GetBinCenter(maxBinFit);
            
            TF1* currentGaus = new TF1("currentGaus", "gaus", x_min_fit, x_max_fit);
            currentGaus->SetParameters(hist->GetMaximum(), hist->GetBinCenter(peakBin), hist->GetRMS());
            
            hist->Fit(currentGaus, "RQN", "", x_min_fit, x_max_fit);
            
            double currentChi2 = currentGaus->GetChisquare();
            int ndf = currentGaus->GetNDF();
            
            if (ndf > 0) {
                double chi2_ndf = currentChi2 / ndf;
                
                // Skip fits with chi2/ndf > 10
                if (chi2_ndf > 10.0) {
                    delete currentGaus;
                    continue;
                }
                
                // Evaluate the fitted function at the left and right bin centers
                double y_left = currentGaus->Eval(x_left_bin_center);
                double y_right = currentGaus->Eval(x_right_bin_center);
                double yValueDiff = TMath::Abs(y_left - y_right);
                double y_max = TMath::Max(y_left, y_right);
                double yValueDiffPercent = (y_max > 0) ? (yValueDiff / y_max * 100.0) : 0.0;
                
                // Skip fits where y-value difference is > 10% of the larger y-value
                if (yValueDiffPercent > 10.0) {
                    delete currentGaus;
                    continue;
                }
                
                // Store the fit result
                FitResult result;
                result.nBinsLeft = nBinsLeft;
                result.nBinsRight = nBinsRight;
                result.totalBins = nBinsLeft + nBinsRight;
                result.chi2_ndf = chi2_ndf;
                result.yValueDiff = yValueDiff;
                result.yValueDiffPercent = yValueDiffPercent;
                result.y_left = y_left;
                result.y_right = y_right;
                result.amplitude = currentGaus->GetParameter(0);
                result.mean = currentGaus->GetParameter(1);
                result.sigma = currentGaus->GetParameter(2);
                result.xMinFit = x_min_fit;
                result.xMaxFit = x_max_fit;
                
                validFits.push_back(result);
                
                // Print the current fit's details
                std::cout << TString::Format("   %7d | %7d | %5d | %8.4f | %6.1f | %6.1f | %8.2f%%", 
                    nBinsLeft, nBinsRight, result.totalBins, chi2_ndf, y_left, y_right, yValueDiffPercent) << std::endl;
            }
            
            delete currentGaus;
        }
    }
    
    if (validFits.empty()) {
        std::cerr << "Error: No valid fits found! (All fits had chi2/ndf > 10 or y-difference > 10%)" << std::endl;
        return;
    }
    
    // Sort fits by: y-symmetry (ascending), then interval size (descending), then chi2 (ascending)
    std::sort(validFits.begin(), validFits.end(), [](const FitResult& a, const FitResult& b) {
        if (TMath::Abs(a.yValueDiffPercent - b.yValueDiffPercent) > 0.01) 
            return a.yValueDiffPercent < b.yValueDiffPercent;  // Better y-symmetry first
        if (a.totalBins != b.totalBins) 
            return a.totalBins > b.totalBins;                 // Larger interval second
        return a.chi2_ndf < b.chi2_ndf;                       // Better chi2 third
    });
    
    // Get the best y-symmetry value
    double bestYSymmetry = validFits[0].yValueDiffPercent;
    
    // Get all fits with the best y-symmetry (within 0.1% tolerance)
    std::vector<FitResult> bestYSymmetryFits;
    for (const auto& fit : validFits) {
        if (fit.yValueDiffPercent <= bestYSymmetry + 0.1) {
            bestYSymmetryFits.push_back(fit);
        }
    }
    
    // Among best y-symmetry fits, find maximum interval size
    int maxSizeAmongBestY = 0;
    for (const auto& fit : bestYSymmetryFits) {
        if (fit.totalBins > maxSizeAmongBestY) {
            maxSizeAmongBestY = fit.totalBins;
        }
    }
    
    // Get fits with best y-symmetry AND maximum size
    std::vector<FitResult> bestYAndSizeFits;
    for (const auto& fit : bestYSymmetryFits) {
        if (fit.totalBins == maxSizeAmongBestY) {
            bestYAndSizeFits.push_back(fit);
        }
    }
    
    // Among those, pick the one with best chi2
    FitResult bestFit = *std::min_element(bestYAndSizeFits.begin(), bestYAndSizeFits.end(),
        [](const FitResult& a, const FitResult& b) {
            return a.chi2_ndf < b.chi2_ndf;
        });
    
    // Store the best fit parameters
    m_xMinFit = bestFit.xMinFit;
    m_xMaxFit = bestFit.xMaxFit;
    m_bestAmplitude = bestFit.amplitude;
    m_bestMean = bestFit.mean;
    m_bestSigma = bestFit.sigma;
    
    // Calculate the final symmetry value for display (keeping original calculation for reference)
    double finalXSymmetryValue = TMath::Abs(m_xMinFit + m_xMaxFit - 2.0 * peakCenter);
    
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "SELECTION SUMMARY:" << std::endl;
    std::cout << "Total valid fits (chi2<10, y-diff<10%): " << validFits.size() << std::endl;
    std::cout << "Best y-symmetry: " << bestYSymmetry << "%" << std::endl;
    std::cout << "Fits with best y-symmetry: " << bestYSymmetryFits.size() << std::endl;
    std::cout << "Max interval size among best y-symmetry: " << maxSizeAmongBestY << " bins" << std::endl;
    std::cout << "Fits with best y-symmetry AND max size: " << bestYAndSizeFits.size() << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "OPTIMAL FIT SELECTED:" << std::endl;
    std::cout << "Bins: " << bestFit.nBinsLeft << " left + " << bestFit.nBinsRight << " right = " 
              << bestFit.totalBins << " total" << std::endl;
    std::cout << "Fit range: [" << m_xMinFit << ", " << m_xMaxFit << "]" << std::endl;
    std::cout << "Chi2/NDF: " << bestFit.chi2_ndf << std::endl;
    std::cout << "Y-values: f(x_L)=" << bestFit.y_left << ", f(x_R)=" << bestFit.y_right << std::endl;
    std::cout << "Y-symmetry: " << bestFit.yValueDiffPercent << "%" << std::endl;
    std::cout << "X-symmetry value (original): " << finalXSymmetryValue << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
}

// Plot method
void PlotOptionsRelRes::Plot(TFile* inputFile) {
    TCanvas* c = new TCanvas("c_relres", "", 1200, 800);
    c->SetBottomMargin(0.2);

    TH1D* hist = (TH1D*)inputFile->Get(m_histName);
    if (!hist) {
        std::cerr << "Error: Histogram " << m_histName << " not found." << std::endl;
        delete c;
        return;
    }

    hist->GetXaxis()->SetTitle(m_xLabel);
    hist->GetYaxis()->SetTitle(m_yLabel);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleOffset(0.9);
    hist->SetLineColor(kGreen + 3);
    hist->SetMarkerColor(kGreen + 3);
    hist->SetMarkerStyle(20);
    hist->Draw("pe");
    
    // Use bins to set the best fit range and get optimized parameters
    SetFitRangeByBins(hist);
    
    std::cout << "peak at: " << hist->GetBinCenter(hist->GetMaximumBin()) << std::endl;
    std::cout << "Chosen fit range: [" << m_xMinFit << ", " << m_xMaxFit << "]" << std::endl;
    
    // Create the final fit function with optimized parameters as starting values
    TF1* gaussianFit = new TF1("gaussianFit", "gaus", m_xMinFit, m_xMaxFit);
    gaussianFit->SetParameters(m_bestAmplitude, m_bestMean, m_bestSigma);
    
    std::cout << "Using optimized parameters as starting values: mean=" << m_bestMean 
              << ", sigma=" << m_bestSigma << ", amplitude=" << m_bestAmplitude << std::endl;
    
    // Fit without parameter limits
    hist->Fit("gaussianFit", "RQ+");

    gaussianFit->SetLineColor(kRed);
    gaussianFit->SetLineWidth(2);
    hist->GetYaxis()->SetRangeUser(0, gaussianFit->GetMaximum() * 1.1);
    gaussianFit->Draw("same");

    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->SetTextColor(kRed);
    latex->DrawLatex(0.2, 0.8, Form("RMS = %.2f", hist->GetRMS()));

    c->SaveAs(m_saveName);
    delete c;
    delete latex;
}

/////////////////////////////////////////////////////////////////

// PlotOptionsBinnedRelRes implementation
PlotOptionsBinnedRelRes::PlotOptionsBinnedRelRes(const TString& histName,
                                                 const char* title,
                                                 const char* xLabel,
                                                 const char* yLabel,
                                                 const std::vector<std::pair<double, double>>& fitRanges,
                                                 const char* saveName,
                                                 const char* binSavePrefix)
    : m_histName(histName),
      m_title(title),
      m_xLabel(xLabel),
      m_yLabel(yLabel),
      m_fitRanges(fitRanges),
      m_xMinFit(0.0),
      m_xMaxFit(0.0),
      m_saveName(saveName),
      m_binSavePrefix(binSavePrefix) {}


void PlotOptionsBinnedRelRes::SetFitRangeByBins(TH1D* hist) {
    // Find the peak bin and its center
    int peakBin = hist->GetMaximumBin();
    double peakCenter = hist->GetBinCenter(peakBin);

    // Structure to store fit results
    struct FitResult {
        int nBinsLeft;
        int nBinsRight;
        int totalBins;
        double chi2_ndf;
        double yValueDiffPercent;
        double xMinFit;
        double xMaxFit;
    };

    std::vector<FitResult> validFits;

    // Print header for the list of fits
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "   nBins_L | nBins_R | Total | Chi2/NDF | |f(x_L) - f(x_R)|% " << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    
    // Loop through different numbers of bins (removed x-symmetry constraint)
    for (int nBinsLeft = 2; nBinsLeft <= 3; ++nBinsLeft) {
        for (int nBinsRight = 2; nBinsRight <= 2; ++nBinsRight) {
            
            int minBinFit = peakBin - nBinsLeft;
            int maxBinFit = peakBin + nBinsRight;
            
            if (minBinFit < 1 || maxBinFit > hist->GetNbinsX()) {
                continue;
            }
            
            double x_min_fit = hist->GetBinLowEdge(minBinFit);
            double x_max_fit = hist->GetBinLowEdge(maxBinFit + 1);
            
            // Get bin centers for y-value evaluation
            double x_left_bin_center = hist->GetBinCenter(minBinFit);
            double x_right_bin_center = hist->GetBinCenter(maxBinFit);
            
            TF1* currentGaus = new TF1("currentGaus", "gaus", x_min_fit, x_max_fit);
            currentGaus->SetParameters(hist->GetMaximum(), hist->GetBinCenter(peakBin), hist->GetRMS());
            
            hist->Fit(currentGaus, "RQN", "", x_min_fit, x_max_fit);
            
            double currentChi2 = currentGaus->GetChisquare();
            int ndf = currentGaus->GetNDF();
            
            if (ndf > 0) {
                double chi2_ndf = currentChi2 / ndf;
                
                // Skip fits with chi2/ndf > 10
                if (chi2_ndf > 10.0) {
                    delete currentGaus;
                    continue;
                }
                
                // Evaluate the fitted function at the left and right bin centers
                double y_left = currentGaus->Eval(x_left_bin_center);
                double y_right = currentGaus->Eval(x_right_bin_center);
                double yValueDiff = TMath::Abs(y_left - y_right);
                double y_max = TMath::Max(y_left, y_right);
                double yValueDiffPercent = (y_max > 0) ? (yValueDiff / y_max * 100.0) : 0.0;
                
                // Skip fits where y-value difference is > 10% of the larger y-value
                if (yValueDiffPercent > 10.0) {
                    delete currentGaus;
                    continue;
                }
                
                // Store the fit result
                FitResult result;
                result.nBinsLeft = nBinsLeft;
                result.nBinsRight = nBinsRight;
                result.totalBins = nBinsLeft + nBinsRight;
                result.chi2_ndf = chi2_ndf;
                result.yValueDiffPercent = yValueDiffPercent;
                result.xMinFit = x_min_fit;
                result.xMaxFit = x_max_fit;
                
                validFits.push_back(result);
                
                // Print the current fit's details
                std::cout << TString::Format("   %7d | %7d | %5d | %8.4f | %8.2f%%", 
                    nBinsLeft, nBinsRight, result.totalBins, chi2_ndf, yValueDiffPercent) << std::endl;
            }
            
            delete currentGaus;
        }
    }
    
    if (validFits.empty()) {
        std::cerr << "Error: No valid fits found! (All fits had chi2/ndf > 10 or y-difference > 10%)" << std::endl;
        // Fallback to a default range if no good fits are found
        m_xMinFit = hist->GetBinCenter(peakBin) - 0.55 *  hist->GetRMS();
        m_xMaxFit = hist->GetBinCenter(peakBin) + 0.45 *  hist->GetRMS(); // empirical asymmetric range
        return;
    }
    
    // Sort fits by: y-symmetry (ascending), then interval size (descending), then chi2 (ascending)
    std::sort(validFits.begin(), validFits.end(), [](const FitResult& a, const FitResult& b) {
        if (TMath::Abs(a.yValueDiffPercent - b.yValueDiffPercent) > 0.01) 
            return a.yValueDiffPercent < b.yValueDiffPercent;  // Better y-symmetry first
        if (a.totalBins != b.totalBins) 
            return a.totalBins > b.totalBins;                 // Larger interval second
        return a.chi2_ndf < b.chi2_ndf;                       // Better chi2 third
    });
    
    // Get the best y-symmetry value
    double bestYSymmetry = validFits[0].yValueDiffPercent;
    
    // Get all fits with the best y-symmetry (within 0.1% tolerance)
    std::vector<FitResult> bestYSymmetryFits;
    for (const auto& fit : validFits) {
        if (fit.yValueDiffPercent <= bestYSymmetry + 0.1) {
            bestYSymmetryFits.push_back(fit);
        }
    }
    
    // Among best y-symmetry fits, find maximum interval size
    int maxSizeAmongBestY = 0;
    for (const auto& fit : bestYSymmetryFits) {
        if (fit.totalBins > maxSizeAmongBestY) {
            maxSizeAmongBestY = fit.totalBins;
        }
    }
    
    // Get fits with best y-symmetry AND maximum size
    std::vector<FitResult> bestYAndSizeFits;
    for (const auto& fit : bestYSymmetryFits) {
        if (fit.totalBins == maxSizeAmongBestY) {
            bestYAndSizeFits.push_back(fit);
        }
    }
    
    // Among those, pick the one with best chi2
    FitResult bestFit = *std::min_element(bestYAndSizeFits.begin(), bestYAndSizeFits.end(),
        [](const FitResult& a, const FitResult& b) {
            return a.chi2_ndf < b.chi2_ndf;
        });
    
    // Store the best fit parameters
    m_xMinFit = bestFit.xMinFit;
    m_xMaxFit = bestFit.xMaxFit;
    
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "SELECTION SUMMARY:" << std::endl;
    std::cout << "Total valid fits (chi2<10, y-diff<10%): " << validFits.size() << std::endl;
    std::cout << "Best y-symmetry: " << bestYSymmetry << "%" << std::endl;
    std::cout << "Fits with best y-symmetry: " << bestYSymmetryFits.size() << std::endl;
    std::cout << "Max interval size among best y-symmetry: " << maxSizeAmongBestY << " bins" << std::endl;
    std::cout << "Fits with best y-symmetry AND max size: " << bestYAndSizeFits.size() << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "OPTIMAL FIT SELECTED:" << std::endl;
    std::cout << "Bins: " << bestFit.nBinsLeft << " left + " << bestFit.nBinsRight << " right = " 
              << bestFit.totalBins << " total" << std::endl;
    std::cout << "Fit range: [" << m_xMinFit << ", " << m_xMaxFit << "]" << std::endl;
    std::cout << "Chi2/NDF: " << bestFit.chi2_ndf << std::endl;
    std::cout << "Y-symmetry: " << bestFit.yValueDiffPercent << "%" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
}

void PlotOptionsBinnedRelRes::Plot(TFile* inputFile) {
    TH2D* h_RelRes_binned = (TH2D*)inputFile->Get(m_histName);
    if (!h_RelRes_binned) {
        std::cerr << "Error: 2D Histogram " << m_histName << " not found." << std::endl;
        return;
    }

    TCanvas* c = new TCanvas("c_binned", "", 1200, 800);
    c->SetLeftMargin(0.15);

    int nbinsX = h_RelRes_binned->GetNbinsX();
    TGraphErrors* g = new TGraphErrors(nbinsX);
    TGraphErrors* g_RMS = new TGraphErrors(nbinsX);

    for (int j = 1; j <= nbinsX; ++j) {
        TH1D* projY = h_RelRes_binned->ProjectionY("", j, j);
        
        if (projY->GetEntries() == 0) {
            delete projY;
            continue;
        }

        // Determine fit range: use provided ranges or automatic detection
        if (m_fitRanges.empty()) {
            // Use the automatic method to set the fit range
            SetFitRangeByBins(projY);
        } else {
            // Use provided fit ranges (bin index j-1 since j starts from 1)
            if (j - 1 < static_cast<int>(m_fitRanges.size())) {
                m_xMinFit = m_fitRanges[j - 1].first;
                m_xMaxFit = m_fitRanges[j - 1].second;
            } else {
                // If not enough ranges provided, fall back to automatic for this bin
                std::cerr << "Warning: Not enough fit ranges provided for bin " << j 
                          << ". Using automatic range detection." << std::endl;
                SetFitRangeByBins(projY);
            }
        }
        
        double _center = h_RelRes_binned->GetXaxis()->GetBinCenter(j);
        
        TF1* gaus = new TF1("gaus", "gaus", m_xMinFit, m_xMaxFit);
        gaus->SetParameters(projY->GetMaximum(), projY->GetMean(), projY->GetRMS());
        projY->Fit(gaus, "RQ");
        
        double mean = gaus->GetParameter(1);
        double sigma = gaus->GetParameter(2);

        TCanvas* c_proj = new TCanvas(Form("c_proj_%s_%d", m_binSavePrefix, j),
                                      Form("Bin Projection %d", j), 800, 600);
        projY->SetTitle(Form("%s Bin: %.1f-%.1f", m_xLabel, 
                             h_RelRes_binned->GetXaxis()->GetBinLowEdge(j), 
                             h_RelRes_binned->GetXaxis()->GetBinUpEdge(j)));
        
        // Disable automatic stats box
        projY->SetStats(0);
        projY->Draw();
        
        double ymax = 1.1 * std::fmax(projY->GetMaximum(), gaus->GetMaximum());
        projY->GetYaxis()->SetRangeUser(0, ymax);
        gaus->SetLineColor(kRed);
        gaus->Draw("same");
        
        // Create custom stats box
        TPaveText* statsBox = new TPaveText(0.6, 0.6, 0.9, 0.9, "NDC");
        statsBox->SetFillColor(kWhite);
        statsBox->SetBorderSize(1);
        statsBox->SetTextAlign(12); // Left aligned
        statsBox->SetTextSize(0.035);
        statsBox->SetTextFont(42);

        statsBox->AddText(Form("Entries: %d", (int)projY->GetEntries()));
        statsBox->AddText(Form("Mean: %.5f", projY->GetMean()));
        statsBox->AddText(Form("RMS: %.5f", projY->GetRMS()));
        statsBox->AddText(Form("#chi^{2} / ndf: %.1f / %d", gaus->GetChisquare(), gaus->GetNDF()));
        statsBox->AddText(Form("Mean: %.5f #pm %.5f", mean, gaus->GetParError(1)));
        statsBox->AddText(Form("Sigma: %.4f #pm %.5f", sigma, gaus->GetParError(2)));
        statsBox->AddText(Form("Fit range: [%.3f, %.3f]", m_xMinFit, m_xMaxFit));
        statsBox->Draw();
        
        c_proj->Update();
        c_proj->SaveAs(Form("figs/%s_bin_%d.png", m_binSavePrefix, j));
        
        c->SetLogx();
        g->SetPoint(j - 1, _center, mean);
        g->SetPointError(j - 1, 0.0, sigma);
        g_RMS->SetPoint(j - 1, _center * 1.05 , projY->GetMean());//small offset to distinguish the two graphs
        g_RMS->SetPointError(j - 1, 0.0, projY->GetRMS());

        delete c_proj;
        delete projY;
        delete gaus;
        delete statsBox;
    }
    
    // The plotting part remains the same.
    
    g_RMS->SetTitle(m_title);
    g_RMS->SetMarkerStyle(20);
    g_RMS->SetMarkerColor(kBlack);
    g_RMS->SetLineColor(kBlack);
    g_RMS->SetLineWidth(2);
    g_RMS->GetXaxis()->SetLimits(5.0, 200);
    g_RMS->Draw("AP");

    g->SetTitle(m_title);
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue + 2);
    g->SetLineColor(kBlue + 2);
    g->SetLineWidth(2);
    // g->GetXaxis()->SetLimits(0, 40);
    g->Draw("PSAME");

    g_RMS->SetTitle(m_title);

    // g->GetYaxis()->SetRangeUser(-0.06, 0.06);

    TLine* line = new TLine(0, 0, g->GetXaxis()->GetXmax(), 0);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();

    // Create a TLegend object
    // Arguments are x1, y1, x2, y2 (normalized coordinates)
    TLegend *legend = new TLegend(0.75, 0.75, 0.9, 0.9);
    // "lp" stands for line and point, indicating how the entry should be represented
    // "e" stands for error bars
    legend->AddEntry(g_RMS, "Histograms"   , "ep");
    legend->AddEntry(g    , "Gaussian Fit" , "ep");
    // Draw the legend on the canvas
    legend->Draw();


    // double y_min = g->GetYaxis()->GetXmin();
    // double y_max = g->GetYaxis()->GetXmax();
    g_RMS->GetYaxis()->SetRangeUser(-0.1,0.1);
    // g_RMS->GetXaxis()->SetRangeUser(1.0,200);
    g->GetXaxis()->SetRangeUser(1.0,200);

    c->Update();
    c->SaveAs(m_saveName);
    delete c;
    delete g;
    delete g_RMS;
    // delete shade;
    delete line;
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

// PlotOptionsResponseMatrix implementation
PlotOptionsResponseMatrix::PlotOptionsResponseMatrix(const TString& histName,
                                                     const char* xLabel,
                                                     const char* yLabel,
                                                     const char* saveName)
    : m_histName(histName),
      m_xLabel(xLabel),
      m_yLabel(yLabel),
      m_saveName(saveName) {}

void PlotOptionsResponseMatrix::Plot(TFile* inputFile) {
    // Retrieve the 2D histogram from the input file
    TH2D* h_matrix_orig = (TH2D*)inputFile->Get(m_histName);
    if (!h_matrix_orig) {
        std::cerr << "Error: 2D Histogram " << m_histName << " not found." << std::endl;
        return;
    }

    // Clone the histogram to store percentages
    // h_matrix_orig->Rebin2D(10, 10);// Rebin the original histogram    // merging 10x10 bins into one
    TH2D* h_matrix_perc = (TH2D*)h_matrix_orig->Clone(Form("%s_percentages", h_matrix_orig->GetName()));
    h_matrix_perc->Reset(); // Clear contents for percentage storage keeping the same binning
    h_matrix_perc->SetTitle("");
    h_matrix_perc->GetZaxis()->SetTitle("Percentage (%)");
    h_matrix_perc->GetZaxis()->SetRangeUser(0, 100);

    // Set plot style and margins
    TCanvas* c = new TCanvas("c_response_matrix", "Response Matrix", 1200, 1000);
    // c->SetLogx();
    // c->SetLogy();
    gStyle->SetPaintTextFormat(".0f"); // Show percentages as integers
    gStyle->SetOptStat(0);             // No stats box
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);
    c->SetTopMargin(0.1);
    c->SetBottomMargin(0.1);
    h_matrix_perc->GetXaxis()->SetRangeUser(4, 300);
    h_matrix_perc->GetYaxis()->SetRangeUser(4, 300);

    // Calculate percentages for each bin
    int nbinsX = h_matrix_orig->GetNbinsX();
    int nbinsY = h_matrix_orig->GetNbinsY();
    std::cout << "Number of bins in X: " << nbinsX << ", Y: " << nbinsY << std::endl;
    for (int ix = 1; ix <= nbinsX; ix++) {
        double total_true = 0.0;
        for (int iy = 1; iy <= nbinsY; iy++) {
            total_true += h_matrix_orig->GetBinContent(ix, iy);
        }
        // std::cout << "Total entries for true bin " << ix << ": " << total_true << std::endl;

        if (total_true > 0) {
            for (int iy = 1; iy <= nbinsY; iy++) {
                double content = h_matrix_orig->GetBinContent(ix, iy);
                double percentage = (content / total_true) * 100.0;

                // std::cout << "Bin (" << x << ", " << y << ") = " << content << " -> " << percentage << "%" << std::endl;
                // if (m_histName.Contains("Corr_Q2_EM")) {
                //  std::cout << "bin content: " << content << ", total_true: " << total_true << ", percentage: " << percentage << std::endl;
                // }
                // std::cout << "Bin (" << x << ", " << y << ") = " << content << " -> " << percentage << "%" << std::endl;
                if (percentage>=1.0){
                    h_matrix_perc->SetBinContent(ix, iy, percentage);
                    // std::cout << "Bin (" << x << ", " << y << ") = " << content << " -> " << percentage << "%" << std::endl;
                }

            }
        }
    }

    // Draw the percentage matrix
    h_matrix_perc->GetXaxis()->SetTitle(m_xLabel);
    h_matrix_perc->GetYaxis()->SetTitle(m_yLabel);
    h_matrix_perc->GetXaxis()->SetTitleOffset(1.3);
    h_matrix_perc->GetYaxis()->SetTitleOffset(1.3);
    gPad->Clear();
    h_matrix_perc->Draw("COLZ"); // COLZ creates a heatmap

    // Draw the percentage values as text on the plot
    TText t;
    t.SetTextSize(0.02);
    t.SetTextAlign(22); // Centered
    for (int iy = 1; iy <= nbinsY; iy++) {
        for (int ix = 1; ix <= nbinsX; ix++) {
            double bin_content = h_matrix_perc->GetBinContent(ix, iy);
            double bin_x = h_matrix_perc->GetXaxis()->GetBinCenter(ix);
            double bin_y = h_matrix_perc->GetYaxis()->GetBinCenter(iy);
            
            // Debug output to match calculation
            // std::cout << "Drawing bin (" << ix << ", " << iy << ") at coords (" << bin_x << "," << bin_y << "): " << bin_content << " -> " << Form("%.0f", bin_content) << std::endl;
            
            if (bin_content > 1.0) { // Only draw non-zero percentages
                t.SetTextColor(kBlack); // adjust the color
                t.DrawText(bin_x, bin_y, Form("%.0f", bin_content));
            }
            else if (bin_content > 0.1) {
                t.SetTextColor(kBlue); // adjust the color
                t.DrawText(bin_x, bin_y, "<1%");
            }
        }
    }

    // Add experiment/simulation text
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.2, 0.92, "#bf{ePIC} Simulation");
    latex.DrawLatex(0.65, 0.92, "#bf{Diff. DIS} #10x100 GeV");

    // gStyle->SetPalette(kStarryNight); // Set color palette
    SetCustomPalette(); // Call the function to set custom palette
    // TColor::InvertPalette();   // Invert colors for better visibility
    c->SetLogx();
    c->SetLogy();
    c->Update();
    c->SaveAs(m_saveName);

    // Clean up
    delete h_matrix_perc;
    delete c;
}