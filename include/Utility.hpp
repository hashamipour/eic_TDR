#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <string> // Needed for std::string
#include <map>
#include <vector>    // Needed for std::map
#include "TROOT.h"

// Declare functions
void SetCustomPalette(const std::string& paletteName = "PetalFlare");
void SetCustomPalette(const int& paletteID);


std::vector<Double_t> GetRoundedLogBins(double min_val, double max_val, int n_bins);
std::vector<Double_t> GetManualQ2Bins();

#endif // UTILITY_HPP
