#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <string> // Needed for std::string
#include <map>
#include <vector>    // Needed for std::map
#include "TROOT.h"
#include <iostream> // Needed for std::cout

// Declare functions
void SetCustomPalette(const std::string& paletteName = "PetalFlare");
void SetCustomPalette(const int& paletteID);


std::vector<Double_t> GetRoundedLogBins(double min_val, double max_val, int n_bins);
std::vector<Double_t> GetManualQ2Bins();

class Logger {
public:
    static void warning(const std::string& message) {
        std::cout << "\033[33mWARNING: \033[0m" << message << std::endl;
    }
    
    static void info(const std::string& message) {
        std::cout << "\033[36mINFO: \033[0m" << message << std::endl;
    }
    
    static void error(const std::string& message) {
        std::cerr << "\033[31mERROR: \033[0m" << message << std::endl;  // Use std::cerr!
    }
    
    static void success(const std::string& message) {
        std::cout << "\033[32mSUCCESS: \033[0m" << message << std::endl;
    }
    
    static void debug(const std::string& message) {
        std::cout << "\033[35mDEBUG: \033[0m" << message << std::endl;
    }
};

#endif // UTILITY_HPP
