#!/bin/bash

# ========================================
# Script to compile and run DDIS analysis
# Usage:
#   ./run_ddis.sh             # Compile & run everything
#   ./run_ddis.sh plot        # Compile & run only plotting
# ========================================

set -e  # Exit on error

if [ "$1" = "plot" ]; then
    echo ">> Only compiling and running plotting code..."
    rm -rf ./figs/*
    g++ -o DDIS_Plots_Q2_OOP DDIS_Plots_Q2_OOP.cpp Plotting.cpp $(root-config --cflags --glibs) -g
    ./DDIS_Plots_Q2_OOP ./DDIS_Skim_Q2_output.root
else
    echo ">> Compiling and running DDIS_Skim_Q2..."
    g++ DDIS_Skim_Q2.cpp -o DDIS_Skim_Q2 $(root-config --cflags --glibs)
    ./DDIS_Skim_Q2 filelist.txt

    echo ">> Cleaning old figures..."
    rm -rf ./figs/*

    echo ">> Compiling and running plotting code..."
    g++ -o DDIS_Plots_Q2_OOP DDIS_Plots_Q2_OOP.cpp Plotting.cpp $(root-config --cflags --glibs) -g
    ./DDIS_Plots_Q2_OOP ./DDIS_Skim_Q2_output.root
fi

