#!/bin/bash

clear
rm -rf ./build/*
rm -rf ./figs/*
make
clear
#./build/ddis_plots_q2_xy ./DDIS_Skim_Q2_output.root
./build/ddis_plots_t ./proton_mandelstam_analysis.root


#rm -rf ./build/*
#make
#./build/ddis_skim_q2_xy data/filelist.txt DDIS_Skim_Q2_output.root
