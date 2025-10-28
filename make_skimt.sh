#!/bin/bash

clear
rm -rf ./build/*
rm proton_mandelstam_analysis.root
make
clear

./build/skim_t ./data/filelist.txt
