#!/bin/sh

#-------------------------------------------------
# README
# Batch file for "Output and Attribute-Based Carbon Regulation Under Uncertainty
# This file will run, in sequence:
# 1. All of the Stata .do files that transform raw data into inputs for the Matlab model
# 2. mainscript.m in Matlab, which runs the model
# 3. LaTex file to build the paper

# To run this file, the user must first clone the OutputBasedRegs repository to a local directory (or download the files, maintaining the directory structure), and download the data files (which I store on dropbox) into a separate directory.
# The root directory definitions below should point to these two directories.

# This script is written to be run from a bash shell. This shell is built-in for Max/Unix/Linux users, and is available through git for Windows if you are on a Windows machine. 
# After navigating to the OutputBasedRegs local directory, I recommend running this file via the command: `bash -x makepaper.sh |& tee makepaper_out.txt`. This command will log output and any error messages to makepaper_out.txt.

# This script is designed to work on a Windows machine with 64-bit Stata SE and Matlab. Mac users or users with other Stata versions may have to modify the commands below (e.g., replace `stataSE-64` with `stata-mp`). 

# Windows users may need to add folder paths for Stata and Matlab to their windows path variable.
# Mac users can add Stata or Matlab to their path via a command such as `sudo ln -s /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp /usr/local/bin`

# You may have to run bash as an administrator and/or give your machine the correct permissions via first running the command `chmod +x /path/to/your/script.sh`

# END OF README
#-------------------------------------------------




# DEFINE ROOT DIRECTORIES FOR REPO AND DROPBOX
# MODIFY THESE TO YOUR FOLDER PATHS
REPODIR="C:/Work/OutputBasedRegs"
DBDIR="C:/Users/Ryan Kellogg/Dropbox/OutputBasedRegs"

# DELETE ALL DROPBOX INTERMEDIATE AND OUTPUT DATA FILES
find "$DBDIR/intdata" -type f -delete
find "$DBDIR/output" -type f -delete


#-------------------------------------------------
# RUN THE STATA BUILD FILES

# Load and clean CPI data
StataSE-64 -e do $REPODIR/code/build/fuelprices/LoadCPI.do

# Load and clean gasoline price data. Calculate variance.
StataSE-64 -e do $REPODIR/code/build/fuelprices/GasolinePriceVariance.do

# Load and clean Henry Hub price data. Calculate variance.
StataSE-64 -e do $REPODIR/code/build/fuelprices/HHPriceVariance.do

# Load and clean natural gas prices for electricity users
StataSE-64 -e do $REPODIR/code/build/fuelprices/NatGasPricesForElec.do

# Calculate B_EE for attribute-based standard application
StataSE-64 -e do $REPODIR/code/build/other/EstimateBee.do

# Calculate lifetime discounted miles per vehicle
StataSE-64 -e do $REPODIR/code/build/other/LifetimeDiscMilesPerVehicle.do

# Load, clean, and export footprint data
StataSE-64 -e do $REPODIR/code/build/other/FootprintAnnual.do

# Load, clean, and export miles traveled data
StataSE-64 -e do $REPODIR/code/build/other/VMTAnnual.do

# Load, clean, and export elec generation data
StataSE-64 -e do $REPODIR/code/build/other/ElecGenAnnual.do

# Delete log files in repo root
rm *.log


#-------------------------------------------------
# RUN THE MATLAB MODEL
# GENERATES RESULTS FIGURES, CSV FILES WITH WELFARE EFFECTS AND REGUALTORY SLOPES, AND FORMATTED .TEX TABLES OF PARAMETER INPUTS
matlab -nosplash -wait -nodesktop -r "cd $REPODIR/code/matlabmodel; mainscript; exit"


#-------------------------------------------------
# COMPILE PAPER
cd paper
pdflatex OutputBasedStandards.tex
pdflatex OutputBasedStandards.tex
pdflatex OutputBasedStandards.tex

# Clean up auxillary files
rm *.aux *.bbl *.blg *.log *.out *.gz

exit
