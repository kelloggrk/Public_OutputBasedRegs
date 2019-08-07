*   Created by Ryan Kellogg on Sunday, 7 July, 2019

/*
FILE DESCRIPTION
Loads in BKZ (AER 2013) spreadsheet with NHTSA VMTs
Calculates lifetime expected discounted VMT
*/


/*********** BASIC SETUP *********************/
clear all
set more off, permanently
capture log close
local fname LifetimeDiscMilesPerVehicle

* Set local git directory and local dropbox directory
* Calling the path file works only if the working directory is nested in the repo
pathutil split "`c(pwd)'"
while "`s(filename)'" != "OutputBasedRegs" && "`s(filename)'" != "outputbasedregs" {
	cd ..
	pathutil split "`c(pwd)'"
}
do "globals.do"

* Location for code and log files
global codedir = "$repodir/code/build/other"
global logdir = "$codedir/logfiles"

* Dropbox folder locations
global rawdir = "$dropbox/rawdata/other"
global outdir = "$dropbox/output/other"

* Create a plain text log file to record output
* Log file has same name as do-file
log using "$logdir/`fname'.txt", replace text



***********************************************************************

* Import the BKZ new car data. VMT weighted by survival probability and car vs truck share
* Each row is vintage (vehicle age)
clear
import excel using "$rawdir/BKZ_AER2013_worksheet.xlsx", sheet("NEW Cars") cellrange(G54:G89)

* Discount rate from Allcott and Wozny (2014)
local disc = 0.062

* Create discounted sum of VMT
gen Ind = _n
gen PV = G / (1+`disc')^(Ind-1)
egen LifeVMT = sum(PV)

* Save
keep LifeVMT
duplicates drop
outfile using "$outdir/LifetimeDiscMiles.csv", comma replace

cap log close
exit, clear
