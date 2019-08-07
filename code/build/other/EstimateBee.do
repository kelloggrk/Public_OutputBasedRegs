**  Created by Ryan Kellogg on September 17, 2016
**  Moved into GH repo on Monday, April 9, 2018


/*
FILE DESCRIPTION
Loads in marginal fuel economy technolgy costs from Natl Research Council (2015)
Estimates B_{ee}

*/

/*********** BASIC SETUP *********************/
clear all
set more off, permanently
capture log close
local fname EstimateBee

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

* Load in cost curves from Excel
import excel using "$rawdir/NAP2015_CostCurves.xlsx", ///
	sheet("Sheet1") firstrow clear

* Simplify variable names
rename Fuel GPM
rename Cost Cost
rename Estimate Type	

* Create index
gen NN = _n
sort Type NN
bysort Type: gen N = _n
drop NN
bysort Type: egen MaxN = max(N)

* Calculate cumulative gpm improvement
gen GPMInc = GPM[_n-1] - GPM
replace GPMInc = 0 if N==1
gen GPMCum = 0
replace GPMCum = GPMCum[_n-1] + GPMInc if N>1

* Calculate MC in $ / gpm
gen MC = Cost / GPMInc if N>1	

* Regressions
regress MC GPMCum if Type=="high" & N>1
regress MC GPMCum if Type=="low" & N>1
regress MC GPMCum if Type=="high" & N>1 & N<MaxN
gen HighEst = _b[GPMCum]
regress MC GPMCum if Type=="low" & N>1 & N<MaxN
gen LowEst = _b[GPMCum]

* Average high and low cost curves
gen Slope = (HighEst + LowEst) / 2
sum Slope

* Export slope to csv
keep Slope
duplicates drop
outfile using "$outdir/Footprint_B_EE.csv", comma replace


cap log close
exit, clear
