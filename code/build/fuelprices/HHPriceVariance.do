* Created by Ryan Kellogg on Sunday, 28 July, 2019


/*
FILE DESCRIPTION
Loads in Henry Hub front month data from EIA
Calculates natural gas price volatility in levels
*/

/*********** BASIC SETUP *********************/
clear all
set more off, permanently
capture log close
local fname HHPriceVariance

* Set local git directory and local dropbox directory
* Calling the path file works only if the working directory is nested in the repo
pathutil split "`c(pwd)'"
while "`s(filename)'" != "OutputBasedRegs" && "`s(filename)'" != "outputbasedregs" {
	cd ..
	pathutil split "`c(pwd)'"
}
do "globals.do"

* Location for code and log files
global codedir = "$repodir/code/build/fuelprices"
global logdir = "$codedir/logfiles"

* Dropbox folder locations
global rawdir = "$dropbox/rawdata/fuelprices"
global intdir = "$dropbox/intdata/fuelprices"
global outdir = "$dropbox/output/fuelprices"

* Create a plain text log file to record output
* Log file has same name as do-file
log using "$logdir/`fname'.txt", replace text



***********************************************************************

* Load and clean the nominal gas price data
import excel "$rawdir/RNGC1m.xls", sheet("Data 1") firstrow cellrange(A2) clear
drop in 1
rename Sourcekey datestring
rename RNGC1 PHH_Nom
gen Date=date(datestring,"DMY")
format Date %dDmCY
drop datestring
gen Month = month(Date)
gen Year = year(Date)
destring PHH_Nom, replace
order Year Month Date PHH_Nom
drop if PHH_Nom==.
tab Year		// check 12 months each year
sort Year Month
tempfile HHNom
save "`HHNom'"

* Merge with CPI data
merge Year Month using "$intdir/CPI.dta"
keep if _merge==3
drop _merge
tab Year		// check 12 months per year

* Get real Jan 2016 prices
gen PHH = PHH_Nom * CPIJan2016 / CPI
label variable PHH "Real HH front month price, Jan 2016 $/mmBtu"
drop PHH_Nom CPI*
sort Date

* Calculate volatility of monthly HH prices at one-month to 15-year horizons
* Then take average. Only use pre-2016 observations
keep if Year<2016
local T = 15 * 12	// 15 year lag
gen Diff`T' = PHH - PHH[_n-`T']
* Loop over differences in moving avg, storing std deviations along the way
gen Ind = _n
gen StdDev = 0
forvalues i = 1/`T' {
	gen D`i' = PHH - PHH[_n-`i']
	qui sum D`i' if Diff`T'~=.
	local SD = r(sd)
	replace StdDev = `SD' if Ind==`i'
}
drop D1-D`T' Diff`T'
* Keep standard deviations by lag length
keep Ind StdDev
keep if Ind<=`T'
drop Ind

* Export average std deviation over 10 years
egen STD = mean(StdDev)
drop StdDev
rename STD StdDev
duplicates drop
outfile using "$outdir/MeanHHStdDeviation_Levels.csv", comma replace


capture log close
exit, clear

