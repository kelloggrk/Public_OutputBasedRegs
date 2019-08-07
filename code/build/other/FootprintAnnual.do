*   Created by Ryan Kellogg on Sunday, April 28, 2019


/*
FILE DESCRIPTION
Loads in sales-weighted vehicle attribute data from Leard, Linn, McConnell 
attributes paper.
Exports dataset of footprints and gasoline prices
*/

/*********** BASIC SETUP *********************/
clear all
set more off, permanently
capture log close
local fname FootprintAnnual

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
global gasdir = "$dropbox/intdata/fuelprices"
global outdir = "$dropbox/output/other"

* Create a plain text log file to record output
* Log file has same name as do-file
log using "$logdir/`fname'.txt", replace text




***********************************************************************

* Load in sales-weighted vehicle attributes from Excel. 
* Data come from Leard, Linn, McConnell attributes paper, RFF-DP-16-04.pdf
import excel using "$rawdir/attributes_means.xlsx", ///
	sheet("Sheet1") firstrow clear

* Keep cars from start of data (1996) through 2010. Constant, non-ABR
* standard of 27.5 mpg during this period
keep model_year fp_cars
keep if model_year<=2010 
rename fp_cars fp
rename model_year year
rename year Year
sort Year

* Merge in gasoline price data (in real $)
tempfile tempgas
save "`tempgas'"
clear
insheet using "$gasdir/AnnualRealGasPrices_MA.csv"
rename v1 Year
rename v2 PGas
sort Year
merge 1:1 Year using "`tempgas'"
keep if _merge==3
drop _merge
sort Year


* Export footprint and gas price
export delimited using "$outdir/AnnualFootprint_PGas.csv", replace

cap log close
exit
