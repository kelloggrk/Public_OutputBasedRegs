*   Created by Ryan Kellogg on Saturday, April 20, 2019


/*
FILE DESCRIPTION
Loads in annual electricity net generation data
Merges with natural gas price data
Exports clean csv to matlab
*/

/*********** BASIC SETUP *********************/
clear all
set more off, permanently
capture log close
local fname ElecGenAnnual

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

* Load in generation data (MWh)
import excel using "$rawdir/annual_generation_state.xls", ///
	cellrange(A2) firstrow clear
	
* Keep US-total generation from all sources
keep if (STATE=="US-TOTAL" | STATE=="US-Total") & TYPEOFPRODUCER=="Total Electric Power Industry"
keep if ENERGYSOURCE=="Total"
drop STATE TYPEOFPRODUCER ENERGYSOURCE F G H I J K
rename GENERATIONMegawatthours gen
rename YEAR Year
sort Year

* Merge in natural gas price data
merge 1:1 Year using "$gasdir/NatGasPrices_Real2016_Annual.dta"
keep if _merge==3
drop _merge

* Export natural gas price and elec gen
export delimited using "$outdir/AnnualElecGen_Pnatgas.csv", replace

* Export 2015 generation and gas price
keep if Year==2015
keep gen PGas
export delimited using "$outdir/ElecGen_Pnatgas_2015.csv", replace



cap log close
exit, clear
