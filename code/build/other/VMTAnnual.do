*   Created by Ryan Kellogg on Friday, 19 April, 2019


/*
FILE DESCRIPTION
Loads in VMT data from St Louis Fed
Then creates annual VMT dataset and merges with real fuel prices
*/

/*********** BASIC SETUP *********************/
clear all
set more off, permanently
capture log close
local fname VMTAnnual

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

* tex single number file folder
global texdir = "$repodir/paper/SingleNumberTex"

* Dropbox folder locations
global rawdir = "$dropbox/rawdata/other"
global gasdir = "$dropbox/intdata/fuelprices"
global outdir = "$dropbox/output/other"

* Create a plain text log file to record output
* Log file has same name as do-file
log using "$logdir/`fname'.txt", replace text



***********************************************************************

* Import the VMT data (in millions of miles)
import delimited using "$rawdir/TRFVOLUSM227NFWA.csv", clear
rename trf vmt

* Collapse data to annual
rename date datestring
gen Date = date(datestring,"YMD")
format Date %dDmCY
drop datestring

gen Year = year(Date)
sort Year
collapse(sum) vmt, by(Year)
tempfile tempall
save "`tempall'"

* Keep data from through 2012 (baseline year for simulation)
keep if Year<=2012


* Merge in gasoline price data (in real $)
tempfile tempvmt
save "`tempvmt'"
clear
insheet using "$gasdir/AnnualRealGasPrices.csv"
rename v1 Year
rename v2 PGas
sort Year
merge 1:1 Year using "`tempvmt'"
keep if _merge==3
drop _merge
sort Year

* Unit root tests for vmt data
tsset Year
reg vmt l.vmt Year
local AR1 = _b[l.vmt]
file open fh using "$texdir/vmt_AR1.tex", write replace text
file write fh %8.2f (`AR1')
file close fh
dfgls vmt


* Export PGas and vmt
export delimited using "$outdir/AnnualVMT_PGas.csv", replace


* Export 2012 vmt
clear
use "`tempall'"
keep if Year==2012
keep vmt
export delimited using "$outdir/VMT2012.csv", replace


cap log close
exit, clear
