**Created by Ryan Kellogg on Saturday, 20 April 2019


/*
FILE DESCRIPTION
Loads in the following from StandardsVsFeebates/FuelPrices/rawdata:
EIA natural gas prices for electricity sector

Loads CPI data and merges, and exports annual dataset

*/

/*********** BASIC SETUP *********************/
clear all
set more off, permanently
capture log close
local fname NatGasPricesForElec

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

* Load and clean the natural gas price data
import excel "$rawdir/N3045US3a.xls", sheet("Data 1") firstrow cellrange(A2) clear
drop in 1
rename Sourcekey datestring
rename N3045US3 PGas_Nom
gen Date=date(datestring,"DMY")
format Date %dDmCY
drop datestring
gen Year = year(Date)
destring PGas_Nom, replace
drop Date
keep if Year<2016 
sort Year
tempfile tempngprice
save "`tempngprice'"

* Load CPI data and merge
clear
use "$intdir/CPI.dta"
sort Year
collapse(mean) CPI*, by(Year)
merge 1:1 Year using "`tempngprice'"
keep if _merge==3
drop _merge

* Obtain real Jan 2016 prices in $mmBtu. Output data.
gen PGas = PGas_Nom * CPIJan2016 / CPI
drop PGas_Nom CPI*
replace PGas = PGas / 1.036		// mmBtu per mcf, using https://www.eia.gov/tools/faqs/faq.php?id=45
label variable PGas "Real nat gas price for elec power, Jan 2016 $/mmBtu"
order Year PGas
keep Year PGas
sort Year
save "$intdir/NatGasPrices_Real2016_Annual.dta", replace


capture log close
exit, clear
