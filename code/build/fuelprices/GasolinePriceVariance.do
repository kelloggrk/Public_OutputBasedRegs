* Created by Ryan Kellogg on Tuesday, April 16, 2019


/*
FILE DESCRIPTION
Loads in retail fuel price data from EIA
Calculates gasoline price volatility in levels
Exports csv of annual real gasoline prices
*/


/*********** BASIC SETUP *********************/
clear all
set more off, permanently
capture log close
local fname GasolinePriceVariance

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

* Load and clean the gas price data
import excel "$rawdir/PET_PRI_GND_DCUS_NUS_M.xls", sheet("Data 1") firstrow cellrange(A2) clear
drop in 1
* Keep all grades all formulations and regular all grades (latter transitions to former)
keep Sourcekey EMM_EPM0_PTE_NUS_DPG EMM_EPMR_PTE_NUS_DPG
rename Sourcekey datestring
rename EMM_EPM0_PTE_NUS_DPG PGas_All
rename EMM_EPMR_PTE_NUS_DPG PGas_Reg
gen Date=date(datestring,"DMY")
format Date %dDmCY
drop datestring
gen Month = month(Date)
gen Year = year(Date)
destring PGas_All PGas_Reg, replace
order Year Month Date PGas_All PGas_Reg
drop if Year==1990		// partial year
replace PGas_All = PGas_Reg if Date<d(01Apr1993)	// before transition to reformulated
drop PGas_Reg
rename PGas_All PGas_Nom
drop if PGas_Nom==.
tab Year		// check 12 months each year; no missing data except Jan 1991
sort Year Month

* Merge with CPI data and convert to real prices
merge Year Month using "$intdir/CPI.dta"
keep if _merge==3
drop _merge
tab Year		// check 12 months per year

* Get real 2012 prices
gen PGas = PGas_Nom * CPI2012 / CPI
label variable PGas "Real US retail gasoline price, all grades all formulations, 2012 $/gal"
drop PGas_Nom CPI*
sort Date
order Date Year Month PGas
sum PGas
sort Date

* Create three year moving average
sort Year Month
gen N = _n
tsset N
tssmooth ma PGas_ma = PGas, window(35 1 0)		// create 3 year moving average
replace PGas_ma = . if N<36						// not all lags available
save "$intdir/GasPrices_Retail_Real2012.dta", replace


* Calculate volatility of monthly gas prices at one-month to ten-year horizons
* Then take average. Use data through 2012
keep if Year<=2012
local T = 10 * 12	// 10 year lag
gen Diff`T' = PGas - PGas[_n-`T']
* Loop over differences in moving avg, storing std deviations along the way
gen Ind = _n
gen StdDev = 0
forvalues i = 1/`T' {
	gen D`i' = PGas - PGas[_n-`i']
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
outfile using "$outdir/MeanGasolinePriceStdDeviation_Levels.csv", comma replace


* Calculate volatility using the moving average
clear
use "$intdir/GasPrices_Retail_Real2012.dta"
keep if Year<=2012
* Create 10 year differences in the moving avg
local T = 10 * 12	// 10 year lag
gen Diff`T' = PGas_ma - PGas_ma[_n-`T']
* Loop over differences in moving avg, storing std deviations along the way
gen Ind = _n
gen StdDev = 0
forvalues i = 1/`T' {
	gen D`i' = PGas_ma - PGas_ma[_n-`i']
	qui sum D`i' if Diff`T'~=.
	local SD = r(sd)
	replace StdDev = `SD' if Ind==`i'
}
drop D1-D`T' Diff`T'
* Keep standard deviations by lag length and export
keep Ind StdDev
keep if Ind<=`T'
drop Ind
outfile using "$outdir/GasolinePriceStdDeviations_Levels_MA.csv", comma replace


* Export csv of annual real gas prices
clear
use "$intdir/GasPrices_Retail_Real2012.dta" 
keep Year PGas
sort Year
collapse(mean) PGas, by(Year)
outfile using "$intdir/AnnualRealGasPrices.csv", comma replace


* Export csv of annual real gas prices, three year MA
clear
use "$intdir/GasPrices_Retail_Real2012.dta" 
keep Year PGas_ma
sort Year
collapse(mean) PGas_ma, by(Year)
outfile using "$intdir/AnnualRealGasPrices_MA.csv", comma replace


* Export 2012 gas price
clear
use "$intdir/GasPrices_Retail_Real2012.dta" 
keep if Year==2012
keep PGas
collapse(mean) PGas
outfile using "$outdir/RealGasPrice2012.csv", comma replace


capture log close

