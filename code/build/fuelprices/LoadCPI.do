**Created by Ryan Kellogg on Sunday, 28 July, 2019


/*
FILE DESCRIPTION
Loads raw CPI data and cleans it, saving stata data file

*/

/*********** BASIC SETUP *********************/
clear all
set more off, permanently
capture log close
local fname LoadCPI

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

* Create a plain text log file to record output
* Log file has same name as do-file
log using "$logdir/`fname'.txt", replace text


***********************************************************************

* Load CPI data
import excel "$rawdir/SeriesReport-20160411224153_e958e7.xlsx", sheet("BLS Data Series") firstrow cellrange(A11) clear
drop HALF1 HALF2
rename Jan CPI1
rename Feb CPI2
rename Mar CPI3
rename Apr CPI4
rename May CPI5
rename Jun CPI6
rename Jul CPI7
rename Aug CPI8
rename Sep CPI9
rename Oct CPI10
rename Nov CPI11
rename Dec CPI12
destring Year CPI*, replace
drop if Year==.
sort Year
reshape long CPI, i(Year) j(Month)

* Define CPI for 2012 and Jan 2016
gen cpi2012 = CPI if Year==2012
egen CPI2012 = mean(cpi2012)
drop cpi2012
gen cpijan2016 = -99
replace cpijan2016 = CPI if Year==2016 & Month==1
egen CPIJan2016 = max(cpijan2016)
drop cpijan2016


* Clean up and save
drop if CPI==.		// months of 2016 with no data
sort Year Month
label variable CPI "CPI, all urban, all goods less energy, not seasonally adjusted"
save "$intdir/CPI.dta", replace



capture log close
exit, clear
