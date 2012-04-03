*==========================================================================*
*                             mmsel.ado   (v2.0)                           *
*--------------------------------------------------------------------------*
*               Adapted for Sample Selection Correction by				   *
*					 Sami Souabni (sami@souabni.com)					   *
*		School of Business and Economics, Swansea University			   *
*			who bears no responsibility for any errors. 				   *
*					  Based on original code by   						   *
*               Mark Bryan, ISER, University of Essex,					   *
*			who bears no responsibility for any errors.                    *
*--------------------------------------------------------------------------*
*  Downloaded from GitHub												   *
*--------------------------------------------------------------------------*
*  Code to simulate (counterfactual) distributions from quantile           *
*  regressions. Based on Machado and Mata (2005). An option to correct     *
*  for sample selection has been added, using as adaptation of the         *
*  procedure described in Albrecht et. al (2009). Multiple options		   *
*  available for different references groups, following Oaxaca (1973),     *
*  Blinder (1973) Oaxaca Ransom (1994) and Jann (2008).                    *
*--------------------------------------------------------------------------*
*  Requires pid variable							                       *
*--------------------------------------------------------------------------*
*  Creates directories: tmp, logs, results, data                           *
*  Log file: results/`filename'(_sel).log								   *
*  Graphs: results/`filename'(_sel).gph                                    *
*  Other logs: logs/gaps                                                   *
*--------------------------------------------------------------------------*
*  Syntax:                                                                 *
*  mmsel varlist(min=2) [if] [in] , GRoup(varname numeric) Reps(integer 200*
*  )Filename(string) [pooled incgrp group1 Method(integer 2)               *
*  ADJust(varname numeric) REDuced(var list numeric min=2) grponlysel       *
* CONSTRaint(string)]								   					   *
*                                                                          *
*==========================================================================*

*Machado and Mata Program, taking into account sample selection

capture program drop mmsel
program define mmsel
version 10.0
syntax varlist(min=2) [if] [in] , GRoup(varname numeric) Filename(string) [pooled incgrp group1 grponlysel Method(integer 2) Reps(integer 200) SIngle(varlist numeric) ADJust(varname numeric) REDuced(varlist numeric min=2) CONSTRaint(string)]
tokenize `varlist'
global lhs "`1'"
macro shift 1
global rhs "`*'"

local lhs $lhs
local rhs $rhs

global group `group'
global adjust `adjust'
global filename `filename'
global pooled `pooled'
global incgrp `incgrp'
global group1 `group'
global method `method'
global reps `reps'
global reduced `reduced'
global constraint `constraint'
global grponlysel `grponlysel'
global single `single'

capture mkdir tmp
capture mkdir logs
capture mkdir results
capture mkdir data

if `"`single'"'!=""&`"`grponlysel'"'!="" {

ge works=0
replace works==1 if `lhs'>0
foreach gen in 0 1 {
if `gen'==0 {
keep if `group'==`gen'

tempfile presingle
save "`presingle'"

	capture program drop single
	capture program drop dsingle

keep works `single'
single works `single', h(0.4)
di ""
di ""
di "BEWARE: If running 64bit Stata, expect an error r(9999) - use 32bit version instead!"
di ""
di ""
use "`presingle'", clear
compress
	matrix gamma = e(b)
	local k_totaal = colsof(gamma)
	matrix cov_step1 = e(V)
	local b1 = gamma[1,1]
	local b0 = gamma[1,colsof(gamma)]

	local i = 0
	tempvar lambda
	tempvar lambda1
	tempvar z
	
	g `z' = `b0'
	local i = 0
	local svar "$svar"
	
	while (`i' < `k_totaal'-1) {
		local i = `i'+ 1
		local b1 = gamma[1,`i']

		tempvar hulp_x
		local name_x = word("`svar'",1)
		g `hulp_x' = `name_x'
		local svar1 : list local(svar) - local(name_x)
		local svar "`svar1'"
		
		quietly replace `z' = `z' + `b1' * `hulp_x'
	}
	replace Ps1 = `z'
	}
	probit works Ps1
	predict ps, xb
	replace Ps1 = normalden(-ps) / (1-normal(-ps))
	else {
	keep if `group'==`gen'
	}
	tempfile single_`gen'
	save "`single_`gen''"
	}
 use "`single_0'"
 append "`single_1'"
 
  keep if works==1

if `"`grponlysel'"'=="" {
drop if Ps1==.
}
else {
drop if Ps1==.&`group'==0
}

foreach var of local rhs {
drop if `var'==.
}

}


if `"`single'"'!=""&`"`grponlysel'"'=="" {
ge works=0
replace works==1 if `lhs'>0

foreach gen in 0 1 {
keep if `group'==`gen'

tempfile presingle
save "`presingle'"

	capture program drop single
	capture program drop dsingle

keep works `single'
single works `single', h(0.4)
di ""
di ""
di "BEWARE: If running 64bit Stata, expect an error r(9999) - use 32bit version instead!"
di ""
di ""
use "`presingle'", clear
compress
	matrix gamma = e(b)
	local k_totaal = colsof(gamma)
	matrix cov_step1 = e(V)
	local b1 = gamma[1,1]
	local b0 = gamma[1,colsof(gamma)]

	local i = 0
	tempvar lambda
	tempvar lambda1
	tempvar z
	
	g `z' = `b0'
	local i = 0
	local svar "$svar"
	
	while (`i' < `k_totaal'-1) {
		local i = `i'+ 1
		local b1 = gamma[1,`i']

		tempvar hulp_x
		local name_x = word("`svar'",1)
		g `hulp_x' = `name_x'
		local svar1 : list local(svar) - local(name_x)
		local svar "`svar1'"
		
		quietly replace `z' = `z' + `b1' * `hulp_x'
	}
	replace Ps1 = `z'
	probit works Ps1
	predict ps, xb
	replace Ps1 = normalden(-ps) / (1-normal(-ps))
	tempfile single_`gen'
	save "`single_`gen''"
	}
use "`single_0'"
append "`single_1'"
 
keep if works==1

if `"`grponlysel'"'=="" {
drop if Ps1==.
}
else {
drop if Ps1==.&`group'==0
}

foreach var of local rhs {
drop if `var'==.
}

}

save data/data_with_ps, replace

*---------------------------------
* Clean up working directory
*---------------------------------
forval i = 1/99 {

capture erase tmp/xfbf`i'.dta
capture erase tmp/xmbm`i'.dta
capture erase tmp/xfbm`i'.dta
capture erase tmp/xfbp`i'.dta
capture erase tmp/xmbp`i'.dta
}

capture erase tmp/xfbf.dta
capture erase tmp/xmbm.dta
capture erase tmp/xfbm.dta
capture erase tmp/xmbp.dta
capture erase tmp/xfbp.dta
clear
estimates clear

use data/data_with_ps, clear
if `"`if'"'!="" {
keep `if'
}

capture ge `adjust'2=`adjust'^2
capture ge `adjust'3=`adjust'2^2

*---------------------------------
* Build matrix of means for each group
*---------------------------------

matrix accum X = `rhs' if `group'==1
matrix X = X["_cons",1...] /* extract totals */
matrix N = X["_cons","_cons"] /* number of obs */
scalar N = N[1, 1] /* number of obs */
matrix xbarm = X / N /* men */

matrix accum X = `rhs' if `group'==0
matrix X = X["_cons",1...] /* extract totals */
matrix N = X["_cons","_cons"] /* number of obs */
scalar N = N[1, 1] /* number of obs */
matrix xbarf = X / N /* women */

*---------------------------------
* At the mean
*---------------------------------
if `"`grponlysel'"'!=""{
if `"`adjust'"'!="" {
	reg `lhs' `rhs' if `group'==1
	matrix bmols = e(b)
	matrix Vmols = e(V)
	quietly reg `lhs' `rhs' `adjust' if `group'==0
	scalar beta1=_b[`adjust']
	ge `lhs'_adj=`lhs'-(beta1*`adjust') if `group'==0
	quietly reg `lhs'_adj `rhs' if `group'==0
	matrix bfols = e(b)
	matrix Vfols = e(V)
	drop `lhs'_adj
	}
		else {
	reg `lhs' `rhs' if `group'==1
	matrix bmols = e(b)
	matrix Vmols = e(V)

	reg `lhs' `rhs' if `group'==0
	matrix bfols = e(b)
	matrix Vfols = e(V)
	}
}
else{
if `"`adjust'"'!="" {
	quietly reg `lhs' `rhs' `adjust' if `group'==1
	scalar beta1=_b[`adjust']
	ge `lhs'_adj=`lhs'-(beta1*`adjust') if `group'==1
	quietly reg `lhs'_adj `rhs' if `group'==1
	matrix bmols = e(b)
	matrix Vmols = e(V)
	drop `lhs'_adj

	quietly reg `lhs' `rhs' `adjust' if `group'==0
	scalar beta1=_b[`adjust']
	ge `lhs'_adj=`lhs'-(beta1*`adjust') if `group'==0
	quietly reg `lhs'_adj `rhs' if `group'==0
	matrix bfols = e(b)
	matrix Vfols = e(V)
	drop `lhs'_adj
	}
	else {
	reg `lhs' `rhs' if `group'==1
	matrix bmols = e(b)
	matrix Vmols = e(V)

	reg `lhs' `rhs' if `group'==0
	matrix bfols = e(b)
	matrix Vfols = e(V)
	}
}


* Predictions, standard errors and differentials 

matrix olsdiffmat = (bmols - bfols) * xbarf' /* gender diff holding characs constant (mean female characs) */
scalar olsdiff=  olsdiffmat[1,1]
matrix V = Vmols + Vfols
matrix v_olsdiff = xbarf*V*xbarf' /* calculate se */
scalar se_olsdiff = v_olsdiff[1,1]
scalar se_olsdiff = sqrt(se_olsdiff)

*---------------------------------
* Reference Wage Structure Selected?
*---------------------------------
if `"`grponlysel'"'!="" {
if `"`pooled'"' != "" {
di "POOLED REFERENCE WAGE STRUCTURE"
di "GROUP=0 SS CORRECTION ONLY"
}
if `"`group1'"' != "" {
di "GROUP=1 REFERENCE STRUCTURE"
di "GROUP=0 SS CORRECTION ONLY"
}
if `"`incgen'"' != "" {
di "POOLED REF INC GROUP DUMMY"
di "GROUP=0 SS CORRECTION ONLY"
}
}
else {
if `"`pooled'"' != "" {
di "POOLED REFERENCE WAGE STRUCTURE"
}
if `"`group1'"' != "" {
di "GROUP=1 REFERENCE STRUCTURE"
}
if `"`incgen'"' != "" {
di "POOLED REF INC GROUP DUMMY"
}
}
*---------------------------------
* Across each quantile
*---------------------------------

* Raw gaps

forval i = 1/19 {
	local q = `i'*0.05
	quietly: qreg `lhs' `group', quantile(`q') nolog 
	scalar raw`i' = _b[`group']
}

*Sex included due to http://www.ssc.wisc.edu/~jmuniz/jann_oaxaca.pdf p.6

	if `"`pooled'"'!="" {
		adjtest m p
		adjtest f p
	}
	
	if `"`incgen'"'!="" {
		adjtest m i
		adjtest f i
	}
	
	if `"`group1'"'!="" {
		adjtest f g
	}
	
	adjtest m m
	adjtest f f


* Calculate wage gaps

set seed 1

set rmsg on /* See how long gengap takes */

if `"`group1'"' != "" {
gengap g
}
if `"`incgen'"'!="" {
gengap i
}
if `"`pooled'"'!="" {
gengap p
} 

set rmsg off
/*
copy tmp/xfbf.dta logs/xfbf.dta, replace /* save simulated women's distribution for future ref */
copy tmp/xmbm.dta logs/xmbm.dta, replace /* save simulated women's distribution for future ref */
copy tmp/xfbm.dta logs/xfbm.dta, replace /* save simulated women's distribution for future ref */
*/
* Compare simulated distribution with actual distribution for women

use tmp/xfbf, clear /* simulated */
su xfbf, de

use data/data_with_ps, clear
if `"`if'"'!="" {
	keep `if'
}
capture ge `adjust'2=`adjust'^2
capture ge `adjust'3=`adjust'^3

su $lhs, de

* Do replications

if `"`group1'"' != "" {
	gengap g
}
if `"`incgen'"'!="" {
	gengap i
}
if `"`pooled'"'!="" {
	gengap p
} 

if `"`group1'"' != "" {
	simulate "gengap g" ovgap1=r(ovgap1) ovgap2=r(ovgap2) ovgap3=r(ovgap3) ovgap4=r(ovgap4) ovgap5=r(ovgap5) ovgap6=r(ovgap6) ovgap7=r(ovgap7) ovgap8=r(ovgap8) /*
	*/ ovgap9=r(ovgap9) ovgap10=r(ovgap10) ovgap11=r(ovgap11) ovgap12=r(ovgap12) ovgap13=r(ovgap13) ovgap14=r(ovgap14) ovgap15=r(ovgap15) ovgap16=r(ovgap16) ovgap17=r(ovgap17) /*
	*/ ovgap18=r(ovgap18) ovgap19=r(ovgap19) expgap1=r(expgap1) expgap2=r(expgap2) expgap3=r(expgap3) expgap4=r(expgap4) expgap5=r(expgap5) expgap6=r(expgap6) expgap7=r(expgap7) expgap8=r(expgap8) /*
	*/ expgap9=r(expgap9) expgap10=r(expgap10) expgap11=r(expgap11) expgap12=r(expgap12) expgap13=r(expgap13) expgap14=r(expgap14) expgap15=r(expgap15) expgap16=r(expgap16) expgap17=r(expgap17) /*
	*/ expgap18=r(expgap18) expgap19=r(expgap19) unexpgap1=r(unexpgap1) unexpgap2=r(unexpgap2) unexpgap3=r(unexpgap3) unexpgap4=r(unexpgap4) unexpgap5=r(unexpgap5) unexpgap6=r(unexpgap6) unexpgap7=r(unexpgap7) unexpgap8=r(unexpgap8) /*
	*/ unexpgap9=r(unexpgap9) unexpgap10=r(unexpgap10) unexpgap11=r(unexpgap11) unexpgap12=r(unexpgap12) unexpgap13=r(unexpgap13) unexpgap14=r(unexpgap14) unexpgap15=r(unexpgap15) unexpgap16=r(unexpgap16) unexpgap17=r(unexpgap17) /*
	*/ unexpgap18=r(unexpgap18) unexpgap19=r(unexpgap19), reps(`reps') dots saving(logs/gaps) replace 
}

if `"`incgen'"'!="" {
	simulate "gengap i" ovgap1=r(ovgap1) ovgap2=r(ovgap2) ovgap3=r(ovgap3) ovgap4=r(ovgap4) ovgap5=r(ovgap5) ovgap6=r(ovgap6) ovgap7=r(ovgap7) ovgap8=r(ovgap8) /*
	*/ ovgap9=r(ovgap9) ovgap10=r(ovgap10) ovgap11=r(ovgap11) ovgap12=r(ovgap12) ovgap13=r(ovgap13) ovgap14=r(ovgap14) ovgap15=r(ovgap15) ovgap16=r(ovgap16) ovgap17=r(ovgap17) /*
	*/ ovgap18=r(ovgap18) ovgap19=r(ovgap19) expgap1=r(expgap1) expgap2=r(expgap2) expgap3=r(expgap3) expgap4=r(expgap4) expgap5=r(expgap5) expgap6=r(expgap6) expgap7=r(expgap7) expgap8=r(expgap8) /*
	*/ expgap9=r(expgap9) expgap10=r(expgap10) expgap11=r(expgap11) expgap12=r(expgap12) expgap13=r(expgap13) expgap14=r(expgap14) expgap15=r(expgap15) expgap16=r(expgap16) expgap17=r(expgap17) /*
	*/ expgap18=r(expgap18) expgap19=r(expgap19) unexpgap1=r(unexpgap1) unexpgap2=r(unexpgap2) unexpgap3=r(unexpgap3) unexpgap4=r(unexpgap4) unexpgap5=r(unexpgap5) unexpgap6=r(unexpgap6) unexpgap7=r(unexpgap7) unexpgap8=r(unexpgap8) /*
	*/ unexpgap9=r(unexpgap9) unexpgap10=r(unexpgap10) unexpgap11=r(unexpgap11) unexpgap12=r(unexpgap12) unexpgap13=r(unexpgap13) unexpgap14=r(unexpgap14) unexpgap15=r(unexpgap15) unexpgap16=r(unexpgap16) unexpgap17=r(unexpgap17) /*
	*/ unexpgap18=r(unexpgap18) unexpgap19=r(unexpgap19), reps(`reps') dots saving(logs/gaps) replace 
}
if `"`pooled'"'!="" {
	simulate "gengap p" ovgap1=r(ovgap1) ovgap2=r(ovgap2) ovgap3=r(ovgap3) ovgap4=r(ovgap4) ovgap5=r(ovgap5) ovgap6=r(ovgap6) ovgap7=r(ovgap7) ovgap8=r(ovgap8) /*
	*/ ovgap9=r(ovgap9) ovgap10=r(ovgap10) ovgap11=r(ovgap11) ovgap12=r(ovgap12) ovgap13=r(ovgap13) ovgap14=r(ovgap14) ovgap15=r(ovgap15) ovgap16=r(ovgap16) ovgap17=r(ovgap17) /*
	*/ ovgap18=r(ovgap18) ovgap19=r(ovgap19) expgap1=r(expgap1) expgap2=r(expgap2) expgap3=r(expgap3) expgap4=r(expgap4) expgap5=r(expgap5) expgap6=r(expgap6) expgap7=r(expgap7) expgap8=r(expgap8) /*
	*/ expgap9=r(expgap9) expgap10=r(expgap10) expgap11=r(expgap11) expgap12=r(expgap12) expgap13=r(expgap13) expgap14=r(expgap14) expgap15=r(expgap15) expgap16=r(expgap16) expgap17=r(expgap17) /*
	*/ expgap18=r(expgap18) expgap19=r(expgap19) unexpgap1=r(unexpgap1) unexpgap2=r(unexpgap2) unexpgap3=r(unexpgap3) unexpgap4=r(unexpgap4) unexpgap5=r(unexpgap5) unexpgap6=r(unexpgap6) unexpgap7=r(unexpgap7) unexpgap8=r(unexpgap8) /*
	*/ unexpgap9=r(unexpgap9) unexpgap10=r(unexpgap10) unexpgap11=r(unexpgap11) unexpgap12=r(unexpgap12) unexpgap13=r(unexpgap13) unexpgap14=r(unexpgap14) unexpgap15=r(unexpgap15) unexpgap16=r(unexpgap16) unexpgap17=r(unexpgap17) /*
	*/ unexpgap18=r(unexpgap18) unexpgap19=r(unexpgap19), reps(`reps') dots saving(logs/gaps) replace 
}

* Predictions, standard errors and differentials 

gen olscoeff1 = olsdiff

forval i = 1/19 {

	su ovgap`i'
	if `i'==1 {
		gen pred1 = r(mean) if _n==`i'
		gen se_pred1 = r(sd) if _n==`i'
		gen rawgap = raw`i' if _n==`i'
	} 
	else {
		replace pred1 = r(mean) if _n==`i'
		replace se_pred1 = r(sd) if _n==`i'
		replace rawgap = raw`i' if _n==`i'
}
	
}

forval i = 1/19 {

	su expgap`i'
	if `i'==1 {
		gen chars1 = r(mean) if _n==`i'
		gen se_chars1 = r(sd) if _n==`i'
		} 
	else {
		replace chars1 = r(mean) if _n==`i'
		replace se_chars1 = r(sd) if _n==`i'
}
	
}



forval i = 1/19 {

	su unexpgap`i'
	if `i'==1 {
		gen coef1 = r(mean) if _n==`i'
		gen se_coef1 = r(sd) if _n==`i'
		gen loconf1 = coef1 - 1.96 * r(sd) if _n==`i'
		gen hiconf1 = coef1 + 1.96 * r(sd) if _n==`i'
	} 
	else {
		replace coef1 = r(mean) if _n==`i'
		replace se_coef1 = r(sd) if _n==`i'
		replace loconf1 = coef1 - 1.96 * r(sd) if _n==`i'
		replace hiconf1 = coef1 + 1.96 * r(sd) if _n==`i'
}
	
}
capture log c
if `"`adjust'"'!="" {
	log using results/`filename'_sel.log, replace
}
else{
	log using results/`filename'.log, replace
}
* Differentials

scalar li olsdiff
scalar li se_olsdiff

list rawgap pred1 se_pred1 chars1 se_chars1 coef1 se_coef1 in 1/20

gen q = _n/20
	replace q = . if q>0.95 /* Quintile number for graph */

log c

if `"`adjust'"'!="" {	
	twoway connected coef1 hiconf1 loconf1 olscoeff1 rawgap chars1 pred1 q, msymbol(o p p p p o o) mcolor(. . . . . gs10 .)  sch(sami) clpattern(solid dash dash dot "-." "." "__.") /*
	*/ xtitle("Quantile") ytitle("Gap") legend(off)/*
	*/ legend(label(1 "Coefficients") label(2 "Coef 95% confidence intervals") label(4 "OLS") label(5 "Raw gap") label(6 "Characteristics") label(7 "Predicted gap") order(1 2 4 5 6 7) position(3)) /*
	*/ xlabel(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9) /*ylabel(-0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6) yscale(r(-0.1 0.6))*//*
	*/ saving(results/`filename'_sel.gph, replace)
}
else {
	twoway connected coef1 hiconf1 loconf1 olscoeff1 rawgap chars1 pred1 q, msymbol(o p p p p o o) mcolor(. . . . . gs10 .)  sch(sami) clpattern(solid dash dash dot "-." "." "__.") /*
	*/ xtitle("Quantile") ytitle("Gap") legend(off)/*
	*/ legend(label(1 "Coefficients") label(2 "Coef 95% confidence intervals") label(4 "OLS") label(5 "Raw gap") label(6 "Characteristics") label(7 "Predicted gap") order(1 2 4 5 6 7) position(3)) /*
	*/ xlabel(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9) /*ylabel(-0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6) yscale(r(-0.1 0.6))*//*
	*/ saving(results/`filename'.gph, replace)
}
*---------------------------------
* Clean up
*---------------------------------

forval i = 1/99 {
capture erase tmp/xfbf`i'.dta
capture erase tmp/xmbm`i'.dta
capture erase tmp/xfbm`i'.dta
capture erase tmp/xmbp`i'.dta
capture erase tmp/xfbp`i'.dta
}

capture erase tmp/xfbf.dta
capture erase tmp/xmbm.dta
capture erase tmp/xfbm.dta
capture erase tmp/xmbp.dta
capture erase tmp/xfbp.dta

clear
estimates clear

end


