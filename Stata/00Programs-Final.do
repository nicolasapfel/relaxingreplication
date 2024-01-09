		**********************************
		*  		    PROGRAMS 		     *
		**********************************

	* Contents of this file:
	* SSCIM
	* SSADA
	* SSCIM-LIML (SSCIM with AR and LIML)
	* BP_adj_* commands, which 
	
************************************

* SSCIM ************************************************************************************
cap program drop sscim
prog sscim, eclass
version 14.2
syntax varlist [if] [in] [aw], endog(varlist) exog(varlist) ssstub(string) [vce(string) c(real 0.1) psif(real 1)]

quietly{
tempfile tempcim
save `tempcim', replace

if "`if'" != "" {
	keep `if'
}

preserve
// Default to vce is robust
if "`vce'"=="" {
  di "vce was empty, default to robust"
  loc vce = "robust"
}

loc all: list varlist | exog

// Weigh data if weight option specified
if "`weight'"!=""{
local wts `weight'`exp'
}

marksample touse
markout `touse' `exog' `endog' 
gettoken lhs varlist : varlist
qui sum `lhs' if `touse' 
tempname obs
sca `obs'=r(N)

// Remove collinearity:
local coll `s(collinear)'
	_rmcoll `varlist' if `touse', ///
		`constan' `coll' 
	local varlist `r(varlist)'
local coll `s(collinear)'
	_rmcoll `exog' if `touse', ///
		`constan' `coll' 
	local exog `r(varlist)'
local coll `s(collinear)'
	_rmcoll `endog' if `touse', ///
		`constan' `coll' 
	local endog `r(varlist)'

//Not more than one endogenous x is allowed
tempname nu_endog
scalar `nu_endog'=`:word count `endog''
if `nu_endog' > 1 {
	noi dis in red "Only one endogenous regressor is allowed."
	use `tempcim', clear
	exit 103 // Too many variables specified
}

//User must state for all variables whether they are endogenous or exogenous
loc all: list varlist | exog

tempname nu_all nu_endog nu_exog
scalar `nu_all'=`:word count `all''
scalar `nu_endog'=`:word count `endog''
scalar `nu_exog'=`:word count `exog''
if (`nu_endog'+`nu_exog')!=`nu_all' {
	noi dis in red "Each variable in `varlist' must be either in option endog() or exog()"
	use `tempcim', clear
	exit 499
}

loc exog_xs: list varlist & exog
loc all_exog: list all-endog

//Remove collinearity from exogenous variables
local coll `s(collinear)'
	_rmcoll `all_exog' if `touse', ///
		`constan' `coll' 
	local all_exog `r(varlist)'

loc ivs: list all_exog-exog_xs

//Remove collinearity from IVs
local coll `s(collinear)'
	_rmcoll `ivs' if `touse', ///
		`constan' `coll' 
	local ivs `r(varlist)'	
	
mat b=J(1,`:word count `ivs'',0)
matname b `ivs',c(.)
// nr of all exogenous, of IVs and of endogenous + controls
tempname nu_moments nu_ivs nu_rhs
scalar `nu_moments'=`:word count `all_exog''
scalar `nu_ivs'=`:word count `ivs''
scalar `nu_rhs'=`:word count `varlist''

sca m=`nu_ivs' //Number of IVs
sca k = `nu_all'

//Model is not identified
if `nu_rhs'>`nu_moments' { // nr of endog (1) + controls > nr of IVs + controls
	noi dis in red "There are more parameters than potential instruments. Model is not at all identified."
	use `tempcim', clear
	exit 481
}

// Number of overidentifying restrictions
tempname nuoverid
scalar `nuoverid'=`nu_ivs'-1

//Stop if there are no overidentifying restrictions
if `nuoverid' < 1 {
	noi dis in red "No overidentifying restrictions."
	use `tempcim', clear
	exit 499
}

use `tempcim', clear
// Initial J 
if "`weight'"!=""{
	ivreg2 `lhs' `exog_xs' (`endog'=`ssstub'*) [`wts'], first `vce' partial(`exog_xs')
}
else{
	ivreg2 `lhs' `exog_xs' (`endog'=`ssstub'*), first `vce' partial(`exog_xs')
}
overid
local Jcf = r(j_oid)
*noi di `Jcf'
qui scalar tau=1-`c'/ln(`obs') // Target significance level
if `Jcf'<invchi2(`nuoverid',scalar(tau)) {
	noi dis in red "Hansen test does not reject in the first place. No evidence for invalid instruments. Use all shift-share products."
	restore
	preserve
	noi di as result  "Post-CIM 2SLS with shares as IVs, classical standard errors:"
	if "`weight'"!=""{
	  noi ivreg2 `lhs' `exog_xs' (`endog'=`ssstub'*) [`wts'], first `vce' partial(`exog_xs')
	}
	else{
	  noi ivreg2 `lhs' `exog_xs' (`endog'=`ssstub'*), first `vce' partial(`exog_xs')
	}
	*overid
	matrix first = e(first)
	loc fs = first[3,1]
	*noi di as result "First-stage F-statistic: " `fs'
	restore
	use `tempcim', clear
	exit 499
}
restore

	dis ""
	noi di as text "{hline 59}"
	noi dis as result "{txt}Confidence Interval Method"
	noi di as text "{hline 59}"

local m = m // Nr of IVs as local

matrix define CoefSE = J(m, 2,.) // Matrix for coefs and SEs, m is number of IVs
*local m = m // Nr of IVs as local

// make matrix with coefficients and standard errors
forval nam = 1/`m' { 
	preserve
	token `ivs', parse("")
	rename ``nam'' tx``nam''
		if "`weight'"!=""{
		 qui ivreg2 `lhs' `exog_xs' `ssstub'* (`endog'=tx``nam'') [`wts'], `vce'  partial(`exog_xs') // Important to control for other variables
		}
		else{
		 qui ivreg2 `lhs' `exog_xs' `ssstub'* (`endog'=tx``nam''), `vce' partial(`exog_xs') 
		}
		mat b = e(b)
		scalar D1 = b[1,1]
		mat V = e(V)
		scalar seD1 = sqrt(V[1,1])
		matrix CoefSE[`nam',1] = D1
		matrix CoefSE[`nam',2] = seD1
		restore
}

preserve // dataset for estimation

local rowsbp = (m*(m-1))/2
quietly drop _all
set obs `rowsbp'
gen breakpoints = .
loc bpind = 1

forval j = 1/`m'{
	local el = `m'
	while `el' > `j' {
		scalar betaj = CoefSE[`j',1]
		scalar betar = CoefSE[`el',1]
		scalar absdist = abs(betaj - betar)
		scalar sej = CoefSE[`j',2]
		scalar ser = CoefSE[`el',2]
		scalar sesum = sej + ser
		scalar bp = absdist/sesum
		quietly replace breakpoints = bp in `bpind'
		local bpind = `bpind' + 1
		loc el = `el' - 1
	}
}

gsort -breakpoints
tempfile test1
save `test1', replace
restore // dataset for estimation

loc mC = m - 1 // mC is number of IVs, mC-1 is number of overidentifying restrictions, will be updated in while-loop
scalar psi = `psif' * sqrt(2.01^2*log(_N))  //initial critical value
loc tpv = `c'/ln(_N) 
noi di as result  " Target P-value: `tpv'"
noi di as text "{hline 59}"
noi di as result  _column(18) "Nr. of invalid IVs:"
loc wiold = ""
unab myvars : `ivs'
loc bprun = 0
loc sle = strlen("`ssstub'") + 1 // length of stub plus one

// Start while-loop
while (`Jcf' >  invchi2(`mC', scalar(tau))){ //Stop when J-test not rejected anymore
		preserve // dataset for estimation
		//Empty matrix for Confidence intervals
		matrix define CImat = J(m, 3,.) 
		matrix rownames CImat = `myvars'
		matrix define votes = J(m, m, 0) 
		//Matrix for votes
		//Fill CImat with upper and lower boundaries of confidence intervals
		forv i = 1/`m'{
			matrix CImat[`i',1] = CoefSE[`i',1] - psi*CoefSE[`i',2]
			matrix CImat[`i',2] = CoefSE[`i',1] + psi*CoefSE[`i',2]
		}

		// Fill in index to identify IVs later on
		forv j = 1/`m'{
			token `ivs', parse("")
			loc strind = substr("``j''",`sle',.) //read number in stub
			matrix CImat[`j',3] = `strind'
		}
		//Order by lower boundary of CIs
		mata: st_matrix("CImat", sort(st_matrix("CImat"), 1)) 
		mat CIt = CImat[1...,3]
		mat CIt = CIt'
		mat votes = votes \ CIt
		//Added number of stubs as last line

		// Fill in votes matrix
		forval j = 1/`m'{
			local el = `j' 
			while `el' >= 1 {
			  scalar lb = CImat[`j',1]
			  scalar ub = CImat[`el',2]
			  if lb < ub {
			    matrix votes[`j',`el'] = 1
			  } 
			  else{
			    matrix votes[`j', `el'] = 0
			  }
			  local el = `el' - 1
			}
		}
		// 1 entered in votes matrix if boundary of CI of IV-estimation from certain IV is smaller than upper boundary of one of the preceding CIs
		//Start comparing with itself - first entry naturally is 1, because lower and upper boundary of same CI
		
		// Change votes to dataset
		drop _all
		svmat votes
		egen rowtotal = rowtotal(_all)
		replace rowtotal = -1 if rowtotal > `m' // Otherwise the stub indices are chosen as max
		egen rtmax = max(rowtotal) // maximal rowsum
		replace rowtotal = rtmax if rowtotal < 0 // stub-row is assigned maximal number, s.t. it is kept
		keep if rowtotal==rtmax
		// Choose the row associated with maximal votes
		//Now we have at least two rows: at least one containing the ones and zeros, second: containing stub indices
		if _N != 2 { // If there are more than two rows with maximal number of ones (one for product indices), i.e. if there is tie
		  drop r*
		  gen Jval = . //J-values will be entered here
		  mkmat _all, matrix(selE)
		  tempfile temp1
		  save `temp1', replace //Save multiple rows and fill in J-values in following
		  local N=_N-1 //Do not count the stubs row - we have N-1 rows with maximal votes
		  forval q = 1/`N' { //For each row associated with max, run HS-test
			use `temp1', clear
			* cap drop i // this still needed?
			qui keep if _n == `q' | _n==_N // Keep row of interest and stub row (last one)
			mkmat _all, matrix(sel) // transform to matrix
			mat sel = sel[1...,1..`m']
			mat sel = sel'
			mata: st_matrix("sel", sort(st_matrix("sel"), 2))
			mat sel = sel'
			restore //This restores data used for estimation
			
			// Locals with valid and invalid IV-names
			loc wv
			loc wi
			local nrowsel = rowsof(sel)
			*noi di `nrowsel'
			forval u = 1/`m' {
			  if sel[1,`u'] == 1 {
			    loc add = sel[`nrowsel',`u']
			    loc add = "`ssstub'`add'"
			    loc wv `wv' `add'
			  }
 			  if sel[1,`u'] == 0 {
			    loc add = sel[`nrowsel',`u']
			    loc add = "`ssstub'`add'"
			    loc wi `wi' `add'
			  }
			}

			preserve // dataset, important to restore here, s.t. it can be restored in next iteration
			if "`weight'"!=""{
			  ivreg2 `lhs' `wi' `exog_xs' (`endog'=`wv') [`wts'], `vce' partial(`exog_xs')
			}
			else{
			  ivreg2 `lhs' `exog_xs' `wi' (`endog'=`wv'), `vce' partial(`exog_xs')
			}
			overid
			local Jc = r(j_oid)

			mat selE[`q',`m' + 1] = `Jc' 
			// in last column of selE, save J-statistic associated with the respective model
		  }	// end for

		drop _all
		qui svmat selE, names(col) // maximal votes rows and names of IV-row with Jval
		egen Jmin = min(Jval)
		qui replace Jval = Jmin if Jval == . // so that names-row also kept
		keep if Jval==Jmin
		drop J* // Jval and Jmin
		mkmat _all, matrix(sel) // to matrix
		}
		else{  // If there is only one line in votes associated with maximum
		  cap drop i //still needed?
		  mkmat _all, matrix(sel)
		}
		restore //data
		
		// Build locals that contain valid and invalid varnames
		loc wv
		loc wi
		mat sel = sel[1...,1..`m']
		mat sel = sel'
		mata: st_matrix("sel", sort(st_matrix("sel"), 2))
		mat sel = sel'
		loc nrval = 0
		loc length_wi = 0
		forval u = 1/`m' {
		  if sel[1,`u'] == 1 {
		    loc add = sel[2,`u'] // share index
		    loc add = "`ssstub'`add'" // complete share name
			loc wv `wv' `add' // valid share local
			loc nrval = `nrval' + 1 // there is one more valid IV
		  }
		  if sel[1,`u'] == 0 {
			loc add = sel[2,`u']
			loc add = "`ssstub'`add'"
			loc wi `wi' `add'
			loc length_wi = `length_wi' + 1
		  }
		}
		
		if "`wi'" != "`wiold'" { // skip if identity of invalid IVs did not change
		local wicount: word count `wi' // number of invalid IVs
		
		// J-statistic, data has been restore above
		  if "`weight'"!=""{
		    noi ivreg2 `lhs' `exog_xs' `wi' (`endog'=`wv') [`wts'], `vce' partial(`exog_xs')
		  }
		  else{
		    noi ivreg2 `lhs' `exog_xs' `wi' (`endog'=`wv'), `vce' partial(`exog_xs')
		  }
		  overid
		  local Jcf = r(j_oid)
		  local pcf = r(p_oid)
		  noi di `pcf'
		*}
		  noi di _continue "`wicount'. "
		  loc mC = `nrval' - 1 // update the degrees of freedom
		  loc wiold = "`wi'" // invalid IVs in this while-iteration
		}
		loc bprun = `bprun' + 1
		preserve
		use `test1', clear
		keep in `bprun'
		levelsof(breakpoints), loc(bpnew)
		sca psi = `bpnew'
		restore

		mat drop sel
		cap mat drop selE
  } // end while
  }
  
  loc psi = psi
  di ""
  di as text "{hline 59}"
  di as result "Final P-value: `pcf'"
  di as result "Final Psi: `psi'"
  di as text "{hline 59}"
  di as result "Valid shares: `wv' "
  di as text "{hline 59}"
  di as result "Invalid shares: `wi'"
  di as text "{hline 59}"
  di as result "Invalid shares: `length_wi'"
  di as text "{hline 59}"
  
  //Generate post-CIM shift-share IV
  di "Post-CIM 2SLS with shares as IVs, classical standard errors:"
  if "`weight'"!=""{
    ivreg2 `lhs' `exog_xs' `wi' (`endog'=`wv') [`weight'=`exp'], first `vce' partial(`exog_xs')
  }
  else {
    ivreg2 `lhs' `exog_xs' `wi' (`endog'=`wv'), first `vce' partial(`exog_xs')
  }
  overid
  local J = r(j_oid)
  di `J'
  use `tempcim', clear
  
  ereturn local wv `wv'
  ereturn local wi `wi'
  ereturn local length_wi `length_wi'
   
end



* SSCIM with Anderson-Rubin test
* 	and LIML

cap program drop sscimliml
prog sscimliml, eclass
version 14.2
syntax varlist [if] [in] [aw], endog(varlist) exog(varlist) ssstub(string) [vce(string) c(real 0.1) psif(real 1)]

quietly{
tempfile tempcim
save `tempcim', replace

if "`if'" != "" {
	keep `if'
}

preserve
// Default to vce is robust
if "`vce'"=="" {
  di "vce was empty, default to robust"
  loc vce = "robust"
}

loc all: list varlist | exog

// Weigh data if weight option specified
if "`weight'"!=""{
local wts `weight'`exp'
}

marksample touse
markout `touse' `exog' `endog' 
gettoken lhs varlist : varlist
qui sum `lhs' if `touse' 
tempname obs
sca `obs'=r(N)

// Remove collinearity:
local coll `s(collinear)'
	_rmcoll `varlist' if `touse', ///
		`constan' `coll' 
	local varlist `r(varlist)'
local coll `s(collinear)'
	_rmcoll `exog' if `touse', ///
		`constan' `coll' 
	local exog `r(varlist)'
local coll `s(collinear)'
	_rmcoll `endog' if `touse', ///
		`constan' `coll' 
	local endog `r(varlist)'

//Not more than one endogenous x is allowed
tempname nu_endog
scalar `nu_endog'=`:word count `endog''
if `nu_endog' > 1 {
	noi dis in red "Only one endogenous regressor is allowed."
	use `tempcim', clear
	exit 103 // Too many variables specified
}

//User must state for all variables whether they are endogenous or exogenous
loc all: list varlist | exog

tempname nu_all nu_endog nu_exog
scalar `nu_all'=`:word count `all''
scalar `nu_endog'=`:word count `endog''
scalar `nu_exog'=`:word count `exog''
if (`nu_endog'+`nu_exog')!=`nu_all' {
	noi dis in red "Each variable in `varlist' must be either in option endog() or exog()"
	use `tempcim', clear
	exit 499
}

loc exog_xs: list varlist & exog
loc all_exog: list all-endog

//Remove collinearity from exogenous variables
local coll `s(collinear)'
	_rmcoll `all_exog' if `touse', ///
		`constan' `coll' 
	local all_exog `r(varlist)'

loc ivs: list all_exog-exog_xs

//Remove collinearity from IVs
local coll `s(collinear)'
	_rmcoll `ivs' if `touse', ///
		`constan' `coll' 
	local ivs `r(varlist)'	
	
mat b=J(1,`:word count `ivs'',0)
matname b `ivs',c(.)
// nr of all exogenous, of IVs and of endogenous + controls
tempname nu_moments nu_ivs nu_rhs
scalar `nu_moments'=`:word count `all_exog''
scalar `nu_ivs'=`:word count `ivs''
scalar `nu_rhs'=`:word count `varlist''

sca m=`nu_ivs' //Number of IVs
sca k = `nu_all'

//Model is not identified
if `nu_rhs'>`nu_moments' { // nr of endog (1) + controls > nr of IVs + controls
	noi dis in red "There are more parameters than potential instruments. Model is not at all identified."
	use `tempcim', clear
	exit 481
}

// Number of overidentifying restrictions
tempname nuoverid
scalar `nuoverid'=`nu_ivs'-1

//Stop if there are no overidentifying restrictions
if `nuoverid' < 1 {
	noi dis in red "No overidentifying restrictions."
	use `tempcim', clear
	exit 499
}

use `tempcim', clear
// Initial J 
ivreg2 `lhs' `exog_xs' (`endog'=`ssstub'*) [`wts'], first liml partial(`exog_xs')
local Jcf = e(arubin)
*di "Initial J"
*noi di `Jcf'

qui scalar tau=1-`c'/ln(`obs') // Target significance level
if `Jcf'<invchi2(`nuoverid',scalar(tau)) {
	noi dis in red "Hansen test does not reject in the first place. No evidence for invalid instruments. Use all shift-share products."
	restore
	preserve
	noi di as result  "Post-CIM 2SLS with shares as IVs, classical standard errors:"
	if "`weight'"!=""{
	  ivreg2 `lhs' `exog_xs' (`endog'=`ssstub'*) [`wts'], first liml partial(`exog_xs')
	}
	else{
	  ivreg2 `lhs' `exog_xs' (`endog'=`ssstub'*), first liml partial(`exog_xs')
	}
    matrix first = e(first)
	loc fs = first[3,1]
	*noi di as result "First-stage F-statistic: " `fs' 
	restore
	use `tempcim', clear
	exit 499
}
restore

	dis ""
	noi di as text "{hline 59}"
	noi dis as result "{txt}Confidence Interval Method"
	noi di as text "{hline 59}"

local m = m // Nr of IVs as local
matrix define CoefSE = J(m, 2,.) // Matrix for coefs and SEs, m is number of IVs
*local m = m // Nr of IVs as local

forval nam = 1/`m' { // make matrix with coefficients and standard errors
	preserve
	token `ivs', parse("")
	rename ``nam'' tx``nam''
		if "`weight'"!=""{
		 cap ivreg2 `lhs' `exog_xs' `ssstub'* (`endog'=tx``nam'') [`wts'], `vce' partial(`exog_xs') // Important to control for other variables
		}
		else{
		 cap ivreg2 `lhs' `exog_xs' `ssstub'* (`endog'=tx``nam''), `vce' partial(`exog_xs')
		}
		mat b = e(b)
		scalar D1 = b[1,1]
		mat V = e(V)
		scalar seD1 = sqrt(V[1,1])
		matrix CoefSE[`nam',1] = D1
		matrix CoefSE[`nam',2] = seD1
		restore
}

preserve // dataset for estimation

local rowsbp = (m*(m-1))/2
quietly drop _all
set obs `rowsbp'
gen breakpoints = .
loc bpind = 1

forval j = 1/`m'{
	local el = `m'
	while `el' > `j' {
		scalar betaj = CoefSE[`j',1]
		scalar betar = CoefSE[`el',1]
		scalar absdist = abs(betaj - betar)
		scalar sej = CoefSE[`j',2]
		scalar ser = CoefSE[`el',2]
		scalar sesum = sej + ser
		scalar bp = absdist/sesum
		quietly replace breakpoints = bp in `bpind'
		local bpind = `bpind' + 1
		loc el = `el' - 1
	}
}

gsort -breakpoints
tempfile test1
save `test1', replace
restore // dataset for estimation

loc mC = m - 1 // mC is number of IVs, mC-1 is number of overidentifying restrictions, will be updated in while-loop
scalar psi = `psif' * sqrt(2.01^2*log(_N))  //initial critical value
loc tpv = `c'/ln(_N) 
noi di as result  " Target P-value: `tpv'"
noi di as text "{hline 59}"
noi di as result _column(18) "Nr. of invalid IVs:"
loc wiold = ""
unab myvars : `ivs'
loc bprun = 0
loc sle = strlen("`ssstub'") + 1 // length of stub plus one

// Start while-loop
while (`Jcf' >  invchi2(`mC', scalar(tau))){ //Stop when J-test not rejected anymore
		preserve // dataset for estimation
		//Empty matrix for Confidence intervals
		matrix define CImat = J(m, 3,.) 
		matrix rownames CImat = `myvars'
		matrix define votes = J(m, m, 0) 
		//Matrix for votes
		//Fill CImat with upper and lower boundaries of confidence intervals
		forv i = 1/`m'{
			matrix CImat[`i',1] = CoefSE[`i',1] - psi*CoefSE[`i',2]
			matrix CImat[`i',2] = CoefSE[`i',1] + psi*CoefSE[`i',2]
		}

		// Fill in index to identify IVs later on
		forv j = 1/`m'{
			token `ivs', parse("")
			loc strind = substr("``j''",`sle',.) //read number in stub
			matrix CImat[`j',3] = `strind'
		}
		//Order by lower boundary of CIs
		mata: st_matrix("CImat", sort(st_matrix("CImat"), 1)) 
		mat CIt = CImat[1...,3]
		mat CIt = CIt'
		mat votes = votes \ CIt
		//Added number of stubs as last line

		// Fill in votes matrix
		forval j = 1/`m'{
			local el = `j' 
			while `el' >= 1 {
			  scalar lb = CImat[`j',1]
			  scalar ub = CImat[`el',2]
			  if lb < ub {
			    matrix votes[`j',`el'] = 1
			  } 
			  else{
			    matrix votes[`j', `el'] = 0
			  }
			  local el = `el' - 1
			}
		}
		// 1 entered in votes matrix if boundary of CI of IV-estimation from certain IV is smaller than upper boundary of one of the preceding CIs
		//Start comparing with itself - first entry naturally is 1, because lower and upper boundary of same CI
		
		// Change votes to dataset
		drop _all
		svmat votes
		egen rowtotal = rowtotal(_all)
		replace rowtotal = -1 if rowtotal > `m' // Otherwise the stub indices are chosen as max
		egen rtmax = max(rowtotal) // maximal rowsum
		replace rowtotal = rtmax if rowtotal < 0 // stub-row is assigned maximal number, s.t. it is kept
		keep if rowtotal==rtmax
		// Choose the row associated with maximal votes
		//Now we have at least two rows: at least one containing the ones and zeros, second: containing stub indices

		if _N != 2 { // If there are more than two rows with maximal number of ones (one for product indices), i.e. if there is tie
		  drop r*
		  gen Jval = . //J-values will be entered here
		  mkmat _all, matrix(selE)
		  tempfile temp1
		  save `temp1', replace //Save multiple rows and fill in J-values in following
		  local N=_N-1 //Do not count the stubs row - we have N-1 rows with maximal votes
		  forval q = 1/`N' { //For each row associated with max, run HS-test
			use `temp1', clear
			qui keep if _n == `q' | _n==_N // Keep row of interest and stub row (last one)
			mkmat _all, matrix(sel) // transform to matrix
			mat sel = sel[1...,1..`m']
			mat sel = sel'
			mata: st_matrix("sel", sort(st_matrix("sel"), 2))
			mat sel = sel'
			restore //This restores data used for estimation
			
			// Locals with valid and invalid IV-names
			loc wv
			loc wi
			local nrowsel = rowsof(sel)
			*noi di `nrowsel'
			forval u = 1/`m' {
			  if sel[1,`u'] == 1 {
			    loc add = sel[`nrowsel',`u']
			    loc add = "`ssstub'`add'"
			    loc wv `wv' `add'
			  }
 			  if sel[1,`u'] == 0 {
			    loc add = sel[`nrowsel',`u']
			    loc add = "`ssstub'`add'"
			    loc wi `wi' `add'
			  }
			}
			
			preserve // dataset, important to restore here, s.t. it can be restored in next iteration
			if "`weight'"!=""{
			  qui ivreg2 `lhs' `exog_xs' `wi' (`endog'=`wv') [`wts'], liml partial(`exog_xs')
			}
			else{
			  qui ivreg2 `lhs' `exog_xs' `wi' (`endog'=`wv'), liml partial(`exog_xs')
			}
			local Jc = e(arubin)  // J-statistic
			mat selE[`q',`m' + 1] = `Jc' 
			// in last column of selE, save J-statistic associated with the respective model
		  }	// end for

		drop _all
		qui svmat selE, names(col) // maximal votes rows and names of IV-row with Jval
		egen Jmin = min(Jval)
		qui replace Jval = Jmin if Jval == . // so that names-row also kept
		keep if Jval==Jmin

		drop J* // Jval and Jmin
		mkmat _all, matrix(sel) // to matrix
		}
		else{  // If there is only one line in votes associated with maximum
		  cap drop i //still needed?
		  mkmat _all, matrix(sel)
		}
		restore //data
		
		// Build locals that contain valid and invalid varnames
		loc wv
		loc wi
		mat sel = sel[1...,1..`m']
		mat sel = sel'
		mata: st_matrix("sel", sort(st_matrix("sel"), 2))
		mat sel = sel'
		loc nrval = 0
		loc length_wi = 0
		forval u = 1/`m' {
		  if sel[1,`u'] == 1 {
		    loc add = sel[2,`u'] // share index
		    loc add = "`ssstub'`add'" // complete share name
			loc wv `wv' `add' // valid share local
			loc nrval = `nrval' + 1 // there is one more valid IV
		  }
		  if sel[1,`u'] == 0 {
			loc add = sel[2,`u']
			loc add = "`ssstub'`add'"
			loc wi `wi' `add'
			loc length_wi = `length_wi' + 1
		  }
		}
		
		if "`wi'" != "`wiold'" { // skip if identity of invalid IVs did not change
		local wicount: word count `wi' // number of invalid IVs
		
		// J-statistic, data has been restore above
		  if "`weight'"!=""{
		    ivreg2 `lhs' `exog_xs' `wi' (`endog'=`wv') [`wts'], liml partial(`exog_xs')
		  }
		  else{
		    ivreg2 `lhs' `exog_xs' `wi' (`endog'=`wv'), liml partial(`exog_xs')
		  }
		  loc Jcf = e(arubin)
		  loc pcf = e(arubinp)
		  di `pcf'
		  di `Jcf'
		*}
		  noi di _continue "`wicount'. "
		  loc mC = `nrval' - 1 // update the degrees of freedom
		  loc wiold = "`wi'" // invalid IVs in this while-iteration
		}
		loc bprun = `bprun' + 1
		preserve
		use `test1', clear
		keep in `bprun'
		levelsof(breakpoints), loc(bpnew)
		sca psi = `bpnew'
		restore

		mat drop sel
		cap mat drop selE
  } // end while
  }
  loc psi = psi
  di ""
  di as text "{hline 59}"
  di as result "Final P-value: `pcf'"
  di as result "Final Psi: `psi'"
  di as text "{hline 59}"
  di as result "Valid shares: `wv' "
  di as text "{hline 59}"
  di as result "Invalid shares: `wi'"
  di as text "{hline 59}"
  
  //Generate post-CIM shift-share IV
  di "Post-CIM 2SLS with shares as IVs, classical standard errors:"
  if "`weight'"!=""{
    ivreg2 `lhs' `exog_xs' `wi' (`endog'=`wv') [`weight'=`exp'], first liml partial(`exog_xs')
  }
  else {
    ivreg2 `lhs' `exog_xs' `wi' (`endog'=`wv'), first liml partial(`exog_xs')
  }
  di e(arubinp)
  use `tempcim', clear
  
  ereturn local wv `wv'
  ereturn local wi `wi'
  ereturn local length_wi `length_wi'
   
end


* SSADA ******************************************************************************************

* This is an adaptation of the Stata program of Helmut Farbmacher "SIVREG: Stata module to perform adaptive Lasso with some invalid instruments"

cap program drop ssada
prog ssada, eclass
version 14.2
syntax varlist [if] [in] [aw], endog(varlist) exog(varlist) id(string) [vce(string) c(real 0.1)]

tempfile tempssada
save `tempssada', replace

preserve
if "`in'" != "" {
	keep `in'
}

if "`if'" != "" {
	keep `if'
}

*preserve
// Default to vce is robust
if "`vce'"=="" {
  di "vce was empty, default to robust"
  loc vce = "robust"
}

loc all: list varlist | exog

// Weigh data if weight option specified
if "`weight'"!=""{
  local wts `weight'`exp'
  loc exp = substr("`exp'",3,.)
  quietly replace `exp' = sqrt(`exp')
foreach i of varlist `all' {
  qui replace `i' = `exp' * `i'
  }
}

marksample touse
markout `touse' `exog' `endog' 
gettoken lhs varlist : varlist
qui sum `lhs' if `touse'
tempname obs
sca `obs'=r(N)

//Remove collinearity:
local coll `s(collinear)'
	_rmcoll `varlist' if `touse', ///
		`constan' `coll' 
	local varlist `r(varlist)'
local coll `s(collinear)'
	_rmcoll `exog' if `touse', ///
		`constan' `coll' 
	local exog `r(varlist)'
local coll `s(collinear)'
	_rmcoll `endog' if `touse', ///
		`constan' `coll' 
	local endog `r(varlist)'	
	
//Not more than one endogenous x is allowed
tempname nu_endog
scalar `nu_endog'=`:word count `endog''
if `nu_endog' > 1 {
	dis in red "Only one endogenous regressor is allowed."
	use `tempssada', clear
	exit 103 // Too many variables specified
}

//User must state for all variables whether they are endogenous or exogenous
loc all: list varlist | exog

tempname nu_all nu_endog nu_exog
scalar `nu_all'=`:word count `all''
scalar `nu_endog'=`:word count `endog''
scalar `nu_exog'=`:word count `exog''
if (`nu_endog'+`nu_exog')!=`nu_all' {
	dis in red "Each variable in `varlist' must be either in option endog() or exog()"
	use `tempssada', clear
	exit 499
}

loc exog_xs: list varlist & exog
loc all_exog: list all-endog

//Remove collinearity from exogenous variables
local coll `s(collinear)'
	_rmcoll `all_exog' if `touse', ///
		`constan' `coll' 
	local all_exog `r(varlist)'

loc ivs: list all_exog-exog_xs

//Remove collinearity from IVs
local coll `s(collinear)'
	_rmcoll `ivs' if `touse', ///
		`constan' `coll' 
	local ivs `r(varlist)'	
	
mat b=J(1,`:word count `ivs'',0)
matname b `ivs',c(.)
tempname nu_moments nu_ivs nu_rhs
scalar `nu_moments'=`:word count `all_exog''
scalar `nu_ivs'=`:word count `ivs''
scalar `nu_rhs'=`:word count `varlist''

sca k = `nu_all'

//Model is not identified
if `nu_rhs'>`nu_moments' {
	dis in red "There are more parameters than potential instruments. Model is not at all identified."
	use `tempssada', clear
	exit 481
}

mat c_cons123=`c'

//Partialling out exog_xs
tempvar lhs_fwl endog_fwl

if "`weight'"!="" {
qui reg `lhs' `exp' `exog_xs', noconstant
qui predict `lhs_fwl', resid
qui reg `endog' `exp' `exog_xs', nocons
qui predict `endog_fwl', resid
}
else {
qui reg `lhs' `exog_xs'
qui predict `lhs_fwl', resid
qui reg `endog' `exog_xs'
qui predict `endog_fwl', resid
}

foreach var of varlist `ivs' {
	tempvar `var'_fwl
	if "`weight'"!="" {
	qui reg `var' `exp' `exog_xs', nocons
	qui predict ``var'_fwl', resid
	}
	else {
	qui reg `var' `exog_xs'
	qui predict ``var'_fwl', resid
	}
	local ivs_fwl `"`ivs_fwl' ``var'_fwl'"'
} // Same for IVs

tempname nuoverid
scalar `nuoverid'=`nu_ivs'-1

//Stop if there are no overidentifying restrictions
if `nuoverid' < 1 {
	dis in red "No overidentifying restrictions."
	use `tempssada', clear
	exit 499
}

tempfile tempmata
save `tempmata', replace
restore

preserve
if "`weight'"!=""{
ivreg2 `lhs' `exog_xs' (`endog'=`ivs') [`weight'=`exp'] `if', first rob partial(`exog_xs')
}
else {
ivreg2 `lhs' `exog_xs' (`endog'=`ivs') `if', first rob partial(`exog_xs')
}
overid
local Jc = r(j_oid)
di `Jc'

use `tempmata', clear
qui scalar tau=1-`c'/ln(`obs')
noi di "needs to be smaller than this in order to stop"
di invchi2(`nuoverid',scalar(tau))
if `Jc' < invchi2(`nuoverid',scalar(tau)) {
	dis in red "Hansen test does not reject in the first place. No evidence for invalid instruments. Use all shift-share products."
	restore
	preserve
	di "Standard 2SLS:"
	if "`weight'"!=""{
	  noi ivreg2 `lhs' `exog_xs' (`endog'=`ivs') [`wts'], first rob partial(`exog_xs')
	}
	else{
	  noi ivreg2 `lhs' `exog_xs' (`endog'=`ivs'), first rob partial(`exog_xs')
	}
	overid
	local pfinal = r(p_oid)
	noi di "Final P-value: `pfinal'"
   	matrix first = e(first)
	loc fs = first[3,1]
	restore
	use `tempssada', clear
	exit 499
}

	dis ""
	di as text "{hline 59}"
	dis "{txt}Adaptive Lasso with some invalid IVs"
	di as text "{hline 59}"

qui putmata yabc=`lhs_fwl' xabc=`endog_fwl' Zabc=(`ivs_fwl'), replace
mata {	
	n=rows(Zabc) // number observations
	m=cols(Zabc) // number IVs
	c_cons=st_matrix("c_cons123") // c-value for cutoff
	adaptive=st_matrix("adaptive123") // Dummy that is one when adaptive Lasso is used
	
	//normalization
	Zabc=Zabc:-mean(Zabc) // elementwise subtraction
	sZabc=sqrt(diagonal(Zabc'Zabc)/n) // standard deviation
	Zabc=Zabc:/sZabc' 
	xabc=xabc:-mean(xabc)
	yabc=yabc:-mean(yabc)
	
	izz=luinv(Zabc'Zabc)
	pihat=izz*Zabc'xabc // gammahat
	redform=izz*Zabc'yabc // Gammahat
	
	*pihat
	*redform
	
	xhat=Zabc*pihat // first-stage fitted values
	Zeta=luinv(xhat'xhat)*xhat'Zabc // coefficient of regression from Z on Dhat
	Zt=Zabc-xhat*Zeta // residuals from said regression -> Ztilde
	
	//calculate 2SLS and 2GMM for initial Hansen statistic
	b2sls=luinv(xabc'Zabc*izz*Zabc'xabc)*(xabc'Zabc*izz*Zabc'yabc)
	u2sls=yabc-xabc*b2sls
	Zu=Zabc:*u2sls // elementwise multiplication
	Wgmm=luinv(Zu'Zu)
	bgmm=luinv(xabc'Zabc*Wgmm*Zabc'xabc)*(xabc'Zabc*Wgmm*Zabc'yabc)
	ghat=Zabc'(yabc-xabc*bgmm)
	J=ghat'Wgmm*ghat
	Jinitial=J
	
	//median estimator to get a consistent estimate for alpha
	ratio=redform:/pihat

	betamed=mm_median(ratio)

	alphamed=redform-pihat*betamed
	vau=1
	W=diag(abs(alphamed):^vau)
	
	Zt=Zt*W		// Ztilde weighted			
	
	// avoid to get negative taus (not possible in chi-square distribution)
	tau=1-c_cons/ln(n)

	if (tau<0) {
		"Note: Your choice of c leads to a negative tau; replaced by default (c=0.1)"
		tau=1-0.1/ln(n)
	}

	//LARS:
	kx=1
	alpha=J(m-1,m,0)
	beta=J(m-1,1,0)
	A=J(1,m,0)
	mu=J(n,1,0)
	ssval=J(1,m,1)
	steps=1	
	mR=m // number of IVs, mR-1 is number of overidentifying restrictions
	
	while (J>invchi2(mR-1,tau)) { //mR - 1 overidentifying restrictions
		cmu=Zt'(yabc-mu)/n
		
		if (kx<2) {
			c=colmaxabs(cmu)
			threshold=c/1000
			ss=(abs(round(cmu,threshold)):==round(c,threshold))'
			A=ss
			j=select((1..cols(ss)), (abs(ss) :== max(abs(ss))))
			Aset=j
		}
		else {	
			c=colmaxabs(select(cmu,ssval'))
			threshold=c/1000
			cmax=c*J(m,1,1)
			ss=(abs(reldif(cmax,abs(cmu))):<threshold)'-A
			A=(abs(reldif(cmax,abs(cmu))):<threshold)'
			j=select((1..cols(ss)), (abs(ss) :== max(abs(ss))))
			Aset=(Aset, j)
		}
		
		sign=sign(cmu)
		signA=select(sign',A)
		SA=diag(signA)
		XA=select(Zt,A)
		XA=XA*SA
		G=XA'XA/n
		iG=luinv(G)
		iota=J(cols(XA),1,1)
		B=1/sqrt(iota'iG*iota)
		w=B*iG*iota	
		u=XA*w
		b=Zt'u/n
		Aset_ordered=sort(Aset',1)'	
		dhat=(signA'):*w
		Dhat=J(m,1,0)
		Asetrow = uniqrows(Aset_ordered')
		Aset_ordered = Asetrow' 
		for(j=1;j<=cols(Aset_ordered)&j<=rows(dhat);j++) {
			jnew=Aset_ordered[j]
			Dhat[jnew]=dhat[j]	
		}
		//check for drops (Lasso modification) and get gammatilde
		gammatilde=.
		if (kx>1) {
			gam=alpha[max((kx-1,1)),.]:/(Dhat')
			drops = gam:<0
			if (drops!=J(1,m,0)) {
				gammatilde=rowmin(abs(select(gam,drops)))			//jump at most to where the first variable crosses zero (this is minimum absolue value of the negative numbers in length)
				if (gammatilde!=.) {
					gamtildej=select((1..m), (abs(gam) :== gammatilde))		//drop minj from selected set A
				}
			}
		}
		//get gammatilde
		ccm=((c:-cmu):/(B:-b))
		ccp=((c:+cmu):/(B:+b))
		ssval=J(1,cols(A),1)-A
		ccm=select(ccm',ssval)'
		ccp=select(ccp',ssval)'
		ccm=abs(ccm)+(ccm:<0)*1000000
		ccp=abs(ccp)+(ccp:<0)*1000000
		ccm=ccm+(ccm:==0)*1000000
		ccp=ccp+(ccp:==0)*1000000
		gammahat=rowmin(colmin((ccm,ccp)))	
		//decision rule b/w gammahat and gammatilde
		if (gammatilde<gammahat) {							//Lasso step
			steps=steps-1
			alpha=(alpha \ J(1,cols(alpha),0))
			beta=(beta \ 0)
			mu=mu+gammatilde*u
			alpha[max((kx,1)),.]=alpha[max((kx-1,1)),.]+gammatilde*Dhat'
			A[.,gamtildej]=0
			ssval[.,gamtildej]=1
			deletej=Aset:-(gamtildej:*J(1,cols(Aset),1))
			Aset=select(Aset,deletej)
		}
		else {												//LARS step
			mu=mu+gammahat*u
			alpha[max((kx,1)),.]=alpha[max((kx-1,1)),.]+gammahat*Dhat'
		}		
		
		//retrieve alpha's if adaptive lasso (Note: W=I(n) if standard (non-adaptive) Lasso)
		alpha[max((kx,1)),.]=alpha[max((kx,1)),.]*W
		
		Zcontrol_abc=select(Zabc,A)
		beta_iv=xhat'(yabc-Zabc*alpha[kx,.]')/(xhat'xhat)
		beta[kx,1]=beta_iv

		lambda_final=c
				
		kx=kx+1
		steps=steps+1
		
		//re-calculate 2SLS and 2GMM for Hansen statistic
		RX=(xabc,Zcontrol_abc) // also need to include the invalid IVs as controls
		b2sls=luinv(RX'Zabc*izz*Zabc'RX)*RX'Zabc*izz*Zabc'yabc // 2sls beta
		u2sls=yabc-RX*b2sls // 2sls residuals
		Zu=Zabc:*u2sls // elementwise multiplication
		
		Wgmm=luinv(Zu'Zu) // inversion -> weight matrix
		bgmm=luinv(RX'Zabc*Wgmm*Zabc'RX)*(RX'Zabc*Wgmm*Zabc'yabc)
		ghat=Zabc'(yabc-RX*bgmm) // Empirical moment
		J=ghat'Wgmm*ghat // HS- J - statistic
		mR=mR-1 // degrees of freedom
		*b2sls and bgmm is already the post-Lasso estimator
		ZuuZ=Zu'Zu
		v2sls=luinv(RX'Zabc*izz*Zabc'RX)*RX'Zabc*izz*ZuuZ*izz*Zabc'RX*luinv(RX'Zabc*izz*Zabc'RX) // post estimator variance
	} // while end - continues until J smaller cutoff
alphaout = alpha[kx-1,.]
ncolalpha = cols(alphaout)
st_matrix("alphaout", alphaout)
st_numscalar("ncolalpha", ncolalpha)
st_matrix("ratio", ratio)
st_numscalar("m", m)
st_numscalar("n", n)
st_numscalar("Jinitial", Jinitial)
st_matrix("Zabc", Zabc)
st_matrix("xabc", xabc)
st_matrix("yabc", yabc)
} // mata end

loc ncolalpha = ncolalpha
loc wv
loc wi
token `ivs', parse("")
// In this loop create wv and wi locals
loc length_wi = 0
forval j = 1(1)`ncolalpha'{
	if alphaout[1,`j']==0 {
		local wv `wv' ``j''
	} 
		else {
		local wi `wi' ``j''
		local length_wi = `length_wi' + 1
	}
}

di "Valid shares: " "`wv' "
di "-----------------------------------------------------------"
di "Invalid shares: " "`wi'"
di "-----------------------------------------------------------"
di "Number of invalid shares: `length_wi'"
di "-----------------------------------------------------------"

// Restore dataset without if-restriction(s)
restore // to situation before weighting

di "Post-adaLasso 2SLS, classical SE:"
if "`weight'"!=""{
noi ivreg2 `lhs' `wi' `exog_xs' (`endog'=`wv') [`weight'=`exp'] `if', first `vce' partial(`exog_xs')
}
else {
noi ivreg2 `lhs' `exog_xs' `wi' (`endog'=`wv') `if', first `vce' partial(`exog_xs')
}
overid
local pfinal = r(p_oid)
noi di "Final P-value: `pfinal'"
matrix first = e(first)
loc fs = first[3,1]

use `tempssada', clear

// return local with products chosen as valid and invalid
ereturn local wv `wv'
ereturn local wi `wi'
ereturn local length_wi `length_wi'
 
end

cap program drop post_tabel
program define post_tabel

syntax varlist, depvar(varname) endog(varname) postmethod(string) method(string) tailn(real) wil(real) signi(real) vname(name)

ivregress `postmethod' `depvar' (`endog'  = year_share*) `varlist' i.stateyear i.city, cl(metarea) // 2SLS
est store `vname'_`method'_`postmethod'_`tailn'
estat firststage
matrix singleres = r(singleresults)
scalar fst = singleres[1,4]
estadd scalar fst
estadd local exiv "`varlist'"
estadd local length_wi "`wil'"
estadd local signi "`signi'"

end

cap program drop orig_tabel
program define orig_tabel

syntax, depvar(varname) endog(varname) postmethod(string) method(string) vname(name)

ivregress `postmethod' `depvar' (`endog'  = year_share*) i.stateyear i.city, cl(metarea) // 2SLS
est store `vname'_`method'_`postmethod'
estat firststage
matrix singleres = r(singleresults)
scalar fst = singleres[1,4]
estadd scalar fst
estadd local length_wi "0"
estadd local signi "-"

end

* Comments:
* This is used in 05BPforR.do
* I construct two IVs for each base year: one with lagged shifts, one with contemporaneous

cap program drop BP_adj_70
cap program drop BP_adj_80
cap program drop BP_adj_b

program define BP_adj_b
preserve
use rawdata/FB, clear

local vars = ""
cap local vars = "`2'"
di "`vars'"

forval sh = 1/19 {
	gen share70_aux_`sh' = share70_`sh'
}

forval sh = 1/19 {
	gen share80_aux_`sh' = share80_`sh'
}

if "`vars'"!="" {
foreach sharei in `3' {
	replace share70_aux_`sharei' = 0
}

foreach shares in `5' {
	replace share80_aux_`shares' = 0
}
} 
else {
 di "Empty!"  
}

merge 1:1 czone survey using data/aggregated_czone.dta, nogen

xtset czone decade

forval val = 1/19 {
	gen prod70_`val' = share70_aux_`val' * shift`val'
}
cap drop ss_1970
egen ss_1970 = rowtotal(prod70_*)

forval val = 1/19 {
	gen prod80_`val' = share80_aux_`val' * shift`val'
}
cap drop ss_1980
egen ss_1980 = rowtotal(prod80_*)

tabulate survey, gen(year)
gen lpop = l.population

ivreg2 `1' `2' `4' (d_immi_std  L.d_immi_std = ss_1970 ss_1980) year* [aw=l.population] if survey>=1990, robust cluster(czone)
underid
scalar fst = `r(j_uid)'
restore

end


program define BP_adj_70
preserve
use rawdata/FB, clear

forval sh = 1/19 {
	gen share70_aux_`sh' = share70_`sh'
	gen shift_aux`sh' = shift`sh'
}

foreach sharei in `3' {
	replace share70_aux_`sharei' = 0
	drop shift_aux`sharei'
}

forval val = 1/19 {
	gen prod70_`val' = share70_aux_`val' * shift`val'
}
cap drop ss_1970
egen ss_1970 = rowtotal(prod70_*)

merge 1:1 czone survey using data/aggregated_czone.dta, nogen

xtset czone decade

tabulate survey, gen(year)
gen lpop = l.population

ivreg2 `1' `2' (d_immi_std = ss_1970) year* [aw=l.population] if survey>=1990, robust cluster(czone)
weakivtest

restore

end

cap program drop shift_cor
program define shift_cor 
preserve
foreach sharei in `1' {
	drop shift`sharei'
}
keep survey shift*
keep in 2/5
reshape long shift, i(survey) j(id)
reshape wide shift, i(id) j(survey)
pwcorr shift*
estpost corr shift*, matrix
restore
end
