*
* 	Version 1, December 2019
* 	Nicolas Apfel
*


cap program drop sscim
prog sscim, eclass
version 10
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
loc exp = substr("`exp'",3,.)
replace `exp' = sqrt(`exp')
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

// Initial J 
qui ivregress gmm `lhs_fwl' (`endog_fwl'=`ivs_fwl'), robust
estat overid, forceweights
loc Jcf = e(J)
qui scalar tau=1-`c'/ln(`obs') // Target significance level
if `Jcf'<invchi2(`nuoverid',scalar(tau)) {
	noi dis in red "Hansen test does not reject in the first place. No evidence for invalid instruments. Use all shift-share products."
	restore
	preserve
	noi di as result  "Post-CIM 2SLS with shares as IVs, classical standard errors:"
	if "`weight'"!=""{
	  noi ivregress 2sls `lhs' `exog_xs' (`endog'=`ssstub'*) [`wts'], first vce(`vce')
	}
	else{
	  noi ivregress 2sls `lhs' `exog_xs' (`endog'=`ssstub'*), first vce(`vce')
	}
    
	qui estat firststage
	matrix singleres = r(singleresults)
	scalar fst = singleres[1,4]
	noi di as result "First-stage F-statistic: " fst 
	restore
	use `tempcim', clear
	exit 499
}
restore

	dis ""
	noi di as text "{hline 59}"
	noi dis as result "{txt}Confidence Interval Method"
	noi di as text "{hline 59}"

matrix define CoefSE = J(m, 2,.) // Matrix for coefs and SEs, m is number of IVs
local m = m // Nr of IVs as local

forval nam = 1/`m' { // make matrix with coefficients and standard errors
	preserve
	token `ivs', parse("")
	rename ``nam'' tx``nam''
		if "`weight'"!=""{
		 ivregress 2sls `lhs' `exog_xs' `ssstub'* (`endog'=tx``nam'') [`wts'], vce(`vce') // Important to control for other variables
		}
		else{
		 ivregress 2sls `lhs' `exog_xs' `ssstub'* (`endog'=tx``nam''), vce(`vce')
		}
		mat b = e(b)
		mat li b
		scalar D1 = b[1,1]
		mat V = e(V)
		mat li V
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
			forval u = 1/`m' {
			  if sel[1,`u'] == 1 {
			    loc add = sel[2,`u']
			    loc add = "`ssstub'`add'"
			    loc wv `wv' `add'
			  }
 			  if sel[1,`u'] == 0 {
			    loc add = sel[2,`u']
			    loc add = "`ssstub'`add'"
			    loc wi `wi' `add'
			  }
			}
			
			preserve // dataset, important to restore here, s.t. it can be restored in next iteration
			if "`weight'"!=""{
			  ivregress gmm `lhs' `exog_xs' `wi' (`endog'=`wv') [`wts'], robust
			}
			else{
			  ivregress gmm `lhs' `exog_xs' `wi' (`endog'=`wv'), robust
			}
			local Jc = e(J)  // J-statistic
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
		  cap drop i 
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
		  }
		}
		
		if "`wi'" != "`wiold'" { // skip if identity of invalid IVs did not change
		local wicount: word count `wi' // number of invalid IVs

		// J-statistic, data has been restore above
		  if "`weight'"!=""{
		    ivregress gmm `lhs' `exog_xs' `wi' (`endog'=`wv') [`wts'], robust
		  }
		  else{
		    ivregress gmm `lhs' `exog_xs' `wi' (`endog'=`wv'), robust
		  }
		  estat overid, forceweights
		  loc Jcf = r(HansenJ)
		  loc pcf = r(p_HansenJ)
		*}
		  noi di _column(26) "`wicount'"
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
    ivregress 2sls `lhs' `exog_xs' `wi' (`endog'=`wv') [`weight'=`exp'], first vce(`vce')
  }
  else {
    ivregress 2sls `lhs' `exog_xs' `wi' (`endog'=`wv'), first vce(`vce')
  }
  
  use `tempcim', clear
  
  ereturn local wv `wv'
  ereturn local wi `wi'
   
end
