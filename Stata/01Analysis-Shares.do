cap log close
* Set path
cd "C:\Users\na1d22\Dropbox\Surrey\IVSelection\Stata" // Adjust accordingly.

* Make programs available
run 00Programs-Final

log using Analysis-Shares.txt, text replace // Adjust the path

* Load the IV-data created before
use rawdata/FB.dta, clear // Creation of this and the following datasets detailed in Basso and Peri (2015), under https://giovanniperi.ucdavis.edu/data-for-basso-and-peri.html

* Merge with rest of the data which contains dependent and independent variables
merge 1:1 czone survey using data/aggregated_czone.dta, nogen // Same holds for this dataset.

* Merge country of origin-specific counts
xtset czone decade

* Create final IV
forval val = 1/19{
	gen prod70_`val' = share70_`val' * shift`val'
}
egen ss_1970 = rowtotal(prod70_*)

preserve
keep czone survey ss_1970
reshape wide ss_1970, i(czone) j(survey)
pwcorr ss_19701990 ss_19702000 ss_19702010
restore

preserve
keep survey shift*
keep in 2/5
reshape long shift, i(survey) j(id)
reshape wide shift, i(id) j(survey)
pwcorr shift*
estpost corr shift*, matrix
eststo corr_all
restore

forval v = 1/19 {
quietly corr share70_`v' share80_`v'
di "&", _c
di round(r(rho),0.001), _c
}

pwcorr share70_* share80_*

* Create year dummies
tabulate survey, gen(year)
gen lpop = l.population
gen ld_immi_std = l.d_immi_std

save data/analysis_data, replace

cap drop sum_shares
egen sum_shares = rowtotal(share70*)

ivreg2 dlweekly_hskill (d_immi_std= share70_*) year* [aw=l.population] if survey>=1990, robust cluster(czone) nocons
weakivtest

ivreg2 dlweekly_hskill (d_immi_std= share70_*) year* [aw=l.population] if survey>=1990, robust cluster(czone) nocons
weakivtest
scalar fst = r(F_eff)
di r(F_eff)
return list

ivreg2 dlweekly_lskill (d_immi_std = ss_1970) year* [aw=l.population] if survey>=1990, robust cluster(czone) nocons
weakivtest
di r(F_eff)

*****************
* MAIN ANALYSIS *
*****************
quietly{
est drop _all
foreach fb in 70 { // Basis year of the IV
	foreach dv in "dlweekly" "dlweekly_hskill" "dlweekly_lskill" { // Dependent variable used
loc dvn
loc fb = `fb'
* Create abbreviations
if "`dv'" == "dlweekly"{
	loc dvn "w"
}
if "`dv'" == "dlweekly_hskill"{
	loc dvn "wh"
}
if "`dv'" == "dlweekly_lskill"{
	loc dvn "wl"
}

use data/analysis_data, replace
* Standard: 1
* All shares are assumed to be valid
**** OLS
xi: reg `dv' d_immi_std i.survey [aw=l.population] if survey>=1990, vce(cluster czone)
est store `dvn'_fb`fb'_1_ols
**** 2SLS
ivreg2 `dv' (d_immi_std= share`fb'_*) year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store `dvn'_fb`fb'_1_tsls
eststo: weakivtest
scalar fst = r(F_eff)
scalar c5TSLS = r(c_TSLS_5)
scalar c10TSLS = r(c_TSLS_10)
scalar c20TSLS = r(c_TSLS_20)
scalar c5LIML = r(c_LIML_5)
scalar c10LIML = r(c_LIML_10)
scalar c20LIML = r(c_LIML_20)
est restore `dvn'_fb`fb'_1_tsls
estadd scalar fst
estadd scalar c5TSLS
estadd scalar c10TSLS
estadd scalar c20TSLS
estadd local nri "0"
estadd local signi "-"
**** LIML
ivregress liml `dv' (d_immi_std= share`fb'_*) year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store `dvn'_fb`fb'_1_liml
estadd scalar fst
estadd scalar c5LIML
estadd scalar c10LIML
estadd scalar c20LIML
estadd local nri "0"
estadd local signi "-"
**** Bartik
ivreg2 `dv' (d_immi_std= ss_1970) year* [aw=l.population]  if survey>=1990, robust cluster(czone)
weakivtest
scalar fst = r(F_eff)
est store `dvn'_fb`fb'_1_ss
scalar c5TSLS = r(c_TSLS_5)
scalar c10TSLS = r(c_TSLS_10)
scalar c20TSLS = r(c_TSLS_20)
est restore `dvn'_fb`fb'_1_ss
estadd scalar fst
estadd scalar c5TSLS
estadd scalar c10TSLS
estadd scalar c20TSLS

foreach cval in 0.1 0.384 0.768 { // Two different
*SSADA: 2
* Selection via aLasso
di "****************************************"
di "*            ADALASSO                  *"
di "****************************************"

loc cabk
if `cval' == 0.1{
loc cabk = 0
loc signi = "0.01302"
}

if `cval' == 0.384{
loc cabk = 1
loc signi = "0.05"
}
if `cval' == 0.768 {
loc cabk = 4
loc signi = "0.1"
}

cap ssada `dv' d_immi_std year* if survey >= 1990 [aw=lpop], endog(d_immi_std) exog(share`fb'_* year*) id(czone survey) vce(cluster czone) c(`cval')
loc wv `e(wv)'
loc wi `e(wi)'
di `wv'
di _rc
if _rc == 0 {
**** 2SLS
ivreg2 `dv' `wi' (d_immi_std= `wv') year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store `dvn'_fb`fb'_2_tsls_`cabk'
weakivtest
scalar fst = r(F_eff)
scalar c5TSLS = r(c_TSLS_5)
scalar c10TSLS = r(c_TSLS_10)
scalar c20TSLS = r(c_TSLS_20)
scalar c5LIML = r(c_LIML_5)
scalar c10LIML = r(c_LIML_10)
scalar c20LIML = r(c_LIML_20)
est restore `dvn'_fb`fb'_2_tsls_`cabk'
estadd scalar fst
estadd scalar c5TSLS
estadd scalar c10TSLS
estadd scalar c20TSLS
local ninv: word count `wi'
estadd local nri "`ninv'"
estadd local signi `signi'
noi di `signi'

**** LIML
ivregress liml `dv' `wi' (d_immi_std= `wv') year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store `dvn'_fb`fb'_2_liml_`cabk'
estadd scalar fst
estadd scalar c5LIML
estadd scalar c10LIML
estadd scalar c20LIML
estadd local nri "`ninv'"
estadd local signi `signi'

loc line "`dv' & aLasso & `signi'"
foreach var of varlist share70_* {
if (strpos("`wi' ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

**** BP-Bartik
loc invind
foreach share in `wi' {
	loc shareind = substr("`share'", 9, .)
	loc invind `invind' `shareind'
}
BP_adj_`fb' "`dv'" "`wi'" "`invind'" "`dvn'_al_hs"
scalar fst = r(F_eff)
est store `dvn'_fb`fb'_2_ss_`cabk'
scalar c5TSLS = r(c_TSLS_5)
scalar c10TSLS = r(c_TSLS_10)
scalar c20TSLS = r(c_TSLS_20)
est restore `dvn'_fb`fb'_2_ss_`cabk'
estadd scalar fst
estadd scalar c5TSLS
estadd scalar c10TSLS
estadd scalar c20TSLS
} 
else {
est restore `dvn'_fb`fb'_1_tsls
est store `dvn'_fb`fb'_2_tsls_`cabk'
estadd local signi `signi', replace
est restore `dvn'_fb`fb'_1_liml
est store `dvn'_fb`fb'_2_liml_`cabk'
estadd local signi `signi', replace
est restore `dvn'_fb`fb'_1_ss
est store `dvn'_fb`fb'_2_ss_`cabk'

loc line = "`dv' & ADA & `signi' & - & - & - & - & - & - & - & - & - & - & - & - & - & - & - & - & - & - & -"
noi di "`line'"
}

*SSCIM: 3
* Selection via CIM
di "****************************************"
di "*              SSCIM                   *"
di "****************************************"
cap sscim `dv' d_immi_std year* if survey >= 1990 [aw=lpop], endog(d_immi_std) exog(share`fb'_* year*) ssstub(share`fb'_) vce("cluster(czone)") c(`cval')
loc wv `e(wv)'
loc wi `e(wi)'
di `wv'
di _rc
if _rc == 0 {
**** 2SLS
ivreg2 `dv' `wi' (d_immi_std= `wv') year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store `dvn'_fb`fb'_3_tsls_`cabk'
weakivtest
scalar fst = r(F_eff)
scalar c5TSLS = r(c_TSLS_5)
scalar c10TSLS = r(c_TSLS_10)
scalar c20TSLS = r(c_TSLS_20)
scalar c5LIML = r(c_LIML_5)
scalar c10LIML = r(c_LIML_10)
scalar c20LIML = r(c_LIML_20)
est restore `dvn'_fb`fb'_3_tsls_`cabk'
estadd scalar fst
estadd scalar c5TSLS
estadd scalar c10TSLS
estadd scalar c20TSLS
local ninv: word count `wi'
estadd local nri "`ninv'"
estadd local signi `signi'
noi di "`ninv'"

**** LIML
ivregress liml `dv' `wi' (d_immi_std= `wv') year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store `dvn'_fb`fb'_3_liml_`cabk'
estadd scalar fst
estadd scalar c5LIML
estadd scalar c10LIML
estadd scalar c20LIML
estadd local nri "`ninv'"
estadd local signi `signi'
**** BP-Bartik
loc invind
foreach share in `wi' {
	loc shareind = substr("`share'", 9, .)
	loc invind `invind' `shareind'
}
noi BP_adj_`fb' "`dv'" "`wi'" "`invind'" "`dvn'_`cabk'_cim_hs"
scalar fst = r(F_eff)
estadd local method "CIM"
estadd local test "HS"
est store `dvn'_fb`fb'_3_ss_`cabk'
scalar c5TSLS = r(c_TSLS_5)
scalar c10TSLS = r(c_TSLS_10)
scalar c20TSLS = r(c_TSLS_20)
est restore `dvn'_fb`fb'_3_ss_`cabk'
estadd scalar fst
estadd scalar c5TSLS
estadd scalar c10TSLS
estadd scalar c20TSLS
} 
else {
est restore `dvn'_fb`fb'_1_tsls
est store `dvn'_fb`fb'_3_tsls_`cabk'
est restore `dvn'_fb`fb'_1_liml
est store `dvn'_fb`fb'_3_liml_`cabk'
estadd local signi `signi', replace
est restore `dvn'_fb`fb'_1_ss
est store `dvn'_fb`fb'_3_ss_`cabk'
} // end last else

loc line "`dv' & CIM & `signi'"
foreach var of varlist share70_* {
if (strpos("`wi' ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

} // end foreach cval
} // end foreach dv
} // end foreach fb
}

* Load Data
use data/analysis_data, clear

*********************************************
* I. SINGLE REGRESSOR Selection via AR-test *
*********************************************
* Selection happens in 04BP-shares.R

* Single regressor case
* Both aLasso and CIM

* 1. aLasso
quietly{
	* dlweekly
foreach var of varlist share70_1  share70_3  share70_5  share70_6  share70_8  share70_9 share70_16 share70_17 share70_19  {
	ren `var' k_`var'
}

ivreg2 dlweekly year3 year4 k_* (d_immi_std = share70_*) [aw=lpop] if survey >=1990, robust cluster(czone)
weakivtest
scalar fst = r(F_eff)
est store w_fb70_2_tsls_ar
scalar c5TSLS = r(c_TSLS_5)
scalar c10TSLS = r(c_TSLS_10)
scalar c20TSLS = r(c_TSLS_20)
scalar c5LIML = r(c_LIML_5)
scalar c10LIML = r(c_LIML_10)
scalar c20LIML = r(c_LIML_20)
est restore w_fb70_2_tsls_ar
estadd scalar fst
estadd scalar c5TSLS
estadd scalar c10TSLS
estadd scalar c20TSLS
estadd local nri "9"
estadd local signi "0.01302"
ivregress liml dlweekly year3 year4 k_* (d_immi_std = share70_*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store w_fb70_2_liml_ar
estadd scalar fst
estadd scalar c5LIML
estadd scalar c10LIML
estadd scalar c20LIML
estadd local nri "9"
estadd local signi "0.01302"

use data/analysis_data, clear
BP_adj_70 dlweekly "share70_1  share70_3  share70_5  share70_6  share70_8  share70_9 share70_16 share70_17 share70_19" "1 3 5 6 8 9 16 17 19" 
scalar fst = r(F_eff)
est store w_fb70_2_ss_ar
scalar c5TSLS = r(c_TSLS_5)
scalar c10TSLS = r(c_TSLS_10)
scalar c20TSLS = r(c_TSLS_20)
est restore w_fb70_2_ss_ar
estadd scalar fst
estadd scalar c5TSLS
estadd scalar c10TSLS
estadd scalar c20TSLS
shift_cor "1 3 5 6 8 9 16 17 19 "
eststo corr_w_al_ar
estadd local method = "AL"
estadd local test = "AR"

use data/analysis_data, clear
loc line "dlweekly & aLasso AR & 0.01302"
foreach var of varlist share70_* {
if (strpos("share70_1  share70_3  share70_5  share70_6  share70_8  share70_9 share70_16 share70_17 share70_19 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

	* dlweekly_hskill
use data/analysis_data, clear
foreach var of varlist share70_2  share70_3  share70_5  share70_6  share70_8  share70_9 share70_13 share70_14 share70_16 share70_17 share70_19 {
	ren `var' k_`var'
}

ivreg2 dlweekly_hskill year3 year4 k_* (d_immi_std = share70_*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store wh_fb70_2_tsls_ar
weakivtest
scalar fst = r(F_eff)
scalar c5TSLS = r(c_TSLS_5)
scalar c10TSLS = r(c_TSLS_10)
scalar c20TSLS = r(c_TSLS_20)
scalar c5LIML = r(c_LIML_5)
scalar c10LIML = r(c_LIML_10)
scalar c20LIML = r(c_LIML_20)
est restore wh_fb70_2_tsls_ar
estadd scalar fst
estadd scalar c5TSLS
estadd scalar c10TSLS
estadd scalar c20TSLS
estadd local nri "11"
estadd local signi "0.01302"
ivregress liml dlweekly_hskill year3 year4 k_* (d_immi_std = share70_*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store wh_fb70_2_liml_ar
estadd scalar fst
estadd scalar c5LIML
estadd scalar c10LIML
estadd scalar c20LIML
estadd local nri "11"
estadd local signi "0.01302"
estat overid, forceweights forcenonrobust

use data/analysis_data, clear
BP_adj_70 dlweekly_hskill "share70_2  share70_3  share70_5  share70_6  share70_8  share70_9 share70_13 share70_14 share70_16 share70_17 share70_19 " "2 3 5 6 8 9 13 14 16 17 19"
scalar fst = r(F_eff)
est store wh_fb70_2_ss_ar
scalar c5TSLS = r(c_TSLS_5)
scalar c10TSLS = r(c_TSLS_10)
scalar c20TSLS = r(c_TSLS_20)
est restore wh_fb70_2_ss_ar
estadd scalar fst
estadd scalar c5TSLS
estadd scalar c10TSLS
estadd scalar c20TSLS
estadd local nri "11"
estadd local signi "0.01302"
shift_cor "2 3 5 6 8 9 13 14 16 17 19"
eststo corr_w_al_ar
estadd local method = "AL"
estadd local test = "AR"

use data/analysis_data, clear
loc line "dlweekly\_hskill & aLasso AR & 0.01302"
foreach var of varlist share70_* {
if (strpos("share70_2  share70_3  share70_5  share70_6  share70_8  share70_9 share70_13 share70_14 share70_16 share70_17 share70_19 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

	* dlweekly_lskill
use data/analysis_data, clear
foreach var of varlist share70_1 share70_5 share70_6 share70_7 share70_8 share70_9 share70_13 share70_17 share70_19 {
	ren `var' k_`var'
}

ivreg2 dlweekly_lskill year3 year4 k_* (d_immi_std = share70_*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store wl_fb70_2_tsls_ar
weakivtest
scalar fst = r(F_eff)
scalar c5TSLS = r(c_TSLS_5)
scalar c10TSLS = r(c_TSLS_10)
scalar c20TSLS = r(c_TSLS_20)
scalar c5LIML = r(c_LIML_5)
scalar c10LIML = r(c_LIML_10)
scalar c20LIML = r(c_LIML_20)
est restore wl_fb70_2_tsls_ar
estadd scalar fst
estadd scalar c5TSLS
estadd scalar c10TSLS
estadd scalar c20TSLS
estadd local nri "9"
estadd local signi "0.01302"
ivregress liml dlweekly_lskill year3 year4 k_* (d_immi_std = share70_*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store wl_fb70_2_liml_ar
estadd scalar fst
estadd scalar c5LIML
estadd scalar c10LIML
estadd scalar c20LIML
estadd local nri "9"
estadd local signi "0.01302"
estat overid, forceweights forcenonrobust

use data/analysis_data, clear
BP_adj_70 dlweekly_lskill "share70_1  share70_5  share70_6  share70_7  share70_8  share70_9 share70_13 share70_17 share70_19 " "1 5 6 7 8 9 13 17 19"
scalar fst = r(F_eff)
est store wl_fb70_2_ss_ar
scalar c5TSLS = r(c_TSLS_5)
scalar c10TSLS = r(c_TSLS_10)
scalar c20TSLS = r(c_TSLS_20)
est restore wl_fb70_2_ss_ar
estadd scalar fst
estadd scalar c5TSLS
estadd scalar c10TSLS
estadd scalar c20TSLS
estadd local nri "9"
estadd local signi "0.01302"
shift_cor "1 5 6 7 8 9 13 17 19"
eststo corr_wl_al_ar
estadd local method = "AL"
estadd local test = "AR"


use data/analysis_data, clear
loc line "dlweekly\_lskill & aLasso AR & 0.01302"
foreach var of varlist share70_* {
if (strpos("share70_1  share70_5  share70_6  share70_7  share70_8  share70_9 share70_13 share70_17 share70_19 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"
}

* 2. Confidence Interval Method

foreach dv in "dlweekly" "dlweekly_hskill" "dlweekly_lskill" { // Dependent variable used
loc dvn
quietly{
* Create abbreviations
if "`dv'" == "dlweekly"{
	loc dvn "w"
}
if "`dv'" == "dlweekly_hskill"{
	loc dvn "wh"
}
if "`dv'" == "dlweekly_lskill"{
	loc dvn "wl"
}

use data/analysis_data, clear

* Selection via CIM
di "****************************************"
di "*              SSCIM                   *"
di "****************************************"
cap sscimliml `dv' d_immi_std year* if survey >= 1990 [aw=lpop], endog(d_immi_std) exog(share70_* year*) ssstub(share70_) vce("cluster(czone)") c(0.1)
loc wv `e(wv)'
loc wi `e(wi)'
di _rc

**** 2SLS
ivreg2 `dv' `wi' (d_immi_std= `wv') year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store `dvn'_fb70_3_tsls_ar
weakivtest
scalar fst = r(F_eff)
scalar c5TSLS = r(c_TSLS_5)
scalar c10TSLS = r(c_TSLS_10)
scalar c20TSLS = r(c_TSLS_20)
scalar c5LIML = r(c_LIML_5)
scalar c10LIML = r(c_LIML_10)
scalar c20LIML = r(c_LIML_20)
est restore `dvn'_fb70_3_tsls_ar
estadd scalar fst
estadd scalar c5TSLS
estadd scalar c10TSLS
estadd scalar c20TSLS
local ninv: word count `wi'
estadd local nri "`ninv'"
estadd local signi "0.01302"

**** BP-Bartik
loc invind
foreach share in `wi' {
	loc shareind = substr("`share'", 9, .)
	loc invind `invind' `shareind'
}
noi BP_adj_70 "`dv'" "`wi'" "`invind'"
scalar fst = r(F_eff)
est store `dvn'_fb70_3_ss_ar
scalar c5TSLS = r(c_TSLS_5)
scalar c10TSLS = r(c_TSLS_10)
scalar c20TSLS = r(c_TSLS_20)
est restore `dvn'_fb70_3_ss_ar
estadd scalar fst
estadd scalar c5TSLS
estadd scalar c10TSLS
estadd scalar c20TSLS
shift_cor "`invind'"
eststo corr_`dvn'_cim_ar
estadd local method = "CIM"
estadd local test = "AR"

**** LIML
ivreg2 `dv' `wi' (d_immi_std= `wv') year* [aw=l.population] if survey>=1990, robust cluster(czone) liml first
est store `dvn'_fb70_3_liml_ar
estadd scalar fst
estadd scalar c5LIML
estadd scalar c10LIML
estadd scalar c20LIML
estadd local nri "`ninv'"
estadd local signi "0.01302"
noi di "`dv'"
noi di "`wi'"

use data/analysis_data, clear
loc line "`dv' & CIM AR & 0.01302"
foreach var of varlist share70_* {
if (strpos("`wi' ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

} 
}

****************************
* II. MULTIPLE REGRESSORS  *
****************************
* w_mult_orig_tsls
* `dv'_mult_`selectionmethod'_`postselectionmethod'

* Selection is made in R - 04BP-shares.R

	* A. dlweekly
	**************
* 1. Original results dlweekly
quietly{
*cap est drop _all
use data/analysis_data, clear

BP_adj_b dlweekly
est store w_mult_orig_ss
estadd scalar fst

ivreg2 dlweekly year3 year4 (d_immi_std L.d_immi_std =share70* share80*) [aw=lpop] if survey >=1990, cluster(czone)
underid
scalar fst = `r(j_uid)'
est store w_mult_orig_tsls
estadd scalar fst
est store w_mult_orig_tsls_short
estadd local nri "0"
estadd local signi "-"
ivreg2 dlweekly year3 year4 (d_immi_std L.d_immi_std =share70* share80*) [aw=lpop] if survey >=1990, robust cluster(czone) liml
underid
scalar fst = `r(j_uid)'
est store w_mult_orig_liml
estadd scalar fst
estadd local nri "0"
estadd local signi "-"

* 2. AL HS selection - dlweekly
est restore w_mult_orig_ss
est store w_mult_hs_ss
est restore w_mult_orig_tsls
est store w_mult_hs_tsls
est restore w_mult_orig_tsls
est store w_mult_hs_tsls_short
estadd local nri "0", replace
estadd local signi "0.1", replace
est restore w_mult_orig_liml
est store w_mult_hs_liml
estadd local nri "0", replace
estadd local signi "0.1", replace

* 3. AL AR selection - dlweekly
use data/analysis_data, clear
foreach var of varlist share80_1  share80_2  share80_4  share80_5 share80_10 share80_13 share80_15 share70_17 share70_19 share80_19 {
	ren `var' k_`var'
}

ivreg2 dlweekly year3 year4 k_* (d_immi_std L.d_immi_std =share70* share80*) [aw=lpop] if survey >=1990, robust cluster(czone)
underid
scalar fst = `r(j_uid)'
est store w_mult_ar_tsls
estadd scalar fst
est store w_mult_ar_tsls_short
estadd local nri "10"
estadd local signi "0.01302"
ivreg2 dlweekly year3 year4 k_* (d_immi_std L.d_immi_std =share70* share80*) [aw=lpop] if survey >=1990, robust cluster(czone) liml
underid
scalar fst = `r(j_uid)'
est store w_mult_ar_liml
estadd scalar fst
estadd local nri "10"
estadd local signi "0.01302"
di `e(arubinp)'

use data/analysis_data, clear
BP_adj_b dlweekly "share70_17 share70_19" "17 19" "share80_1  share80_2  share80_4  share80_5 share80_10 share80_13 share80_15 share80_19 " "1 2 4 5 10 13 15 19"
est store w_mult_ar_ss
estadd scalar fst

use data/analysis_data, clear
loc line "dlweekly & aLasso AR & 1970"
foreach var of varlist share70_* {
if (strpos("share70_17 share70_19 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

use data/analysis_data, clear
loc line " &  & 1980"
foreach var of varlist share80_* {
if (strpos("share80_1 share80_2 share80_4 share80_5 share80_10 share80_13 share80_15 share80_19 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

* 4. Original results dlweekly_hskill
use data/analysis_data, clear

BP_adj_b dlweekly_hskill
est store wh_mult_orig_ss
estadd scalar fst

ivreg2 dlweekly_hskill year3 year4 (d_immi_std L.d_immi_std =share70* share80*) [aw=lpop] if survey >=1990, cluster(czone)
underid
scalar fst = `r(j_uid)'
est store wh_mult_orig_tsls
estadd scalar fst
est store wh_mult_orig_tsls_short
estadd local nri "0"
estadd local signi "-"
ivreg2 dlweekly_hskill year3 year4 (d_immi_std L.d_immi_std =share70* share80*) [aw=lpop] if survey >=1990, robust cluster(czone) liml
underid
scalar fst = `r(j_uid)'
est store wh_mult_orig_liml
estadd scalar fst
estadd local nri "0"
estadd local signi "-"

* 5. AL HS selection - dlweekly_hskill
use data/analysis_data, clear
foreach var of varlist share70_5 share80_5 {
	ren `var' k_`var'
}

ivreg2 dlweekly_hskill year3 year4 k_* (d_immi_std L.d_immi_std =share70* share80*) [aw=lpop] if survey >=1990, robust cluster(czone)
underid
scalar fst = `r(j_uid)'
est store wh_mult_hs_tsls
estadd scalar fst
est store wh_mult_hs_tsls_short
estadd local nri "2"
estadd local signi "0.1"
ivreg2 dlweekly_hskill year3 year4 k_* (d_immi_std L.d_immi_std =share70* share80*) [aw=lpop] if survey >=1990, robust cluster(czone) liml
underid
scalar fst = `r(j_uid)'
est store wh_mult_hs_liml
estadd scalar fst
estadd local nri "2"
estadd local signi "0.1"

use data/analysis_data, clear
BP_adj_b dlweekly_hskill "share70_5" "5" "share80_5" "5"
est store wh_mult_hs_ss
estadd scalar fst

use data/analysis_data, clear
loc line "dlweekly\_hskill & aLasso HS & 1970"
foreach var of varlist share70_* {
if (strpos("share70_5 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

use data/analysis_data, clear
loc line " &  & 1980"
foreach var of varlist share80_* {
if (strpos("share80_5 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

* 6. AL AR selection - dlweekly_hskill
use data/analysis_data, clear
foreach var of varlist share80_1  share80_2  share80_3  share80_4  share70_5  share80_5 share80_11 share80_13 share80_15 share70_17 share70_19 {
	ren `var' k_`var'
}

ivreg2 dlweekly_hskill year3 year4 k_* (d_immi_std L.d_immi_std =share70* share80*) [aw=lpop] if survey >=1990, robust cluster(czone)
underid
scalar fst = `r(j_uid)'
est store wh_mult_ar_tsls
estadd scalar fst
est store wh_mult_ar_tsls_short
estadd local nri "11"
estadd local signi "0.01302"
ivreg2 dlweekly_hskill year3 year4 k_* (d_immi_std L.d_immi_std =share70* share80*) [aw=lpop] if survey >=1990, robust cluster(czone) liml
underid
scalar fst = `r(j_uid)'
est store wh_mult_ar_liml
estadd scalar fst
estadd local nri "11"
estadd local signi "0.01302"
di `e(arubinp)'

use data/analysis_data, clear
BP_adj_b dlweekly_hskill "share70_5 share70_17 share70_19" "5 17 19" "share80_1 share80_2 share80_3 share80_4 share80_5 share80_11 share80_13 share80_15" "1 2 4 5 11 13 15"
est store wh_mult_ar_ss
estadd scalar fst

use data/analysis_data, clear
loc line "dlweekly_\hskill & aLasso AR & 1970"
foreach var of varlist share70_* {
if (strpos("share70_5 share70_17 share70_19 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

use data/analysis_data, clear
loc line " &  & 1980"
foreach var of varlist share80_* {
if (strpos("share80_1 share80_2 share80_3 share80_4 share80_5 share80_11 share80_13 share80_15 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

* 7. Original results dlweekly_lskill
use data/analysis_data, clear

BP_adj_b dlweekly_lskill
est store wl_mult_orig_ss
estadd scalar fst

ivreg2 dlweekly_lskill year3 year4 (d_immi_std L.d_immi_std =share70* share80*) [aw=lpop] if survey >=1990, cluster(czone)
underid
scalar fst = `r(j_uid)'
est store wl_mult_orig_tsls
estadd scalar fst
est store wl_mult_orig_tsls_short
estadd local nri "0"
estadd local signi "-"
ivreg2 dlweekly_lskill year3 year4 (d_immi_std L.d_immi_std =share70* share80*) [aw=lpop] if survey >=1990, robust cluster(czone) liml
underid
scalar fst = `r(j_uid)'
est store wl_mult_orig_liml
estadd scalar fst
estadd local nri "0"
estadd local signi "-"

* 8. AL HS selection - dlweekly_lskill
est restore wl_mult_orig_ss
est store wl_mult_hs_ss
est restore wl_mult_orig_tsls
est store wl_mult_hs_tsls
est store wl_mult_hs_tsls_short
estadd local nri "0", replace
estadd local signi "0.1", replace
est restore wl_mult_orig_liml
est store wl_mult_hs_liml
estadd local nri "0", replace
estadd local signi "0.1", replace

* 9. AL AR selection - dlweekly_lskill
use data/analysis_data, clear
foreach var of varlist share80_1 share70_2 share80_2 share80_4 share80_5 share80_6 share80_16 share70_17 share70_19 {
	ren `var' k_`var'
}

ivreg2 dlweekly_lskill year3 year4 k_* (d_immi_std L.d_immi_std =share70* share80*) [aw=lpop] if survey >=1990, robust cluster(czone)
underid
scalar fst = `r(j_uid)'
est store wl_mult_ar_tsls
estadd scalar fst
est store wl_mult_ar_tsls_short
estadd local nri "9"
estadd local signi "0.01302"
ivreg2 dlweekly_lskill year3 year4 k_* (d_immi_std L.d_immi_std =share70* share80*) [aw=lpop] if survey >=1990, robust cluster(czone) liml
underid
scalar fst = `r(j_uid)'
est store wl_mult_ar_liml
estadd scalar fst
estadd local nri "9"
estadd local signi "0.01302"
di `e(arubinp)'

use data/analysis_data, clear
BP_adj_b dlweekly_lskill "share70_2 share70_17 share70_19" "2 17 19" "share80_2 share80_4 share80_5 share80_6 share80_16" "2 4 5 6 16"
est store wl_mult_ar_ss
estadd scalar fst

use data/analysis_data, clear
loc line "dlweekly_\lskill & aLasso AR & 1970"
foreach var of varlist share70_* {
if (strpos("share70_2 share70_17 share70_19 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

use data/analysis_data, clear
loc line " &  & 1980"
foreach var of varlist share80_* {
if (strpos("share80_1 share80_2 share80_4 share80_5 share80_6 share80_16 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

}

****************************
* 	   CREATE TABLES	   *
****************************
* Including these subtables in the order displayed here creates the tables in the paper

*All results: For dissertation
*Wages
esttab w_fb70_1_ss w_fb70_2_ss_1 w_fb70_3_ss_1 w_fb70_2_ss_4 w_fb70_3_ss_4 w_fb70_2_ss_ar  w_fb70_3_ss_ar using tables/w_ss_all.tex, keep(d_immi_std) scalar("fst F") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth ) mtitles("Standard" "AL (HS)" "CIM (HS)" "AL (HS)" "CIM (HS)" "AL (AR)" "CIM (AR)") ///
	coeflabels(d_immi_std "\toprule \multicolumn{8}{c}{\textit{Panel A}: Change in average log weekly wages} \\ \midrule SSIV") nolines star(+ 0.10 * 0.05 ** 0.01 *** 0.001) nonotes booktabs nogaps replace
label variable d_immi_std "2SLS"
esttab w_fb70_1_tsls w_fb70_2_tsls_1 w_fb70_3_tsls_1 w_fb70_2_tsls_4 w_fb70_3_tsls_4 w_fb70_2_tsls_ar  w_fb70_3_tsls_ar using tables/w_tsls_all.tex, keep(d_immi_std) scalar("fst F") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth ) nomtitles nonumbers nolines nonotes posthead(\midrule) star(+ 0.10 * 0.05 ** 0.01 *** 0.001) booktabs nogaps replace 
label variable d_immi_std "LIML"
esttab w_fb70_1_liml w_fb70_2_liml_1 w_fb70_3_liml_1 w_fb70_2_liml_4 w_fb70_3_liml_4 w_fb70_2_liml_ar  w_fb70_3_liml_ar using tables/w_liml_all.tex, keep(d_immi_std) scalar("fst F" "nri \midrule \# inv" "signi Sign.") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth ) nomtitles star(+ 0.10 * 0.05 ** 0.01 *** 0.001) nonumbers nolines nonotes posthead(\midrule) booktabs nogaps replace 
	
	* Wages, HS
esttab wh_fb70_1_ss wh_fb70_2_ss_1 wh_fb70_3_ss_1 wh_fb70_2_ss_4 wh_fb70_3_ss_4 wh_fb70_2_ss_ar  wh_fb70_3_ss_ar using tables/wh_ss_all.tex, keep(d_immi_std) scalar("fst F") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth ) nomtitles ///
	coeflabels(d_immi_std "\midrule \multicolumn{8}{c}{\textit{Panel B}: Change in avg. log weekly wages of high-skilled} \\ \midrule SSIV") nolines star(+ 0.10 * 0.05 ** 0.01 *** 0.001) nonumbers nonotes booktabs nogaps replace
label variable d_immi_std "2SLS"
esttab wh_fb70_1_tsls wh_fb70_2_tsls_1 wh_fb70_3_tsls_1 wh_fb70_2_tsls_4 wh_fb70_3_tsls_4 wh_fb70_2_tsls_ar  wh_fb70_3_tsls_ar using tables/wh_tsls_all.tex, keep(d_immi_std) scalar("fst F") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth ) nomtitles star(+ 0.10 * 0.05 ** 0.01 *** 0.001) nonumbers nolines nonotes posthead(\midrule) booktabs nogaps replace 
label variable d_immi_std "LIML"
esttab wh_fb70_1_liml wh_fb70_2_liml_1 wh_fb70_3_liml_1 wh_fb70_2_liml_4 wh_fb70_3_liml_4 wh_fb70_2_liml_ar  wh_fb70_3_liml_ar using tables/wh_liml_all.tex, keep(d_immi_std) scalar("fst F" "nri \midrule \# inv" "signi Sign.") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth ) nomtitles star(+ 0.10 * 0.05 ** 0.01 *** 0.001) nonumbers nolines nonotes posthead(\midrule) booktabs nogaps replace 

	* Wages, LS
esttab wl_fb70_1_ss wl_fb70_2_ss_1 wl_fb70_3_ss_1 wl_fb70_2_ss_4 wl_fb70_3_ss_4 wl_fb70_2_ss_ar  wl_fb70_3_ss_ar using tables/wl_ss_all.tex, keep(d_immi_std) scalar("fst F") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth ) nomtitles ///
	star(+ 0.10 * 0.05 ** 0.01 *** 0.001) coeflabels(d_immi_std "\midrule \multicolumn{8}{c}{\textit{Panel C}: Change in avg. log weekly wages of low-skilled} \\ \midrule SSIV") nolines nonumbers nonotes booktabs nogaps replace
label variable d_immi_std "2SLS"
esttab wl_fb70_1_tsls wl_fb70_2_tsls_1 wl_fb70_3_tsls_1 wl_fb70_2_tsls_4 wl_fb70_3_tsls_4 wl_fb70_2_tsls_ar  wl_fb70_3_tsls_ar using tables/wl_tsls_all.tex, keep(d_immi_std) scalar("fst F") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth ) nomtitles star(+ 0.10 * 0.05 ** 0.01 *** 0.001) nonumbers nolines nonotes posthead(\midrule) booktabs nogaps replace 
label variable d_immi_std "LIML"
esttab wl_fb70_1_liml wl_fb70_2_liml_1 wl_fb70_3_liml_1 wl_fb70_2_liml_4 wl_fb70_3_liml_4 wl_fb70_2_liml_ar  wl_fb70_3_liml_ar using tables/wl_liml_all.tex, keep(d_immi_std) scalar("fst F" "nri \midrule \# inv" "signi Sign.") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth ) nomtitles star(+ 0.10 * 0.05 ** 0.01 *** 0.001) nonumbers nolines nonotes posthead(\midrule) ///
	postfoot(\bottomrule  ///
	\end{tabular*} }) booktabs nogaps replace
	
*For paper:
* Wages
esttab w_fb70_1_ss w_fb70_2_ss_4 w_fb70_3_ss_4 w_fb70_2_ss_ar  w_fb70_3_ss_ar using tables/w_ss.tex, keep(d_immi_std) scalar("fst F" "c5TSLS $\tau=0.05$" "c10TSLS $\tau=0.1$" "c20TSLS $\tau=0.2$") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth) mtitles("Standard" "AL HS" "CIM HS" "AL AR" "CIM AR") nonumbers ///
	coeflabels(d_immi_std "\toprule \multicolumn{6}{c}{\textit{Panel A}: Change in average log weekly wages} \\ \midrule SSIV") nolines star nonotes booktabs nogaps replace
label variable d_immi_std "2SLS"
esttab w_fb70_1_tsls w_fb70_2_tsls_4 w_fb70_3_tsls_4 w_fb70_2_tsls_ar  w_fb70_3_tsls_ar using tables/w_tsls.tex, keep(d_immi_std) scalar("fst F" "c5TSLS $\tau=0.05$" "c10TSLS $\tau=0.1$" "c20TSLS $\tau=0.2$") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth) nomtitles nonumbers nolines nonotes posthead(\midrule) star booktabs nogaps replace 
label variable d_immi_std "LIML"
esttab w_fb70_1_liml w_fb70_2_liml_4 w_fb70_3_liml_4 w_fb70_2_liml_ar  w_fb70_3_liml_ar using tables/w_liml.tex, keep(d_immi_std) scalar("fst F" "c5LIML $\tau=0.05$" "c10LIML $\tau=0.1$" "c20LIML $\tau=0.2$" "nri \midrule \# inv") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth) nomtitles star nonumbers nolines nonotes posthead(\midrule) booktabs nogaps replace 
	
	* Wages, HS
esttab wh_fb70_1_ss wh_fb70_2_ss_4 wh_fb70_3_ss_4 wh_fb70_2_ss_ar  wh_fb70_3_ss_ar using tables/wh_ss.tex, keep(d_immi_std) scalar("fst F" "c5TSLS $\tau=0.05$" "c10TSLS $\tau=0.1$" "c20TSLS $\tau=0.2$") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth) nomtitles ///
	coeflabels(d_immi_std "\midrule \multicolumn{6}{c}{\textit{Panel B}: Change in avg. log weekly wages of high-skilled} \\ \midrule SSIV") nolines star nonumbers nonotes booktabs nogaps replace
label variable d_immi_std "2SLS"
esttab wh_fb70_1_tsls wh_fb70_2_tsls_4 wh_fb70_3_tsls_4 wh_fb70_2_tsls_ar wh_fb70_3_tsls_ar using tables/wh_tsls.tex, keep(d_immi_std) scalar("fst F") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth) nomtitles star nonumbers nolines nonotes posthead(\midrule) booktabs nogaps replace 
label variable d_immi_std "LIML"
esttab wh_fb70_1_liml wh_fb70_2_liml_4 wh_fb70_3_liml_4 wh_fb70_2_liml_ar  wh_fb70_3_liml_ar using tables/wh_liml.tex, keep(d_immi_std) scalar("fst F" "c5LIML $\tau=0.05$" "c10LIML $\tau=0.1$" "c20LIML $\tau=0.2$" "nri \midrule \# inv") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth) nomtitles star nonumbers nolines nonotes posthead(\midrule) booktabs nogaps replace 

	* Wages, LS
esttab wl_fb70_1_ss wl_fb70_2_ss_4 wl_fb70_3_ss_4 wl_fb70_2_ss_ar  wl_fb70_3_ss_ar using tables/wl_ss.tex, keep(d_immi_std) scalar("fst F" "c5TSLS $\tau=0.05$" "c10TSLS $\tau=0.1$" "c20TSLS $\tau=0.2$") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth) nomtitles ///
	star coeflabels(d_immi_std "\midrule \multicolumn{6}{c}{\textit{Panel C}: Change in avg. log weekly wages of low-skilled} \\ \midrule SSIV") nolines nonumbers nonotes booktabs nogaps replace
label variable d_immi_std "2SLS"
esttab wl_fb70_1_tsls wl_fb70_2_tsls_4 wl_fb70_3_tsls_4 wl_fb70_2_tsls_ar  wl_fb70_3_tsls_ar using tables/wl_tsls.tex, keep(d_immi_std) scalar("fst F") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth) nomtitles star nonumbers nolines nonotes posthead(\midrule) booktabs nogaps replace 
label variable d_immi_std "LIML"
esttab wl_fb70_1_liml wl_fb70_2_liml_4 wl_fb70_3_liml_4 wl_fb70_2_liml_ar  wl_fb70_3_liml_ar using tables/wl_liml.tex, keep(d_immi_std) scalar("fst F" "c5LIML $\tau=0.05$" "c10LIML $\tau=0.1$" "c20LIML $\tau=0.2$" "nri \midrule \# inv" "signi Sign.") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth) nomtitles star nonumbers nolines nonotes posthead(\midrule) ///
	postfoot(\bottomrule  ///
	\end{tabular*} }) booktabs nogaps replace
	
* JRSSA Version
* Wages
label variable d_immi_std "TSLS"
esttab w_fb70_1_tsls w_fb70_2_tsls_4 w_fb70_3_tsls_4 w_fb70_2_tsls_ar  w_fb70_3_tsls_ar using tables/w_tsls_short.tex, keep(d_immi_std) scalar("fst F" "c5TSLS $\tau=0.05$" "c10TSLS $\tau=0.1$" "c20TSLS $\tau=0.2$" "nri \midrule \# inv") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth) mtitles("Standard" "AL HS" "CIM HS" "AL AR" "CIM AR") coeflabels(d_immi_std "\toprule \multicolumn{6}{c}{\textit{Panel A}: Change in average log weekly wages}   \\ \midrule TSLS") nonumbers nolines nonotes star booktabs nogaps replace 
	
	* Wages, HS
label variable d_immi_std "TSLS"
esttab wh_fb70_1_tsls wh_fb70_2_tsls_4 wh_fb70_3_tsls_4 wh_fb70_2_tsls_ar  wh_fb70_3_tsls_ar using tables/wh_tsls_short.tex, keep(d_immi_std) scalar("fst F" "c5TSLS $\tau=0.05$" "c10TSLS $\tau=0.1$" "c20TSLS $\tau=0.2$" "nri \midrule \# inv") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth) coeflabels(d_immi_std "\toprule \multicolumn{6}{c}{\textit{Panel B}: Change in avg. log weekly wages of high-skilled}   \\ \midrule TSLS") nonumbers nomtitles nolines nonotes star booktabs nogaps replace 

	* Wages, LS
*est store wl_fb70_2_tsls_4
*estadd local signi "0.1", replace
label variable d_immi_std "TSLS"
esttab wl_fb70_1_tsls wl_fb70_2_tsls_4 wl_fb70_3_tsls_4 wl_fb70_2_tsls_ar wl_fb70_3_tsls_ar using tables/wl_tsls_short.tex, keep(d_immi_std) scalar("fst F" "c5TSLS $\tau=0.05$" "c10TSLS $\tau=0.1$" "c20TSLS $\tau=0.2$" "nri \midrule \# inv" "signi Sign.") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm} C{2cm}) width(\textwidth) coeflabels(d_immi_std "\toprule \multicolumn{6}{c}{\textit{Panel C}: Change in avg. log weekly wages of low-skilled}   \\ \midrule TSLS") nonumbers postfoot(\bottomrule  ///
	\end{tabular*} }) nomtitles nolines nonotes star booktabs nogaps replace 
	
	
* Multiple table

* SSIV
esttab w_mult_orig_ss w_mult_hs_ss w_mult_ar_ss wh_mult_orig_ss wh_mult_hs_ss wh_mult_ar_ss wl_mult_orig_ss wl_mult_hs_ss wl_mult_ar_ss using tables/ss_mult.tex, keep(d_immi_std L.d_immi_std) scalar("fst CD") se not label noobs ///
	alignment(C{1.11cm} C{1.11cm} C{1.38cm} | C{1.09cm}  C{1.27cm} C{1.4cm} | C{1.37cm} C{1.37cm} C{1.25cm} C{1cm}) star(+ 0.10 * 0.05 ** 0.01 *** 0.001) width(\textwidth + \tabcolsep) mtitles("Stand." "AL(HS)" "AL(AR)" "Stand." "AL(HS)" "AL(AR)" "Stand." "AL(HS)" "AL(AR)") ///
	coeflabels(d_immi_std "\toprule \multicolumn{1}{l}{DV:} & \multicolumn{3}{c}{Wages} & \multicolumn{3}{c}{High-skilled} & \multicolumn{3}{c}{Low-skilled} \\ \midrule SSIV $\Delta immi_t$" L.d_immi_std "\hspace{0.7cm} $\Delta immi_{t-10}$") nolines star nonotes booktabs nogaps replace substitute(\_ _)
* TSLS
esttab w_mult_orig_tsls w_mult_hs_tsls w_mult_ar_tsls wh_mult_orig_tsls wh_mult_hs_tsls wh_mult_ar_tsls wl_mult_orig_tsls wl_mult_hs_tsls wl_mult_ar_tsls using tables/tsls_mult.tex, keep(d_immi_std L.d_immi_std) scalar("fst CD") se not label noobs ///
	alignment(C{1.11cm} C{1.11cm} C{1.38cm} | C{1.09cm}  C{1.27cm} C{1.4cm} | C{1.37cm} C{1.37cm} C{1.25cm} C{1cm} ) star(+ 0.10 * 0.05 ** 0.01 *** 0.001) width(\textwidth + \tabcolsep) nomtitles nonumbers ///
	coeflabels(d_immi_std "\midrule TSLS $\Delta immi_t$" L.d_immi_std "\hspace{0.7cm} $\Delta immi_{t-10}$") nolines star nonotes booktabs nogaps replace substitute(\_ _)
* LIML
esttab w_mult_orig_liml w_mult_hs_liml w_mult_ar_liml wh_mult_orig_liml wh_mult_hs_liml wh_mult_ar_liml wl_mult_orig_liml wl_mult_hs_liml wl_mult_ar_liml using tables/liml_mult.tex, keep(d_immi_std L.d_immi_std) scalar("fst CD" "nri \midrule \# inv" "signi Sign.") se not label noobs ///
	alignment(C{1.11cm} C{1.11cm} C{1.38cm} | C{1.09cm}  C{1.27cm} C{1.4cm} | C{1.37cm} C{1.37cm} C{1.25cm} C{1cm} ) star(+ 0.10 * 0.05 ** 0.01 *** 0.001) width(\textwidth + \tabcolsep) nomtitles nonumbers ///
	coeflabels(d_immi_std "\midrule LIML $\Delta immi_t$" L.d_immi_std "\hspace{0.7cm} $\Delta immi_{t-10}$") nolines star nonotes postfoot(\bottomrule  ///
	\end{tabular*} }) booktabs nogaps replace substitute(\_ _)

* TSLS short
esttab w_mult_orig_tsls_short w_mult_hs_tsls_short w_mult_ar_tsls_short wh_mult_orig_tsls_short wh_mult_hs_tsls_short wh_mult_ar_tsls_short wl_mult_orig_tsls_short wl_mult_hs_tsls_short wl_mult_ar_tsls_short using tables/tsls_mult_short.tex, keep(d_immi_std L.d_immi_std) scalar("fst CD" "nri \midrule \# inv" "signi Sign.") se not label noobs ///
	alignment(C{1.11cm} C{1.11cm} C{1.18cm} | C{1.09cm}  C{1.27cm} C{1.16cm} | C{1.37cm} C{1.37cm} C{1.18cm} C{1cm} ) star(+ 0.10 * 0.05 ** 0.01 *** 0.001) width(\textwidth+\tabcolsep) mtitles("Stand." "AL(HS)" "AL(AR)" "Stand." "AL(HS)" "AL(AR)" "Stand." "AL(HS)" "AL(AR)") nonumbers ///
	coeflabels(d_immi_std "\toprule \multicolumn{1}{l}{DV:} & \multicolumn{3}{c}{Wages} & \multicolumn{3}{c}{High-skilled} & \multicolumn{3}{c}{Low-skilled} \\ \midrule TSLS $\Delta immi_t$" L.d_immi_std "\hspace{0.7cm} $\Delta immi_{t-10}$") nolines star nonotes postfoot(\bottomrule  ///
	\end{tabular*} }) booktabs nogaps replace substitute(\_ _)
		
log close
