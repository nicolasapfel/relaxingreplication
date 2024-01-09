cd "C:\Users\na1d22\Dropbox\Surrey\IVSelection\Stata" // Adjust accordingly.

use "data/analysis_data", clear // Adjust this path.
ren survey year

ren ss_1970 shs1
replace shs1 = shs1/1000

merge m:1 year using data/brd, nogen
drop prod*
* Create BRD IV
forval val = 1/19{
	gen prod70_`val' = share70_`val' * brd`val'
}
egen shs2 = rowtotal(prod70_*)
replace shs2 = shs2/1000

merge m:1 year using data/Nonstate, nogen
drop prod*
* Create Nonstate IV
forval val = 1/19{
	gen prod70_`val' = share70_`val' * nonstate`val'
}
egen shs3 = rowtotal(prod70_*)
replace shs3 = shs3/1000

merge m:1 year using data/Onesided, nogen
drop prod*
* Create Onesided violence IV
forval val = 1/19{
	gen prod70_`val' = share70_`val' * onesided`val'
}
egen shs4 = rowtotal(prod70_*)

merge m:1 year using data/Population, nogen
drop prod*
* Create Population IV
forval val = 1/19{
	gen prod70_`val' = share70_`val' * Pop`val'
}
egen shs5 = rowtotal(prod70_*)
replace shs5 = shs5/1000000000

merge m:1 year using data/FreedomHouse, nogen
drop prod*
* Create Civil liberties IV
forval val = 1/19{
	gen prod70_`val' = share70_`val' * CL`val'
}
egen shs6 = rowtotal(prod70_*)

drop prod*
* Create Political rights IV
forval val = 1/19{
	gen prod70_`val' = share70_`val' * PR`val'
}
egen shs7 = rowtotal(prod70_*)

drop prod*
* Create FH Status IV
forval val = 1/19{
	gen prod70_`val' = share70_`val' * Status`val'
}
egen shs8 = rowtotal(prod70_*)

merge m:1 year using data/PolityV, nogen
drop prod*
* Create Polity IV
forval val = 1/19{
	gen prod70_`val' = share70_`val' * polity`val'
}
egen shs9 = rowtotal(prod70_*)

merge m:1 year using data/FreedomPress, nogen

drop prod*
* Create Press Fredom Status IV
forval val = 1/19{
	gen prod70_`val' = share70_`val' * PStatus`val'
}
egen shs10 = rowtotal(prod70_*)

drop prod*
* Create Press Fredom Score IV
forval val = 1/19{
	gen prod70_`val' = share70_`val' * PScore`val'
}
egen shs11 = rowtotal(prod70_*)

ren year survey
xtset czone decade

save data/analysis_shifts_data, replace
use data/analysis_shifts_data, clear

*****************
* MAIN ANALYSIS *
*****************
quietly{
*est drop _al
foreach dv in "dlweekly" "dlweekly_hskill" "dlweekly_lskill" { // Dependent variable used
loc dvn

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

use data/analysis_shifts_data, replace
* Standard: 1
* All shares are assumed to be valid
**** OLS
xi: reg `dv' d_immi_std i.survey [aw=l.population] if survey>=1990, vce(cluster czone)
est store `dvn'_shs_sta_ols
**** 2SLS
ivregress 2sls `dv' (d_immi_std= shs*) year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store `dvn'_shs_sta_tsls
estat firststage
matrix singleres = r(singleresults)
scalar fst = singleres[1,4]
estadd scalar fst
estadd local nri "0"
**** LIML
ivregress liml `dv' (d_immi_std= shs*) year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store `dvn'_shs_sta_liml
estat firststage
matrix singleres = r(singleresults)
scalar fst = singleres[1,4]
estadd scalar fst
estadd local nri "0"
estadd local signi "-"

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

cap ssada `dv' d_immi_std year* if survey >= 1990 [aw=lpop], endog(d_immi_std) exog(shs* year*) id(czone survey) vce(cluster czone) c(`cval')
loc wv `e(wv)'
loc wi `e(wi)'
di `wv'
di _rc
if _rc == 0 {
**** 2SLS
ivregress 2sls `dv' `wi' (d_immi_std= `wv') year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store `dvn'_shs_ada_tsls_`cabk'
estat firststage
matrix singleres = r(singleresults)
scalar fst = singleres[1,4]
estadd scalar fst
estadd local nri "`ninv'"
local ninv: word count `wi'
**** LIML
ivregress liml `dv' `wi' (d_immi_std= `wv') year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store `dvn'_shs_ada_liml_`cabk'
estat firststage
matrix singleres = r(singleresults)
scalar fst = singleres[1,4]
estadd scalar fst
estadd local nri "`ninv'"
estadd local signi `signi'

loc line "`dv' & ADA & `signi'"
foreach var of varlist shs* {
if (strpos("`wi' ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

} 
else {
est restore `dvn'_shs_sta_tsls
est store `dvn'_shs_ada_tsls_`cabk'
est restore `dvn'_shs_sta_liml
est store `dvn'_shs_ada_liml_`cabk'
estadd local signi `signi', replace

noi di "`dv' & aLasso & `signi' & -& -& -& -& -& -& -& -& -& -& -\\"
}

*SSCIM: 3
* Selection via CIM
di "****************************************"
di "*              SSCIM                   *"
di "****************************************"
cap sscim `dv' d_immi_std year* if survey >= 1990 [aw=lpop], endog(d_immi_std) exog(shs* year*) ssstub(shs) vce("cluster(czone)") c(`cval')
loc wv `e(wv)'
loc wi `e(wi)'
di `wv'
di _rc
if _rc == 0 {
**** 2SLS
ivregress 2sls `dv' `wi' (d_immi_std= `wv') year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store `dvn'_shs_cim_tsls_`cabk'
estat firststage
matrix singleres = r(singleresults)
scalar fst = singleres[1,4]
estadd scalar fst
estadd local nri "`ninv'"
local ninv: word count `wi'
**** LIML
ivregress liml `dv' `wi' (d_immi_std= `wv') year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store `dvn'_shs_cim_liml_`cabk'
estat firststage
matrix singleres = r(singleresults)
scalar fst = singleres[1,4]
estadd scalar fst
estadd local nri "`ninv'"
estadd local signi `signi'

} 
else {
est restore `dvn'_shs_sta_tsls
est store `dvn'_shs_cim_tsls_`cabk'
est restore `dvn'_shs_sta_liml
est store `dvn'_shs_cim_liml_`cabk'
estadd local signi `signi', replace
} // end last else

loc line "`dv' & CIM & `signi'"
foreach var of varlist shs* {
if (strpos("`wi' ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'\\"

} // end foreach cval
} // end foreach dv
}

* AR Analyses
* Following selection analyses for ALasso AR can be found in IVSelection/R/Code/BP-shifts
* Load Data
use data/analysis_shifts_data, clear

*********************************************
* I. SINGLE REGRESSOR Selection via AR-test *
*********************************************
* Single regressor case
* Both aLasso and CIM

* 1. aLasso
quietly{
	* dlweekly
foreach var of varlist shs1 shs2 shs3 shs6 shs7 shs9 shs11 {
	ren `var' k_`var'
}

ivregress 2sls dlweekly year3 year4 k_* (d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store w_shs_ada_tsls_ar
estat firststage
matrix singleres = r(singleresults)
scalar fst = singleres[1,4]
estadd scalar fst
ivregress liml dlweekly year3 year4 k_* (d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store w_shs_ada_liml_ar
estat firststage
matrix singleres = r(singleresults)
scalar fst = singleres[1,4]
estadd scalar fst
estadd local nri "7"
estadd local signi "0.01302"

use data/analysis_shifts_data, clear
loc line "dlweekly & ADA AR & 0.01302"
foreach var of varlist shs* {
if (strpos("shs1  shs2  shs3  shs6  shs7  shs9 shs11 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'\\"

	* dlweekly_hskill
use data/analysis_shifts_data, clear
foreach var of varlist shs1 shs2 shs3 shs4 shs6 shs7 shs8 shs9 shs11 {
	ren `var' k_`var'
}

ivregress 2sls dlweekly_hskill year3 year4 k_* (d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store wh_shs_ada_tsls_ar
estat firststage
matrix singleres = r(singleresults)
scalar fst = singleres[1,4]
estadd scalar fst
ivregress liml dlweekly_hskill year3 year4 k_* (d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store wh_shs_ada_liml_ar
estat firststage
matrix singleres = r(singleresults)
scalar fst = singleres[1,4]
estadd scalar fst
estadd local nri "9"
estadd local signi "0.01302"
estat overid, forceweights forcenonrobust

use data/analysis_shifts_data, clear
loc line "dlweekly\_hskill & ADA AR & 0.01302"
foreach var of varlist shs* {
if (strpos("shs1 shs2 shs3 shs4 shs6 shs7 shs8 shs9 shs11 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'\\"

	* dlweekly_lskill
use data/analysis_shifts_data, clear
foreach var of varlist shs1 shs2 shs3 shs5 shs6 shs7 shs10 shs11 {
	ren `var' k_`var'
}

ivregress 2sls dlweekly_lskill year3 year4 k_* (d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store wl_shs_ada_tsls_ar
estat firststage
matrix singleres = r(singleresults)
scalar fst = singleres[1,4]
estadd scalar fst
ivregress liml dlweekly_lskill year3 year4 k_* (d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store wl_shs_ada_liml_ar
estat firststage
matrix singleres = r(singleresults)
scalar fst = singleres[1,4]
estadd scalar fst
estadd local nri "8"
estadd local signi "0.01302"
estat overid, forceweights forcenonrobust

use data/analysis_shifts_data, clear
loc line "dlweekly\_lskill & ADA AR & 0.01302"
foreach var of varlist shs* {
if (strpos("shs1 shs2 shs3 shs5 shs6 shs7 shs10 shs11 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"
*}

* 2. Confidence Interval Method

foreach dv in "dlweekly" "dlweekly_hskill" "dlweekly_lskill" { // Dependent variable used
loc dvn
*quietly{
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

use data/analysis_shifts_data, clear

* Selection via CIM
di "****************************************"
di "*              SSCIM                   *"
di "****************************************"
sscimliml `dv' d_immi_std year* if survey >= 1990 [aw=lpop], endog(d_immi_std) exog(shs* year*) ssstub(shs) vce("cluster(czone)") c(0.1)
loc wv `e(wv)'
loc wi `e(wi)'
di _rc

**** 2SLS
ivregress 2sls `dv' `wi' (d_immi_std= `wv') year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store `dvn'_shs_cim_tsls_ar
estat firststage
matrix singleres = r(singleresults)
scalar fst = singleres[1,4]
estadd scalar fst
local ninv: word count `wi'

**** LIML
ivreg2 `dv' `wi' (d_immi_std= `wv') year* [aw=l.population] if survey>=1990, robust cluster(czone) liml first
est store `dvn'_shs_cim_liml_ar
matrix singleres = e(first)
scalar fst = singleres[4,1]
estadd scalar fst
estadd local nri "`ninv'"
estadd local signi "0.01302"
di "`dv'"
di "`wi'"

use data/analysis_shifts_data, clear
loc line "`dv' & CIM AR & 0.01302"
foreach var of varlist shs* {
if (strpos("`wi' ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'\\"

} 
}


*Table:
* Wages
esttab w_shs_sta_tsls w_shs_ada_tsls_1 w_shs_cim_tsls_1 w_shs_ada_tsls_4 w_shs_cim_tsls_4 w_shs_ada_tsls_ar w_shs_cim_tsls_ar using tables/w_shs_tsls.tex, keep(d_immi_std) scalar("fst F") se not label noobs ///
	alignment(p{2cm} p{2cm} p{2cm} p{2cm} p{2cm} p{2cm}) width(\textwidth) mtitles("Standard" "AL HS" "CIM HS" "AL HS" "CIM HS" "AL AR" "CIM AR") nonumbers ///
	coeflabels(d_immi_std "\toprule \multicolumn{6}{c}{\textit{Panel A}: Change in average log weekly wages} \\ \midrule 2SLS") nolines star nonotes booktabs nogaps replace
label variable d_immi_std "LIML"
esttab w_shs_sta_liml w_shs_ada_liml_1 w_shs_cim_liml_1 w_shs_ada_liml_4 w_shs_cim_liml_4 w_shs_ada_liml_ar w_shs_cim_liml_ar using tables/w_shs_liml.tex, keep(d_immi_std) scalar("fst F" "nri \midrule \# inv" "signi Sign.") se not label noobs ///
	alignment(p{2cm} p{2cm} p{2cm} p{2cm} p{2cm} p{2cm}) width(\textwidth) nomtitles star nonumbers nolines nonotes posthead(\midrule) booktabs nogaps replace 
	
* Wages, HS
esttab wh_shs_sta_tsls wh_shs_ada_tsls_1 wh_shs_cim_tsls_1 wh_shs_ada_tsls_4 wh_shs_cim_tsls_4 wh_shs_ada_tsls_ar wh_shs_cim_tsls_ar using tables/wh_shs_tsls.tex, keep(d_immi_std) scalar("fst F") se not label noobs ///
	alignment(p{2cm} p{2cm} p{2cm} p{2cm} p{2cm} p{2cm}) width(\textwidth) nomtitles nonumbers ///
	coeflabels(d_immi_std "\midrule \multicolumn{6}{c}{\textit{Panel B}: Change in avg. log weekly wages of high-skilled} \\ \midrule 2SLS") nolines star nonotes booktabs nogaps replace
label variable d_immi_std "LIML"
esttab wh_shs_sta_liml wh_shs_ada_liml_1 wh_shs_cim_liml_1 wh_shs_ada_liml_4 wh_shs_cim_liml_4 wh_shs_ada_liml_ar wh_shs_cim_liml_ar using tables/wh_shs_liml.tex, keep(d_immi_std) scalar("fst F" "nri \midrule \# inv" "signi Sign.") se not label noobs ///
	alignment(p{2cm} p{2cm} p{2cm} p{2cm} p{2cm} p{2cm}) width(\textwidth) nomtitles star nonumbers nolines nonotes posthead(\midrule) booktabs nogaps replace 

* Wages, LS
esttab wl_shs_sta_tsls wl_shs_ada_tsls_1 wl_shs_cim_tsls_1 wl_shs_ada_tsls_4 wl_shs_cim_tsls_4 wl_shs_ada_tsls_ar wl_shs_cim_tsls_ar using tables/wl_shs_tsls.tex, keep(d_immi_std) scalar("fst F") se not label noobs ///
	alignment(p{2cm} p{2cm} p{2cm} p{2cm} p{2cm} p{2cm}) width(\textwidth) nomtitles nonumbers ///
	coeflabels(d_immi_std "\midrule \multicolumn{6}{c}{\textit{Panel C}: Change in avg. log weekly wages of low-skilled} \\ \midrule 2SLS") nolines star nonotes booktabs nogaps replace
label variable d_immi_std "LIML"
esttab wl_shs_sta_liml wl_shs_ada_liml_1 wl_shs_cim_liml_1 wl_shs_ada_liml_4 wl_shs_cim_liml_4 wl_shs_ada_liml_ar wl_shs_cim_liml_ar using tables/wl_shs_liml.tex, keep(d_immi_std) scalar("fst F" "nri \midrule \# inv" "signi Sign.") se not label noobs ///
	alignment(p{2cm} p{2cm} p{2cm} p{2cm} p{2cm} p{2cm}) width(\textwidth) nomtitles star nonumbers nolines nonotes posthead(\midrule) postfoot(\bottomrule  ///
	\end{tabular*}})booktabs nogaps replace 

* WITH LAGGED IMMIGRATION
* Selection can be found in IVSelection/R/Code/BP-shifts.R

* Standard
	* Wages
**** 2SLS
ivregress 2sls dlweekly (d_immi_std L.d_immi_std = shs*) year* [aw=l.population] if survey>=1990, robust cluster(czone)
underid
scalar fst = `r(j_uid)'
est store w_shsmult_sta_tsls
estadd scalar fst
estadd local nri "0"
**** LIML
ivregress liml dlweekly (d_immi_std L.d_immi_std = shs*) year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store w_shsmult_sta_liml
estadd scalar fst
estadd local nri "0"
estadd local signi "-"

	* Wages, HS
**** 2SLS
ivregress 2sls dlweekly_hskill (d_immi_std L.d_immi_std = shs*) year* [aw=l.population] if survey>=1990, robust cluster(czone)
underid
scalar fst = `r(j_uid)'
est store wh_shsmult_sta_tsls
estadd scalar fst
estadd local nri "0"
**** LIML
ivregress liml dlweekly_hskill (d_immi_std L.d_immi_std = shs*) year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store wh_shsmult_sta_liml
estadd scalar fst
estadd local nri "0"
estadd local signi "-"

	* Wages, LS
**** 2SLS
ivregress 2sls dlweekly_lskill (d_immi_std L.d_immi_std = shs*) year* [aw=l.population] if survey>=1990, robust cluster(czone)
underid
scalar fst = `r(j_uid)'
est store wl_shsmult_sta_tsls
estadd scalar fst
estadd local nri "0"
**** LIML
ivregress liml dlweekly_lskill (d_immi_std L.d_immi_std = shs*) year* [aw=l.population] if survey>=1990, robust cluster(czone)
est store wl_shsmult_sta_liml
estadd scalar fst
estadd local nri "0"
estadd local signi "-"

* AdaLasso
	*Wages
use data/analysis_shifts_data, clear
foreach var of varlist shs2 shs7 shs11 {
	ren `var' k_`var'
}

ivregress 2sls dlweekly year3 year4 k_* (d_immi_std L.d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
underid
scalar fst = `r(j_uid)'
est store w_shsmult_ada_tsls_4
estadd scalar fst
ivregress liml dlweekly year3 year4 k_* (d_immi_std L.d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store w_shsmult_ada_liml_4
estadd scalar fst
estadd local nri "3"
estadd local signi "0.1"

use data/analysis_shifts_data, clear
loc line "dlweekly & aLasso HS & 0.1"
foreach var of varlist shs* {
if (strpos("shs2 shs7 shs11 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

	*Wages, HS
use data/analysis_shifts_data, clear
foreach var of varlist shs2 shs7 {
	ren `var' k_`var'
}

ivregress 2sls dlweekly_hskill year3 year4 k_* (d_immi_std L.d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store wh_shsmult_ada_tsls_4
underid
scalar fst = `r(j_uid)'
estadd scalar fst
ivregress liml dlweekly_hskill year3 year4 k_* (d_immi_std L.d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store wh_shsmult_ada_liml_4
estadd scalar fst
estadd local nri "2"
estadd local signi "0.1"

use data/analysis_shifts_data, clear
loc line "dlweekly\\_hskill & aLasso HS & 0.1"
foreach var of varlist shs* {
if (strpos("shs2 shs7 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

	*Wages, LS
use data/analysis_shifts_data, clear
foreach var of varlist shs2 shs7 shs11 {
	ren `var' k_`var'
}

ivregress 2sls dlweekly_lskill year3 year4 k_* (d_immi_std L.d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store wl_shsmult_ada_tsls_4
underid
scalar fst = `r(j_uid)'
estadd scalar fst
ivregress liml dlweekly_lskill year3 year4 k_* (d_immi_std L.d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store wl_shsmult_ada_liml_4
estadd scalar fst
estadd local nri "3"
estadd local signi "0.1"

use data/analysis_shifts_data, clear
loc line "dlweekly\\_lskill & aLasso HS & 0.1"
foreach var of varlist shs* {
if (strpos("shs2 shs7 shs11 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

* AdaLasso AR
	*Wages
use data/analysis_shifts_data, clear
foreach var of varlist shs2 shs3 shs4 shs6 shs7 shs8 shs10 shs11 {
	ren `var' k_`var'
}

ivregress 2sls dlweekly year3 year4 k_* (d_immi_std L.d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
underid
scalar fst = `r(j_uid)'
est store w_shsmult_ada_tsls_ar
estadd scalar fst
ivregress liml dlweekly year3 year4 k_* (d_immi_std L.d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store w_shsmult_ada_liml_ar
estadd scalar fst
estadd local nri "8"
estadd local signi "0.01302"

use data/analysis_shifts_data, clear
loc line "dlweekly & aLasso AR & 0.01302"
foreach var of varlist shs* {
if (strpos("shs2 shs3 shs4 shs6 shs7 shs8 shs10 shs11 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

	*Wages, HS
use data/analysis_shifts_data, clear
foreach var of varlist shs2 shs3 shs4 shs6 shs7 shs8 shs11 {
	ren `var' k_`var'
}

ivregress 2sls dlweekly_hskill year3 year4 k_* (d_immi_std L.d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
underid
scalar fst = `r(j_uid)'
est store wh_shsmult_ada_tsls_ar
estadd scalar fst
ivregress liml dlweekly_hskill year3 year4 k_* (d_immi_std L.d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store wh_shsmult_ada_liml_ar
estadd scalar fst
estadd local nri "7"
estadd local signi "0.01302"

use data/analysis_shifts_data, clear
loc line "dlweekly\\_hskill & aLasso AR & 0.01302"
foreach var of varlist shs* {
if (strpos("shs2 shs3 shs4 shs6 shs7 shs8 shs11 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

	*Wages, LS
use data/analysis_shifts_data, clear
foreach var of varlist shs2 shs3 shs5 shs6 shs7 shs10 shs11 {
	ren `var' k_`var'
}

ivregress 2sls dlweekly_lskill year3 year4 k_* (d_immi_std L.d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
underid
scalar fst = `r(j_uid)'
est store wl_shsmult_ada_tsls_ar
estadd scalar fst
ivregress liml dlweekly_lskill year3 year4 k_* (d_immi_std L.d_immi_std = shs*) [aw=lpop] if survey >=1990, robust cluster(czone)
est store wl_shsmult_ada_liml_ar
estadd scalar fst
estadd local nri "7"
estadd local signi "0.01302"

use data/analysis_shifts_data, clear
loc line "dlweekly\\_lskill & aLasso AR & 0.01302"
foreach var of varlist shs* {
if (strpos("shs2 shs3 shs5 shs6 shs7 shs10 shs11 ", "`var' ") != 0) {
loc line = "`line' & x"
} 
else {
loc line = "`line' & -"
}
}
noi di "`line'"

* Multiple table

* 2SLS
esttab w_shsmult_sta_tsls w_shsmult_ada_tsls_4 w_shsmult_ada_tsls_ar wh_shsmult_sta_tsls wh_shsmult_ada_tsls_4 wh_shsmult_ada_tsls_ar wl_shsmult_sta_tsls wl_shsmult_ada_tsls_4 wl_shsmult_ada_tsls_ar using tables/tsls_shs_mult.tex, keep(d_immi_std L.d_immi_std) scalar("fst J") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm}  | C{2cm} C{2cm} C{2cm} | C{2cm} C{2cm} C{2cm} C{2cm} ) star(+ 0.10 * 0.05 ** 0.01 *** 0.001) width(\textwidth + 4\tabcolsep) mtitles("Standard" "AL (HS)" "AL (AR)" "Standard" "AL (HS)" "AL (AR)" "Standard" "AL (HS)" "AL (AR)") ///
	coeflabels(d_immi_std "\toprule \multicolumn{1}{l}{DV:} & \multicolumn{3}{c}{Wages} & \multicolumn{3}{c}{High-skilled} & \multicolumn{3}{c}{Low-skilled} \\ \midrule 2SLS $\Delta immi_t$" L.d_immi_std "\hspace{0.7cm} $\Delta immi_{t-10}$") nolines star nonotes booktabs nogaps replace substitute(\_ _)
* LIML
esttab w_shsmult_sta_liml w_shsmult_ada_liml_4 w_shsmult_ada_liml_ar wh_shsmult_sta_liml wh_shsmult_ada_liml_4 wh_shsmult_ada_liml_ar wl_shsmult_sta_liml wl_shsmult_ada_liml_4 wl_shsmult_ada_liml_ar using tables/liml_shs_mult.tex, keep(d_immi_std L.d_immi_std) scalar("fst J" "nri \midrule \# inv" "signi Sign.") se not label noobs ///
	alignment(C{2cm} C{2cm} C{2cm}  | C{2cm} C{2cm} C{2cm} | C{2cm} C{2cm} C{2cm} C{2cm} ) star(+ 0.10 * 0.05 ** 0.01 *** 0.001) width(\textwidth + 4\tabcolsep) nomtitles nonumbers ///
	coeflabels(d_immi_std "\midrule LIML $\Delta immi_t$" L.d_immi_std "\hspace{0.7cm} $\Delta immi_{t-10}$") nolines star nonotes postfoot(\bottomrule  ///
	\end{tabular*} }) booktabs nogaps replace substitute(\_ _)
	
	
