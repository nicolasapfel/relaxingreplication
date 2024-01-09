cd "C:\Users\na1d22\Dropbox\Surrey\IVSelection\Stata"

* Battle related deaths
use "rawdata\ucdp-brd-conf-201-dta\ucdp-brd-conf-201.dta", clear // Data from Uppsala Conflict Data Program.

set varabbrev off

keep location_inc year bd_best 
collapse (sum) bd_best, by(location_inc year)
sort location_inc
drop if year > 2010
drop if location_inc == "United States of America"
drop if location_inc == "Afghanistan, United Kingdom, United States of America"
drop if location_inc == "Australia, Iraq, United Kingdom, United States of America"

gen ygroups = .
replace ygroups = 1990 if year <= 1990
replace ygroups = 2000 if year > 1990 & year <= 2000
replace ygroups = 2010 if year > 2000 & year <= 2010
collapse (sum) bd_best, by(location_inc ygroups)
ren ygroups year

gen bplgroup = .

replace bplgroup = 3 if location_inc == "Mexico"
replace bplgroup = 5 if location_inc == "United Kingdom"
replace bplgroup = 6 if location_inc == "China"
replace bplgroup = 9 if location_inc == "Philippines"
replace bplgroup = 11 if location_inc == "India"
replace bplgroup = 13 if location_inc == "Algeria"
replace bplgroup = 13 if location_inc == "Angola"
replace bplgroup = 19 if location_inc == "Afghanistan"
replace bplgroup = 13 if location_inc == "Algeria"
replace bplgroup = 19 if location_inc == "Azerbaijan"
replace bplgroup = 19 if location_inc == "Bangladesh"
replace bplgroup = 17 if location_inc == "Bosnia-Herzegovina"
replace bplgroup = 13 if location_inc == "Burundi"
replace bplgroup = 19 if location_inc == "Cambodia (Kampuchea)"
replace bplgroup = 13 if location_inc == "Cameroon, Nigeria"
replace bplgroup = 13 if location_inc == "Central African Republic"
replace bplgroup = 13 if location_inc == "Chad"
replace bplgroup = 2 if location_inc == "Colombia"
replace bplgroup = 13 if location_inc == "Comoros"
replace bplgroup = 13 if location_inc == "Congo"
replace bplgroup = 17 if location_inc == "Croatia"
replace bplgroup = 13 if location_inc == "DR Congo (Zaire)"
replace bplgroup = 13 if location_inc == "Djibouti"
replace bplgroup = 13 if location_inc == "Djibouti, Eritrea"
replace bplgroup = 2 if location_inc == "Ecuador, Peru"
replace bplgroup = 13 if location_inc == "Egypt"
replace bplgroup = 2 if location_inc == "El Salvador"
replace bplgroup = 13 if location_inc == "Eritrea"
replace bplgroup = 13 if location_inc == "Eritrea, Ethiopia"
replace bplgroup = 13 if location_inc == "Ethiopia"
replace bplgroup = 19 if location_inc == "Georgia"
replace bplgroup = 2 if location_inc == "Guatemala"
replace bplgroup = 13 if location_inc == "Guinea"
replace bplgroup = 13 if location_inc == "Guinea-Bissau"
replace bplgroup = 2 if location_inc == "Haiti"
replace bplgroup = 11 if location_inc == "India, Pakistan"
replace bplgroup = 19 if location_inc == "Indonesia"
replace bplgroup = 19 if location_inc == "Iran"
replace bplgroup = 19 if location_inc == "Iraq"
replace bplgroup = 19 if location_inc == "Iraq, Kuwait"
replace bplgroup = 19 if location_inc == "Israel"
replace bplgroup = 13 if location_inc == "Ivory Coast"
replace bplgroup = 19 if location_inc == "Laos"
replace bplgroup = 19 if location_inc == "Lebanon"
replace bplgroup = 13 if location_inc == "Lesotho"
replace bplgroup = 13 if location_inc == "Liberia"
replace bplgroup = 13 if location_inc == "Mali"
replace bplgroup = 13 if location_inc == "Mauritania"
replace bplgroup = 17 if location_inc == "Moldova"
replace bplgroup = 13 if location_inc == "Morocco"
replace bplgroup = 13 if location_inc == "Mozambique"
replace bplgroup = 19 if location_inc == "Myanmar (Burma)"
replace bplgroup = 19 if location_inc == "Nepal"
replace bplgroup = 2 if location_inc == "Nicaragua"
replace bplgroup = 13 if location_inc == "Niger"
replace bplgroup = 13 if location_inc == "Nigeria"
replace bplgroup = 17 if location_inc == "North Macedonia"
replace bplgroup = 19 if location_inc == "Pakistan"
replace bplgroup = 2 if location_inc == "Panama"
replace bplgroup = 2 if location_inc == "Panama, United States of America"
replace bplgroup = 2 if location_inc == "Papua New Guinea"
replace bplgroup = 2 if location_inc == "Paraguay"
replace bplgroup = 2 if location_inc == "Peru"
replace bplgroup = 17 if location_inc == "Romania"
replace bplgroup = 17 if location_inc == "Russia (Soviet Union)"
replace bplgroup = 13 if location_inc == "Rwanda"
replace bplgroup = 13 if location_inc == "Senegal"
replace bplgroup = 17 if location_inc == "Serbia (Yugoslavia)"
replace bplgroup = 13 if location_inc == "Sierra Leone"
replace bplgroup = 13 if location_inc == "Somalia"
replace bplgroup = 12 if location_inc == "Spain"
replace bplgroup = 19 if location_inc == "Sri Lanka"
replace bplgroup = 13 if location_inc == "Sudan"
replace bplgroup = 19 if location_inc == "Tajikistan"
replace bplgroup = 19 if location_inc == "Thailand"
replace bplgroup = 13 if location_inc == "Trinidad and Tobago"
replace bplgroup = 19 if location_inc == "Turkey"
replace bplgroup = 13 if location_inc == "Uganda"
replace bplgroup = 19 if location_inc == "Uzbekistan"
replace bplgroup = 2 if location_inc == "Venezuela"
replace bplgroup = 13 if location_inc == "Yemen (North Yemen)"

collapse (sum) bd_best, by(bplgroup year)

ren bd_best brd
reshape wide brd, i(year) j(bplgroup)

forval val = 1/19 {
	cap confirm variable brd`val'
	if _rc {
		gen brd`val'=0
	}
	else {
		replace brd`val'=0 if brd`val'==.
	}
}

save data/brd, replace
