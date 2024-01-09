cd "C:\Users\na1d22\Dropbox\Surrey\IVSelection\Stata"

use "rawdata\ucdp-nonstate-201-dta\ucdp-nonstate-201.dta", clear

keep location year best_fatality_estimate

collapse (sum) best_fatality_estimate, by(location year)
sort location
drop if year > 2010

gen ygroups = .
replace ygroups = 1990 if year <= 1990
replace ygroups = 2000 if year > 1990 & year <= 2000
replace ygroups = 2010 if year > 2000 & year <= 2010
collapse (sum) best_fatality_estimate, by(location ygroups)
ren ygroups year


gen bplgroup = .

replace bplgroup = 3 if location == "Mexico"
replace bplgroup = 5 if location == "United Kingdom"
replace bplgroup = 6 if location == "China"
replace bplgroup = 9 if location == "Philippines"
replace bplgroup = 11 if location == "India"
replace bplgroup = 13 if location == "Algeria"
replace bplgroup = 13 if location == "Angola"
replace bplgroup = 19 if location == "Afghanistan"
replace bplgroup = 13 if location == "Algeria"
replace bplgroup = 19 if location == "Azerbaijan"
replace bplgroup = 19 if location == "Bangladesh"
replace bplgroup = 2 if location == "Bolivia"
replace bplgroup = 17 if location == "Bosnia-Herzegovina"
replace bplgroup = 2 if location == "Brazil"
replace bplgroup = 13 if location == "Burundi"
replace bplgroup = 13 if location == "Burundi, Tanzania"
replace bplgroup = 19 if location == "Cambodia (Kampuchea)"
replace bplgroup = 13 if location == "Cameroon"
replace bplgroup = 13 if location == "Cameroon, Nigeria"
replace bplgroup = 1 if location == "Canada"
replace bplgroup = 13 if location == "Central African Republic"
replace bplgroup = 13 if location == "Chad"
replace bplgroup = 2 if location == "Colombia"
replace bplgroup = 13 if location == "Comoros"
replace bplgroup = 13 if location == "Congo"
replace bplgroup = 17 if location == "Croatia"
replace bplgroup = 13 if location == "DR Congo (Zaire)"
replace bplgroup = 13 if location == "Djibouti"
replace bplgroup = 13 if location == "Djibouti, Eritrea"
replace bplgroup = 2 if location == "Ecuador, Peru"
replace bplgroup = 2 if location == "Ecuador"
replace bplgroup = 2 if location == "Ecuador, Colombia"
replace bplgroup = 13 if location == "Egypt"
replace bplgroup = 2 if location == "El Salvador"
replace bplgroup = 13 if location == "Eritrea"
replace bplgroup = 13 if location == "Eritrea, Ethiopia"
replace bplgroup = 13 if location == "Ethiopia, Kenya"
replace bplgroup = 13 if location == "Ethiopia"
replace bplgroup = 19 if location == "Georgia"
replace bplgroup = 13 if location == "Ghana"
replace bplgroup = 2 if location == "Guatemala"
replace bplgroup = 2 if location == "Guatemala, Mexico"
replace bplgroup = 13 if location == "Guinea"
replace bplgroup = 13 if location == "Guinea-Bissau"
replace bplgroup = 2 if location == "Haiti"
replace bplgroup = 2 if location == "Honduras"
replace bplgroup = 19 if location == "India, Pakistan"
replace bplgroup = 19 if location == "Indonesia"
replace bplgroup = 19 if location == "Iran"
replace bplgroup = 19 if location == "Iraq"
replace bplgroup = 19 if location == "Iraq, Kuwait"
replace bplgroup = 19 if location == "Iraq, Turkey"
replace bplgroup = 19 if location == "Israel"
replace bplgroup = 13 if location == "Ivory Coast"
replace bplgroup = 2 if location == "Jamaica"
replace bplgroup = 13 if location == "Kenya"
replace bplgroup = 13 if location == "Kenya, Ethiopia"
replace bplgroup = 13 if location == "Kenya, Somalia"
replace bplgroup = 13 if location == "Kenya, Uganda"
replace bplgroup = 19 if location == "Kyrgyzstan"
replace bplgroup = 19 if location == "Laos"
replace bplgroup = 19 if location == "Lebanon"
replace bplgroup = 13 if location == "Lesotho"
replace bplgroup = 13 if location == "Liberia"
replace bplgroup = 13 if location == "Liberia, Ivory Coast"
replace bplgroup = 13 if location == "Liberia, Sierra Leone"
replace bplgroup = 13 if location == "Madagascar (Malagasy)"
replace bplgroup = 13 if location == "Mali"
replace bplgroup = 13 if location == "Mali, Niger"
replace bplgroup = 13 if location == "Mauritania"
replace bplgroup = 1 if location == "Mexico, Canada"
replace bplgroup = 2 if location == "Mexico, Honduras"
replace bplgroup = 17 if location == "Moldova"
replace bplgroup = 13 if location == "Morocco"
replace bplgroup = 13 if location == "Mozambique"
replace bplgroup = 19 if location == "Myanmar (Burma)"
replace bplgroup = 19 if location == "Myanmar (Burma), Thailand"
replace bplgroup = 19 if location == "Nepal"
replace bplgroup = 2 if location == "Nicaragua"
replace bplgroup = 13 if location == "Niger"
replace bplgroup = 13 if location == "Nigeria"
replace bplgroup = 17 if location == "North Macedonia"
replace bplgroup = 19 if location == "Pakistan"
replace bplgroup = 2 if location == "Panama"
replace bplgroup = 2 if location == "Panama, United States of America"
replace bplgroup = 2 if location == "Papua New Guinea"
replace bplgroup = 2 if location == "Paraguay"
replace bplgroup = 2 if location == "Peru"
replace bplgroup = 17 if location == "Romania"
replace bplgroup = 17 if location == "Russia (Soviet Union)"
replace bplgroup = 13 if location == "Rwanda"
replace bplgroup = 13 if location == "Senegal"
replace bplgroup = 13 if location == "Senegal, Mauritania"
replace bplgroup = 17 if location == "Serbia (Yugoslavia)"
replace bplgroup = 13 if location == "Sierra Leone"
replace bplgroup = 13 if location == "Somalia"
replace bplgroup = 13 if location == "Somalia, Djibouti"
replace bplgroup = 13 if location == "South Africa"
replace bplgroup = 12 if location == "Spain"
replace bplgroup = 19 if location == "Sri Lanka"
replace bplgroup = 13 if location == "Sudan"
replace bplgroup = 13 if location == "Sudan, Chad"
replace bplgroup = 19 if location == "Tajikistan"
replace bplgroup = 19 if location == "Thailand"
replace bplgroup = 13 if location == "Togo"
replace bplgroup = 13 if location == "Trinidad and Tobago"
replace bplgroup = 19 if location == "Turkey"
replace bplgroup = 13 if location == "Uganda"
replace bplgroup = 19 if location == "Uzbekistan"
replace bplgroup = 2 if location == "Venezuela"
replace bplgroup = 13 if location == "Yemen (North Yemen)"

collapse (sum) best_fatality_estimate, by(bplgroup year)

ren best_fatality_estimate nonstate
reshape wide nonstate, i(year) j(bplgroup)

forval val = 1/19 {
	cap confirm variable nonstate`val'
	if _rc {
		gen nonstate`val'=0
	}
	else {
		replace nonstate`val'=0 if nonstate`val'==.
	}
}

save data/Nonstate, replace
