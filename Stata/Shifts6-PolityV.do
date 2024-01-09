cd "C:\Users\na1d22\Dropbox\Surrey\IVSelection\Stata"

set varabbrev off

import excel "C:rawdata\p5v2018.xls", sheet("p5v2018") cellrange(E1:v17549) firstrow clear

keep country year polity2
ren country countryname
keep if year == 1990 | year == 2000 | year == 2010

gen bplgroup=.

drop if countryname == "Germany East" | countryname == "Germany West" | countryname == "United States" | countryname == "Bosnia" 

replace bplgroup = 18 if countryname == "Latvia" | countryname == "Lithuania" | countryname == "Estonia"

replace bplgroup = 10 if countryname == "Vietnam"
replace bplgroup = 3 if countryname == "Mexico"
replace bplgroup = 5 if countryname == "United Kingdom"
replace bplgroup = 5 if countryname == "Pakistan, United Kingdom"
replace bplgroup = 6 if countryname == "China"
replace bplgroup = 9 if countryname == "Philippines"
replace bplgroup = 11 if countryname == "India"
replace bplgroup = 13 if countryname == "Algeria"
replace bplgroup = 13 if countryname == "Algeria, Mali, Mauritania"
replace bplgroup = 13 if countryname == "Algeria, Mauritania"
replace bplgroup = 13 if countryname == "Angola"
replace bplgroup = 13 if countryname == "Angola, Congo"
replace bplgroup = 13 if countryname == "Angola, Namibia"
replace bplgroup = 13 if countryname == "Angola, Namibia, Zambia"
replace bplgroup = 19 if countryname == "Afghanistan"
replace bplgroup = 19 if countryname == "Afghanistan, Pakistan"
replace bplgroup = 19 if countryname == "Afghanistan, Russia (Soviet Union)"
replace bplgroup = 13 if countryname == "Algeria"
replace bplgroup = 2 if countryname == "Argentina"
replace bplgroup = 19 if countryname == "Azerbaijan"
replace bplgroup = 19 if countryname == "Bangladesh"
replace bplgroup = 11 if countryname == "Bangladesh, India"
replace bplgroup = 11 if countryname == "Bangladesh, Myanmar (Burma)"
replace bplgroup = 11 if countryname == "Bangladesh, Myanmar (Burma), Thailand"
replace bplgroup = 19 if countryname == "Bangladesh, India, Pakistan"
replace bplgroup = 19 if countryname == "Bhutan, India"
replace bplgroup = 2 if countryname == "Bolivia"
replace bplgroup = 17 if countryname == "Bosnia-Herzegovina"
replace bplgroup = 17 if countryname == "Bosnia-Herzegovina, Croatia"
replace bplgroup = 17 if countryname == "Bosnia-Herzegovina, Serbia (Yugoslavia)"
replace bplgroup = 13 if countryname == "Botswana, South Africa"
replace bplgroup = 2 if countryname == "Brazil"
replace bplgroup = 13 if countryname == "Burundi"
replace bplgroup = 13 if countryname == "Burundi, DR Congo (Zaire)"
replace bplgroup = 13 if countryname == "Burundi, DR Congo (Zaire), Rwanda"
replace bplgroup = 13 if countryname == "Burundi, Rwanda"
replace bplgroup = 13 if countryname == "Burundi, Tanzania"
replace bplgroup = 19 if countryname == "Cambodia (Kampuchea)"
replace bplgroup = 13 if countryname == "Cameroon"
replace bplgroup = 13 if countryname == "Cameroon, Nigeria"
replace bplgroup = 1 if countryname == "Canada"
replace bplgroup = 13 if countryname == "Central African Republic"
replace bplgroup = 13 if countryname == "Central African Republic, Chad"
replace bplgroup = 13 if countryname == "Central African Republic, Chad, Sudan"
replace bplgroup = 13 if countryname == "Central African Republic, Sudan"
replace bplgroup = 13 if countryname == "Central African Republic, DR Congo (Zaire), Sudan"
replace bplgroup = 13 if countryname == "Chad"
replace bplgroup = 13 if countryname == "Chad, Sudan"
replace bplgroup = 2 if countryname == "Colombia"
replace bplgroup = 13 if countryname == "Comoros"
replace bplgroup = 13 if countryname == "Congo"
replace bplgroup = 17 if countryname == "Croatia"
replace bplgroup = 13 if countryname == "DR Congo (Zaire)"
replace bplgroup = 13 if countryname == "DR Congo (Zaire), "
replace bplgroup = 13 if countryname == "DR Congo (Zaire), Rwanda"
replace bplgroup = 13 if countryname == "DR Congo (Zaire), Rwanda, Uganda"
replace bplgroup = 13 if countryname == "DR Congo (Zaire), Sudan"
replace bplgroup = 13 if countryname == "DR Congo (Zaire), Sudan, Uganda"
replace bplgroup = 13 if countryname == "DR Congo (Zaire), Tanzania"
replace bplgroup = 13 if countryname == "DR Congo (Zaire), Uganda"
replace bplgroup = 13 if countryname == "DR Congo (Zaire), Zambia"
replace bplgroup = 13 if countryname == "Djibouti"
replace bplgroup = 13 if countryname == "Djibouti, Eritrea"
replace bplgroup = 2 if countryname == "Ecuador, Peru"
replace bplgroup = 2 if countryname == "Ecuador"
replace bplgroup = 2 if countryname == "Ecuador, Colombia"
replace bplgroup = 13 if countryname == "Egypt"
replace bplgroup = 2 if countryname == "El Salvador"
replace bplgroup = 13 if countryname == "Eritrea"
replace bplgroup = 13 if countryname == "Eritrea, Ethiopia"
replace bplgroup = 13 if countryname == "Ethiopia, Kenya"
replace bplgroup = 13 if countryname == "Ethiopia, Sudan"
replace bplgroup = 13 if countryname == "Ethiopia, Senegal"
replace bplgroup = 13 if countryname == "Ethiopia, Somalia"
replace bplgroup = 13 if countryname == "Ethiopia"
replace bplgroup = 13 if countryname == "Gambia, Senegal"
replace bplgroup = 19 if countryname == "Georgia"
replace bplgroup = 13 if countryname == "Ghana"
replace bplgroup = 2 if countryname == "Guatemala"
replace bplgroup = 2 if countryname == "Guatemala, Mexico"
replace bplgroup = 13 if countryname == "Guinea"
replace bplgroup = 13 if countryname == "Guinea, Liberia"
replace bplgroup = 13 if countryname == "Guinea, Sierra Leone"
replace bplgroup = 13 if countryname == "Guinea-Bissau"
replace bplgroup = 2 if countryname == "Guyana"
replace bplgroup = 2 if countryname == "Haiti"
replace bplgroup = 2 if countryname == "Honduras"
replace bplgroup = 11 if countryname == "India, Nepal"
replace bplgroup = 19 if countryname == "India, Pakistan"
replace bplgroup = 19 if countryname == "Bangladesh, India, Pakistan"
replace bplgroup = 19 if countryname == "Indonesia"
replace bplgroup = 19 if countryname == "Iran"
replace bplgroup = 19 if countryname == "Iraq"
replace bplgroup = 19 if countryname == "Iraq, Jordan"
replace bplgroup = 19 if countryname == "Iraq, Kuwait"
replace bplgroup = 19 if countryname == "Iraq, Lebanon"
replace bplgroup = 19 if countryname == "Iraq, Kuwait, Turkey"
replace bplgroup = 19 if countryname == "Iraq, Turkey"
replace bplgroup = 19 if countryname == "Israel"
replace bplgroup = 19 if countryname == "Israel, Lebanon"
replace bplgroup = 13 if countryname == "Ivory Coast"
replace bplgroup = 13 if countryname == "Ivory Coast, Liberia"
replace bplgroup = 13 if countryname == "Ivory Coast, Liberia, Sierra Leone"
replace bplgroup = 2 if countryname == "Jamaica"
replace bplgroup = 13 if countryname == "Kenya"
replace bplgroup = 13 if countryname == "Kenya, Ethiopia"
replace bplgroup = 13 if countryname == "Kenya, Somalia"
replace bplgroup = 13 if countryname == "Kenya, Uganda"
replace bplgroup = 13 if countryname == "Kingdom of eSwatini (Swaziland), Namibia, South Africa"
replace bplgroup = 19 if countryname == "Kyrgyzstan"
replace bplgroup = 19 if countryname == "Laos"
replace bplgroup = 19 if countryname == "Lebanon"
replace bplgroup = 13 if countryname == "Lesotho"
replace bplgroup = 13 if countryname == "Liberia"
replace bplgroup = 13 if countryname == "Liberia, Ivory Coast"
replace bplgroup = 13 if countryname == "Liberia, Sierra Leone"
replace bplgroup = 13 if countryname == "Madagascar (Malagasy)"
replace bplgroup = 13 if countryname == "Mali"
replace bplgroup = 13 if countryname == "Mali, Niger"
replace bplgroup = 13 if countryname == "Mauritania"
replace bplgroup = 13 if countryname == "Mauritania, Senegal"
replace bplgroup = 1 if countryname == "Mexico, Canada"
replace bplgroup = 2 if countryname == "Mexico, Honduras"
replace bplgroup = 17 if countryname == "Moldova"
replace bplgroup = 13 if countryname == "Morocco"
replace bplgroup = 13 if countryname == "Mozambique"
replace bplgroup = 13 if countryname == "Mozambique, Zambia, Zimbabwe (Rhodesia)"
replace bplgroup = 19 if countryname == "Myanmar (Burma)"
replace bplgroup = 19 if countryname == "Myanmar (Burma), Thailand"
replace bplgroup = 19 if countryname == "Nepal"
replace bplgroup = 2 if countryname == "Nicaragua"
replace bplgroup = 13 if countryname == "Niger"
replace bplgroup = 13 if countryname == "Nigeria"
replace bplgroup = 13 if countryname == "Nigeria, Sierra Leone"
replace bplgroup = 17 if countryname == "North Macedonia"
replace bplgroup = 17 if countryname == "North Macedonia, Serbia (Yugoslavia)"
replace bplgroup = 19 if countryname == "Pakistan"
replace bplgroup = 2 if countryname == "Panama"
replace bplgroup = 2 if countryname == "Panama, United States of America"
replace bplgroup = 2 if countryname == "Papua New Guinea"
replace bplgroup = 2 if countryname == "Paraguay"
replace bplgroup = 2 if countryname == "Peru"
replace bplgroup = 17 if countryname == "Romania"
replace bplgroup = 17 if countryname == "Russia (Soviet Union)"
replace bplgroup = 13 if countryname == "Rwanda"
replace bplgroup = 19 if countryname == "Saudi Arabia"
replace bplgroup = 19 if countryname == "Saudi Arabia, Turkey"
replace bplgroup = 13 if countryname == "Senegal"
replace bplgroup = 13 if countryname == "Senegal, Mauritania"
replace bplgroup = 17 if countryname == "Serbia (Yugoslavia)"
replace bplgroup = 13 if countryname == "Sierra Leone"
replace bplgroup = 13 if countryname == "Somalia"
replace bplgroup = 13 if countryname == "Somalia, Djibouti"
replace bplgroup = 13 if countryname == "Somalia, Uganda"
replace bplgroup = 13 if countryname == "South Africa"
replace bplgroup = 12 if countryname == "Spain"
replace bplgroup = 19 if countryname == "Sri Lanka"
replace bplgroup = 13 if countryname == "Sudan"
replace bplgroup = 13 if countryname == "Sudan, Chad"
replace bplgroup = 13 if countryname == "Sudan, Uganda"
replace bplgroup = 19 if countryname == "Tajikistan"
replace bplgroup = 13 if countryname == "Tanzania"
replace bplgroup = 19 if countryname == "Thailand"
replace bplgroup = 13 if countryname == "Togo"
replace bplgroup = 13 if countryname == "Trinidad and Tobago"
replace bplgroup = 19 if countryname == "Turkey"
replace bplgroup = 13 if countryname == "Uganda"
replace bplgroup = 19 if countryname == "Uzbekistan"
replace bplgroup = 2 if countryname == "Venezuela"
replace bplgroup = 13 if countryname == "Yemen (North Yemen)"
replace bplgroup = 13 if countryname == "Zimbabwe (Rhodesia)"
replace bplgroup = 13 if countryname == "The Gambia"
replace bplgroup = 13 if countryname == "Congo (Brazzaville)"
replace bplgroup = 13 if countryname == "Congo (Kinshasa)"
replace bplgroup = 13 if countryname == "Congo Kinshasa"
replace bplgroup = 13 if countryname == "Congo Brazzaville"
replace bplgroup = 13 if countryname == "Congo-Brazzaville"

replace bplgroup = 8 if countryname == "Korea, Rep."
replace bplgroup = 8 if countryname == "South Korea"
replace bplgroup = 8 if countryname == "Korea South"
replace bplgroup = 7 if countryname == "Japan"
replace bplgroup = 2 if countryname == "Aruba"
replace bplgroup = 2 if countryname == "American Samoa"
replace bplgroup = 2 if countryname == "Antigua and Barbuda"
replace bplgroup = 2 if countryname == "Bahamas, The"
replace bplgroup = 2 if countryname == "Barbados"
replace bplgroup = 2 if countryname == "Suriname"
replace bplgroup = 2 if countryname == "Uruguay"
replace bplgroup = 2 if countryname == "Venezuela, RB"
replace bplgroup = 19 if countryname == "Armenia"
replace bplgroup = 19 if countryname == "Bahrain"
replace bplgroup = 19 if countryname == "Brunei Darussalam"
replace bplgroup = 19 if countryname == "Brunei"
replace bplgroup = 19 if countryname == "Bhutan"
replace bplgroup = 19 if countryname == "Cyprus"
replace bplgroup = 19 if countryname == "Hong Kong SAR, China"
replace bplgroup = 19 if countryname == "Iran, Islamic Rep."
replace bplgroup = 19 if countryname == "Jordan"
replace bplgroup = 19 if countryname == "Kazakhstan"
replace bplgroup = 19 if countryname == "Kyrgyz Republic"
replace bplgroup = 19 if countryname == "Cambodia"
replace bplgroup = 19 if countryname == "Lao PDR"
replace bplgroup = 19 if countryname == "Kuwait"
replace bplgroup = 19 if countryname == "Macao SAR, China"
replace bplgroup = 19 if countryname == "Myanmar"
replace bplgroup = 19 if countryname == "Mongolia"
replace bplgroup = 19 if countryname == "Malaysia"
replace bplgroup = 19 if countryname == "North Korea"
replace bplgroup = 19 if countryname == "Korea, Dem. Peoples Rep."
replace bplgroup = 19 if countryname == "Korea North"
replace bplgroup = 19 if countryname == "West Bank and Gaza"
replace bplgroup = 19 if countryname == "Qatar"
replace bplgroup = 19 if countryname == "South Asia"
replace bplgroup = 19 if countryname == "Singapore"
replace bplgroup = 19 if countryname == "Timor-Leste"
replace bplgroup = 19 if countryname == "Timor Leste"
replace bplgroup = 19 if countryname == "Yemen"
replace bplgroup = 19 if countryname == "Yemen North"
replace bplgroup = 19 if countryname == "Yemen South"
replace bplgroup = 19 if countryname == "Yemen, Rep."
replace bplgroup = 19 if countryname == "United Arab Emirates"
replace bplgroup = 19 if countryname == "UAE"
replace bplgroup = 19 if countryname == "Taiwan"
replace bplgroup = 19 if countryname == "Turkmenistan"
replace bplgroup = 19 if countryname == "Syrian Arab Republic"
replace bplgroup = 19 if countryname == "Syria"
replace bplgroup = 19 if countryname == "Guam"
replace bplgroup = 19 if countryname == "Maldives"
replace bplgroup = 19 if countryname == "Oman"


* Oceania
replace bplgroup = 14 if countryname == "Australia"
replace bplgroup = 14 if countryname == "Fiji"
replace bplgroup = 14 if countryname == "Micronesia, Fed. Sts."
replace bplgroup = 14 if countryname == "Kiribati"
replace bplgroup = 14 if countryname == "Marshall Islands"
replace bplgroup = 14 if countryname == "New Caledonia"
replace bplgroup = 14 if countryname == "Nauru"
replace bplgroup = 14 if countryname == "New Zealand"
replace bplgroup = 14 if countryname == "Palau"
replace bplgroup = 14 if countryname == "Pacific island small states"
replace bplgroup = 14 if countryname == "Solomon Islands"
replace bplgroup = 14 if countryname == "Vanuatu"
replace bplgroup = 14 if countryname == "Tonga"
replace bplgroup = 14 if countryname == "Tuvalu"
replace bplgroup = 14 if countryname == "Samoa"

* Central and Eastern Europe
replace bplgroup = 17 if countryname == "Austria"
replace bplgroup = 17 if countryname == "Bulgaria"
replace bplgroup = 17 if countryname == "Bosnia and Herzegovina"
replace bplgroup = 17 if countryname == "Kosovo"
replace bplgroup = 17 if countryname == "Czech Republic"
replace bplgroup = 17 if countryname == "Poland"
replace bplgroup = 17 if countryname == "Serbia"
replace bplgroup = 17 if countryname == "Hungary"
replace bplgroup = 17 if countryname == "Slovenia"
replace bplgroup = 17 if countryname == "Slovak Republic"
replace bplgroup = 17 if countryname == "Germany"
replace bplgroup = 17 if countryname == "Belarus"
replace bplgroup = 17 if countryname == "Ukraine"
replace bplgroup = 17 if countryname == "Russian Federation"
replace bplgroup = 17 if countryname == "Macedonia"
replace bplgroup = 17 if countryname == "Montenegro"
replace bplgroup = 17 if countryname == "USSR"
replace bplgroup = 17 if countryname == "Yugoslavia"
replace bplgroup = 17 if countryname == "Czechoslovakia"
replace bplgroup = 17 if countryname == "Russia"


* Rest of Americas
replace bplgroup = 2 if countryname == "Belize"
replace bplgroup = 2 if countryname == "Bermuda"
replace bplgroup = 2 if countryname == "Costa Rica"
replace bplgroup = 2 if countryname == "Caribbean small states"
replace bplgroup = 2 if countryname == "Cuba"
replace bplgroup = 2 if countryname == "Curacao"
replace bplgroup = 2 if countryname == "Cayman Islands"
replace bplgroup = 2 if countryname == "Dominica"
replace bplgroup = 2 if countryname == "Dominican Republic"
replace bplgroup = 2 if countryname == "St. Kitts and Nevis"
replace bplgroup = 2 if countryname == "St. Lucia"
replace bplgroup = 2 if countryname == "St. Martin (French part)"
replace bplgroup = 2 if countryname == "Puerto Rico"
replace bplgroup = 2 if countryname == "St. Vincent and the Grenadines"
replace bplgroup = 2 if countryname == "British Virgin Islands"
replace bplgroup = 2 if countryname == "Virgin Islands (U.S.)"
replace bplgroup = 2 if countryname == "Turks and Caicos Islands"
replace bplgroup = 2 if countryname == "Chile"
replace bplgroup = 2 if countryname == "Grenada"

* Africa
replace bplgroup = 13 if countryname == "Benin"
replace bplgroup = 13 if countryname == "Burkina Faso"
replace bplgroup = 13 if countryname == "Botswana"
replace bplgroup = 13 if countryname == "Cote d'Ivoire"
replace bplgroup = 13 if countryname == "Cote D'Ivoire"
replace bplgroup = 13 if countryname == "Congo, Dem. Rep."
replace bplgroup = 13 if countryname == "Congo, Rep."
replace bplgroup = 13 if countryname == "Cabo Verde"
replace bplgroup = 13 if countryname == "Cape Verde"
replace bplgroup = 13 if countryname == "Gambia, The"
replace bplgroup = 13 if countryname == "Equatorial Guinea"
replace bplgroup = 13 if countryname == "Madagascar"
replace bplgroup = 13 if countryname == "Mauritius"
replace bplgroup = 13 if countryname == "Malawi"
replace bplgroup = 13 if countryname == "Namibia"
replace bplgroup = 13 if countryname == "Zambia"
replace bplgroup = 13 if countryname == "Zimbabwe"
replace bplgroup = 13 if countryname == "South Sudan"
replace bplgroup = 13 if countryname == "Eswatini"
replace bplgroup = 13 if countryname == "Tunisia"
replace bplgroup = 13 if countryname == "Seychelles"
replace bplgroup = 13 if countryname == "Sao Tome and Principe"
replace bplgroup = 13 if countryname == "Egypt, Arab Rep."
replace bplgroup = 13 if countryname == "Gabon"
replace bplgroup = 13 if countryname == "Libya"
replace bplgroup = 13 if countryname == "Gambia"
replace bplgroup = 13 if countryname == "Swaziland"

* Other
replace bplgroup = 15 if countryname == "Other small states"
replace bplgroup = 15 if countryname == "Small states"
replace bplgroup = 15 if countryname == "French Polynesia"
replace bplgroup = 15 if countryname == "Bahamas"

* Baltics
replace bplgroup = 18 if countryname == "Latvia"
replace bplgroup = 18 if countryname == "Lithuania"
replace bplgroup = 18 if countryname == "Estonia"

* Mediterranean countries
replace bplgroup = 16 if countryname == "Andorra"
replace bplgroup = 16 if countryname == "Albania"
replace bplgroup = 16 if countryname == "Greece"
replace bplgroup = 16 if countryname == "Italy"
replace bplgroup = 16 if countryname == "Malta"
replace bplgroup = 16 if countryname == "Portugal"
replace bplgroup = 16 if countryname == "San Marino"
replace bplgroup = 16 if countryname == "Gibraltar"
replace bplgroup = 16 if countryname == ""

* United Kingdom and Ireland
replace bplgroup = 5 if countryname == "Channel Islands"
replace bplgroup = 5 if countryname == "Isle of Man"
replace bplgroup = 5 if countryname == "Ireland"

* Northern
replace bplgroup = 4 if countryname == "Denmark"
replace bplgroup = 4 if countryname == "Finland"
replace bplgroup = 4 if countryname == "Iceland"
replace bplgroup = 4 if countryname == "Norway"
replace bplgroup = 4 if countryname == "Sweden"
replace bplgroup = 4 if countryname == "Greenland"

replace bplgroup = 12 if bplgroup == .

replace polity2 = -8 if countryname == "Afghanistan" & year == 2010
local new = _N +1
set obs `new'
replace countryname = "Armenia" in 470
replace year = 1990 in 470
replace bplgroup = 19 in 470
replace polity2 = 6 in 470

local new = _N +1
set obs `new'
replace countryname = "Azerbaijan" in 471
replace year = 1990 in 471
replace bplgroup = 19 in 471
replace polity2 = -3 in 471

local new = _N +1
set obs `new'
replace countryname = "Belarus" in 472
replace year = 1990 in 472
replace bplgroup = 17 in 472
replace polity2 = 6 in 472

local new = _N +1
set obs `new'
replace countryname = "Croatia" in 473
replace year = 1990 in 473
replace bplgroup = 17 in 473
replace polity2 = -3 in 473

local new = _N +1
set obs `new'
replace countryname = "Estonia" in 474
replace year = 1990 in 474
replace bplgroup = 18 in 474
replace polity2 = 6 in 474

drop if countryname == "Timor Leste"

local new = _N +1
set obs `new'
replace countryname = "Georgia" in 474
replace year = 1990 in 474
replace bplgroup = 17 in 474
replace polity2 = 3 in 474

drop if countryname == "Kosovo"
drop if countryname == "Yemen North"
drop if countryname == "Yemen South"
replace polity2 = -9 if countryname == "Kuwait" & year == 1990

local new = _N +1
set obs `new'
replace countryname = "Kyrgyzstan" in 472
replace year = 1990 in 472
replace bplgroup = 19 in 472
replace polity2 = -3 in 472

local new = _N +1
set obs `new'
replace countryname = "Kazakhstan" in 473
replace year = 1990 in 473
replace bplgroup = 19 in 473
replace polity2 = -3 in 473

local new = _N +1
set obs `new'
replace countryname = "Latvia" in 474
replace year = 1990 in 474
replace bplgroup = 18 in 474
replace polity2 = 8 in 474

drop if countryname == "Lebanon"

local new = _N +1
set obs `new'
replace countryname = "Lithuania" in 472
replace year = 1990 in 472
replace bplgroup = 18 in 472
replace polity2 = 10 in 472

local new = _N +1
set obs `new'
replace countryname = "Macedonia" in 473
replace year = 1990 in 473
replace bplgroup = 17 in 473
replace polity2 = 6 in 473

local new = _N +1
set obs `new'
replace countryname = "Moldova" in 474
replace year = 1990 in 474
replace bplgroup = 17 in 474
replace polity2 = 5 in 474

drop if countryname == "Montenegro"

local new = _N +1
set obs `new'
replace countryname = "Russia" in 474
replace year = 1990 in 474
replace bplgroup = 17 in 474
replace polity2 = 0 in 474

collapse (mean) polity2, by(bplgroup year)
ren polity2 polity
replace polity = polity *(-1)

reshape wide polity, i(year) j(bplgroup)

gen polity15 = 0

save data/PolityV, replace
