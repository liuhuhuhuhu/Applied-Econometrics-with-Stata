********************************************************************************
* Project: Discrete Choice, Matching, Difference-in-Differences, and RDD
* File:    04_discrete_choice_causal_inference.do
* Purpose: Implement conditional logit models, propensity-score matching,
*          difference-in-differences/event-study analysis, and regression
*          discontinuity designs in Stata.
*
* Expected input files:
*   data/energy.dta
*   data/lalonde.dta
*   data/health.dta
*
* Required user-written commands:
*   asclogit
*   psmatch2
*   coefplot
*
* Run this file from the project root directory.
********************************************************************************

version 18.0
clear all
set more off
set linesize 120
capture log close _all

*===============================================================================
* 0. Project configuration
*===============================================================================

global PROJECT_ROOT "`c(pwd)'"
global DATA_DIR    "$PROJECT_ROOT/data"
global OUTPUT_DIR  "$PROJECT_ROOT/output"

capture mkdir "$OUTPUT_DIR"

local energy_data  "$DATA_DIR/energy.dta"
local lalonde_data "$DATA_DIR/lalonde.dta"
local health_data  "$DATA_DIR/health.dta"

log using "$OUTPUT_DIR/discrete_choice_causal_inference.log", ///
    text replace name(mainlog)

* Install required user-written packages only if unavailable.
foreach cmd in asclogit psmatch2 coefplot {
    capture which `cmd'
    if _rc {
        display as text "Installing required package: `cmd'"
        quietly ssc install `cmd', replace
    }
}

*===============================================================================
* 1. Alternative-specific conditional logit model
*===============================================================================

capture confirm file "`energy_data'"
if _rc {
    display as error "Required file not found: `energy_data'"
    log close mainlog
    exit 601
}

use "`energy_data'", clear

local energy_vars depvar ic oc idcase idalt income
foreach var of local energy_vars {
    capture confirm variable `var'
    if _rc {
        display as error "Required variable not found in energy dataset: `var'"
        log close mainlog
        exit 111
    }
}

* Baseline conditional logit model.
asclogit depvar ic oc, ///
    case(idcase) ///
    alternatives(idalt) ///
    basealternative(1)

estimates store clogit_baseline

* Predicted choice probabilities.
predict double phat, pr
label variable phat "Predicted choice probability"

* Store coefficients.
scalar b_ic = _b[ic]
scalar b_oc = _b[oc]

* Own average partial effects:
* effect of increasing the cost of alternative j on Pr(choosing j)
generate double ape_ic_own = b_ic * phat * (1 - phat)
generate double ape_oc_own = b_oc * phat * (1 - phat)

label variable ape_ic_own "Own partial effect of installation cost"
label variable ape_oc_own "Own partial effect of operating cost"

preserve
collapse (mean) phat depvar ape_ic_own ape_oc_own, by(idalt)
sort idalt
list, clean noobs
export delimited using "$OUTPUT_DIR/conditional_logit_average_effects.csv", replace
restore

*===============================================================================
* 2. Conditional logit with income interactions
*===============================================================================

generate double ic_income = ic * income
generate double oc_income = oc * income

label variable ic_income "Installation cost x income"
label variable oc_income "Operating cost x income"

asclogit depvar ic oc ic_income oc_income, ///
    case(idcase) ///
    alternatives(idalt) ///
    basealternative(1)

estimates store clogit_income_interactions

*===============================================================================
* 3. Lalonde data: group construction and descriptive balance
*===============================================================================

capture confirm file "`lalonde_data'"
if _rc {
    display as error "Required file not found: `lalonde_data'"
    log close mainlog
    exit 601
}

use "`lalonde_data'", clear

local lalonde_vars ///
    sample treated age educ black married nodegree hisp kids18 kidmiss re74 re75 re78

foreach var of local lalonde_vars {
    capture confirm variable `var'
    if _rc {
        display as error "Required variable not found in Lalonde dataset: `var'"
        log close mainlog
        exit 111
    }
}

tabulate sample
tabulate treated

* Sample coding:
* sample = 1: NSW experimental sample
* sample = 2: CPS comparison sample

generate byte group = .
replace group = 1 if sample == 1 & treated == 1
replace group = 2 if sample == 1 & treated == 0
replace group = 3 if sample == 2

label define group_lbl ///
    1 "NSW Treatment" ///
    2 "NSW Control" ///
    3 "CPS"
label values group group_lbl
label variable group "Analysis group"

tabulate group

tabstat age educ black married nodegree hisp kids18 kidmiss re74 re75 re78, ///
    by(group) ///
    statistics(mean sd) ///
    columns(statistics)

*===============================================================================
* 4. Propensity-score matching
*===============================================================================

use "`lalonde_data'", clear

* Keep NSW treated observations and CPS comparison observations.
keep if (sample == 1 & treated == 1) | sample == 2

generate byte treatment = 0
replace treatment = 1 if sample == 1 & treated == 1
label variable treatment "NSW treatment versus CPS comparison"

tabulate treatment

* Estimate ATT for 1978 earnings.
psmatch2 treatment ///
    age educ black married nodegree hisp kids18 kidmiss re74 re75, ///
    outcome(re78) ///
    neighbor(1) ///
    common ///
    ate

* Covariate-balance diagnostics.
pstest age educ black married nodegree hisp kids18 kidmiss re74 re75, both graph

graph export "$OUTPUT_DIR/propensity_score_balance.png", ///
    replace width(2000)

*===============================================================================
* 5. Difference-in-differences analysis
*===============================================================================

capture confirm variable year
if _rc {
    display as error "Variable 'year' is required for the DiD section."
    log close mainlog
    exit 111
}

capture confirm variable ada
if _rc {
    display as error "Variable 'ada' is required for the DiD section."
    log close mainlog
    exit 111
}

capture confirm variable disabled
if _rc {
    display as error "Variable 'disabled' is required for the DiD section."
    log close mainlog
    exit 111
}

capture confirm variable wkswork1
if _rc {
    display as error "Variable 'wkswork1' is required for the DiD section."
    log close mainlog
    exit 111
}

capture confirm variable weekly_earn
if _rc {
    display as error "Variable 'weekly_earn' is required for the DiD section."
    log close mainlog
    exit 111
}

* Post-policy period and treatment interaction.
generate byte post = (year >= 1992)
generate byte ada_post = ada * post

label variable post "Post-ADA period"
label variable ada_post "ADA treatment x post period"

* Baseline DiD specifications.
regress wkswork1 ada i.year, vce(robust)
estimates store did_weeks_baseline

regress weekly_earn ada i.year, vce(robust)
estimates store did_wage_baseline

* Extended DiD with demographic controls and state-clustered standard errors.
local controls i.age white black hispanic lths hsgrad somecol

regress wkswork1 ada_post ada i.year `controls' i.region, ///
    vce(cluster statefip)
estimates store did_weeks_controls

regress weekly_earn ada_post ada i.year `controls' i.region, ///
    vce(cluster statefip)
estimates store did_wage_controls

*===============================================================================
* 6. Event-study analysis
*===============================================================================

capture drop ada1988 ada1989 ada1990 ada1992 ada1993 ///
    ada1994 ada1995 ada1996 ada1997

foreach y in 1988 1989 1990 1992 1993 1994 1995 1996 1997 {
    generate byte ada`y' = disabled * (year == `y')
    label variable ada`y' "Disabled x year `y'"
}

* Base year is 1991.
regress wkswork1 ///
    ada1988 ada1989 ada1990 ///
    ada1992 ada1993 ada1994 ada1995 ada1996 ada1997 ///
    disabled i.year, ///
    vce(cluster statefip)

estimates store weeks_event

regress weekly_earn ///
    ada1988 ada1989 ada1990 ///
    ada1992 ada1993 ada1994 ada1995 ada1996 ada1997 ///
    disabled i.year, ///
    vce(cluster statefip)

estimates store wage_event

coefplot weeks_event, ///
    keep(ada1988 ada1989 ada1990 ada1992 ada1993 ada1994 ada1995 ada1996 ada1997) ///
    vertical ///
    yline(0) ///
    xline(4, lpattern(dash)) ///
    title("Event Study: ADA Effect on Weeks Worked") ///
    ytitle("Effect relative to 1991") ///
    xtitle("Year") ///
    coeflabels( ///
        ada1988 = "1988" ///
        ada1989 = "1989" ///
        ada1990 = "1990" ///
        ada1992 = "1992" ///
        ada1993 = "1993" ///
        ada1994 = "1994" ///
        ada1995 = "1995" ///
        ada1996 = "1996" ///
        ada1997 = "1997" ///
    ) ///
    graphregion(color(white)) ///
    name(event_weeks, replace)

graph export "$OUTPUT_DIR/event_study_weeks.png", ///
    replace width(2000)

coefplot wage_event, ///
    keep(ada1988 ada1989 ada1990 ada1992 ada1993 ada1994 ada1995 ada1996 ada1997) ///
    vertical ///
    yline(0) ///
    xline(4, lpattern(dash)) ///
    title("Event Study: ADA Effect on Weekly Earnings") ///
    ytitle("Effect relative to 1991") ///
    xtitle("Year") ///
    coeflabels( ///
        ada1988 = "1988" ///
        ada1989 = "1989" ///
        ada1990 = "1990" ///
        ada1992 = "1992" ///
        ada1993 = "1993" ///
        ada1994 = "1994" ///
        ada1995 = "1995" ///
        ada1996 = "1996" ///
        ada1997 = "1997" ///
    ) ///
    graphregion(color(white)) ///
    name(event_wage, replace)

graph export "$OUTPUT_DIR/event_study_wage.png", ///
    replace width(2000)

*===============================================================================
* 7. Health-insurance coverage and age-65 discontinuity
*===============================================================================

capture confirm file "`health_data'"
if _rc {
    display as error "Required file not found: `health_data'"
    log close mainlog
    exit 601
}

use "`health_data'", clear

local health_vars ///
    age4 medicare medicaid private anyinsurance sawdr inhosp ///
    female hispanic white black bnh emp dropout hs somecoll college ///
    region year

foreach var of local health_vars {
    capture confirm variable `var'
    if _rc {
        display as error "Required variable not found in health dataset: `var'"
        log close mainlog
        exit 111
    }
}

tabstat medicare medicaid private anyinsurance, ///
    by(d65) ///
    statistics(mean sd n)

preserve
collapse (mean) medicare private anyinsurance, by(age4)

twoway ///
    (connected anyinsurance age4, sort) ///
    (connected medicare age4, sort) ///
    (connected private age4, sort), ///
    xline(65, lpattern(dash)) ///
    title("Insurance Coverage by Age") ///
    ytitle("Coverage Rate") ///
    xtitle("Age") ///
    legend(order(1 "Any insurance" 2 "Medicare" 3 "Private")) ///
    graphregion(color(white)) ///
    name(insurance_age, replace)

graph export "$OUTPUT_DIR/insurance_coverage_by_age.png", ///
    replace width(2000)

restore

* Reduced-form outcome regressions around age 65.
regress sawdr ///
    anyinsurance female hispanic white black bnh emp ///
    dropout hs somecoll college age age2 ///
    i.region##i.year, ///
    vce(robust)

estimates store doctor_visit_controls

regress inhosp ///
    anyinsurance female hispanic white black bnh emp ///
    dropout hs somecoll college age age2 ///
    i.region##i.year, ///
    vce(robust)

estimates store hospital_admission_controls

*===============================================================================
* 8. Regression discontinuity at age 65
*===============================================================================

use "`health_data'", clear

* Running variable centered at age 65.
generate double age = age4 - 65
generate double age2 = age^2

* Treatment indicator at the Medicare eligibility threshold.
generate byte d65 = (age4 >= 65)

* Interactions for different slopes on either side of the cutoff.
generate double age65 = age * d65
generate double age2d65 = age2 * d65

label variable age "Age centered at 65"
label variable d65 "Age 65 or older"
label variable age65 "Centered age x age-65 indicator"
label variable age2d65 "Centered age squared x age-65 indicator"

* Quadratic RDD specifications.
regress sawdr d65 age age2 age65 age2d65, vce(robust)
estimates store rdd_doctor_visit

regress inhosp d65 age age2 age65 age2d65, vce(robust)
estimates store rdd_hospital_admission

*===============================================================================
* 9. Save derived health data and finish
*===============================================================================

save "$OUTPUT_DIR/health_with_rdd_variables.dta", replace

display as text "Analysis completed successfully."
display as text "Outputs saved to: $OUTPUT_DIR"

log close mainlog
