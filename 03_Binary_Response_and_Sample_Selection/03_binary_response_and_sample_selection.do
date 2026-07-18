********************************************************************************
* Project: Binary Response Models and Sample Selection
* File:    03_binary_response_and_sample_selection.do
* Purpose: Estimate linear probability and probit models for arrest outcomes,
*          evaluate nonlinear effects and predictive performance, and implement
*          Heckman sample-selection corrections for women's wage equations.
*
* Expected input files:
*   data/arrest_data.dta
*   data/mroz.dta
*
* NOTE:
*   If your local data files use different names, update the two file names in
*   Section 0 below. Run this do-file from the project folder.
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

local arrest_data "$DATA_DIR/arrest_data.dta"
local mroz_data   "$DATA_DIR/mroz.dta"

log using "$OUTPUT_DIR/binary_response_and_sample_selection.log", ///
    text replace name(mainlog)

*===============================================================================
* 1. Linear probability model for arrest outcomes
*===============================================================================

capture confirm file "`arrest_data'"
if _rc {
    display as error "Required file not found: `arrest_data'"
    display as error "Place the arrest dataset in the data/ folder or update local arrest_data."
    log close mainlog
    exit 601
}

use "`arrest_data'", clear

local arrest_vars narr86 pcnv avgsen tottime ptime86 inc86 black hispan born60
foreach var of local arrest_vars {
    capture confirm variable `var'
    if _rc {
        display as error "Required variable not found in arrest dataset: `var'"
        log close mainlog
        exit 111
    }
}

* Binary indicator equal to one if the individual was arrested in 1986.
generate byte arr86 = (narr86 > 0) if !missing(narr86)
label variable arr86 "Any arrest in 1986"

* Baseline linear probability model.
regress arr86 pcnv avgsen tottime ptime86 inc86 black hispan born60
estimates store lpm_nonrobust

* Heteroskedasticity-robust linear probability model.
regress arr86 pcnv avgsen tottime ptime86 inc86 black hispan born60, ///
    vce(robust)
estimates store lpm_robust

* Test whether average sentence length and total prior prison time are jointly
* significant predictors of arrest probability.
test avgsen tottime

*===============================================================================
* 2. Probit model, discrete changes, and prediction
*===============================================================================

probit arr86 pcnv avgsen tottime ptime86 inc86 black hispan born60
estimates store probit_baseline

* Evaluate predicted arrest probabilities at representative covariate values
* while changing prior conviction probability from 0.25 to 0.75.
margins, at( ///
    pcnv  = (0.25 0.75) ///
    avgsen = 0.6322936 ///
    tottime = 0.8387523 ///
    inc86 = 54.96705 ///
    ptime86 = 0.387156 ///
    black = 1 ///
    hispan = 0 ///
    born60 = 1 ///
) post

* Discrete change in predicted probability between pcnv = 0.25 and 0.75.
lincom _b[2._at] - _b[1._at]

* Re-estimate the baseline probit because margins, post replaces e(b).
probit arr86 pcnv avgsen tottime ptime86 inc86 black hispan born60
estimates store probit_for_prediction

* Predicted arrest probabilities.
predict double phat if e(sample), pr
label variable phat "Predicted probability of arrest"

* Classification using a 0.50 cutoff.
generate byte arr86_hat = (phat > 0.50) if !missing(phat)
label variable arr86_hat "Predicted arrest status, cutoff 0.50"

* Classification table with row percentages.
tabulate arr86 arr86_hat, row

* Additional predictive-performance measures.
quietly count if arr86 == arr86_hat & !missing(arr86, arr86_hat)
local correct = r(N)

quietly count if !missing(arr86, arr86_hat)
local classified = r(N)

display as result "Classification accuracy at 0.50 cutoff: " ///
    %9.4f (`correct' / `classified')

*===============================================================================
* 3. Nonlinear probit specification
*===============================================================================

* Squared terms for prior conviction probability, prior prison time, and income.
generate double pcnv2    = pcnv^2
generate double ptime862 = ptime86^2
generate double inc862   = inc86^2

label variable pcnv2    "Squared prior conviction probability"
label variable ptime862 "Squared prior prison time"
label variable inc862   "Squared 1986 income"

probit arr86 ///
    pcnv pcnv2 ///
    avgsen tottime ///
    ptime86 ptime862 ///
    inc86 inc862 ///
    black hispan born60

estimates store probit_quadratic

* Joint significance of the nonlinear terms.
test pcnv2 ptime862 inc862

* Turning point for pcnv in the latent-index equation:
* -beta_pcnv / (2*beta_pcnv2)
nlcom (pcnv_turning_point: -_b[pcnv] / (2 * _b[pcnv2]))

* Average marginal effect of pcnv evaluated over values from 0 to 1.
margins, dydx(pcnv) at(pcnv = (0(0.1)1))

marginsplot, ///
    title("Marginal Effect of Prior Conviction Probability") ///
    xtitle("Prior conviction probability (pcnv)") ///
    ytitle("Marginal effect on arrest probability") ///
    graphregion(color(white)) ///
    name(pcnv_marginal_effect, replace)

graph export "$OUTPUT_DIR/pcnv_marginal_effect.png", ///
    replace width(2000)

*===============================================================================
* 4. Heckman sample-selection model: MROZ application
*===============================================================================

capture confirm file "`mroz_data'"
if _rc {
    display as error "Required file not found: `mroz_data'"
    display as error "Place the MROZ dataset in the data/ folder or update local mroz_data."
    log close mainlog
    exit 601
}

use "`mroz_data'", clear

local mroz_vars ///
    inlf lwage exper expersq educ age kidslt6 kidsge6 ///
    nwifeinc motheduc fatheduc huseduc

foreach var of local mroz_vars {
    capture confirm variable `var'
    if _rc {
        display as error "Required variable not found in MROZ dataset: `var'"
        log close mainlog
        exit 111
    }
}

describe `mroz_vars'
summarize `mroz_vars'

* Create squared labor-market experience if it is not already available.
capture confirm variable expersq
if _rc {
    generate double expersq = exper^2
}
label variable expersq "Squared labor-market experience"

*-------------------------------------------------------------------------------
* 4.1 Manual Heckman two-step procedure
*-------------------------------------------------------------------------------

* Selection equation: labor-force participation.
probit inlf ///
    educ age kidslt6 kidsge6 nwifeinc motheduc fatheduc huseduc

estimates store selection_probit

* Linear index from the first-stage probit.
predict double xb_selection if e(sample), xb
label variable xb_selection "Selection-equation linear prediction"

* Inverse Mills ratio for observations selected into employment.
generate double imr = normalden(xb_selection) / normal(xb_selection) ///
    if inlf == 1 & !missing(xb_selection)
label variable imr "Inverse Mills ratio"

* Wage equation corrected for nonrandom labor-force participation.
regress lwage exper expersq imr if inlf == 1, vce(robust)
estimates store heckman_manual

* Test for sample-selection bias.
test imr

*-------------------------------------------------------------------------------
* 4.2 Stata built-in Heckman two-step estimator
*-------------------------------------------------------------------------------

heckman lwage exper expersq, ///
    select(inlf = ///
        educ age kidslt6 kidsge6 nwifeinc motheduc fatheduc huseduc ///
    ) ///
    twostep

estimates store heckman_twostep

*-------------------------------------------------------------------------------
* 4.3 Full-information maximum-likelihood Heckman estimator
*-------------------------------------------------------------------------------

heckman lwage exper expersq, ///
    select(inlf = ///
        educ age kidslt6 kidsge6 nwifeinc motheduc fatheduc huseduc ///
    )

estimates store heckman_mle

*-------------------------------------------------------------------------------
* 4.4 Naive selected-sample OLS benchmark
*-------------------------------------------------------------------------------

regress lwage exper expersq if inlf == 1, vce(robust)
estimates store ols_selected

*===============================================================================
* 5. Save derived datasets
*===============================================================================

save "$OUTPUT_DIR/mroz_with_selection_terms.dta", replace

display as text "Analysis completed successfully."
display as text "Outputs saved to: $OUTPUT_DIR"

log close mainlog
