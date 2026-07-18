*******************************************************************************
* Demand System Estimation in Stata
* Linear and Exact Almost Ideal Demand System (AIDS)
*
* This script:
*   1. Imports and validates product-level market data
*   2. Constructs expenditures, budget shares, and price indices
*   3. Estimates LA-AIDS models using OLS, IV/2SLS, and GMM
*   4. Estimates an exact AIDS specification iteratively
*   5. Computes representative elasticities
*
* Expected input:
*   data/raw/products.xlsx
*   worksheet: products
*
* Required variables:
*   market, product, servings_sold, price_per_serving, price_instrument
*
* Output:
*   output/demand_system_estimation.log
*
* Notes:
*   - Run this file from the module root:
*       01_Demand_System_Estimation/
*   - The raw dataset is not included in the public repository unless its
*     redistribution is permitted.
*******************************************************************************

version 18.0
clear all
set more off
capture log close _all

* -----------------------------------------------------------------------------
* 0. Project paths
* -----------------------------------------------------------------------------

local root "."
local data_raw "`root'/data/raw"
local output   "`root'/output"

capture mkdir "`root'/output"
log using "`output'/demand_system_estimation.log", text replace

* -----------------------------------------------------------------------------
* 1. Import and validate data
* -----------------------------------------------------------------------------

capture confirm file "`data_raw'/products.xlsx"
if _rc {
    di as error "Input file not found: `data_raw'/products.xlsx"
    di as error "Place products.xlsx in data/raw/ and rerun the script."
    exit 601
}

import excel "`data_raw'/products.xlsx", ///
    sheet("products") firstrow clear

local required_vars ///
    market product servings_sold price_per_serving price_instrument

foreach var of local required_vars {
    capture confirm variable `var'
    if _rc {
        di as error "Required variable is missing: `var'"
        exit 111
    }
}

* Convert quantity to numeric when imported as a string.
capture confirm numeric variable servings_sold
if _rc {
    destring servings_sold, replace ignore(",") force
}

* Remove observations that cannot be used in logarithmic transformations.
drop if missing(market, product, servings_sold, ///
    price_per_serving, price_instrument)
drop if servings_sold < 0
drop if price_per_serving <= 0
drop if price_instrument <= 0

* Parse market identifiers.
gen str3 city    = substr(market, 1, 3)
gen str2 quarter = substr(market, 4, 2)

encode product, gen(product_id)
egen market_id = group(market), label

* Basic data checks.
isid market product
tab city
tab quarter
tab product

* -----------------------------------------------------------------------------
* 2. Construct expenditures and budget shares
* -----------------------------------------------------------------------------

gen double expenditure = price_per_serving * servings_sold
bysort market: egen double total_expenditure = total(expenditure)

assert total_expenditure > 0

gen double budget_share = expenditure / total_expenditure
label variable budget_share "Product expenditure share within market"

bysort market: egen double budget_share_sum = total(budget_share)
assert abs(budget_share_sum - 1) < 1e-8

summarize budget_share, detail

* Log transformations.
gen double ln_price             = ln(price_per_serving)
gen double ln_price_iv          = ln(price_instrument)
gen double ln_total_expenditure = ln(total_expenditure)

* Number of products in each market.
bysort market: gen int products_per_market = _N
quietly summarize products_per_market, meanonly
scalar Kbar = r(mean)

* -----------------------------------------------------------------------------
* 3. Stone price index and LA-AIDS regressors
* -----------------------------------------------------------------------------

gen double weighted_ln_price = budget_share * ln_price
bysort market: egen double ln_stone_index = total(weighted_ln_price)

gen double ln_real_expenditure = ///
    ln_total_expenditure - ln_stone_index

* Restricted price term:
* sum_{k != i} ln(p_k) - (K - 1)ln(p_i)
bysort market: egen double sum_ln_price = total(ln_price)
gen double sum_other_ln_price = sum_ln_price - ln_price
gen double price_term = ///
    sum_other_ln_price - (products_per_market - 1) * ln_price

* Construct the analogous excluded instrument.
bysort market: egen double sum_ln_price_iv = total(ln_price_iv)
gen double sum_other_ln_price_iv = sum_ln_price_iv - ln_price_iv
gen double price_term_iv = ///
    sum_other_ln_price_iv - (products_per_market - 1) * ln_price_iv

* Alternative normalized restricted terms used for the IV diagnostic.
gen double mean_ln_price_others = ///
    sum_other_ln_price / (products_per_market - 1)
gen double restricted_ln_price = ///
    ln_price - mean_ln_price_others

gen double mean_ln_price_iv_others = ///
    sum_other_ln_price_iv / (products_per_market - 1)
gen double restricted_ln_price_iv = ///
    ln_price_iv - mean_ln_price_iv_others

* -----------------------------------------------------------------------------
* 4. LA-AIDS estimation using OLS
* -----------------------------------------------------------------------------

* Unrestricted price specification with product fixed effects.
regress budget_share ///
    ln_price sum_other_ln_price ln_real_expenditure ///
    i.product_id, vce(cluster market_id)
estimates store laaids_unrestricted_ols

* Restricted specification.
regress budget_share ///
    price_term ln_real_expenditure ///
    i.product_id, vce(cluster market_id)
estimates store laaids_restricted_ols

scalar gamma_cross_ols = _b[price_term]
scalar gamma_own_ols   = -(Kbar - 1) * gamma_cross_ols
scalar beta_ols        = _b[ln_real_expenditure]

display as text "Restricted OLS estimates"
display as result "gamma_cross = " gamma_cross_ols
display as result "gamma_own   = " gamma_own_ols
display as result "beta        = " beta_ols

* -----------------------------------------------------------------------------
* 5. Instrument relevance and 2SLS estimation
* -----------------------------------------------------------------------------

* First-stage relevance diagnostic.
regress restricted_ln_price ///
    restricted_ln_price_iv i.product_id, ///
    vce(cluster market_id)

test restricted_ln_price_iv
scalar first_stage_F = r(F)

display as text "First-stage excluded-instrument F statistic"
display as result first_stage_F

* OLS benchmark using the normalized restricted regressor.
regress budget_share ///
    restricted_ln_price ln_real_expenditure ///
    i.product_id, vce(cluster market_id)
estimates store normalized_ols

scalar theta_ols = _b[restricted_ln_price]

* 2SLS model.
ivregress 2sls budget_share ///
    (restricted_ln_price = restricted_ln_price_iv) ///
    ln_real_expenditure i.product_id, ///
    vce(cluster market_id)
estimates store laaids_2sls

scalar theta_2sls = _b[restricted_ln_price]
scalar beta_2sls  = _b[ln_real_expenditure]

* Convert normalized coefficients into implied own- and cross-price terms.
scalar gamma_cross_normalized_ols = -theta_ols / Kbar
scalar gamma_own_normalized_ols   = ///
    (Kbar - 1) * theta_ols / Kbar

scalar gamma_cross_2sls = -theta_2sls / Kbar
scalar gamma_own_2sls   = ///
    (Kbar - 1) * theta_2sls / Kbar

display as text "Normalized OLS implied coefficients"
display as result "gamma_cross = " gamma_cross_normalized_ols
display as result "gamma_own   = " gamma_own_normalized_ols

display as text "2SLS implied coefficients"
display as result "gamma_cross = " gamma_cross_2sls
display as result "gamma_own   = " gamma_own_2sls
display as result "beta        = " beta_2sls

* -----------------------------------------------------------------------------
* 6. GMM estimation
* -----------------------------------------------------------------------------

ivregress gmm budget_share ///
    (price_term = price_term_iv) ///
    ln_real_expenditure i.product_id, ///
    vce(cluster market_id)
estimates store laaids_gmm

scalar gamma_cross_gmm = _b[price_term]
scalar gamma_own_gmm   = -(Kbar - 1) * gamma_cross_gmm
scalar beta_gmm        = _b[ln_real_expenditure]

display as text "GMM estimates"
display as result "gamma_cross = " gamma_cross_gmm
display as result "gamma_own   = " gamma_own_gmm
display as result "beta        = " beta_gmm

* -----------------------------------------------------------------------------
* 7. Exact AIDS estimation by iteration
* -----------------------------------------------------------------------------

* Approximate alpha_i with average product budget shares.
bysort product_id: egen double alpha_i = mean(budget_share)

* Use restricted LA-AIDS estimates as starting values.
quietly regress budget_share ///
    price_term ln_real_expenditure ///
    i.product_id, vce(cluster market_id)

scalar gamma_cross_exact = _b[price_term]
scalar beta_exact        = _b[ln_real_expenditure]
scalar gamma_own_exact   = ///
    -(Kbar - 1) * gamma_cross_exact

local tolerance = 1e-4
local max_iter  = 200
local iteration = 0
local difference = 1

while (`difference' > `tolerance') & (`iteration' < `max_iter') {

    local ++iteration

    scalar gamma_cross_old = gamma_cross_exact
    scalar beta_old        = beta_exact

    bysort market: egen double sum_alpha_lnp = ///
        total(alpha_i * ln_price)
    bysort market: egen double sum_lnp = ///
        total(ln_price)
    bysort market: egen double sum_lnp_sq = ///
        total(ln_price^2)

    gen double ln_exact_price_index = ///
        sum_alpha_lnp + ///
        0.5 * ( ///
            gamma_own_exact * sum_lnp_sq + ///
            gamma_cross_exact * (sum_lnp^2 - sum_lnp_sq) ///
        )

    gen double ln_exact_real_expenditure = ///
        ln_total_expenditure - ln_exact_price_index

    quietly regress budget_share ///
        price_term ln_exact_real_expenditure ///
        i.product_id, vce(cluster market_id)

    scalar gamma_cross_exact = _b[price_term]
    scalar beta_exact        = _b[ln_exact_real_expenditure]
    scalar gamma_own_exact   = ///
        -(Kbar - 1) * gamma_cross_exact

    scalar parameter_difference = max( ///
        abs(gamma_cross_exact - gamma_cross_old), ///
        abs(beta_exact - beta_old) ///
    )

    local difference = parameter_difference

    drop ln_exact_price_index ln_exact_real_expenditure ///
        sum_alpha_lnp sum_lnp sum_lnp_sq
}

display as text "Exact AIDS iteration completed"
display as result "Iterations  = `iteration'"
display as result "Final diff  = `difference'"
display as result "gamma_cross = " gamma_cross_exact
display as result "gamma_own   = " gamma_own_exact
display as result "beta        = " beta_exact

if (`iteration' == `max_iter') & (`difference' > `tolerance') {
    display as error ///
        "Warning: the exact AIDS estimator reached the iteration limit."
}

* -----------------------------------------------------------------------------
* 8. Representative elasticity calculations
* -----------------------------------------------------------------------------

quietly summarize budget_share, meanonly
scalar mean_budget_share = r(mean)

* Elasticities evaluated at the average budget share.
scalar uncompensated_own_price_elasticity = ///
    -1 + gamma_own_gmm / mean_budget_share - beta_gmm

scalar expenditure_elasticity = ///
    1 + beta_gmm / mean_budget_share

scalar compensated_own_price_elasticity = ///
    uncompensated_own_price_elasticity + ///
    expenditure_elasticity * mean_budget_share

display as text "Representative elasticities evaluated at mean budget share"
display as result "Uncompensated own-price elasticity = " ///
    uncompensated_own_price_elasticity
display as result "Expenditure elasticity            = " ///
    expenditure_elasticity
display as result "Compensated own-price elasticity   = " ///
    compensated_own_price_elasticity

* -----------------------------------------------------------------------------
* 9. Close log
* -----------------------------------------------------------------------------

log close
display as text "Demand-system estimation completed successfully."
