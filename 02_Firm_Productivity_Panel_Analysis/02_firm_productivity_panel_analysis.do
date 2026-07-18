********************************************************************************
* Project: Firm Productivity, Panel Estimation, and Export Dynamics
* File:    02_firm_productivity_panel_analysis.do
* Purpose: Estimate production functions using pooled OLS, fixed effects,
*          random effects, and System GMM; construct firm-level TFP; and
*          analyze export participation dynamics.
*
* Required input:
*   data/firm_panel.dta
*
* Expected variables:
*   plant   - firm/plant identifier
*   year    - year
*   lnv     - log output or value added
*   lnk     - log capital
*   lnl     - log labor
*   esales  - export sales
*
* User-written command:
*   xtabond2
*
* Reproducibility:
*   Run this file from the project root directory.
********************************************************************************

version 18.0
clear all
set more off
set linesize 120
capture log close _all

*===============================================================================
* 0. Project paths and setup
*===============================================================================

* Run from the repository module root:
* 02_Firm_Productivity_Panel_Analysis/
global PROJECT_ROOT "`c(pwd)'"
global DATA_DIR    "$PROJECT_ROOT/data"
global OUTPUT_DIR  "$PROJECT_ROOT/output"

capture mkdir "$OUTPUT_DIR"

local data_file "$DATA_DIR/firm_panel.dta"

capture confirm file "`data_file'"
if _rc {
    display as error "Required data file not found:"
    display as error "`data_file'"
    display as error "Place firm_panel.dta in the data/ directory and rerun."
    exit 601
}

log using "$OUTPUT_DIR/firm_productivity_panel_analysis.log", ///
    text replace name(mainlog)

* Install required user-written packages only when unavailable.
capture which xtabond2
if _rc {
    display as text "Installing required package: xtabond2"
    quietly ssc install xtabond2, replace
}

use "`data_file'", clear

* Confirm required variables before estimation.
local required_vars plant year lnv lnk lnl esales
foreach var of local required_vars {
    capture confirm variable `var'
    if _rc {
        display as error "Required variable not found: `var'"
        log close mainlog
        exit 111
    }
}

* Declare the panel structure.
isid plant year
xtset plant year

*===============================================================================
* 1. Panel structure and descriptive analysis
*===============================================================================

xtdescribe

* Number of firms, time periods, observed firm-years, and panel coverage.
egen byte tag_plant = tag(plant)
quietly count if tag_plant
local N = r(N)

egen byte tag_year = tag(year)
quietly count if tag_year
local T = r(N)

local observed = _N
local possible = `N' * `T'

display as result "Number of firms (N):                 " %10.0f `N'
display as result "Number of time periods (T):          " %10.0f `T'
display as result "Observed firm-year records:          " %10.0f `observed'
display as result "Possible firm-year records (N x T):  " %10.0f `possible'
display as result "Panel coverage ratio:                " %10.4f ///
    (`observed' / `possible')

drop tag_plant tag_year

* Export participation indicator.
generate byte exporter = (esales > 0) if !missing(esales)
label define exporter_lbl 0 "Non-exporter" 1 "Exporter"
label values exporter exporter_lbl
label variable exporter "Positive export sales"

* Summary statistics by export status.
by exporter, sort: summarize lnv lnk lnl

* Raw exporter productivity premium:
* mean(log output | exporter) - mean(log output | non-exporter).
quietly summarize lnv if exporter == 1
local mean_lnv_exporter = r(mean)

quietly summarize lnv if exporter == 0
local mean_lnv_nonexporter = r(mean)

display as result "Raw exporter productivity premium in lnv: " %9.4f ///
    (`mean_lnv_exporter' - `mean_lnv_nonexporter')

* Capital intensity: log(K/L) = log(K) - log(L).
generate double ln_kl = lnk - lnl
label variable ln_kl "Log capital-labor ratio"
by exporter, sort: summarize ln_kl

* Within- and between-firm variance decomposition.
xtsum lnv lnk lnl

* Productivity distributions by export status.
twoway ///
    (kdensity lnv if exporter == 1, lpattern(solid)) ///
    (kdensity lnv if exporter == 0, lpattern(dash)), ///
    legend(order(1 "Exporters" 2 "Non-exporters") rows(1)) ///
    title("Distribution of Firm Productivity by Export Status") ///
    subtitle("Kernel density estimates of log output") ///
    xtitle("Log output (lnv)") ///
    ytitle("Density") ///
    graphregion(color(white)) ///
    name(productivity_density, replace)

graph export "$OUTPUT_DIR/productivity_distribution_by_export_status.png", ///
    replace width(2000)

*===============================================================================
* 2. Pooled OLS production-function estimates
*===============================================================================

* Baseline Cobb-Douglas specification:
* lnv_it = alpha + beta_K lnk_it + beta_L lnl_it + error_it
regress lnv lnk lnl, vce(robust)
estimates store pooled_ols

* Test constant returns to scale: H0: beta_K + beta_L = 1.
lincom lnk + lnl
test lnk + lnl = 1

* Add export status as an observable firm characteristic.
regress lnv lnk lnl exporter, vce(robust)
estimates store pooled_ols_exporter

lincom lnk + lnl
test lnk + lnl = 1

*===============================================================================
* 3. Fixed-effects estimates
*===============================================================================

* Firm fixed effects with firm-clustered standard errors.
xtreg lnv lnk lnl, fe vce(cluster plant)
estimates store firm_fe

lincom lnk + lnl
test lnk + lnl = 1

* Two-way fixed effects: firm effects plus year effects.
xtreg lnv lnk lnl i.year, fe vce(cluster plant)
estimates store two_way_fe

lincom lnk + lnl
test lnk + lnl = 1

*===============================================================================
* 4. Random-effects estimates and model selection
*===============================================================================

* Random-effects model with firm-clustered standard errors for reporting.
xtreg lnv lnk lnl, re vce(cluster plant)
estimates store random_effects

lincom lnk + lnl
test lnk + lnl = 1

display as result "Estimated firm-effect variance (sigma_u^2): " ///
    %9.4f (e(sigma_u)^2)
display as result "Estimated idiosyncratic variance (sigma_e^2): " ///
    %9.4f (e(sigma_e)^2)

* Classical Hausman comparison requires compatible covariance estimators.
quietly xtreg lnv lnk lnl, fe
estimates store hausman_fe

quietly xtreg lnv lnk lnl, re
estimates store hausman_re

hausman hausman_fe hausman_re, sigmamore

*===============================================================================
* 5. Dynamic production function: two-step System GMM
*===============================================================================

* Dynamic model with collapsed instruments to limit instrument proliferation.
* The specification follows the original analysis:
* - lagged output and current/lagged inputs treated as GMM-style variables
* - year indicators entered as standard instruments
* - two-step Windmeijer-corrected robust inference with small-sample adjustment
xtabond2 ///
    lnv L.lnv lnk lnl L.lnk L.lnl i.year, ///
    gmmstyle(L.lnv lnk lnl L.lnk L.lnl, lag(2 .) collapse eq(diff)) ///
    gmmstyle(L.lnv lnk lnl L.lnk L.lnl, lag(2 .) collapse eq(level)) ///
    ivstyle(i.year, eq(level)) ///
    twostep robust small

estimates store system_gmm

* Check implied dynamic restrictions:
* delta_K = -rho * beta_K
* delta_L = -rho * beta_L
nlcom ///
    (implied_deltaK: -_b[L.lnv] * _b[lnk]) ///
    (implied_deltaL: -_b[L.lnv] * _b[lnl])

display as result "Estimated coefficient on lagged capital (L.lnk): " ///
    %9.4f _b[L.lnk]
display as result "Estimated coefficient on lagged labor (L.lnl):   " ///
    %9.4f _b[L.lnl]

* Short-run returns to scale.
lincom lnk + lnl
lincom lnk + lnl - 1

* Store preferred System GMM input coefficients for TFP construction.
scalar beta_k_gmm = _b[lnk]
scalar beta_l_gmm = _b[lnl]

display as result "System GMM capital coefficient used for TFP: " ///
    %9.6f scalar(beta_k_gmm)
display as result "System GMM labor coefficient used for TFP:   " ///
    %9.6f scalar(beta_l_gmm)

*===============================================================================
* 6. Firm-level TFP and export dynamics
*===============================================================================

* Construct estimated firm-year TFP as the production-function residual:
* tfp_hat_it = lnv_it - beta_K*lnk_it - beta_L*lnl_it
generate double tfp_hat = ///
    lnv - scalar(beta_k_gmm) * lnk - scalar(beta_l_gmm) * lnl
label variable tfp_hat "Estimated total factor productivity"

* Lagged export status captures persistence in export participation.
generate byte L_exporter = L.exporter
label variable L_exporter "Lagged export participation"

* Linear probability model with year effects and firm-clustered standard errors.
regress exporter tfp_hat L_exporter i.year, vce(cluster plant)
estimates store export_dynamics_lpm

*===============================================================================
* 7. Save derived data and close log
*===============================================================================

compress
save "$OUTPUT_DIR/firm_panel_with_derived_variables.dta", replace

display as text "Analysis completed successfully."
display as text "Outputs saved to: $OUTPUT_DIR"

log close mainlog
