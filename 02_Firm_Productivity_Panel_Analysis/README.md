# Firm Productivity, Panel Estimation, and Export Dynamics

## Overview

This project analyzes firm-level productivity and export participation using an unbalanced panel dataset covering 2,199 manufacturing plants from 1981 to 1991.

The analysis estimates a Cobb-Douglas production function using several panel-data methods and evaluates how estimator choice affects the estimated output elasticities of capital and labor. It also constructs a firm-level total factor productivity measure and examines the relationship between productivity and export participation.

## Research Questions

This project addresses the following questions:

1. Are exporting firms larger, more capital intensive, and more productive than non-exporting firms?
2. How do pooled OLS, fixed effects, and random effects estimates differ?
3. Does unobserved firm heterogeneity affect production-function estimates?
4. Are the assumptions underlying the random-effects model supported by the data?
5. Can System GMM address dynamic endogeneity in the production function?
6. Are more productive firms more likely to participate in export markets?

## Data

The dataset contains:

- 10,060 firm-year observations
- 2,199 manufacturing plants
- Annual observations from 1981 to 1991
- An unbalanced panel structure

Key variables include:

- `lnv`: log value added
- `lnk`: log capital
- `lnl`: log labor
- `esales`: export sales
- `plant`: firm identifier
- `year`: time identifier

The original dataset is not included because it was provided for academic use and may be subject to redistribution restrictions.

## Econometric Methods

The project implements:

- Descriptive panel-data analysis
- Kernel density estimation
- Pooled OLS with robust standard errors
- Firm fixed effects
- Two-way fixed effects
- Random effects
- Hausman specification test
- Returns-to-scale hypothesis tests
- Two-step System GMM
- Arellano-Bond serial-correlation tests
- Hansen and Difference-in-Hansen tests
- Total factor productivity construction
- Linear probability model for export participation

## Main Findings

Exporting firms have higher average value added, capital, labor, and capital intensity than non-exporting firms. The raw exporter productivity premium is approximately 1.55 log points.

The fixed-effects estimates are smaller than the pooled OLS estimates, indicating that pooled OLS attributes part of persistent firm productivity differences to capital and labor inputs.

The Hausman test strongly rejects the random-effects assumption that firm-specific effects are uncorrelated with the regressors. The fixed-effects estimator is therefore preferred.

Although System GMM is used to address dynamic endogeneity, the reported specification fails the AR(2), Hansen, and Difference-in-Hansen diagnostic tests. These results highlight the importance of evaluating instrument validity rather than interpreting GMM estimates mechanically.

Estimated firm productivity is positively associated with export participation, while lagged exporter status shows strong persistence in firms' export behavior.

## Repository Files

- `02_firm_productivity_panel_analysis.do`  
  Reproducible Stata workflow covering data validation, descriptive analysis, panel estimation, System GMM, TFP construction, and export dynamics.

- `Firm_Productivity_Panel_Analysis_Report.pdf`  
  Portfolio-style research report summarizing the methodology, results, diagnostics, interpretation, and limitations.

## How to Run

1. Download or clone this repository.
2. Place `firm_panel.dta` in a local `data` folder.
3. Open Stata and set the working directory to this project folder.
4. Run:

```stata
do 02_firm_productivity_panel_analysis.do
```

The script requires the user-written Stata package `xtabond2`. If it is not already installed, the script installs it automatically.

## Software

- Stata 18
- `xtabond2`

## Limitations

The firm fixed-effects estimator removes time-invariant firm heterogeneity but does not fully address simultaneity between current input choices and time-varying productivity shocks.

The System GMM specification fails several diagnostic tests, so its coefficient estimates should not be interpreted as preferred causal estimates. Future work could explore alternative lag structures, smaller instrument sets, difference GMM, or proxy-variable production-function estimators.

