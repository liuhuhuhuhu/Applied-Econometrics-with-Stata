# Binary Response Models and Sample Selection

## Overview

This project applies binary-response and sample-selection methods to two empirical settings:

1. Modeling arrest outcomes using linear probability and probit models.
2. Estimating women's wage equations while correcting for nonrandom labor-force participation using the Heckman selection model.

The project demonstrates how nonlinear probability models and selection corrections can improve empirical analysis when ordinary least squares is not well suited to the data-generating process.

## Research Questions

This project addresses the following questions:

1. How does the probability of conviction affect the probability of arrest?
2. How do linear probability and probit estimates differ?
3. Does the effect of conviction probability vary across its distribution?
4. How well does the probit model classify arrest outcomes?
5. Does nonrandom labor-force participation create selection bias in wage estimation?
6. How do selected-sample OLS and Heckman-corrected wage estimates differ?

## Empirical Applications

### 1. Arrest Outcomes

The first application models whether an individual was arrested in 1986.

The dependent variable is:

- `arr86`: indicator equal to 1 if the individual was arrested in 1986

Key explanatory variables include:

- `pcnv`: probability of conviction
- `avgsen`: average sentence length
- `tottime`: total prior prison time
- `ptime86`: prison time in 1986
- `inc86`: income in 1986
- `black`: race indicator
- `hispan`: ethnicity indicator
- `born60`: birth-cohort indicator

The sample contains 2,725 observations.

### 2. Women's Wage and Labor-Force Participation

The second application uses the MROZ dataset to estimate a wage equation when wages are observed only for women participating in the labor force.

The wage equation includes:

- `lwage`: log hourly wage
- `exper`: labor-market experience
- `expersq`: squared experience

The labor-force participation equation includes:

- `educ`
- `age`
- `kidslt6`
- `kidsge6`
- `nwifeinc`
- `motheduc`
- `fatheduc`
- `huseduc`

The sample contains 753 women, including 428 labor-force participants with observed wages.

## Econometric Methods

The project implements:

- Linear Probability Model
- Heteroskedasticity-robust standard errors
- Joint significance tests
- Probit estimation
- Predicted probabilities
- Discrete probability changes
- Marginal effects
- Classification analysis
- Quadratic probit specifications
- Turning-point calculations
- Manual Heckman two-step correction
- Stata built-in Heckman two-step estimator
- Heckman maximum-likelihood estimation
- Selected-sample OLS benchmark

## Main Findings

### Arrest Model

The robust linear probability model estimates a coefficient of approximately `-0.154` on the probability of conviction.

Increasing the probability of conviction from 0.25 to 0.75 is associated with:

- A reduction of approximately 7.7 percentage points in the linear probability model
- A reduction of approximately 10.2 percentage points in the probit model

Both models therefore suggest that a higher probability of conviction is associated with a lower probability of arrest.

The variables `avgsen` and `tottime` are not jointly statistically significant.

### Predictive Performance

Using a 0.50 probability cutoff, the probit model achieves approximately 72.7% overall classification accuracy.

However:

- 96.6% of non-arrests are correctly predicted
- Only 10.3% of arrests are correctly predicted

The high overall accuracy is therefore largely driven by the majority non-arrest class. This illustrates why overall accuracy can be misleading when the dependent variable is imbalanced.

### Nonlinear Effects

Quadratic terms for conviction probability, prison time, and income are statistically significant.

The estimated relationship between conviction probability and the latent arrest index is concave, with a turning point of approximately `0.126`.

Above this value, the marginal effect of conviction probability becomes negative.

### Heckman Selection Model

The Heckman correction produces an inverse Mills ratio coefficient of approximately `-0.409`, which is statistically significant.

This provides evidence that labor-force participation is nonrandom and that estimating the wage equation only for working women produces sample-selection bias.

The estimated return to experience is:

- Approximately `0.0476` under selected-sample OLS
- Approximately `0.0424` under the Heckman two-step estimator

The lower Heckman estimate suggests that selected-sample OLS overstates the return to experience.

## Repository Files

- `03_binary_response_and_sample_selection.do`  
  Reproducible Stata workflow covering linear probability models, probit estimation, marginal effects, nonlinear specifications, prediction, classification, and Heckman selection correction.

- `Binary_Response_and_Sample_Selection_Report.pdf`  
  Portfolio-style research report summarizing the theoretical framework, empirical methods, results, interpretation, and limitations.

## How to Run

1. Download or clone this repository.
2. Place the required datasets in a local `data` folder.
3. Confirm the dataset names used near the beginning of the do-file.
4. Open Stata and set the working directory to this project folder.
5. Run:

```stata
do 03_binary_response_and_sample_selection.do
```

The script is currently configured to look for:

```text
data/arrest_data.dta
data/mroz.dta
```

If your local files use different names, update the corresponding file paths near the beginning of the do-file.


## Data Availability

The original datasets are not included in this repository because they were provided for academic use and may be subject to redistribution restrictions.

Users can reproduce the workflow using legally obtained copies of the corresponding datasets and the variable structures documented above.

## Limitations

The linear probability model may generate fitted probabilities outside the unit interval and requires heteroskedasticity-robust inference.

The probit classification results are sensitive to the selected probability threshold. Measures such as precision, recall, specificity, and ROC performance may provide a more complete evaluation than overall accuracy.

The Heckman model relies on distributional assumptions and is more credible when the selection equation contains valid exclusion restrictions that affect participation but do not directly affect wages.

The empirical relationships should therefore be interpreted as conditional associations rather than definitive causal effects.

## Software

- Stata 18

