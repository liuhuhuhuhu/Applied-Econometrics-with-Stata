# Discrete Choice and Causal Inference

## Overview

This project applies several econometric methods to consumer choice and policy evaluation problems.

The analysis consists of four empirical components:

1. Household heating-system choice using an alternative-specific conditional logit model.
2. Job-training evaluation using propensity-score matching.
3. Labor-market effects of the Americans with Disabilities Act using difference-in-differences and event-study methods.
4. Health-insurance and health-care utilization changes at age 65 using regression discontinuity designs.

The project emphasizes both estimation and identification. In addition to reporting coefficients, it evaluates the assumptions required for each method and discusses the limitations of interpreting the results causally.

## Research Questions

This project addresses the following questions:

1. How do installation and operating costs affect household heating-system choices?
2. Does household income change sensitivity to energy-system costs?
3. Can propensity-score matching recover the effect of job training when the comparison group differs substantially from the treated group?
4. How did the Americans with Disabilities Act affect weeks worked and weekly earnings?
5. Do the event-study estimates support the parallel-trends assumption?
6. How does Medicare eligibility at age 65 affect insurance coverage and health-care utilization?

## Empirical Applications

### 1. Household Energy-System Choice

The first application estimates an alternative-specific conditional logit model for household heating-system choice.

The utility of household `i` from alternative `j` depends on:

- `ic`: installation cost
- `oc`: annual operating cost
- Alternative-specific constants
- Interactions between cost and household income

The dataset contains 900 households choosing among five heating-system alternatives.

### 2. Job-Training Evaluation

The second application uses the Lalonde dataset to compare:

- NSW treatment participants
- NSW experimental controls
- CPS nonexperimental controls

The analysis evaluates covariate balance and estimates the effect of job training on 1978 earnings using propensity-score matching.

### 3. Americans with Disabilities Act

The third application evaluates the labor-market effects of the ADA using its implementation in 1992 as a policy intervention.

The treatment group consists of workers with disabilities, while workers without disabilities form the comparison group.

Outcomes include:

- Total weeks worked
- Weekly earnings

The analysis uses:

- Difference-in-differences
- Demographic and regional controls
- State-clustered standard errors
- Event-study specifications

### 4. Medicare Eligibility at Age 65

The final application studies the age-65 Medicare eligibility threshold.

The analysis examines discontinuities in:

- Medicare coverage
- Private insurance
- Any insurance coverage
- Doctor visits
- Hospital admissions

Age centered at 65 is used as the running variable in quadratic regression discontinuity specifications.

## Econometric Methods

The project implements:

- Alternative-specific conditional logit
- Predicted choice probabilities
- Own average partial effects
- Cost-income interaction models
- Covariate balance analysis
- Propensity-score matching
- Average treatment effect on the treated
- Difference-in-differences
- Event-study estimation
- Cluster-robust standard errors
- Regression discontinuity design
- Polynomial running-variable specifications
- Graphical policy analysis

## Main Findings

### Household Energy Choice

Installation and operating costs both reduce the probability that a household chooses a heating system.

The estimated coefficients are approximately:

- Installation cost: `-0.00153`
- Operating cost: `-0.00699`

The operating-cost coefficient is about 4.6 times larger in absolute value than the installation-cost coefficient. This suggests that households place greater weight on recurring annual costs than on one-time installation expenses.

Predicted choice shares closely match observed market shares because the model includes alternative-specific constants.

When cost-income interactions are added, operating cost remains negative and statistically significant, while the interaction terms are not significant. The data therefore provide limited evidence that income systematically changes price sensitivity.

### Job Training and Propensity-Score Matching

The NSW treatment and experimental control groups are similar across observed characteristics, reflecting randomized assignment.

The CPS comparison group differs substantially from the NSW participants in age, education, race, marital status, and pretreatment earnings.

Propensity-score matching produces an estimated treatment effect on 1978 earnings of approximately:

```text
ATET = $494
```

However, the estimate is not statistically significant and is imprecise.

This result illustrates the difficulty of using nonexperimental comparison groups when common support is limited and unobserved differences may remain after matching.

### ADA Difference-in-Differences

The baseline estimates indicate negative post-1992 changes for workers with disabilities relative to workers without disabilities.

The estimated effects are approximately:

- `-26.8` weeks worked
- `-$130` in weekly earnings

These estimates suggest that the ADA may have been associated with reduced labor-market outcomes for disabled workers.

However, the magnitudes should be interpreted cautiously because the validity of the difference-in-differences design depends on parallel trends and the absence of other group-specific shocks.

The event-study specifications are used to evaluate pre-policy trends and examine how effects evolve over time.

### Medicare Eligibility at Age 65

Insurance coverage changes sharply at the age-65 eligibility threshold.

Medicare coverage rises substantially, while the composition of private and public insurance changes.

Regression discontinuity models are used to estimate whether eligibility also affects doctor visits and hospital admissions.

The design provides a stronger causal framework than simple age comparisons because identification relies on observations close to the eligibility cutoff.

## Repository Files

- `04_discrete_choice_causal_inference.do`  
  Reproducible Stata workflow covering conditional logit estimation, propensity-score matching, difference-in-differences, event studies, and regression discontinuity analysis.

- `Discrete_Choice_and_Causal_Inference_Report.pdf`  
  Portfolio-style research report summarizing the empirical methods, results, interpretation, identification assumptions, and limitations.

## How to Run

1. Download or clone this repository.
2. Place the required datasets in a local `data` folder.
3. Confirm the dataset names near the beginning of the do-file.
4. Open Stata and set the working directory to this project folder.
5. Run:

```stata
do 04_discrete_choice_causal_inference.do
```

The script is configured to look for:

```text
data/energy.dta
data/lalonde.dta
data/health.dta
```

The ADA analysis may require a separate labor-market dataset containing variables such as:

```text
year
disabled
ada
wkswork1
weekly_earn
statefip
```

Update the corresponding `use` statement if that section uses a separate file in your local project.

## Required Stata Packages

The script uses the following user-written commands:

- `asclogit`
- `psmatch2`
- `coefplot`

If they are not already installed, the do-file attempts to install them automatically.


## Data Availability

The original datasets are not included because they were provided for academic use and may be subject to redistribution restrictions.

The code can be reproduced using legally obtained datasets with the variable structures documented above.

## Identification and Limitations

The conditional logit model relies on the independence of irrelevant alternatives assumption. This assumption may be restrictive when heating alternatives share unobserved characteristics.

Propensity-score matching removes differences in observed covariates but does not address selection on unobserved characteristics.

Difference-in-differences requires parallel trends between treatment and comparison groups in the absence of the policy. Event-study estimates should be examined carefully before assigning a causal interpretation.

Regression discontinuity estimates are local to individuals near age 65 and depend on smooth potential outcomes and the absence of manipulation around the threshold.

The reported estimates should therefore be interpreted in light of the assumptions specific to each design.

## Software

- Stata 18

