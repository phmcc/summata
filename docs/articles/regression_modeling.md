# Regression Modeling

Regression analysis quantifies the relationship between a response
variable and one or more predictors while accounting for potential
confounding. The standard analytical workflow proceeds in two phases:
*univariable* (unadjusted) analysis, in which each predictor is examined
independently; and *multivariable* (adjusted) analysis, in which effects
conditional on other predictors are examined collectively. Comparing
these estimates found in these analyses reveals the extent of
confounding and identifies independent predictors of a given outcome.

The `summata` package provides three principal functions for regression
analysis:

| Function | Purpose |
|:---|:---|
| [`uniscreen()`](https://phmcc.github.io/summata/reference/uniscreen.md) | Univariable screening of multiple predictors |
| [`fit()`](https://phmcc.github.io/summata/reference/fit.md) | Single univariable/multivariable model |
| [`fullfit()`](https://phmcc.github.io/summata/reference/fullfit.md) | Combined univariable and multivariable analysis |

All functions support a wide range of model types including linear
models, generalized linear models (logistic, Poisson, Gaussian, Gamma,
negative binomial), Cox proportional hazards models, and mixed-effects
models, with consistent syntax and formatted output.

As with other `summata` functions, these functions adhere to the
standard calling convention:

``` r
fullfit(data, outcome, predictors, model_type, ...)
```

where `data` is the dataset, `outcome` is the outcome variable,
`predictors` is a vector of predicting variables, and `model_type` is
the type of model to be implemented. This vignette demonstrates the
various capabilities of these functions using the included sample
dataset.

------------------------------------------------------------------------

## Preliminaries

The examples in this vignette use the `clintrial` dataset included with
`summata`:

``` r
library(summata)
library(survival)

data(clintrial)
data(clintrial_labels)
```

The `clintrial` dataset simulates a multisite oncology clinical trial
with realistic outcomes and site-level clustering. Variables suitable
for different regression approaches include:

| Variable | Description | Type | Application |
|:---|:---|:---|:---|
| `pain_score` | Postoperative Pain Score (0-10) | Continuous | Linear regression |
| `readmission_30d` | 30-Day Readmission | Binary | Logistic regression |
| `any_complication` | Any Complication | Binary | Logistic regression |
| `los_days` | Length of Hospital Stay (days) | Continuous (positive) | Linear or Gamma regression |
| `recovery_days` | Days to Functional Recovery | Continuous (positive) | Linear or Gamma regression |
| `ae_count` | Adverse Event Count | Count (overdispersed) | Negative binomial / quasipoisson |
| `fu_count` | Follow-Up Visit Count | Count (equidispersed) | Poisson regression |
| `pfs_months`, `pfs_status` | Progression-Free Survival Time (months) | Time-to-event | Cox regression |
| `os_months`, `os_status` | Overall Survival Time (months) | Time-to-event | Cox regression |
| `site` | Study Site | Clustering variable | Mixed-effects models |

------------------------------------------------------------------------

## Supported Model Types

The following table summarizes all supported model types (`model_type`)
and families/links (`family`). For most estimates, effect measures are
displayed by default; set `exponentiate = FALSE` to display raw
coefficients.

| `model_type` | `family` | Outcome | Effect Measure | Package |
|:---|:---|:---|:---|:---|
| `"lm"` | — | Continuous | *β* coefficient | `stats` |
| `"glm"` | `"binomial"` | Binary | Odds ratio | `stats` |
| `"glm"` | `"quasibinomial"` | Binary (overdispersed) | Odds ratio | `stats` |
| `"glm"` | `"poisson"` | Count | Rate ratio | `stats` |
| `"glm"` | `"quasipoisson"` | Count (overdispersed) | Rate ratio | `stats` |
| `"glm"` | `"gaussian"` | Continuous | *β* coefficient | `stats` |
| `"glm"` | `"Gamma"` | Positive continuous | Ratio (with log link) | `stats` |
| `"glm"` | `"inverse.gaussian"` | Positive continuous | Coefficient | `stats` |
| `"negbin"` | — | Count (overdispersed) | Rate ratio | `MASS` |
| `"coxph"` | — | Time-to-event | Hazard ratio | `survival` |
| `"clogit"` | — | Matched binary | Odds ratio | `survival` |
| `"lmer"` | — | Continuous (clustered) | *β* coefficient | `lme4` |
| `"glmer"` | `"binomial"` | Binary (clustered) | Odds ratio | `lme4` |
| `"glmer"` | `"poisson"` | Count (clustered) | Rate ratio | `lme4` |
| `"coxme"` | — | Time-to-event (clustered) | Hazard ratio | `coxme` |

------------------------------------------------------------------------

## Univariable Screening

The
[`uniscreen()`](https://phmcc.github.io/summata/reference/uniscreen.md)
function fits separate regression models for each predictor, combining
results into a single table. This approach identifies candidate
predictors for multivariable modeling.

### **Example 1:** Binary Outcome Screen (Logistic Regression)

For binary outcomes, use `model_type = "glm"` (logistic regression is
the default family). Here we examine predictors of 30-day hospital
readmission:

``` r
screening_vars <- c("age", "sex", "race", "bmi", "smoking", 
                    "diabetes", "stage", "ecog", "treatment")

example1 <- uniscreen(
  data = clintrial,
  outcome = "readmission_30d",
  predictors = screening_vars,
  model_type = "glm",
  labels = clintrial_labels
)

example1
#> 
#> Univariable Screening Results
#> Outcome: readmission_30d
#> Model Type: Logistic
#> Predictors Screened: 9
#> Significant (p < 0.05): 7
#> 
#>                    Variable   Group      n Events       OR (95% CI) p-value
#>                      <char>  <char> <char> <char>            <char>  <char>
#>  1:             Age (years)       -    850    427  1.03 (1.02-1.04) < 0.001
#>  2:                     Sex  Female    450    210         reference       -
#>  3:                            Male    400    217  1.36 (1.03-1.78)   0.027
#>  4:                    Race   White    598    298         reference       -
#>  5:                           Black    126     66  1.11 (0.75-1.63)   0.603
#>  6:                           Asian     93     51  1.22 (0.79-1.90)   0.370
#>  7:                           Other     33     12  0.58 (0.28-1.19)   0.136
#>  8: Body Mass Index (kg/m²)       -    838    418  1.00 (0.98-1.03)   0.787
#>  9:          Smoking Status   Never    337    156         reference       -
#> 10:                          Former    311    141  0.96 (0.71-1.31)   0.808
#> 11:                         Current    185    119  2.09 (1.45-3.03) < 0.001
#> 12:                Diabetes      No    637    301         reference       -
#> 13:                             Yes    197    116  1.60 (1.16-2.21)   0.004
#> 14:           Disease Stage       I    211     89         reference       -
#> 15:                              II    263    126  1.26 (0.88-1.82)   0.213
#> 16:                             III    241    133  1.69 (1.16-2.45)   0.006
#> 17:                              IV    132     77  1.92 (1.23-2.98)   0.004
#> 18: ECOG Performance Status       0    265    108         reference       -
#> 19:                               1    302    145  1.34 (0.96-1.87)   0.083
#> 20:                               2    238    139  2.04 (1.43-2.91) < 0.001
#> 21:                               3     37     30 6.23 (2.64-14.70) < 0.001
#> 22:         Treatment Group Control    196     91         reference       -
#> 23:                          Drug A    292    127  0.89 (0.62-1.28)   0.523
#> 24:                          Drug B    362    209  1.58 (1.11-2.24)   0.011
#>                    Variable   Group      n Events       OR (95% CI) p-value
#>                      <char>  <char> <char> <char>            <char>  <char>
```

Each entry represents a separate univariable model. For categorical
predictors, each non-reference level appears as a separate row.

### **Example 2:** Filtering Screen by *p*-value

The `p_threshold` parameter retains only predictors meeting a
significance criterion:

``` r
example2 <- uniscreen(
  data = clintrial,
  outcome = "readmission_30d",
  predictors = screening_vars,
  model_type = "glm",
  p_threshold = 0.01,
  labels = clintrial_labels
)

example2
#> 
#> Univariable Screening Results
#> Outcome: readmission_30d
#> Model Type: Logistic
#> Predictors Screened: 9
#> Significant (p < 0.01): 5
#> 
#>                    Variable   Group      n Events       OR (95% CI) p-value
#>                      <char>  <char> <char> <char>            <char>  <char>
#>  1:             Age (years)       -    850    427  1.03 (1.02-1.04) < 0.001
#>  2:                     Sex  Female    450    210         reference       -
#>  3:                            Male    400    217  1.36 (1.03-1.78)   0.027
#>  4:                    Race   White    598    298         reference       -
#>  5:                           Black    126     66  1.11 (0.75-1.63)   0.603
#>  6:                           Asian     93     51  1.22 (0.79-1.90)   0.370
#>  7:                           Other     33     12  0.58 (0.28-1.19)   0.136
#>  8: Body Mass Index (kg/m²)       -    838    418  1.00 (0.98-1.03)   0.787
#>  9:          Smoking Status   Never    337    156         reference       -
#> 10:                          Former    311    141  0.96 (0.71-1.31)   0.808
#> 11:                         Current    185    119  2.09 (1.45-3.03) < 0.001
#> 12:                Diabetes      No    637    301         reference       -
#> 13:                             Yes    197    116  1.60 (1.16-2.21)   0.004
#> 14:           Disease Stage       I    211     89         reference       -
#> 15:                              II    263    126  1.26 (0.88-1.82)   0.213
#> 16:                             III    241    133  1.69 (1.16-2.45)   0.006
#> 17:                              IV    132     77  1.92 (1.23-2.98)   0.004
#> 18: ECOG Performance Status       0    265    108         reference       -
#> 19:                               1    302    145  1.34 (0.96-1.87)   0.083
#> 20:                               2    238    139  2.04 (1.43-2.91) < 0.001
#> 21:                               3     37     30 6.23 (2.64-14.70) < 0.001
#> 22:         Treatment Group Control    196     91         reference       -
#> 23:                          Drug A    292    127  0.89 (0.62-1.28)   0.523
#> 24:                          Drug B    362    209  1.58 (1.11-2.24)   0.011
#>                    Variable   Group      n Events       OR (95% CI) p-value
#>                      <char>  <char> <char> <char>            <char>  <char>
```

### **Example 3:** Continuous Outcome Screen (Linear Regression)

For continuous outcomes, specify `model_type = "lm"`:

``` r
example4 <- uniscreen(
  data = clintrial,
  outcome = "los_days",
  predictors = c("age", "sex", "stage", "diabetes", "ecog"),
  model_type = "lm",
  labels = clintrial_labels
)

example4
#> 
#> Univariable Screening Results
#> Outcome: los_days
#> Model Type: Linear
#> Predictors Screened: 5
#> Significant (p < 0.05): 5
#> 
#>                    Variable  Group      n Coefficient (95% CI) p-value
#>                      <char> <char> <char>               <char>  <char>
#>  1:             Age (years)      -    830     0.14 (0.11-0.16) < 0.001
#>  2:                     Sex Female    450            reference       -
#>  3:                           Male    400     0.86 (0.20-1.52)   0.010
#>  4:           Disease Stage      I    211            reference       -
#>  5:                             II    263     1.16 (0.30-2.02)   0.009
#>  6:                            III    241     2.72 (1.84-3.60) < 0.001
#>  7:                             IV    132     2.96 (1.92-4.01) < 0.001
#>  8:                Diabetes     No    637            reference       -
#>  9:                            Yes    197     2.59 (1.83-3.35) < 0.001
#> 10: ECOG Performance Status      0    265            reference       -
#> 11:                              1    302     2.12 (1.35-2.88) < 0.001
#> 12:                              2    238     3.78 (2.96-4.59) < 0.001
#> 13:                              3     37     4.90 (3.32-6.48) < 0.001
```

### **Example 4:** Time-to-Event Outcome Screen (Cox Regression)

For survival outcomes, specify the outcome using
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) notation and set
`model_type = "coxph"`:

``` r
example3 <- uniscreen(
  data = clintrial,
  outcome = "Surv(os_months, os_status)",
  predictors = c("age", "sex", "treatment", "stage", "ecog"),
  model_type = "coxph",
  labels = clintrial_labels
)

example3
#> 
#> Univariable Screening Results
#> Outcome: Surv(os_months, os_status)
#> Model Type: Cox PH
#> Predictors Screened: 5
#> Significant (p < 0.05): 5
#> 
#>                    Variable   Group      n Events      HR (95% CI) p-value
#>                      <char>  <char> <char> <char>           <char>  <char>
#>  1:             Age (years)       -    850    609 1.03 (1.03-1.04) < 0.001
#>  2:                     Sex  Female    450    298        reference       -
#>  3:                            Male    400    311 1.30 (1.11-1.53)   0.001
#>  4:         Treatment Group Control    196    151        reference       -
#>  5:                          Drug A    292    184 0.64 (0.52-0.80) < 0.001
#>  6:                          Drug B    362    274 0.94 (0.77-1.15)   0.567
#>  7:           Disease Stage       I    211    127        reference       -
#>  8:                              II    263    172 1.12 (0.89-1.41)   0.337
#>  9:                             III    241    186 1.69 (1.35-2.11) < 0.001
#> 10:                              IV    132    121 3.18 (2.47-4.09) < 0.001
#> 11: ECOG Performance Status       0    265    159        reference       -
#> 12:                               1    302    212 1.36 (1.11-1.67)   0.003
#> 13:                               2    238    194 1.86 (1.51-2.29) < 0.001
#> 14:                               3     37     37 3.06 (2.13-4.38) < 0.001
```

### **Example 5:** Retaining Model Objects

Setting `keep_models = TRUE` stores the fitted model objects for
diagnostics:

``` r
example5 <- uniscreen(
  data = clintrial,
  outcome = "readmission_30d",
  predictors = c("age", "sex", "stage"),
  model_type = "glm",
  keep_models = TRUE
)

# Access individual models
models <- attr(example5, "models")
names(models)
#> [1] "age"   "sex"   "stage"

# Examine a specific model
summary(models[["age"]])
#> 
#> Call:
#> stats::glm(formula = formula, family = family, data = data, model = keep_models, 
#>     x = FALSE, y = TRUE)
#> 
#> Coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -1.780821   0.369774  -4.816 1.46e-06 ***
#> age          0.029855   0.006056   4.930 8.22e-07 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 1178.3  on 849  degrees of freedom
#> Residual deviance: 1152.9  on 848  degrees of freedom
#> AIC: 1156.9
#> 
#> Number of Fisher Scoring iterations: 4
```

------------------------------------------------------------------------

## Multivariable Modeling

The [`fit()`](https://phmcc.github.io/summata/reference/fit.md) function
estimates a single regression model with multiple predictors, producing
adjusted effect estimates. It is effectively a wrapper for existing
model functions in R and other supported packages that enforces
`summata` syntax for more consistent function calling.

### **Example 6:** Logistic Regression

For binary outcomes with multiple predictors:

``` r
example6 <- fit(
  data = clintrial,
  outcome = "readmission_30d",
  predictors = c("age", "sex", "treatment", "stage", "diabetes"),
  model_type = "glm",
  labels = clintrial_labels
)

example6
#> 
#> Multivariable Logistic Model
#> Formula: readmission_30d ~ age + sex + treatment + stage + diabetes
#> n = 834, Events = 417
#> 
#>            Variable   Group      n Events     aOR (95% CI) p-value
#>              <char>  <char> <char> <char>           <char>  <char>
#>  1:     Age (years)       -    834    417 1.03 (1.02-1.05) < 0.001
#>  2:             Sex  Female    444    207        reference       -
#>  3:                    Male    390    210 1.38 (1.04-1.83)   0.028
#>  4: Treatment Group Control    191     88        reference       -
#>  5:                  Drug A    288    126 0.85 (0.58-1.25)   0.409
#>  6:                  Drug B    355    203 1.40 (0.97-2.02)   0.075
#>  7:   Disease Stage       I    207     86        reference       -
#>  8:                      II    261    124 1.25 (0.86-1.84)   0.242
#>  9:                     III    237    131 1.77 (1.19-2.62)   0.005
#> 10:                      IV    129     76 2.05 (1.29-3.24)   0.002
#> 11:        Diabetes      No    637    301        reference       -
#> 12:                     Yes    197    116 1.67 (1.19-2.34)   0.003
```

### **Example 8:** Linear Regression

For continuous outcomes:

``` r
example8 <- fit(
  data = clintrial,
  outcome = "los_days",
  predictors = c("age", "sex", "stage", "ecog"),
  model_type = "lm",
  labels = clintrial_labels
)

example8
#> 
#> Multivariable Linear Model
#> Formula: los_days ~ age + sex + stage + ecog
#> n = 822
#> 
#>                    Variable  Group      n Adj. Coefficient (95% CI) p-value
#>                      <char> <char> <char>                    <char>  <char>
#>  1:             Age (years)      -    822          0.14 (0.11-0.16) < 0.001
#>  2:                     Sex Female    440                 reference       -
#>  3:                           Male    382          0.90 (0.33-1.46)   0.002
#>  4:           Disease Stage      I    206                 reference       -
#>  5:                             II    257          1.21 (0.46-1.97)   0.002
#>  6:                            III    235          2.92 (2.15-3.70) < 0.001
#>  7:                             IV    124          3.18 (2.26-4.11) < 0.001
#>  8: ECOG Performance Status      0    261                 reference       -
#>  9:                              1    295          2.15 (1.46-2.83) < 0.001
#> 10:                              2    229          3.63 (2.89-4.36) < 0.001
#> 11:                              3     37          4.93 (3.50-6.35) < 0.001
```

### **Example 7:** Cox Regression

For time-to-event outcomes:

``` r
example7 <- fit(
  data = clintrial,
  outcome = "Surv(os_months, os_status)",
  predictors = c("age", "sex", "treatment", "stage"),
  model_type = "coxph",
  labels = clintrial_labels
)

example7
#> 
#> Multivariable Cox PH Model
#> Formula: Surv(os_months, os_status) ~ age + sex + treatment + stage
#> n = 847, Events = 606
#> 
#>            Variable   Group      n Events     aHR (95% CI) p-value
#>              <char>  <char> <char> <char>           <char>  <char>
#>  1:     Age (years)       -    847    606 1.04 (1.03-1.04) < 0.001
#>  2:             Sex  Female    450    298        reference       -
#>  3:                    Male    400    311 1.33 (1.13-1.56) < 0.001
#>  4: Treatment Group Control    196    151        reference       -
#>  5:                  Drug A    292    184 0.56 (0.45-0.70) < 0.001
#>  6:                  Drug B    362    274 0.83 (0.67-1.01)   0.062
#>  7:   Disease Stage       I    211    127        reference       -
#>  8:                      II    263    172 1.16 (0.92-1.46)   0.214
#>  9:                     III    241    186 1.90 (1.51-2.38) < 0.001
#> 10:                      IV    132    121 3.57 (2.77-4.60) < 0.001
```

### **Example 9:** Poisson Regression

For count outcomes where variance approximately equals the mean, use
`model_type = "glm"` with `family = "poisson"`. The `fu_count` variable
represents the number of follow-up clinic visits, which is equidispersed
(suitable for standard Poisson):

``` r
example9 <- fit(
  data = clintrial,
  outcome = "fu_count",
  predictors = c("age", "stage", "treatment", "surgery"),
  model_type = "glm",
  family = "poisson",
  labels = clintrial_labels
)

example9
#> 
#> Multivariable Poisson Model
#> Formula: fu_count ~ age + stage + treatment + surgery
#> n = 839, Events = 5517
#> 
#>               Variable   Group      n Events     aRR (95% CI) p-value
#>                 <char>  <char> <char> <char>           <char>  <char>
#>  1:        Age (years)       -    839  5,517 1.00 (1.00-1.00)   0.038
#>  2:      Disease Stage       I    207  1,238        reference       -
#>  3:                         II    261  1,638 1.05 (0.98-1.13)   0.169
#>  4:                        III    240  1,638 1.14 (1.06-1.23) < 0.001
#>  5:                         IV    131  1,003 1.33 (1.21-1.45) < 0.001
#>  6:    Treatment Group Control    191  1,169        reference       -
#>  7:                     Drug A    290  1,910 1.06 (0.99-1.14)   0.101
#>  8:                     Drug B    358  2,438 1.12 (1.04-1.20)   0.002
#>  9: Surgical Resection      No    476  3,091        reference       -
#> 10:                        Yes    363  2,426 1.09 (1.03-1.16)   0.003
```

Output displays rate ratios (RR). A rate ratio of 1.10 indicates a 10%
higher event rate compared to the reference group. For overdispersed
counts (variance \> mean), consider negative binomial or quasipoisson
regression instead (see Example 18).

### **Example 10:** Toggling Reference Rows for Factors

By default, reference rows are shown in output tables
(`reference_rows = TRUE`). Setting this parameter to `FALSE` removes
these reference rows:

``` r
example10 <- fit(
  data = clintrial,
  outcome = "readmission_30d",
  predictors = c("sex", "stage", "treatment"),
  model_type = "glm",
  reference_rows = FALSE,
  labels = clintrial_labels
)

example10
#> 
#> Multivariable Logistic Model
#> Formula: readmission_30d ~ sex + stage + treatment
#> n = 398, Events = 216
#> 
#>           Variable  Group      n Events     aOR (95% CI) p-value
#>             <char> <char> <char> <char>           <char>  <char>
#> 1:             Sex   Male    398    216 1.38 (1.05-1.82)   0.022
#> 2:   Disease Stage     II    263    126 1.23 (0.85-1.78)   0.270
#> 3:                    III    241    133 1.59 (1.09-2.33)   0.016
#> 4:                     IV    132     77 1.91 (1.22-2.99)   0.004
#> 5: Treatment Group Drug A    292    127 0.88 (0.61-1.28)   0.510
#> 6:                 Drug B    361    208 1.50 (1.05-2.14)   0.027
```

### **Example 11:** Confidence Level

The `conf_level` parameter adjusts the confidence interval width:

``` r
example11 <- fit(
  data = clintrial,
  outcome = "readmission_30d",
  predictors = c("age", "sex", "stage"),
  model_type = "glm",
  conf_level = 0.90
)

example11
#> 
#> Multivariable Logistic Model
#> Formula: readmission_30d ~ age + sex + stage
#> n = 847, Events = 425
#> 
#>    Variable  Group      n Events     aOR (90% CI) p-value
#>      <char> <char> <char> <char>           <char>  <char>
#> 1:      age      -    847    425 1.03 (1.02-1.04) < 0.001
#> 2:      sex Female    449    209        reference       -
#> 3:            Male    398    216 1.41 (1.12-1.79)   0.014
#> 4:    stage      I    211     89        reference       -
#> 5:              II    263    126 1.28 (0.94-1.75)   0.196
#> 6:             III    241    133 1.83 (1.32-2.51)   0.002
#> 7:              IV    132     77 2.03 (1.39-2.96)   0.002
```

### **Example 12:** Raw Coefficients

For logistic and Cox models, set `exponentiate = FALSE` to display
log-scale coefficients (*β* rather than *e^(β)*):

``` r
example12 <- fit(
  data = clintrial,
  outcome = "readmission_30d",
  predictors = c("age", "sex", "stage"),
  model_type = "glm",
  exponentiate = FALSE
)

example12
#> 
#> Multivariable Logistic Model
#> Formula: readmission_30d ~ age + sex + stage
#> n = 847, Events = 425
#> 
#>    Variable  Group      n Events Adj. Coefficient (95% CI) p-value
#>      <char> <char> <char> <char>                    <char>  <char>
#> 1:      age      -    847    425       0.03 (0.02 to 0.04) < 0.001
#> 2:      sex Female    449    209                 reference       -
#> 3:            Male    398    216       0.35 (0.07 to 0.62)   0.014
#> 4:    stage      I    211     89                 reference       -
#> 5:              II    263    126      0.25 (-0.13 to 0.62)   0.196
#> 6:             III    241    133       0.60 (0.22 to 0.98)   0.002
#> 7:              IV    132     77       0.71 (0.26 to 1.16)   0.002
```

------------------------------------------------------------------------

## Combined Analysis

The [`fullfit()`](https://phmcc.github.io/summata/reference/fullfit.md)
function integrates univariable screening and multivariable modeling
into a single workflow, producing a combined table with both unadjusted
and adjusted estimates.

This function extends
[`fit()`](https://phmcc.github.io/summata/reference/fit.md) with an
additional parameter (`method`) for variable selection:

``` r
fullfit(data, outcome, predictors, model_type, method, ...)
```

For the `method` parameter, the following options are available. Setting
`method` to one of the following options controls the predictors
entering the multivariable model:

| Method | Description |
|:---|:---|
| `"screen"` | Only predictors meeting the *p*-threshold in univariable analysis are used in the multivariable model (default) |
| `"all"` | All predictors are used in both univariable and multivariable analyses |
| `"custom"` | Multivariable predictors are explicitly specified |

### **Example 13:** Screening-Based Selection

The default `method = "screen"` approach specifies that only predictors
with univariable *p*-value below `p_threshold` are subsequently used in
the multivariable analysis:

``` r
example13 <- fullfit(
  data = clintrial,
  outcome = "readmission_30d",
  predictors = c("age", "sex", "bmi", "smoking", "diabetes",
                 "stage", "treatment"),
  model_type = "glm",
  method = "screen",
  p_threshold = 0.05,
  labels = clintrial_labels
)

example13
#> 
#> Fullfit Analysis Results
#> Outcome: readmission_30d
#> Model Type: glm
#> Method: screen (p < 0.05)
#> Predictors Screened: 7
#> Multivariable Predictors: 6
#> 
#>                    Variable   Group      n Events      OR (95% CI)   Uni p     aOR (95% CI) Multi p
#>                      <char>  <char> <char> <char>           <char>  <char>           <char>  <char>
#>  1:             Age (years)       -    850    427 1.03 (1.02-1.04) < 0.001 1.04 (1.02-1.05) < 0.001
#>  2:                     Sex  Female    450    210        reference       -        reference       -
#>  3:                            Male    400    217 1.36 (1.03-1.78)   0.027 1.39 (1.04-1.86)   0.024
#>  4: Body Mass Index (kg/m²)       -    838    418 1.00 (0.98-1.03)   0.787                -       -
#>  5:          Smoking Status   Never    337    156        reference       -        reference       -
#>  6:                          Former    311    141 0.96 (0.71-1.31)   0.808 1.00 (0.72-1.39)   0.988
#>  7:                         Current    185    119 2.09 (1.45-3.03) < 0.001 2.44 (1.65-3.59) < 0.001
#>  8:                Diabetes      No    637    301        reference       -        reference       -
#>  9:                             Yes    197    116 1.60 (1.16-2.21)   0.004 1.80 (1.28-2.53) < 0.001
#> 10:           Disease Stage       I    211     89        reference       -        reference       -
#> 11:                              II    263    126 1.26 (0.88-1.82)   0.213 1.25 (0.85-1.85)   0.251
#> 12:                             III    241    133 1.69 (1.16-2.45)   0.006 1.76 (1.18-2.63)   0.005
#> 13:                              IV    132     77 1.92 (1.23-2.98)   0.004 2.03 (1.27-3.25)   0.003
#> 14:         Treatment Group Control    196     91        reference       -        reference       -
#> 15:                          Drug A    292    127 0.89 (0.62-1.28)   0.523 0.84 (0.57-1.24)   0.380
#> 16:                          Drug B    362    209 1.58 (1.11-2.24)   0.011 1.39 (0.96-2.03)   0.083
```

### **Example 14:** All Predictors

The `method = "all"` approach includes the same predictors in both
univariable and multivariable analysis:

``` r
example14 <- fullfit(
  data = clintrial,
  outcome = "any_complication",
  predictors = c("age", "sex", "treatment", "stage"),
  model_type = "glm",
  method = "all",
  labels = clintrial_labels
)

example14
#> 
#> Fullfit Analysis Results
#> Outcome: any_complication
#> Model Type: glm
#> Method: all
#> Predictors Screened: 4
#> Multivariable Predictors: 4
#> 
#>            Variable   Group      n Events      OR (95% CI)  Uni p     aOR (95% CI) Multi p
#>              <char>  <char> <char> <char>           <char> <char>           <char>  <char>
#>  1:     Age (years)       -    850    480 1.01 (1.00-1.02)  0.220 1.01 (1.00-1.02)   0.244
#>  2:             Sex  Female    450    241        reference      -        reference       -
#>  3:                    Male    400    239 1.29 (0.98-1.69)  0.069 1.28 (0.97-1.69)   0.076
#>  4: Treatment Group Control    196    111        reference      -        reference       -
#>  5:                  Drug A    292    143 0.73 (0.51-1.06)  0.097 0.74 (0.51-1.07)   0.113
#>  6:                  Drug B    362    226 1.27 (0.89-1.81)  0.182 1.24 (0.87-1.78)   0.236
#>  7:   Disease Stage       I    211    120        reference      -        reference       -
#>  8:                      II    263    133 0.78 (0.54-1.12)  0.172 0.76 (0.52-1.10)   0.140
#>  9:                     III    241    147 1.19 (0.81-1.73)  0.374 1.14 (0.78-1.68)   0.490
#> 10:                      IV    132     77 1.06 (0.68-1.65)  0.790 1.05 (0.68-1.65)   0.815
```

### **Example 15:** Custom Selection

The `method = "custom"` approach allows explicit specification of
multivariable predictors via the `multi_predictors` argument:

``` r
example15 <- fullfit(
  data = clintrial,
  outcome = "icu_admission",
  predictors = c("age", "sex", "bmi", "smoking", "stage", "treatment"),
  model_type = "glm",
  method = "custom",
  multi_predictors = c("age", "sex", "stage", "treatment"),
  labels = clintrial_labels
)

example15
#> 
#> Fullfit Analysis Results
#> Outcome: icu_admission
#> Model Type: glm
#> Method: custom
#> Predictors Screened: 6
#> Multivariable Predictors: 4
#> 
#>                    Variable   Group      n Events      OR (95% CI)   Uni p     aOR (95% CI) Multi p
#>                      <char>  <char> <char> <char>           <char>  <char>           <char>  <char>
#>  1:             Age (years)       -    850    320 1.02 (1.01-1.03) < 0.001 1.02 (1.01-1.03)   0.001
#>  2:                     Sex  Female    450    173        reference       -        reference       -
#>  3:                            Male    400    147 0.93 (0.70-1.23)   0.611 0.94 (0.71-1.25)   0.680
#>  4: Body Mass Index (kg/m²)       -    838    315 0.99 (0.97-1.02)   0.717                -       -
#>  5:          Smoking Status   Never    337    123        reference       -                -       -
#>  6:                          Former    311    120 1.09 (0.80-1.50)   0.584                -       -
#>  7:                         Current    185     69 1.03 (0.71-1.50)   0.856                -       -
#>  8:           Disease Stage       I    211     76        reference       -        reference       -
#>  9:                              II    263     85 0.85 (0.58-1.24)   0.398 0.85 (0.58-1.24)   0.397
#> 10:                             III    241    105 1.37 (0.94-2.00)   0.103 1.37 (0.93-2.01)   0.113
#> 11:                              IV    132     51 1.12 (0.71-1.75)   0.625 1.12 (0.71-1.76)   0.631
#> 12:         Treatment Group Control    196     64        reference       -        reference       -
#> 13:                          Drug A    292    111 1.26 (0.86-1.85)   0.226 1.26 (0.85-1.85)   0.245
#> 14:                          Drug B    362    145 1.38 (0.96-1.99)   0.085 1.31 (0.90-1.90)   0.159
```

### **Example 16:** Controlling Output Columns

The `columns` parameter controls which results are displayed:

``` r
# Univariable only
example16a <- fullfit(
  data = clintrial,
  outcome = "wound_infection",
  predictors = c("age", "sex", "stage"),
  model_type = "glm",
  columns = "uni"
)

example16a
#> 
#> Fullfit Analysis Results
#> Outcome: wound_infection
#> Model Type: glm
#> Method: screen (p < 0.05)
#> Predictors Screened: 3
#> 
#>    Variable  Group      n Events      OR (95% CI) p-value
#>      <char> <char> <char> <char>           <char>  <char>
#> 1:      age      -    850    167 0.99 (0.98-1.01)   0.447
#> 2:      sex Female    450     89        reference       -
#> 3:            Male    400     78 0.98 (0.70-1.38)   0.919
#> 4:    stage      I    211     51        reference       -
#> 5:              II    263     52 0.77 (0.50-1.20)   0.249
#> 6:             III    241     40 0.62 (0.39-0.99)   0.046
#> 7:              IV    132     23 0.66 (0.38-1.15)   0.141

# Multivariable only
example16b <- fullfit(
  data = clintrial,
  outcome = "wound_infection",
  predictors = c("age", "sex", "stage"),
  model_type = "glm",
  columns = "multi"
)

example16b
#> 
#> Fullfit Analysis Results
#> Outcome: wound_infection
#> Model Type: glm
#> Method: screen (p < 0.05)
#> Predictors Screened: 3
#> Multivariable Predictors: 1
#> 
#>    Variable  Group      n Events      OR (95% CI) p-value
#>      <char> <char> <char> <char>           <char>  <char>
#> 1:    stage      I    211     51        reference       -
#> 2:              II    263     52 0.77 (0.50-1.20)   0.249
#> 3:             III    241     40 0.62 (0.39-0.99)   0.046
#> 4:              IV    132     23 0.66 (0.38-1.15)   0.141
```

### **Example 17:** Survival Analysis

Output tables can use Cox regression for survival outcomes by specifying
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) notation and
`model_type = "coxph"`:

``` r
example17 <- fullfit(
  data = clintrial,
  outcome = "Surv(os_months, os_status)",
  predictors = c("age", "sex", "treatment", "stage", "ecog"),
  model_type = "coxph",
  method = "screen",
  p_threshold = 0.10,
  labels = clintrial_labels
)

example17
#> 
#> Fullfit Analysis Results
#> Outcome: Surv(os_months, os_status)
#> Model Type: coxph
#> Method: screen (p < 0.1)
#> Predictors Screened: 5
#> Multivariable Predictors: 5
#> 
#>                    Variable   Group      n Events      HR (95% CI)   Uni p     aHR (95% CI) Multi p
#>                      <char>  <char> <char> <char>           <char>  <char>           <char>  <char>
#>  1:             Age (years)       -    850    609 1.03 (1.03-1.04) < 0.001 1.04 (1.03-1.04) < 0.001
#>  2:                     Sex  Female    450    298        reference       -        reference       -
#>  3:                            Male    400    311 1.30 (1.11-1.53)   0.001 1.29 (1.10-1.52)   0.002
#>  4:         Treatment Group Control    196    151        reference       -        reference       -
#>  5:                          Drug A    292    184 0.64 (0.52-0.80) < 0.001 0.54 (0.44-0.68) < 0.001
#>  6:                          Drug B    362    274 0.94 (0.77-1.15)   0.567 0.84 (0.69-1.03)   0.095
#>  7:           Disease Stage       I    211    127        reference       -        reference       -
#>  8:                              II    263    172 1.12 (0.89-1.41)   0.337 1.13 (0.90-1.43)   0.288
#>  9:                             III    241    186 1.69 (1.35-2.11) < 0.001 1.93 (1.53-2.43) < 0.001
#> 10:                              IV    132    121 3.18 (2.47-4.09) < 0.001 3.87 (2.98-5.01) < 0.001
#> 11: ECOG Performance Status       0    265    159        reference       -        reference       -
#> 12:                               1    302    212 1.36 (1.11-1.67)   0.003 1.51 (1.22-1.85) < 0.001
#> 13:                               2    238    194 1.86 (1.51-2.29) < 0.001 1.93 (1.56-2.39) < 0.001
#> 14:                               3     37     37 3.06 (2.13-4.38) < 0.001 3.33 (2.31-4.80) < 0.001
```

------------------------------------------------------------------------

## Extended Model Types

The `summata` package supports several additional model types for
specialized analyses. These require external packages that are suggested
dependencies.

### **Example 18:** Negative Binomial

Negative binomial models are appropriate for overdispersed count data
where the variance exceeds the mean—a common violation of the Poisson
assumption. This model type requires the `MASS` package. The `ae_count`
variable in `clintrial` is generated with overdispersion, making it
ideal for this demonstration.

``` r
example18 <- fit(
  data = clintrial,
  outcome = "ae_count",
  predictors = c("age", "sex", "diabetes", "treatment"),
  model_type = "negbin",
  labels = clintrial_labels
)

example18
#> 
#> Multivariable Negative Binomial Model
#> Formula: ae_count ~ age + sex + diabetes + treatment
#> n = 824, Events = 4506
#> 
#>           Variable   Group      n Events     aRR (95% CI) p-value
#>             <char>  <char> <char> <char>           <char>  <char>
#> 1:     Age (years)       -    824  4,506 1.01 (1.01-1.01) < 0.001
#> 2:             Sex  Female    440  2,248        reference       -
#> 3:                    Male    384  2,258 1.10 (0.99-1.23)   0.072
#> 4:        Diabetes      No    630  2,998        reference       -
#> 5:                     Yes    194  1,508 1.56 (1.38-1.76) < 0.001
#> 6: Treatment Group Control    189    826        reference       -
#> 7:                  Drug A    287  1,221 0.96 (0.83-1.12)   0.635
#> 8:                  Drug B    348  2,459 1.52 (1.32-1.75) < 0.001
```

### **Example 19:** Gamma Regression

Gamma regression is appropriate for positive, continuous, right-skewed
outcomes. In the `clintrial` dataset, length of hospital stay
(`los_days`) fits this description. Note that the default link for Gamma
is the *inverse* link, which produces coefficients on a
difficult-to-interpret scale. The log link is generally preferred for
interpretability—exponentiated coefficients then represent
multiplicative effects on the mean.

``` r
example19 <- fit(
  data = clintrial,
  outcome = "los_days",
  predictors = c("age", "ecog", "stage", "treatment"),
  model_type = "glm",
  family = Gamma(link = "log"),
  labels = clintrial_labels
)

example19
#> 
#> Multivariable Gamma Model
#> Formula: los_days ~ age + ecog + stage + treatment
#> n = 822
#> 
#>                    Variable   Group      n Adj. Coefficient (95% CI) p-value
#>                      <char>  <char> <char>                    <char>  <char>
#>  1:             Age (years)       -    822          1.01 (1.01-1.01) < 0.001
#>  2: ECOG Performance Status       0    261                 reference       -
#>  3:                               1    295          1.14 (1.10-1.18) < 0.001
#>  4:                               2    229          1.21 (1.17-1.25) < 0.001
#>  5:                               3     37          1.29 (1.21-1.38) < 0.001
#>  6:           Disease Stage       I    206                 reference       -
#>  7:                              II    257          1.06 (1.02-1.10)   0.002
#>  8:                             III    235          1.13 (1.08-1.17) < 0.001
#>  9:                              IV    124          1.15 (1.10-1.20) < 0.001
#> 10:         Treatment Group Control    190                 reference       -
#> 11:                          Drug A    287          0.95 (0.92-0.99)   0.007
#> 12:                          Drug B    345          1.13 (1.09-1.17) < 0.001
```

### **Example 20:** Random-Intercept Model

Linear mixed-effects models (LMMs) account for hierarchical or clustered
data structures. The `clintrial` dataset includes a `site` variable
representing the study site, making it suitable for demonstrating random
effects. This model type requires the `lme4` package.

Include random intercepts for study site using `(1|site)` notation in
the predictors:

``` r
example20 <- fit(
  data = clintrial,
  outcome = "los_days",
  predictors = c("age", "sex", "treatment", "stage", "(1|site)"),
  model_type = "lmer",
  labels = clintrial_labels
)

example20
#> 
#> Multivariable Linear Mixed Model
#> Formula: los_days ~ age + sex + treatment + stage + (1|site)
#> n = 827
#> 
#>            Variable   Group      n Adj. Coefficient (95% CI) p-value
#>              <char>  <char> <char>                    <char>  <char>
#>  1:     Age (years)       -    827       0.14 (0.11 to 0.16) < 0.001
#>  2:             Sex  Female    441                 reference       -
#>  3:                    Male    386       0.82 (0.28 to 1.37)   0.003
#>  4: Treatment Group Control    190                 reference       -
#>  5:                  Drug A    288    -0.94 (-1.68 to -0.20)   0.013
#>  6:                  Drug B    349       2.41 (1.70 to 3.13) < 0.001
#>  7:   Disease Stage       I    207                 reference       -
#>  8:                      II    259       1.25 (0.51 to 1.98) < 0.001
#>  9:                     III    235       2.43 (1.67 to 3.19) < 0.001
#> 10:                      IV    126       2.96 (2.07 to 3.85) < 0.001
```

### **Example 21:** GLMM for Binary Outcome

Generalized linear mixed-effects models (GLMMs) extend mixed-effects
models to non-normal outcomes (binary, count). This model type also
requires the `lme4` package.

Model 30-day readmission with site-level random effects:

``` r
example21 <- fit(
  data = clintrial,
  outcome = "readmission_30d",
  predictors = c("age", "sex", "diabetes", "treatment", "(1|site)"),
  model_type = "glmer",
  family = "binomial",
  labels = clintrial_labels
)

example21
#> 
#> Multivariable glmerMod Model
#> Formula: readmission_30d ~ age + sex + diabetes + treatment + (1|site)
#> n = 834, Events = 417
#> 
#>           Variable   Group      n Events     aOR (95% CI) p-value
#>             <char>  <char> <char> <char>           <char>  <char>
#> 1:     Age (years)       -    834    417 1.03 (1.02-1.05) < 0.001
#> 2:             Sex  Female    444    207        reference       -
#> 3:                    Male    390    210 1.36 (1.02-1.81)   0.036
#> 4:        Diabetes      No    637    301        reference       -
#> 5:                     Yes    197    116 1.65 (1.17-2.33)   0.004
#> 6: Treatment Group Control    191     88        reference       -
#> 7:                  Drug A    288    126 0.87 (0.59-1.28)   0.471
#> 8:                  Drug B    355    203 1.58 (1.09-2.28)   0.016
```

### **Example 22:** GLMM for Count Outcome

For count outcomes with clustering, use `fu_count` (equidispersed) with
site-level random effects. Standard Poisson GLMMs assume equidispersion:

``` r
example22 <- fit(
  data = clintrial,
  outcome = "fu_count",
  predictors = c("age", "stage", "treatment", "(1|site)"),
  model_type = "glmer",
  family = "poisson",
  labels = clintrial_labels
)

example22
#> 
#> Multivariable glmerMod Model
#> Formula: fu_count ~ age + stage + treatment + (1|site)
#> n = 839
#> 
#>           Variable   Group      n Events     aRR (95% CI) p-value
#>             <char>  <char> <char> <char>           <char>  <char>
#> 1:     Age (years)       -    839   <NA> 1.00 (0.99-1.00)   0.005
#> 2:   Disease Stage       I    207  1,238        reference       -
#> 3:                      II    261  1,638 1.05 (0.97-1.13)   0.233
#> 4:                     III    240  1,638 1.12 (1.04-1.21)   0.002
#> 5:                      IV    131  1,003 1.27 (1.17-1.38) < 0.001
#> 6: Treatment Group Control    191  1,169        reference       -
#> 7:                  Drug A    290  1,910 1.07 (1.00-1.16)   0.053
#> 8:                  Drug B    358  2,438 1.11 (1.03-1.19)   0.005
```

### **Example 23:** Cox Mixed-Effects Model

Cox mixed-effects models account for within-cluster correlation in
survival outcomes. This is useful when patients are nested within sites
or physicians and outcomes may be correlated within clusters. This model
type requires the `coxme` package. Model overall survival with
site-level random effects:

``` r
example23 <- fit(
  data = clintrial,
  outcome = "Surv(os_months, os_status)",
  predictors = c("age", "sex", "treatment", "stage", "(1|site)"),
  model_type = "coxme",
  labels = clintrial_labels
)

example23
#> 
#> Multivariable Mixed Effects Cox Model
#> Formula: Surv(os_months, os_status) ~ age + sex + treatment + stage + (1|site)
#> n = 847, Events = 606
#> 
#>            Variable   Group      n Events     aHR (95% CI) p-value
#>              <char>  <char> <char> <char>           <char>  <char>
#>  1:     Age (years)       -    847    606 1.04 (1.03-1.05) < 0.001
#>  2:             Sex  Female    450    298        reference       -
#>  3:                    Male    400    311 1.36 (1.16-1.60) < 0.001
#>  4: Treatment Group Control    196    151        reference       -
#>  5:                  Drug A    292    184 0.57 (0.46-0.72) < 0.001
#>  6:                  Drug B    362    274 0.93 (0.76-1.14)   0.492
#>  7:   Disease Stage       I    211    127        reference       -
#>  8:                      II    263    172 1.17 (0.93-1.48)   0.178
#>  9:                     III    241    186 2.04 (1.62-2.58) < 0.001
#> 10:                      IV    132    121 4.09 (3.16-5.30) < 0.001
```

### **Example 24:** Quasibinomial for Overdispersed Binary Data

When binary or count data exhibit overdispersion (residual deviance \>\>
residual degrees of freedom), quasi-likelihood models provide more
appropriate standard errors.

``` r
example24 <- fit(
  data = clintrial,
  outcome = "any_complication",
  predictors = c("age", "sex", "diabetes", "stage"),
  model_type = "glm",
  family = "quasibinomial",
  labels = clintrial_labels
)

example24
#> 
#> Multivariable Quasi-Binomial Model
#> Formula: any_complication ~ age + sex + diabetes + stage
#> n = 834, Events = 468
#> 
#>         Variable  Group      n Events     aOR (95% CI) p-value
#>           <char> <char> <char> <char>           <char>  <char>
#> 1:   Age (years)      -    834    468 1.01 (1.00-1.02)   0.155
#> 2:           Sex Female    444    236        reference       -
#> 3:                 Male    390    232 1.33 (1.00-1.76)   0.048
#> 4:      Diabetes     No    637    335        reference       -
#> 5:                  Yes    197    133 1.93 (1.37-2.71) < 0.001
#> 6: Disease Stage      I    207    117        reference       -
#> 7:                   II    261    131 0.75 (0.52-1.09)   0.132
#> 8:                  III    237    144 1.22 (0.83-1.79)   0.317
#> 9:                   IV    129     76 1.07 (0.68-1.69)   0.764
```

### **Example 25:** Conditional Logistic Model for Matched Data

Conditional logistic regression is used for matched case-control studies
where cases and controls are matched within strata. The stratification
variable should be specified using the `strata` parameter. This model
type requires the `survival` package.

``` r
# Matched case-control dataset (100 pairs)
set.seed(456)
n_pairs <- 100
matched_data <- do.call(rbind, lapply(1:n_pairs, function(i) {
  base_bmi <- rnorm(1, 27, 4)
  data.frame(
    match_id = i,
    case = c(1, 0),
    smoking = factor(c(rbinom(1, 1, 0.45), rbinom(1, 1, 0.30)), 
                     levels = c(0, 1), labels = c("No", "Yes")),
    diabetes = factor(c(rbinom(1, 1, 0.30), rbinom(1, 1, 0.20)), 
                      levels = c(0, 1), labels = c("No", "Yes")),
    bmi = c(base_bmi + rnorm(1, 1.5, 2), base_bmi + rnorm(1, -0.5, 2))
  )
}))

example25 <- fit(
  data = matched_data,
  outcome = "case",
  predictors = c("smoking", "diabetes", "bmi"),
  model_type = "clogit",
  strata = "match_id"
)

example25
#> 
#> Multivariable Conditional Logistic Model
#> Formula: case ~ smoking + diabetes + bmi + strata( match_id )
#> n = 124, Events = 100
#> 
#>    Variable  Group      n Events     aHR (95% CI) p-value
#>      <char> <char> <char> <char>           <char>  <char>
#> 1:  smoking     No    124    100        reference       -
#> 2:             Yes     76    100 3.46 (1.48-8.12)   0.004
#> 3: diabetes     No    140    100        reference       -
#> 4:             Yes     60    100 0.90 (0.40-2.03)   0.795
#> 5:      bmi      -    200    100 1.61 (1.30-2.00) < 0.001
```

The output shows odds ratios conditional on the matching strata.

------------------------------------------------------------------------

## Exporting Results

Regression tables can be exported to various formats:

``` r
# Microsoft Word
table2docx(
  table = example13,
  file = file.path(tempdir(), "Table2_Regression.docx"),
  caption = "Table 2. Univariable and Multivariable Analysis"
)

# PDF
table2pdf(
  table = example13,
  file = file.path(tempdir(), "Table2_Regression.pdf"),
  caption = "Table 2. Regression Results"
)
```

See the [Table
Export](https://phmcc.github.io/summata/articles/table_export.md)
vignette for comprehensive documentation.

------------------------------------------------------------------------

## Best Practices

### Variable Selection Strategy

1.  **Exploratory phase**: Use
    [`uniscreen()`](https://phmcc.github.io/summata/reference/uniscreen.md)
    with liberal *p*-threshold (e.g., 0.20)
2.  **Model building**: Use
    [`fullfit()`](https://phmcc.github.io/summata/reference/fullfit.md)
    with `method = "screen"` or `method = "custom"`
3.  **Confirmatory analysis**: Use
    [`fit()`](https://phmcc.github.io/summata/reference/fit.md) with
    pre-specified predictors
4.  **Sensitivity analysis**: Compare specifications with
    [`compfit()`](https://phmcc.github.io/summata/reference/compfit.md)
    (*see* [Model
    Comparison](https://phmcc.github.io/summata/articles/model_comparison.md))

### Choosing the Right Model

| Outcome Type  | Characteristics          | Recommended Model          |
|:--------------|:-------------------------|:---------------------------|
| Binary        | Independent observations | `glm` with `binomial`      |
| Binary        | Overdispersed            | `glm` with `quasibinomial` |
| Binary        | Clustered/hierarchical   | `glmer` with `binomial`    |
| Binary        | Matched case-control     | `clogit`                   |
| Count         | Mean ≈ variance          | `glm` with `poisson`       |
| Count         | Variance \> mean         | `negbin` or `quasipoisson` |
| Count         | Clustered                | `glmer` with `poisson`     |
| Continuous    | Symmetric, normal errors | `lm`                       |
| Continuous    | Positive, right-skewed   | `glm` with `Gamma`         |
| Continuous    | Clustered                | `lmer`                     |
| Time-to-event | Independent              | `coxph`                    |
| Time-to-event | Clustered                | `coxme`                    |

### Interpreting Results

The comparison between univariable and multivariable estimates reveals
confounding:

- **Large differences**: Substantial confounding present
- **Similar estimates**: Association robust to adjustment
- **Sign reversal**: Simpson’s paradox; investigate carefully

### Sample Size Considerations

Adequate events per predictor are required for stable estimation:

- Logistic/Cox models: ≥10 events per predictor (rule of thumb)
- Linear models: ≥10–20 observations per predictor
- Mixed-effects models: ≥5–10 clusters recommended
- Use `method = "screen"` to reduce predictors when sample size is
  limited

------------------------------------------------------------------------

## Common Issues

### Missing Data

Regression functions use complete-case analysis by default:

``` r
# Check missingness
sapply(clintrial[, c("age", "sex", "stage")], function(x) sum(is.na(x)))

# Create complete-case dataset explicitly
complete_data <- na.omit(clintrial[, c("readmission_30d", "age", "sex", "stage")])
```

### Factor Reference Levels

Ensure reference categories are set appropriately:

``` r
# Set specific reference level
clintrial$stage <- relevel(factor(clintrial$stage), ref = "I")
```

### Convergence Issues

For models that fail to converge:

``` r
# Access model for diagnostics
result <- fit(data, outcome, predictors, model_type = "glm")
model <- attr(result, "model")

# Check convergence
model$converged

# Large coefficients may indicate separation
coef(model)
```

### Overdispersion Diagnostics

To check for overdispersion in Poisson or binomial models:

``` r
# Fit Poisson model
result <- fit(data, outcome, predictors, model_type = "glm", family = "poisson")
model <- attr(result, "model")

# Dispersion estimate (should be ~1 for no overdispersion)
sum(residuals(model, type = "pearson")^2) / model$df.residual
```

If the dispersion is substantially greater than 1, consider using
`quasipoisson`, `quasibinomial`, or `negbin` instead.

### Mixed-Effects Model Convergence

Mixed-effects models can be sensitive to starting values and
optimization:

``` r
# Increase iterations
result <- fit(
  data = clintrial,
  outcome = "readmission_30d",
  predictors = c("age", "treatment", "(1|site)"),
  model_type = "glmer",
  family = "binomial",
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)
```

### Multinomial and Ordinal Outcomes

The `summata` package currently supports binary outcomes (via logistic
regression) but not multinomial or ordinal outcomes with more than two
categories. For categorical outcomes with more than two levels, use
dedicated packages such as
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) for
multinomial regression or
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) for ordinal
regression.

------------------------------------------------------------------------

## Further Reading

- [Descriptive
  Tables](https://phmcc.github.io/summata/articles/descriptive_tables.md):
  [`desctable()`](https://phmcc.github.io/summata/reference/desctable.md)
  for baseline characteristics
- [Model
  Comparison](https://phmcc.github.io/summata/articles/model_comparison.md):
  [`compfit()`](https://phmcc.github.io/summata/reference/compfit.md)
  for comparing models
- [Table
  Export](https://phmcc.github.io/summata/articles/table_export.md):
  Export to PDF, Word, and other formats
- [Forest
  Plots](https://phmcc.github.io/summata/articles/forest_plots.md):
  Visualization of regression results
- [Multivariate
  Regression](https://phmcc.github.io/summata/articles/multivariate_regression.md):
  [`multifit()`](https://phmcc.github.io/summata/reference/multifit.md)
  for multi-outcome analysis
- [Advanced
  Workflows](https://phmcc.github.io/summata/articles/advanced_workflows.md):
  Interactions and mixed-effects models
