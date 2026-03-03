# Multivariate Regression Analysis

Performs regression analyses of a single predictor (exposure) across
multiple outcomes. This function is designed for studies where a single
exposure variable is tested against multiple endpoints, such as
complication screening, biomarker associations, or phenome-wide
association studies. Returns publication-ready formatted results with
optional covariate adjustment. Supports interactions, mixed-effects
models, stratification, and clustered standard errors.

## Usage

``` r
multifit(
  data,
  outcomes,
  predictor,
  covariates = NULL,
  interactions = NULL,
  random = NULL,
  strata = NULL,
  cluster = NULL,
  model_type = "glm",
  family = "binomial",
  columns = "adjusted",
  p_threshold = 1,
  conf_level = 0.95,
  show_n = TRUE,
  show_events = TRUE,
  digits = 2,
  p_digits = 3,
  labels = NULL,
  predictor_label = NULL,
  include_predictor = TRUE,
  keep_models = FALSE,
  exponentiate = NULL,
  parallel = TRUE,
  n_cores = NULL,
  number_format = NULL,
  verbose = NULL,
  ...
)
```

## Arguments

- data:

  Data frame or data.table containing the analysis dataset. The function
  automatically converts data frames to data.tables for efficient
  processing.

- outcomes:

  Character vector of outcome variable names to analyze. Each outcome is
  tested in its own model with the predictor. For time-to-event
  analysis, use [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
  syntax for the outcome variable (*e.g.,*
  `c("Surv(time1, status1)", "Surv(time2, status2)")`).

- predictor:

  Character string specifying the predictor (exposure) variable name.
  This variable is tested against each outcome. Can be continuous or
  categorical (factor).

- covariates:

  Optional character vector of covariate variable names to include in
  adjusted models. When specified, models are fit as
  `outcome ~ predictor + covariate1 + covariate2 + ...`, and only the
  predictor effect is reported. Default is `NULL` (unadjusted models).

- interactions:

  Optional character vector of interaction terms to include in adjusted
  models, using colon notation (*e.g.,* `c("predictor:sex")`).
  Interactions involving the predictor will have their effects extracted
  and reported. Default is `NULL`.

- random:

  Optional character string specifying random effects formula for mixed
  effects models (*e.g.,* `"(1|hospital)"` or `"(1|site/patient)"`).
  Required when `model_type` is `"glmer"`, `"lmer"`, or `"coxme"` unless
  random effects are included in the `covariates` vector. Alternatively,
  random effects can be included directly in the `covariates` vector
  using the same syntax (*e.g.,*
  `covariates = c("age", "sex", "(1|site)")`). Default is `NULL`.

- strata:

  Optional character string naming the stratification variable for Cox
  or conditional logistic models. Creates separate baseline hazards for
  each stratum. Default is `NULL`.

- cluster:

  Optional character string naming the clustering variable for Cox
  models. Computes robust clustered standard errors. Default is `NULL`.

- model_type:

  Character string specifying the type of regression model to fit.
  Options include:

  - `"glm"` - Generalized linear model (default). Supports multiple
    distributions via the `family` parameter including logistic,
    Poisson, Gamma, Gaussian, and quasi-likelihood models.

  - `"lm"` - Linear regression for continuous outcomes with normally
    distributed errors.

  - `"coxph"` - Cox proportional hazards model for time-to-event
    survival analysis. Requires
    [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) outcome
    syntax.

  - `"clogit"` - Conditional logistic regression for matched
    case-control studies.

  - `"negbin"` - Negative binomial regression for overdispersed count
    data (requires MASS package). Estimates an additional dispersion
    parameter compared to Poisson regression.

  - `"glmer"` - Generalized linear mixed-effects model for hierarchical
    or clustered data with non-normal outcomes (requires lme4 package
    and `random` parameter).

  - `"lmer"` - Linear mixed-effects model for hierarchical or clustered
    data with continuous outcomes (requires lme4 package and `random`
    parameter).

  - `"coxme"` - Cox mixed-effects model for clustered survival data
    (requires coxme package and `random` parameter).

- family:

  For GLM and GLMER models, specifies the error distribution and link
  function. Can be a character string, a family function, or a family
  object. Ignored for non-GLM/GLMER models.

  **Binary/Binomial outcomes:**

  - `"binomial"` or
    [`binomial()`](https://rdrr.io/r/stats/family.html) - Logistic
    regression for binary outcomes (0/1, TRUE/FALSE). Returns odds
    ratios (OR). Default.

  - `"quasibinomial"` or
    [`quasibinomial()`](https://rdrr.io/r/stats/family.html) - Logistic
    regression with overdispersion. Use when residual deviance \>\>
    residual df.

  - `binomial(link = "probit")` - Probit regression (normal CDF link).

  - `binomial(link = "cloglog")` - Complementary log-log link for
    asymmetric binary outcomes.

  **Count outcomes:**

  - `"poisson"` or [`poisson()`](https://rdrr.io/r/stats/family.html) -
    Poisson regression for count data. Returns rate ratios (RR). Assumes
    mean = variance.

  - `"quasipoisson"` or
    [`quasipoisson()`](https://rdrr.io/r/stats/family.html) - Poisson
    regression with overdispersion. Use when variance \> mean.

  **Continuous outcomes:**

  - `"gaussian"` or
    [`gaussian()`](https://rdrr.io/r/stats/family.html) -
    Normal/Gaussian distribution for continuous outcomes. Equivalent to
    linear regression.

  - `gaussian(link = "log")` - Log-linear model for positive continuous
    outcomes. Returns multiplicative effects.

  **Positive continuous outcomes:**

  - `"Gamma"` or [`Gamma()`](https://rdrr.io/r/stats/family.html) -
    Gamma distribution for positive, right-skewed continuous data
    (*e.g.,* costs, lengths of stay). Default log link.

  - `Gamma(link = "inverse")` - Gamma with inverse (canonical) link.

  - `Gamma(link = "identity")` - Gamma with identity link for additive
    effects on positive outcomes.

  - `"inverse.gaussian"` or
    [`inverse.gaussian()`](https://rdrr.io/r/stats/family.html) -
    Inverse Gaussian for positive, highly right-skewed data.

  For negative binomial regression (overdispersed counts), use
  `model_type = "negbin"` instead of the `family` parameter.

  See [`family`](https://rdrr.io/r/stats/family.html) for additional
  details and options.

- columns:

  Character string specifying which result columns to display when both
  unadjusted and adjusted models are fit (*i.e.,* when `covariates` is
  specified):

  - `"adjusted"` - Show only adjusted (covariate-controlled) results
    \[default\]

  - `"unadjusted"` - Show only unadjusted (crude) results

  - `"both"` - Show both unadjusted and adjusted results side-by-side

  Ignored when `covariates = NULL`.

- p_threshold:

  Numeric value between 0 and 1 specifying a *p*-value threshold for
  filtering results. Only outcomes with *p*-value less than or equal to
  the threshold are included in the output. Default is 1 (no filtering;
  all outcomes returned).

- conf_level:

  Numeric confidence level for confidence intervals. Must be between 0
  and 1. Default is 0.95 (95% confidence intervals).

- show_n:

  Logical. If `TRUE`, includes the sample size column in the output
  table. Default is `TRUE`.

- show_events:

  Logical. If `TRUE`, includes the events column in the output table
  (relevant for survival and logistic regression). Default is `TRUE`.

- digits:

  Integer specifying the number of decimal places for effect estimates
  (OR, HR, RR, coefficients). Default is 2.

- p_digits:

  Integer specifying the number of decimal places for *p*-values. Values
  smaller than `10^(-p_digits)` are displayed as `"< 0.001"` (for
  `p_digits = 3`), `"< 0.0001"` (for `p_digits = 4`), etc. Default is 3.

- labels:

  Named character vector or list providing custom display labels for
  variables. Can include labels for outcomes, predictors, and
  covariates. Names should match variable names, values are the display
  labels. Labels are applied to: (1) outcome names in the Outcome
  column, (2) predictor variable name when displayed, and (3) variable
  names in formatted interaction terms. Variables not in `labels` use
  their original names. Default is `NULL`.

- predictor_label:

  Optional character string providing a custom display label for the
  predictor variable. Takes precedence over `labels` for the predictor.
  Default is `NULL` (uses label from `labels` or original name).

- include_predictor:

  Logical. If `TRUE` (default), includes the predictor column in the
  output table showing which level of a factor predictor is being
  compared. If `FALSE`, omits the predictor column, which may be useful
  when the predictor information will be explained in a table caption or
  figure legend.

- keep_models:

  Logical. If `TRUE`, stores all fitted model objects in the output as
  an attribute. Models are accessible via `attr(result, "models")`.
  Default is `FALSE`.

- exponentiate:

  Logical. Whether to exponentiate coefficients (display OR/HR/RR
  instead of log odds/log hazards). Default is `NULL`, which
  automatically displays raw coefficients for linear models and
  exponentiates for logistic, Poisson, and Cox models.

- parallel:

  Logical. If `TRUE` (default), fits models in parallel for improved
  performance with many outcomes.

- n_cores:

  Integer specifying the number of CPU cores to use for parallel
  processing. Default is `NULL` (auto-detect: uses `detectCores() - 1`).
  Ignored when `parallel = FALSE`.

- number_format:

  Character string or two-element character vector controlling thousand
  and decimal separators in formatted output. Named presets:

  - `"us"` - Comma thousands, period decimal: `1,234.56` \[default\]

  - `"eu"` - Period thousands, comma decimal: `1.234,56`

  - `"space"` - Thin-space thousands, period decimal: `1 234.56` (SI/ISO
    31-0)

  - `"none"` - No thousands separator: `1234.56`

  Or provide a custom two-element vector `c(big.mark, decimal.mark)`,
  *e.g.*, `c("'", ".")` for Swiss-style: `1'234.56`.

  When `NULL` (default), uses
  `getOption("summata.number_format", "us")`. Set the global option once
  per session to avoid passing this argument repeatedly:

          options(summata.number_format = "eu")
        

- verbose:

  Logical. If `TRUE`, displays model fitting warnings (*e.g.,* singular
  fit, convergence issues). If `FALSE` (default), routine fitting
  messages are suppressed while unexpected warnings are preserved. When
  `NULL`, uses `getOption("summata.verbose", FALSE)`.

- ...:

  Additional arguments passed to the underlying model fitting functions.

## Value

A data.table with S3 class `"multifit_result"` containing formatted
multivariate regression results. The table structure includes:

- Outcome:

  Character. Outcome variable name or custom label

- Predictor:

  Character. For factor predictors: formatted as "Variable (Level)"
  showing the level being compared to reference. For binary variables
  where the non-reference level is an affirmative value (Yes, 1, True,
  Present, Positive, +), shows just "Variable". For continuous
  predictors: the variable name. For interactions: the formatted
  interaction term (*e.g.,* "Treatment (Drug A) × Sex (Male)")

- n:

  Integer. Sample size used in the model (if `show_n = TRUE`)

- Events:

  Integer. Number of events (if `show_events = TRUE`)

- OR/HR/RR/Coefficient (95% CI):

  Character. Unadjusted effect estimate with CI (if
  `columns = "unadjusted"` or `"both"`)

- aOR/aHR/aRR/Adj. Coefficient (95% CI):

  Character. Adjusted effect estimate with CI (if `columns = "adjusted"`
  or `"both"`)

- Uni *p* / Multi *p* / *p*-value:

  Character. Formatted *p*-value(s). Column names depend on `columns`
  setting

The returned object includes the following attributes accessible via
[`attr()`](https://rdrr.io/r/base/attr.html):

- raw_data:

  data.table. Unformatted numeric results with separate columns for
  effect estimates, standard errors, confidence intervals, and
  *p*-values. Suitable for custom analysis or visualization

- models:

  list (if `keep_models = TRUE`). Named list of fitted model objects,
  with outcome names as list names. Each element contains `$unadjusted`
  and/or `$adjusted` models depending on settings

- predictor:

  Character. The predictor variable name

- outcomes:

  Character vector. The outcome variable names

- covariates:

  Character vector or `NULL`. The covariate variable names

- interactions:

  Character vector or `NULL`. The interaction terms

- random:

  Character or `NULL`. The random effects formula

- strata:

  Character or `NULL`. The stratification variable

- cluster:

  Character or `NULL`. The clustering variable

- model_type:

  Character. The regression model type used

- columns:

  Character. Which columns were displayed

- analysis_type:

  Character. `"multi_outcome"` to identify analysis type

- significant:

  Character vector. Names of outcomes with *p* \< 0.05 for the predictor
  (uses adjusted *p*-values when available)

## Details

**Analysis Approach:**

The function implements a multivariate (multi-outcome) screening
workflow that inverts the typical regression paradigm:

1.  For each outcome in `outcomes`, fits a separate model with the
    predictor as the main exposure

2.  If `covariates` specified, fits adjusted model:
    `outcome ~ predictor + covariates + interactions`

3.  Extracts only the predictor effect(s) from each model, ignoring
    covariate coefficients

4.  Combines results into a single table for comparison across outcomes

5.  Optionally filters by *p*-value threshold

This is conceptually opposite to
[`uniscreen()`](https://phmcc.github.io/summata/reference/uniscreen.md),
which tests multiple predictors against a single outcome. Use
`multifit()` when you have one exposure of interest and want to screen
across multiple endpoints.

**When to Use Multivariate Regression Analysis:**

- **Complication screening**: Test one exposure (*e.g.,* operative time,
  BMI, biomarker level) against multiple postoperative complications

- **Treatment effects**: Test one treatment against multiple efficacy
  and safety endpoints simultaneously

- **Biomarker studies**: Test one biomarker against multiple clinical
  outcomes to understand its prognostic value

- **Phenome-wide association studies (PheWAS)**: Test genetic variants
  or exposures against many phenotypes

- **Risk factor profiling**: Understand how one risk factor relates to a
  spectrum of outcomes

**Handling Categorical Predictors:**

When the predictor is a factor variable with multiple levels:

- Each non-reference level gets its own row for each outcome

- Reference category is determined by factor level ordering

- The Predictor column shows "Variable (Level)" format (*e.g.,*
  "Treatment (Drug A)", "Treatment (Drug B)")

- For binary variables with affirmative non-reference levels (Yes, 1,
  True, Present, Positive, +), shows just "Variable" (*e.g.,* "Diabetes"
  instead of "Diabetes (Yes)")

- Effect estimates compare each level to the reference

**Adjusted vs. Unadjusted Results:**

When `covariates` is specified, the function fits both models but only
extracts predictor effects:

- `columns = "adjusted"`: Reports only covariate-adjusted effects.
  Column labeled "aOR/aHR," *etc.*

- `columns = "unadjusted"`: Reports only crude effects. Column labeled
  "OR/HR," *etc.*

- `columns = "both"`: Reports both side-by-side. Useful for identifying
  confounding (large change in effect) or independent effects (similar
  estimates)

**Interaction Terms:**

When `interactions` includes terms involving the predictor:

- Main effect of predictor is always reported

- Interaction effects are extracted and displayed with formatted names

- Format: `Variable (Level) × Variable (Level)` using multiplication
  sign notation

- Useful for testing effect modification (*e.g.,* does treatment effect
  differ by sex?)

**Mixed-Effects Models:**

For clustered or hierarchical data (*e.g.,* patients within hospitals):

- Use `model_type = "glmer"` with `random = "(1|cluster)"` for random
  intercept models

- Nested random effects: `random = "(1|site/patient)"`

- Crossed random effects: `random = "(1|site) + (1|doctor)"`

- For survival outcomes, use `model_type = "coxme"`

**Stratification and Clustering (Cox models):**

For Cox proportional hazards models:

- `strata`: Creates separate baseline hazards for each stratum level.
  Use when hazards are non-proportional across strata but stratum
  effects do not need to be estimated

- `cluster`: Computes robust (sandwich) standard errors accounting for
  within-cluster correlation. Alternative to mixed effects when only
  robust SEs are needed

**Filtering based on *p*-value:**

The `p_threshold` parameter filters results after fitting all models:

- Only outcomes with *p* less than or equal to the threshold are
  retained in output

- For factor predictors, outcome is kept if any level is significant

- Useful for focusing on significant associations in exploratory
  analyses

- Default is 1 (no filtering) - recommended for confirmatory analyses

**Outcome Homogeneity:**

All outcomes in a single `multifit()` call should be of the same type
(all binary, all continuous, or all survival). Mixing outcome types
produces tables with incompatible effect measures (*e.g.,* odds ratios
alongside regression coefficients), which can mislead readers. The
function validates outcome compatibility and issues a warning when mixed
types are detected.

For analyses involving multiple outcome types, run separate `multifit()`
calls for each type:

    # Binary outcomes
    binary_results <- multifit(data, outcomes = c("death", "readmission"),
                               predictor = "treatment", model_type = "glm")

    # Continuous outcomes
    continuous_results <- multifit(data, outcomes = c("los_days", "cost"),
                                   predictor = "treatment", model_type = "lm")

**Effect Measures by Model Type:**

- **Logistic** (`model_type = "glm"`, `family = "binomial"`): Odds
  ratios (OR/aOR)

- **Cox** (`model_type = "coxph"`): Hazard ratios (HR/aHR)

- **Poisson** (`model_type = "glm"`, `family = "poisson"`): Rate ratios
  (RR/aRR)

- **Linear** (`model_type = "lm"`): Coefficient estimates

- **Mixed effects**: Same as fixed-effects counterparts

**Memory and Performance:**

- `parallel = TRUE` (default) uses multiple cores for faster fitting

- `keep_models = FALSE` (default) discards model objects to save memory

- For many outcomes, parallel processing provides substantial speedup

- Set `keep_models = TRUE` only when you need model diagnostics

## See also

[`uniscreen`](https://phmcc.github.io/summata/reference/uniscreen.md)
for screening multiple predictors against one outcome,
[`multiforest`](https://phmcc.github.io/summata/reference/multiforest.md)
for creating forest plots from multifit results,
[`fit`](https://phmcc.github.io/summata/reference/fit.md) for
single-outcome regression with full coefficient output,
[`fullfit`](https://phmcc.github.io/summata/reference/fullfit.md) for
complete univariable-to-multivariable workflow

Other regression functions:
[`compfit()`](https://phmcc.github.io/summata/reference/compfit.md),
[`fit()`](https://phmcc.github.io/summata/reference/fit.md),
[`fullfit()`](https://phmcc.github.io/summata/reference/fullfit.md),
[`print.compfit_result()`](https://phmcc.github.io/summata/reference/print.compfit_result.md),
[`print.fit_result()`](https://phmcc.github.io/summata/reference/print.fit_result.md),
[`print.fullfit_result()`](https://phmcc.github.io/summata/reference/print.fullfit_result.md),
[`print.multifit_result()`](https://phmcc.github.io/summata/reference/print.multifit_result.md),
[`print.uniscreen_result()`](https://phmcc.github.io/summata/reference/print.uniscreen_result.md),
[`uniscreen()`](https://phmcc.github.io/summata/reference/uniscreen.md)

## Examples

``` r
# Load example data
data(clintrial)
data(clintrial_labels)

# Example 1: Basic multivariate analysis (unadjusted)
# Test treatment effect on multiple binary outcomes
result1 <- multifit(
    data = clintrial,
    outcomes = c("surgery", "pfs_status", "os_status"),
    predictor = "treatment",
    labels = clintrial_labels,
    parallel = FALSE
)
print(result1)
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 3
#> Model Type: glm
#> Display: unadjusted
#> 
#>                       Outcome                Predictor      n Events      OR (95% CI) p-value
#>                        <char>                   <char> <char> <char>           <char>  <char>
#> 1:         Surgical Resection Treatment Group (Drug A)    292    173 1.58 (1.10-2.27)   0.014
#> 2:         Surgical Resection Treatment Group (Drug B)    362    103 0.43 (0.30-0.62) < 0.001
#> 3: Progression or Death Event Treatment Group (Drug A)    292    227 0.44 (0.26-0.74)   0.002
#> 4: Progression or Death Event Treatment Group (Drug B)    362    322 1.02 (0.59-1.77)   0.950
#> 5:                Death Event Treatment Group (Drug A)    292    184 0.51 (0.34-0.76)   0.001
#> 6:                Death Event Treatment Group (Drug B)    362    274 0.93 (0.62-1.40)   0.721
# Shows odds ratios comparing Drug A and Drug B to Control

# \donttest{

# Example 2: Adjusted analysis with covariates
# Adjust for age, sex, and disease stage
result2 <- multifit(
    data = clintrial,
    outcomes = c("surgery", "pfs_status", "os_status"),
    predictor = "treatment",
    covariates = c("age", "sex", "stage"),
    labels = clintrial_labels,
    parallel = FALSE
)
print(result2)
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 3
#> Model Type: glm
#> Covariates: age, sex, stage
#> Display: adjusted
#> 
#>                       Outcome                Predictor      n Events     aOR (95% CI) p-value
#>                        <char>                   <char> <char> <char>           <char>  <char>
#> 1:         Surgical Resection Treatment Group (Drug A)    292    173 1.84 (1.23-2.75)   0.003
#> 2:         Surgical Resection Treatment Group (Drug B)    361    102 0.45 (0.30-0.66) < 0.001
#> 3: Progression or Death Event Treatment Group (Drug A)    292    227 0.36 (0.21-0.63) < 0.001
#> 4: Progression or Death Event Treatment Group (Drug B)    361    321 0.77 (0.43-1.38)   0.374
#> 5:                Death Event Treatment Group (Drug A)    292    184 0.44 (0.28-0.68) < 0.001
#> 6:                Death Event Treatment Group (Drug B)    361    273 0.74 (0.48-1.16)   0.192
# Shows adjusted odds ratios (aOR)

# Example 3: Compare unadjusted and adjusted results
result3 <- multifit(
    data = clintrial,
    outcomes = c("surgery", "pfs_status", "os_status"),
    predictor = "treatment",
    covariates = c("age", "sex", "stage"),
    columns = "both",
    labels = clintrial_labels,
    parallel = FALSE
)
print(result3)
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 3
#> Model Type: glm
#> Covariates: age, sex, stage
#> Display: both
#> 
#>                       Outcome                Predictor      n Events      OR (95% CI)   Uni p     aOR (95% CI) Multi p
#>                        <char>                   <char> <char> <char>           <char>  <char>           <char>  <char>
#> 1:                Death Event Treatment Group (Drug A)    292    184 0.51 (0.34-0.76)   0.001 0.44 (0.28-0.68) < 0.001
#> 2:                Death Event Treatment Group (Drug B)    362    274 0.93 (0.62-1.40)   0.721 0.74 (0.48-1.16)   0.192
#> 3: Progression or Death Event Treatment Group (Drug A)    292    227 0.44 (0.26-0.74)   0.002 0.36 (0.21-0.63) < 0.001
#> 4: Progression or Death Event Treatment Group (Drug B)    362    322 1.02 (0.59-1.77)   0.950 0.77 (0.43-1.38)   0.374
#> 5:         Surgical Resection Treatment Group (Drug A)    292    173 1.58 (1.10-2.27)   0.014 1.84 (1.23-2.75)   0.003
#> 6:         Surgical Resection Treatment Group (Drug B)    362    103 0.43 (0.30-0.62) < 0.001 0.45 (0.30-0.66) < 0.001
# Useful for identifying confounding effects

# Example 4: Continuous predictor across outcomes
# Test age effect on multiple outcomes
result4 <- multifit(
    data = clintrial,
    outcomes = c("surgery", "pfs_status", "os_status"),
    predictor = "age",
    covariates = c("sex", "treatment", "stage"),
    labels = clintrial_labels,
    parallel = FALSE
)
print(result4)
#> 
#> Multivariate Analysis Results
#> Predictor: age
#> Outcomes: 3
#> Model Type: glm
#> Covariates: sex, treatment, stage
#> Display: adjusted
#> 
#>                       Outcome      n Events     aOR (95% CI) p-value
#>                        <char> <char> <char>           <char>  <char>
#> 1:         Surgical Resection    847    368 0.96 (0.94-0.97) < 0.001
#> 2: Progression or Death Event    847    721 1.02 (1.01-1.04)   0.006
#> 3:                Death Event    847    606 1.05 (1.04-1.07) < 0.001
# One row per outcome for continuous predictor

# Example 5: Cox regression for survival outcomes
library(survival)
cox_result <- multifit(
    data = clintrial,
    outcomes = c("Surv(pfs_months, pfs_status)", 
                 "Surv(os_months, os_status)"),
    predictor = "treatment",
    covariates = c("age", "sex", "stage"),
    model_type = "coxph",
    labels = clintrial_labels,
    parallel = FALSE
)
print(cox_result)
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 2
#> Model Type: coxph
#> Covariates: age, sex, stage
#> Display: adjusted
#> 
#>                               Outcome                Predictor      n Events     aHR (95% CI) p-value
#>                                <char>                   <char> <char> <char>           <char>  <char>
#> 1: Progression-Free Survival (months) Treatment Group (Drug A)    721    721 0.54 (0.44-0.66) < 0.001
#> 2: Progression-Free Survival (months) Treatment Group (Drug B)    721    721 0.82 (0.68-0.98)   0.033
#> 3:          Overall Survival (months) Treatment Group (Drug A)    606    606 0.56 (0.45-0.70) < 0.001
#> 4:          Overall Survival (months) Treatment Group (Drug B)    606    606 0.83 (0.67-1.01)   0.062
# Returns hazard ratios (HR/aHR)

# Example 6: Cox with stratification by site
cox_strat <- multifit(
    data = clintrial,
    outcomes = c("Surv(os_months, os_status)"),
    predictor = "treatment",
    covariates = c("age", "sex"),
    strata = "site",
    model_type = "coxph",
    labels = clintrial_labels,
    parallel = FALSE
)
print(cox_strat)
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 1
#> Model Type: coxph
#> Covariates: age, sex
#> Strata: site
#> Display: adjusted
#> 
#>                      Outcome                Predictor      n Events     aHR (95% CI) p-value
#>                       <char>                   <char> <char> <char>           <char>  <char>
#> 1: Overall Survival (months) Treatment Group (Drug A)    609    609 0.62 (0.50-0.77) < 0.001
#> 2: Overall Survival (months) Treatment Group (Drug B)    609    609 1.00 (0.81-1.22)   0.976

# Example 7: Cox with clustered standard errors
cox_cluster <- multifit(
    data = clintrial,
    outcomes = c("Surv(os_months, os_status)"),
    predictor = "treatment",
    covariates = c("age", "sex", "stage"),
    cluster = "site",
    model_type = "coxph",
    labels = clintrial_labels,
    parallel = FALSE
)
print(cox_cluster)
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 1
#> Model Type: coxph
#> Covariates: age, sex, stage
#> Cluster: site
#> Display: adjusted
#> 
#>                      Outcome                Predictor      n Events     aHR (95% CI) p-value
#>                       <char>                   <char> <char> <char>           <char>  <char>
#> 1: Overall Survival (months) Treatment Group (Drug A)    606    606 0.56 (0.45-0.70) < 0.001
#> 2: Overall Survival (months) Treatment Group (Drug B)    606    606 0.83 (0.67-1.01)   0.167

# Example 8: Interaction between predictor and covariate
# Test if treatment effect differs by sex
result_int <- multifit(
    data = clintrial,
    outcomes = c("surgery", "os_status"),
    predictor = "treatment",
    covariates = c("age", "sex", "stage"),
    interactions = c("treatment:sex"),
    labels = clintrial_labels,
    parallel = FALSE
)
print(result_int)
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 2
#> Model Type: glm
#> Covariates: age, sex, stage
#> Interactions: treatment:sex
#> Display: adjusted
#> 
#>               Outcome                             Predictor      n Events     aOR (95% CI) p-value
#>                <char>                                <char> <char> <char>           <char>  <char>
#> 1: Surgical Resection              Treatment Group (Drug A)    292    173 2.09 (1.20-3.62)   0.009
#> 2: Surgical Resection              Treatment Group (Drug B)    361    102 0.50 (0.29-0.87)   0.013
#> 3: Surgical Resection Treatment Group (Drug A) × Sex (Male)      -      - 0.77 (0.34-1.71)   0.515
#> 4: Surgical Resection Treatment Group (Drug B) × Sex (Male)      -      - 0.79 (0.36-1.73)   0.554
#> 5:        Death Event              Treatment Group (Drug A)    292    184 0.39 (0.22-0.70)   0.002
#> 6:        Death Event              Treatment Group (Drug B)    361    273 0.61 (0.34-1.10)   0.099
#> 7:        Death Event Treatment Group (Drug A) × Sex (Male)      -      - 1.29 (0.53-3.18)   0.575
#> 8:        Death Event Treatment Group (Drug B) × Sex (Male)      -      - 1.59 (0.65-3.91)   0.308
# Shows main effects and interaction terms with × notation

# Example 9: Linear model for continuous outcomes
linear_result <- multifit(
    data = clintrial,
    outcomes = c("los_days", "biomarker_x"),
    predictor = "treatment",
    covariates = c("age", "sex"),
    model_type = "lm",
    labels = clintrial_labels,
    parallel = FALSE
)
print(linear_result)
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 2
#> Model Type: lm
#> Covariates: age, sex
#> Display: adjusted
#> 
#>                           Outcome                Predictor      n Adj. Coefficient (95% CI) p-value
#>                            <char>                   <char> <char>                    <char>  <char>
#> 1: Length of Hospital Stay (days) Treatment Group (Drug A)    288    -0.80 (-1.58 to -0.02)   0.045
#> 2: Length of Hospital Stay (days) Treatment Group (Drug B)    350       2.66 (1.91 to 3.41) < 0.001
#> 3:            Biomarker X (ng/mL) Treatment Group (Drug A)    290     -0.39 (-0.94 to 0.16)   0.164
#> 4:            Biomarker X (ng/mL) Treatment Group (Drug B)    358     -0.01 (-0.53 to 0.52)   0.983
# Returns coefficient estimates, not ratios

# Example 10: Poisson regression for equidispersed count outcomes
# fu_count has variance ~= mean, appropriate for standard Poisson
poisson_result <- multifit(
    data = clintrial,
    outcomes = c("fu_count"),
    predictor = "treatment",
    covariates = c("age", "stage"),
    model_type = "glm",
    family = "poisson",
    labels = clintrial_labels,
    parallel = FALSE
)
print(poisson_result)
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 1
#> Model Type: glm
#> Covariates: age, stage
#> Display: adjusted
#> 
#>                  Outcome                Predictor      n Events     aRR (95% CI) p-value
#>                   <char>                   <char> <char> <char>           <char>  <char>
#> 1: Follow-Up Visit Count Treatment Group (Drug A)    290  1,910 1.08 (1.00-1.16)   0.052
#> 2: Follow-Up Visit Count Treatment Group (Drug B)    358  2,438 1.11 (1.03-1.19)   0.005
# Returns rate ratios (RR)
# For overdispersed counts (ae_count), use model_type = "negbin" instead

# Example 11: Filter to significant results only
sig_results <- multifit(
    data = clintrial,
    outcomes = c("surgery", "pfs_status", "os_status"),
    predictor = "stage",
    p_threshold = 0.05,
    labels = clintrial_labels,
    parallel = FALSE
)
print(sig_results)
#> 
#> Multivariate Analysis Results
#> Predictor: stage
#> Outcomes: 3
#> Model Type: glm
#> Display: unadjusted
#> 
#>                       Outcome           Predictor      n Events        OR (95% CI) p-value
#>                        <char>              <char> <char> <char>             <char>  <char>
#> 1:         Surgical Resection Disease Stage (III)    241     91   0.42 (0.29-0.61) < 0.001
#> 2:         Surgical Resection  Disease Stage (IV)    132     19   0.12 (0.07-0.20) < 0.001
#> 3: Progression or Death Event  Disease Stage (II)    263    223   2.37 (1.52-3.71) < 0.001
#> 4: Progression or Death Event Disease Stage (III)    241    222   4.97 (2.86-8.65) < 0.001
#> 5: Progression or Death Event  Disease Stage (IV)    132    128 13.62 (4.82-38.46) < 0.001
#> 6:                Death Event Disease Stage (III)    241    186   2.24 (1.49-3.36) < 0.001
#> 7:                Death Event  Disease Stage (IV)    132    121  7.28 (3.70-14.30) < 0.001
# Only outcomes with significant associations shown

# Example 12: Custom outcome labels
result_labeled <- multifit(
    data = clintrial,
    outcomes = c("surgery", "pfs_status", "os_status"),
    predictor = "treatment",
    labels = c(
        surgery = "Surgical Resection",
        pfs_status = "Disease Progression",
        os_status = "Death",
        treatment = "Treatment Group"
    ),
    parallel = FALSE
)
print(result_labeled)
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 3
#> Model Type: glm
#> Display: unadjusted
#> 
#>                Outcome                Predictor      n Events      OR (95% CI) p-value
#>                 <char>                   <char> <char> <char>           <char>  <char>
#> 1:  Surgical Resection Treatment Group (Drug A)    292    173 1.58 (1.10-2.27)   0.014
#> 2:  Surgical Resection Treatment Group (Drug B)    362    103 0.43 (0.30-0.62) < 0.001
#> 3: Disease Progression Treatment Group (Drug A)    292    227 0.44 (0.26-0.74)   0.002
#> 4: Disease Progression Treatment Group (Drug B)    362    322 1.02 (0.59-1.77)   0.950
#> 5:               Death Treatment Group (Drug A)    292    184 0.51 (0.34-0.76)   0.001
#> 6:               Death Treatment Group (Drug B)    362    274 0.93 (0.62-1.40)   0.721

# Example 13: Keep models for diagnostics
result_models <- multifit(
    data = clintrial,
    outcomes = c("surgery", "os_status"),
    predictor = "treatment",
    covariates = c("age", "sex"),
    keep_models = TRUE,
    parallel = FALSE
)

# Access stored models
models <- attr(result_models, "models")
names(models)
#> [1] "surgery"   "os_status"

# Get adjusted model for surgery outcome
surgery_model <- models$surgery$adjusted
summary(surgery_model)
#> 
#> Call:
#> stats::glm(formula = adj_formula, family = family, data = .data, 
#>     model = keep_models, x = FALSE, y = TRUE)
#> 
#> Coefficients:
#>                 Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)      2.20552    0.41603   5.301 1.15e-07 ***
#> treatmentDrug A  0.49773    0.19099   2.606  0.00916 ** 
#> treatmentDrug B -0.83000    0.18854  -4.402 1.07e-05 ***
#> age             -0.03765    0.00647  -5.820 5.89e-09 ***
#> sexMale         -0.12698    0.14757  -0.861  0.38950    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 1164.1  on 849  degrees of freedom
#> Residual deviance: 1061.7  on 845  degrees of freedom
#> AIC: 1071.7
#> 
#> Number of Fisher Scoring iterations: 4
#> 

# Example 14: Access raw numeric data
result <- multifit(
    data = clintrial,
    outcomes = c("surgery", "os_status"),
    predictor = "age",
    parallel = FALSE
)

# Get unformatted results for custom analysis
raw_data <- attr(result, "raw_data")
print(raw_data)
#>      outcome predictor  group     n events coefficient          se exp_coef  ci_lower  ci_upper statistic      p_value effect_type adjusted
#>       <char>    <char> <char> <int>  <num>       <num>       <num>    <num>     <num>     <num>     <num>        <num>      <char>   <lgcl>
#> 1:   surgery       age      -   850    370 -0.03608738 0.006211760 0.964556 0.9528839 0.9763711 -5.809525 6.265023e-09          OR    FALSE
#> 2: os_status       age      -   850    609  0.04695122 0.006987363 1.048071 1.0338154 1.0625229  6.719447 1.824154e-11          OR    FALSE
# Contains exp_coef, ci_lower, ci_upper, p_value, \emph{etc.}

# Example 15: Hide sample size and event columns
result_minimal <- multifit(
    data = clintrial,
    outcomes = c("surgery", "os_status"),
    predictor = "treatment",
    show_n = FALSE,
    show_events = FALSE,
    parallel = FALSE
)
print(result_minimal)
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 2
#> Model Type: glm
#> Display: unadjusted
#> 
#>      Outcome          Predictor      OR (95% CI) p-value
#>       <char>             <char>           <char>  <char>
#> 1:   surgery treatment (Drug A) 1.58 (1.10-2.27)   0.014
#> 2:   surgery treatment (Drug B) 0.43 (0.30-0.62) < 0.001
#> 3: os_status treatment (Drug A) 0.51 (0.34-0.76)   0.001
#> 4: os_status treatment (Drug B) 0.93 (0.62-1.40)   0.721

# Example 16: Customize decimal places
result_digits <- multifit(
    data = clintrial,
    outcomes = c("surgery", "os_status"),
    predictor = "age",
    digits = 3,
    p_digits = 4,
    parallel = FALSE
)
print(result_digits)
#> 
#> Multivariate Analysis Results
#> Predictor: age
#> Outcomes: 2
#> Model Type: glm
#> Display: unadjusted
#> 
#>      Outcome      n Events         OR (95% CI)  p-value
#>       <char> <char> <char>              <char>   <char>
#> 1:   surgery    850    370 0.965 (0.953-0.976) < 0.0001
#> 2: os_status    850    609 1.048 (1.034-1.063) < 0.0001

# Example 17: Force coefficient display (no exponentiation)
result_coef <- multifit(
    data = clintrial,
    outcomes = c("surgery"),
    predictor = "age",
    exponentiate = FALSE,
    parallel = FALSE
)
print(result_coef)
#> 
#> Multivariate Analysis Results
#> Predictor: age
#> Outcomes: 1
#> Model Type: glm
#> Display: unadjusted
#> 
#>    Outcome      n Events      OR (95% CI) p-value
#>     <char> <char> <char>           <char>  <char>
#> 1: surgery    850    370 0.96 (0.95-0.98) < 0.001

# Example 18: Complete publication workflow
final_table <- multifit(
    data = clintrial,
    outcomes = c("surgery", "pfs_status", "os_status"),
    predictor = "treatment",
    covariates = c("age", "sex", "stage", "grade"),
    columns = "both",
    labels = clintrial_labels,
    digits = 2,
    p_digits = 3,
    parallel = FALSE
)
print(final_table)
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 3
#> Model Type: glm
#> Covariates: age, sex, stage, grade
#> Display: both
#> 
#>                       Outcome                Predictor      n Events      OR (95% CI)   Uni p     aOR (95% CI) Multi p
#>                        <char>                   <char> <char> <char>           <char>  <char>           <char>  <char>
#> 1:                Death Event Treatment Group (Drug A)    292    184 0.51 (0.34-0.76)   0.001 0.41 (0.26-0.65) < 0.001
#> 2:                Death Event Treatment Group (Drug B)    362    274 0.93 (0.62-1.40)   0.721 0.69 (0.44-1.08)   0.105
#> 3: Progression or Death Event Treatment Group (Drug A)    292    227 0.44 (0.26-0.74)   0.002 0.34 (0.20-0.60) < 0.001
#> 4: Progression or Death Event Treatment Group (Drug B)    362    322 1.02 (0.59-1.77)   0.950 0.73 (0.40-1.31)   0.291
#> 5:         Surgical Resection Treatment Group (Drug A)    292    173 1.58 (1.10-2.27)   0.014 1.78 (1.19-2.67)   0.005
#> 6:         Surgical Resection Treatment Group (Drug B)    362    103 0.43 (0.30-0.62) < 0.001 0.43 (0.29-0.64) < 0.001

# Example 19: Gamma regression for positive continuous outcomes
gamma_result <- multifit(
    data = clintrial,
    outcomes = c("los_days", "recovery_days"),
    predictor = "treatment",
    covariates = c("age", "surgery"),
    model_type = "glm",
    family = Gamma(link = "log"),
    labels = clintrial_labels,
    parallel = FALSE
)
print(gamma_result)
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 2
#> Model Type: glm
#> Covariates: age, surgery
#> Display: adjusted
#> 
#>                           Outcome                Predictor      n Adj. Coefficient (95% CI) p-value
#>                            <char>                   <char> <char>                    <char>  <char>
#> 1: Length of Hospital Stay (days) Treatment Group (Drug A)    288    -0.06 (-0.10 to -0.02)       -
#> 2: Length of Hospital Stay (days) Treatment Group (Drug B)    350       0.15 (0.11 to 0.19)       -
#> 3:    Days to Functional Recovery Treatment Group (Drug A)    288    -0.10 (-0.17 to -0.03)       -
#> 4:    Days to Functional Recovery Treatment Group (Drug B)    352       0.17 (0.10 to 0.23)       -
# Returns multiplicative effects on positive continuous data

# Example 20: Quasipoisson for overdispersed counts
quasi_result <- multifit(
    data = clintrial,
    outcomes = c("ae_count"),
    predictor = "treatment",
    covariates = c("age", "diabetes"),
    model_type = "glm",
    family = "quasipoisson",
    labels = clintrial_labels,
    parallel = FALSE
)
print(quasi_result)
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 1
#> Model Type: glm
#> Covariates: age, diabetes
#> Display: adjusted
#> 
#>                Outcome                Predictor      n Events     aRR (95% CI) p-value
#>                 <char>                   <char> <char> <char>           <char>  <char>
#> 1: Adverse Event Count Treatment Group (Drug A)    287  1,221 0.95 (0.81-1.12)       -
#> 2: Adverse Event Count Treatment Group (Drug B)    348  2,459 1.55 (1.34-1.80)       -
# Adjusts standard errors for overdispersion

# Example 21: Generalized linear mixed effects (GLMER)
# Test treatment across outcomes with site clustering
if (requireNamespace("lme4", quietly = TRUE)) {
    glmer_result <- suppressWarnings(multifit(
        data = clintrial,
        outcomes = c("surgery", "pfs_status", "os_status"),
        predictor = "treatment",
        covariates = c("age", "sex"),
        random = "(1|site)",
        model_type = "glmer",
        family = "binomial",
        labels = clintrial_labels,
        parallel = FALSE
    ))
    print(glmer_result)
}
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 3
#> Model Type: glmer
#> Covariates: age, sex
#> Random Effects: (1|site)
#> Display: adjusted
#> 
#>                       Outcome                Predictor      n Events     aOR (95% CI) p-value
#>                        <char>                   <char> <char> <char>           <char>  <char>
#> 1:         Surgical Resection Treatment Group (Drug A)    292    173 1.64 (1.13-2.40)   0.010
#> 2:         Surgical Resection Treatment Group (Drug B)    362    103 0.44 (0.30-0.63) < 0.001
#> 3: Progression or Death Event Treatment Group (Drug A)    292    227 0.41 (0.24-0.71)   0.001
#> 4: Progression or Death Event Treatment Group (Drug B)    362    322 1.01 (0.58-1.78)   0.971
#> 5:                Death Event Treatment Group (Drug A)    292    184 0.44 (0.28-0.69) < 0.001
#> 6:                Death Event Treatment Group (Drug B)    362    274 0.88 (0.56-1.38)   0.571

# Example 22: Cox mixed effects with random site effects
if (requireNamespace("coxme", quietly = TRUE)) {
    coxme_result <- multifit(
        data = clintrial,
        outcomes = c("Surv(pfs_months, pfs_status)", 
                     "Surv(os_months, os_status)"),
        predictor = "treatment",
        covariates = c("age", "sex", "stage"),
        random = "(1|site)",
        model_type = "coxme",
        labels = clintrial_labels,
        parallel = FALSE
    )
    print(coxme_result)
}
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 2
#> Model Type: coxme
#> Covariates: age, sex, stage
#> Random Effects: (1|site)
#> Display: adjusted
#> 
#>                               Outcome                Predictor      n Events     aHR (95% CI) p-value
#>                                <char>                   <char> <char> <char>           <char>  <char>
#> 1: Progression-Free Survival (months) Treatment Group (Drug A)    721    721 0.53 (0.44-0.66) < 0.001
#> 2: Progression-Free Survival (months) Treatment Group (Drug B)    721    721 0.85 (0.71-1.03)   0.099
#> 3:          Overall Survival (months) Treatment Group (Drug A)    606    606 0.57 (0.46-0.72) < 0.001
#> 4:          Overall Survival (months) Treatment Group (Drug B)    606    606 0.93 (0.76-1.14)   0.492

# Example 23: Multiple interactions across outcomes
multi_int <- multifit(
    data = clintrial,
    outcomes = c("surgery", "pfs_status", "os_status"),
    predictor = "treatment",
    covariates = c("age", "sex", "stage"),
    interactions = c("treatment:stage", "treatment:sex"),
    labels = clintrial_labels,
    parallel = FALSE
)
print(multi_int)
#> 
#> Multivariate Analysis Results
#> Predictor: treatment
#> Outcomes: 3
#> Model Type: glm
#> Covariates: age, sex, stage
#> Interactions: treatment:stage, treatment:sex
#> Display: adjusted
#> 
#>                        Outcome                                      Predictor      n Events          aOR (95% CI) p-value
#>                         <char>                                         <char> <char> <char>                <char>  <char>
#>  1:         Surgical Resection                       Treatment Group (Drug A)    292    173      2.29 (0.99-5.25)   0.052
#>  2:         Surgical Resection                       Treatment Group (Drug B)    361    102      0.49 (0.22-1.11)   0.087
#>  3:         Surgical Resection  Treatment Group (Drug A) × Disease Stage (II)      -      -      0.99 (0.36-2.72)   0.978
#>  4:         Surgical Resection  Treatment Group (Drug B) × Disease Stage (II)      -      -      0.99 (0.37-2.61)   0.977
#>  5:         Surgical Resection Treatment Group (Drug A) × Disease Stage (III)      -      -      0.98 (0.32-2.96)   0.973
#>  6:         Surgical Resection Treatment Group (Drug B) × Disease Stage (III)      -      -      0.98 (0.34-2.82)   0.969
#>  7:         Surgical Resection  Treatment Group (Drug A) × Disease Stage (IV)      -      -      0.55 (0.13-2.35)   0.417
#>  8:         Surgical Resection  Treatment Group (Drug B) × Disease Stage (IV)      -      -      1.53 (0.35-6.75)   0.577
#>  9:         Surgical Resection          Treatment Group (Drug A) × Sex (Male)      -      -      0.75 (0.33-1.69)   0.487
#> 10:         Surgical Resection          Treatment Group (Drug B) × Sex (Male)      -      -      0.81 (0.37-1.76)   0.587
#> 11: Progression or Death Event                       Treatment Group (Drug A)    292    227      0.41 (0.16-1.06)   0.066
#> 12: Progression or Death Event                       Treatment Group (Drug B)    361    321      0.45 (0.17-1.25)   0.127
#> 13: Progression or Death Event  Treatment Group (Drug A) × Disease Stage (II)      -      -      1.55 (0.47-5.14)   0.472
#> 14: Progression or Death Event  Treatment Group (Drug B) × Disease Stage (II)      -      -     3.97 (1.10-14.31)   0.035
#> 15: Progression or Death Event Treatment Group (Drug A) × Disease Stage (III)      -      -      0.29 (0.03-2.74)   0.281
#> 16: Progression or Death Event Treatment Group (Drug B) × Disease Stage (III)      -      -     1.47 (0.14-15.30)   0.745
#> 17: Progression or Death Event  Treatment Group (Drug A) × Disease Stage (IV)      -      -     1.37 (0.12-16.05)   0.804
#> 18: Progression or Death Event  Treatment Group (Drug B) × Disease Stage (IV)      -      - 3691058.18 (0.00-Inf)   0.977
#> 19: Progression or Death Event          Treatment Group (Drug A) × Sex (Male)      -      -      0.74 (0.24-2.27)   0.602
#> 20: Progression or Death Event          Treatment Group (Drug B) × Sex (Male)      -      -      0.68 (0.21-2.25)   0.529
#> 21:                Death Event                       Treatment Group (Drug A)    292    184      0.65 (0.29-1.48)   0.303
#> 22:                Death Event                       Treatment Group (Drug B)    361    273      0.45 (0.19-1.07)   0.071
#> 23:                Death Event  Treatment Group (Drug A) × Disease Stage (II)      -      -      0.68 (0.24-1.90)   0.460
#> 24:                Death Event  Treatment Group (Drug B) × Disease Stage (II)      -      -      1.43 (0.50-4.05)   0.504
#> 25:                Death Event Treatment Group (Drug A) × Disease Stage (III)      -      -      0.35 (0.10-1.23)   0.100
#> 26:                Death Event Treatment Group (Drug B) × Disease Stage (III)      -      -      1.08 (0.30-3.81)   0.909
#> 27:                Death Event  Treatment Group (Drug A) × Disease Stage (IV)      -      -      0.21 (0.02-1.98)   0.171
#> 28:                Death Event  Treatment Group (Drug B) × Disease Stage (IV)      -      -     3.32 (0.18-61.53)   0.421
#> 29:                Death Event          Treatment Group (Drug A) × Sex (Male)      -      -      1.25 (0.51-3.08)   0.624
#> 30:                Death Event          Treatment Group (Drug B) × Sex (Male)      -      -      1.78 (0.71-4.45)   0.218
#>                        Outcome                                      Predictor      n Events          aOR (95% CI) p-value
#>                         <char>                                         <char> <char> <char>                <char>  <char>
# Shows how treatment effects vary by stage and sex across outcomes

# }
```
