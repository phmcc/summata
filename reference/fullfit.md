# Complete Regression Analysis Workflow

Executes a comprehensive regression analysis pipeline that combines
univariable screening, automatic/manual variable selection, and
multivariable modeling in a single function call. This function is
designed to streamline the complete analytical workflow from initial
exploration to final adjusted models, with publication-ready formatted
output showing both univariable and multivariable results side-by-side
if desired.

## Usage

``` r
fullfit(
  data,
  outcome,
  predictors,
  method = "screen",
  multi_predictors = NULL,
  p_threshold = 0.05,
  columns = "both",
  model_type = "glm",
  family = "binomial",
  random = NULL,
  conf_level = 0.95,
  reference_rows = TRUE,
  show_n = TRUE,
  show_events = TRUE,
  digits = 2,
  p_digits = 3,
  labels = NULL,
  metrics = "both",
  return_type = "table",
  keep_models = FALSE,
  exponentiate = NULL,
  conf_method = NULL,
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

- outcome:

  Character string specifying the outcome variable name. For
  time-to-event analysis, use
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) syntax for the
  outcome variable (*e.g.,* `"Surv(os_months, os_status)"`).

- predictors:

  Character vector of predictor variable names to analyze. All
  predictors are tested in univariable models. The subset included in
  the multivariable model depends on the `method` parameter.

- method:

  Character string specifying the variable selection strategy:

  - `"screen"` - Automatic selection based on univariable *p*-value
    threshold. Only predictors with p \\\le\\ `p_threshold` in
    univariable analysis are included in the multivariable model
    \[default\]

  - `"all"` - Include all predictors in both univariable and
    multivariable analyses (no selection)

  - `"custom"` - Manual selection. All predictors in univariable
    analysis, but only those specified in `multi_predictors` are
    included in multivariable model

- multi_predictors:

  Character vector of predictors to include in the multivariable model
  when `method = "custom"`. Required when using custom selection.
  Ignored for other methods. Default is `NULL`.

- p_threshold:

  Numeric *p*-value threshold for automatic variable selection when
  `method = "screen"`. Predictors with univariable *p*-value less than
  or equal to the threshold are included in multivariable model. Common
  values: 0.05 (strict), 0.10 (moderate), 0.20 (liberal screening).
  Default is 0.05. Ignored for other methods.

- columns:

  Character string specifying which result columns to display:

  - `"both"` - Show both univariable and multivariable results
    side-by-side \[default\]

  - `"uni"` - Show only univariable results

  - `"multi"` - Show only multivariable results

- model_type:

  Character string specifying the regression model type:

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
    (*e.g.,* costs, lengths of stay). When passed as a string, resolves
    to log link for interpretable multiplicative effects.

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

- random:

  Character string specifying the random-effects formula for
  mixed-effects models (`glmer`, `lmer`, `coxme`). Use standard
  lme4/coxme syntax, *e.g.,* `"(1|site)"` for random intercepts by site,
  `"(1|site/patient)"` for nested random effects. Required when
  `model_type` is a mixed-effects model type unless random effects are
  included in the `predictors` vector. Alternatively, random effects can
  be included directly in the `predictors` vector using the same syntax
  (*e.g.,* `predictors = c("age", "sex", "(1|site)")`), though they will
  not be screened as predictors. Default is `NULL`.

- conf_level:

  Numeric confidence level for confidence intervals. Must be between 0
  and 1. Default is 0.95 (95% CI).

- reference_rows:

  Logical. If `TRUE`, adds rows for reference categories of factor
  variables with baseline values. Default is `TRUE`.

- show_n:

  Logical. If `TRUE`, includes sample size columns. Default is `TRUE`.

- show_events:

  Logical. If `TRUE`, includes events columns (survival and logistic
  models). Default is `TRUE`.

- digits:

  Integer specifying decimal places for effect estimates. Default is 2.

- p_digits:

  Integer specifying the number of decimal places for *p*-values. Values
  smaller than `10^(-p_digits)` are displayed as `"< 0.001"` (for
  `p_digits = 3`), `"< 0.0001"` (for `p_digits = 4`), etc. Default is 3.

- labels:

  Named character vector or list providing custom display labels for
  variables. Names should match variable names, values are display
  labels. Default is `NULL`.

- metrics:

  Character specification for which statistics to display:

  - `"both"` - Show effect estimates with CI and *p*-values \[default\]

  - `"effect"` - Show only effect estimates with CI

  - `"p"` - Show only *p*-values

  Can also be a character vector: `c("effect", "p")` is equivalent to
  `"both"`.

- return_type:

  Character string specifying what to return:

  - `"table"` - Return formatted results table only \[default\]

  - `"model"` - Return multivariable model object only

  - `"both"` - Return list with both table and model

- keep_models:

  Logical. If `TRUE`, stores univariable model objects in the output.
  Can consume significant memory for many predictors. Default is
  `FALSE`.

- exponentiate:

  Logical. Whether to exponentiate coefficients. Default is `NULL`,
  which automatically exponentiates for logistic, Poisson, and Cox
  models, and displays raw coefficients for linear models.

- conf_method:

  Character string controlling the confidence interval method. If `NULL`
  (default), uses `getOption("summata.conf_method", "profile")`.

  - `"profile"` - Profile likelihood intervals for GLM and negative
    binomial models (via
    [`MASS::confint.glm()`](https://rdrr.io/pkg/MASS/man/confint.html)),
    exact *t*-distribution intervals for linear models. Falls back to
    Wald on profiling failure. Quasi-likelihood families always use Wald
    (no true likelihood).

  - `"wald"` - Wald intervals (coefficient \\\pm\\ *z* \\\times\\ SE)
    for all model types. Faster but less accurate near boundary
    conditions or with small subgroups.

  Cox and mixed-effects models use Wald intervals regardless of this
  setting. Set globally with `options(summata.conf_method = "wald")` to
  use Wald throughout a session. Note: when `method = "screen"` and
  `columns = "multi"`, the internal screening pass always uses Wald
  since only *p*-values are needed for variable selection.

- parallel:

  Logical. If `TRUE` (default), fits univariable models in parallel
  using multiple CPU cores for improved performance.

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

  Additional arguments passed to model fitting functions (*e.g.,*
  `weights`, `subset`, `na.action`).

## Value

Depends on `return_type` parameter:

When `return_type = "table"` (default): A data.table with S3 class
`"fullfit_result"` containing:

- Variable:

  Character. Predictor name or custom label

- Group:

  Character. Category level for factors, empty for continuous

- n/n_group:

  Integer. Sample sizes (if `show_n = TRUE`). For variables included in
  the multivariable model, reflects the complete-case sample size from
  the fitted model (listwise deletion across all included predictors).
  For variables not selected into the multivariable model, reflects the
  per-variable sample size from the univariable analysis. This follows
  STROBE guideline item 12, which recommends reporting the number of
  participants included at each stage of analysis.

- events/events_group:

  Integer. Event counts (if `show_events = TRUE`). Same complete-case
  convention as `n`: multivariable rows show events from the fitted
  model, univariable-only rows show per-variable counts.

- OR/HR/RR/Coefficient (95% CI):

  Character. Unadjusted effect (if `columns` includes "uni" and
  `metrics` includes "effect")

- Uni p:

  Character. Univariable *p*-value (if `columns` includes "uni" and
  `metrics` includes "p")

- aOR/aHR/aRR/Adj. Coefficient (95% CI):

  Character. Adjusted effect (if `columns` includes "multi" and
  `metrics` includes "effect")

- Multi p:

  Character. Multivariable *p*-value (if `columns` includes "multi" and
  `metrics` includes "p")

When `return_type = "model"`: The fitted multivariable model object
(glm, lm, coxph, *etc.*).

When `return_type = "both"`: A list with two elements:

- table:

  The formatted results data.table

- model:

  The fitted multivariable model object

The table includes the following attributes:

- outcome:

  Character. The outcome variable name

- model_type:

  Character. The regression model type

- method:

  Character. The variable selection method used

- columns:

  Character. Which columns were displayed

- model:

  The multivariable model object (if fitted)

- uni_results:

  The complete univariable screening results

- n_multi:

  Integer. Number of predictors in multivariable model

- screened:

  Character vector. Names of predictors that passed univariable
  screening at the specified *p*-value threshold

- significant:

  Character vector. Names of variables with p \< 0.05 in the
  multivariable model (or univariable if multivariable was not fitted)

## Details

**Analysis Workflow:**

The function implements a complete regression analysis pipeline:

1.  **Univariable screening**: Fits separate models for each predictor
    (outcome ~ predictor). Each predictor is tested independently to
    understand crude associations.

2.  **Variable selection**: Based on the `method` parameter:

    - `"screen"`: Automatically selects predictors with univariable p
      \\\le\\ `p_threshold`

    - `"all"`: Includes all predictors (no selection)

    - `"custom"`: Uses predictors specified in `multi_predictors`

3.  **Multivariable modeling**: Fits a single model with selected
    predictors (outcome ~ predictor1 + predictor2 + ...). Estimates are
    adjusted for all other variables in the model.

4.  **Output formatting**: Combines results into publication-ready table
    with appropriate effect measures and formatting.

**Variable Selection Strategies:**

*"Screen" Method* (`method = "screen"`):

- Uses *p*-value threshold for automatic selection

- Liberal thresholds (*e.g.,* 0.20) cast a wide net to avoid missing
  important predictors

- Stricter thresholds (*e.g.,* 0.05) focus on strongly associated
  predictors

- Helps reduce overfitting and multicollinearity

- Common in exploratory analyses and when sample size is limited

*"All" Method* (`method = "all"`):

- No variable selection - includes all predictors

- Appropriate when all variables are theoretically important

- Risk of overfitting with many predictors relative to sample size

- Useful for confirmatory analyses with pre-specified models

*"Custom" Method* (`method = "custom"`):

- Manual selection based on subject matter knowledge

- Runs univariable analysis for all predictors (for comparison)

- Includes only specified predictors in multivariable model

- Ideal for theory-driven model building

- Allows comparison of unadjusted vs adjusted effects for all variables

**Interpreting Results:**

When `columns = "both"` (default), tables show:

- **Univariable columns**: Crude associations, unadjusted for other
  variables. Labeled as "OR/HR/RR/Coefficient (95% CI)" and "Uni *p*"

- **Multivariable columns**: Adjusted associations, accounting for all
  other predictors in the model. Labeled as "aOR/aHR/aRR/Adj.
  Coefficient (95% CI)" and "Multi *p*" ("a" = adjusted)

- Variables not meeting selection criteria show "-" in multivariable
  columns

Comparing univariable and multivariable results helps identify:

- **Confounding**: Large changes in effect estimates

- **Independent effects**: Similar univariable and multivariable
  estimates

- **Mediation**: Attenuated effects in multivariable model

- **Suppression**: Effects that emerge only after adjustment

**Sample Size Considerations:**

Rule of thumb for multivariable models:

- **Logistic regression**: \\\ge\\ 10 events per predictor variable

- **Cox regression**: \\\ge\\ 10 events per predictor variable

- **Linear regression**: \\\ge\\ 10-20 observations per predictor

Use screening methods to reduce predictor count when these ratios are
not met.

## See also

[`uniscreen`](https://phmcc.github.io/summata/reference/uniscreen.md)
for univariable screening only,
[`fit`](https://phmcc.github.io/summata/reference/fit.md) for fitting a
single multivariable model,
[`compfit`](https://phmcc.github.io/summata/reference/compfit.md) for
comparing multiple models,
[`desctable`](https://phmcc.github.io/summata/reference/desctable.md)
for descriptive statistics

Other regression functions:
[`compfit()`](https://phmcc.github.io/summata/reference/compfit.md),
[`fit()`](https://phmcc.github.io/summata/reference/fit.md),
[`multifit()`](https://phmcc.github.io/summata/reference/multifit.md),
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

# Example 1: Basic screening with p < 0.05 threshold
result1 <- fullfit(
    data = clintrial,
    outcome = "os_status",
    predictors = c("age", "sex", "bmi", "smoking",
                   "hypertension", "diabetes",
                   "treatment", "stage"),
    method = "screen",
    p_threshold = 0.05,
    labels = clintrial_labels
)
#> Running univariable analysis...
#> Fitting multivariable model with 5 predictors...
print(result1)
#> 
#> Fullfit Analysis Results
#> Outcome: os_status
#> Model Type: glm
#> Method: screen (p < 0.05)
#> Predictors Screened: 8
#> Multivariable Predictors: 5
#> 
#>                    Variable   Group      n Events       OR (95% CI)   Uni p      aOR (95% CI) Multi p
#>                      <char>  <char> <char> <char>            <char>  <char>            <char>  <char>
#>  1:             Age (years)       -    833    594  1.05 (1.03-1.06) < 0.001  1.05 (1.04-1.07) < 0.001
#>  2:                     Sex  Female    443    292         reference       -         reference       -
#>  3:                            Male    390    302  1.78 (1.31-2.43) < 0.001  1.98 (1.41-2.78) < 0.001
#>  4: Body Mass Index (kg/m²)       -    838    599  1.01 (0.98-1.05)   0.347                 -       -
#>  5:          Smoking Status   Never    337    248         reference       -         reference       -
#>  6:                          Former    311    203  0.67 (0.48-0.94)   0.022  0.72 (0.50-1.05)   0.087
#>  7:                         Current    185    143  1.22 (0.81-1.87)   0.351  1.38 (0.88-2.19)   0.162
#>  8:            Hypertension      No    504    354         reference       -                 -       -
#>  9:                             Yes    331    242  1.15 (0.85-1.57)   0.369                 -       -
#> 10:                Diabetes      No    637    457         reference       -                 -       -
#> 11:                             Yes    197    138  0.92 (0.65-1.31)   0.646                 -       -
#> 12:         Treatment Group Control    191    146         reference       -         reference       -
#> 13:                          Drug A    288    181  0.51 (0.34-0.76)   0.001  0.44 (0.28-0.68) < 0.001
#> 14:                          Drug B    354    267  0.93 (0.61-1.39)   0.721  0.77 (0.49-1.20)   0.251
#> 15:           Disease Stage       I    207    125         reference       -         reference       -
#> 16:                              II    261    170  1.25 (0.86-1.82)   0.243  1.26 (0.84-1.90)   0.263
#> 17:                             III    237    182  2.24 (1.49-3.38) < 0.001  2.54 (1.64-3.99) < 0.001
#> 18:                              IV    128    117 7.28 (3.85-15.04) < 0.001 8.61 (4.40-18.33) < 0.001
# Shows both univariable and multivariable results
# Only significant univariable predictors in multivariable model

# \donttest{

# Example 2: Include all predictors (no selection)
result2 <- fullfit(
    data = clintrial,
    outcome = "os_status",
    predictors = c("age", "sex", "treatment", "stage"),
    method = "all",
    labels = clintrial_labels
)
#> Running univariable analysis...
#> Fitting multivariable model with 4 predictors...
print(result2)
#> 
#> Fullfit Analysis Results
#> Outcome: os_status
#> Model Type: glm
#> Method: all
#> Predictors Screened: 4
#> Multivariable Predictors: 4
#> 
#>            Variable   Group      n Events       OR (95% CI)   Uni p      aOR (95% CI) Multi p
#>              <char>  <char> <char> <char>            <char>  <char>            <char>  <char>
#>  1:     Age (years)       -    847    606  1.05 (1.03-1.06) < 0.001  1.05 (1.04-1.07) < 0.001
#>  2:             Sex  Female    449    297         reference       -         reference       -
#>  3:                    Male    398    309  1.78 (1.31-2.43) < 0.001  2.02 (1.45-2.83) < 0.001
#>  4: Treatment Group Control    194    149         reference       -         reference       -
#>  5:                  Drug A    292    184  0.51 (0.34-0.76)   0.001  0.44 (0.28-0.68) < 0.001
#>  6:                  Drug B    361    273  0.93 (0.61-1.39)   0.721  0.74 (0.47-1.16)   0.192
#>  7:   Disease Stage       I    211    127         reference       -         reference       -
#>  8:                      II    263    172  1.25 (0.86-1.82)   0.243  1.32 (0.88-1.97)   0.181
#>  9:                     III    241    186  2.24 (1.49-3.38) < 0.001  2.70 (1.74-4.21) < 0.001
#> 10:                      IV    132    121 7.28 (3.85-15.04) < 0.001 9.08 (4.66-19.28) < 0.001

# Example 3: Custom variable selection
result3 <- fullfit(
    data = clintrial,
    outcome = "os_status",
    predictors = c("age", "sex", "bmi", "smoking", "treatment", "stage"),
    method = "custom",
    multi_predictors = c("age", "treatment", "stage"),
    labels = clintrial_labels
)
#> Running univariable analysis...
#> Fitting multivariable model with 3 predictors...
print(result3)
#> 
#> Fullfit Analysis Results
#> Outcome: os_status
#> Model Type: glm
#> Method: custom
#> Predictors Screened: 6
#> Multivariable Predictors: 3
#> 
#>                    Variable   Group      n Events       OR (95% CI)   Uni p      aOR (95% CI) Multi p
#>                      <char>  <char> <char> <char>            <char>  <char>            <char>  <char>
#>  1:             Age (years)       -    847    606  1.05 (1.03-1.06) < 0.001  1.05 (1.04-1.07) < 0.001
#>  2:                     Sex  Female    450    298         reference       -                 -       -
#>  3:                            Male    400    311  1.78 (1.31-2.43) < 0.001                 -       -
#>  4: Body Mass Index (kg/m²)       -    838    599  1.01 (0.98-1.05)   0.347                 -       -
#>  5:          Smoking Status   Never    337    248         reference       -                 -       -
#>  6:                          Former    311    203  0.67 (0.48-0.94)   0.022                 -       -
#>  7:                         Current    185    143  1.22 (0.81-1.87)   0.351                 -       -
#>  8:         Treatment Group Control    194    149         reference       -         reference       -
#>  9:                          Drug A    292    184  0.51 (0.34-0.76)   0.001  0.44 (0.28-0.68) < 0.001
#> 10:                          Drug B    361    273  0.93 (0.61-1.39)   0.721  0.78 (0.50-1.20)   0.260
#> 11:           Disease Stage       I    211    127         reference       -         reference       -
#> 12:                              II    263    172  1.25 (0.86-1.82)   0.243  1.31 (0.88-1.95)   0.181
#> 13:                             III    241    186  2.24 (1.49-3.38) < 0.001  2.56 (1.66-3.96) < 0.001
#> 14:                              IV    132    121 7.28 (3.85-15.04) < 0.001 8.43 (4.36-17.76) < 0.001
# Univariable for all, multivariable for selected only

# Example 4: Cox regression with screening
library(survival)
cox_result <- fullfit(
    data = clintrial,
    outcome = "Surv(os_months, os_status)",
    predictors = c("age", "sex", "treatment", "stage"),
    model_type = "coxph",
    method = "screen",
    p_threshold = 0.10,
    labels = clintrial_labels
)
#> Running univariable analysis...
#> Fitting multivariable model with 4 predictors...
print(cox_result)
#> 
#> Fullfit Analysis Results
#> Outcome: Surv(os_months, os_status)
#> Model Type: coxph
#> Method: screen (p < 0.1)
#> Predictors Screened: 4
#> Multivariable Predictors: 4
#> 
#>            Variable   Group      n Events      HR (95% CI)   Uni p     aHR (95% CI) Multi p
#>              <char>  <char> <char> <char>           <char>  <char>           <char>  <char>
#>  1:     Age (years)       -    847    606 1.03 (1.03-1.04) < 0.001 1.04 (1.03-1.04) < 0.001
#>  2:             Sex  Female    450    298        reference       -        reference       -
#>  3:                    Male    400    311 1.30 (1.11-1.53)   0.001 1.33 (1.13-1.56) < 0.001
#>  4: Treatment Group Control    196    151        reference       -        reference       -
#>  5:                  Drug A    292    184 0.64 (0.52-0.80) < 0.001 0.56 (0.45-0.70) < 0.001
#>  6:                  Drug B    362    274 0.94 (0.77-1.15)   0.567 0.83 (0.67-1.01)   0.062
#>  7:   Disease Stage       I    211    127        reference       -        reference       -
#>  8:                      II    263    172 1.12 (0.89-1.41)   0.337 1.16 (0.92-1.46)   0.214
#>  9:                     III    241    186 1.69 (1.35-2.11) < 0.001 1.90 (1.51-2.38) < 0.001
#> 10:                      IV    132    121 3.18 (2.47-4.09) < 0.001 3.57 (2.77-4.60) < 0.001

# Example 5: Linear regression without screening
linear_result <- fullfit(
    data = clintrial,
    outcome = "bmi",
    predictors = c("age", "sex", "smoking", "creatinine"),
    model_type = "lm",
    method = "all",
    labels = clintrial_labels
)
#> Running univariable analysis...
#> Fitting multivariable model with 4 predictors...
print(linear_result)
#> 
#> Fullfit Analysis Results
#> Outcome: bmi
#> Model Type: lm
#> Method: all
#> Predictors Screened: 4
#> Multivariable Predictors: 4
#> 
#>                       Variable   Group      n  Coefficient (95% CI)  Uni p Adj. Coefficient (95% CI) Multi p
#>                         <char>  <char> <char>                <char> <char>                    <char>  <char>
#> 1:                 Age (years)       -    833 -0.02 (-0.05 to 0.01)  0.140     -0.02 (-0.05 to 0.01)   0.138
#> 2:                         Sex  Female    443             reference      -                 reference       -
#> 3:                                Male    390 -0.36 (-1.03 to 0.31)  0.296     -0.34 (-1.01 to 0.34)   0.328
#> 4:              Smoking Status   Never    337             reference      -                 reference       -
#> 5:                              Former    311  0.22 (-0.55 to 0.98)  0.578      0.15 (-0.61 to 0.92)   0.697
#> 6:                             Current    185  0.31 (-0.57 to 1.20)  0.490      0.26 (-0.62 to 1.15)   0.558
#> 7: Baseline Creatinine (mg/dL)       -    833  0.52 (-0.61 to 1.65)  0.367      0.55 (-0.58 to 1.68)   0.340

# Example 6: Poisson regression for count outcomes
poisson_result <- fullfit(
    data = clintrial,
    outcome = "fu_count",
    predictors = c("age", "stage", "treatment", "surgery"),
    model_type = "glm",
    family = "poisson",
    method = "all",
    labels = clintrial_labels
)
#> Running univariable analysis...
#> Fitting multivariable model with 4 predictors...
print(poisson_result)
#> 
#> Fullfit Analysis Results
#> Outcome: fu_count
#> Model Type: glm
#> Method: all
#> Predictors Screened: 4
#> Multivariable Predictors: 4
#> 
#>               Variable   Group      n Events      RR (95% CI)   Uni p     aRR (95% CI) Multi p
#>                 <char>  <char> <char> <char>           <char>  <char>           <char>  <char>
#>  1:        Age (years)       -    839  5,517 1.00 (0.99-1.00)   0.005 1.00 (1.00-1.00)   0.038
#>  2:      Disease Stage       I    207  1,238        reference       -        reference       -
#>  3:                         II    261  1,638 1.05 (0.97-1.13)   0.201 1.05 (0.98-1.13)   0.169
#>  4:                        III    240  1,638 1.14 (1.06-1.23) < 0.001 1.14 (1.06-1.23) < 0.001
#>  5:                         IV    131  1,003 1.28 (1.18-1.39) < 0.001 1.33 (1.21-1.45) < 0.001
#>  6:    Treatment Group Control    191  1,169        reference       -        reference       -
#>  7:                     Drug A    290  1,910 1.08 (1.00-1.16)   0.045 1.06 (0.99-1.14)   0.101
#>  8:                     Drug B    358  2,438 1.11 (1.04-1.19)   0.003 1.12 (1.04-1.20)   0.002
#>  9: Surgical Resection      No    476  3,091        reference       -        reference       -
#> 10:                        Yes    363  2,426 1.03 (0.97-1.08)   0.322 1.09 (1.03-1.16)   0.003

# Example 7: Show only multivariable results
multi_only <- fullfit(
    data = clintrial,
    outcome = "os_status",
    predictors = c("age", "sex", "treatment", "stage"),
    method = "all",
    columns = "multi",
    labels = clintrial_labels
)
#> Fitting multivariable model with 4 predictors...
print(multi_only)
#> 
#> Fullfit Analysis Results
#> Outcome: os_status
#> Model Type: glm
#> Method: all
#> Predictors Screened: 4
#> Multivariable Predictors: 4
#> 
#>            Variable   Group      n Events      aOR (95% CI) p-value
#>              <char>  <char> <char> <char>            <char>  <char>
#>  1:     Age (years)       -    847    606  1.05 (1.04-1.07) < 0.001
#>  2:             Sex  Female    449    297         reference       -
#>  3:                    Male    398    309  2.02 (1.45-2.83) < 0.001
#>  4: Treatment Group Control    194    149         reference       -
#>  5:                  Drug A    292    184  0.44 (0.28-0.68) < 0.001
#>  6:                  Drug B    361    273  0.74 (0.47-1.16)   0.192
#>  7:   Disease Stage       I    211    127         reference       -
#>  8:                      II    263    172  1.32 (0.88-1.97)   0.181
#>  9:                     III    241    186  2.70 (1.74-4.21) < 0.001
#> 10:                      IV    132    121 9.08 (4.66-19.28) < 0.001

# Example 8: Return both table and model object
both <- fullfit(
    data = clintrial,
    outcome = "os_status",
    predictors = c("age", "sex", "treatment", "stage"),
    method = "all",
    return_type = "both"
)
#> Running univariable analysis...
#> Fitting multivariable model with 4 predictors...
print(both$table)
#> 
#> Fullfit Analysis Results
#> Outcome: os_status
#> Model Type: glm
#> Method: all
#> Predictors Screened: 4
#> Multivariable Predictors: 4
#> 
#>      Variable   Group      n Events       OR (95% CI)   Uni p      aOR (95% CI) Multi p
#>        <char>  <char> <char> <char>            <char>  <char>            <char>  <char>
#>  1:       age       -    847    606  1.05 (1.03-1.06) < 0.001  1.05 (1.04-1.07) < 0.001
#>  2:       sex  Female    449    297         reference       -         reference       -
#>  3:              Male    398    309  1.78 (1.31-2.43) < 0.001  2.02 (1.45-2.83) < 0.001
#>  4: treatment Control    194    149         reference       -         reference       -
#>  5:            Drug A    292    184  0.51 (0.34-0.76)   0.001  0.44 (0.28-0.68) < 0.001
#>  6:            Drug B    361    273  0.93 (0.61-1.39)   0.721  0.74 (0.47-1.16)   0.192
#>  7:     stage       I    211    127         reference       -         reference       -
#>  8:                II    263    172  1.25 (0.86-1.82)   0.243  1.32 (0.88-1.97)   0.181
#>  9:               III    241    186  2.24 (1.49-3.38) < 0.001  2.70 (1.74-4.21) < 0.001
#> 10:                IV    132    121 7.28 (3.85-15.04) < 0.001 9.08 (4.66-19.28) < 0.001
summary(both$model)
#> 
#> Call:
#> stats::glm(formula = formula, family = family, data = data)
#> 
#> Coefficients:
#>                  Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)     -2.646316   0.491675  -5.382 7.36e-08 ***
#> age              0.052666   0.007447   7.072 1.53e-12 ***
#> sexMale          0.702430   0.170410   4.122 3.76e-05 ***
#> treatmentDrug A -0.829478   0.227358  -3.648 0.000264 ***
#> treatmentDrug B -0.297335   0.227823  -1.305 0.191854    
#> stageII          0.274497   0.205201   1.338 0.180994    
#> stageIII         0.992064   0.224514   4.419 9.93e-06 ***
#> stageIV          2.205770   0.359692   6.132 8.66e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 1011.63  on 846  degrees of freedom
#> Residual deviance:  869.99  on 839  degrees of freedom
#>   (3 observations deleted due to missingness)
#> AIC: 885.99
#> 
#> Number of Fisher Scoring iterations: 5
#> 

# Example 9: Keep univariable models for diagnostics
with_models <- fullfit(
    data = clintrial,
    outcome = "os_status",
    predictors = c("age", "bmi", "creatinine"),
    keep_models = TRUE
)
#> Running univariable analysis...
#> Fitting multivariable model with 1 predictors...
uni_results <- attr(with_models, "uni_results")
uni_models <- attr(uni_results, "models")
summary(uni_models[["age"]])
#> 
#> Call:
#> stats::glm(formula = formula, family = family, data = data, model = keep_models, 
#>     x = FALSE, y = TRUE)
#> 
#> Coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -1.825668   0.408920  -4.465 8.02e-06 ***
#> age          0.046951   0.006987   6.719 1.82e-11 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 1013.6  on 849  degrees of freedom
#> Residual deviance:  964.4  on 848  degrees of freedom
#> AIC: 968.4
#> 
#> Number of Fisher Scoring iterations: 4
#> 

# Example 10: Linear mixed effects with site clustering
if (requireNamespace("lme4", quietly = TRUE)) {
    lmer_result <- fullfit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "treatment", "surgery", "stage"),
        random = "(1|site)",
        model_type = "lmer",
        method = "all",
        labels = clintrial_labels
    )
    print(lmer_result)
}
#> Running univariable analysis...
#> Fitting multivariable model with 4 predictors...
#> 
#> Fullfit Analysis Results
#> Outcome: los_days
#> Model Type: lmer
#> Method: all
#> Predictors Screened: 4
#> Multivariable Predictors: 4
#> 
#>               Variable   Group      n  Coefficient (95% CI)   Uni p Adj. Coefficient (95% CI) Multi p
#>                 <char>  <char> <char>                <char>  <char>                    <char>  <char>
#>  1:        Age (years)       -    827   0.14 (0.11 to 0.16) < 0.001       0.16 (0.14 to 0.19) < 0.001
#>  2:    Treatment Group Control    190             reference       -                 reference       -
#>  3:                     Drug A    288 -0.74 (-1.55 to 0.07)   0.074    -1.37 (-2.07 to -0.68) < 0.001
#>  4:                     Drug B    349   2.90 (2.12 to 3.68) < 0.001       2.95 (2.27 to 3.63) < 0.001
#>  5: Surgical Resection      No    459             reference       -                 reference       -
#>  6:                        Yes    368  0.62 (-0.03 to 1.28)   0.060       3.28 (2.70 to 3.85) < 0.001
#>  7:      Disease Stage       I    207             reference       -                 reference       -
#>  8:                         II    259   1.30 (0.46 to 2.15)   0.003       1.50 (0.81 to 2.19) < 0.001
#>  9:                        III    235   2.68 (1.82 to 3.54) < 0.001       3.00 (2.28 to 3.72) < 0.001
#> 10:                         IV    126   3.04 (2.02 to 4.06) < 0.001       4.35 (3.48 to 5.22) < 0.001

# }
```
