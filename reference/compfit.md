# Compare Multiple Regression Models

Fits multiple regression models and provides a comprehensive comparison
table with model quality metrics, convergence diagnostics, and selection
guidance. Computes a composite score combining multiple quality metrics
to facilitate rapid model comparison and selection.

## Usage

``` r
compfit(
  data,
  outcome,
  model_list,
  model_names = NULL,
  interactions_list = NULL,
  random = NULL,
  model_type = "auto",
  family = "binomial",
  conf_level = 0.95,
  p_digits = 3,
  include_coefficients = FALSE,
  scoring_weights = NULL,
  labels = NULL,
  number_format = NULL,
  verbose = NULL,
  ...
)
```

## Arguments

- data:

  Data frame or data.table containing the dataset.

- outcome:

  Character string specifying the outcome variable. For survival
  analysis, use [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
  syntax (*e.g.,* `"Surv(time, status)"`).

- model_list:

  List of character vectors, each containing predictor names for one
  model. Can also be a single character vector to auto-generate nested
  models.

- model_names:

  Character vector of names for each model. If `NULL`, uses "Model 1",
  "Model 2", *etc.* Default is `NULL`.

- interactions_list:

  List of character vectors specifying interaction terms for each model.
  Each element corresponds to one model in model_list. Use `NULL` for
  models without interactions. Use colon notation for interactions
  (*e.g.,* `c("age:treatment")`). If `NULL`, no interactions are added
  to any model. Default is `NULL`.

- random:

  Character string specifying the random-effects formula for
  mixed-effects models (`glmer`, `lmer`, `coxme`). Use standard
  `lme4`/`coxme` syntax, *e.g.,* `"(1|site)"` for random intercepts by
  site. This random effects formula is applied to all models in the
  comparison. Alternatively, random effects can be included directly in
  the predictor vectors within `model_list` using the same syntax, which
  allows different random effects structures across models. Default is
  `NULL`.

- model_type:

  Character string specifying model type. If `"auto"`, detects based on
  outcome. Options include:

  - `"auto"` - Automatically detect based on outcome type (default)

  - `"glm"` - Generalized linear model. Supports logistic, Poisson,
    Gamma, Gaussian via `family` parameter.

  - `"lm"` - Linear regression for continuous outcomes

  - `"coxph"` - Cox proportional hazards for survival analysis

  - `"negbin"` - Negative binomial regression for overdispersed counts
    (requires MASS package)

  - `"lmer"` - Mixed-effects linear regression for clustered continuous
    outcomes

  - `"glmer"` - Mixed-effects logistic regression for clustered
    categorical outcomes

  - `"coxme"` - Mixed-effects Cox regression for clustered time-to-event
    outcomes

- family:

  For GLM and GLMER models, specifies the error distribution and link
  function. Common options include:

  - `"binomial"` - Logistic regression for binary outcomes (default)

  - `"poisson"` - Poisson regression for count data

  - `"quasibinomial"` - Logistic with overdispersion

  - `"quasipoisson"` - Poisson with overdispersion

  - `"gaussian"` - Normal distribution (linear regression via GLM)

  - `"Gamma"` - Gamma for positive continuous data

  - `"inverse.gaussian"` - For positive, highly skewed data

  For negative binomial, use `model_type = "negbin"` instead. See
  [`family`](https://rdrr.io/r/stats/family.html) for all options.

- conf_level:

  Numeric confidence level for intervals. Default is 0.95.

- p_digits:

  Integer specifying the number of decimal places for *p*-values. Values
  smaller than `10^(-p_digits)` are displayed as `"< 0.001"` (for
  `p_digits = 3`), `"< 0.0001"` (for `p_digits = 4`), etc. Default is 3.

- include_coefficients:

  Logical. If TRUE, includes a second table with coefficient estimates.
  Default is FALSE.

- scoring_weights:

  Named list of scoring weights. Each weight should be between 0 and 1,
  and they should sum to 1. Available metrics depend on model type. If
  `NULL`, uses sensible defaults. See Details for available metrics.

- labels:

  Named character vector providing custom display labels for variables.
  Default is `NULL`.

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

  Additional arguments passed to model fitting functions.

## Value

A data.table with class "compfit_result" containing:

- Model:

  Model name/identifier

- CMS:

  Composite Model Score for model selection (higher is better)

- N:

  Sample size

- Events:

  Number of events (for survival/logistic)

- Predictors:

  Number of predictors

- Converged:

  Whether model converged properly

- AIC:

  Akaike Information Criterion

- BIC:

  Bayesian Information Criterion

- *R*\\^2\\ / Pseudo-*R*\\^2\\:

  McFadden pseudo-R-squared (GLM)

- Concordance:

  C-statistic (logistic/survival)

- Brier Score:

  Brier accuracy score (logistic)

- Global *p*:

  Overall model *p*-value

Attributes include:

- models:

  List of fitted model objects

- coefficients:

  Coefficient comparison table (if requested)

- best_model:

  Name of recommended model

## Details

This function fits all specified models and computes comprehensive
quality metrics for comparison. It generates a Composite Model Score
(CMS) that combines multiple metrics: lower AIC/BIC (information
criteria), higher concordance (discrimination), and model convergence
status.

For GLMs, McFadden's pseudo-R-squared is calculated as 1 -
(logLik/logLik_null). For survival models, the global *p*-value comes
from the log-rank test.

Models that fail to converge are flagged and penalized in the composite
score.

**Interaction Terms:**

When `interactions_list` is provided, each element specifies the
interaction terms for the corresponding model in `model_list`. This is
particularly useful for testing whether adding interactions improves
model fit:

- Use `NULL` for models without interactions

- Specify interactions using colon notation:
  `c("age:treatment", "sex:stage")`

- Main effects for all variables in interactions must be in the
  predictor list

- Common pattern: Compare main effects model vs model with interactions

Scoring weights can be customized based on model type:

- GLM: `"convergence"`, `"aic"`, `"concordance"`, `"pseudo_r2"`,
  `"brier"`

- Cox: `"convergence"`, `"aic"`, `"concordance"`, `"global_p"`

- Linear: `"convergence"`, `"aic"`, `"pseudo_r2"`, `"rmse"`

Default weights emphasize discrimination (concordance) and model fit
(AIC).

The composite score is designed as a tool to quickly rank models by
their quality metrics. It should be used alongside traditional model
selection criteria rather than as a definitive model selection method.

## See also

[`fit`](https://phmcc.codeberg.page/summata/reference/fit.md) for
individual model fitting,
[`fullfit`](https://phmcc.codeberg.page/summata/reference/fullfit.md)
for automated variable selection,
[`tablesave`](https://phmcc.codeberg.page/summata/reference/tablesave.md)
for exporting tables

Other regression functions:
[`fit()`](https://phmcc.codeberg.page/summata/reference/fit.md),
[`fullfit()`](https://phmcc.codeberg.page/summata/reference/fullfit.md),
[`multifit()`](https://phmcc.codeberg.page/summata/reference/multifit.md),
[`print.compfit_result()`](https://phmcc.codeberg.page/summata/reference/print.compfit_result.md),
[`print.fit_result()`](https://phmcc.codeberg.page/summata/reference/print.fit_result.md),
[`print.fullfit_result()`](https://phmcc.codeberg.page/summata/reference/print.fullfit_result.md),
[`print.multifit_result()`](https://phmcc.codeberg.page/summata/reference/print.multifit_result.md),
[`print.uniscreen_result()`](https://phmcc.codeberg.page/summata/reference/print.uniscreen_result.md),
[`uniscreen()`](https://phmcc.codeberg.page/summata/reference/uniscreen.md)

## Examples

``` r
# Load example data
data(clintrial)
data(clintrial_labels)

# Information criteria assume a common sample, so the candidate predictors
# are restricted to complete cases before any comparison
comparison_data <- na.omit(
    clintrial[, c("readmission_30d", "os_months", "os_status", "age", "sex",
                  "smoking", "diabetes", "stage", "ecog", "grade",
                  "treatment")]
)

# Example 1: Compare nested logistic regression models
models <- list(
    base = c("age", "sex"),
    clinical = c("age", "sex", "smoking", "diabetes"),
    full = c("age", "sex", "smoking", "diabetes", "stage", "ecog")
)

comparison <- compfit(
    data = comparison_data,
    outcome = "readmission_30d",
    model_list = models,
    model_names = c("Base", "Clinical", "Full")
)
#> Auto-detected binary outcome, using logistic regression
#> Fitting Base with 2 predictors...
#> Fitting Clinical with 4 predictors...
#> Fitting Full with 6 predictors...
comparison
#> 
#> Model Comparison Results
#> Outcome: readmission_30d
#> Model Type: glm
#> 
#> CMS Weights:
#>   Convergence: 15%
#>   AIC: 25%
#>   Concordance: 40%
#>   Pseudo-R²: 15%
#>   Brier score: 5%
#> 
#> Recommended Model: Full (CMS: 66.9)
#> 
#> Models ranked by selection score:
#>       Model   CMS     N Events Predictors Converged    AIC    BIC Pseudo-R² Concordance Brier Score Global p
#>      <char> <num> <int>  <num>      <int>    <char>  <num>  <num>     <num>       <num>       <num>   <char>
#> 1:     Full  66.9   833    416          6       Yes 1069.1 1125.8     0.095       0.692       0.220  < 0.001
#> 2: Clinical  49.3   833    416          4       Yes 1102.1 1130.4     0.056       0.656       0.231  < 0.001
#> 3:     Base  32.1   833    416          2       Yes 1131.1 1145.3     0.026       0.602       0.241  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor

# \donttest{

# Example 2: Compare Cox survival models
surv_models <- list(
    simple = c("age", "sex"),
    clinical = c("age", "sex", "stage", "grade")
)

surv_comparison <- compfit(
    data = comparison_data,
    outcome = "Surv(os_months, os_status)",
    model_list = surv_models,
    model_type = "coxph"
)
#> Fitting simple with 2 predictors...
#> Fitting clinical with 4 predictors...
surv_comparison
#> 
#> Model Comparison Results
#> Outcome: Surv(os_months, os_status)
#> Model Type: coxph
#> 
#> CMS Weights:
#>   Convergence: 15%
#>   AIC: 30%
#>   Concordance: 40%
#>   Global p-value: 15%
#> 
#> Recommended Model: clinical (CMS: 74.1)
#> 
#> Models ranked by selection score:
#>       Model   CMS     N Events Predictors Converged    AIC    BIC Pseudo-R² Concordance Global p
#>      <char> <num> <int>  <num>      <int>    <char>  <num>  <num>     <num>       <num>   <char>
#> 1: clinical  74.1   833    594          4       Yes 7147.4 7178.1     0.218       0.678  < 0.001
#> 2:   simple  38.8   833    594          2       Yes 7257.4 7266.2     0.097       0.612  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor

# Example 3: Test the effect of adding an interaction term.
# Treatment efficacy differs by disease stage in these data, so the
# interaction is expected to improve the fit.
interaction_models <- list(
    main = c("age", "treatment", "stage"),
    interact = c("age", "treatment", "stage")
)

interaction_comp <- compfit(
    data = comparison_data,
    outcome = "Surv(os_months, os_status)",
    model_list = interaction_models,
    model_type = "coxph",
    model_names = c("Main Effects", "With Interaction"),
    interactions_list = list(
        NULL,
        c("treatment:stage")
    )
)
#> Fitting Main Effects with 3 predictors...
#> Fitting With Interaction with 3 predictors + 1 interaction...
interaction_comp
#> 
#> Model Comparison Results
#> Outcome: Surv(os_months, os_status)
#> Model Type: coxph
#> 
#> CMS Weights:
#>   Convergence: 15%
#>   AIC: 30%
#>   Concordance: 40%
#>   Global p-value: 15%
#> 
#> Recommended Model: Main Effects (CMS: 73.9)
#> 
#> Models ranked by selection score:
#>               Model   CMS     N Events Predictors Converged    AIC    BIC Pseudo-R² Concordance Global p
#>              <char> <num> <int>  <num>      <int>    <char>  <num>  <num>     <num>       <num>   <char>
#> 1:     Main Effects  73.9   833    594          3       Yes 7144.4 7170.7     0.219       0.676  < 0.001
#> 2: With Interaction  44.3   833    594          3       Yes 7150.1 7202.7     0.225       0.680  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor

# Example 4: Include coefficient comparison table
detailed <- compfit(
    data = comparison_data,
    outcome = "readmission_30d",
    model_list = models,
    include_coefficients = TRUE,
    labels = clintrial_labels
)
#> Auto-detected binary outcome, using logistic regression
#> Fitting base with 2 predictors...
#> Fitting clinical with 4 predictors...
#> Fitting full with 6 predictors...

# Access coefficient table
coef_table <- attr(detailed, "coefficients")
coef_table
#>        Model                Variable   Group      n Events      aOR (95% CI) p-value
#>       <char>                  <char>  <char> <char> <char>            <char>  <char>
#>  1:     base             Age (years)       -    833    416  1.03 (1.02-1.04) < 0.001
#>  2:     base                     Sex  Female    443    206         reference       -
#>  3:     base                            Male    390    210  1.36 (1.03-1.79)   0.031
#>  4: clinical             Age (years)       -    833    416  1.03 (1.02-1.05) < 0.001
#>  5: clinical                     Sex  Female    443    206         reference       -
#>  6: clinical                            Male    390    210  1.38 (1.04-1.83)   0.027
#>  7: clinical          Smoking Status   Never    337    156         reference       -
#>  8: clinical                          Former    311    141  1.03 (0.75-1.42)   0.862
#>  9: clinical                         Current    185    119  2.44 (1.67-3.59) < 0.001
#> 10: clinical                Diabetes      No    636    300         reference       -
#> 11: clinical                             Yes    197    116  1.81 (1.29-2.54) < 0.001
#> 12:     full             Age (years)       -    833    416  1.04 (1.02-1.05) < 0.001
#> 13:     full                     Sex  Female    443    206         reference       -
#> 14:     full                            Male    390    210  1.41 (1.05-1.89)   0.021
#> 15:     full          Smoking Status   Never    337    156         reference       -
#> 16:     full                          Former    311    141  1.04 (0.75-1.44)   0.833
#> 17:     full                         Current    185    119  2.43 (1.65-3.62) < 0.001
#> 18:     full                Diabetes      No    636    300         reference       -
#> 19:     full                             Yes    197    116  1.82 (1.29-2.58) < 0.001
#> 20:     full           Disease Stage       I    207     86         reference       -
#> 21:     full                              II    261    124  1.34 (0.90-1.99)   0.149
#> 22:     full                             III    237    131  1.95 (1.31-2.92)   0.001
#> 23:     full                              IV    128     75  2.26 (1.41-3.65) < 0.001
#> 24:     full ECOG Performance Status       0    263    107         reference       -
#> 25:     full                               1    298    143  1.37 (0.97-1.95)   0.078
#> 26:     full                               2    235    136  2.06 (1.41-3.00) < 0.001
#> 27:     full                               3     37     30 6.80 (2.91-17.98) < 0.001
#>        Model                Variable   Group      n Events      aOR (95% CI) p-value
#>       <char>                  <char>  <char> <char> <char>            <char>  <char>

# Example 5: Access fitted model objects
fitted_models <- attr(comparison, "models")
names(fitted_models)
#> [1] "Base"     "Clinical" "Full"    

# Example 6: Get best model recommendation
best <- attr(comparison, "best_model")
cat("Recommended model:", best, "\n")
#> Recommended model: Full 

# }
```
