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

- Pseudo-*R*\\^2\\:

  McFadden's pseudo-R-squared (GLM)

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

[`fit`](https://phmcc.github.io/summata/reference/fit.md) for individual
model fitting,
[`fullfit`](https://phmcc.github.io/summata/reference/fullfit.md) for
automated variable selection,
[`table2pdf`](https://phmcc.github.io/summata/reference/table2pdf.md)
for exporting results

Other regression functions:
[`fit()`](https://phmcc.github.io/summata/reference/fit.md),
[`fullfit()`](https://phmcc.github.io/summata/reference/fullfit.md),
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

# Example 1: Compare nested logistic regression models
models <- list(
    base = c("age", "sex"),
    clinical = c("age", "sex", "smoking", "diabetes"),
    full = c("age", "sex", "smoking", "diabetes", "stage", "ecog")
)

comparison <- compfit(
    data = clintrial,
    outcome = "os_status",
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
#> Outcome: os_status
#> Model Type: glm
#> 
#> CMS Weights:
#>   Convergence: 15%
#>   AIC: 25%
#>   Concordance: 40%
#>   Pseudo-R²: 15%
#>   Brier score: 5%
#> 
#> Recommended Model: Full (CMS: 73.2)
#> 
#> Models ranked by selection score:
#>       Model   CMS     N Events Predictors Converged   AIC   BIC Pseudo-R² Concordance Brier Score Global p
#>      <char> <num> <int>  <num>      <int>    <char> <num> <num>     <num>       <num>       <num>   <char>
#> 1:     Full  73.2   833    594          6   Suspect 836.0 892.7     0.187       0.784       0.162  < 0.001
#> 2: Clinical  43.1   833    594          4       Yes 940.4 968.7     0.070       0.680       0.187  < 0.001
#> 3:     Base  39.4   850    609          2       Yes 955.3 969.5     0.064       0.675       0.188  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor

# \donttest{
# Example 2: Compare Cox survival models
library(survival)
surv_models <- list(
    simple = c("age", "sex"),
    clinical = c("age", "sex", "stage", "grade")
)

surv_comparison <- compfit(
    data = clintrial,
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
#> 1: clinical  74.1   838    598          4       Yes 7201.9 7232.6     0.219       0.678  < 0.001
#> 2:   simple  38.9   850    609          2       Yes 7458.5 7467.3     0.101       0.613  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor

# Example 3: Test effect of adding interaction terms
interaction_models <- list(
    main = c("age", "treatment", "sex"),
    interact = c("age", "treatment", "sex")
)

interaction_comp <- compfit(
    data = clintrial,
    outcome = "os_status",
    model_list = interaction_models,
    model_names = c("Main Effects", "With Interaction"),
    interactions_list = list(
        NULL,
        c("treatment:sex")
    )
)
#> Auto-detected binary outcome, using logistic regression
#> Fitting Main Effects with 3 predictors...
#> Fitting With Interaction with 3 predictors + 1 interaction...
interaction_comp
#> 
#> Model Comparison Results
#> Outcome: os_status
#> Model Type: glm
#> 
#> CMS Weights:
#>   Convergence: 15%
#>   AIC: 25%
#>   Concordance: 40%
#>   Pseudo-R²: 15%
#>   Brier score: 5%
#> 
#> Recommended Model: Main Effects (CMS: 66.7)
#> 
#> Models ranked by selection score:
#>               Model   CMS     N Events Predictors Converged   AIC   BIC Pseudo-R² Concordance Brier Score Global p
#>              <char> <num> <int>  <num>      <int>    <char> <num> <num>     <num>       <num>       <num>   <char>
#> 1:     Main Effects  66.7   850    609          3       Yes 943.4 967.2     0.079       0.697       0.185  < 0.001
#> 2: With Interaction  41.7   850    609          3       Yes 947.2 980.4     0.079       0.697       0.184  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor

# Example 4: Include coefficient comparison table
detailed <- compfit(
    data = clintrial,
    outcome = "os_status",
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
#>        Model                Variable   Group      n Events           aOR (95% CI) p-value
#>       <char>                  <char>  <char> <char> <char>                 <char>  <char>
#>  1:     base             Age (years)       -    850    609       1.05 (1.03-1.06) < 0.001
#>  2:     base                     Sex  Female    450    298              reference       -
#>  3:     base                            Male    400    311       1.86 (1.36-2.55) < 0.001
#>  4: clinical             Age (years)       -    833    594       1.05 (1.03-1.06) < 0.001
#>  5: clinical                     Sex  Female    443    292              reference       -
#>  6: clinical                            Male    390    302       1.83 (1.33-2.52) < 0.001
#>  7: clinical          Smoking Status   Never    337    248              reference       -
#>  8: clinical                          Former    311    203       0.74 (0.52-1.05)   0.089
#>  9: clinical                         Current    185    143       1.38 (0.89-2.13)   0.151
#> 10: clinical                Diabetes      No    636    456              reference       -
#> 11: clinical                             Yes    197    138       0.99 (0.69-1.43)   0.974
#> 12:     full             Age (years)       -    833    594       1.06 (1.04-1.07) < 0.001
#> 13:     full                     Sex  Female    443    292              reference       -
#> 14:     full                            Male    390    302       2.00 (1.42-2.83) < 0.001
#> 15:     full          Smoking Status   Never    337    248              reference       -
#> 16:     full                          Former    311    203       0.72 (0.49-1.05)   0.089
#> 17:     full                         Current    185    143       1.32 (0.83-2.10)   0.244
#> 18:     full                Diabetes      No    636    456              reference       -
#> 19:     full                             Yes    197    138       0.93 (0.62-1.39)   0.719
#> 20:     full           Disease Stage       I    207    125              reference       -
#> 21:     full                              II    261    170       1.32 (0.87-2.01)   0.197
#> 22:     full                             III    237    182       2.66 (1.69-4.18) < 0.001
#> 23:     full                              IV    128    117      9.69 (4.75-19.77) < 0.001
#> 24:     full ECOG Performance Status       0    263    158              reference       -
#> 25:     full                               1    298    208       1.70 (1.15-2.50)   0.007
#> 26:     full                               2    235    191       3.35 (2.13-5.25) < 0.001
#> 27:     full                               3     37     37 36022880.19 (0.00-Inf)   0.977
#>        Model                Variable   Group      n Events           aOR (95% CI) p-value
#>       <char>                  <char>  <char> <char> <char>                 <char>  <char>

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
