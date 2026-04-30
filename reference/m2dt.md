# Convert Model to Data Table

Extracts coefficients, confidence intervals, and comprehensive model
statistics from fitted regression models and converts them to a
standardized data.table format suitable for further analysis or
publication. This is a core utility function frequently used internally
by other summata regression functions, although it can be used as a
standalone function as well.

## Usage

``` r
m2dt(
  data,
  model,
  conf_level = 0.95,
  keep_qc_stats = TRUE,
  include_intercept = TRUE,
  terms_to_exclude = NULL,
  reference_rows = TRUE,
  reference_label = "reference",
  skip_counts = FALSE,
  conf_method = NULL
)
```

## Arguments

- data:

  Data frame or data.table containing the dataset used to fit the model.
  Required for computing group-level sample sizes and event counts.

- model:

  Fitted model object. Supported classes include:

  - `glm` - Generalized linear models (logistic, Poisson, *etc.*)

  - `lm` - Linear models

  - `coxph` - Cox proportional hazards models

  - `clogit` - Conditional logistic regression

  - `coxme` - Mixed effects Cox models

  - `lmerMod` - Linear mixed effects models

  - `glmerMod` - Generalized linear mixed effects models

- conf_level:

  Numeric confidence level for confidence intervals. Must be between 0
  and 1. Default is 0.95 (95% CI).

- keep_qc_stats:

  Logical. If `TRUE`, includes model quality statistics such as AIC,
  BIC, *R*\\^2\\, concordance, and model fit tests. These appear as
  additional columns in the output. Default is `TRUE`.

- include_intercept:

  Logical. If `TRUE`, includes the model intercept in output. If
  `FALSE`, removes the intercept row from results. Useful for creating
  cleaner presentation tables. Default is `TRUE`.

- terms_to_exclude:

  Character vector of term names to exclude from output. Useful for
  removing specific unwanted parameters (*e.g.,* nuisance variables,
  spline terms). Default is `NULL`. Note: If
  `include_intercept = FALSE`, "(Intercept)" is automatically added to
  this list.

- reference_rows:

  Logical. If `TRUE`, adds rows for reference categories of factor
  variables with appropriate labels and baseline values (OR/HR = 1,
  Coefficient = 0). This makes tables more complete and easier to
  interpret. Default is `TRUE`.

- reference_label:

  Character string used to label reference category rows in the output.
  Appears in the `reference` column. Default is `"reference"`.

- skip_counts:

  Logical. If `TRUE`, skips computation of group-level sample sizes and
  event counts (faster but less informative). Default is `FALSE`.

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
  use Wald throughout a session.

## Value

A `data.table` containing extracted model information with the following
standard columns:

- model_scope:

  Character. Either "Univariable" (unadjusted model with single
  predictor) or "Multivariable" (adjusted model with multiple
  predictors)

- model_type:

  Character. Type of regression (*e.g.,* "Logistic", "Linear", "Cox PH",
  "Poisson", *etc.*)

- variable:

  Character. Variable name (for factor variables, the base variable name
  without the level)

- group:

  Character. Group/level name for factor variables; empty string for
  continuous variables

- n:

  Integer. Total sample size used in the model

- n_group:

  Integer. Sample size for this specific variable level (factor
  variables only)

- events:

  Integer. Total number of events in the model (for survival and
  logistic models)

- events_group:

  Integer. Number of events for this specific variable level (for
  survival and logistic models with factor variables)

- coefficient:

  Numeric. Raw regression coefficient (log odds, log hazard, *etc.*)

- se:

  Numeric. Standard error of the coefficient

- OR/HR/RR/Coefficient:

  Numeric. Effect estimate - column name depends on model type:

  - `OR` for logistic regression (odds ratio)

  - `HR` for Cox models (hazard ratio)

  - `RR` for Poisson regression (rate/risk ratio)

  - `Coefficient` for linear models or other GLMs

- ci_lower:

  Numeric. Lower bound of confidence interval for effect estimate

- ci_upper:

  Numeric. Upper bound of confidence interval for effect estimate

- statistic:

  Numeric. Test statistic (z-value for GLM/Cox, t-value for LM)

- p_value:

  Numeric. *p*-value for coefficient test

- sig:

  Character. Significance markers: `***` (p \< 0.001), `**` (*p* \<
  0.01), `*` (*p* \< 0.05), `.` (*p* \< 0.10).

- sig_binary:

  Logical. Binary indicator: `TRUE` if *p* \< 0.05, `FALSE` otherwise

- reference:

  Character. Contains `reference_label` for reference category rows when
  `reference_rows = TRUE`, empty string otherwise

## Details

This function is the core extraction utility used by
[`fit()`](https://phmcc.github.io/summata/reference/fit.md) and other
regression functions. It handles the complexities of different model
classes and provides a consistent output format suitable for tables and
forest plots.

**Model Type Detection:** The function automatically detects model type
and applies appropriate:

- Effect measure naming (OR, HR, RR, Coefficient)

- Confidence interval calculation (see below)

- Event counting for binary/survival outcomes

**Confidence Interval Methods:** The CI method is selected per model
class using [`stats::confint()`](https://rdrr.io/r/stats/confint.html)
dispatch:

- **GLM/negative binomial**: Profile likelihood via
  [`MASS::confint.glm()`](https://rdrr.io/pkg/MASS/man/confint.html),
  except quasi-families which use Wald

- **Linear models**: Exact *t*-distribution via
  [`confint.lm()`](https://rdrr.io/r/stats/confint.html)

- **Cox PH**: Wald intervals (coefficient \\\pm\\ *z* \\\times\\ SE)

- **Mixed-effects models**: Wald intervals

Falls back to Wald on profiling failure.

**Mixed Effects Models:** For lme4 models (glmer, lmer), the function
extracts fixed effects only. Random effects variance components are not
included in the output table, as they represent clustering structure
rather than predictor effects.

## See also

[`fit`](https://phmcc.github.io/summata/reference/fit.md) for the main
regression interface,
[`glmforest`](https://phmcc.github.io/summata/reference/glmforest.md),
[`coxforest`](https://phmcc.github.io/summata/reference/coxforest.md),
[`lmforest`](https://phmcc.github.io/summata/reference/lmforest.md) for
forest plot visualization

## Examples

``` r
# Load example data
data(clintrial)

# Example 1: Extract from logistic regression
glm_model <- glm(os_status ~ age + sex + treatment, 
                 data = clintrial, family = binomial)

glm_result <- m2dt(clintrial, glm_model)
glm_result
#>      model_scope model_type             term     n events coefficient          se        coef  coef_lower  coef_upper  exp_coef  exp_lower exp_upper        OR   ci_lower  ci_upper  statistic      p_value      AIC      BIC deviance
#>           <char>     <char>           <char> <int>  <num>       <num>       <num>       <num>       <num>       <num>     <num>      <num>     <num>     <num>      <num>     <num>      <num>        <num>    <num>    <num>    <num>
#> 1: Multivariable   Logistic      (Intercept)   850    609 -1.87542261 0.447739150 -1.87542261 -2.76262522 -1.00520266 0.1532902 0.06312583 0.3659705 0.1532902 0.06312583 0.3659705 -4.1886501 2.806187e-05 943.4251 967.1512 933.4251
#> 2: Multivariable   Logistic              age   850    609  0.04905071 0.007172847  0.04905071  0.03520894  0.06335862 1.0502736 1.03583611 1.0654089 1.0502736 1.03583611 1.0654089  6.8383874 8.008956e-12 943.4251 967.1512 933.4251
#> 3: Multivariable   Logistic        sexFemale   450    298  0.00000000          NA  0.00000000          NA          NA 1.0000000         NA        NA 1.0000000         NA        NA         NA           NA 943.4251 967.1512 933.4251
#> 4: Multivariable   Logistic          sexMale   400    311  0.60465005 0.163104011  0.60465005  0.28714086  0.92708853 1.8306115 1.33261191 2.5271408 1.8306115 1.33261191 2.5271408  3.7071439 2.096098e-04 943.4251 967.1512 933.4251
#> 5: Multivariable   Logistic treatmentControl   196    151  0.00000000          NA  0.00000000          NA          NA 1.0000000         NA        NA 1.0000000         NA        NA         NA           NA 943.4251 967.1512 933.4251
#> 6: Multivariable   Logistic  treatmentDrug A   292    184 -0.74004132 0.217839394 -0.74004132 -1.17367172 -0.31839399 0.4770942 0.30922945 0.7273162 0.4770942 0.30922945 0.7273162 -3.3971877 6.808224e-04 943.4251 967.1512 933.4251
#> 7: Multivariable   Logistic  treatmentDrug B   362    274 -0.14959091 0.217516344 -0.14959091 -0.58168093  0.27246175 0.8610602 0.55895801 1.3131932 0.8610602 0.55895801 1.3131932 -0.6877226 4.916275e-01 943.4251 967.1512 933.4251
#>    null_deviance df_residual c_statistic hoslem_chi2  hoslem_p    variable   group n_group events_group reference    sig sig_binary
#>            <num>       <int>       <num>       <num>     <num>      <char>  <char>   <num>        <num>    <char> <char>     <lgcl>
#> 1:      1013.635         845   0.6968638    10.04661 0.2617695 (Intercept)              NA           NA              ***       TRUE
#> 2:      1013.635         845   0.6968638    10.04661 0.2617695         age              NA           NA              ***       TRUE
#> 3:      1013.635          NA   0.6968638          NA        NA         sex  Female     450          298 reference             FALSE
#> 4:      1013.635         845   0.6968638    10.04661 0.2617695         sex    Male     400          311              ***       TRUE
#> 5:      1013.635          NA   0.6968638          NA        NA   treatment Control     196          151 reference             FALSE
#> 6:      1013.635         845   0.6968638    10.04661 0.2617695   treatment  Drug A     292          184              ***       TRUE
#> 7:      1013.635         845   0.6968638    10.04661 0.2617695   treatment  Drug B     362          274                       FALSE

# \donttest{

# Example 2: Extract from linear model
lm_model <- lm(los_days ~ age + sex + surgery, data = clintrial)

lm_result <- m2dt(clintrial, lm_model)
lm_result
#>      model_scope model_type        term     n events coefficient         se       coef coef_lower coef_upper   exp_coef exp_lower  exp_upper Coefficient  ci_lower   ci_upper statistic      p_value        R2    adj_R2      AIC      BIC
#>           <char>     <char>      <char> <int>  <num>       <num>      <num>      <num>      <num>      <num>      <num>     <num>      <num>       <num>     <num>      <num>     <num>        <num>     <num>     <num>    <num>    <num>
#> 1: Multivariable     Linear (Intercept)   830     NA  10.1480443 0.87333008 10.1480443  8.4338370 11.8622517 10.1480443 8.4338370 11.8622517  10.1480443 8.4338370 11.8622517 11.619941 5.128976e-29 0.1392244 0.1360981 4859.439 4883.046
#> 2: Multivariable     Linear         age   830     NA   0.1487493 0.01344942  0.1487493  0.1223502  0.1751483  0.1487493 0.1223502  0.1751483   0.1487493 0.1223502  0.1751483 11.059904 1.302192e-26 0.1392244 0.1360981 4859.439 4883.046
#> 3: Multivariable     Linear   sexFemale   442     NA   0.0000000         NA  0.0000000         NA         NA  1.0000000        NA         NA   0.0000000        NA         NA        NA           NA        NA        NA 4859.439 4883.046
#> 4: Multivariable     Linear     sexMale   388     NA   0.9071448 0.31347883  0.9071448  0.2918360  1.5224536  0.9071448 0.2918360  1.5224536   0.9071448 0.2918360  1.5224536  2.893799 3.906286e-03 0.1392244 0.1360981 4859.439 4883.046
#> 5: Multivariable     Linear   surgeryNo   460     NA   0.0000000         NA  0.0000000         NA         NA  1.0000000        NA         NA   0.0000000        NA         NA        NA           NA        NA        NA 4859.439 4883.046
#> 6: Multivariable     Linear  surgeryYes   370     NA   1.3070747 0.32101589  1.3070747  0.6769718  1.9371776  1.3070747 0.6769718  1.9371776   1.3070747 0.6769718  1.9371776  4.071682 5.116365e-05 0.1392244 0.1360981 4859.439 4883.046
#>       sigma df_residual    variable  group n_group events_group reference    sig sig_binary
#>       <num>       <int>      <char> <char>   <num>        <num>    <char> <char>     <lgcl>
#> 1: 4.503368         826 (Intercept)             NA           NA              ***       TRUE
#> 2: 4.503368         826         age             NA           NA              ***       TRUE
#> 3:       NA          NA         sex Female     442           NA reference             FALSE
#> 4: 4.503368         826         sex   Male     388           NA               **       TRUE
#> 5:       NA          NA     surgery     No     460           NA reference             FALSE
#> 6: 4.503368         826     surgery    Yes     370           NA              ***       TRUE

# Example 3: Cox proportional hazards model
library(survival)
cox_model <- coxph(Surv(os_months, os_status) ~ age + sex + stage,
                   data = clintrial)

cox_result <- m2dt(clintrial, cox_model)
cox_result
#>      model_scope model_type      term     n events coefficient          se       coef coef_lower coef_upper exp_coef exp_lower exp_upper       HR  ci_lower ci_upper statistic      p_value concordance concordance_se       rsq   rsq_max
#>           <char>     <char>    <char> <int>  <num>       <num>       <num>      <num>      <num>      <num>    <num>     <num>     <num>    <num>     <num>    <num>     <num>        <num>       <num>          <num>     <num>     <num>
#> 1: Multivariable     Cox PH       age   847    606  0.03461725 0.003728664 0.03461725  0.0273092  0.0419253 1.035223 1.0276855  1.042817 1.035223 1.0276855 1.042817  9.284090 1.630945e-20   0.6691193     0.01092972 0.2039982 0.9998581
#> 2: Multivariable     Cox PH sexFemale   450    298  0.00000000          NA 0.00000000         NA         NA 1.000000        NA        NA 1.000000        NA       NA        NA           NA          NA             NA 0.2039982        NA
#> 3: Multivariable     Cox PH   sexMale   400    311  0.30621038 0.081605529 0.30621038  0.1462665  0.4661543 1.358268 1.1575046  1.593853 1.358268 1.1575046 1.593853  3.752324 1.752029e-04   0.6691193     0.01092972 0.2039982 0.9998581
#> 4: Multivariable     Cox PH    stageI   211    127  0.00000000          NA 0.00000000         NA         NA 1.000000        NA        NA 1.000000        NA       NA        NA           NA          NA             NA 0.2039982        NA
#> 5: Multivariable     Cox PH   stageII   263    172  0.12276723 0.117008958 0.12276723 -0.1065661  0.3521006 1.130621 0.8989156  1.422052 1.130621 0.8989156 1.422052  1.049212 2.940804e-01   0.6691193     0.01092972 0.2039982 0.9998581
#> 6: Multivariable     Cox PH  stageIII   241    186  0.59770106 0.115604112 0.59770106  0.3711212  0.8242810 1.817935 1.4493587  2.280241 1.817935 1.4493587 2.280241  5.170241 2.337929e-07   0.6691193     0.01092972 0.2039982 0.9998581
#> 7: Multivariable     Cox PH   stageIV   132    121  1.22555424 0.129235162 1.22555424  0.9722580  1.4788505 3.406053 2.6439076  4.387899 3.406053 2.6439076 4.387899  9.483133 2.467603e-21   0.6691193     0.01092972 0.2039982 0.9998581
#>    likelihood_ratio_test likelihood_ratio_df likelihood_ratio_p wald_test wald_df       wald_p score_test score_df      score_p variable  group n_group events_group reference    sig sig_binary
#>                    <num>               <num>              <num>     <num>   <num>        <num>      <num>    <num>        <num>   <char> <char>   <num>        <num>    <char> <char>     <lgcl>
#> 1:              193.2463                   5       7.903413e-40    198.43       5 6.155672e-41   207.1364        5 8.440913e-43      age             NA           NA              ***       TRUE
#> 2:                    NA                  NA                 NA        NA      NA           NA         NA       NA           NA      sex Female     450          298 reference             FALSE
#> 3:              193.2463                   5       7.903413e-40    198.43       5 6.155672e-41   207.1364        5 8.440913e-43      sex   Male     400          311              ***       TRUE
#> 4:                    NA                  NA                 NA        NA      NA           NA         NA       NA           NA    stage      I     211          127 reference             FALSE
#> 5:              193.2463                   5       7.903413e-40    198.43       5 6.155672e-41   207.1364        5 8.440913e-43    stage     II     263          172                       FALSE
#> 6:              193.2463                   5       7.903413e-40    198.43       5 6.155672e-41   207.1364        5 8.440913e-43    stage    III     241          186              ***       TRUE
#> 7:              193.2463                   5       7.903413e-40    198.43       5 6.155672e-41   207.1364        5 8.440913e-43    stage     IV     132          121              ***       TRUE

# Example 4: Exclude intercept for cleaner tables
clean_result <- m2dt(clintrial, glm_model, include_intercept = FALSE)
clean_result
#>      model_scope model_type             term     n events coefficient          se        coef  coef_lower  coef_upper  exp_coef exp_lower exp_upper        OR  ci_lower  ci_upper  statistic      p_value      AIC      BIC deviance
#>           <char>     <char>           <char> <int>  <num>       <num>       <num>       <num>       <num>       <num>     <num>     <num>     <num>     <num>     <num>     <num>      <num>        <num>    <num>    <num>    <num>
#> 1: Multivariable   Logistic              age   850    609  0.04905071 0.007172847  0.04905071  0.03520894  0.06335862 1.0502736 1.0358361 1.0654089 1.0502736 1.0358361 1.0654089  6.8383874 8.008956e-12 943.4251 967.1512 933.4251
#> 2: Multivariable   Logistic        sexFemale   450    298  0.00000000          NA  0.00000000          NA          NA 1.0000000        NA        NA 1.0000000        NA        NA         NA           NA 943.4251 967.1512 933.4251
#> 3: Multivariable   Logistic          sexMale   400    311  0.60465005 0.163104011  0.60465005  0.28714086  0.92708853 1.8306115 1.3326119 2.5271408 1.8306115 1.3326119 2.5271408  3.7071439 2.096098e-04 943.4251 967.1512 933.4251
#> 4: Multivariable   Logistic treatmentControl   196    151  0.00000000          NA  0.00000000          NA          NA 1.0000000        NA        NA 1.0000000        NA        NA         NA           NA 943.4251 967.1512 933.4251
#> 5: Multivariable   Logistic  treatmentDrug A   292    184 -0.74004132 0.217839394 -0.74004132 -1.17367172 -0.31839399 0.4770942 0.3092294 0.7273162 0.4770942 0.3092294 0.7273162 -3.3971877 6.808224e-04 943.4251 967.1512 933.4251
#> 6: Multivariable   Logistic  treatmentDrug B   362    274 -0.14959091 0.217516344 -0.14959091 -0.58168093  0.27246175 0.8610602 0.5589580 1.3131932 0.8610602 0.5589580 1.3131932 -0.6877226 4.916275e-01 943.4251 967.1512 933.4251
#>    null_deviance df_residual c_statistic hoslem_chi2  hoslem_p  variable   group n_group events_group reference    sig sig_binary
#>            <num>       <int>       <num>       <num>     <num>    <char>  <char>   <num>        <num>    <char> <char>     <lgcl>
#> 1:      1013.635         845   0.6968638    10.04661 0.2617695       age              NA           NA              ***       TRUE
#> 2:      1013.635          NA   0.6968638          NA        NA       sex  Female     450          298 reference             FALSE
#> 3:      1013.635         845   0.6968638    10.04661 0.2617695       sex    Male     400          311              ***       TRUE
#> 4:      1013.635          NA   0.6968638          NA        NA treatment Control     196          151 reference             FALSE
#> 5:      1013.635         845   0.6968638    10.04661 0.2617695 treatment  Drug A     292          184              ***       TRUE
#> 6:      1013.635         845   0.6968638    10.04661 0.2617695 treatment  Drug B     362          274                       FALSE

# Example 5: Change confidence level
result_90ci <- m2dt(clintrial, glm_model, conf_level = 0.90)
result_90ci
#>      model_scope model_type             term     n events coefficient          se        coef coef_lower coef_upper  exp_coef  exp_lower exp_upper        OR   ci_lower  ci_upper  statistic      p_value      AIC      BIC deviance
#>           <char>     <char>           <char> <int>  <num>       <num>       <num>       <num>      <num>      <num>     <num>      <num>     <num>     <num>      <num>     <num>      <num>        <num>    <num>    <num>    <num>
#> 1: Multivariable   Logistic      (Intercept)   850    609 -1.87542261 0.447739150 -1.87542261 -2.6185570 -1.1442562 0.1532902 0.07290799 0.3184607 0.1532902 0.07290799 0.3184607 -4.1886501 2.806187e-05 943.4251 967.1512 933.4251
#> 2: Multivariable   Logistic              age   850    609  0.04905071 0.007172847  0.04905071  0.0374070  0.0610228 1.0502736 1.03811545 1.0629232 1.0502736 1.03811545 1.0629232  6.8383874 8.008956e-12 943.4251 967.1512 933.4251
#> 3: Multivariable   Logistic        sexFemale   450    298  0.00000000          NA  0.00000000         NA         NA 1.0000000         NA        NA 1.0000000         NA        NA         NA           NA 943.4251 967.1512 933.4251
#> 4: Multivariable   Logistic          sexMale   400    311  0.60465005 0.163104011  0.60465005  0.3379278  0.8748420 1.8306115 1.40203927 2.3984964 1.8306115 1.40203927 2.3984964  3.7071439 2.096098e-04 943.4251 967.1512 933.4251
#> 5: Multivariable   Logistic treatmentControl   196    151  0.00000000          NA  0.00000000         NA         NA 1.0000000         NA        NA 1.0000000         NA        NA         NA           NA 943.4251 967.1512 933.4251
#> 6: Multivariable   Logistic  treatmentDrug A   292    184 -0.74004132 0.217839394 -0.74004132 -1.1029725 -0.3855404 0.4770942 0.33188310 0.6800830 0.4770942 0.33188310 0.6800830 -3.3971877 6.808224e-04 943.4251 967.1512 933.4251
#> 7: Multivariable   Logistic  treatmentDrug B   362    274 -0.14959091 0.217516344 -0.14959091 -0.5113443  0.2051031 0.8610602 0.59968888 1.2276517 0.8610602 0.59968888 1.2276517 -0.6877226 4.916275e-01 943.4251 967.1512 933.4251
#>    null_deviance df_residual c_statistic hoslem_chi2  hoslem_p    variable   group n_group events_group reference    sig sig_binary
#>            <num>       <int>       <num>       <num>     <num>      <char>  <char>   <num>        <num>    <char> <char>     <lgcl>
#> 1:      1013.635         845   0.6968638    10.04661 0.2617695 (Intercept)              NA           NA              ***       TRUE
#> 2:      1013.635         845   0.6968638    10.04661 0.2617695         age              NA           NA              ***       TRUE
#> 3:      1013.635          NA   0.6968638          NA        NA         sex  Female     450          298 reference             FALSE
#> 4:      1013.635         845   0.6968638    10.04661 0.2617695         sex    Male     400          311              ***       TRUE
#> 5:      1013.635          NA   0.6968638          NA        NA   treatment Control     196          151 reference             FALSE
#> 6:      1013.635         845   0.6968638    10.04661 0.2617695   treatment  Drug A     292          184              ***       TRUE
#> 7:      1013.635         845   0.6968638    10.04661 0.2617695   treatment  Drug B     362          274                       FALSE

# }
```
