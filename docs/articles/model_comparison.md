# Model Comparison

Model selection is a fundamental problem in applied statistics. Given a
set of candidate models, it is important to determine which
specification best balances goodness-of-fit against parsimony.
Traditional approaches rely on information criteria (AIC, BIC),
discrimination metrics (C-statistic), and hypothesis tests for nested
models.

The [`compfit()`](https://phmcc.github.io/summata/reference/compfit.md)
function synthesizes these metrics into a weighted Composite Model Score
(CMS) to facilitate systematic comparison between models. Like other
functions in this package, it follows the standard `summata` input
paradigm, with an additional parameter to allow for multiple models:

``` r
compfit(data, outcome, model_list, model_type, ...)
```

where `data` is the dataset, `outcome` is the common endpoint variable,
`model_list` is a named list of predictor vectors (each representing a
different model to be compared), and `model_type` is the type of model
to be implemented. This vignette demonstrates the various capabilities
of this function using the included sample dataset.

------------------------------------------------------------------------

## Preliminaries

The examples in this vignette use the `clintrial` dataset included with
`summata`:

``` r
library(summata)

data(clintrial)
data(clintrial_labels)
```

------------------------------------------------------------------------

## Model Quality Metrics

The [`compfit()`](https://phmcc.github.io/summata/reference/compfit.md)
function uses several metrics for comparing models. For mathematical
details, see the [Statistical
Foundations](https://phmcc.github.io/summata/articles/statistical_foundations.md)
article.

| Metric      | Interpretation                                         | Better |
|:------------|:-------------------------------------------------------|:-------|
| AIC         | Information criterion balancing fit and complexity     | Lower  |
| BIC         | Information criterion with stronger complexity penalty | Lower  |
| Concordance | Discrimination ability (0.5 = random, 1.0 = perfect)   | Higher |
| Pseudo-*R*² | Proportion of variation explained (GLMs)               | Higher |
| Brier Score | Prediction accuracy for binary outcomes                | Lower  |

------------------------------------------------------------------------

## Basic Usage

### **Example 1:** Nested Model Comparison

Compare models with increasing complexity:

``` r
example1 <- compfit(
  data = clintrial,
  outcome = "surgery",
  model_list = list(
    "Demographics" = c("age", "sex"),
    "Plus Stage" = c("age", "sex", "stage", "ecog"),
    "Full Model" = c("age", "sex", "stage", "ecog", "treatment", "smoking")
  ),
  model_type = "glm"
)

example1
#> 
#> Model Comparison Results
#> Outcome: surgery
#> Model Type: glm
#> 
#> CMS Weights:
#>   Convergence: 15%
#>   AIC: 25%
#>   Concordance: 40%
#>   Pseudo-R²: 15%
#>   Brier score: 5%
#> 
#> Recommended Model: Full Model (CMS: 81.2)
#> 
#> Models ranked by selection score:
#>           Model   CMS     N Events Predictors Converged    AIC    BIC Pseudo-R² Concordance Brier Score Global p
#>          <char> <num> <int>  <num>      <int>    <char>  <num>  <num>     <num>       <num>       <num>   <char>
#> 1:   Full Model  81.2   833    362          6       Yes  900.0  961.4     0.234       0.811       0.175  < 0.001
#> 2:   Plus Stage  68.9   842    368          4       Yes  969.9 1012.5     0.175       0.773       0.192  < 0.001
#> 3: Demographics  33.6   850    370          2       Yes 1132.9 1147.1     0.032       0.617       0.235  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor
```

### **Example 2:** Linear Regression Models

For continuous outcomes, use `model_type = "lm"`:

``` r
example2 <- compfit(
  data = clintrial,
  outcome = "los_days",
  model_list = list(
    "Simple" = c("age", "sex"),
    "Disease" = c("age", "sex", "stage", "ecog"),
    "Treatment" = c("age", "sex", "stage", "ecog", "surgery", "treatment")
  ),
  model_type = "lm"
)

example2
#> 
#> Model Comparison Results
#> Outcome: los_days
#> Model Type: lm
#> 
#> CMS Weights:
#>   Convergence: 15%
#>   AIC: 25%
#>   Pseudo-R²: 45%
#>   RMSE: 15%
#> 
#> Recommended Model: Treatment (CMS: 72.2)
#> 
#> Models ranked by selection score:
#>        Model   CMS     N Events Predictors Converged    AIC    BIC Pseudo-R² Concordance Global p
#>       <char> <num> <int>  <num>      <int>    <char>  <num>  <num>     <num>       <num>   <char>
#> 1: Treatment  72.2   822     NA          6       Yes 4297.3 4358.6     0.550          NA  < 0.001
#> 2:   Disease  44.2   822     NA          4       Yes 4670.3 4717.5     0.287          NA  < 0.001
#> 3:    Simple  28.0   830     NA          2       Yes 4873.9 4892.8     0.122          NA  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor
```

### **Example 3:** Cox Regression Models

For time-to-event outcomes, use `model_type = "coxph"`:

``` r
example3 <- compfit(
  data = clintrial,
  outcome = "Surv(os_months, os_status)",
  model_list = list(
    "Unadjusted" = c("treatment"),
    "Demographics" = c("treatment", "age", "sex"),
    "Full" = c("treatment", "age", "sex", "stage", "ecog")
  ),
  model_type = "coxph"
)

example3
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
#> Recommended Model: Full (CMS: 75.8)
#> 
#> Models ranked by selection score:
#>           Model   CMS     N Events Predictors Converged    AIC    BIC Pseudo-R² Concordance Global p
#>          <char> <num> <int>  <num>      <int>    <char>  <num>  <num>     <num>       <num>   <char>
#> 1:         Full  75.8   196    151          5       Yes 7188.4 7232.4     0.283       0.699  < 0.001
#> 2: Demographics  48.0   196    151          3       Yes 7438.4 7456.0     0.126       0.629  < 0.001
#> 3:   Unadjusted  33.8   196    151          1       Yes 7527.1 7535.9     0.026       0.549  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor
```

### **Example 4:** Count Models

For count outcomes, use `model_type = "glm"` with a count model family
(e.g., `family = "poisson"`):

``` r
example4 <- compfit(
  data = clintrial,
  outcome = "fu_count",
  model_list = list(
    "Minimal" = c("age", "treatment"),
    "Clinical" = c("age", "treatment", "stage", "ecog"),
    "Full" = c("age", "treatment", "stage", "ecog", "surgery", "diabetes")
  ),
  model_type = "glm",
  family = "poisson",
  labels = clintrial_labels
)

example4
#> 
#> Model Comparison Results
#> Outcome: fu_count
#> Model Type: glm
#> 
#> CMS Weights:
#>   Convergence: 15%
#>   AIC: 25%
#>   Concordance: 40%
#>   Pseudo-R²: 15%
#>   Brier score: 5%
#> 
#> Recommended Model: Full (CMS: 67.5)
#> 
#> Models ranked by selection score:
#>       Model   CMS     N Events Predictors Converged    AIC    BIC Pseudo-R² Concordance Brier Score Global p
#>      <char> <num> <int>  <num>      <int>    <char>  <num>  <num>     <num>       <num>       <num>   <char>
#> 1:     Full  67.5   826   5410          6       Yes 3956.9 4013.5        NA          NA          NA     <NA>
#> 2: Clinical  57.0   834   5482          4       Yes 4009.2 4056.4        NA          NA          NA     <NA>
#> 3:  Minimal  42.5   842   5533          2       Yes 4081.4 4100.4        NA          NA          NA     <NA>
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor
```

------------------------------------------------------------------------

## Interaction Testing

A key application of
[`compfit()`](https://phmcc.github.io/summata/reference/compfit.md) is
testing whether interaction terms improve model fit.

### **Example 5:** Single Interaction

Compare a main-effects model to one with an interaction term:

``` r
example5 <- compfit(
  data = clintrial,
  outcome = "surgery",
  model_list = list(
    "Main Effects" = c("age", "treatment", "sex"),
    "With Interaction" = c("age", "treatment", "sex")
  ),
  interactions_list = list(
    NULL,
    c("sex:treatment")
  ),
  model_type = "glm"
)

example5
#> 
#> Model Comparison Results
#> Outcome: surgery
#> Model Type: glm
#> 
#> CMS Weights:
#>   Convergence: 15%
#>   AIC: 25%
#>   Concordance: 40%
#>   Pseudo-R²: 15%
#>   Brier score: 5%
#> 
#> Recommended Model: Main Effects (CMS: 67.1)
#> 
#> Models ranked by selection score:
#>               Model   CMS     N Events Predictors Converged    AIC    BIC Pseudo-R² Concordance Brier Score Global p
#>              <char> <num> <int>  <num>      <int>    <char>  <num>  <num>     <num>       <num>       <num>   <char>
#> 1:     Main Effects  67.1   850    370          3       Yes 1071.7 1095.5     0.088       0.697       0.218  < 0.001
#> 2: With Interaction  42.1   850    370          3       Yes 1075.6 1108.8     0.088       0.697       0.218  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor
```

### **Example 6:** Multiple Interactions

Compare different interaction hypotheses:

``` r
example6 <- compfit(
  data = clintrial,
  outcome = "Surv(os_months, os_status)",
  model_list = list(
    "Main Effects" = c("age", "treatment", "sex", "stage"),
    "Age × Treatment" = c("age", "treatment", "sex", "stage"),
    "Sex × Treatment" = c("age", "treatment", "sex", "stage"),
    "Both" = c("age", "treatment", "sex", "stage")
  ),
  interactions_list = list(
    NULL,
    c("age:treatment"),
    c("sex:treatment"),
    c("age:treatment", "sex:treatment")
  ),
  model_type = "coxph"
)

example6
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
#> Recommended Model: Main Effects (CMS: 74.2)
#> 
#> Models ranked by selection score:
#>              Model   CMS     N Events Predictors Converged    AIC    BIC Pseudo-R² Concordance Global p
#>             <char> <num> <int>  <num>      <int>    <char>  <num>  <num>     <num>       <num>   <char>
#> 1:    Main Effects  74.2   847    606          4       Yes 7295.8 7326.7     0.232        0.68  < 0.001
#> 2: Age × Treatment  60.2   847    606          4       Yes 7298.8 7338.4     0.232        0.68  < 0.001
#> 3: Sex × Treatment  57.8   847    606          4       Yes 7299.3 7339.0     0.232        0.68  < 0.001
#> 4:            Both  44.3   847    606          4       Yes 7302.2 7350.7     0.233        0.68  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor
```

------------------------------------------------------------------------

## Accessing Detailed Results

### **Example 7:** Coefficient Comparison

Setting `include_coefficients = TRUE` generates a table comparing
coefficients across models:

``` r
example7 <- compfit(
  data = clintrial,
  outcome = "surgery",
  model_list = list(
    "Model A" = c("age", "sex"),
    "Model B" = c("age", "sex", "stage"),
    "Model C" = c("age", "sex", "stage", "treatment")
  ),
  model_type = "glm",
  include_coefficients = TRUE,
  labels = clintrial_labels
)

# Main comparison
example7
#> 
#> Model Comparison Results
#> Outcome: surgery
#> Model Type: glm
#> 
#> CMS Weights:
#>   Convergence: 15%
#>   AIC: 25%
#>   Concordance: 40%
#>   Pseudo-R²: 15%
#>   Brier score: 5%
#> 
#> Recommended Model: Model C (CMS: 75.4)
#> 
#> Models ranked by selection score:
#>      Model   CMS     N Events Predictors Converged    AIC    BIC Pseudo-R² Concordance Brier Score Global p
#>     <char> <num> <int>  <num>      <int>    <char>  <num>  <num>     <num>       <num>       <num>   <char>
#> 1: Model C  75.4   847    368          4       Yes  984.4 1022.4     0.165       0.765       0.195  < 0.001
#> 2: Model B  59.4   847    368          3       Yes 1045.0 1073.5     0.109       0.719       0.211  < 0.001
#> 3: Model A  33.6   850    370          2       Yes 1132.9 1147.1     0.032       0.617       0.235  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor
#> 
#> Note: Coefficient comparison available via attr(result, 'coefficients')

# Coefficient table
coef_table <- attr(example7, "coefficients")
coef_table
#>       Model        Variable   Group      n Events     aOR (95% CI) p-value
#>      <char>          <char>  <char> <char> <char>           <char>  <char>
#>  1: Model A     Age (years)       -    850    370 0.96 (0.95-0.98) < 0.001
#>  2: Model A             Sex  Female    450    204        reference       -
#>  3: Model A                    Male    400    166 0.85 (0.64-1.12)   0.252
#>  4: Model B     Age (years)       -    847    368 0.96 (0.95-0.97) < 0.001
#>  5: Model B             Sex  Female    449    203        reference       -
#>  6: Model B                    Male    398    165 0.79 (0.59-1.06)   0.123
#>  7: Model B   Disease Stage       I    211    125        reference       -
#>  8: Model B                      II    263    133 0.68 (0.47-0.99)   0.046
#>  9: Model B                     III    241     91 0.37 (0.25-0.55) < 0.001
#> 10: Model B                      IV    132     19 0.10 (0.06-0.18) < 0.001
#> 11: Model C     Age (years)       -    847    368 0.96 (0.94-0.97) < 0.001
#> 12: Model C             Sex  Female    449    203        reference       -
#> 13: Model C                    Male    398    165 0.84 (0.62-1.14)   0.263
#> 14: Model C   Disease Stage       I    211    125        reference       -
#> 15: Model C                      II    263    133 0.69 (0.47-1.03)   0.070
#> 16: Model C                     III    241     91 0.41 (0.27-0.62) < 0.001
#> 17: Model C                      IV    132     19 0.09 (0.05-0.17) < 0.001
#> 18: Model C Treatment Group Control    194     93        reference       -
#> 19: Model C                  Drug A    292    173 1.84 (1.23-2.75)   0.003
#> 20: Model C                  Drug B    361    102 0.45 (0.30-0.66) < 0.001
#>       Model        Variable   Group      n Events     aOR (95% CI) p-value
#>      <char>          <char>  <char> <char> <char>           <char>  <char>
```

### **Example 8:** Fitted Model Objects

The underlying model objects are stored as attributes for further
analysis:

``` r
models <- attr(example7, "models")
names(models)
#> [1] "Model A" "Model B" "Model C"

# Examine a specific model
summary(models[["Model C"]])
#> 
#> Call:
#> stats::glm(formula = formula, family = family, data = data)
#> 
#> Coefficients:
#>                  Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)      3.195791   0.470598   6.791 1.11e-11 ***
#> age             -0.043364   0.006879  -6.304 2.91e-10 ***
#> sexMale         -0.174895   0.156321  -1.119  0.26322    
#> stageII         -0.364326   0.201174  -1.811  0.07014 .  
#> stageIII        -0.892278   0.209902  -4.251 2.13e-05 ***
#> stageIV         -2.384440   0.302499  -7.882 3.21e-15 ***
#> treatmentDrug A  0.609610   0.205036   2.973  0.00295 ** 
#> treatmentDrug B -0.805077   0.200419  -4.017 5.89e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 1159.60  on 846  degrees of freedom
#> Residual deviance:  968.44  on 839  degrees of freedom
#>   (3 observations deleted due to missingness)
#> AIC: 984.44
#> 
#> Number of Fisher Scoring iterations: 4
```

### **Example 9:** Recommended Model

The best-performing model is identified automatically based on the
Composite Model Score (CMS):

``` r
example9 <- compfit(
  data = clintrial,
  outcome = "Surv(os_months, os_status)",
  model_list = list(
    "Minimal" = c("treatment"),
    "Standard" = c("treatment", "age", "sex", "stage"),
    "Extended" = c("treatment", "age", "sex", "stage", "ecog", "grade")
  ),
  model_type = "coxph"
)

recommended <- attr(example9, "best_model")
cat("Recommended model:", recommended, "\n")
#> Recommended model: Extended
```

------------------------------------------------------------------------

## The Composite Model Score (CMS)

The Composite Model Score (CMS) is a metric designed to facilitate rapid
model comparison. It combines multiple quality indicators into a single
value ranging from 0 to 100.

### Score Interpretation

| Range  | Interpretation |
|:-------|:---------------|
| 85–100 | Excellent      |
| 75–84  | Very Good      |
| 65–74  | Good           |
| 55–64  | Fair           |
| \< 55  | Poor           |

### Default Weights by Model Type

The score components and their default weights vary by model type.

**Linear regression:**

| Component     | Weight |
|:--------------|:-------|
| Convergence   | 15%    |
| AIC           | 25%    |
| Adjusted *R²* | 45%    |
| RMSE          | 15%    |

**Logistic regression:**

| Component   | Weight |
|:------------|:-------|
| Convergence | 15%    |
| AIC         | 25%    |
| Concordance | 40%    |
| Pseudo-*R²* | 15%    |
| Brier Score | 5%     |

**Cox regression:**

| Component   | Weight |
|:------------|:-------|
| Convergence | 15%    |
| AIC         | 30%    |
| Concordance | 40%    |
| Global *p*  | 15%    |

### **Example 10:** Custom Weights

Modify the CMS scoring weights to emphasize specific metrics:

``` r
example10 <- compfit(
  data = clintrial,
  outcome = "surgery",
  model_list = list(
    "Simple" = c("age", "sex"),
    "Standard" = c("age", "sex", "stage"),
    "Full" = c("age", "sex", "stage", "treatment", "ecog")
  ),
  model_type = "glm",
  scoring_weights = list(
    convergence = 0.05,
    aic = 0.20,
    concordance = 0.60,
    pseudo_r2 = 0.10,
    brier = 0.05
  )
)

example10
#> 
#> Model Comparison Results
#> Outcome: surgery
#> Model Type: glm
#> 
#> CMS Weights:
#>   Convergence: 5%
#>   AIC: 20%
#>   Concordance: 60%
#>   Pseudo-R²: 10%
#>   Brier score: 5%
#> 
#> Recommended Model: Full (CMS: 79.5)
#> 
#> Models ranked by selection score:
#>       Model   CMS     N Events Predictors Converged    AIC    BIC Pseudo-R² Concordance Brier Score Global p
#>      <char> <num> <int>  <num>      <int>    <char>  <num>  <num>     <num>       <num>       <num>   <char>
#> 1:     Full  79.5   842    368          5       Yes  904.8  956.9     0.235       0.811       0.175  < 0.001
#> 2: Standard  53.7   847    368          3       Yes 1045.0 1073.5     0.109       0.719       0.211  < 0.001
#> 3:   Simple  31.8   850    370          2       Yes 1132.9 1147.1     0.032       0.617       0.235  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor
```

------------------------------------------------------------------------

## Application Scenarios

### **Scenario 1:** Confounder Assessment

Evaluate the impact of covariate adjustment on effect estimates:

``` r
scenario1 <- compfit(
  data = clintrial,
  outcome = "Surv(os_months, os_status)",
  model_list = list(
    "Crude" = c("treatment"),
    "Age-Sex Adjusted" = c("treatment", "age", "sex"),
    "Fully Adjusted" = c("treatment", "age", "sex", "stage", "ecog")
  ),
  model_type = "coxph",
  include_coefficients = TRUE,
  labels = clintrial_labels
)

scenario1
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
#> Recommended Model: Fully Adjusted (CMS: 75.8)
#> 
#> Models ranked by selection score:
#>               Model   CMS     N Events Predictors Converged    AIC    BIC Pseudo-R² Concordance Global p
#>              <char> <num> <int>  <num>      <int>    <char>  <num>  <num>     <num>       <num>   <char>
#> 1:   Fully Adjusted  75.8   196    151          5       Yes 7188.4 7232.4     0.283       0.699  < 0.001
#> 2: Age-Sex Adjusted  48.0   196    151          3       Yes 7438.4 7456.0     0.126       0.629  < 0.001
#> 3:            Crude  33.8   196    151          1       Yes 7527.1 7535.9     0.026       0.549  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor
#> 
#> Note: Coefficient comparison available via attr(result, 'coefficients')

# Compare treatment effect across models
attr(scenario1, "coefficients")
#>                Model                Variable   Group      n Events      HR (95% CI) p-value     aHR (95% CI)
#>               <char>                  <char>  <char> <char> <char>           <char>  <char>           <char>
#>  1:            Crude         Treatment Group Control    196    151        reference       -             <NA>
#>  2:            Crude                          Drug A    292    184 0.64 (0.52-0.80) < 0.001             <NA>
#>  3:            Crude                          Drug B    362    274 0.94 (0.77-1.15)   0.567             <NA>
#>  4: Age-Sex Adjusted         Treatment Group Control    196    151             <NA>       -        reference
#>  5: Age-Sex Adjusted                          Drug A    292    184             <NA> < 0.001 0.62 (0.50-0.77)
#>  6: Age-Sex Adjusted                          Drug B    362    274             <NA>   0.324 0.90 (0.74-1.10)
#>  7: Age-Sex Adjusted             Age (years)       -    850    609             <NA> < 0.001 1.03 (1.03-1.04)
#>  8: Age-Sex Adjusted                     Sex  Female    450    298             <NA>       -        reference
#>  9: Age-Sex Adjusted                            Male    400    311             <NA>   0.004 1.26 (1.08-1.48)
#> 10:   Fully Adjusted         Treatment Group Control    196    151             <NA>       -        reference
#> 11:   Fully Adjusted                          Drug A    292    184             <NA> < 0.001 0.54 (0.44-0.68)
#> 12:   Fully Adjusted                          Drug B    362    274             <NA>   0.095 0.84 (0.69-1.03)
#> 13:   Fully Adjusted             Age (years)       -    842    602             <NA> < 0.001 1.04 (1.03-1.04)
#> 14:   Fully Adjusted                     Sex  Female    450    298             <NA>       -        reference
#> 15:   Fully Adjusted                            Male    400    311             <NA>   0.002 1.29 (1.10-1.52)
#> 16:   Fully Adjusted           Disease Stage       I    211    127             <NA>       -        reference
#> 17:   Fully Adjusted                              II    263    172             <NA>   0.288 1.13 (0.90-1.43)
#> 18:   Fully Adjusted                             III    241    186             <NA> < 0.001 1.93 (1.53-2.43)
#> 19:   Fully Adjusted                              IV    132    121             <NA> < 0.001 3.87 (2.98-5.01)
#> 20:   Fully Adjusted ECOG Performance Status       0    265    159             <NA>       -        reference
#> 21:   Fully Adjusted                               1    302    212             <NA> < 0.001 1.51 (1.22-1.85)
#> 22:   Fully Adjusted                               2    238    194             <NA> < 0.001 1.93 (1.56-2.39)
#> 23:   Fully Adjusted                               3     37     37             <NA> < 0.001 3.33 (2.31-4.80)
#>                Model                Variable   Group      n Events      HR (95% CI) p-value     aHR (95% CI)
#>               <char>                  <char>  <char> <char> <char>           <char>  <char>           <char>
```

### **Scenario 2:** Variable Selection Validation

Compare automated versus theory-driven selection:

``` r
# Identify candidates via screening
screening <- uniscreen(
  data = clintrial,
  outcome = "surgery",
  predictors = c("age", "sex", "bmi", "smoking", "diabetes",
                 "hypertension", "stage", "ecog", "treatment"),
  model_type = "glm",
  p_threshold = 0.10
)

# Extract significant predictors
sig_vars <- attr(screening, "significant")

scenario2 <- compfit(
  data = clintrial,
  outcome = "surgery",
  model_list = list(
    "Theory-Driven" = c("age", "sex", "stage", "treatment"),
    "Data-Driven" = sig_vars,
    "Combined" = unique(c("age", "sex", "stage", "treatment", sig_vars))
  ),
  model_type = "glm"
)

scenario2
#> 
#> Model Comparison Results
#> Outcome: surgery
#> Model Type: glm
#> 
#> CMS Weights:
#>   Convergence: 15%
#>   AIC: 25%
#>   Concordance: 40%
#>   Pseudo-R²: 15%
#>   Brier score: 5%
#> 
#> Recommended Model: Data-Driven (CMS: 81.3)
#> 
#> Models ranked by selection score:
#>            Model   CMS     N Events Predictors Converged   AIC    BIC Pseudo-R² Concordance Brier Score Global p
#>           <char> <num> <int>  <num>      <int>    <char> <num>  <num>     <num>       <num>       <num>   <char>
#> 1:   Data-Driven  81.3   842    368          4       Yes 903.2  950.5     0.235       0.811       0.176  < 0.001
#> 2:      Combined  80.8   842    368          5       Yes 904.8  956.9     0.235       0.811       0.175  < 0.001
#> 3: Theory-Driven  50.4   847    368          4       Yes 984.4 1022.4     0.165       0.765       0.195  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor
```

### **Scenario 3:** Parsimony Assessment

Test whether additional predictors meaningfully improve fit:

``` r
scenario3 <- compfit(
  data = clintrial,
  outcome = "los_days",
  model_list = list(
    "3 Predictors" = c("age", "surgery", "ecog"),
    "5 Predictors" = c("age", "surgery", "ecog", "stage", "treatment"),
    "8 Predictors" = c("age", "surgery", "ecog", "stage", "treatment",
                       "sex", "smoking", "diabetes")
  ),
  model_type = "lm",
  labels = clintrial_labels
)

scenario3
#> 
#> Model Comparison Results
#> Outcome: los_days
#> Model Type: lm
#> 
#> CMS Weights:
#>   Convergence: 15%
#>   AIC: 25%
#>   Pseudo-R²: 45%
#>   RMSE: 15%
#> 
#> Recommended Model: 8 Predictors (CMS: 75.1)
#> 
#> Models ranked by selection score:
#>           Model   CMS     N Events Predictors Converged    AIC    BIC Pseudo-R² Concordance Global p
#>          <char> <num> <int>  <num>      <int>    <char>  <num>  <num>     <num>       <num>   <char>
#> 1: 8 Predictors  75.1   813     NA          8       Yes 4135.8 4211.0     0.613          NA  < 0.001
#> 2: 5 Predictors  64.1   822     NA          5       Yes 4310.6 4367.1     0.542          NA  < 0.001
#> 3: 3 Predictors  34.1   822     NA          3       Yes 4696.4 4729.4     0.258          NA  < 0.001
#> 
#> CMS interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor
```

------------------------------------------------------------------------

## Exporting Results

Comparison tables can be exported to various formats:

``` r
# Main comparison table
table2docx(
  table = example1,
  file = "Model_Comparison.docx",
  caption = "Table 3. Model Comparison Results"
)

# Coefficient table
table2docx(
  table = attr(example6, "coefficients"),
  file = "Coefficient_Comparison.docx",
  caption = "Table S1. Coefficient Estimates Across Models"
)

# PDF with landscape orientation for wide tables
table2pdf(
  table = example1,
  file = "Model_Comparison.pdf",
  caption = "Model Comparison",
  orientation = "landscape"
)
```

------------------------------------------------------------------------

## Best Practices

### When to Compare Models

1.  **Covariate selection**: Determine appropriate adjustment level
2.  **Interaction testing**: Evaluate effect modification
3.  **Parsimony**: Balance complexity against fit
4.  **Sensitivity analysis**: Compare different specifications

### Interpreting Results

1.  **Multiple metrics**: The Composite Model Score (CMS) provides a
    summary, but examine individual metrics
2.  **Practical significance**: Statistical improvement may not
    translate to meaningful differences
3.  **Sample size**: Complex models require larger samples
4.  **Convergence**: Non-converged models should be interpreted
    cautiously

### Limitations

1.  The Composite Model Score (CMS) is a heuristic, not a formal
    statistical test
2.  Comparisons assume models are fit on identical observations
3.  Information criteria are most meaningful for nested models
4.  Small score differences may not be practically significant

------------------------------------------------------------------------

## Common Issues

### Non-Convergence

Check convergence status and consider simplifying non-converging models:

``` r
# Check convergence status
comparison[, .(Model, Converged)]

# For non-converging models:
# 1. Reduce complexity
# 2. Check for separation (logistic)
# 3. Examine predictor correlations
```

When encountering non-converging models, consider reducing complexity,
examining predictor correlations, and checking for perfect separation
(primarily in logistic models).

### Differing Sample Sizes

Models with missing data may have different sample sizes, which
complicates comparison:

``` r
# Check sample sizes
comparison[, .(Model, N, Events)]

# Use complete cases for fair comparison
complete_data <- na.omit(data[, relevant_vars, with = FALSE])
```

### Interpreting Close Scores

When scores are similar, prefer parsimony and consider domain knowledge:

``` r
# Examine individual metrics
comparison[, .(Model, `Composite Model Score (CMS)`, AIC, Concordance)]

# Prefer parsimony when scores are close
# Consider interpretability
```

------------------------------------------------------------------------

## Further Reading

- [Statistical
  Foundations](https://phmcc.github.io/summata/articles/statistical_foundations.md):
  Mathematical details for model quality metrics
- [Descriptive
  Tables](https://phmcc.github.io/summata/articles/descriptive_tables.md):
  [`desctable()`](https://phmcc.github.io/summata/reference/desctable.md)
  for baseline characteristics
- [Regression
  Modeling](https://phmcc.github.io/summata/articles/regression_modeling.md):
  [`fit()`](https://phmcc.github.io/summata/reference/fit.md),
  [`uniscreen()`](https://phmcc.github.io/summata/reference/uniscreen.md),
  and
  [`fullfit()`](https://phmcc.github.io/summata/reference/fullfit.md)
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
