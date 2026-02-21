# Perform statistical tests for continuous variables

Conducts hypothesis tests comparing continuous variables across groups.
Supports t-tests, Wilcoxon tests, ANOVA, and Kruskal-Wallis tests with
automatic selection based on number of groups.

## Usage

``` r
perform_continuous_test(var_vec, grp_vec, test_type, stat_type)
```

## Arguments

- var_vec:

  Numeric vector of the continuous variable.

- grp_vec:

  Factor or character vector defining groups.

- test_type:

  Character string: "parametric", "nonparametric", or "auto".

- stat_type:

  Character string indicating primary statistic being tested.

## Value

Numeric *p*-value from the hypothesis test.
