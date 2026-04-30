# Perform survival comparison test

Performs statistical test comparing survival curves across groups.

## Usage

``` r
perform_survival_test(surv_obj, group_var, test_type)
```

## Arguments

- surv_obj:

  Survival object created by Surv().

- group_var:

  Vector of group assignments.

- test_type:

  Character string specifying test type.

## Value

List with test statistic, *p*-value, and test type.
