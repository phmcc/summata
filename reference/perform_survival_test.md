# Perform survival comparison test

Performs statistical test comparing survival curves across groups. The
log-rank test is computed by
[`survival::survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html).
The Gehan-Breslow (`"wilcoxon"`), Tarone-Ware (`"tarone"`), and
Peto-Peto (`"petopeto"`) tests weight each event time by the number at
risk, its square root, and the Peto-Peto estimate of the pooled survival
function, respectively (see
[`weighted_logrank_test()`](https://phmcc.codefloe.page/summata/reference/weighted_logrank_test.md)).
[`survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html) offers
only the weights \\\hat{S}(t-)^\rho\\ of the Fleming-Harrington tests,
of which \\\rho = 1\\ approximates the Peto-Peto test and none gives the
other two. Groups without expected events, such as a group whose
observations are all censored before the first event, are left out of
the test, as
[`survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html) leaves
them out, and the degrees of freedom count the groups that remain.

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

List with test statistic, degrees of freedom, *p*-value, and test type.
The statistic, degrees of freedom, and *p*-value are `NA` when fewer
than two groups have expected events.
