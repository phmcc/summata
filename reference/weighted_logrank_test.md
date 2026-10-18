# Weighted log-rank test

Compares survival curves across groups with a weighted log-rank
statistic. At each distinct event time \\t_j\\, with \\n_j\\ subjects at
risk and \\d_j\\ events in the pooled sample, the observed and expected
numbers of events in each group are compared with weight, with
hypergeometric covariance. The statistic \\v' V^{-1} v\\, computed after
removing one group, has a \\\chi^2\\ distribution with one degree of
freedom fewer than the number of groups. These are the definitions of
the SAS LIFETEST and Stata `sts test` procedures; with unit weights the
statistic equals that of
[`survival::survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html),
and times are aligned with
[`survival::aeqSurv()`](https://rdrr.io/pkg/survival/man/aeqSurv.html)
as [`survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html) aligns
them. Groups with no subjects at risk at any event time are excluded, as
in [`survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html); when
fewer than two groups remain, the statistic and *p*-value are `NA`.

## Usage

``` r
weighted_logrank_test(surv_obj, group_var, weights)
```

## Arguments

- surv_obj:

  Survival object created by
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html), with
  right-censored data.

- group_var:

  Vector of group assignments of the same length.

- weights:

  Function of the numbers at risk and of events at the ordered distinct
  event times of the pooled sample, returning the weight of each event
  time.

## Value

List with the components of a `survdiff` object: the number of subjects
in each group (`n`), the weighted observed and expected numbers of
events in each group (`obs`, `exp`), their covariance matrix (`var`),
the test statistic (`chisq`), and the *p*-value (`pvalue`), with its
degrees of freedom (`df`).
