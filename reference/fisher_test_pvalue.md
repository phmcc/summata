# Fisher's exact test with a Monte Carlo fallback

Computes the *p*-value of Fisher's exact test. For tables larger than 2
\\\times\\ 2,
[`stats::fisher.test()`](https://rdrr.io/r/stats/fisher.test.html) uses
a network algorithm with a bounded workspace (`workspace = 200000` by
default), which is exhausted by many tables with a few hundred
observations or more. A larger workspace is not a practical remedy,
because computation time grows steeply and the algorithm still fails for
many such tables. When the exact computation fails, the *p*-value is
estimated by Monte Carlo simulation conditional on the margins, as
[`fisher.test()`](https://rdrr.io/r/stats/fisher.test.html) itself
recommends, and a warning reports the substitution.

## Usage

``` r
fisher_test_pvalue(tab, args = list(), var_name = NULL)
```

## Arguments

- tab:

  Contingency table without missing-value rows or columns.

- args:

  Named list of arguments for
  [`fisher.test()`](https://rdrr.io/r/stats/fisher.test.html). With
  `simulate.p.value = TRUE`, the simulation is used directly. A supplied
  `B` sets the number of simulated tables for the fallback, which
  otherwise uses 10,000.

- var_name:

  Character string naming the variable, used in the warning.

## Value

Numeric *p*-value.
