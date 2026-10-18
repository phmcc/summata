# Wilcoxon rank sum test p-value

Computes the two-sided *p*-value of the Wilcoxon rank sum test for two
groups. When either group has 50 or more observations,
[`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html) uses
the normal approximation with continuity correction and tie correction
by default; that computation is reproduced here with ranks and tie
counts from
[`average_ranks()`](https://phmcc.codefloe.page/summata/reference/average_ranks.md),
because [`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html)
ranks with [`rank()`](https://rdrr.io/r/base/rank.html) and counts ties
with [`table()`](https://rdrr.io/r/base/table.html), which converts
every rank to a character string; on large samples these two steps
dominate its running time. Smaller samples, for which
[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html) computes the
exact distribution by default, and calls with additional test arguments
are passed to
[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html) itself, as
are data without variation, whose *p*-value differs between versions of
R.

## Usage

``` r
rank_sum_pvalue(x, g, test_args = list())
```

## Arguments

- x:

  Numeric vector without missing values.

- g:

  Grouping vector of the same length without missing values, with
  exactly two groups.

- test_args:

  Named list of additional arguments for
  [`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html).

## Value

Numeric *p*-value.

## Details

[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html) ranks the
values rounded to the significant digits of its argument `digits.rank`,
which defaults to 7 in recent versions of R and is absent from earlier
ones, which do not round. The values are rounded as the installed
[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html) rounds them
by default, so that the ties, and the *p*-value, are those it finds.
