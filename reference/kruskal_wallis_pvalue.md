# Kruskal-Wallis rank sum test p-value

Computes the *p*-value of the Kruskal-Wallis rank sum test with the
statistic, tie correction, and degrees of freedom of
[`stats::kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html),
using ranks and tie counts from
[`average_ranks()`](https://phmcc.codefloe.page/summata/reference/average_ranks.md).
[`kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html) ranks with
[`rank()`](https://rdrr.io/r/base/rank.html) and counts ties with
[`table()`](https://rdrr.io/r/base/table.html), which converts every
value to a character string; on large samples these two steps account
for nearly all of its running time.

## Usage

``` r
kruskal_wallis_pvalue(x, g)
```

## Arguments

- x:

  Numeric vector without missing values.

- g:

  Grouping vector of the same length without missing values.

## Value

Numeric *p*-value.
