# Average ranks with tie counts

Ranks a numeric vector with ties assigned their average rank, as
[`rank()`](https://rdrr.io/r/base/rank.html) does by default, from a
single sort with
[`data.table::frank()`](https://rdrr.io/pkg/data.table/man/frank.html).
The \\j\\th smallest distinct value, shared by \\t_j\\ observations,
occupies sorted positions \\c\_{j-1} + 1, \dots, c_j\\, where \\c\\ is
the cumulative sum of the tie counts, so its average rank is \\c_j -
(t_j - 1)/2\\. Every quantity is an integer or half-integer, so the
ranks are exact.

## Usage

``` r
average_ranks(x)
```

## Arguments

- x:

  Numeric vector without missing values.

## Value

Numeric vector of ranks, with attribute `"ties"` holding the number of
observations at each distinct value in increasing order.
