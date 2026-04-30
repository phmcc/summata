# Determine CI separator for forest plot text annotations

Returns the appropriate separator string between CI lower and upper
bounds in forest plot annotations. Considers whether values may be
negative and the current locale's decimal mark.

## Usage

``` r
forest_ci_separator(has_negatives, marks = NULL)
```

## Arguments

- has_negatives:

  Logical indicating whether any CI bounds are negative.

- marks:

  Optional list with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.codeberg.page/summata/reference/resolve_number_marks.md).

## Value

Character string separator.
