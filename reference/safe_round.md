# Safe rounding that handles `NULL` and NA

Rounds numeric values while gracefully handling `NULL`, empty, and NA
inputs.

## Usage

``` r
safe_round(x, digits)
```

## Arguments

- x:

  Numeric value to round.

- digits:

  Integer number of decimal places.

## Value

Rounded numeric value, or NA_real\_ if input is `NULL`/NA/empty.
