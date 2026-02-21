# Perform statistical tests for categorical variables

Conducts chi-square or Fisher's exact tests for categorical variables
across groups. Automatically selects Fisher's exact test for small
expected frequencies.

## Usage

``` r
perform_categorical_test(tab, test_type)
```

## Arguments

- tab:

  Contingency table (matrix or table object).

- test_type:

  Character string: "chisq" for chi-square, "fisher" for Fisher's exact,
  or "auto" for automatic selection.

## Value

Numeric *p*-value from the hypothesis test.
