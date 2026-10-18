# Perform statistical tests for categorical variables

Conducts chi-square or Fisher's exact tests for categorical variables
across groups. Automatically selects Fisher's exact test for small
expected frequencies. The test compares the observed levels of the
variable across the defined groups. The missing-value row and column
that `table(useNA = "ifany")` adds are removed, as are levels and groups
ithout observations: an empty row or column contributes nothing to
either test, but its zero expected counts would otherwise select the
exact test under `"auto"` and make the chi-square statistic undefined.

## Usage

``` r
perform_categorical_test(tab, test_type, test_args = list(), var_name = NULL)
```

## Arguments

- tab:

  Contingency table (matrix or table object).

- test_type:

  Character string: "chisq" for chi-square, "fisher" for Fisher's exact,
  or "auto" for automatic selection.

- test_args:

  Named list of additional test arguments from `...`; each test receives
  only the arguments it accepts (see
  [`filter_test_args()`](https://phmcc.codefloe.page/summata/reference/filter_test_args.md)).

- var_name:

  Character string naming the variable, used in warnings.

## Value

Numeric *p*-value from the hypothesis test.
