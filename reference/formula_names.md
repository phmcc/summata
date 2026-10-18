# Write variable names as a model formula reads them

The variables of
[`fit()`](https://phmcc.codefloe.page/summata/reference/fit.md),
[`uniscreen()`](https://phmcc.codefloe.page/summata/reference/uniscreen.md),
[`multifit()`](https://phmcc.codefloe.page/summata/reference/multifit.md),
and
[`fullfit()`](https://phmcc.codefloe.page/summata/reference/fullfit.md)
are given as column names, from which the model formula is written as
text. A column name that is not a syntactic name, such as
`"tumor stage"`, is read in a formula as code unless it is quoted with
backticks. A term that is a column name, each component of an
interaction term (`"a:b"`) whose components are all column names, and
each argument of a
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) outcome that is a
column name are quoted where the name requires it. A column name may
also be given with backticks. Other terms, such as function calls and
random effects, are kept as given, as is a term in which no name
changes.

## Usage

``` r
formula_names(terms, var_names, quote = TRUE)
```

## Arguments

- terms:

  Character vector of terms, or `NULL`.

- var_names:

  Character vector of column names.

- quote:

  Logical. If `TRUE` (default), names are written as a formula reads
  them. If `FALSE`, the backticks of column names are removed instead,
  so that the terms match the columns of the data;
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) outcomes are
  then kept as given.

## Value

Character vector of terms, or `NULL` when `terms` is `NULL`.
