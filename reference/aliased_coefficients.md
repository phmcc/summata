# Names of the coefficients that a model could not estimate

Returns the names of the aliased coefficients of a model, those whose
columns of the model matrix are linear combinations of other columns.
[`stats::coef()`](https://rdrr.io/r/stats/coef.html) reports them as
`NA` for `lm`, `glm`, negative binomial and `coxph` models; lme4 drops
their columns from the model matrix before fitting and records them in
its `"col.dropped"` attribute. The backticks that enclose non-syntactic
variable names in coefficient names are removed.

## Usage

``` r
aliased_coefficients(model)
```

## Arguments

- model:

  Fitted model object.

## Value

Character vector of coefficient names, empty when every coefficient was
estimated.
