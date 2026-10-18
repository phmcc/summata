# Confidence intervals for model coefficients

Computes confidence intervals for the named coefficients by the method
that `conf_method` selects. With `"profile"`,
[`stats::confint()`](https://rdrr.io/r/stats/confint.html) gives profile
likelihood intervals for GLM and negative binomial models, exact *t*
intervals for linear models, and the model's own intervals for other
classes (Wald intervals for Cox models); if it fails, as profiling can
under separation, Wald intervals are returned instead. With `"wald"`,
and for quasi-likelihood families, which have no likelihood to profile,
Wald intervals are returned.
[`m2dt()`](https://phmcc.codefloe.page/summata/reference/m2dt.md) and
the forest plot functions share this function so that a table and a plot
of the same model show the same intervals.

## Usage

``` r
coef_confint(model, parm, conf_level, conf_method)
```

## Arguments

- model:

  Fitted model object.

- parm:

  Character vector of coefficient names.

- conf_level:

  Numeric confidence level.

- conf_method:

  Character string, `"profile"` or `"wald"`.

## Value

A matrix with one row per element of `parm` and two columns (lower and
upper bounds).
