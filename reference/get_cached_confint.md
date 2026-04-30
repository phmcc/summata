# Retrieve confidence intervals with cache support

Returns confidence intervals for a fitted model, using a cached result
when available from upstream
[`fit()`](https://phmcc.codeberg.page/summata/reference/fit.md) or
[`m2dt()`](https://phmcc.codeberg.page/summata/reference/m2dt.md) calls
to avoid redundant computation. This is particularly beneficial for GLM
and negative binomial models where
[`confint()`](https://rdrr.io/r/stats/confint.html) performs profile
likelihood profiling, which can be expensive for models with many
parameters.

## Usage

``` r
get_cached_confint(model, conf_level = 0.95)
```

## Arguments

- model:

  Fitted model object. If the model carries a `"cached_confint"`
  attribute (set by
  [`fit()`](https://phmcc.codeberg.page/summata/reference/fit.md) during
  table generation) and the cached confidence level matches
  `conf_level`, the cached result is returned directly.

- conf_level:

  Numeric confidence level. Must match the cached level for the cache to
  be used.

## Value

A matrix with one row per model coefficient and two columns (lower and
upper bounds), as returned by
[`stats::confint()`](https://rdrr.io/r/stats/confint.html).
