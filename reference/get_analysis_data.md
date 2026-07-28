# Restrict data to the observations used in model fitting

Group-level sample sizes must be counted over the observations that
entered the model, not over every row supplied by the user. Several
fitted objects carry no model frame from which those observations can be
read directly:
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) does
not store one by default, coxme models never store one, and the
memory-conserving fits performed internally by
[`uniscreen()`](https://phmcc.codeberg.page/summata/reference/uniscreen.md)
and
[`multifit()`](https://phmcc.codeberg.page/summata/reference/multifit.md)
pass `model = FALSE`. In these cases the supplied data still contains
the observations dropped during fitting, and counting them inflates
group sizes whenever the response or a covariate is missing.

## Usage

``` r
get_analysis_data(model, model_class, data)
```

## Arguments

- model:

  Fitted model object.

- model_class:

  Character string of the model class.

- data:

  Data frame or data.table supplied to the model call.

## Value

Data.table restricted to the observations used in fitting. The input is
returned unchanged when the analysis rows cannot be determined.

## Details

Complete cases are determined from the variables named in the model
formula, reproducing the default `na.omit` behavior of the modeling
functions. A complete-case restriction can only remove observations that
the model itself could not have used, so the subset is retained even
when it does not reproduce the fitted sample size exactly, as may happen
when a fit used `subset` or case weights.
