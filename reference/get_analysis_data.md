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
