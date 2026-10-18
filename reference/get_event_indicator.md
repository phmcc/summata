# Code the events of a survival model

Evaluates the event argument of the model's
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) response (see
[`surv_event_expression()`](https://phmcc.codefloe.page/summata/reference/surv_event_expression.md))
on the rows of `data` and codes it as
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) codes it: 1 for
an event and 0 for a censored observation, whether the status is coded 0
and 1, `FALSE` and `TRUE`, or 1 and 2, or is given as an expression such
as `status == 2`. A factor is coded by its second level. Summing the
variable named in the response instead counted a status coded 1 and 2
twice for each event and once for each censored observation, and reading
the second variable of the formula took the stop time of a
counting-process response and the first predictor of a
[`clogit()`](https://rdrr.io/pkg/survival/man/clogit.html) model.

## Usage

``` r
get_event_indicator(model, data)
```

## Arguments

- model:

  Fitted `coxph`, `clogit`, or `coxme` model.

- data:

  Data frame in which the event argument is evaluated, such as the
  observations analyzed or the observations supplied.

## Value

Numeric vector of 0 and 1 (missing where the event is), one per row of
`data`, or `NULL` when the events cannot be read.
