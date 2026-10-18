# Code the events of a `Surv()` response

Code the events of a
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) response

## Usage

``` r
surv_event_indicator(response, data, env = parent.frame())
```

## Arguments

- response:

  Left-hand side of a model formula, as a call, or a character string
  such as `"Surv(os_months, os_status)"`.

- data:

  Data frame in which the event argument is evaluated.

- env:

  Environment in which variables not in `data` are found.

## Value

See
[`get_event_indicator()`](https://phmcc.codefloe.page/summata/reference/get_event_indicator.md).
