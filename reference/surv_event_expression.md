# Read the event argument of a `Surv()` response

The event of a [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
response is its event argument: the second argument of
`Surv(time, event)`, the third of `Surv(time, time2, event)`, or the
argument named `event`.
[`survival::clogit()`](https://rdrr.io/pkg/survival/man/clogit.html)
writes its response as `Surv(rep(1, N), case)`, whose event is `case`.
The argument may be an expression, such as `status == 2`.

## Usage

``` r
surv_event_expression(response)
```

## Arguments

- response:

  Left-hand side of a model formula, as a call, or a character string
  such as `"Surv(os_months, os_status)"`.

## Value

The expression of the event argument, or `NULL` when the response is not
a call to [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) with an
event.
