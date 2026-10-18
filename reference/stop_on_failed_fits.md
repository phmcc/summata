# Stop when any model could not be fitted

Checks the results of a run over predictors or outcomes and stops with
an error that names every element whose fit failed, with its message. An
element has failed when its fit raised an error (captured by
[`capture_fit_error()`](https://phmcc.codefloe.page/summata/reference/capture_fit_error.md)),
when a forked worker failed outside the fitting function
([`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
returns a `"try-error"` object), or when a worker returned no result.

## Usage

``` r
stop_on_failed_fits(results, inputs, what)
```

## Arguments

- results:

  List of results, one per element of `inputs`.

- inputs:

  Character vector naming the elements.

- what:

  Character string naming one element in the message, such as
  `"predictor"` or `"outcome"`.

## Value

`results`, invisibly, when every element was fitted.
