# Return model-fitting errors as values

Wraps a function of one argument so that an error is returned as an
object of class `"summata_fit_error"` holding the condition message
instead of being signaled. Sequential and parallel runs then deliver
every failure in the same form to
[`stop_on_failed_fits()`](https://phmcc.codefloe.page/summata/reference/stop_on_failed_fits.md).

## Usage

``` r
capture_fit_error(fun)
```

## Arguments

- fun:

  Function of one argument.

## Value

Function of one argument.
