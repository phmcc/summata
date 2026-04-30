# Print method for fit results

Displays a summary header with model scope (Univariable/Multivariable),
model type, formula, sample size, and event count before printing the
formatted results table.

## Usage

``` r
# S3 method for class 'fit_result'
print(x, ...)
```

## Arguments

- x:

  fit_result object.

- ...:

  Additional arguments passed to print methods.

## Value

Invisibly returns the input object `x`. Called for its side effect of
printing a formatted summary to the console.

## See also

Other regression functions:
[`compfit()`](https://phmcc.github.io/summata/reference/compfit.md),
[`fit()`](https://phmcc.github.io/summata/reference/fit.md),
[`fullfit()`](https://phmcc.github.io/summata/reference/fullfit.md),
[`multifit()`](https://phmcc.github.io/summata/reference/multifit.md),
[`print.compfit_result()`](https://phmcc.github.io/summata/reference/print.compfit_result.md),
[`print.fullfit_result()`](https://phmcc.github.io/summata/reference/print.fullfit_result.md),
[`print.multifit_result()`](https://phmcc.github.io/summata/reference/print.multifit_result.md),
[`print.uniscreen_result()`](https://phmcc.github.io/summata/reference/print.uniscreen_result.md),
[`uniscreen()`](https://phmcc.github.io/summata/reference/uniscreen.md)
