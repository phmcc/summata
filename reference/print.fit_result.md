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
[`compfit()`](https://phmcc.codeberg.page/summata/reference/compfit.md),
[`fit()`](https://phmcc.codeberg.page/summata/reference/fit.md),
[`fullfit()`](https://phmcc.codeberg.page/summata/reference/fullfit.md),
[`multifit()`](https://phmcc.codeberg.page/summata/reference/multifit.md),
[`print.compfit_result()`](https://phmcc.codeberg.page/summata/reference/print.compfit_result.md),
[`print.fullfit_result()`](https://phmcc.codeberg.page/summata/reference/print.fullfit_result.md),
[`print.multifit_result()`](https://phmcc.codeberg.page/summata/reference/print.multifit_result.md),
[`print.uniscreen_result()`](https://phmcc.codeberg.page/summata/reference/print.uniscreen_result.md),
[`uniscreen()`](https://phmcc.codeberg.page/summata/reference/uniscreen.md)
