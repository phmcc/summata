# Fit a model with selective warning suppression

Wraps model fitting expressions to suppress routine warnings from
mixed-effects and GLM fitting (*e.g.,* singular fits, convergence
messages, separation warnings) while allowing unexpected warnings
through. When `verbose = TRUE`, all warnings are displayed.

## Usage

``` r
quiet_fit(expr, verbose = FALSE)
```

## Arguments

- expr:

  An unevaluated expression (model fitting call) to execute.

- verbose:

  Logical. If `TRUE`, all warnings are shown. If `FALSE`, routine
  fitting warnings are suppressed. Default is `FALSE`.

## Value

The result of evaluating `expr`.
