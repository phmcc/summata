# Add the levels of logical predictors to the factor levels of a model

R codes a logical predictor as a factor with the levels `FALSE` and
`TRUE`, naming its coefficient by the second (`olderTRUE`) with `FALSE`
as the reference, but records levels only for factors and character
variables (`xlevels`). With the two levels added, tables and forest
plots show a logical predictor by level, with a reference row for
`FALSE`, as they show a factor, rather than as a single term named by
its coefficient.

## Usage

``` r
add_logical_levels(xlevels, classes)
```

## Arguments

- xlevels:

  Named list of factor levels, or `NULL`.

- classes:

  Named character vector of the classes of the predictors, as in the
  `"dataClasses"` attribute of a terms object.

## Value

`xlevels` with the levels `"FALSE"` and `"TRUE"` for each logical
predictor it does not already contain.
