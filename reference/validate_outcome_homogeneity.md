# Validate outcome homogeneity for multifit

Checks whether all outcomes in a multifit analysis are compatible with
the specified model type. Issues a warning when outcomes appear to be of
mixed types (*e.g.,* binary and continuous outcomes in the same
analysis), which would produce tables with incompatible effect measures.

## Usage

``` r
validate_outcome_homogeneity(data, outcomes, model_type, family = "binomial")
```

## Arguments

- data:

  Data.table containing the analysis data.

- outcomes:

  Character vector of outcome variable names.

- model_type:

  Character string specifying the model type.

- family:

  Character string specifying the GLM family (for glm/glmer).

## Value

Invisible `NULL`. Issues warnings if problems are detected.
