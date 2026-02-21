# Get factor levels from model (works with S3 and S4)

Extracts factor level information from fitted model objects. Handles
both S3 models (glm, lm, coxph) via xlevels slot and S4 models (lme4)
via the model frame.

## Usage

``` r
get_model_xlevels(model)
```

## Arguments

- model:

  Fitted model object (S3 or S4).

## Value

Named list of factor levels, or `NULL` if no factors present.
