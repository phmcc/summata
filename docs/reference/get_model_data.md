# Get data from model object (works with S3 and S4)

Retrieves the original data used to fit a model. Checks multiple
locations including model attributes, \$data slot, \$model slot, and
@frame for S4.

## Usage

``` r
get_model_data(model)
```

## Arguments

- model:

  Fitted model object (S3 or S4).

## Value

Data frame or data.table used to fit the model, or `NULL` if
unavailable.
