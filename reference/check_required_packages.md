# Check required packages for model type

Verifies that necessary packages are installed for the specified model
type. Stops with informative error if required packages are missing.

## Usage

``` r
check_required_packages(model_type)
```

## Arguments

- model_type:

  Character string indicating model type.

## Value

`NULL` (invisibly). Stops execution if packages missing.
