# Detect outcome type from data

Automatically determines whether an outcome variable is binary,
continuous, or count-based by examining the data values. Used for
automatic model type selection and validation. Binary outcomes have 2
unique values, continuous have many values or non-integers, counts have
integers \>= 0.

## Usage

``` r
detect_outcome_type(data, outcome)
```

## Arguments

- data:

  Data frame or data.table containing the outcome variable.

- outcome:

  Character string naming the outcome variable.

## Value

Character string: "binary", "continuous", "count", or "unknown".
