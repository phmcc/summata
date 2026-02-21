# Combine coefficient tables from multiple models

Merges coefficient tables from multiple fitted models into a single
data.table with a Model identifier column.

## Usage

``` r
combine_coefficient_tables(coef_list, model_names)
```

## Arguments

- coef_list:

  List of data.tables containing coefficient information.

- model_names:

  Character vector of model names corresponding to coef_list.

## Value

Combined data.table with Model column, or `NULL` if empty.
