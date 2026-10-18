# Coefficients that belong to a variable whose name extends a factor's name

Identifies the coefficients whose names begin with the name of a factor
but belong to another variable of the model whose name begins with it as
well, such as the continuous "smoking_years" or the coefficient
"smokerTRUE" of a logical variable "smoker" beside a factor "smoking". A
coefficient whose remainder after the factor's name is one of the
factor's levels is kept for the factor.

## Usage

``` r
claimed_by_longer_term(terms, var, var_levels, labels)
```

## Arguments

- terms:

  Character vector of coefficient names, without backticks.

- var:

  Character string. Name of the factor.

- var_levels:

  Character vector of the levels of the factor.

- labels:

  Character vector of the term labels of the other variables of the
  model, without backticks.

## Value

Logical vector, `TRUE` for the coefficients that belong to a variable
with a longer name.
