# Parse term into variable and group

Splits coefficient term names into base variable names and factor
levels. For example, "sexMale" becomes variable="sex" and group="Male".
Handles interaction terms and continuous variables appropriately.

## Usage

``` r
parse_term(terms, xlevels = NULL, model = NULL)
```

## Arguments

- terms:

  Character vector of coefficient term names.

- xlevels:

  Named list of factor levels from the model.

- model:

  Optional model object, from which the terms of the model are read and
  the factor levels of coxme models are reconstructed.

## Value

Data.table with 'variable' and 'group' columns.

## Details

A term that belongs to another variable of the model whose name begins
with a factor's name, such as a continuous variable "smoking_years" or
the coefficient "smokerTRUE" of a logical variable "smoker" beside a
factor "smoking", is not taken for a level of that factor (see
[`claimed_by_longer_term()`](https://phmcc.codefloe.page/summata/reference/claimed_by_longer_term.md)).
