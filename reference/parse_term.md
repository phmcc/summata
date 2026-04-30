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

  Optional model object for extracting factor info from coxme models.

## Value

Data.table with 'variable' and 'group' columns.
