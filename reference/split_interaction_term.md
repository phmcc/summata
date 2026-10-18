# Split an interaction coefficient into variables and levels

Splits the name of an interaction coefficient (*e.g.,* "treatmentDrug
A:stageII") at each colon and matches each component to the variable it
names. A component of a factor is the variable name followed by one of
its levels; a component of a continuous variable is the variable name
alone.

## Usage

``` r
split_interaction_term(term, var_names, xlevels = NULL)
```

## Arguments

- term:

  Character string. Name of an interaction coefficient.

- var_names:

  Character vector. Names of the variables to match.

- xlevels:

  Named list of factor levels, or `NULL`.

## Value

List with one element per component, each a list with elements `part`
(the component), `var` (the matched variable, or `NA`) and `level` (the
level, or `NA` for a continuous variable or an unmatched component).

## Details

The backticks that enclose non-syntactic variable names in coefficient
names (*e.g.,* "\`high bp\`Yes") are removed first. Names are compared
literally, longest first, and a variable is preferred whose levels
include the remainder of the component, so that "sexualityHetero" is
assigned to `sexuality` rather than to `sex`. A level that contains a
colon (*e.g.,* "08:00") is kept in one component when it is among the
levels in `xlevels`.
