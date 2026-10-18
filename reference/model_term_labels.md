# Term labels of the fixed effects of a model

Returns the labels of the fixed-effect terms of a model as R writes them
in coefficient names, with non-syntactic variable names enclosed in
backticks (*e.g.,* "\`high bp\`") and calls written out (*e.g.,*
"factor(stage)").

## Usage

``` r
model_term_labels(model)
```

## Arguments

- model:

  Fitted model object.

## Value

Character vector of term labels, or `character(0)` when the terms cannot
be recovered.
