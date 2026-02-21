# Format interaction term for display

Converts R's internal interaction term format (*e.g.,* "treatmentDrug
A:stageII") to a more readable format (*e.g.,* "Treatment (Drug A) ×
Stage (II)").

## Usage

``` r
format_interaction_term(term, labels = NULL)
```

## Arguments

- term:

  Character string of the interaction term from model coefficients.

- labels:

  Optional named vector of labels for variable names.

## Value

Formatted interaction term string.
