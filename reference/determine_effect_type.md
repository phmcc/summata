# Determine Effect Type Label

Identifies the appropriate effect measure label (OR, HR, RR,
Coefficient, aOR, aHR, aRR, Adj. Coefficient) based on model type,
exponentiation setting, and whether the estimate is adjusted
(multivariable) or unadjusted (univariable).

## Usage

``` r
determine_effect_type(uni_raw, multi_raw, exponentiate, adjusted = FALSE)
```

## Arguments

- uni_raw:

  Raw univariable data.table containing coefficient columns. Used to
  detect effect type when `adjusted = FALSE`.

- multi_raw:

  Raw multivariable data.table containing coefficient columns. Used to
  detect effect type when `adjusted = TRUE`.

- exponentiate:

  Logical or `NULL` controlling label selection:

  `TRUE`

  :   Force exponentiated labels (OR, HR, RR, Exp(Coef))

  `FALSE`

  :   Force coefficient labels (Coefficient, Adj. Coefficient)

  `NULL`

  :   Auto-detect from column names in raw data

- adjusted:

  Logical. If `TRUE`, returns adjusted effect labels (aOR, aHR, aRR,
  Adj. Coefficient) for multivariable results. If `FALSE`, returns
  unadjusted labels (OR, HR, RR, Coefficient) for univariable results.

## Value

Character string with the effect measure label:
