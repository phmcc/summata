# Finalize Column Names for Display

Renames internal column names (uni_effect, uni_p, multi_effect, multi_p)
to publication-ready display names with appropriate effect measure
labels (OR, HR, RR, aOR, aHR, aRR, Coefficient, *etc.*).

## Usage

``` r
finalize_column_names(
  result,
  uni_raw,
  multi_raw,
  exponentiate,
  columns,
  metrics,
  conf_level = 0.95
)
```

## Arguments

- result:

  Data.table with columns to rename. Expected to contain some
  combination of: uni_effect, uni_p, multi_effect, multi_p.

- uni_raw:

  Raw univariable data.table used to determine effect type by checking
  for presence of OR, HR, RR, or Coefficient columns.

- multi_raw:

  Raw multivariable data.table used to determine adjusted effect type.

- exponentiate:

  Logical or `NULL`. If `TRUE`, forces exponentiated labels (OR, HR,
  RR). If `FALSE`, uses "Coefficient". If `NULL`, auto-detects from raw
  data columns.

- columns:

  Character string indicating which columns are present (`"both"`,
  `"uni"`, or `"multi"`). Used for context but renaming is based on
  actual column presence.

- metrics:

  Character vector of metrics being displayed (`"effect"`, `"p"`, or
  both). Used for context.

## Value

The input data.table with columns renamed
