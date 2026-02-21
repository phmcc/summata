# Apply zebra stripes with proper variable group detection for indented tables

Applies alternating background colors to variable groups in flextable
objects. Handles both indented tables (detects groups by leading
whitespace) and non-indented tables (uses pre-identified groups).

## Usage

``` r
apply_zebra_stripes_ft(ft, df, var_groups)
```

## Arguments

- ft:

  flextable object.

- df:

  The source data.table used to create the flextable.

- var_groups:

  List of row index vectors for variable groups.

## Value

Flextable object with zebra stripe formatting applied.
