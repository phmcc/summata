# Determine alignment for exported tables

Creates column alignment string for LaTeX tables. Variable and Group
columns are left-aligned; all others are centered.

## Usage

``` r
determine_alignment(df)
```

## Arguments

- df:

  Data.frame or data.table to determine alignment for.

## Value

Character string with alignment codes (*e.g.,* "rlcc").
