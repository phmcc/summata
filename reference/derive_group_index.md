# Locate the observations of each group

Finds the row positions of each group once, so that the summaries of
every variable can subset them instead of comparing the full grouping
vector with each group again.

## Usage

``` r
derive_group_index(grp_vec, groups)
```

## Arguments

- grp_vec:

  Grouping vector.

- groups:

  Vector of groups, as returned by
  [`derive_group_levels()`](https://phmcc.codefloe.page/summata/reference/derive_group_levels.md).

## Value

List with one integer vector of row positions per group, in the order of
`groups`. Observations with a missing group appear in none.
