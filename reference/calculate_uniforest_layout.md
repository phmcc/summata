# Calculate table layout for uniforest plots

Internal function to determine column positions and widths for forest
plot table section. Positions are calculated in the same units as the
data (log scale for OR/HR/RR, linear for coefficients).

## Usage

``` r
calculate_uniforest_layout(
  to_show_exp_clean,
  show_n,
  show_events,
  indent_groups,
  table_width,
  center_padding,
  effect_abbrev,
  font_size,
  log_scale,
  rangeb,
  ci_pct = 95
)
```

## Arguments

- to_show_exp_clean:

  Data.table with formatted data for plotting.

- show_n:

  Logical whether to include n column.

- show_events:

  Logical whether to include events column.

- indent_groups:

  Logical whether levels are indented.

- table_width:

  Proportion of width for table.

- center_padding:

  Padding between table and forest.

- effect_abbrev:

  Effect type abbreviation.

- font_size:

  Font size multiplier.

- log_scale:

  Logical whether using log scale.

- rangeb:

  Numeric vector with plot range bounds (in data units).

## Value

List with column positions, widths, and layout parameters.
