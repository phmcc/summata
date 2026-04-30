# Calculate table layout for forest plots

Computes column widths and positions for the table portion of a forest
plot. Determines spacing based on content width, font size, and desired
table/forest proportion. Returns positions in log-scale units for plot
coordinate system.

## Usage

``` r
calculate_forest_layout(
  to_show_exp_clean,
  show_n,
  show_events,
  indent_groups,
  condense_table,
  effect_label,
  ref_label,
  font_size,
  table_width = 0.6,
  rangeb,
  center_padding
)
```

## Arguments

- to_show_exp_clean:

  Data.table with formatted display data for the plot.

- show_n:

  Logical whether to include sample size column.

- show_events:

  Logical whether to include events column.

- indent_groups:

  Logical whether groups are indented (affects level column).

- condense_table:

  Logical whether table is condensed (affects level column).

- effect_label:

  Character string describing effect measure type.

- ref_label:

  Character string label for reference categories.

- font_size:

  Numeric font size for width calculations.

- table_width:

  Numeric proportion of total width for table (0-1).

- rangeb:

  Numeric vector of length 2 with plot x-axis range.

- center_padding:

  Numeric additional padding for effect column.

## Value

List with table_width, forest_width, positions, rangeplot_start,
total_width, and effect_abbrev components.
