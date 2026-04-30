# Create Forest Plot for Multivariate Regression

Generates a publication-ready forest plot from a
[`multifit()`](https://phmcc.codeberg.page/summata/reference/multifit.md)
output object. The plot displays effect estimates (OR, HR, RR, or
coefficients) with confidence intervals across multiple outcomes,
organized by outcome with the predictor levels shown for each.

## Usage

``` r
multiforest(
  x,
  title = "Multivariate Analysis",
  effect_label = NULL,
  column = "adjusted",
  digits = 2,
  p_digits = 3,
  conf_level = 0.95,
  font_size = 1,
  annot_size = 3.88,
  header_size = 5.82,
  title_size = 23.28,
  plot_width = NULL,
  plot_height = NULL,
  table_width = 0.6,
  show_n = TRUE,
  show_events = NULL,
  show_predictor = NULL,
  covariates_footer = TRUE,
  indent_predictor = FALSE,
  bold_variables = TRUE,
  center_padding = 4,
  zebra_stripes = TRUE,
  color = NULL,
  null_line = NULL,
  log_scale = NULL,
  labels = NULL,
  units = "in",
  number_format = NULL
)
```

## Arguments

- x:

  Multifit result object (data.table with class attributes from
  [`multifit()`](https://phmcc.codeberg.page/summata/reference/multifit.md)).

- title:

  Character string specifying the plot title. Default is
  `"Multivariate Analysis"`. Use descriptive titles for publication.

- effect_label:

  Character string for the effect measure label on the forest plot axis.
  Default is `NULL`, which auto-detects based on model type (*e.g.,*
  "Odds Ratio", "Hazard Ratio", "Rate Ratio", "Estimate").

- column:

  Character string specifying which results to plot when
  [`multifit()`](https://phmcc.codeberg.page/summata/reference/multifit.md)
  was called with `columns = "both"`. Options are `"adjusted"` (default)
  or `"unadjusted"`. Ignored when the `multifit` result contains only
  one column type.

- digits:

  Integer specifying the number of decimal places for effect estimates
  and confidence intervals. Default is 2.

- p_digits:

  Integer specifying the number of decimal places for *p*-values. Values
  smaller than `10^(-p_digits)` are displayed as `"< 0.001"` (for
  `p_digits = 3`), `"< 0.0001"` (for `p_digits = 4`), etc. Default is 3.

- conf_level:

  Numeric confidence level for confidence intervals. Must be between 0
  and 1. Default is 0.95 (95% confidence intervals). The CI percentage
  is automatically displayed in column headers (*e.g.,* "90% CI" when
  `conf_level = 0.90`).

- font_size:

  Numeric multiplier controlling the base font size for all text
  elements. Default is 1.0.

- annot_size:

  Numeric value controlling the relative font size for data annotations.
  Default is 3.88.

- header_size:

  Numeric value controlling the relative font size for column headers.
  Default is 5.82.

- title_size:

  Numeric value controlling the relative font size for the main plot
  title. Default is 23.28.

- plot_width:

  Numeric value specifying the intended output width in specified
  `units`. Used for optimizing layout. Default is `NULL` (automatic).
  Recommended: 10-16 inches.

- plot_height:

  Numeric value specifying the intended output height in specified
  `units`. Default is `NULL` (automatic based on number of rows).

- table_width:

  Numeric value between 0 and 1 specifying the proportion of total plot
  width allocated to the data table. Default is 0.6 (60% table, 40%
  forest plot).

- show_n:

  Logical. If `TRUE`, includes a column showing sample sizes. Default is
  `TRUE`.

- show_events:

  Logical. If `TRUE`, includes a column showing the number of events for
  each row. Default is `TRUE` for binomial and survival models, `FALSE`
  for linear models.

- show_predictor:

  Logical. If `TRUE`, includes the "Predictor" column showing which
  level of a factor predictor is being compared. If `FALSE`, omits the
  column (useful when predictor info is in the caption). Default is
  `NULL`, which uses the `include_predictor` setting from
  [`multifit()`](https://phmcc.codeberg.page/summata/reference/multifit.md)
  if available, otherwise `TRUE`.

- covariates_footer:

  Logical. If `TRUE` (default), displays a footer listing the covariates
  used in adjusted models. Covariate names are formatted using the
  `labels` parameter if provided. Only shown when displaying adjusted
  results.

- indent_predictor:

  Logical. If `TRUE`, indents predictor levels under outcome names,
  creating a hierarchical display. If `FALSE` (default), shows outcome
  and predictor level in separate columns.

- bold_variables:

  Logical. If `TRUE`, variable names are displayed in bold. If `FALSE`
  (default), variable names are displayed in plain text.

- center_padding:

  Numeric value specifying horizontal spacing between table and forest
  plot. Default is 4.

- zebra_stripes:

  Logical. If `TRUE`, applies alternating gray background shading to
  different outcomes for improved readability. Default is `TRUE`.

- color:

  Character string specifying the color for point estimates in the
  forest plot. Default is `NULL`, which auto-selects based on effect
  type: purple (`"#8A61D8"`) for hazard ratios (Cox), teal (`"#4BA6B6"`)
  for odds ratios (logistic), blue (`"#3F87EE"`) for rate/risk ratios
  (Poisson, Gamma, *etc.*), and green (`"#5A8F5A"`) for coefficients
  (linear models). Use hex codes or R color names for custom colors.

- null_line:

  Numeric value for the reference line position. Default is `NULL`,
  which uses 1 for ratio measures (OR, HR, RR) and 0 for coefficients.

- log_scale:

  Logical. If `TRUE`, uses log scale for the x-axis. Default is `NULL`,
  which auto-detects (`TRUE` for OR/HR/RR, `FALSE` for coefficients).

- labels:

  Named character vector providing custom display labels for outcomes
  and variables. Applied to outcome names in the plot. Default is `NULL`
  (uses labels already applied in multifit, or original names).

- units:

  Character string specifying units for plot dimensions: `"in"`
  (inches), `"cm"`, or `"mm"`. Default is `"in"`.

- number_format:

  Character string or two-element character vector controlling thousand
  and decimal separators in formatted output. Named presets:

  - `"us"` - Comma thousands, period decimal: `1,234.56` \[default\]

  - `"eu"` - Period thousands, comma decimal: `1.234,56`

  - `"space"` - Thin-space thousands, period decimal: `1 234.56` (SI/ISO
    31-0)

  - `"none"` - No thousands separator: `1234.56`

  Or provide a custom two-element vector `c(big.mark, decimal.mark)`,
  *e.g.*, `c("'", ".")` for Swiss-style: `1'234.56`.

  When `NULL` (default), uses
  `getOption("summata.number_format", "us")`. Set the global option once
  per session to avoid passing this argument repeatedly:

          options(summata.number_format = "eu")
        

## Value

A `ggplot` object containing the complete forest plot. The plot can be:

- Displayed directly: `print(plot)`

- Saved to file: `ggsave("forest.pdf", plot, width = 12, height = 8)`

- Further customized with ggplot2 functions

The returned object includes an attribute `"rec_dims"` accessible via
`attr(plot, "rec_dims")`, which is a list containing:

- width:

  Numeric. Recommended plot width in specified units

- height:

  Numeric. Recommended plot height in specified units

These recommendations are automatically calculated based on the number
of variables, text sizes, and layout parameters, and are printed to
console if `plot_width` or `plot_height` are not specified.

## Details

**Plot Layout:**

The forest plot is organized with outcomes as grouping headers and
predictor levels (or interaction terms) as rows within each outcome.
This provides a clear visual comparison of how a single predictor
affects multiple outcomes.

1.  **Title**: Centered at top

2.  **Data Table** (left): Contains:

    - Outcome column (or grouped headers)

    - Predictor/Group column

    - n: Sample sizes (optional)

    - Events: Event counts (optional, for applicable models)

    - Effect (95% CI); *p*-value

3.  **Forest Plot** (right):

    - Point estimates (squares)

    - 95% confidence intervals

    - Reference line at null value (1 or 0)

    - Log scale for ratio measures

**Data Source:**

The function extracts effect estimates directly from the multifit output
object's `raw_data` attribute, which contains the numeric values needed
for plotting. This approach is efficient and ensures consistency with
the formatted table output.

## See also

[`autoforest`](https://phmcc.codeberg.page/summata/reference/autoforest.md)
for automatic model detection,
[`multifit`](https://phmcc.codeberg.page/summata/reference/multifit.md)
for multi-outcome regression analysis,
[`glmforest`](https://phmcc.codeberg.page/summata/reference/glmforest.md)
for single GLM forest plots,
[`coxforest`](https://phmcc.codeberg.page/summata/reference/coxforest.md)
for single Cox model forest plots,
[`lmforest`](https://phmcc.codeberg.page/summata/reference/lmforest.md)
for single linear model forest plots,
[`uniforest`](https://phmcc.codeberg.page/summata/reference/uniforest.md)
for univariable screening forest plots

Other visualization functions:
[`autoforest()`](https://phmcc.codeberg.page/summata/reference/autoforest.md),
[`coxforest()`](https://phmcc.codeberg.page/summata/reference/coxforest.md),
[`glmforest()`](https://phmcc.codeberg.page/summata/reference/glmforest.md),
[`lmforest()`](https://phmcc.codeberg.page/summata/reference/lmforest.md),
[`uniforest()`](https://phmcc.codeberg.page/summata/reference/uniforest.md)

## Examples

``` r
data(clintrial)
data(clintrial_labels)
library(survival)

# Create example multifit result
result <- multifit(
    data = clintrial,
    outcomes = c("surgery", "pfs_status", "os_status"),
    predictor = "treatment",
    covariates = c("age", "sex", "stage"),
    parallel = FALSE
)

# Example 1: Basic multivariate forest plot
p <- multiforest(result)
#> Recommended plot dimensions: width = 14.6 in, height = 5.0 in

# \donttest{

old_width <- options(width = 180)

# Example 2: With custom title and labels
plot2 <- multiforest(
    result,
    title = "Treatment Effects Across Clinical Outcomes",
    labels = clintrial_labels
)
#> Recommended plot dimensions: width = 16.7 in, height = 5.0 in

# Example 3: Customize appearance
plot3 <- multiforest(
    result,
    color = "#E74C3C",
    zebra_stripes = TRUE,
    labels = clintrial_labels
)
#> Recommended plot dimensions: width = 16.7 in, height = 5.0 in

# Example 4: Save with recommended dimensions
dims <- attr(plot3, "rec_dims")
ggplot2::ggsave(file.path(tempdir(), "multioutcome_forest.pdf"),
                plot3, width = dims$width, height = dims$height)

options(old_width)

# }
```
