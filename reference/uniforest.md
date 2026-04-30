# Create Forest Plot for Univariable Screening

Generates a publication-ready forest plot from a
[`uniscreen()`](https://phmcc.codeberg.page/summata/reference/uniscreen.md)
output object. The plot displays effect estimates (OR, HR, RR, or
coefficients) with confidence intervals for each predictor tested in
univariable analysis against a single outcome.

## Usage

``` r
uniforest(
  x,
  title = "Univariable Screening",
  effect_label = NULL,
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
  indent_groups = FALSE,
  condense_table = FALSE,
  bold_variables = FALSE,
  center_padding = 4,
  zebra_stripes = TRUE,
  color = NULL,
  null_line = NULL,
  log_scale = NULL,
  labels = NULL,
  show_footer = TRUE,
  units = "in",
  number_format = NULL
)
```

## Arguments

- x:

  Univariable screen result object (data.table with class attributes
  from
  [`uniscreen()`](https://phmcc.codeberg.page/summata/reference/uniscreen.md)).

- title:

  Character string specifying the plot title. Default is
  `"Univariable Screening"`. Use descriptive titles for publication.

- effect_label:

  Character string for the effect measure label on the forest plot axis.
  Default is `NULL`, which auto-detects based on model type (*e.g.,*
  "Odds Ratio", "Hazard Ratio", "Rate Ratio", "Coefficient").

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
  `conf_level = 0.90`). Note: This parameter affects display only; the
  underlying CIs come from the uniscreen result.

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
  each row. Default is `NULL`, which auto-detects based on model type
  (`TRUE` for binomial/survival, `FALSE` for linear).

- indent_groups:

  Logical. If `TRUE`, indents factor levels under their parent variable
  name, creating a hierarchical display. If `FALSE` (default), shows
  variable and level in separate columns.

- condense_table:

  Logical. If `TRUE`, condenses binary categorical variables into single
  rows by showing only the non-reference category. Automatically sets
  `indent_groups = TRUE`. Useful for tables with many binary variables.
  Default is `FALSE`.

- bold_variables:

  Logical. If `TRUE`, variable names are displayed in bold. If `FALSE`
  (default), variable names are displayed in plain text.

- center_padding:

  Numeric value specifying horizontal spacing between table and forest
  plot. Default is 4.

- zebra_stripes:

  Logical. If `TRUE`, applies alternating gray background shading to
  different variables for improved readability. Default is `TRUE`.

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

  Named character vector providing custom display labels for variables.
  Applied to predictor names in the plot. Default is `NULL` (uses
  original variable names).

- show_footer:

  Logical. If `TRUE`, displays a footer with the outcome variable name.
  Default is `TRUE`.

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

The forest plot displays univariable (unadjusted) associations between
each predictor and the outcome. This is useful for:

- Visualizing results of initial variable screening

- Identifying potential predictors for multivariable modeling

- Presenting crude associations alongside adjusted results

- Quick visual assessment of effect sizes and significance

The plot automatically handles:

- Different effect types (OR, HR, RR, coefficients) with appropriate
  axis scaling (log vs linear)

- Factor variables with multiple levels (grouped under variable name)

- Continuous variables (single row per predictor)

- Reference categories for categorical variables

## See also

[`autoforest`](https://phmcc.codeberg.page/summata/reference/autoforest.md)
for automatic model detection,
[`uniscreen`](https://phmcc.codeberg.page/summata/reference/uniscreen.md)
for generating univariable screening results,
[`multiforest`](https://phmcc.codeberg.page/summata/reference/multiforest.md)
for multi-outcome forest plots,
[`coxforest`](https://phmcc.codeberg.page/summata/reference/coxforest.md),
[`glmforest`](https://phmcc.codeberg.page/summata/reference/glmforest.md),
[`lmforest`](https://phmcc.codeberg.page/summata/reference/lmforest.md)
for single-model forest plots

Other visualization functions:
[`autoforest()`](https://phmcc.codeberg.page/summata/reference/autoforest.md),
[`coxforest()`](https://phmcc.codeberg.page/summata/reference/coxforest.md),
[`glmforest()`](https://phmcc.codeberg.page/summata/reference/glmforest.md),
[`lmforest()`](https://phmcc.codeberg.page/summata/reference/lmforest.md),
[`multiforest()`](https://phmcc.codeberg.page/summata/reference/multiforest.md)

## Examples

``` r
data(clintrial)
data(clintrial_labels)

# Create example uniscreen result
uni_results <- uniscreen(
    data = clintrial,
    outcome = "os_status",
    predictors = c("age", "sex", "smoking", "treatment", "stage"),
    labels = clintrial_labels,
    parallel = FALSE
)

# Example 1: Basic univariable forest plot
p <- uniforest(uni_results, title = "Univariable Associations with Mortality")
#> Recommended plot dimensions: width = 11.4 in, height = 6.2 in

# \donttest{

old_width <- options(width = 180)

# Example 2: Survival analysis
library(survival)
surv_results <- uniscreen(
    data = clintrial,
    outcome = "Surv(os_months, os_status)",
    predictors = c("age", "sex", "treatment", "stage"),
    model_type = "coxph",
    labels = clintrial_labels,
    parallel = FALSE
)

p2 <- uniforest(surv_results, title = "Univariable Survival Analysis")
#> Recommended plot dimensions: width = 11.2 in, height = 5.5 in

# Example 3: Linear regression
lm_results <- uniscreen(
    data = clintrial,
    outcome = "los_days",
    predictors = c("age", "sex", "surgery", "diabetes"),
    model_type = "lm",
    labels = clintrial_labels,
    parallel = FALSE
)

p3 <- uniforest(lm_results, title = "Predictors of Length of Stay")
#> Recommended plot dimensions: width = 10.0 in, height = 5.0 in

# Example 4: Customize appearance
p4 <- uniforest(
    uni_results,
    title = "Crude Associations with Mortality",
    color = "#E74C3C",
    indent_groups = TRUE,
    zebra_stripes = TRUE,
    bold_variables = TRUE
)
#> Recommended plot dimensions: width = 10.2 in, height = 7.2 in

# Example 5: Save with recommended dimensions
dims <- attr(p4, "rec_dims")
ggplot2::ggsave(file.path(tempdir(), "univariable_forest.pdf"),
                p4, width = dims$width, height = dims$height)

options(old_width)

# }
```
