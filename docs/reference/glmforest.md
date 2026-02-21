# Create Forest Plot for Generalized Linear Models

Generates a publication-ready forest plot that combines a formatted data
table with a graphical representation of effect estimates (odds ratios,
risk ratios, or coefficients) from a generalized linear model. The plot
integrates variable names, group levels, sample sizes, effect estimates
with confidence intervals, *p*-values, and model diagnostics in a single
comprehensive visualization designed for manuscripts and presentations.

## Usage

``` r
glmforest(
  x,
  data = NULL,
  title = "Generalized Linear Model",
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
  show_events = TRUE,
  indent_groups = FALSE,
  condense_table = FALSE,
  bold_variables = FALSE,
  center_padding = 4,
  zebra_stripes = TRUE,
  ref_label = "reference",
  labels = NULL,
  color = NULL,
  exponentiate = NULL,
  qc_footer = TRUE,
  units = "in",
  number_format = NULL
)
```

## Arguments

- x:

  Either a fitted GLM object (class `glm` or `glmerMod`), a `fit_result`
  object from
  [`fit()`](https://phmcc.github.io/summata/reference/fit.md), or a
  `fullfit_result` object from
  [`fullfit()`](https://phmcc.github.io/summata/reference/fullfit.md).
  When a `fit_result` or `fullfit_result` is provided, the model, data,
  and labels are automatically extracted.

- data:

  Data frame or data.table containing the original data used to fit the
  model. If `NULL` (default) and `x` is a model, the function attempts
  to extract data from the model object. If `x` is a `fit_result`, data
  is extracted automatically. Providing data explicitly is recommended
  when passing a model directly.

- title:

  Character string specifying the plot title displayed at the top.
  Default is `"Generalized Linear Model"`. Use descriptive titles like
  "Risk Factors for Disease Outcome" for publication.

- effect_label:

  Character string for the effect measure label on the forest plot axis.
  If `NULL` (default), automatically determined based on model family
  and link function: "Odds Ratio" for logistic regression
  (`family = binomial, link = logit`), "Risk Ratio" for log-link models,
  "Exp(Coefficient)" for other exponential families, or "Coefficient"
  for identity link.

- digits:

  Integer specifying the number of decimal places for effect estimates
  and confidence intervals in the data table. Default is 2.

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
  elements. Values \> 1 increase all fonts proportionally, values \< 1
  decrease them. Default is 1.0. Useful for adjusting readability across
  different output sizes.

- annot_size:

  Numeric value controlling the relative font size for data annotations
  (variable names, values in table cells). Default is 3.88. Adjust
  relative to `font_size`.

- header_size:

  Numeric value controlling the relative font size for column headers
  ("Variable", "Group", "n", *etc.*). Default is 5.82. Headers are
  typically larger than annotations for hierarchy.

- title_size:

  Numeric value controlling the relative font size for the main plot
  title. Default is 23.28. The title is typically the largest text
  element.

- plot_width:

  Numeric value specifying the intended output width in specified
  `units`. Used for optimizing layout and text sizing. Default is `NULL`
  (automatic). Recommended: 10-16 inches for standard publications.

- plot_height:

  Numeric value specifying the intended output height in specified
  `units`. Default is `NULL` (automatic based on number of rows). The
  function provides recommendations if not specified.

- table_width:

  Numeric value between 0 and 1 specifying the proportion of total plot
  width allocated to the data table (left side). The forest plot
  occupies `1 - table_width`. Default is 0.6 (60% table, 40% forest).
  Increase for longer variable names, decrease to emphasize the forest
  plot.

- show_n:

  Logical. If `TRUE`, includes a column showing group-specific sample
  sizes for categorical variables and total sample size for continuous
  variables. Default is `TRUE`.

- show_events:

  Logical. If `TRUE`, includes a column showing the number of events for
  each group. Relevant for logistic regression (number of cases) and
  other binary outcomes. Default is `TRUE`.

- indent_groups:

  Logical. If `TRUE`, indents factor levels under their parent variable
  name, creating a hierarchical visual structure. When `TRUE`, the
  "Group" column is hidden. Default is `FALSE`.

- condense_table:

  Logical. If `TRUE`, condenses binary categorical variables into single
  rows by showing only the non-reference category. Automatically sets
  `indent_groups = TRUE`. Useful for tables with many binary variables.
  Default is `FALSE`.

- bold_variables:

  Logical. If `TRUE`, variable names are displayed in bold. If `FALSE`
  (default), variable names are displayed in plain text.

- center_padding:

  Numeric value specifying the horizontal spacing (in character units)
  between the data table and forest plot. Increase for more separation,
  decrease to fit more content. Default is 4.

- zebra_stripes:

  Logical. If `TRUE`, applies alternating gray background shading to
  different variables (not rows) to improve visual grouping and
  readability. Default is `TRUE`.

- ref_label:

  Character string to display for reference categories of factor
  variables. Typically shown in place of effect estimates. Default is
  `"reference"`. Common alternatives: "ref", "1.00 (ref)".

- labels:

  Named character vector or list providing custom display labels for
  variables. Names should match variable names in the model, values are
  the labels to display. Example:
  `c(age = "Age (years)", bmi = "Body Mass Index")`. Default is `NULL`
  (use original variable names).

- color:

  Character string specifying the color for effect estimate point
  markers in the forest plot. Use hex codes or R color names. Default is
  `NULL`, which auto-selects based on effect type: `"#4BA6B6"` (teal)
  for odds ratios (binomial/quasibinomial with logit link), `"#3F87EE"`
  (blue) for rate/risk ratios (Poisson, Gamma, inverse

  Gaussian with log link), and `"#5A8F5A"` (green) for coefficients
  (Gaussian/identity link). This scheme matches
  [`uniforest()`](https://phmcc.github.io/summata/reference/uniforest.md)
  and
  [`multiforest()`](https://phmcc.github.io/summata/reference/multiforest.md).
  Choose colors that contrast well with black error bars.

- exponentiate:

  Logical. If `TRUE`, exponentiates coefficients to display odds ratios,
  risk ratios, *etc.* If `FALSE`, shows raw coefficients. Default is
  `NULL`, which automatically exponentiates for logit, log, and cloglog
  links, and shows raw coefficients for identity link.

- qc_footer:

  Logical. If `TRUE`, displays model quality control statistics in the
  footer (observations analyzed, model family, deviance, pseudo-*R*²,
  AIC). Default is `TRUE`.

- units:

  Character string specifying the units for plot dimensions. Options:
  `"in"` (inches), `"cm"` (centimeters), `"mm"` (millimeters). Default
  is `"in"`. Affects interpretation of `plot_width` and `plot_height`.

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

**Plot Components:**

The forest plot consists of several integrated components:

1.  **Title**: Centered at top, describes the analysis

2.  **Data Table** (left side): Contains columns for:

    - Variable: Predictor names (or custom labels)

    - Group: Factor levels (optional, hidden when indenting)

    - n: Sample sizes by group (optional)

    - Events: Event counts by group (optional)

    - Effect (95% CI); *p*-value: Formatted estimates with *p*-values

3.  **Forest Plot** (right side): Graphical display with:

    - Point estimates (squares sized by sample size)

    - 95% confidence intervals (error bars)

    - Reference line (at OR/RR = 1 or coefficient = 0)

    - Log scale for odds/risk ratios

    - Labeled axis

4.  **Model Statistics** (footer): Summary of:

    - Observations analyzed (with percentage of total data)

    - Model family (Binomial, Poisson, *etc.*)

    - Deviance statistics

    - Pseudo-*R*² (McFadden)

    - AIC

**Automatic Effect Measure Selection:**

When `effect_label = NULL` and `exponentiate = NULL`, the function
intelligently selects the appropriate effect measure:

- **Logistic regression** (`family = binomial(link = "logit")`): Odds
  Ratios (OR)

- **Log-link models** (`link = "log"`): Risk Ratios (RR) or Rate Ratios

- **Other exponential families**: exp(coefficient)

- **Identity link**: Raw coefficients

**Reference Categories:**

For factor variables, the first level (determined by factor ordering or
alphabetically for character variables) serves as the reference
category:

- Displayed with the `ref_label` instead of an estimate

- No confidence interval or *p*-value shown

- Visually aligned with other categories

- When `condense_table = TRUE`, reference-only variables may be omitted
  entirely

**Layout Optimization:**

The function automatically optimizes layout based on content:

- Calculates appropriate axis ranges to accommodate all confidence
  intervals

- Selects meaningful tick marks on log or linear scales

- Sizes point markers proportional to sample size (larger = more data)

- Adjusts table width based on variable name lengths when
  `table_width = NULL`

- Recommends overall dimensions based on number of rows

**Visual Grouping Options:**

Three display modes are available:

1.  **Standard** (`indent_groups = FALSE`, `condense_table = FALSE`):
    Separate "Variable" and "Group" columns, all categories shown

2.  **Indented** (`indent_groups = TRUE`, `condense_table = FALSE`):
    Hierarchical display with groups indented under variables

3.  **Condensed** (`condense_table = TRUE`): Binary variables shown in
    single rows, automatically indented

**Zebra Striping:**

When `zebra_stripes = TRUE`, alternating variables (not individual rows)
receive light gray backgrounds. This helps visually group all levels of
a factor variable together, making the plot easier to read especially
with many multi-level factors.

**Model Statistics Display:**

The footer shows key diagnostic information:

- **Observations analyzed**: Total N and percentage of original data
  (accounting for missing values)

- **Null/Residual Deviance**: Model fit improvement

- **Pseudo-*R*²**: McFadden *R*² = 1 - (log L_1 / log L_2)

- **AIC**: For model comparison (lower is better)

For logistic regression, concordance (C-statistic/AUC) may also be
displayed if available.

**Saving Plots:**

Use
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
with recommended dimensions:

      p <- glmforest(model, data)
      dims <- attr(p, "rec_dims")
      ggsave("forest.pdf", p, width = dims$width, height = dims$height)

Or specify custom dimensions:

    ggsave("forest.png", p, width = 12, height = 8, dpi = 300)

## See also

[`autoforest`](https://phmcc.github.io/summata/reference/autoforest.md)
for automatic model detection,
[`coxforest`](https://phmcc.github.io/summata/reference/coxforest.md)
for Cox proportional hazards forest plots,
[`lmforest`](https://phmcc.github.io/summata/reference/lmforest.md) for
linear model forest plots,
[`uniforest`](https://phmcc.github.io/summata/reference/uniforest.md)
for univariable screening forest plots,
[`multiforest`](https://phmcc.github.io/summata/reference/multiforest.md)
for multi-outcome forest plots,
[`glm`](https://rdrr.io/r/stats/glm.html) for fitting GLMs,
[`fit`](https://phmcc.github.io/summata/reference/fit.md) for regression
modeling

Other visualization functions:
[`autoforest()`](https://phmcc.github.io/summata/reference/autoforest.md),
[`coxforest()`](https://phmcc.github.io/summata/reference/coxforest.md),
[`lmforest()`](https://phmcc.github.io/summata/reference/lmforest.md),
[`multiforest()`](https://phmcc.github.io/summata/reference/multiforest.md),
[`uniforest()`](https://phmcc.github.io/summata/reference/uniforest.md)

## Examples

``` r
if (FALSE) {
# Load example data
data(clintrial)
data(clintrial_labels)

# Example 1: Basic logistic regression forest plot
model1 <- glm(os_status ~ age + sex + bmi + treatment,
              data = clintrial,
              family = binomial)

plot1 <- glmforest(model1, data = clintrial)
print(plot1)

  options(width = 180)
# Example 2: With custom variable labels
plot2 <- glmforest(
    x = model1,
    data = clintrial,
    title = "Risk Factors for Mortality",
    labels = clintrial_labels
)
print(plot2)

# Example 3: Customize appearance
plot3 <- glmforest(
    x = model1,
    data = clintrial,
    title = "Adjusted Odds Ratios",
    color = "#D62728",  # Red points
    font_size = 1.2,    # Larger text
    zebra_stripes = FALSE,
    labels = clintrial_labels
)
print(plot3)

# Example 4: Indented layout for hierarchical view
plot4 <- glmforest(
    x = model1,
    data = clintrial,
    indent_groups = TRUE,
    labels = clintrial_labels
)
print(plot4)
# Group column hidden, levels indented under variables

# Example 5: Condensed layout for many binary variables
model5 <- glm(os_status ~ age + sex + smoking + hypertension + 
                  diabetes + surgery,
              data = clintrial,
              family = binomial)

plot5 <- glmforest(
    x = model5,
    data = clintrial,
    condense_table = TRUE,
    labels = clintrial_labels
)
print(plot5)
# Binary variables shown in single rows

# Example 6: Hide sample size and events columns
plot6 <- glmforest(
    x = model1,
    data = clintrial,
    show_n = FALSE,
    show_events = FALSE,
    labels = clintrial_labels
)
print(plot6)

# Example 7: Adjust table/forest proportions
plot7 <- glmforest(
    x = model1,
    data = clintrial,
    table_width = 0.7,  # More space for table
    labels = clintrial_labels
)
print(plot7)

# Example 8: Get and use recommended dimensions
plot8 <- glmforest(model1, data = clintrial)

dims <- attr(plot8, "rec_dims")
cat("Recommended: ", dims$width, "x", dims$height, "inches\n")

# Save with recommended dimensions
# ggsave("forest.pdf", plot8, width = dims$width, height = dims$height)

# Example 9: Specify exact output dimensions
plot9 <- glmforest(
    x = model1,
    data = clintrial,
    plot_width = 14,
    plot_height = 10,
    labels = clintrial_labels
)
# No dimension recommendations printed

# Example 10: Use different units (centimeters)
plot10 <- glmforest(
    x = model1,
    data = clintrial,
    plot_width = 30,  # 30 cm
    plot_height = 20,  # 20 cm
    units = "cm",
    labels = clintrial_labels
)

# Example 11: Poisson regression for count data
model11 <- glm(ae_count ~ age + treatment + diabetes + surgery,
               data = clintrial,
               family = poisson)

plot11 <- glmforest(
    x = model11,
    data = clintrial,
    title = "Rate Ratios for Adverse Events",
    labels = clintrial_labels
)
print(plot11)
# Shows rate ratios with blue color for count models

# Example 12: Gamma regression for positive continuous outcomes
model12 <- glm(los_days ~ age + treatment + surgery + stage,
               data = clintrial,
               family = Gamma(link = "log"))

plot12 <- glmforest(
    x = model12,
    data = clintrial,
    title = "Multiplicative Effects on Length of Stay",
    labels = clintrial_labels
)
print(plot12)
# Shows ratios with blue color (same as other ratio models)

# Example 13: Force display of raw coefficients
plot13 <- glmforest(
    x = model1,
    data = clintrial,
    exponentiate = FALSE,  # Show log odds
    effect_label = "Log Odds",
    labels = clintrial_labels
)
print(plot13)
# Reference line at 0 instead of 1

# Example 14: Custom reference label
plot14 <- glmforest(
    x = model1,
    data = clintrial,
    ref_label = "1.00 (ref)",
    labels = clintrial_labels
)
print(plot14)

# Example 15: Adjust font sizes for presentations
plot15 <- glmforest(
    x = model1,
    data = clintrial,
    font_size = 1.5,      # 50% larger
    title_size = 30,      # Larger title
    labels = clintrial_labels
)
print(plot15)

# Example 16: Complete publication-ready plot
final_model <- glm(
    os_status ~ age + sex + bmi + smoking + hypertension + 
        diabetes + treatment + stage,
    data = clintrial,
    family = binomial
)

final_plot <- glmforest(
    x = final_model,
    data = clintrial,
    title = "Multivariable Logistic Regression: Risk Factors for Mortality",
    labels = clintrial_labels,
    indent_groups = TRUE,
    zebra_stripes = TRUE,
    color = "#4BA6B6",
    font_size = 1.0,
    digits = 2
)

# Save for publication
dims <- attr(final_plot, "rec_dims")
# ggsave("figure1.pdf", final_plot, 
#        width = dims$width, height = dims$height)
# ggsave("figure1.png", final_plot,
#        width = dims$width, height = dims$height, dpi = 300)
}
```
