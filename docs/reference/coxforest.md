# Create Forest Plot for Cox Proportional Hazards Models

Generates a publication-ready forest plot that combines a formatted data
table with a graphical representation of hazard ratios from a Cox
proportional hazards survival model. The plot integrates variable names,
group levels, sample sizes, event counts, hazard ratios with confidence
intervals, *p*-values, and model diagnostics in a single comprehensive
visualization designed for manuscripts and presentations.

## Usage

``` r
coxforest(
  x,
  data = NULL,
  title = "Cox Proportional Hazards Model",
  effect_label = "Hazard Ratio",
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
  color = "#8A61D8",
  qc_footer = TRUE,
  units = "in",
  number_format = NULL
)
```

## Arguments

- x:

  Either a fitted Cox model object (class `coxph` or `coxme`), a
  `fit_result` object from
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
  Default is `"Cox Proportional Hazards Model"`. Use descriptive titles
  for publication.

- effect_label:

  Character string for the effect measure label on the forest plot axis.
  Default is `"Hazard Ratio"`. Alternatives might include "HR" for
  space-constrained plots.

- digits:

  Integer specifying the number of decimal places for hazard ratios and
  confidence intervals. Default is 2.

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

  Logical. If `TRUE`, includes a column showing group-specific sample
  sizes. Default is `TRUE`.

- show_events:

  Logical. If `TRUE`, includes a column showing the number of events
  (deaths, failures) for each group. Critical for survival analysis
  interpretation. Default is `TRUE`.

- indent_groups:

  Logical. If `TRUE`, indents factor levels under their parent variable,
  creating hierarchical structure. The "Group" column is hidden when
  `TRUE`. Default is `FALSE`.

- condense_table:

  Logical. If `TRUE`, condenses binary variables to show only the
  non-reference level (*e.g.,* "Yes" for Yes/No variables), with the
  variable name indicating the displayed level. This creates a more
  compact table for models with many binary predictors. Default is
  `FALSE`.

- bold_variables:

  Logical. If `TRUE`, variable names are displayed in bold. If `FALSE`
  (default), variable names are displayed in plain text.

- center_padding:

  Numeric value specifying horizontal spacing between table and forest
  plot. Default is 4.

- zebra_stripes:

  Logical. If `TRUE`, applies alternating gray background shading to
  different variables for improved readability. Default is `TRUE`.

- ref_label:

  Character string to display for reference categories. Default is
  `"reference"`.

- labels:

  Named character vector providing custom display labels for variables.
  Example: `c(age = "Age (years)", stage = "Disease Stage")`. Default is
  `NULL`.

- color:

  Character string specifying the color for hazard ratio point estimates
  in the forest plot. Default is `"#8A61D8"` (purple). Use hex codes or
  R color names.

- qc_footer:

  Logical. If `TRUE`, displays model quality control statistics in the
  footer (events analyzed, global log-rank *p*-value, concordance index,
  AIC). Default is `TRUE`.

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

**Survival-Specific Features:**

The Cox forest plot includes several survival analysis-specific
components:

- **Event counts**: Number of events (deaths, failures) shown for each
  predictor category, critical for assessing statistical power

- **Hazard ratios**: Always exponentiated coefficients (never raw),
  interpreted as the multiplicative change in hazard

- **Log scale**: Forest plot uses log scale for HR (reference line at 1)

- **Model diagnostics**: Includes concordance (C-index), global log-rank
  test *p*-value, and AIC

**Plot Components:**

1.  **Title**: Centered at top

2.  **Data Table** (left): Contains:

    - Variable and Group columns

    - n: Sample sizes by group

    - Events: Event counts by group (critical for survival)

    - aHR (95% CI); *p*-value: Adjusted hazard ratios with CIs and
      *p*-values

3.  **Forest Plot** (right):

    - Point estimates (squares sized by sample size)

    - 95% confidence intervals

    - Reference line at HR = 1

    - Log scale for hazard ratios

4.  **Model Statistics** (footer):

    - Events analyzed (with percentage of total)

    - Global log-rank test *p*-value

    - Concordance (C-index) with standard error

    - AIC

**Interpreting Hazard Ratios:**

- **HR = 1**: No effect on hazard (reference)

- **HR \> 1**: Increased hazard (worse survival)

- **HR \< 1**: Decreased hazard (better survival)

- Example: HR = 2.0 means twice the hazard of the event at any time

**Event Counts:**

The "Events" column is particularly important in survival analysis:

- Indicates the number of actual events (not censored observations) in
  each group

- Essential for assessing statistical power

- Categories with very few events may have unreliable HR estimates

- The footer shows total events analyzed and percentage of all events in
  the original data

**Concordance (C-index):**

The concordance statistic displayed in the footer indicates
discrimination:

- Range: 0.5 to 1.0

- 0.5 = random prediction (coin flip)

- 0.7-0.8 = acceptable discrimination

- \> 0.8 = excellent discrimination

- Standard error provided for confidence interval calculation

**Global Log-Rank Test:**

The global *p*-value tests the null hypothesis that all coefficients are
zero:

- Significant *p*-value (\< 0.05) indicates the model as a whole
  predicts survival

- Non-significant global test doesn't preclude significant individual
  predictors

- Based on the score (log-rank) test

**Stratification and Clustering:**

If the model includes stratification
([`strata()`](https://rdrr.io/pkg/survival/man/strata.html)) or
clustering
([`cluster()`](https://rdrr.io/pkg/survival/man/cluster.html)):

- Stratified variables are not shown in the forest plot (they don't have
  HRs)

- Clustering affects standard errors but not point estimates

- Both are handled automatically by the function

**Proportional Hazards Assumption:**

The forest plot assumes proportional hazards (constant HR over time).
Users should verify this assumption using:

- `cox.zph(model)` for testing

- Stratification for variables violating the assumption

- Time-dependent coefficients if needed

## See also

[`autoforest`](https://phmcc.github.io/summata/reference/autoforest.md)
for automatic model detection,
[`glmforest`](https://phmcc.github.io/summata/reference/glmforest.md)
for logistic/GLM forest plots,
[`lmforest`](https://phmcc.github.io/summata/reference/lmforest.md) for
linear model forest plots,
[`uniforest`](https://phmcc.github.io/summata/reference/uniforest.md)
for univariable screening forest plots,
[`multiforest`](https://phmcc.github.io/summata/reference/multiforest.md)
for multi-outcome forest plots,
[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) for fitting Cox
models, [`fit`](https://phmcc.github.io/summata/reference/fit.md) for
regression modeling

Other visualization functions:
[`autoforest()`](https://phmcc.github.io/summata/reference/autoforest.md),
[`glmforest()`](https://phmcc.github.io/summata/reference/glmforest.md),
[`lmforest()`](https://phmcc.github.io/summata/reference/lmforest.md),
[`multiforest()`](https://phmcc.github.io/summata/reference/multiforest.md),
[`uniforest()`](https://phmcc.github.io/summata/reference/uniforest.md)

## Examples

``` r
if (FALSE) {
# Load example data
data(clintrial)
data(clintrial_labels)
library(survival)

# Example 1: Basic Cox model forest plot
model1 <- coxph(Surv(os_months, os_status) ~ age + sex + treatment,
                data = clintrial)

plot1 <- coxforest(model1, data = clintrial)
print(plot1)

  options(width = 180)
# Example 2: With custom labels and title
plot2 <- coxforest(
    x = model1,
    data = clintrial,
    title = "Prognostic Factors for Overall Survival",
    labels = clintrial_labels
)
print(plot2)

# Example 3: Comprehensive multivariable model
model3 <- coxph(
    Surv(os_months, os_status) ~ age + sex + bmi + smoking + 
        treatment + stage + grade,
    data = clintrial
)

plot3 <- coxforest(
    x = model3,
    data = clintrial,
    labels = clintrial_labels,
    indent_groups = TRUE
)
print(plot3)

# Example 4: Condensed layout for many binary predictors
model4 <- coxph(
    Surv(os_months, os_status) ~ age + sex + smoking + 
        hypertension + diabetes + surgery,
    data = clintrial
)

plot4 <- coxforest(
    x = model4,
    data = clintrial,
    condense_table = TRUE,
    labels = clintrial_labels
)
print(plot4)

# Example 5: Custom color scheme
plot5 <- coxforest(
    x = model1,
    data = clintrial,
    color = "#E74C3C",  # Red
    zebra_stripes = FALSE,
    labels = clintrial_labels
)
print(plot5)

# Example 6: Hide sample sizes, show only events
plot6 <- coxforest(
    x = model1,
    data = clintrial,
    show_n = FALSE,
    show_events = TRUE,
    labels = clintrial_labels
)
print(plot6)

# Example 7: Adjust table width
plot7 <- coxforest(
    x = model3,
    data = clintrial,
    table_width = 0.65,  # More space for long variable names
    labels = clintrial_labels
)
print(plot7)

# Example 8: Specify exact dimensions
plot8 <- coxforest(
    x = model1,
    data = clintrial,
    plot_width = 14,
    plot_height = 8,
    labels = clintrial_labels
)

# Example 9: Use recommended dimensions for saving
plot9 <- coxforest(model1, data = clintrial)
dims <- attr(plot9, "rec_dims")

# Save to PDF
# ggsave("survival_forest.pdf", plot9,
#        width = dims$width, height = dims$height)

# Example 10: Different units (centimeters)
plot10 <- coxforest(
    x = model1,
    data = clintrial,
    plot_width = 35,
    plot_height = 25,
    units = "cm",
    labels = clintrial_labels
)

# Example 11: Stratified Cox model
# Stratification variable (site) won't appear in plot
model11 <- coxph(
    Surv(os_months, os_status) ~ age + sex + treatment + strata(site),
    data = clintrial
)

plot11 <- coxforest(
    x = model11,
    data = clintrial,
    title = "Stratified by Study Site",
    labels = clintrial_labels
)
print(plot11)

# Example 12: With clustering for robust SE
model12 <- coxph(
    Surv(os_months, os_status) ~ age + sex + treatment + cluster(site),
    data = clintrial
)

plot12 <- coxforest(
    x = model12,
    data = clintrial,
    title = "Clustered by Study Site",
    labels = clintrial_labels
)
print(plot12)
# Standard errors account for site clustering

# Example 13: Custom reference label
plot13 <- coxforest(
    x = model1,
    data = clintrial,
    ref_label = "1.00 (ref)",
    labels = clintrial_labels
)
print(plot13)

# Example 14: Presentation-sized fonts
plot14 <- coxforest(
    x = model1,
    data = clintrial,
    font_size = 1.4,
    title_size = 28,
    labels = clintrial_labels
)
print(plot14)

# Example 15: Publication-ready final plot
final_model <- coxph(
    Surv(os_months, os_status) ~ age + sex + bmi + smoking + 
        hypertension + diabetes + ecog + treatment + stage + grade,
    data = clintrial
)

final_plot <- coxforest(
    x = final_model,
    data = clintrial,
    title = "Multivariable Cox Regression: Prognostic Factors for Overall Survival",
    labels = clintrial_labels,
    indent_groups = TRUE,
    zebra_stripes = TRUE,
    show_n = TRUE,
    show_events = TRUE,
    color = "#8A61D8"
)

# Check model fit
summary(final_model)

# Verify proportional hazards assumption
# cox.zph(final_model)

# Save for publication
dims <- attr(final_plot, "rec_dims")
# ggsave("figure2_survival.pdf", final_plot,
#        width = dims$width, height = dims$height)
# ggsave("figure2_survival.png", final_plot,
#        width = dims$width, height = dims$height, dpi = 300)
}
```
