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
  [`fit()`](https://phmcc.codeberg.page/summata/reference/fit.md), or a
  `fullfit_result` object from
  [`fullfit()`](https://phmcc.codeberg.page/summata/reference/fullfit.md).
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

[`autoforest`](https://phmcc.codeberg.page/summata/reference/autoforest.md)
for automatic model detection,
[`glmforest`](https://phmcc.codeberg.page/summata/reference/glmforest.md)
for logistic/GLM forest plots,
[`lmforest`](https://phmcc.codeberg.page/summata/reference/lmforest.md)
for linear model forest plots,
[`uniforest`](https://phmcc.codeberg.page/summata/reference/uniforest.md)
for univariable screening forest plots,
[`multiforest`](https://phmcc.codeberg.page/summata/reference/multiforest.md)
for multi-outcome forest plots,
[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) for fitting Cox
models, [`fit`](https://phmcc.codeberg.page/summata/reference/fit.md)
for regression modeling

Other visualization functions:
[`autoforest()`](https://phmcc.codeberg.page/summata/reference/autoforest.md),
[`glmforest()`](https://phmcc.codeberg.page/summata/reference/glmforest.md),
[`lmforest()`](https://phmcc.codeberg.page/summata/reference/lmforest.md),
[`multiforest()`](https://phmcc.codeberg.page/summata/reference/multiforest.md),
[`uniforest()`](https://phmcc.codeberg.page/summata/reference/uniforest.md)

## Examples

``` r
data(clintrial)
data(clintrial_labels)
library(survival)

# Create example model
model1 <- coxph(
    survival::Surv(os_months, os_status) ~ age + sex + treatment,
    data = clintrial)

# Example 1: Basic Cox model forest plot
p <- coxforest(model1, data = clintrial)
#> Recommended plot dimensions: width = 12.9 in, height = 5.0 in

# \donttest{

old_width <- options(width = 180)

# Example 2: With custom labels and title
plot2 <- coxforest(
    x = model1,
    data = clintrial,
    title = "Prognostic Factors for Overall Survival",
    labels = clintrial_labels
)
#> Recommended plot dimensions: width = 13.7 in, height = 5.0 in

# Example 3: Comprehensive model with indented layout
model3 <- coxph(
    Surv(os_months, os_status) ~ age + sex + bmi + smoking + 
        treatment + stage + grade,
    data = clintrial
)

plot3 <- coxforest(
    x = model3,
    data = clintrial,
    labels = clintrial_labels,
    indent_groups = TRUE,
    zebra_stripes = TRUE
)
#> Recommended plot dimensions: width = 13.7 in, height = 8.5 in

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
#> Recommended plot dimensions: width = 12.2 in, height = 5.2 in

# Example 5: Stratified Cox model
model5 <- coxph(
    Surv(os_months, os_status) ~ age + sex + treatment + strata(site),
    data = clintrial
)

plot5 <- coxforest(
    x = model5,
    data = clintrial,
    title = "Stratified by Study Site",
    labels = clintrial_labels
)
#> Recommended plot dimensions: width = 13.7 in, height = 5.0 in

# Example 6: Save with recommended dimensions
dims <- attr(plot5, "rec_dims")
ggplot2::ggsave(file.path(tempdir(), "survival_forest.pdf"),
                plot5, width = dims$width, height = dims$height)

options(old_width)

# }
```
