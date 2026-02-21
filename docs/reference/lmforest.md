# Create Forest Plot for Linear Models

Generates a publication-ready forest plot that combines a formatted data
table with a graphical representation of regression coefficients from a
linear model. The plot integrates variable names, group levels, sample
sizes, coefficients with confidence intervals, *p*-values, and model
diagnostics (*R*², *F*-statistic, AIC) in a single comprehensive
visualization designed for manuscripts and presentations.

## Usage

``` r
lmforest(
  x,
  data = NULL,
  title = "Linear Model",
  effect_label = "Coefficient",
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
  indent_groups = FALSE,
  condense_table = FALSE,
  bold_variables = FALSE,
  center_padding = 4,
  zebra_stripes = TRUE,
  ref_label = "reference",
  labels = NULL,
  units = "in",
  color = "#5A8F5A",
  qc_footer = TRUE,
  number_format = NULL
)
```

## Arguments

- x:

  Either a fitted linear model object (class `lm` or `lmerMod`), a
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
  Default is `"Linear Model"`.

- effect_label:

  Character string for the effect measure label on the forest plot axis.
  Default is `"Coefficient"`.

- digits:

  Integer specifying the number of decimal places for coefficients and
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
  `units`. Default is `NULL` (automatic).

- plot_height:

  Numeric value specifying the intended output height in specified
  `units`. Default is `NULL` (automatic).

- table_width:

  Numeric value between 0 and 1 specifying the proportion of total plot
  width allocated to the data table. Default is 0.6.

- show_n:

  Logical. If `TRUE`, includes a column showing group-specific sample
  sizes. Default is `TRUE`.

- indent_groups:

  Logical. If `TRUE`, indents factor levels under their parent variable
  name, creating hierarchical structure. The "Group" column is hidden
  when `TRUE`. Default is `FALSE`.

- condense_table:

  Logical. If `TRUE`, condenses binary categorical variables into single
  rows. Automatically sets `indent_groups = TRUE`. Default is `FALSE`.

- bold_variables:

  Logical. If `TRUE`, variable names are displayed in bold. If `FALSE`
  (default), variable names are displayed in plain text.

- center_padding:

  Numeric value specifying horizontal spacing between table and forest
  plot. Default is 4.

- zebra_stripes:

  Logical. If `TRUE`, applies alternating gray background shading to
  different variables. Default is `TRUE`.

- ref_label:

  Character string to display for reference categories of factor
  variables. Default is `"reference"`.

- labels:

  Named character vector providing custom display labels for variables.
  Example: `c(age = "Age (years)", height = "Height (cm)")`. Default is
  `NULL`.

- units:

  Character string specifying units for plot dimensions: `"in"`
  (inches), `"cm"`, or `"mm"`. Default is `"in"`.

- color:

  Character string specifying the color for coefficient point estimates
  in the forest plot. Default is `"#5A8F5A"` (green). Use hex codes or R
  color names.

- qc_footer:

  Logical. If `TRUE`, displays model quality control statistics in the
  footer (observations analyzed, *R*², adjusted *R*², *F*-statistic,
  AIC). Default is `TRUE`.

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

**Linear Model-Specific Features:**

The linear model forest plot differs from logistic and Cox plots in
several ways:

- **Coefficients**: Raw regression coefficients shown (not
  exponentiated)

- **Reference line**: At coefficient = 0 (not at 1)

- **Linear scale**: Forest plot uses linear scale (not log scale)

- **No events column**: Only sample sizes shown (no event counts)

- ***R*² statistics**: Model fit assessed by *R*² and adjusted *R*²

- ***F*-test**: Overall model significance from *F*-statistic

**Plot Components:**

1.  **Title**: Centered at top

2.  **Data Table** (left): Contains:

    - Variable: Predictor names

    - Group: Factor levels (if applicable)

    - *n*: Sample sizes by group

    - Coefficient (95% CI); *p*-value: Raw coefficients with CIs and
      *p*-values

3.  **Forest Plot** (right):

    - Point estimates (squares sized by sample size)

    - 95% confidence intervals (error bars)

    - Reference line at coefficient = 0

    - Linear scale

4.  **Model Statistics** (footer):

    - Observations analyzed (with percentage of total data)

    - *R*² and adjusted *R*²

    - *F*-statistic with degrees of freedom and *p*-value

    - AIC

**Interpreting Coefficients:**

Linear regression coefficients represent the change in the outcome
variable for a one-unit change in the predictor:

- **Continuous predictors**: Coefficient = change in Y per unit of X

- **Binary predictors**: Coefficient = difference in Y between groups

- **Factor predictors**: Coefficients = differences from reference
  category

- **Sign matters**: Positive = increase in Y, Negative = decrease in Y

- **Zero crossing**: CI crossing zero suggests no significant effect

Example: If the coefficient for "age" is 0.50 when predicting BMI, BMI
increases by 0.50 kg/m² for each additional year of age.

**Model Fit Statistics:**

The footer displays key diagnostics:

- ***R*²**: Proportion of variance explained (0 to 1)

  - 0.0-0.3: Weak explanatory power

  - 0.3-0.5: Moderate

  - 0.5-0.7: Good

  - \> 0.7: Strong (rare in social/biological sciences)

- **Adjusted *R*²**: *R*² penalized for number of predictors

  - Always ≤ *R*²

  - Preferred for model comparison

  - Accounts for model complexity

- ***F*-statistic**: Tests null hypothesis that all coefficients = 0

  - Degrees of freedom: df1 = \# predictors, df2 = \# observations - \#
    predictors - 1

  - Significant *p*-value indicates model explains variance better than
    intercept-only

- **AIC**: For model comparison (lower is better)

**Assumptions:**

Linear regression assumes:

1.  Linearity of relationships

2.  Independence of observations

3.  Homoscedasticity (constant variance)

4.  Normality of residuals

5.  No multicollinearity

Check assumptions using:

- `plot(model)` for diagnostic plots

- `car::vif(model)` for multicollinearity

- `lmtest::bptest(model)` for heteroscedasticity

- `shapiro.test(residuals(model))` for normality

**Reference Categories:**

For factor variables:

- First level is the reference (coefficient = 0)

- Other levels show difference from reference

- Reference displayed with `ref_label`

- Relevel factors before modeling if needed:
  `factor(x, levels = c("desired_ref", ...))`

**Sample Size Reporting:**

The "*n*" column shows:

- For continuous variables: Total observations with non-missing data

- For factor variables: Number of observations in each category

- Footer shows total observations analyzed and percentage of original
  data (accounting for missing values)

## See also

[`autoforest`](https://phmcc.github.io/summata/reference/autoforest.md)
for automatic model detection,
[`glmforest`](https://phmcc.github.io/summata/reference/glmforest.md)
for logistic/GLM forest plots,
[`coxforest`](https://phmcc.github.io/summata/reference/coxforest.md)
for Cox model forest plots,
[`uniforest`](https://phmcc.github.io/summata/reference/uniforest.md)
for univariable screening forest plots,
[`multiforest`](https://phmcc.github.io/summata/reference/multiforest.md)
for multi-outcome forest plots, [`lm`](https://rdrr.io/r/stats/lm.html)
for fitting linear models,
[`fit`](https://phmcc.github.io/summata/reference/fit.md) for regression
modeling

Other visualization functions:
[`autoforest()`](https://phmcc.github.io/summata/reference/autoforest.md),
[`coxforest()`](https://phmcc.github.io/summata/reference/coxforest.md),
[`glmforest()`](https://phmcc.github.io/summata/reference/glmforest.md),
[`multiforest()`](https://phmcc.github.io/summata/reference/multiforest.md),
[`uniforest()`](https://phmcc.github.io/summata/reference/uniforest.md)

## Examples

``` r
if (FALSE) {
# Load example data
data(clintrial)
data(clintrial_labels)

# Example 1: Basic linear model forest plot
model1 <- lm(bmi ~ age + sex + smoking,
             data = clintrial)

plot1 <- lmforest(model1, data = clintrial)
print(plot1)

  options(width = 180)
# Example 2: With custom labels and title
plot2 <- lmforest(
    x = model1,
    data = clintrial,
    title = "Predictors of Body Mass Index",
    effect_label = "Change in BMI (kg/m²)",
    labels = clintrial_labels
)
print(plot2)

# Example 3: More comprehensive model
model3 <- lm(
    bmi ~ age + sex + smoking + hypertension + diabetes + creatinine,
    data = clintrial
)

plot3 <- lmforest(
    x = model3,
    data = clintrial,
    labels = clintrial_labels,
    indent_groups = TRUE
)
print(plot3)

# Example 4: Check model diagnostics
summary(model3)
plot(model3)  # Diagnostic plots

# Example 5: Condensed layout
model5 <- lm(
    bmi ~ age + sex + smoking + hypertension + diabetes,
    data = clintrial
)

plot5 <- lmforest(
    x = model5,
    data = clintrial,
    condense_table = TRUE,
    labels = clintrial_labels
)
print(plot5)

# Example 6: Custom color
plot6 <- lmforest(
    x = model1,
    data = clintrial,
    color = "#2ECC71",  # Bright green
    labels = clintrial_labels
)
print(plot6)

# Example 7: Hide sample sizes
plot7 <- lmforest(
    x = model1,
    data = clintrial,
    show_n = FALSE,
    labels = clintrial_labels
)
print(plot7)

# Example 8: Adjust table width
plot8 <- lmforest(
    x = model3,
    data = clintrial,
    table_width = 0.55,
    labels = clintrial_labels
)
print(plot8)

# Example 9: Get recommended dimensions
plot9 <- lmforest(model1, data = clintrial)
dims <- attr(plot9, "rec_dims")
cat("Recommended:", dims$width, "x", dims$height, dims$units, "\n")

# Example 10: Specify exact dimensions
plot10 <- lmforest(
    x = model1,
    data = clintrial,
    plot_width = 12,
    plot_height = 7,
    labels = clintrial_labels
)

# Example 11: Different units
plot11 <- lmforest(
    x = model1,
    data = clintrial,
    plot_width = 30,
    plot_height = 20,
    units = "cm",
    labels = clintrial_labels
)

# Example 12: No zebra stripes
plot12 <- lmforest(
    x = model1,
    data = clintrial,
    zebra_stripes = FALSE,
    labels = clintrial_labels
)
print(plot12)

# Example 13: Custom reference label
plot13 <- lmforest(
    x = model1,
    data = clintrial,
    ref_label = "0.00 (ref)",
    labels = clintrial_labels
)
print(plot13)

# Example 14: Larger fonts for presentation
plot14 <- lmforest(
    x = model1,
    data = clintrial,
    font_size = 1.3,
    title_size = 26,
    labels = clintrial_labels
)
print(plot14)

# Example 15: More decimal places
plot15 <- lmforest(
    x = model1,
    data = clintrial,
    digits = 3,
    labels = clintrial_labels
)
print(plot15)

# Example 16: Hemoglobin as outcome
model16 <- lm(
    hemoglobin ~ age + sex + bmi + smoking + creatinine,
    data = clintrial
)

plot16 <- lmforest(
    x = model16,
    data = clintrial,
    title = "Predictors of Baseline Hemoglobin",
    effect_label = "Change in Hemoglobin (g/dL)",
    labels = clintrial_labels
)
print(plot16)

# Check assumptions
par(mfrow = c(2, 2))
plot(model16)
par(mfrow = c(1, 1))

# Example 17: Publication-ready final plot
final_model <- lm(
    bmi ~ age + sex + race + smoking + hypertension + 
        diabetes + creatinine + hemoglobin,
    data = clintrial
)

# Check model fit
summary(final_model)

final_plot <- lmforest(
    x = final_model,
    data = clintrial,
    title = "Multivariable Linear Regression: Predictors of Body Mass Index",
    effect_label = "Change in BMI (kg/m²)",
    labels = clintrial_labels,
    indent_groups = TRUE,
    zebra_stripes = TRUE,
    show_n = TRUE,
    color = "#5A8F5A",
    digits = 2
)

# Save for publication
dims <- attr(final_plot, "rec_dims")
# ggsave("figure3_linear.pdf", final_plot,
#        width = dims$width, height = dims$height)
# ggsave("figure3_linear.png", final_plot,
#        width = dims$width, height = dims$height, dpi = 300)
}
```
