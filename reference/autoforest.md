# Create Forest Plot with Automatic Model Detection

A convenience wrapper function that automatically detects the input type
and routes to the appropriate specialized forest plot function. This
eliminates the need to remember which forest function to call for
different model types or analysis objects, making it ideal for
exploratory analysis and rapid prototyping.

## Usage

``` r
autoforest(x, data = NULL, title = NULL, ...)
```

## Arguments

- x:

  One of the following:

  - A fitted model object: `glm`, `lm`, `coxph`, *etc.*

  - A `fit_result` object from
    [`fit()`](https://phmcc.codeberg.page/summata/reference/fit.md)

  - A `fullfit_result` object from
    [`fullfit()`](https://phmcc.codeberg.page/summata/reference/fullfit.md)

  - A `uniscreen_result` object from
    [`uniscreen()`](https://phmcc.codeberg.page/summata/reference/uniscreen.md)

  - A `multifit_result` object from
    [`multifit()`](https://phmcc.codeberg.page/summata/reference/multifit.md)

- data:

  Data frame or data.table containing the original data. Required when
  `x` is a raw model object. Ignored when `x` is a result object from
  [`fit()`](https://phmcc.codeberg.page/summata/reference/fit.md),
  [`fullfit()`](https://phmcc.codeberg.page/summata/reference/fullfit.md),
  [`uniscreen()`](https://phmcc.codeberg.page/summata/reference/uniscreen.md),
  or
  [`multifit()`](https://phmcc.codeberg.page/summata/reference/multifit.md)
  since these contain embedded data.

- title:

  Character string for plot title. If `NULL` (default), an appropriate
  title is generated based on the detected model type:

  - Cox models: "Cox Proportional Hazards Model"

  - Logistic regression: "Logistic Regression Model"

  - Poisson regression: "Poisson Regression Model"

  - Linear regression: "Linear Regression Model"

  - Uniscreen results: "Univariable \[Type\] Screening"

  - Multifit results: "Multivariate \[Type\] Analysis"

- ...:

  Additional arguments passed to the specific forest plot function.
  Common arguments include:

  labels

  :   Named character vector for variable labels

  digits

  :   Number of decimal places for estimates (default 2)

  p_digits

  :   Number of decimal places for *p*-values (default 3)

  conf_level

  :   Confidence level for intervals (default 0.95)

  show_n

  :   Logical, show sample sizes (default `TRUE`)

  show_events

  :   Logical, show event counts (default `TRUE` for survival/binomial)

  qc_footnote

  :   Logical, show model QC stats in footer (default `TRUE`)

  zebra_stripes

  :   Logical, alternating row shading (default `TRUE`)

  indent_groups

  :   Logical, indent factor levels (default `FALSE`)

  color

  :   Color for point estimates

  table_width

  :   Proportion of width for table (default 0.6)

  plot_width, plot_height

  :   Explicit dimensions

  units

  :   Dimension units: `"in"`, `"cm"`, or `"mm"`

  See the documentation for the specific forest function for all
  available options.

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

This function provides a convenient wrapper around the specialized
forest plot functions, automatically routing to the appropriate function
based on the model class or result type. All parameters are passed
through to the underlying function, so the full range of options remains
available.

For model-specific advanced features, individual forest functions may be
called directly.

**Automatic Detection Logic:**

The function uses the following priority order for detection:

1.  **uniscreen results**: Detected by class `"uniscreen_result"` or
    presence of attributes `outcome`, `predictors`, `model_type`, and
    `model_scope = "Univariable"`. Routes to
    [`uniforest()`](https://phmcc.codeberg.page/summata/reference/uniforest.md).

2.  **multifit results**: Detected by presence of attributes
    `predictor`, `outcomes`, `model_type`, and `raw_data`. Routes to
    [`multiforest()`](https://phmcc.codeberg.page/summata/reference/multiforest.md).

3.  **Cox models**: Classes `coxph` or `clogit`. Routes to
    [`coxforest()`](https://phmcc.codeberg.page/summata/reference/coxforest.md).

4.  **GLM models**: Class `glm`. Routes to
    [`glmforest()`](https://phmcc.codeberg.page/summata/reference/glmforest.md).

5.  **Linear models**: Class `lm` (but not `glm`). Routes to
    [`lmforest()`](https://phmcc.codeberg.page/summata/reference/lmforest.md).

## See also

[`glmforest`](https://phmcc.codeberg.page/summata/reference/glmforest.md)
for GLM forest plots,
[`coxforest`](https://phmcc.codeberg.page/summata/reference/coxforest.md)
for Cox model forest plots,
[`lmforest`](https://phmcc.codeberg.page/summata/reference/lmforest.md)
for linear model forest plots,
[`uniforest`](https://phmcc.codeberg.page/summata/reference/uniforest.md)
for univariable screening forest plots,
[`multiforest`](https://phmcc.codeberg.page/summata/reference/multiforest.md)
for multi-outcome forest plots,
[`fit`](https://phmcc.codeberg.page/summata/reference/fit.md) for
single-model regression,
[`fullfit`](https://phmcc.codeberg.page/summata/reference/fullfit.md)
for combined univariable/multivariable regression,
[`uniscreen`](https://phmcc.codeberg.page/summata/reference/uniscreen.md)
for univariable screening,
[`multifit`](https://phmcc.codeberg.page/summata/reference/multifit.md)
for multi-outcome analysis

Other visualization functions:
[`coxforest()`](https://phmcc.codeberg.page/summata/reference/coxforest.md),
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
glm_model <- glm(surgery ~ age + sex + bmi + smoking,
                 family = binomial, data = clintrial)

# Example 1: Logistic regression model
p <- autoforest(glm_model, data = clintrial)
#> Recommended plot dimensions: width = 13.3 in, height = 5.0 in
# Automatically detects GLM and routes to glmforest()

# \donttest{

# Example 2: Cox proportional hazards model
cox_model <- coxph(Surv(os_months, os_status) ~ age + sex + treatment + stage,
                   data = clintrial)

plot2 <- autoforest(cox_model, data = clintrial)
#> Recommended plot dimensions: width = 12.9 in, height = 5.5 in
# Automatically detects coxph and routes to coxforest()

# Example 3: Linear regression model
lm_model <- lm(biomarker_x ~ age + sex + bmi + treatment, data = clintrial)

plot3 <- autoforest(lm_model, data = clintrial)
#> Recommended plot dimensions: width = 13.0 in, height = 5.0 in
# Automatically detects lm and routes to lmforest()

# Example 4: With custom labels and formatting options
plot4 <- autoforest(
    cox_model,
    data = clintrial,
    labels = clintrial_labels,
    title = "Prognostic Factors for Overall Survival",
    zebra_stripes = TRUE,
    indent_groups = TRUE
)
#> Recommended plot dimensions: width = 11.8 in, height = 6.2 in

# Example 5: From fit() result - data and labels extracted automatically
fit_result <- fit(
    data = clintrial,
    outcome = "surgery",
    predictors = c("age", "sex", "bmi", "treatment"),
    labels = clintrial_labels
)

plot5 <- autoforest(fit_result)
#> Recommended plot dimensions: width = 13.5 in, height = 5.0 in
# No need to pass data or labels - extracted from fit_result

# Save with recommended dimensions
dims <- attr(plot5, "rec_dims")
ggplot2::ggsave(file.path(tempdir(), "forest.pdf"),
                plot5, width = dims$width, height = dims$height)

# }
```
