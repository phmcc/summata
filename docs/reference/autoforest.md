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
    [`fit()`](https://phmcc.github.io/summata/reference/fit.md)

  - A `fullfit_result` object from
    [`fullfit()`](https://phmcc.github.io/summata/reference/fullfit.md)

  - A `uniscreen_result` object from
    [`uniscreen()`](https://phmcc.github.io/summata/reference/uniscreen.md)

  - A `multifit_result` object from
    [`multifit()`](https://phmcc.github.io/summata/reference/multifit.md)

- data:

  Data frame or data.table containing the original data. Required when
  `x` is a raw model object. Ignored when `x` is a result object from
  [`fit()`](https://phmcc.github.io/summata/reference/fit.md),
  [`fullfit()`](https://phmcc.github.io/summata/reference/fullfit.md),
  [`uniscreen()`](https://phmcc.github.io/summata/reference/uniscreen.md),
  or
  [`multifit()`](https://phmcc.github.io/summata/reference/multifit.md)
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
    [`uniforest()`](https://phmcc.github.io/summata/reference/uniforest.md).

2.  **multifit results**: Detected by presence of attributes
    `predictor`, `outcomes`, `model_type`, and `raw_data`. Routes to
    [`multiforest()`](https://phmcc.github.io/summata/reference/multiforest.md).

3.  **Cox models**: Classes `coxph` or `clogit`. Routes to
    [`coxforest()`](https://phmcc.github.io/summata/reference/coxforest.md).

4.  **GLM models**: Class `glm`. Routes to
    [`glmforest()`](https://phmcc.github.io/summata/reference/glmforest.md).

5.  **Linear models**: Class `lm` (but not `glm`). Routes to
    [`lmforest()`](https://phmcc.github.io/summata/reference/lmforest.md).

## See also

[`glmforest`](https://phmcc.github.io/summata/reference/glmforest.md)
for GLM forest plots,
[`coxforest`](https://phmcc.github.io/summata/reference/coxforest.md)
for Cox model forest plots,
[`lmforest`](https://phmcc.github.io/summata/reference/lmforest.md) for
linear model forest plots,
[`uniforest`](https://phmcc.github.io/summata/reference/uniforest.md)
for univariable screening forest plots,
[`multiforest`](https://phmcc.github.io/summata/reference/multiforest.md)
for multi-outcome forest plots,
[`fit`](https://phmcc.github.io/summata/reference/fit.md) for
single-model regression,
[`fullfit`](https://phmcc.github.io/summata/reference/fullfit.md) for
combined univariable/multivariable regression,
[`uniscreen`](https://phmcc.github.io/summata/reference/uniscreen.md)
for univariable screening,
[`multifit`](https://phmcc.github.io/summata/reference/multifit.md) for
multi-outcome analysis

Other visualization functions:
[`coxforest()`](https://phmcc.github.io/summata/reference/coxforest.md),
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

# Example 1: Logistic regression model
glm_model <- glm(surgery ~ age + sex + bmi + smoking,
                 family = binomial, data = clintrial)

plot1 <- autoforest(glm_model, data = clintrial)
print(plot1)
# Automatically detects GLM and routes to glmforest()

# Example 2: Cox proportional hazards model
cox_model <- coxph(Surv(os_months, os_status) ~ age + sex + treatment + stage,
                   data = clintrial)

plot2 <- autoforest(cox_model, data = clintrial)
print(plot2)
# Automatically detects coxph and routes to coxforest()

# Example 3: Linear regression model
lm_model <- lm(biomarker_x ~ age + sex + bmi + treatment, data = clintrial)

plot3 <- autoforest(lm_model, data = clintrial)
print(plot3)
# Automatically detects lm and routes to lmforest()

# Example 4: With custom labels
plot4 <- autoforest(
    glm_model,
    data = clintrial,
    labels = clintrial_labels,
    title = "Risk Factors for Surgical Intervention"
)
print(plot4)

# Example 5: Pass additional formatting options
plot5 <- autoforest(
    cox_model,
    data = clintrial,
    labels = clintrial_labels,
    zebra_stripes = TRUE,
    indent_groups = TRUE,
    color = "#E74C3C"
)
print(plot5)

# Example 6: From fit() result - data and labels extracted automatically
fit_result <- fit(
    data = clintrial,
    outcome = "surgery",
    predictors = c("age", "sex", "bmi", "treatment"),
    labels = clintrial_labels
)

plot6 <- autoforest(fit_result)
print(plot6)
# No need to pass data or labels - extracted from fit_result

# Example 7: From fullfit() result
ff_result <- fullfit(
    data = clintrial,
    outcome = "surgery",
    predictors = c("age", "sex", "bmi", "smoking", "treatment"),
    labels = clintrial_labels
)

plot7 <- autoforest(ff_result)
print(plot7)

# Example 8: From uniscreen() result
uni_result <- uniscreen(
    data = clintrial,
    outcome = "surgery",
    predictors = c("age", "sex", "bmi", "smoking", "treatment", "stage"),
    labels = clintrial_labels
)

plot8 <- autoforest(uni_result)
print(plot8)
# Automatically detects uniscreen_result and routes to uniforest()

# Example 9: From multifit() result
mf_result <- multifit(
    data = clintrial,
    outcomes = c("surgery", "any_complication", "readmission_30d"),
    predictor = "treatment",
    labels = clintrial_labels
)

plot9 <- autoforest(mf_result)
print(plot9)
# Automatically detects multifit result and routes to multiforest()

# Example 10: Survival uniscreen
surv_uni <- uniscreen(
    data = clintrial,
    outcome = "Surv(os_months, os_status)",
    predictors = c("age", "sex", "treatment", "stage", "grade"),
    model_type = "coxph",
    labels = clintrial_labels
)

plot10 <- autoforest(surv_uni)
print(plot10)
# Title automatically set to "Univariable Survival Analysis Screening"

# Example 11: Override auto-generated title
plot11 <- autoforest(
    surv_uni,
    title = "Univariable Predictors of Overall Survival"
)
print(plot11)

# Example 12: Save with recommended dimensions
plot12 <- autoforest(cox_model, data = clintrial, labels = clintrial_labels)
dims <- attr(plot12, "rec_dims")

# Save to file
# ggsave("forest_plot.pdf", plot12, width = dims$width, height = dims$height)
# ggsave("forest_plot.png", plot12, width = dims$width, height = dims$height, dpi = 300)

# Example 13: Poisson regression
clintrial$event_count <- rpois(nrow(clintrial), lambda = 3)
pois_model <- glm(event_count ~ age + sex + treatment,
                  family = poisson, data = clintrial)

plot13 <- autoforest(pois_model, data = clintrial)
print(plot13)
# Automatically detects Poisson GLM, title set to "Poisson Regression Model"

# Example 14: Mixed-Effects uniscreen
if (requireNamespace("lme4", quietly = TRUE)) {
    me_uni <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex", "treatment"),
        model_type = "glmer",
        random = "(1 | site)",
        labels = clintrial_labels
    )
    
    plot14 <- autoforest(me_uni)
    print(plot14)
    # Title: "Univariable Mixed-Effects Logistic Screening"
}

# Example 15: Quick comparison workflow
# Fit multiple model types and visualize each
models <- list(
    logistic = glm(surgery ~ age + sex + treatment, binomial, clintrial),
    survival = coxph(Surv(os_months, os_status) ~ age + sex + treatment, clintrial),
    linear = lm(biomarker_x ~ age + sex + treatment, clintrial)
)

# autoforest handles each appropriately
plots <- lapply(models, function(m) autoforest(m, data = clintrial))

# Each plot uses the correct forest function automatically
print(plots$logistic)
print(plots$survival)
print(plots$linear)
}
```
