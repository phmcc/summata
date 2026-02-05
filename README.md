# <span class="pkg-name">summata</span> <a href="https://phmcc.github.io/summata/"><img src="man/figures/summata.svg" align="right" height="139" alt="summata website" /></a>

<!-- badges: start -->
[![R-CMD-check](https://github.com/phmcc/summata/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/phmcc/summata/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/phmcc/summata/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/phmcc/summata/actions/workflows/test-coverage.yaml)
[![lint](https://github.com/phmcc/summata/actions/workflows/lint.yaml/badge.svg)](https://github.com/phmcc/summata/actions/workflows/lint.yaml)
[![Codecov test coverage](https://codecov.io/gh/phmcc/summata/branch/main/graph/badge.svg)](https://codecov.io/gh/phmcc/summata?branch=main)
[![CRAN status](https://www.r-pkg.org/badges/version/summata)](https://CRAN.R-project.org/package=summata)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

> ***summata*** | /suːˈmɑːtə/ | *Latin, n. pl. of* summātum*, gerundive of* summāre*: those that have been summarized*
>
> Concise, publication-ready statistical summaries.

## Overview

The `summata` package provides a comprehensive framework for generating summary tables and visualizations from statistical analyses. Built on `data.table` for computational efficiency, it streamlines the workflow from descriptive statistics through regression modeling to final output—all using a unified interface with standardized, presentation-ready results.

<img src="man/figures/README_coxforest.png" alt="Cox regression forest plot" width="100%">

## Installation

Install this package from GitHub (stable) or Codeberg (development):

```r
# Stable release
devtools::install_github("phmcc/summata")

# Development version
devtools::install_git("https://codeberg.org/phmcc/summata.git")
```

## Package Composition

### Design Principles

The architecture of `summata` reflects three guiding principles:

1. **Consistent syntax.** All modeling functions share a common signature: data first, followed by the outcome or variable of interest, then additional covariates. This convention facilitates both pipe-based workflows and pedagogical clarity.

2. **Transparent computation.** Functions attach their underlying model objects and raw numerical results as attributes, permitting verification of computations or extension of analyses beyond the formatted output.
   
3. **Separation of concerns.** Analysis, formatting, and export constitute distinct operations, allowing each stage to be modified independently.

These principles manifest in the standard calling convention:

``` r
result <- fit(data, variable, c("covar1", "covar2", ..., "covarN"), ...)
```

where `data` is the dataset, `variable` is the variable of interest (dependent variable, endpoint/outcome, stratification/grouping variable, *etc.*), and `c("covar1", "covar2", ..., "covarN")`is a vector of covariates (independent variables, model predictors, compared/explanatory factors, *etc.*). The result is a formatted table for export, with readily accessible attributes for further analysis.

``` r
# The formatted table
print(result)

# The underlying model object
model <- attr(result, "model")

# Raw numerical results
raw <- attr(result, "raw_data")
```

### Functional Reference

This package provides a variety of functions for different stages of statistical analysis:

#### Descriptive analysis

Tables to provide quick summary data and comparison tests between different groups.

| Function | Purpose |
|:---------|:--------|
| `desctable()` | Descriptive statistics with stratification and hypothesis testing |
| `survtable()` | Survival probability estimates at specified time points |

#### Predictive analysis

Fitted univariable and multivariable regression results for predictive modeling.

| Function | Purpose |
|:---------|:--------|
| `uniscreen()` | Systematic univariable analysis across multiple predictors |
| `fit()` | Single regression model with formatted coefficient extraction |
| `fullfit()` | Integrated univariable screening with multivariable regression |
| `compfit()` | Nested model comparison with composite scoring |
| `multifit()` | Multivariate analysis with a single predictor evaluated against multiple outcomes |

#### Table export

Export of finalized tables to various commonly used formats.

| Function | Format | Dependencies |
|:---------|:-------|:-----------|
| `autoexport()` | Auto-detect from file extension | Varies |
| `table2pdf()` | PDF | `xtable`, LaTeX distribution |
| `table2html()` | HTML | `xtable` |
| `table2tex()` | LaTeX source | `xtable` |
| `table2docx()` | Microsoft Word | `officer`, `flextable` |
| `table2pptx()` | Microsoft PowerPoint | `officer`, `flextable` |
| `table2rtf()` | Rich Text Format | `officer`, `flextable` |

#### Data visualization

Generation of publication-ready forest plot graphics to summarize regression models.

| Function | Application |
|:---------|:------------|
| `autoforest()` | Automatic model class detection |
| `lmforest()` | Linear models |
| `glmforest()` | Generalized linear models |
| `coxforest()` | Proportional hazards models |
| `uniforest()` | Univariable screening results |
| `multiforest()` | Multivariate analysis results |

### Supported Model Classes

The following regression models are currently supported by `summata`. Specify the model using the `model_type` parameter in the appropriate regression function (`uniscreen()`, `fit()`, `fullfit()`, `compfit()`, or `multifit()`):

| Model Class | `model_type` | Function | Effect Measure |
|:------------|:----------|:---------|:---------------|
| Linear regression | `lm` | `stats::lm()` | *β* coefficient |
| Logistic regression | `glm`, `family = "binomial"` | `stats::glm()` | Odds ratio |
| Logistic (overdispersed) | `glm`, `family = "quasibinomial"` | `stats::glm()` | Odds ratio |
| Poisson regression | `glm`, `family = "poisson"` | `stats::glm()` | Rate ratio |
| Poisson (overdispersed) | `glm`, `family = "quasipoisson"` | `stats::glm()` | Rate ratio |
| Gaussian (via GLM) | `glm`, `family = "gaussian"` | `stats::glm()` | *β* coefficient |
| Gamma regression | `glm`, `family = "Gamma"` | `stats::glm()` | Ratio* |
| Inverse Gaussian | `glm`, `family = "inverse.gaussian"` | `stats::glm()` | Ratio* |
| Cox proportional hazards | `coxph` | `survival::coxph()` | Hazard ratio |
| Conditional logistic | `clogit` | `survival::clogit()` | Odds ratio |
| Negative binomial | `negbin` | `MASS::glm.nb()` | Rate ratio |
| Linear mixed effects | `lmer` | `lme4::lmer()` | *β* coefficient |
| Generalized linear mixed effects | `glmer` | `lme4::glmer()` | Odds/rate ratio |
| Cox mixed effects | `coxme` | `coxme::coxme()` | Hazard ratio |

*with log link; coefficient with identity link

## Comparison with Related Packages

The R ecosystem includes several established packages for regression table generation. The following comparison identifies areas of overlap and distinction:

| Capability | summata | gtsummary | finalfit | arsenal |
|:-----------|:-------:|:---------:|:--------:|:-------:|
| Descriptive statistics | ✓ | ✓ | ✓ | ✓ |
| Survival summaries | ✓ | ✓ | ◐ | — |
| Univariable screening | ✓ | ✓ | ✓ | ✓ |
| Multivariable regression | ✓ | ✓ | ✓ | ✓ |
| Multi-format export | ✓ | ✓ | ✓ | ✓ |
| Integrated forest plots | ✓ | ◐ | ✓ | — |
| Model comparison | ✓ | ✓ | ◐ | — |
| Mixed-effects models | ✓ | ◐ | ◐ | — |
| Multivariate analysis | ✓ | — | — | — |

<sub>✓ Full support | ◐ Partial support | — Not available</sub>

A detailed feature comparison is available in the [package documentation](articles/benchmarks.html).

## Illustrative Example

The `clintrial` dataset included with this package provides simulated clinical trial data comprising patient identifiers, baseline characteristics, therapeutic interventions, short-term outcomes, and long-term survival statistics. The following example demonstrates how `summata` functions can be used to analyze perioperative factors affecting 30-day hospital readmission.

### **Step 0:** Data Preparation

Prior to analysis, load the dataset, apply labels, and define predictors:

``` r
library(summata)

# Load example data ≥
data("clintrial")
data("clintrial_labels")

# Define candidate predictors
predictors <- c("age", "sex", "race", "ethnicity", "bmi", "smoking",
                "hypertension", "diabetes", "ecog", "creatinine",
                "hemoglobin", "biomarker_x", "biomarker_y", "grade",
                "treatment", "surgery", "los_days", "stage")
```

### **Step 1:** Descriptive Statistics

Use the `desctable()` function to generate summary statistics with stratification by a grouping variable (in this case, 30-day readmission):

``` r
table1 <- desctable(
    data = clintrial,
    by = "readmission_30d",
    variables = predictors,
    labels = clintrial_labels
)

table2pdf(table1, "table1.pdf",
          caption = "\\textbf{Table 1} - Baseline characteristics by 30-day readmission status",
          paper = "auto",
          condense_table = TRUE,
          dark_header = TRUE,
          zebra_stripes = TRUE
)
```

<img src="man/figures/README_desctable.png" alt="Descriptive statistics table" width="100%">

### **Step 2:** Regression Analysis

Perform an integrated univariable-to-multivariable regression workflow using the `fullfit()` function:

``` r
table2 <- fullfit(
    data = clintrial,
    outcome = "readmission_30d",
    predictors = predictors,
    method = "screen",
    p_threshold = 0.05,
    model_type = "glm",
    labels = clintrial_labels
)

table2pdf(table2, "table2.pdf",
          caption = "\\textbf{Table 2} - Predictors of 30-day readmission, multivariable Cox regression with univariable screen",
          paper = "auto",
          condense_table = TRUE,
          dark_header = TRUE,
          zebra_stripes = TRUE
)
```

<img src="man/figures/README_fullfit.png" alt="Regression table" width="100%">

### **Step 3:** Forest Plot

Finally, generate a forest plot to provide a graphical representation of effect estimates using the `glmforest()` function:

``` r
forest_30d <- glmforest(table2,
                        title = "Factors Affecting 30-Day Readmission",
                        labels = clintrial_labels,
                        indent_groups = TRUE
                        )

ggsave("forest_30d.pdf", forest_30d,
       width = attr(forest_30d, "recommended_dims")$width,
       height = attr(forest_30d, "recommended_dims")$height, 
       units = "in")
```

<img src="man/figures/README_glmforest.png" alt="Forest plot" width="100%">

## Development

### Repository

- **Primary development**: [codeberg.org/phmcc/summata](https://codeberg.org/phmcc/summata)
- **GitHub releases**: [github.com/phmcc/summata](https://github.com/phmcc/summata)

### Contributing

Bug reports and feature requests may be submitted via the [issue tracker](https://github.com/phmcc/summata/issues). Contributions are welcome; please consult the contributing guidelines prior to submitting pull requests.

## Acknowledgments

The design of `summata` draws inspiration from several existing packages:

- **finalfit** (Harrison) — Regression workflow concepts  
- **gtsummary** (Sjoberg et al.) — Table generation architecture  
- **arsenal** (Heinzen et al.) — Descriptive statistics methodology  
- **data.table** (Dowle & Srinivasan) — High-performance data operations

## License

GPL ≥ 3.0

## Citation

``` r
citation("summata")

To cite summata in publications, use:

  McClelland PH (2026). _summata: Publication-Ready Summary Tables and Forest Plots_. R package version 1.0.0, <https://phmcc.github.io/summata/>.

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {summata: Publication-Ready Summary Tables and Forest Plots},
    author = {Paul Hsin-ti McClelland},
    year = {2026},
    note = {R package version 1.0.0},
    url = {https://phmcc.github.io/summata/},
  }
```

## Further Resources

- **Function documentation**: `?function_name` or the [reference index](reference/index.html)
- **Vignettes**: `vignette("summata")` or [online articles](articles/index.html)
- **Issue tracker**: [GitHub Issues](https://github.com/phmcc/summata/issues)

---

<sub>The `summata` package is under active development. The API may change prior to CRAN submission.</sub>
