# Feature Comparison

This article compares `summata` with other R packages for regression
table generation and summary statistics.

------------------------------------------------------------------------

## Similar Packages

The R ecosystem includes several well-established packages for creating
publication-ready regression tables. The following comparison identifies
areas of overlap and distinction between `summata` and its alternatives.

| Package | Primary Focus |
|:---|:---|
| **summata** | Fast regression workflows, forest plots, multivariate regression |
| **gtsummary** | Comprehensive table generation, `gt` ecosystem, maximum flexibility |
| **finalfit** | Clinical research, missing data handling, bootstrap simulations |
| **arsenal** | Large-scale summaries, SAS-like output |
| **stargazer** | Econometrics, LaTeX output, academic journal formatting |
| **tableone** | Simple Table 1 generation, SMD calculations |
| **compareGroups** | Bivariate analysis, clinical epidemiology |

------------------------------------------------------------------------

### Core Feature Matrix

| Feature | summata | gtsummary | finalfit | arsenal | stargazer | tableone | compareGroups |
|:---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| Descriptive Tables | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Stratified Summaries | ✓ | ✓ | ✓ | ✓ | — | ✓ | ✓ |
| Univariable Screening | ✓ | ✓ | ✓ | ✓ | — | — | ✓ |
| Multivariable Workflow | ✓ | ✓ | ✓ | ✓ | ✓ | — | — |
| Multivariate Regression | ✓ | — | — | — | — | — | — |
| Model Comparison | ✓ | ✓ | ◐ | — | ✓ | — | — |
| Mixed-Effects Models | ✓ | ✓ | ◐ | — | ◐ | — | — |
| Cox/Survival Models | ✓ | ✓ | ✓ | ✓ | ✓ | — | ✓ |
| Interaction Formatting | ✓ | ◐ | ◐ | ◐ | — | — | — |
| Forest Plots | ✓ | ◐ | ✓ | — | — | — | — |
| Table Merge/Stack | — | ✓ | ◐ | — | — | — | — |
| Export (Word/PDF) | ✓ | ✓ | ✓ | ✓ | ✓ | ◐ | ✓ |
| Variable Labels | ✓ | ✓ | ✓ | ✓ | ◐ | ✓ | ✓ |

**Legend:** ✓ Full support \| ◐ Partial support \| — Not available

------------------------------------------------------------------------

### Feature Definitions

| Feature | Description |
|:---|:---|
| **Descriptive Tables** | Summary statistics tables (mean, SD, median, IQR, *n*, %) |
| **Stratified Summaries** | Table 1-style summaries stratified by group with *p*-values |
| **Univariable Screening** | Test multiple predictors against one outcome (crude associations) |
| **Multivariable Workflow** | Combined univariable + multivariable analysis in one table |
| **Multivariate Regression** | Test one predictor across multiple outcomes |
| **Model Comparison** | Compare multiple models side-by-side with fit statistics |
| **Mixed-Effects Models** | Support for `lmer`/`glmer`/`coxme` random effects models |
| **Cox/Survival Models** | Cox proportional hazards and survival analysis |
| **Interaction Formatting** | Native support for interaction terms with formatted output |
| **Forest Plots** | Integrated forest plot visualization from regression results |
| **Table Merge/Stack** | Combine separate tables horizontally or vertically |
| **Export (Word/PDF)** | Direct export to publication formats |
| **Variable Labels** | Apply custom labels to variables in output |

------------------------------------------------------------------------

## Unique Strengths of *summata*

The following features distinguish `summata` from comparable packages:

### Multivariate Regression

The
[`multifit()`](https://phmcc.github.io/summata/reference/multifit.md)
and
[`multiforest()`](https://phmcc.github.io/summata/reference/multiforest.md)
functions implement an inverted screening paradigm: testing a single
predictor across multiple outcomes simultaneously. This workflow is
common in epidemiological and clinical research but not directly
supported by other packages.

``` r
# Test treatment effect across multiple outcomes
result <- multifit(
  data = clintrial,
  outcomes = c("surgery", "pfs_status", "os_status"),
  predictor = "treatment",
  covariates = c("age", "sex", "stage")
)
```

### Integrated Forest Plots

Customizable forest plots are generated directly from analysis results
without intermediate steps:

``` r
# From univariable screening
screen_result <- uniscreen(data, outcome, predictors)
uniforest(screen_result)

# From multivariate regression
multi_result <- multifit(data, outcomes, predictor)
multiforest(multi_result)
```

### Performance Optimization

Built on `data.table` for computational efficiency, `summata` provides
competitive performance across all benchmarked workflows. With
`conf_method = "wald"`, it is the fastest option tested for GLM-based
regression tables and univariable screening. With the default profile
likelihood CIs, performance is comparable to other packages that use the
same CI method (`finalfit`, `broom`). See the
[Benchmarks](https://phmcc.github.io/summata/articles/benchmarks.md)
article for detailed comparisons.

### Unified API

All modeling functions share consistent syntax across model types:

``` r
# Same interface for different model types
fit(data, outcome, predictors, model_type = "glm")
fit(data, outcome, predictors, model_type = "coxph")
fit(data, outcome, predictors, model_type = "lmer", random = "(1|site)")
```

### Complete Mixed-Effects Support

Full support for `coxme` (mixed-effects Cox models) alongside `lmer` and
`glmer`, which is limited or absent in other packages.

------------------------------------------------------------------------

## Additional Resources

- [Benchmarks](https://phmcc.github.io/summata/articles/benchmarks.md) —
  Performance comparisons
- [gtsummary documentation](https://www.danieldsjoberg.com/gtsummary/)
- [finalfit documentation](https://finalfit.org/)
- [arsenal documentation](https://mayoverse.github.io/arsenal/)
