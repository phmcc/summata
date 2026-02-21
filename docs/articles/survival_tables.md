# Survival Tables

Survival analysis requires specialized summary tables that report
time-to-event outcomes in formats appropriate for longitudinal research.
While
[`desctable()`](https://phmcc.github.io/summata/reference/desctable.md)
includes basic survival summaries (median with 95% CI), detailed
survival analysis often requires more comprehensive reporting: survival
probabilities at specified time points, multiple quantiles, and group
comparisons with appropriate statistical tests.

The
[`survtable()`](https://phmcc.github.io/summata/reference/survtable.md)
function generates publication-ready survival tables with flexible
output options. It uses the familiar
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) syntax for
specifying survival outcomes and adheres to the standard `summata`
calling convention:

``` r
survtable(data, outcome, by, times, probs, ...)
```

where `data` is the dataset, `outcome` specifies the survival endpoint
using [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) notation,
`by` defines the grouping variable, `times` specifies landmark time
points, and `probs` specifies survival quantiles to report.

------------------------------------------------------------------------

## Preliminaries

The examples in this vignette use the `clintrial` dataset included with
`summata`:

``` r
library(summata)
library(survival)

data(clintrial)
data(clintrial_labels)
```

The `clintrial` dataset contains 850 observations with time-to-event
variables suitable for demonstrating survival summaries.

------------------------------------------------------------------------

## Landmark Survival Estimates

Landmark survival estimates report the probability of survival at
specific time points (*e.g.*, 1-year, 2-year survival rates).

### **Example 1:** Basic Landmark Analysis

The following example reports survival probabilities at 12, 24, and 36
months (with the default median survival reporting included):

``` r
example1 <- survtable(
  data = clintrial,
  outcome = "Surv(os_months, os_status)",
  by = "treatment",
  times = c(12, 24, 36),
  time_unit = "months"
)

example1
#> 
#> Survival Summary Table
#> Outcome: Surv(os_months, os_status)
#> Stratified by: treatment
#> Time points: 12, 24, 36 months
#> Quantiles: 50%
#> Statistic: Survival probability
#> Test: Log-rank (p = < 0.001)
#> 
#>    Variable    12 months    24 months    36 months  Median (95% CI) p-value
#>      <char>       <char>       <char>       <char>           <char>  <char>
#> 1:    Total 59% (56-62%) 46% (43-49%) 36% (33-40%) 19.4 (16.2-23.4) < 0.001
#> 2:  Control 55% (49-63%) 37% (31-44%) 26% (20-33%) 14.7 (10.5-19.2)        
#> 3:   Drug A 68% (63-74%) 56% (51-62%) 48% (42-54%) 33.6 (24.5-42.2)        
#> 4:   Drug B 54% (49-59%) 43% (38-48%) 33% (28-38%) 14.7 (11.0-21.8)
```

### **Example 2:** Remove Median Survival Column

By default, a column showing median survival is displayed
(`probs = 0.5`). To remove it, set `probs = NULL`.

``` r
example2 <- survtable(
  data = clintrial,
  outcome = "Surv(os_months, os_status)",
  by = "treatment",
  times = c(12, 24),
  probs = NULL,
  time_unit = "months",
  total = FALSE
)

example2
#> 
#> Survival Summary Table
#> Outcome: Surv(os_months, os_status)
#> Stratified by: treatment
#> Time points: 12, 24 months
#> Statistic: Survival probability
#> Test: Log-rank (p = < 0.001)
#> 
#>    Variable    12 months    24 months p-value
#>      <char>       <char>       <char>  <char>
#> 1:  Control 55% (49-63%) 37% (31-44%) < 0.001
#> 2:   Drug A 68% (63-74%) 56% (51-62%)        
#> 3:   Drug B 54% (49-59%) 43% (38-48%)
```

------------------------------------------------------------------------

## Survival Quantiles

Survival quantiles report the time at which a specified proportion of
subjects have experienced the event. The median survival time (50th
percentile) is the most common, but other quantiles provide additional
context.

### **Example 3:** Survival Quartiles

Report multiple survival quantiles by specifying the `probs` parameter:

``` r
example3 <- survtable(
  data = clintrial,
  outcome = "Surv(os_months, os_status)",
  by = "stage",
  times = NULL,
  probs = c(0.25, 0.5, 0.75),
  labels = clintrial_labels
)

example3
#> 
#> Survival Summary Table
#> Outcome: Surv(os_months, os_status)
#> Stratified by: stage
#> Quantiles: 25%, 50%, 75%
#> Statistic: Survival probability
#> Test: Log-rank (p = < 0.001)
#> 
#>    Variable 75% Survival Time (95% CI)  Median (95% CI) 25% Survival Time (95% CI) p-value
#>      <char>                     <char>           <char>                     <char>  <char>
#> 1:    Total              5.1 (4.3-6.1) 19.5 (16.1-23.5)                         NR < 0.001
#> 2:        I             9.0 (7.4-12.5) 33.7 (22.4-45.2)                         NR        
#> 3:       II             7.7 (6.1-12.3) 32.2 (24.8-42.0)                         NR        
#> 4:      III              3.9 (3.0-5.2) 13.7 (10.4-22.3)             44.7 (35.8-NR)        
#> 5:       IV              2.0 (1.4-2.8)    5.9 (4.3-8.2)           15.6 (10.8-26.7)
```

------------------------------------------------------------------------

## Multiple Endpoints

Studies often include multiple time-to-event outcomes, such as
progression-free survival (PFS) and overall survival (OS). The
[`survtable()`](https://phmcc.github.io/summata/reference/survtable.md)
function handles multiple endpoints in a single call.

### **Example 4:** PFS and OS Comparison

Pass a vector of outcomes to compare multiple survival endpoints:

``` r
example4 <- survtable(
  data = clintrial,
  outcome = c("Surv(pfs_months, pfs_status)", "Surv(os_months, os_status)"),
  by = "treatment",
  times = c(12, 24),
  probs = 0.5,
  time_unit = "months",
  total = FALSE,
  labels = c(
    "Surv(pfs_months, pfs_status)" = "Progression-Free Survival",
    "Surv(os_months, os_status)" = "Overall Survival"
  )
)

example4
#> 
#> Survival Summary Table
#> Outcomes: 2
#> Stratified by: treatment
#> Time points: 12, 24 months
#> Quantiles: 50%
#> Statistic: Survival probability
#> Test: Log-rank
#> 
#>                     Variable   Group    12 months    24 months  Median (95% CI) p-value
#>                       <char>  <char>       <char>       <char>           <char>  <char>
#> 1: Progression-Free Survival Control 29% (23-36%)  14% (9-20%)    4.7 (3.6-6.9) < 0.001
#> 2:                            Drug A 44% (39-51%) 30% (25-36%)   9.0 (7.5-12.4)        
#> 3:                            Drug B 30% (26-36%) 18% (14-22%)    4.7 (3.7-6.2)        
#> 4:          Overall Survival Control 55% (49-63%) 37% (31-44%) 14.7 (10.5-19.2) < 0.001
#> 5:                            Drug A 68% (63-74%) 56% (51-62%) 33.6 (24.5-42.2)        
#> 6:                            Drug B 54% (49-59%) 43% (38-48%) 14.7 (11.0-21.8)
```

------------------------------------------------------------------------

## Output Variations

### **Example 5:** Cumulative Event Rates

For competing risks analyses or when reporting event rates rather than
survival probabilities, use `type = "risk"` to display cumulative
incidence (1 − survival):

``` r
example5 <- survtable(
  data = clintrial,
  outcome = "Surv(os_months, os_status)",
  by = "treatment",
  times = c(12, 24, 36),
  type = "risk",
  time_unit = "months"
)

example5
#> 
#> Survival Summary Table
#> Outcome: Surv(os_months, os_status)
#> Stratified by: treatment
#> Time points: 12, 24, 36 months
#> Quantiles: 50%
#> Statistic: Cumulative incidence
#> Test: Log-rank (p = < 0.001)
#> 
#>    Variable    12 months    24 months    36 months  Median (95% CI) p-value
#>      <char>       <char>       <char>       <char>           <char>  <char>
#> 1:    Total 41% (38-44%) 54% (51-57%) 64% (60-67%) 19.4 (16.2-23.4) < 0.001
#> 2:  Control 45% (37-51%) 63% (56-69%) 74% (67-80%) 14.7 (10.5-19.2)        
#> 3:   Drug A 32% (26-37%) 44% (38-49%) 52% (46-58%) 33.6 (24.5-42.2)        
#> 4:   Drug B 46% (41-51%) 57% (52-62%) 67% (62-72%) 14.7 (11.0-21.8)
```

### **Example 6:** Including At-Risk Counts

The number at risk at each time point provides context for the precision
of survival estimates:

``` r
example6 <- survtable(
  data = clintrial,
  outcome = "Surv(os_months, os_status)",
  by = "treatment",
  times = c(12, 24),
  stats = c("survival", "ci", "n_risk"),
  time_unit = "months",
  total = FALSE
)

example6
#> 
#> Survival Summary Table
#> Outcome: Surv(os_months, os_status)
#> Stratified by: treatment
#> Time points: 12, 24 months
#> Quantiles: 50%
#> Statistic: Survival probability
#> Test: Log-rank (p = < 0.001)
#> 
#>    Variable              12 months              24 months  Median (95% CI) p-value
#>      <char>                 <char>                 <char>           <char>  <char>
#> 1:  Control 55% (49-63%) [n = 106]  37% (31-44%) [n = 70] 14.7 (10.5-19.2) < 0.001
#> 2:   Drug A 68% (63-74%) [n = 192] 56% (51-62%) [n = 154] 33.6 (24.5-42.2)        
#> 3:   Drug B 54% (49-59%) [n = 190] 43% (38-48%) [n = 145] 14.7 (11.0-21.8)
```

------------------------------------------------------------------------

## Pairing with Kaplan-Meier Curves

Survival tables are often presented alongside Kaplan-Meier curves in
publications. While `summata` focuses on tabular output, the `survminer`
and `ggsurvfit` packages provide excellent options for survival curves.

### **Example 7:** Coordinated Figure-Table Output

The following workflow produces a matched figure-table pair suitable for
publication:

``` r
library(ggsurvfit)
library(survival)

# Fit the survival model
km_fit <- survfit(Surv(os_months, os_status) ~ treatment, data = clintrial)

# Create Kaplan-Meier plot with risk table
ggsurvfit(km_fit) +
  add_confidence_interval() +
  add_risktable() +
  add_quantile(y_value = 0.5, linetype = "dashed") +
  scale_ggsurvfit() +
  labs(
    title = "Overall Survival by Treatment",
    x = "Time (months)",
    y = "Survival Probability"
  ) +
  theme_minimal()

# Generate the companion table
surv_table <- survtable(
  data = clintrial,
  outcome = "Surv(os_months, os_status)",
  by = "treatment",
  times = c(12, 24, 36),
  probs = 0.5,
  time_unit = "months"
)

# Export both for publication
ggsave("km_curve.pdf", width = 8, height = 6)
table2pdf(surv_table, "survival_table.pdf", 
          caption = "Table 2. Survival Estimates by Treatment Group")
```

This workflow ensures that the Kaplan-Meier curve and survival table
report consistent time points and groupings.

------------------------------------------------------------------------

## Exporting Tables

Survival tables can be exported to various formats using the standard
`summata` export functions. See the [Table
Export](https://phmcc.github.io/summata/articles/table_export.md)
vignette for comprehensive documentation.

``` r
# Microsoft Word
table2docx(
  table = example1,
  file = "SurvivalTable.docx",
  caption = "Table 2. Survival Estimates by Treatment Group"
)

# PDF (requires LaTeX)
table2pdf(
  table = example1,
  file = "SurvivalTable.pdf",
  caption = "Table 2. Survival Estimates by Treatment Group"
)
```

------------------------------------------------------------------------

## Best Practices

### Time Point Selection

When selecting landmark time points, consider the following:

1.  Use clinically meaningful intervals (e.g., 1-year, 2-year for
    oncology)
2.  Ensure adequate follow-up at selected time points
3.  Match time points across treatment arms for comparability
4.  Consider the median follow-up when selecting the latest time point

### Reporting Recommendations

1.  Always report confidence intervals alongside point estimates
2.  Include the number at risk when space permits
3.  Specify the time unit in table headers or footnotes
4.  Use consistent decimal precision across estimates

------------------------------------------------------------------------

## Further Reading

- [Descriptive
  Tables](https://phmcc.github.io/summata/articles/descriptive_tables.md):
  [`desctable()`](https://phmcc.github.io/summata/reference/desctable.md)
  for baseline characteristics
- [Regression
  Modeling](https://phmcc.github.io/summata/articles/regression_modeling.md):
  [`fit()`](https://phmcc.github.io/summata/reference/fit.md),
  [`uniscreen()`](https://phmcc.github.io/summata/reference/uniscreen.md),
  and
  [`fullfit()`](https://phmcc.github.io/summata/reference/fullfit.md)
- [Model
  Comparison](https://phmcc.github.io/summata/articles/model_comparison.md):
  [`compfit()`](https://phmcc.github.io/summata/reference/compfit.md)
  for comparing models
- [Table
  Export](https://phmcc.github.io/summata/articles/table_export.md):
  Export to PDF, Word, and other formats
- [Forest
  Plots](https://phmcc.github.io/summata/articles/forest_plots.md):
  Visualization of regression results
- [Multivariate
  Regression](https://phmcc.github.io/summata/articles/multivariate_regression.md):
  [`multifit()`](https://phmcc.github.io/summata/reference/multifit.md)
  for multi-outcome analysis
- [Advanced
  Workflows](https://phmcc.github.io/summata/articles/advanced_workflows.md):
  Interactions and mixed-effects models
