# Descriptive Tables

Descriptive statistics provide the foundation for any quantitative
analysis. Before estimating relationships between variables, it is
essential to first characterize the distribution of each variable and
assess balance across comparison groups. A well-constructed descriptive
table—often designated “Table 1” in published research—accomplishes
three objectives: it summarizes the central tendency and dispersion of
continuous variables, tabulates the frequency distribution of
categorical variables, and tests for systematic differences between
groups.

The
[`desctable()`](https://phmcc.github.io/summata/reference/desctable.md)
function generates publication-ready descriptive tables with automatic
detection of variable types, appropriate summary statistics, and
optional hypothesis testing. It adheres to the to the standard `summata`
calling convention:

``` r
desctable(data, by, variables, ...)
```

where `data` is the dataset, `by` specifies the grouping variable
(optional), and `variables` lists the variables to summarize. This
vignette demonstrates the function’s capabilities using the included
sample dataset.

------------------------------------------------------------------------

## Preliminaries

The examples in this vignette use the `clintrial` dataset included with
`summata`:

``` r
library(summata)

data(clintrial)
data(clintrial_labels)
```

The `clintrial` dataset contains 850 observations with continuous,
categorical, and time-to-event variables suitable for demonstrating
descriptive statistics. The `clintrial_labels` vector provides
human-readable labels for display.

------------------------------------------------------------------------

## Summary Statistics and Tests

The
[`desctable()`](https://phmcc.github.io/summata/reference/desctable.md)
function automatically selects appropriate summary statistics and
hypothesis tests based on variable type:

| Variable Type              | Summary Statistic | Two Groups   | Three+ Groups  |
|:---------------------------|:------------------|:-------------|:---------------|
| Continuous (parametric)    | Mean ± SD         | *t*-test     | ANOVA          |
| Continuous (nonparametric) | Median \[IQR\]    | Wilcoxon     | Kruskal–Wallis |
| Categorical                | *n* (%)           | χ² or Fisher | χ² or Fisher   |
| Time-to-event              | Median (95% CI)   | Log-rank     | Log-rank       |

For categorical variables, Fisher exact test is used when any expected
cell count falls below 5. For continuous variables, the test selection
follows the displayed statistic: parametric tests are used with
mean-based statistics, nonparametric tests with median-based statistics.

------------------------------------------------------------------------

## Basic Usage

### **Example 1:** Grouped Descriptive Table

The most common use-case for descriptive tables is comparing
characteristics across groups:

``` r
desc_vars <- c("age", "sex", "race", "bmi", "stage", "ecog", 
               "Surv(os_months, os_status)")

example1 <- desctable(
  data = clintrial, 
  by = "treatment", 
  variables = desc_vars,
  labels = clintrial_labels
)

example1
#>                      Variable           Group            Total          Control           Drug A           Drug B p-value
#>                        <char>          <char>           <char>           <char>           <char>           <char>  <char>
#>  1:                         N                              850              196              292              362        
#>  2:               Age (years)    Median [IQR] 60.0 [52.0-67.0] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0]   0.465
#>  3:                       Sex          Female      450 (52.9%)      100 (51.0%)      164 (56.2%)      186 (51.4%)   0.394
#>  4:                                      Male      400 (47.1%)       96 (49.0%)      128 (43.8%)      176 (48.6%)        
#>  5:                      Race           White      598 (70.4%)      147 (75.0%)      198 (67.8%)      253 (69.9%)   0.794
#>  6:                                     Black      126 (14.8%)       25 (12.8%)       47 (16.1%)       54 (14.9%)        
#>  7:                                     Asian       93 (10.9%)        17 (8.7%)       35 (12.0%)       41 (11.3%)        
#>  8:                                     Other        33 (3.9%)         7 (3.6%)        12 (4.1%)        14 (3.9%)        
#>  9:   Body Mass Index (kg/m²)    Median [IQR] 28.0 [24.7-31.5] 28.0 [24.4-31.9] 28.2 [24.9-31.3] 27.8 [24.7-31.5]   0.831
#> 10:             Disease Stage               I      211 (24.9%)       59 (30.4%)       78 (26.7%)       74 (20.5%)   0.007
#> 11:                                        II      263 (31.1%)       65 (33.5%)       93 (31.8%)      105 (29.1%)        
#> 12:                                       III      241 (28.5%)       39 (20.1%)       75 (25.7%)      127 (35.2%)        
#> 13:                                        IV      132 (15.6%)       31 (16.0%)       46 (15.8%)       55 (15.2%)        
#> 14:   ECOG Performance Status               0      265 (31.5%)       56 (28.9%)       83 (28.5%)      126 (35.3%)   0.150
#> 15:                                         1      302 (35.9%)       79 (40.7%)      115 (39.5%)      108 (30.3%)        
#> 16:                                         2      238 (28.3%)       51 (26.3%)       79 (27.1%)      108 (30.3%)        
#> 17:                                         3        37 (4.4%)         8 (4.1%)        14 (4.8%)        15 (4.2%)        
#> 18: Overall Survival (months) Median (95% CI) 19.4 (16.2-23.4) 14.7 (10.5-19.2) 33.6 (24.5-42.2) 14.7 (11.0-21.8) < 0.001
```

The output includes a “Total” column by default, showing overall
statistics alongside group-specific values.

### **Example 2:** Ungrouped Summary Statistics

Omitting the `by` argument produces overall summary statistics without
group comparisons:

``` r
example2 <- desctable(
  data = clintrial,
  variables = c("age", "bmi", "sex", "stage"),
  labels = clintrial_labels
)

example2
#>                   Variable        Group            Total
#>                     <char>       <char>           <char>
#> 1:                       N                           850
#> 2:             Age (years) Median [IQR] 60.0 [52.0-67.0]
#> 3: Body Mass Index (kg/m²) Median [IQR] 28.0 [24.7-31.5]
#> 4:                     Sex       Female      450 (52.9%)
#> 5:                                 Male      400 (47.1%)
#> 6:           Disease Stage            I      211 (24.9%)
#> 7:                                   II      263 (31.1%)
#> 8:                                  III      241 (28.5%)
#> 9:                                   IV      132 (15.6%)
```

------------------------------------------------------------------------

## Customizing Summary Statistics

The default summary statistics can be customized for both continuous and
categorical variables.

### **Example 3:** Continuous Variables

The `stats_continuous` parameter controls how continuous variables are
summarized:

| Value            | Output Format              |
|:-----------------|:---------------------------|
| `"mean_sd"`      | Mean ± SD                  |
| `"median_iqr"`   | Median \[Q1–Q3\] (default) |
| `"median_range"` | Median (min–max)           |

``` r
example3 <- desctable(
  data = clintrial,
  by = "treatment",
  variables = c("age", "bmi", "los_days"),
  stats_continuous = c("mean_sd", "median_iqr", "median_range"),
  labels = clintrial_labels
)

example3
#>                           Variable          Group            Total          Control           Drug A           Drug B p-value
#>                             <char>         <char>           <char>           <char>           <char>           <char>  <char>
#>  1:                              N                             850              196              292              362        
#>  2:                    Age (years)      Mean ± SD      60.0 ± 11.8      59.2 ± 11.7      59.9 ± 12.1      60.5 ± 11.8   0.472
#>  3:                                  Median [IQR] 60.0 [52.0-67.0] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0]        
#>  4:                                Median (Range) 60.0 (18.0-90.0) 59.0 (26.0-90.0) 60.0 (18.0-90.0) 61.0 (24.0-90.0)        
#>  5:        Body Mass Index (kg/m²)      Mean ± SD       28.1 ± 4.9       28.0 ± 5.0       28.2 ± 4.9       28.1 ± 5.0   0.880
#>  6:                                  Median [IQR] 28.0 [24.7-31.5] 28.0 [24.4-31.9] 28.2 [24.9-31.3] 27.8 [24.7-31.5]        
#>  7:                                Median (Range) 28.0 (15.0-42.6) 28.0 (17.4-40.8) 28.2 (15.5-40.9) 27.8 (15.0-42.6)        
#>  8: Length of Hospital Stay (days)      Mean ± SD       20.1 ± 4.8       19.1 ± 4.5       18.4 ± 4.3       21.9 ± 4.8 < 0.001
#>  9:                                  Median [IQR] 19.9 [16.5-23.3] 19.0 [15.9-22.2] 18.1 [15.2-21.4] 21.7 [18.8-25.4]        
#> 10:                                Median (Range)  19.9 (6.5-36.8)  19.0 (6.5-29.0)  18.1 (8.3-32.2)  21.7 (8.4-36.8)
```

### **Example 4:** Categorical Variables

The `stats_categorical` parameter controls categorical variable display:

| Value         | Output Format     |
|:--------------|:------------------|
| `"n_percent"` | *n* (%) (default) |
| `"n"`         | *n* only          |
| `"percent"`   | % only            |

``` r
example4 <- desctable(
  data = clintrial,
  by = "treatment",
  variables = c("sex", "stage", "ecog"),
  stats_categorical = "percent",
  labels = clintrial_labels
)

example4
#>                    Variable  Group  Total Control Drug A Drug B p-value
#>                      <char> <char> <char>  <char> <char> <char>  <char>
#>  1:                       N           850     196    292    362        
#>  2:                     Sex Female  52.9%   51.0%  56.2%  51.4%   0.394
#>  3:                           Male  47.1%   49.0%  43.8%  48.6%        
#>  4:           Disease Stage      I  24.9%   30.4%  26.7%  20.5%   0.007
#>  5:                             II  31.1%   33.5%  31.8%  29.1%        
#>  6:                            III  28.5%   20.1%  25.7%  35.2%        
#>  7:                             IV  15.6%   16.0%  15.8%  15.2%        
#>  8: ECOG Performance Status      0  31.5%   28.9%  28.5%  35.3%   0.150
#>  9:                              1  35.9%   40.7%  39.5%  30.3%        
#> 10:                              2  28.3%   26.3%  27.1%  30.3%        
#> 11:                              3   4.4%    4.1%   4.8%   4.2%
```

### **Example 5:** Numeric Precision

Control decimal places with `digits` (for statistics) and `p_digits`
(for *p*-values):

``` r
example5 <- desctable(
  data = clintrial,
  by = "treatment",
  variables = c("age", "bmi", "sex"),
  digits = 2,
  p_digits = 4,
  test = TRUE,
  labels = clintrial_labels
)

example5
#>                   Variable        Group               Total             Control              Drug A              Drug B p-value
#>                     <char>       <char>              <char>              <char>              <char>              <char>  <char>
#> 1:                       N                              850                 196                 292                 362        
#> 2:             Age (years) Median [IQR] 60.00 [52.00-67.00] 59.00 [52.00-66.50] 60.00 [51.00-68.00] 61.00 [53.00-68.00]  0.4654
#> 3: Body Mass Index (kg/m²) Median [IQR] 28.00 [24.70-31.50] 28.00 [24.40-31.90] 28.25 [24.90-31.30] 27.80 [24.70-31.50]  0.8314
#> 4:                     Sex       Female         450 (52.9%)         100 (51.0%)         164 (56.2%)         186 (51.4%)  0.3943
#> 5:                                 Male         400 (47.1%)          96 (49.0%)         128 (43.8%)         176 (48.6%)
```

------------------------------------------------------------------------

## Statistical Testing

When comparing groups, hypothesis tests assess whether observed
differences are statistically significant.

### **Example 6:** Disabling Automatic Test Selection

By default, automatic hypothesis testing based on the summary statistic
is displayed (`test = TRUE`). Setting `test = FALSE` disables this
functionality:

``` r
example6 <- desctable(
  data = clintrial,
  by = "treatment",
  variables = c("age", "bmi", "sex", "stage"),
  test = FALSE,
  labels = clintrial_labels
)

example6
#>                   Variable        Group            Total          Control           Drug A           Drug B
#>                     <char>       <char>           <char>           <char>           <char>           <char>
#> 1:                       N                           850              196              292              362
#> 2:             Age (years) Median [IQR] 60.0 [52.0-67.0] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0]
#> 3: Body Mass Index (kg/m²) Median [IQR] 28.0 [24.7-31.5] 28.0 [24.4-31.9] 28.2 [24.9-31.3] 27.8 [24.7-31.5]
#> 4:                     Sex       Female      450 (52.9%)      100 (51.0%)      164 (56.2%)      186 (51.4%)
#> 5:                                 Male      400 (47.1%)       96 (49.0%)      128 (43.8%)      176 (48.6%)
#> 6:           Disease Stage            I      211 (24.9%)       59 (30.4%)       78 (26.7%)       74 (20.5%)
#> 7:                                   II      263 (31.1%)       65 (33.5%)       93 (31.8%)      105 (29.1%)
#> 8:                                  III      241 (28.5%)       39 (20.1%)       75 (25.7%)      127 (35.2%)
#> 9:                                   IV      132 (15.6%)       31 (16.0%)       46 (15.8%)       55 (15.2%)
```

### **Example 7:** Specifying Tests Manually

Override automatic selection with `test_continuous` and
`test_categorical`. Available test specifications include:

**Continuous** (`test_continuous`):

- `"auto"`: Automatic selection (default)
- `"t"`: Student *t*-test
- `"wrs"`: Wilcoxon rank-sum test
- `"aov"`: One-way ANOVA
- `"kwt"`: Kruskal–Wallis test

**Categorical** (`test_categorical`):

- `"auto"`: Automatic selection (default)
- `"chisq"`: Pearson χ² test
- `"fisher"`: Fisher exact test

The following example forces parametric tests for continuous variables:

``` r
example7a <- desctable(
  data = clintrial,
  by = "treatment",
  variables = c("age", "bmi", "los_days"),
  test = TRUE,
  test_continuous = "aov",  # ANOVA
  labels = clintrial_labels
)

example7a
#>                          Variable        Group            Total          Control           Drug A           Drug B p-value
#>                            <char>       <char>           <char>           <char>           <char>           <char>  <char>
#> 1:                              N                           850              196              292              362        
#> 2:                    Age (years) Median [IQR] 60.0 [52.0-67.0] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0]   0.472
#> 3:        Body Mass Index (kg/m²) Median [IQR] 28.0 [24.7-31.5] 28.0 [24.4-31.9] 28.2 [24.9-31.3] 27.8 [24.7-31.5]   0.880
#> 4: Length of Hospital Stay (days) Median [IQR] 19.9 [16.5-23.3] 19.0 [15.9-22.2] 18.1 [15.2-21.4] 21.7 [18.8-25.4] < 0.001
```

This example forces the Fisher exact test for categorical variables:

``` r
example7b <- desctable(
  data = clintrial,
  by = "treatment",
  variables = c("sex", "stage"),
  test = TRUE,
  test_categorical = "fisher",
  labels = clintrial_labels
)

example7b
#>         Variable  Group       Total     Control      Drug A      Drug B p-value
#>           <char> <char>      <char>      <char>      <char>      <char>  <char>
#> 1:             N                850         196         292         362        
#> 2:           Sex Female 450 (52.9%) 100 (51.0%) 164 (56.2%) 186 (51.4%)   0.394
#> 3:                 Male 400 (47.1%)  96 (49.0%) 128 (43.8%) 176 (48.6%)        
#> 4: Disease Stage      I 211 (24.9%)  59 (30.4%)  78 (26.7%)  74 (20.5%)        
#> 5:                   II 263 (31.1%)  65 (33.5%)  93 (31.8%) 105 (29.1%)        
#> 6:                  III 241 (28.5%)  39 (20.1%)  75 (25.7%) 127 (35.2%)        
#> 7:                   IV 132 (15.6%)  31 (16.0%)  46 (15.8%)  55 (15.2%)
```

------------------------------------------------------------------------

## Handling Missing Data

Missing values require special consideration in descriptive tables.
Options control whether missing values are displayed and how percentages
are calculated.

### **Example 8:** Including Missing Values

By default, missing values are excluded from calculations. Set
`na_include = TRUE` to display them as a separate category:

``` r
example8 <- desctable(
  data = clintrial,
  by = "treatment",
  variables = c("smoking", "diabetes"),
  na_include = TRUE,
  labels = clintrial_labels
)

example8
#>          Variable   Group       Total     Control      Drug A      Drug B p-value
#>            <char>  <char>      <char>      <char>      <char>      <char>  <char>
#> 1:              N                 850         196         292         362        
#> 2: Smoking Status   Never 337 (40.5%)  83 (43.5%) 123 (42.7%) 131 (37.0%)   0.450
#> 3:                 Former 311 (37.3%)  66 (34.6%) 101 (35.1%) 144 (40.7%)        
#> 4:                Current 185 (22.2%)  42 (22.0%)  64 (22.2%)  79 (22.3%)        
#> 5:                Unknown          17           5           4           8        
#> 6:       Diabetes      No 637 (76.4%) 152 (79.6%) 218 (75.7%) 267 (75.2%)   0.490
#> 7:                    Yes 197 (23.6%)  39 (20.4%)  70 (24.3%)  88 (24.8%)        
#> 8:                Unknown          16           5           4           7
```

### **Example 9:** Missing Value Denominators

The `na_percent` parameter controls whether missing values are included
in percentage calculations:

``` r
# Percentages exclude missing (denominator = non-missing)
example9a <- desctable(
  data = clintrial,
  by = "treatment",
  variables = c("smoking"),
  na_include = TRUE,
  na_percent = FALSE,
  labels = clintrial_labels
)

example9a
#>          Variable   Group       Total    Control      Drug A      Drug B p-value
#>            <char>  <char>      <char>     <char>      <char>      <char>  <char>
#> 1:              N                 850        196         292         362        
#> 2: Smoking Status   Never 337 (40.5%) 83 (43.5%) 123 (42.7%) 131 (37.0%)   0.450
#> 3:                 Former 311 (37.3%) 66 (34.6%) 101 (35.1%) 144 (40.7%)        
#> 4:                Current 185 (22.2%) 42 (22.0%)  64 (22.2%)  79 (22.3%)        
#> 5:                Unknown          17          5           4           8

# Percentages include missing (denominator = total)
example9b <- desctable(
  data = clintrial,
  by = "treatment",
  variables = c("smoking"),
  na_include = TRUE,
  na_percent = TRUE,
  labels = clintrial_labels
)

example9b
#>          Variable   Group       Total    Control      Drug A      Drug B p-value
#>            <char>  <char>      <char>     <char>      <char>      <char>  <char>
#> 1:              N                 850        196         292         362        
#> 2: Smoking Status   Never 337 (39.6%) 83 (42.3%) 123 (42.1%) 131 (36.2%)   0.450
#> 3:                 Former 311 (36.6%) 66 (33.7%) 101 (34.6%) 144 (39.8%)        
#> 4:                Current 185 (21.8%) 42 (21.4%)  64 (21.9%)  79 (21.8%)        
#> 5:                Unknown   17 (2.0%)   5 (2.6%)    4 (1.4%)    8 (2.2%)
```

### **Example 10:** Custom Missing Value Label

The label for missing values can be customized using the `na_label`
parameter:

``` r
example10 <- desctable(
  data = clintrial,
  by = "treatment",
  variables = c("smoking"),
  na_include = TRUE,
  na_label = "Not Reported",
  labels = clintrial_labels
)

example10
#>          Variable        Group       Total    Control      Drug A      Drug B p-value
#>            <char>       <char>      <char>     <char>      <char>      <char>  <char>
#> 1:              N                      850        196         292         362        
#> 2: Smoking Status        Never 337 (40.5%) 83 (43.5%) 123 (42.7%) 131 (37.0%)   0.450
#> 3:                      Former 311 (37.3%) 66 (34.6%) 101 (35.1%) 144 (40.7%)        
#> 4:                     Current 185 (22.2%) 42 (22.0%)  64 (22.2%)  79 (22.3%)        
#> 5:                Not Reported          17          5           4           8
```

------------------------------------------------------------------------

## Total Column Options

The total column provides overall statistics alongside group-specific
values.

### **Example 11:** Total Column Configuration

The `total` parameter controls the presence and position of the total
column:

| Value             | Effect                       |
|:------------------|:-----------------------------|
| `TRUE`, `"first"` | Total column first (default) |
| `"last"`          | Total column last            |
| `FALSE`           | No total column              |

``` r
# Total column in last position
example11a <- desctable(
  data = clintrial,
  by = "treatment",
  variables = c("age", "sex", "stage"),
  total = "last",
  labels = clintrial_labels
)

example11a
#>         Variable        Group          Control           Drug A           Drug B            Total p-value
#>           <char>       <char>           <char>           <char>           <char>           <char>  <char>
#> 1:             N                           196              292              362              850        
#> 2:   Age (years) Median [IQR] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0] 60.0 [52.0-67.0]   0.465
#> 3:           Sex       Female      100 (51.0%)      164 (56.2%)      186 (51.4%)      450 (52.9%)   0.394
#> 4:                       Male       96 (49.0%)      128 (43.8%)      176 (48.6%)      400 (47.1%)        
#> 5: Disease Stage            I       59 (30.4%)       78 (26.7%)       74 (20.5%)      211 (24.9%)   0.007
#> 6:                         II       65 (33.5%)       93 (31.8%)      105 (29.1%)      263 (31.1%)        
#> 7:                        III       39 (20.1%)       75 (25.7%)      127 (35.2%)      241 (28.5%)        
#> 8:                         IV       31 (16.0%)       46 (15.8%)       55 (15.2%)      132 (15.6%)

# No total column
example11b <- desctable(
  data = clintrial,
  by = "treatment",
  variables = c("age", "sex", "stage"),
  total = FALSE,
  labels = clintrial_labels
)

example11b
#>         Variable        Group          Control           Drug A           Drug B p-value
#>           <char>       <char>           <char>           <char>           <char>  <char>
#> 1:             N                           196              292              362        
#> 2:   Age (years) Median [IQR] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0]   0.465
#> 3:           Sex       Female      100 (51.0%)      164 (56.2%)      186 (51.4%)   0.394
#> 4:                       Male       96 (49.0%)      128 (43.8%)      176 (48.6%)        
#> 5: Disease Stage            I       59 (30.4%)       78 (26.7%)       74 (20.5%)   0.007
#> 6:                         II       65 (33.5%)       93 (31.8%)      105 (29.1%)        
#> 7:                        III       39 (20.1%)       75 (25.7%)      127 (35.2%)        
#> 8:                         IV       31 (16.0%)       46 (15.8%)       55 (15.2%)
```

------------------------------------------------------------------------

## Complete Example

The following demonstrates a comprehensive descriptive table suitable
for publication:

``` r
table1 <- desctable(
  data = clintrial,
  by = "treatment",
  variables = c(
    "age", "sex", "race", "ethnicity", "bmi",
    "smoking", "diabetes", "hypertension",
    "stage", "grade", "ecog",
    "Surv(os_months, os_status)"
  ),
  labels = clintrial_labels,
  stats_continuous = "mean_sd",
  stats_categorical = "n_percent",
  test = TRUE,
  total = TRUE,
  digits = 1,
  p_digits = 3
)

table1
#>                      Variable                     Group            Total          Control           Drug A           Drug B p-value
#>                        <char>                    <char>           <char>           <char>           <char>           <char>  <char>
#>  1:                         N                                        850              196              292              362        
#>  2:               Age (years)                 Mean ± SD      60.0 ± 11.8      59.2 ± 11.7      59.9 ± 12.1      60.5 ± 11.8   0.472
#>  3:                       Sex                    Female      450 (52.9%)      100 (51.0%)      164 (56.2%)      186 (51.4%)   0.394
#>  4:                                                Male      400 (47.1%)       96 (49.0%)      128 (43.8%)      176 (48.6%)        
#>  5:                      Race                     White      598 (70.4%)      147 (75.0%)      198 (67.8%)      253 (69.9%)   0.794
#>  6:                                               Black      126 (14.8%)       25 (12.8%)       47 (16.1%)       54 (14.9%)        
#>  7:                                               Asian       93 (10.9%)        17 (8.7%)       35 (12.0%)       41 (11.3%)        
#>  8:                                               Other        33 (3.9%)         7 (3.6%)        12 (4.1%)        14 (3.9%)        
#>  9:                 Ethnicity              Non-Hispanic      744 (87.5%)      172 (87.8%)      253 (86.6%)      319 (88.1%)   0.846
#> 10:                                            Hispanic      106 (12.5%)       24 (12.2%)       39 (13.4%)       43 (11.9%)        
#> 11:   Body Mass Index (kg/m²)                 Mean ± SD       28.1 ± 4.9       28.0 ± 5.0       28.2 ± 4.9       28.1 ± 5.0   0.880
#> 12:            Smoking Status                     Never      337 (40.5%)       83 (43.5%)      123 (42.7%)      131 (37.0%)   0.450
#> 13:                                              Former      311 (37.3%)       66 (34.6%)      101 (35.1%)      144 (40.7%)        
#> 14:                                             Current      185 (22.2%)       42 (22.0%)       64 (22.2%)       79 (22.3%)        
#> 15:                  Diabetes                        No      637 (76.4%)      152 (79.6%)      218 (75.7%)      267 (75.2%)   0.490
#> 16:                                                 Yes      197 (23.6%)       39 (20.4%)       70 (24.3%)       88 (24.8%)        
#> 17:              Hypertension                        No      504 (60.4%)      108 (56.2%)      179 (62.2%)      217 (61.1%)   0.401
#> 18:                                                 Yes      331 (39.6%)       84 (43.8%)      109 (37.8%)      138 (38.9%)        
#> 19:             Disease Stage                         I      211 (24.9%)       59 (30.4%)       78 (26.7%)       74 (20.5%)   0.007
#> 20:                                                  II      263 (31.1%)       65 (33.5%)       93 (31.8%)      105 (29.1%)        
#> 21:                                                 III      241 (28.5%)       39 (20.1%)       75 (25.7%)      127 (35.2%)        
#> 22:                                                  IV      132 (15.6%)       31 (16.0%)       46 (15.8%)       55 (15.2%)        
#> 23:               Tumor Grade       Well-differentiated      153 (18.2%)       42 (21.8%)       54 (18.7%)       57 (15.9%)   0.331
#> 24:                           Moderately differentiated      412 (49.0%)       96 (49.7%)      143 (49.5%)      173 (48.3%)        
#> 25:                               Poorly differentiated      275 (32.7%)       55 (28.5%)       92 (31.8%)      128 (35.8%)        
#> 26:   ECOG Performance Status                         0      265 (31.5%)       56 (28.9%)       83 (28.5%)      126 (35.3%)   0.150
#> 27:                                                   1      302 (35.9%)       79 (40.7%)      115 (39.5%)      108 (30.3%)        
#> 28:                                                   2      238 (28.3%)       51 (26.3%)       79 (27.1%)      108 (30.3%)        
#> 29:                                                   3        37 (4.4%)         8 (4.1%)        14 (4.8%)        15 (4.2%)        
#> 30: Overall Survival (months)           Median (95% CI) 19.4 (16.2-23.4) 14.7 (10.5-19.2) 33.6 (24.5-42.2) 14.7 (11.0-21.8) < 0.001
#>                      Variable                     Group            Total          Control           Drug A           Drug B p-value
#>                        <char>                    <char>           <char>           <char>           <char>           <char>  <char>
```

### Accessing Raw Data

The underlying numeric values are stored as an attribute for
programmatic access:

``` r
raw_data <- attr(table1, "raw_data")
head(raw_data)
#>       Variable   Group stat_type     Total Total_sd Total_n   Control Control_sd Control_n    Drug A Drug A_sd Drug A_n    Drug B Drug B_sd Drug B_n   p_value Total_total
#>         <char>  <char>    <char>     <num>    <num>   <int>     <num>      <num>     <int>     <num>     <num>    <int>     <num>     <num>    <int>     <num>       <int>
#> 1: Age (years) mean_sd   mean_sd  59.96588 11.84967     850  59.18367   11.69821       196  59.87671  12.07559      292  60.46133  11.75488      362 0.4720805          NA
#> 2:         Sex  Female  category 450.00000       NA      NA 100.00000         NA        NA 164.00000        NA       NA 186.00000        NA       NA 0.3942586         850
#> 3:                Male  category 400.00000       NA      NA  96.00000         NA        NA 128.00000        NA       NA 176.00000        NA       NA        NA         850
#> 4:        Race   White  category 598.00000       NA      NA 147.00000         NA        NA 198.00000        NA       NA 253.00000        NA       NA 0.7939277         850
#> 5:               Black  category 126.00000       NA      NA  25.00000         NA        NA  47.00000        NA       NA  54.00000        NA       NA        NA         850
#> 6:               Asian  category  93.00000       NA      NA  17.00000         NA        NA  35.00000        NA       NA  41.00000        NA       NA        NA         850
#>    Control_total Drug A_total Drug B_total Total_ci_lower Total_ci_upper Control_ci_lower Control_ci_upper Drug A_ci_lower Drug A_ci_upper Drug B_ci_lower Drug B_ci_upper
#>            <int>        <int>        <int>          <num>          <num>            <num>            <num>           <num>           <num>           <num>           <num>
#> 1:            NA           NA           NA             NA             NA               NA               NA              NA              NA              NA              NA
#> 2:           196          292          362             NA             NA               NA               NA              NA              NA              NA              NA
#> 3:           196          292          362             NA             NA               NA               NA              NA              NA              NA              NA
#> 4:           196          292          362             NA             NA               NA               NA              NA              NA              NA              NA
#> 5:           196          292          362             NA             NA               NA               NA              NA              NA              NA              NA
#> 6:           196          292          362             NA             NA               NA               NA              NA              NA              NA              NA
```

------------------------------------------------------------------------

## Survival Summary Tables

For detailed survival analysis—including landmark survival estimates,
survival quantiles, and multiple endpoints—see the dedicated [Survival
Tables](https://phmcc.github.io/summata/articles/survival_tables.md)
vignette. The
[`survtable()`](https://phmcc.github.io/summata/reference/survtable.md)
function provides comprehensive options for reporting time-to-event
outcomes.

------------------------------------------------------------------------

## Exporting Tables

Descriptive tables can be exported to various formats. See the [Table
Export](https://phmcc.github.io/summata/articles/table_export.md)
vignette for comprehensive documentation.

``` r
# Microsoft Word
table2docx(
  table = table1,
  file = "Table1.docx",
  caption = "Table 1. Baseline Characteristics by Group"
)

# PDF (requires LaTeX)
table2pdf(
  table = table1,
  file = "Table1.pdf",
  caption = "Table 1. Baseline Characteristics by Group"
)

# HTML
table2html(
  table = table1,
  file = "Table1.html",
  caption = "Table 1. Baseline Characteristics by Group"
)
```

------------------------------------------------------------------------

## Best Practices

### Variable Selection

1.  Include all relevant baseline characteristics
2.  Order variables logically (typically chronologically and by domain)
3.  Exclude the grouping variable from the variables list

### Statistical Considerations

1.  Use automatic test selection unless there is specific justification
    otherwise
2.  Report exact *p*-values when possible; very small values display as
    “\< 0.001” (or to preferred degree of precision)
3.  Consider multiple comparison adjustments when testing many variables
4.  For skewed or non-normally distributed continuous variables, report
    using nonparametric statistical procedures (e.g., median with IQR)
    rather than parametric ones (e.g., mean with SD)

### Formatting Recommendations

1.  Use consistent decimal precision within variable types
2.  Include units in variable labels (e.g., “Age (years)”)
3.  Include a total column for context
4.  Use landscape orientation for tables with many columns

------------------------------------------------------------------------

## Common Issues

### Empty Cells in Categorical Variables

When a factor level has zero observations in a group, ensure all levels
are explicitly defined:

``` r
data$stage <- factor(data$stage, levels = c("I", "II", "III", "IV"))
```

### Skewed Continuous Variables

For highly skewed distributions, use median and IQR:

``` r
desctable(data, by, variables, stats_continuous = "median_iqr")
```

### Large Tables

For tables with many variables, consider splitting by category or using
landscape orientation for export:

``` r
table2pdf(table, "table1.pdf", orientation = "landscape", font_size = 8)
```

------------------------------------------------------------------------

## Further Reading

- [Survival
  Tables](https://phmcc.github.io/summata/articles/survival_tables.md):
  [`survtable()`](https://phmcc.github.io/summata/reference/survtable.md)
  for time-to-event summaries
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
