# Create Publication-Ready Descriptive Statistics Tables

Generates comprehensive descriptive statistics tables with automatic
variable type detection, group comparisons, and appropriate statistical
testing. This function is designed to create "Table 1"-style summaries
commonly used in clinical and epidemiological research, with full
support for continuous, categorical, and time-to-event variables.

## Usage

``` r
desctable(
  data,
  by = NULL,
  variables,
  stats_continuous = c("median_iqr"),
  stats_categorical = "n_percent",
  digits = 1,
  p_digits = 3,
  conf_level = 0.95,
  p_per_stat = FALSE,
  na_include = FALSE,
  na_label = "Unknown",
  na_percent = FALSE,
  test = TRUE,
  test_continuous = "auto",
  test_categorical = "auto",
  total = TRUE,
  total_label = "Total",
  labels = NULL,
  number_format = NULL,
  ...
)
```

## Arguments

- data:

  Data frame or data.table containing the dataset to summarize.
  Automatically converted to a data.table for efficient processing.

- by:

  Character string specifying the column name of the grouping variable
  for stratified analysis (*e.g.*, treatment arm, exposure status). When
  `NULL` (default), produces overall summaries only without group
  comparisons or statistical tests.

- variables:

  Character vector of variable names to summarize. Can include standard
  column names for continuous or categorical variables, and survival
  expressions using
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) syntax (*e.g.*,
  `"Surv(os_months, os_status)"`). Variables are processed in the order
  provided.

- stats_continuous:

  Character vector specifying which statistics to compute for continuous
  variables. Multiple values create separate rows for each variable.
  Options:

  - `"mean_sd"` - Mean ± standard deviation

  - `"median_iqr"` - Median \[interquartile range\]

  - `"median_range"` - Median (minimum–maximum)

  - `"range"` - Minimum–maximum only

  Default is `"median_iqr"`.

- stats_categorical:

  Character string specifying the format for categorical variable
  summaries:

  - `"n"` - Count only

  - `"percent"` - Percentage only

  - `"n_percent"` - Count (percentage) \[default\]

- digits:

  Integer specifying the number of decimal places for continuous
  statistics. Default is 1.

- p_digits:

  Integer specifying the number of decimal places for *p*-values. Values
  smaller than `10^(-p_digits)` are displayed as `"< 0.001"` (for
  `p_digits = 3`), `"< 0.0001"` (for `p_digits = 4`), etc. Default is 3.

- conf_level:

  Numeric confidence level for confidence intervals in survival variable
  summaries (median survival time with CI). Must be between 0 and 1.
  Default is 0.95 (95% confidence intervals).

- p_per_stat:

  Logical. If `TRUE`, displays *p*-values on each row (per statistic)
  rather than only on the first row of each variable. Useful when
  different statistics within a variable warrant separate significance
  testing. Default is `FALSE`.

- na_include:

  Logical. If `TRUE`, missing values (NAs) are displayed as a separate
  category/row for each variable. If `FALSE`, missing values are
  silently excluded from calculations. Default is `FALSE`.

- na_label:

  Character string used to label the missing values row when
  `na_include = TRUE`. Default is `"Unknown"`.

- na_percent:

  Logical. Controls how percentages are calculated for categorical
  variables when `na_include = TRUE`:

  - If `TRUE`, percentages include NAs in the denominator (all
    categories sum to 100%)

  - If `FALSE`, percentages exclude NAs from the denominator
    (non-missing categories sum to 100%, missing shown separately)

  Only affects categorical variables. Default is `FALSE`.

- test:

  Logical. If `TRUE`, performs appropriate statistical tests for group
  comparisons and adds a *p*-value column. Requires `by` to be
  specified. Tests are automatically selected based on variable type and
  test parameters. Default is `TRUE`.

- test_continuous:

  Character string specifying the statistical test for continuous
  variables:

  - `"auto"` - Automatic selection: *t*-test/ANOVA for means,
    Wilcoxon/Kruskal-Wallis for medians \[default\]

  - `"t"` - Independent samples *t*-test (2 groups only)

  - `"aov"` - One-way ANOVA (2+ groups)

  - `"wrs"` - Wilcoxon rank-sum test (2 groups only)

  - `"kwt"` - Kruskal-Wallis test (2+ groups)

- test_categorical:

  Character string specifying the statistical test for categorical
  variables:

  - `"auto"` - Automatic selection: Fisher exact test if any expected
    cell frequency \< 5, otherwise χ² test \[default\]

  - `"fisher"` - Fisher exact test

  - `"chisq"` - χ² test

- total:

  Logical or character string controlling the total column:

  - `TRUE` or `"first"` - Include total column as first column after
    Variable/Group \[default\]

  - `"last"` - Include total column as last column before *p*-value

  - `FALSE` - Exclude total column

- total_label:

  Character string for the total column header. Default is `"Total"`.

- labels:

  Named character vector or list providing custom display labels for
  variables. Names should match variable names (or
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) expressions),
  values are the display labels. Variables not in `labels` use their
  original names. Can also label the grouping variable specified in
  `by`. Default is `NULL`.

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
        

- ...:

  Additional arguments passed to the underlying statistical test
  functions (*e.g.*, `var.equal = TRUE` for *t*-tests,
  `simulate.p.value = TRUE` for Fisher test).

## Value

A data.table with S3 class `"desctable"` containing formatted
descriptive statistics. The table structure includes:

- Variable:

  Variable name or label (from `labels`)

- Group:

  For continuous variables: statistic type (*e.g.*, "Mean ± SD", "Median
  \[IQR\]"). For categorical variables: category level. Empty for
  variable name rows.

- Total:

  Statistics for the total sample (if `total = TRUE`)

- Group columns:

  Statistics for each group level (when `by` is specified). Column names
  match group levels.

- *p*-value:

  Formatted *p*-values from statistical tests (when `test = TRUE` and
  `by` is specified)

The first row always shows sample sizes for each column. All numeric
output (counts, statistics, *p*-values) respects the `number_format`
setting for locale-appropriate formatting.

The returned object includes the following attributes accessible via
[`attr()`](https://rdrr.io/r/base/attr.html):

- raw_data:

  A data.table containing unformatted numeric values suitable for
  further statistical analysis or custom formatting. Includes additional
  columns for standard deviations, quartiles, etc.

- by_variable:

  The grouping variable name used (value of `by`)

- variables:

  The variables analyzed (value of `variables`)

## Details

**Variable Type Detection:**

The function automatically detects variable types and applies
appropriate summaries:

- **Continuous**: Numeric variables (integer or double) receive
  statistics specified in `stats_continuous`

- **Categorical**: Character, factor, or logical variables receive
  frequency counts and percentages

- **Time-to-Event**: Variables specified as `Surv(time, event)` display
  median survival with confidence intervals (level controlled by
  `conf_level`)

**Statistical Testing:**

When `test = TRUE` and `by` is specified:

- **Continuous with "auto"**: Parametric tests (*t*-test, ANOVA) for
  mean-based statistics; non-parametric tests (Wilcoxon, Kruskal-Wallis)
  for median-based statistics

- **Categorical with "auto"**: Fisher exact test when any expected cell
  frequency \< 5; χ² test otherwise

- **Survival**: Log-rank test for comparing survival curves

- **Range statistics**: No *p*-value computed (ranges are descriptive)

**Missing Data Handling:**

Missing values are handled differently by variable type:

- **Continuous**: NAs excluded from calculations; optionally shown as
  count when `na_include = TRUE`

- **Categorical**: NAs can be included as a category when
  `na_include = TRUE`. The `na_percent` parameter controls whether
  percentages are calculated with or without NAs in the denominator

- **Survival**: NAs in time or event excluded from analysis

**Formatting Conventions:**

All numeric output respects the `number_format` parameter. Separators
within ranges and confidence intervals adapt automatically to avoid
ambiguity:

- Mean ± SD: `"45.2 ± 12.3"` (US) or `"45,2 ± 12,3"` (EU)

- Median \[IQR\]: `"38.0 [28.0-52.0]"` (US) or `"38,0 [28,0–52,0]"` (EU,
  en-dash separator)

- Range: `"18.0-75.0"` (positive, US), `"-5.0 to 10.0"` (when bounds are
  negative)

- Survival: `"24.5 (21.2-28.9)"` (US) or `"24,5 (21,2–28,9)"` (EU)

- Counts ≥ 1000: `"1,234"` (US) or `"1.234"` (EU)

- *p*-values: `"< 0.001"` (US) or `"< 0,001"` (EU)

## See also

[`survtable`](https://phmcc.github.io/summata/reference/survtable.md)
for detailed survival summary tables,
[`fit`](https://phmcc.github.io/summata/reference/fit.md) for regression
modeling,
[`table2pdf`](https://phmcc.github.io/summata/reference/table2pdf.md)
for PDF export,
[`table2docx`](https://phmcc.github.io/summata/reference/table2docx.md)
for Word export,
[`table2html`](https://phmcc.github.io/summata/reference/table2html.md)
for HTML export

Other descriptive functions:
[`print.survtable()`](https://phmcc.github.io/summata/reference/print.survtable.md),
[`survtable()`](https://phmcc.github.io/summata/reference/survtable.md)

## Examples

``` r
# Load example clinical trial data
data(clintrial)

# Example 1: Basic descriptive table without grouping
desctable(clintrial,
        variables = c("age", "sex", "bmi"))
#>    Variable        Group            Total
#>      <char>       <char>           <char>
#> 1:        N                           850
#> 2:      age Median [IQR] 60.0 [52.0-67.0]
#> 3:      sex       Female      450 (52.9%)
#> 4:                  Male      400 (47.1%)
#> 5:      bmi Median [IQR] 28.0 [24.7-31.5]


# \donttest{
# Example 2: Grouped comparison with default tests
desctable(clintrial,
        by = "treatment",
        variables = c("age", "sex", "race", "bmi"))
#>    Variable        Group            Total          Control           Drug A           Drug B p-value
#>      <char>       <char>           <char>           <char>           <char>           <char>  <char>
#> 1:        N                           850              196              292              362        
#> 2:      age Median [IQR] 60.0 [52.0-67.0] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0]   0.465
#> 3:      sex       Female      450 (52.9%)      100 (51.0%)      164 (56.2%)      186 (51.4%)   0.394
#> 4:                  Male      400 (47.1%)       96 (49.0%)      128 (43.8%)      176 (48.6%)        
#> 5:     race        White      598 (70.4%)      147 (75.0%)      198 (67.8%)      253 (69.9%)   0.794
#> 6:                 Black      126 (14.8%)       25 (12.8%)       47 (16.1%)       54 (14.9%)        
#> 7:                 Asian       93 (10.9%)        17 (8.7%)       35 (12.0%)       41 (11.3%)        
#> 8:                 Other        33 (3.9%)         7 (3.6%)        12 (4.1%)        14 (3.9%)        
#> 9:      bmi Median [IQR] 28.0 [24.7-31.5] 28.0 [24.4-31.9] 28.2 [24.9-31.3] 27.8 [24.7-31.5]   0.831

# Example 3: Customize continuous statistics
desctable(clintrial,
        by = "treatment",
        variables = c("age", "bmi", "creatinine"),
        stats_continuous = c("median_iqr", "range"))
#>      Variable        Group            Total          Control           Drug A           Drug B p-value
#>        <char>       <char>           <char>           <char>           <char>           <char>  <char>
#> 1:          N                           850              196              292              362        
#> 2:        age Median [IQR] 60.0 [52.0-67.0] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0]   0.465
#> 3:                   Range        18.0-90.0        26.0-90.0        18.0-90.0        24.0-90.0        
#> 4:        bmi Median [IQR] 28.0 [24.7-31.5] 28.0 [24.4-31.9] 28.2 [24.9-31.3] 27.8 [24.7-31.5]   0.831
#> 5:                   Range        15.0-42.6        17.4-40.8        15.5-40.9        15.0-42.6        
#> 6: creatinine Median [IQR]    1.0 [0.8-1.2]    1.0 [0.8-1.2]    1.0 [0.8-1.2]    1.0 [0.8-1.2]   0.483
#> 7:                   Range          0.3-2.2          0.5-1.8          0.3-2.2          0.4-2.1        

# Example 4: Change categorical display format
desctable(clintrial,
        by = "treatment",
        variables = c("sex", "race", "smoking"),
        stats_categorical = "n")  # Show counts only
#>     Variable   Group  Total Control Drug A Drug B p-value
#>       <char>  <char> <char>  <char> <char> <char>  <char>
#>  1:        N            850     196    292    362        
#>  2:      sex  Female    450     100    164    186   0.394
#>  3:             Male    400      96    128    176        
#>  4:     race   White    598     147    198    253   0.794
#>  5:            Black    126      25     47     54        
#>  6:            Asian     93      17     35     41        
#>  7:            Other     33       7     12     14        
#>  8:  smoking   Never    337      83    123    131   0.450
#>  9:           Former    311      66    101    144        
#> 10:          Current    185      42     64     79        

# Example 5: Include missing values
desctable(clintrial,
        by = "treatment",
        variables = c("age", "smoking", "hypertension"),
        na_include = TRUE,
        na_label = "Missing")
#>        Variable        Group            Total          Control           Drug A           Drug B p-value
#>          <char>       <char>           <char>           <char>           <char>           <char>  <char>
#> 1:            N                           850              196              292              362        
#> 2:          age Median [IQR] 60.0 [52.0-67.0] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0]   0.465
#> 3:      smoking        Never      337 (40.5%)       83 (43.5%)      123 (42.7%)      131 (37.0%)   0.450
#> 4:                    Former      311 (37.3%)       66 (34.6%)      101 (35.1%)      144 (40.7%)        
#> 5:                   Current      185 (22.2%)       42 (22.0%)       64 (22.2%)       79 (22.3%)        
#> 6:                   Missing               17                5                4                8        
#> 7: hypertension           No      504 (60.4%)      108 (56.2%)      179 (62.2%)      217 (61.1%)   0.401
#> 8:                       Yes      331 (39.6%)       84 (43.8%)      109 (37.8%)      138 (38.9%)        
#> 9:                   Missing               15                4                4                7        

# Example 6: Disable statistical testing
desctable(clintrial,
        by = "treatment",
        variables = c("age", "sex", "bmi"),
        test = FALSE)
#>    Variable        Group            Total          Control           Drug A           Drug B
#>      <char>       <char>           <char>           <char>           <char>           <char>
#> 1:        N                           850              196              292              362
#> 2:      age Median [IQR] 60.0 [52.0-67.0] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0]
#> 3:      sex       Female      450 (52.9%)      100 (51.0%)      164 (56.2%)      186 (51.4%)
#> 4:                  Male      400 (47.1%)       96 (49.0%)      128 (43.8%)      176 (48.6%)
#> 5:      bmi Median [IQR] 28.0 [24.7-31.5] 28.0 [24.4-31.9] 28.2 [24.9-31.3] 27.8 [24.7-31.5]

# Example 7: Force specific tests
desctable(clintrial,
        by = "surgery",
        variables = c("age", "sex"),
        test_continuous = "t",      # t-test instead of auto
        test_categorical = "fisher") # Fisher test instead of auto
#>    Variable        Group            Total               No              Yes p-value
#>      <char>       <char>           <char>           <char>           <char>  <char>
#> 1:        N                           850              480              370        
#> 2:      age Median [IQR] 60.0 [52.0-67.0] 62.0 [54.0-70.0] 57.0 [49.0-66.0] < 0.001
#> 3:      sex       Female      450 (52.9%)      246 (51.2%)      204 (55.1%)   0.268
#> 4:                  Male      400 (47.1%)      234 (48.8%)      166 (44.9%)        

# Example 8: Adjust decimal places
desctable(clintrial,
        by = "treatment",
        variables = c("age", "bmi"),
        digits = 2,    # 2 decimals for continuous
        p_digits = 4)  # 4 decimals for p-values
#>    Variable        Group               Total             Control              Drug A              Drug B p-value
#>      <char>       <char>              <char>              <char>              <char>              <char>  <char>
#> 1:        N                              850                 196                 292                 362        
#> 2:      age Median [IQR] 60.00 [52.00-67.00] 59.00 [52.00-66.50] 60.00 [51.00-68.00] 61.00 [53.00-68.00]  0.4654
#> 3:      bmi Median [IQR] 28.00 [24.70-31.50] 28.00 [24.40-31.90] 28.25 [24.90-31.30] 27.80 [24.70-31.50]  0.8314

# Example 9: Custom variable labels
labels <- c(
    age = "Age (years)",
    sex = "Sex",
    bmi = "Body Mass Index (kg/m\u00b2)",
    treatment = "Treatment Arm"
)

desctable(clintrial,
        by = "treatment",
        variables = c("age", "sex", "bmi"),
        labels = labels)
#>                   Variable        Group            Total          Control           Drug A           Drug B p-value
#>                     <char>       <char>           <char>           <char>           <char>           <char>  <char>
#> 1:                       N                           850              196              292              362        
#> 2:             Age (years) Median [IQR] 60.0 [52.0-67.0] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0]   0.465
#> 3:                     Sex       Female      450 (52.9%)      100 (51.0%)      164 (56.2%)      186 (51.4%)   0.394
#> 4:                                 Male      400 (47.1%)       96 (49.0%)      128 (43.8%)      176 (48.6%)        
#> 5: Body Mass Index (kg/m²) Median [IQR] 28.0 [24.7-31.5] 28.0 [24.4-31.9] 28.2 [24.9-31.3] 27.8 [24.7-31.5]   0.831

# Example 10: Position total column last
desctable(clintrial,
        by = "treatment",
        variables = c("age", "sex"),
        total = "last")
#>    Variable        Group          Control           Drug A           Drug B            Total p-value
#>      <char>       <char>           <char>           <char>           <char>           <char>  <char>
#> 1:        N                           196              292              362              850        
#> 2:      age Median [IQR] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0] 60.0 [52.0-67.0]   0.465
#> 3:      sex       Female      100 (51.0%)      164 (56.2%)      186 (51.4%)      450 (52.9%)   0.394
#> 4:                  Male       96 (49.0%)      128 (43.8%)      176 (48.6%)      400 (47.1%)        

# Example 11: Exclude total column
desctable(clintrial,
        by = "treatment",
        variables = c("age", "sex"),
        total = FALSE)
#>    Variable        Group          Control           Drug A           Drug B p-value
#>      <char>       <char>           <char>           <char>           <char>  <char>
#> 1:        N                           196              292              362        
#> 2:      age Median [IQR] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0]   0.465
#> 3:      sex       Female      100 (51.0%)      164 (56.2%)      186 (51.4%)   0.394
#> 4:                  Male       96 (49.0%)      128 (43.8%)      176 (48.6%)        

# Example 12: Survival analysis
desctable(clintrial,
        by = "treatment",
        variables = "Surv(os_months, os_status)")
#>                      Variable           Group            Total          Control           Drug A           Drug B p-value
#>                        <char>          <char>           <char>           <char>           <char>           <char>  <char>
#> 1:                          N                              850              196              292              362        
#> 2: Surv(os_months, os_status) Median (95% CI) 19.4 (16.2-23.4) 14.7 (10.5-19.2) 33.6 (24.5-42.2) 14.7 (11.0-21.8) < 0.001

# Example 13: Multiple survival endpoints
desctable(clintrial,
        by = "treatment",
        variables = c(
            "Surv(pfs_months, pfs_status)",
            "Surv(os_months, os_status)"
        ),
        labels = c(
            "Surv(pfs_months, pfs_status)" = "Progression-Free Survival",
            "Surv(os_months, os_status)" = "Overall Survival"
        ))
#>                     Variable           Group            Total          Control           Drug A           Drug B p-value
#>                       <char>          <char>           <char>           <char>           <char>           <char>  <char>
#> 1:                         N                              850              196              292              362        
#> 2: Progression-Free Survival Median (95% CI)    6.0 (5.0-7.0)    4.7 (3.6-6.9)   9.0 (7.5-12.4)    4.7 (3.7-6.2) < 0.001
#> 3:          Overall Survival Median (95% CI) 19.4 (16.2-23.4) 14.7 (10.5-19.2) 33.6 (24.5-42.2) 14.7 (11.0-21.8) < 0.001

# Example 14: Mixed variable types
desctable(clintrial,
        by = "treatment",
        variables = c(
            "age", "sex", "race",           # Demographics
            "bmi", "creatinine",            # Labs
            "smoking", "hypertension",      # Risk factors
            "Surv(os_months, os_status)"    # Survival
        ))
#>                       Variable           Group            Total          Control           Drug A           Drug B p-value
#>                         <char>          <char>           <char>           <char>           <char>           <char>  <char>
#>  1:                          N                              850              196              292              362        
#>  2:                        age    Median [IQR] 60.0 [52.0-67.0] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0]   0.465
#>  3:                        sex          Female      450 (52.9%)      100 (51.0%)      164 (56.2%)      186 (51.4%)   0.394
#>  4:                                       Male      400 (47.1%)       96 (49.0%)      128 (43.8%)      176 (48.6%)        
#>  5:                       race           White      598 (70.4%)      147 (75.0%)      198 (67.8%)      253 (69.9%)   0.794
#>  6:                                      Black      126 (14.8%)       25 (12.8%)       47 (16.1%)       54 (14.9%)        
#>  7:                                      Asian       93 (10.9%)        17 (8.7%)       35 (12.0%)       41 (11.3%)        
#>  8:                                      Other        33 (3.9%)         7 (3.6%)        12 (4.1%)        14 (3.9%)        
#>  9:                        bmi    Median [IQR] 28.0 [24.7-31.5] 28.0 [24.4-31.9] 28.2 [24.9-31.3] 27.8 [24.7-31.5]   0.831
#> 10:                 creatinine    Median [IQR]    1.0 [0.8-1.2]    1.0 [0.8-1.2]    1.0 [0.8-1.2]    1.0 [0.8-1.2]   0.483
#> 11:                    smoking           Never      337 (40.5%)       83 (43.5%)      123 (42.7%)      131 (37.0%)   0.450
#> 12:                                     Former      311 (37.3%)       66 (34.6%)      101 (35.1%)      144 (40.7%)        
#> 13:                                    Current      185 (22.2%)       42 (22.0%)       64 (22.2%)       79 (22.3%)        
#> 14:               hypertension              No      504 (60.4%)      108 (56.2%)      179 (62.2%)      217 (61.1%)   0.401
#> 15:                                        Yes      331 (39.6%)       84 (43.8%)      109 (37.8%)      138 (38.9%)        
#> 16: Surv(os_months, os_status) Median (95% CI) 19.4 (16.2-23.4) 14.7 (10.5-19.2) 33.6 (24.5-42.2) 14.7 (11.0-21.8) < 0.001

# Example 15: Export table
table1 <- desctable(clintrial,
                  by = "treatment",
                  variables = c("age", "sex", "bmi"))

# Can export directly to PDF/LaTeX/HTML for publication
# table2pdf(table1, "table1.pdf")
# table2docx(table1, "table1.docx")

# Example 16: Three or more groups
desctable(clintrial,
        by = "stage",  # Assuming stage has 3+ levels
        variables = c("age", "sex", "bmi"))
#>    Variable        Group            Total                I               II              III               IV p-value
#>      <char>       <char>           <char>           <char>           <char>           <char>           <char>  <char>
#> 1:        N                           850              211              263              241              132        
#> 2:      age Median [IQR] 60.0 [52.0-67.0] 61.0 [53.0-67.0] 61.0 [51.0-70.0] 58.0 [51.0-67.0] 60.0 [53.0-66.0]   0.481
#> 3:      sex       Female      449 (53.0%)      110 (52.1%)      131 (49.8%)      134 (55.6%)       74 (56.1%)   0.515
#> 4:                  Male      398 (47.0%)      101 (47.9%)      132 (50.2%)      107 (44.4%)       58 (43.9%)        
#> 5:      bmi Median [IQR] 28.0 [24.7-31.5] 27.9 [24.8-30.9] 28.1 [24.2-31.5] 27.8 [24.7-31.5] 28.5 [25.1-32.1]   0.595
# Automatically uses ANOVA/Kruskal-Wallis and chi-squared

# Example 17: Access raw unformatted data
result <- desctable(clintrial,
                  by = "treatment",
                  variables = c("age", "bmi"))
raw_data <- attr(result, "raw_data")
print(raw_data)
#>    Variable      Group  stat_type Total Total_q1 Total_q3 Total_n Control Control_q1 Control_q3 Control_n Drug A Drug A_q1 Drug A_q3 Drug A_n Drug B Drug B_q1 Drug B_q3 Drug B_n
#>      <char>     <char>     <char> <num>    <num>    <num>   <int>   <num>      <num>      <num>     <int>  <num>     <num>     <num>    <int>  <num>     <num>     <num>    <int>
#> 1:      age median_iqr median_iqr    60     52.0     67.0     850      59       52.0       66.5       196  60.00      51.0      68.0      292   61.0      53.0      68.0      362
#> 2:      bmi median_iqr median_iqr    28     24.7     31.5     838      28       24.4       31.9       194  28.25      24.9      31.3      288   27.8      24.7      31.5      356
#>      p_value
#>        <num>
#> 1: 0.4653885
#> 2: 0.8314435
# Raw data includes unformatted numbers, SDs, quartiles, etc.

# Example 18: Check which grouping variable was used
result <- desctable(clintrial,
                  by = "treatment",
                  variables = c("age", "sex"))
attr(result, "by_variable")  # "treatment"
#> [1] "treatment"

# Example 19: NA percentage calculation options
# Include NAs in percentage denominator (all sum to 100%)
desctable(clintrial,
        by = "treatment",
        variables = "smoking",
        na_include = TRUE,
        na_percent = TRUE)
#>    Variable   Group       Total    Control      Drug A      Drug B p-value
#>      <char>  <char>      <char>     <char>      <char>      <char>  <char>
#> 1:        N                 850        196         292         362        
#> 2:  smoking   Never 337 (39.6%) 83 (42.3%) 123 (42.1%) 131 (36.2%)   0.450
#> 3:           Former 311 (36.6%) 66 (33.7%) 101 (34.6%) 144 (39.8%)        
#> 4:          Current 185 (21.8%) 42 (21.4%)  64 (21.9%)  79 (21.8%)        
#> 5:          Unknown   17 (2.0%)   5 (2.6%)    4 (1.4%)    8 (2.2%)        

# Exclude NAs from denominator (non-missing sum to 100%)
desctable(clintrial,
        by = "treatment",
        variables = "smoking",
        na_include = TRUE,
        na_percent = FALSE)
#>    Variable   Group       Total    Control      Drug A      Drug B p-value
#>      <char>  <char>      <char>     <char>      <char>      <char>  <char>
#> 1:        N                 850        196         292         362        
#> 2:  smoking   Never 337 (40.5%) 83 (43.5%) 123 (42.7%) 131 (37.0%)   0.450
#> 3:           Former 311 (37.3%) 66 (34.6%) 101 (35.1%) 144 (40.7%)        
#> 4:          Current 185 (22.2%) 42 (22.0%)  64 (22.2%)  79 (22.3%)        
#> 5:          Unknown          17          5           4           8        

# Example 20: Passing additional test arguments
# Equal variance t-test
desctable(clintrial,
        by = "sex",
        variables = "age",
        test_continuous = "t",
        var.equal = TRUE)
#>    Variable        Group            Total           Female             Male p-value
#>      <char>       <char>           <char>           <char>           <char>  <char>
#> 1:        N                           850              450              400        
#> 2:      age Median [IQR] 60.0 [52.0-67.0] 60.0 [52.0-67.0] 60.0 [51.5-69.0]   0.951

# Example 21: European number formatting
desctable(clintrial,
        by = "treatment",
        variables = c("age", "sex", "bmi"),
        number_format = "eu")
#>    Variable        Group            Total          Control           Drug A           Drug B p-value
#>      <char>       <char>           <char>           <char>           <char>           <char>  <char>
#> 1:        N                           850              196              292              362        
#> 2:      age Median [IQR] 60,0 [52,0–67,0] 59,0 [52,0–66,5] 60,0 [51,0–68,0] 61,0 [53,0–68,0]   0,465
#> 3:      sex       Female      450 (52,9%)      100 (51,0%)      164 (56,2%)      186 (51,4%)   0,394
#> 4:                  Male      400 (47,1%)       96 (49,0%)      128 (43,8%)      176 (48,6%)        
#> 5:      bmi Median [IQR] 28,0 [24,7–31,5] 28,0 [24,4–31,9] 28,2 [24,9–31,3] 27,8 [24,7–31,5]   0,831

# Example 22: Complete Table 1 for publication
table1 <- desctable(
    data = clintrial,
    by = "treatment",
    variables = c(
        "age", "sex", "race", "ethnicity", "bmi",
        "smoking", "hypertension", "diabetes",
        "ecog", "creatinine", "hemoglobin",
        "site", "stage", "grade",
        "Surv(os_months, os_status)"
    ),
    labels = clintrial_labels,
    stats_continuous = c("median_iqr", "range"),
    total = TRUE,
    na_include = FALSE
)
print(table1)
#>                        Variable                     Group            Total          Control           Drug A           Drug B p-value
#>                          <char>                    <char>           <char>           <char>           <char>           <char>  <char>
#>  1:                           N                                        850              196              292              362        
#>  2:                 Age (years)              Median [IQR] 60.0 [52.0-67.0] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0]   0.465
#>  3:                                                 Range        18.0-90.0        26.0-90.0        18.0-90.0        24.0-90.0        
#>  4:                         Sex                    Female      450 (52.9%)      100 (51.0%)      164 (56.2%)      186 (51.4%)   0.394
#>  5:                                                  Male      400 (47.1%)       96 (49.0%)      128 (43.8%)      176 (48.6%)        
#>  6:                        Race                     White      598 (70.4%)      147 (75.0%)      198 (67.8%)      253 (69.9%)   0.794
#>  7:                                                 Black      126 (14.8%)       25 (12.8%)       47 (16.1%)       54 (14.9%)        
#>  8:                                                 Asian       93 (10.9%)        17 (8.7%)       35 (12.0%)       41 (11.3%)        
#>  9:                                                 Other        33 (3.9%)         7 (3.6%)        12 (4.1%)        14 (3.9%)        
#> 10:                   Ethnicity              Non-Hispanic      744 (87.5%)      172 (87.8%)      253 (86.6%)      319 (88.1%)   0.846
#> 11:                                              Hispanic      106 (12.5%)       24 (12.2%)       39 (13.4%)       43 (11.9%)        
#> 12:     Body Mass Index (kg/m²)              Median [IQR] 28.0 [24.7-31.5] 28.0 [24.4-31.9] 28.2 [24.9-31.3] 27.8 [24.7-31.5]   0.831
#> 13:                                                 Range        15.0-42.6        17.4-40.8        15.5-40.9        15.0-42.6        
#> 14:              Smoking Status                     Never      337 (40.5%)       83 (43.5%)      123 (42.7%)      131 (37.0%)   0.450
#> 15:                                                Former      311 (37.3%)       66 (34.6%)      101 (35.1%)      144 (40.7%)        
#> 16:                                               Current      185 (22.2%)       42 (22.0%)       64 (22.2%)       79 (22.3%)        
#> 17:                Hypertension                        No      504 (60.4%)      108 (56.2%)      179 (62.2%)      217 (61.1%)   0.401
#> 18:                                                   Yes      331 (39.6%)       84 (43.8%)      109 (37.8%)      138 (38.9%)        
#> 19:                    Diabetes                        No      637 (76.4%)      152 (79.6%)      218 (75.7%)      267 (75.2%)   0.490
#> 20:                                                   Yes      197 (23.6%)       39 (20.4%)       70 (24.3%)       88 (24.8%)        
#> 21:     ECOG Performance Status                         0      265 (31.5%)       56 (28.9%)       83 (28.5%)      126 (35.3%)   0.150
#> 22:                                                     1      302 (35.9%)       79 (40.7%)      115 (39.5%)      108 (30.3%)        
#> 23:                                                     2      238 (28.3%)       51 (26.3%)       79 (27.1%)      108 (30.3%)        
#> 24:                                                     3        37 (4.4%)         8 (4.1%)        14 (4.8%)        15 (4.2%)        
#> 25: Baseline Creatinine (mg/dL)              Median [IQR]    1.0 [0.8-1.2]    1.0 [0.8-1.2]    1.0 [0.8-1.2]    1.0 [0.8-1.2]   0.483
#> 26:                                                 Range          0.3-2.2          0.5-1.8          0.3-2.2          0.4-2.1        
#> 27:  Baseline Hemoglobin (g/dL)              Median [IQR] 13.0 [11.7-14.2] 13.0 [11.4-14.0] 12.8 [11.7-14.1] 13.1 [11.8-14.5]   0.129
#> 28:                                                 Range         7.5-18.0         8.0-17.6         7.5-18.0         8.1-18.0        
#> 29:                  Study Site                Site Alpha        76 (8.9%)        18 (9.2%)        23 (7.9%)        35 (9.7%)   0.252
#> 30:                                             Site Beta        80 (9.4%)        13 (6.6%)        27 (9.2%)       40 (11.0%)        
#> 31:                                            Site Gamma        72 (8.5%)       20 (10.2%)        15 (5.1%)       37 (10.2%)        
#> 32:                                            Site Delta       89 (10.5%)       22 (11.2%)        27 (9.2%)       40 (11.0%)        
#> 33:                                          Site Epsilon       86 (10.1%)       20 (10.2%)       33 (11.3%)        33 (9.1%)        
#> 34:                                             Site Zeta       93 (10.9%)        16 (8.2%)       42 (14.4%)        35 (9.7%)        
#> 35:                                              Site Eta       89 (10.5%)        19 (9.7%)        29 (9.9%)       41 (11.3%)        
#> 36:                                            Site Theta        75 (8.8%)        18 (9.2%)       32 (11.0%)        25 (6.9%)        
#> 37:                                             Site Iota      102 (12.0%)       30 (15.3%)       34 (11.6%)       38 (10.5%)        
#> 38:                                            Site Kappa       88 (10.4%)       20 (10.2%)       30 (10.3%)       38 (10.5%)        
#> 39:               Disease Stage                         I      211 (24.9%)       59 (30.4%)       78 (26.7%)       74 (20.5%)   0.007
#> 40:                                                    II      263 (31.1%)       65 (33.5%)       93 (31.8%)      105 (29.1%)        
#> 41:                                                   III      241 (28.5%)       39 (20.1%)       75 (25.7%)      127 (35.2%)        
#> 42:                                                    IV      132 (15.6%)       31 (16.0%)       46 (15.8%)       55 (15.2%)        
#> 43:                 Tumor Grade       Well-differentiated      153 (18.2%)       42 (21.8%)       54 (18.7%)       57 (15.9%)   0.331
#> 44:                             Moderately differentiated      412 (49.0%)       96 (49.7%)      143 (49.5%)      173 (48.3%)        
#> 45:                                 Poorly differentiated      275 (32.7%)       55 (28.5%)       92 (31.8%)      128 (35.8%)        
#> 46:   Overall Survival (months)           Median (95% CI) 19.4 (16.2-23.4) 14.7 (10.5-19.2) 33.6 (24.5-42.2) 14.7 (11.0-21.8) < 0.001
#>                        Variable                     Group            Total          Control           Drug A           Drug B p-value
#>                          <char>                    <char>           <char>           <char>           <char>           <char>  <char>
# }
```
