# Simulated Clinical Trial Dataset

A simulated dataset from a hypothetical multi-center oncology clinical
trial comparing two experimental drugs against control. Designed to
demonstrate the full capabilities of descriptive and regression analysis
functions.

## Usage

``` r
clintrial
```

## Format

A data frame with 850 observations and 32 variables:

- patient_id:

  Unique patient identifier (character)

- age:

  Age at enrollment in years (numeric: 18-90)

- sex:

  Biological sex (factor: Female, Male)

- race:

  Self-reported race (factor: White, Black, Asian, Other)

- ethnicity:

  Hispanic ethnicity (factor: Non-Hispanic, Hispanic)

- bmi:

  Body mass index in kg/m\\^2\\ (numeric)

- smoking:

  Smoking history (factor: Never, Former, Current)

- hypertension:

  Hypertension diagnosis (factor: No, Yes)

- diabetes:

  Diabetes diagnosis (factor: No, Yes)

- ecog:

  ECOG performance status (factor: 0, 1, 2, 3)

- creatinine:

  Baseline creatinine in mg/dL (numeric)

- hemoglobin:

  Baseline hemoglobin in g/dL (numeric)

- biomarker_x:

  Serum biomarker A in ng/mL (numeric)

- biomarker_y:

  Serum biomarker B in U/L (numeric)

- site:

  Enrolling site (factor: Site Alpha through Site Kappa)

- grade:

  Tumor grade (factor: Well/Moderately/Poorly differentiated)

- stage:

  Disease stage at diagnosis (factor: I, II, III, IV)

- treatment:

  Randomized treatment (factor: Control, Drug A, Drug B)

- surgery:

  Surgical resection (factor: No, Yes)

- any_complication:

  Any post-operative complication (factor: No, Yes)

- wound_infection:

  Post-operative wound infection (factor: No, Yes)

- icu_admission:

  ICU admission required (factor: No, Yes)

- readmission_30d:

  Hospital readmission within 30 days (factor: No, Yes)

- pain_score:

  Pain score at discharge (numeric: 0-10)

- recovery_days:

  Days to functional recovery (numeric)

- los_days:

  Hospital length of stay in days (numeric)

- ae_count:

  Adverse event count (integer). Overdispersed count suitable for
  negative binomial or quasipoisson regression.

- fu_count:

  Follow-up visit count (integer). Equidispersed count suitable for
  standard Poisson regression.

- pfs_months:

  Progression-Free Survival Time (months)

- pfs_status:

  Progression or Death Event

- os_months:

  Overall survival time in months (numeric)

- os_status:

  Death indicator (numeric: 0=censored, 1=death)

## Source

Simulated data for demonstration purposes

## Details

This dataset includes realistic correlations between variables: -
Survival is worse with higher stage, ECOG, age, and biomarker_x -
Treatment effects show Drug B \> Drug A \> Control - `ae_count` is
overdispersed (variance \> mean) for negative binomial demos -
`fu_count` is equidispersed (variance \\\approx\\ mean) for Poisson
demos - Approximately 2% of values are missing at random - Median
follow-up is approximately 30 months

## See also

Other sample data:
[`clintrial_labels`](https://phmcc.github.io/summata/reference/clintrial_labels.md)

## Examples

``` r
data(clintrial)
data(clintrial_labels)

# Descriptive statistics by treatment arm
desctable(clintrial,
        by = "treatment", 
        variables = c("age", "sex", "stage", "ecog", 
                     "biomarker_x", "Surv(os_months, os_status)"),
        labels = clintrial_labels)
#>                      Variable           Group            Total          Control           Drug A           Drug B p-value
#>                        <char>          <char>           <char>           <char>           <char>           <char>  <char>
#>  1:                         N                              850              196              292              362        
#>  2:               Age (years)    Median [IQR] 60.0 [52.0-67.0] 59.0 [52.0-66.5] 60.0 [51.0-68.0] 61.0 [53.0-68.0]   0.465
#>  3:                       Sex          Female      450 (52.9%)      100 (51.0%)      164 (56.2%)      186 (51.4%)   0.394
#>  4:                                      Male      400 (47.1%)       96 (49.0%)      128 (43.8%)      176 (48.6%)        
#>  5:             Disease Stage               I      211 (24.9%)       59 (30.4%)       78 (26.7%)       74 (20.5%)   0.007
#>  6:                                        II      263 (31.1%)       65 (33.5%)       93 (31.8%)      105 (29.1%)        
#>  7:                                       III      241 (28.5%)       39 (20.1%)       75 (25.7%)      127 (35.2%)        
#>  8:                                        IV      132 (15.6%)       31 (16.0%)       46 (15.8%)       55 (15.2%)        
#>  9:   ECOG Performance Status               0      265 (31.5%)       56 (28.9%)       83 (28.5%)      126 (35.3%)   0.150
#> 10:                                         1      302 (35.9%)       79 (40.7%)      115 (39.5%)      108 (30.3%)        
#> 11:                                         2      238 (28.3%)       51 (26.3%)       79 (27.1%)      108 (30.3%)        
#> 12:                                         3        37 (4.4%)         8 (4.1%)        14 (4.8%)        15 (4.2%)        
#> 13:       Biomarker X (ng/mL)    Median [IQR]    5.4 [3.9-7.5]    5.3 [4.1-7.2]    5.2 [3.8-7.3]    5.5 [3.8-7.7]   0.318
#> 14: Overall Survival (months) Median (95% CI) 19.4 (16.2-23.4) 14.7 (10.5-19.2) 33.6 (24.5-42.2) 14.7 (11.0-21.8) < 0.001

# Poisson regression for equidispersed counts
fit(clintrial,
    outcome = "fu_count",
    predictors = c("age", "stage", "treatment"),
    model_type = "glm",
    family = "poisson",
    labels = clintrial_labels)
#> 
#> Multivariable Poisson Model
#> Formula: fu_count ~ age + stage + treatment
#> n = 839, Events = 5517
#> 
#>           Variable   Group      n Events     aRR (95% CI) p-value
#>             <char>  <char> <char> <char>           <char>  <char>
#> 1:     Age (years)       -    839  5,517 1.00 (0.99-1.00)   0.005
#> 2:   Disease Stage       I    207  1,238        reference       -
#> 3:                      II    261  1,638 1.05 (0.97-1.13)   0.234
#> 4:                     III    240  1,638 1.12 (1.04-1.21)   0.002
#> 5:                      IV    131  1,003 1.27 (1.17-1.38) < 0.001
#> 6: Treatment Group Control    191  1,169        reference       -
#> 7:                  Drug A    290  1,910 1.08 (1.00-1.16)   0.052
#> 8:                  Drug B    358  2,438 1.11 (1.03-1.19)   0.005

# Negative binomial for overdispersed counts
fit(clintrial,
    outcome = "ae_count",
    predictors = c("age", "treatment", "diabetes"),
    model_type = "negbin",
    labels = clintrial_labels)
#> 
#> Multivariable Negative Binomial Model
#> Formula: ae_count ~ age + treatment + diabetes
#> n = 824, Events = 4506
#> 
#>           Variable   Group      n Events     aRR (95% CI) p-value
#>             <char>  <char> <char> <char>           <char>  <char>
#> 1:     Age (years)       -    824  4,506 1.01 (1.01-1.01) < 0.001
#> 2: Treatment Group Control    189    826        reference       -
#> 3:                  Drug A    287  1,221 0.96 (0.82-1.11)   0.585
#> 4:                  Drug B    348  2,459 1.52 (1.32-1.75) < 0.001
#> 5:        Diabetes      No    630  2,998        reference       -
#> 6:                     Yes    194  1,508 1.56 (1.38-1.77) < 0.001

# Complete analysis pipeline
fullfit(clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "sex", "stage", "grade", "ecog",
                      "smoking", "biomarker_x", "biomarker_y", "treatment"),
        method = "screen",
        p_threshold = 0.20,
        model_type = "coxph",
        labels = clintrial_labels)
#> Running univariable analysis...
#> Fitting multivariable model with 8 predictors...
#> 
#> Fullfit Analysis Results
#> Outcome: Surv(os_months, os_status)
#> Model Type: coxph
#> Method: screen (p < 0.2)
#> Predictors Screened: 9
#> Multivariable Predictors: 8
#> 
#>                    Variable                     Group      n Events      HR (95% CI)   Uni p     aHR (95% CI) Multi p
#>                      <char>                    <char> <char> <char>           <char>  <char>           <char>  <char>
#>  1:             Age (years)                         -    850    609 1.03 (1.03-1.04) < 0.001 1.03 (1.03-1.04) < 0.001
#>  2:                     Sex                    Female    450    298        reference       -        reference       -
#>  3:                                              Male    400    311 1.30 (1.11-1.53)   0.001 1.33 (1.13-1.57) < 0.001
#>  4:           Disease Stage                         I    211    127        reference       -        reference       -
#>  5:                                                II    263    172 1.12 (0.89-1.41)   0.337 1.10 (0.87-1.39)   0.413
#>  6:                                               III    241    186 1.69 (1.35-2.11) < 0.001 1.81 (1.43-2.28) < 0.001
#>  7:                                                IV    132    121 3.18 (2.47-4.09) < 0.001 3.76 (2.89-4.89) < 0.001
#>  8:             Tumor Grade       Well-differentiated    153     95        reference       -        reference       -
#>  9:                         Moderately differentiated    412    297 1.36 (1.08-1.71)   0.010 1.52 (1.20-1.93) < 0.001
#> 10:                             Poorly differentiated    275    208 1.62 (1.27-2.06) < 0.001 1.92 (1.49-2.46) < 0.001
#> 11: ECOG Performance Status                         0    265    159        reference       -        reference       -
#> 12:                                                 1    302    212 1.36 (1.11-1.67)   0.003 1.52 (1.23-1.87) < 0.001
#> 13:                                                 2    238    194 1.86 (1.51-2.29) < 0.001 2.06 (1.66-2.55) < 0.001
#> 14:                                                 3     37     37 3.06 (2.13-4.38) < 0.001 3.43 (2.37-4.97) < 0.001
#> 15:          Smoking Status                     Never    337    248        reference       -        reference       -
#> 16:                                            Former    311    203 0.84 (0.70-1.02)   0.074 0.94 (0.78-1.14)   0.546
#> 17:                                           Current    185    143 1.19 (0.97-1.46)   0.103 1.29 (1.05-1.60)   0.017
#> 18:     Biomarker X (ng/mL)                         -    842    602 1.11 (1.08-1.13) < 0.001 1.08 (1.05-1.11) < 0.001
#> 19:       Biomarker Y (U/L)                         -    842    601 1.00 (1.00-1.00)   0.934                -       -
#> 20:         Treatment Group                   Control    196    151        reference       -        reference       -
#> 21:                                            Drug A    292    184 0.64 (0.52-0.80) < 0.001 0.53 (0.42-0.66) < 0.001
#> 22:                                            Drug B    362    274 0.94 (0.77-1.15)   0.567 0.77 (0.63-0.95)   0.014
#>                    Variable                     Group      n Events      HR (95% CI)   Uni p     aHR (95% CI) Multi p
#>                      <char>                    <char> <char> <char>           <char>  <char>           <char>  <char>
        
```
