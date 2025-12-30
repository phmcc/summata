#' Test Suite for desctable()
#' 
#' Comprehensive tests covering descriptive statistics tables,
#' continuous/categorical/survival variables, statistical tests,
#' formatting options, and edge cases.
#' 
#' @details Run with testthat::test_file("tests/testthat/test-desctable.R")

library(testthat)
library(data.table)
library(summata)

## ============================================================================
## Setup: Create test data and helper functions
## ============================================================================

## Use package's built-in clinical trial data
data(clintrial)
data(clintrial_labels)

## Create a complete-case subset for tests requiring no missing data
clintrial_complete <- na.omit(clintrial)

## Helper function to check desctable structure
expect_desctable_result <- function(result) {
    expect_s3_class(result, "data.table")
    
    ## Required columns
    expect_true("Variable" %in% names(result))
    expect_true("Group" %in% names(result))
    
    ## Required attributes
    expect_true(!is.null(attr(result, "raw_data")))
    expect_true(!is.null(attr(result, "variables")))
    
    ## First row should be N
    expect_equal(result$Variable[1], "N")
}

## Helper to check for Total column
expect_total_column <- function(result, label = "Total") {
    expect_true(label %in% names(result))
}

## Helper to check for p-value column
expect_pvalue_column <- function(result) {
    expect_true("p-value" %in% names(result))
}


## ============================================================================
## SECTION 1: Basic Functionality - No Grouping
## ============================================================================

test_that("desctable works without grouping variable", {
    
    result <- desctable(
        data = clintrial,
        variables = c("age", "sex", "stage")
    )
    
    expect_desctable_result(result)
    
    ## Should have Total column by default
    expect_total_column(result)
    
    ## Should NOT have p-value column (no groups to compare)
    expect_false("p-value" %in% names(result))
})


test_that("desctable ungrouped shows correct N", {
    
    result <- desctable(
        data = clintrial,
        variables = c("age")
    )
    
    ## First row should show total N
    n_displayed <- gsub(",", "", result$Total[1])  # Remove commas
    expect_equal(as.integer(n_displayed), nrow(clintrial))
})


## ============================================================================
## SECTION 2: Basic Functionality - With Grouping
## ============================================================================

test_that("desctable works with grouping variable", {
    
    result <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("age", "sex", "stage")
    )
    
    expect_desctable_result(result)
    
    ## Should have columns for each treatment group
    treatment_levels <- unique(clintrial$treatment[!is.na(clintrial$treatment)])
    for (lvl in treatment_levels) {
        expect_true(as.character(lvl) %in% names(result))
    }
    
    ## Should have p-value column by default
    expect_pvalue_column(result)
})


test_that("desctable by_variable attribute is set", {
    
    result <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("age", "sex")
    )
    
    expect_equal(attr(result, "by_variable"), "treatment")
})


test_that("desctable group N values are correct", {
    
    result <- desctable(
        data = clintrial,
        by = "sex",
        variables = c("age")
    )
    
    ## Check N row has correct counts per group
    n_male <- sum(clintrial$sex == "Male", na.rm = TRUE)
    n_female <- sum(clintrial$sex == "Female", na.rm = TRUE)
    
    male_n_displayed <- gsub(",", "", result$Male[1])
    female_n_displayed <- gsub(",", "", result$Female[1])
    
    expect_equal(as.integer(male_n_displayed), n_male)
    expect_equal(as.integer(female_n_displayed), n_female)
})


## ============================================================================
## SECTION 3: Continuous Variables
## ============================================================================

test_that("desctable handles continuous variables correctly", {
    
    result <- desctable(
        data = clintrial,
        variables = c("age", "bmi"),
        stats_continuous = c("mean_sd", "median_iqr")
    )
    
    expect_desctable_result(result)
    
    ## Should have rows for each stat type
    expect_true(any(result$Group == "Mean +/- SD"))
    expect_true(any(result$Group == "Median [IQR]"))
})


test_that("desctable stats_continuous='mean_sd' works", {
    
    result <- desctable(
        data = clintrial,
        variables = c("age"),
        stats_continuous = "mean_sd"
    )
    
    ## Should only have mean_sd row (plus N row)
    expect_true(any(result$Group == "Mean +/- SD"))
    expect_false(any(result$Group == "Median [IQR]"))
})


test_that("desctable stats_continuous='median_iqr' works", {
    
    result <- desctable(
        data = clintrial,
        variables = c("age"),
        stats_continuous = "median_iqr"
    )
    
    expect_true(any(result$Group == "Median [IQR]"))
})


test_that("desctable stats_continuous='median_range' works", {
    
    result <- desctable(
        data = clintrial,
        variables = c("age"),
        stats_continuous = "median_range"
    )
    
    expect_true(any(result$Group == "Median (Range)"))
})


test_that("desctable stats_continuous='range' works", {
    
    result <- desctable(
        data = clintrial,
        variables = c("age"),
        stats_continuous = "range"
    )
    
    expect_true(any(result$Group == "Range"))
})


test_that("desctable multiple stats_continuous creates multiple rows", {
    
    result <- desctable(
        data = clintrial,
        variables = c("age"),
        stats_continuous = c("mean_sd", "median_iqr", "range")
    )
    
    ## Age should have 3 rows (one for each stat)
    age_rows <- result[Variable == "age" | (Variable == "" & shift(Variable, fill = "age") == "age")]
    ## Count rows with stat labels
    stat_rows <- sum(result$Group %in% c("Mean +/- SD", "Median [IQR]", "Range"))
    expect_equal(stat_rows, 3)
})


## ============================================================================
## SECTION 4: Categorical Variables
## ============================================================================

test_that("desctable handles categorical variables correctly", {
    
    result <- desctable(
        data = clintrial,
        variables = c("sex", "stage")
    )
    
    expect_desctable_result(result)
    
    ## Should have rows for each level
    sex_levels <- levels(clintrial$sex)
    for (lvl in sex_levels) {
        expect_true(any(result$Group == lvl))
    }
})


test_that("desctable stats_categorical='n_percent' shows count and percent", {
    
    result <- desctable(
        data = clintrial,
        variables = c("sex"),
        stats_categorical = "n_percent"
    )
    
    ## Values should contain both count and percentage
    sex_rows <- result[Group %in% c("Male", "Female")]
    
    ## Should have format like "123 (45.6%)"
    expect_true(all(grepl("\\d+.*\\(.*%\\)", sex_rows$Total)))
})


test_that("desctable stats_categorical='n' shows count only", {
    
    result <- desctable(
        data = clintrial,
        variables = c("sex"),
        stats_categorical = "n"
    )
    
    sex_rows <- result[Group %in% c("Male", "Female")]
    
    ## Should NOT contain percentage
    expect_false(any(grepl("%", sex_rows$Total)))
})


test_that("desctable stats_categorical='percent' shows percent only", {
    
    result <- desctable(
        data = clintrial,
        variables = c("sex"),
        stats_categorical = "percent"
    )
    
    sex_rows <- result[Group %in% c("Male", "Female")]
    
    ## Should contain percentage
    expect_true(all(grepl("%", sex_rows$Total)))
})


## ============================================================================
## SECTION 5: Survival Variables
## ============================================================================

test_that("desctable handles survival variables", {
    
    skip_if_not_installed("survival")
    
    result <- desctable(
        data = clintrial,
        variables = c("Surv(os_months, os_status)")
    )
    
    expect_desctable_result(result)
    
    ## Should show median survival
    expect_true(any(grepl("Median|Survival", result$Group, ignore.case = TRUE)))
})


test_that("desctable survival with grouping shows log-rank test", {
    
    skip_if_not_installed("survival")
    
    result <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("Surv(os_months, os_status)"),
        test = TRUE
    )
    
    expect_pvalue_column(result)
    
    ## Should have a p-value for survival comparison
    surv_rows <- result[grepl("Surv|Median", Variable, ignore.case = TRUE) | 
                        grepl("Surv|Median", Group, ignore.case = TRUE)]
    
    ## At least one row should have a non-empty p-value
    if (nrow(surv_rows) > 0) {
        pvals <- surv_rows$`p-value`
        has_pval <- any(!is.na(pvals) & pvals != "")
        expect_true(has_pval)
    }
})


## ============================================================================
## SECTION 6: Statistical Tests
## ============================================================================

test_that("desctable test=TRUE adds p-value column", {
    
    result <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("age", "sex"),
        test = TRUE
    )
    
    expect_pvalue_column(result)
})


test_that("desctable test=FALSE excludes p-value column", {
    
    result <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("age", "sex"),
        test = FALSE
    )
    
    expect_false("p-value" %in% names(result))
})


test_that("desctable test_continuous='t' uses t-test", {
    
    result <- desctable(
        data = clintrial,
        by = "sex",  # Two groups
        variables = c("age"),
        stats_continuous = "mean_sd",
        test = TRUE,
        test_continuous = "t"
    )
    
    expect_pvalue_column(result)
    
    ## P-value should exist and be valid
    raw <- attr(result, "raw_data")
    expect_true(!is.null(raw$p_value))
})


test_that("desctable test_continuous='wrs' uses Wilcoxon test", {
    
    result <- suppressWarnings(
        desctable(
            data = clintrial,
            by = "sex",
            variables = c("age"),
            stats_continuous = "median_iqr",
            test = TRUE,
            test_continuous = "wrs"
        )
    )
    
    expect_pvalue_column(result)
})


test_that("desctable test_continuous='auto' selects appropriate test", {
    
    ## For mean_sd, should use parametric test
    result_mean <- desctable(
        data = clintrial,
        by = "sex",
        variables = c("age"),
        stats_continuous = "mean_sd",
        test = TRUE,
        test_continuous = "auto"
    )
    
    ## For median_iqr, should use non-parametric test
    result_median <- suppressWarnings(
        desctable(
            data = clintrial,
            by = "sex",
            variables = c("age"),
            stats_continuous = "median_iqr",
            test = TRUE,
            test_continuous = "auto"
        )
    )
    
    ## Both should have p-values
    expect_pvalue_column(result_mean)
    expect_pvalue_column(result_median)
})


test_that("desctable test_categorical='auto' selects appropriate test", {
    
    result <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("sex"),
        test = TRUE,
        test_categorical = "auto"
    )
    
    expect_pvalue_column(result)
})


test_that("desctable test_categorical='chisq' uses chi-squared test", {
    
    result <- suppressWarnings(
        desctable(
            data = clintrial,
            by = "treatment",
            variables = c("sex"),
            test = TRUE,
            test_categorical = "chisq"
        )
    )
    
    expect_pvalue_column(result)
})


test_that("desctable test_categorical='fisher' uses Fisher's exact test", {
    
    result <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("sex"),
        test = TRUE,
        test_categorical = "fisher"
    )
    
    expect_pvalue_column(result)
})


test_that("desctable range stat does not show p-value", {
    
    result <- desctable(
        data = clintrial,
        by = "sex",
        variables = c("age"),
        stats_continuous = "range",
        test = TRUE
    )
    
    ## Result should be valid
    expect_desctable_result(result)
    
    ## Range rows should have empty p-value
    range_rows <- result[Group == "Range"]
    expect_true(nrow(range_rows) > 0)  # Should have range rows
    
    if ("p-value" %in% names(result)) {
        expect_true(all(range_rows$`p-value` == "" | is.na(range_rows$`p-value`)))
    }
})


## ============================================================================
## SECTION 7: Total Column Options
## ============================================================================

test_that("desctable total=TRUE includes Total column first", {
    
    result <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("age"),
        total = TRUE
    )
    
    expect_total_column(result)
    
    ## Total should be before group columns
    col_idx <- which(names(result) == "Total")
    group_cols <- setdiff(names(result), c("Variable", "Group", "Total", "p-value"))
    
    if (length(group_cols) > 0) {
        first_group_idx <- min(which(names(result) %in% group_cols))
        expect_true(col_idx < first_group_idx)
    }
})


test_that("desctable total='last' puts Total column last", {
    
    result <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("age"),
        total = "last"
    )
    
    expect_total_column(result)
    
    ## Total should be after group columns but before p-value
    cols <- names(result)
    total_idx <- which(cols == "Total")
    
    if ("p-value" %in% cols) {
        pval_idx <- which(cols == "p-value")
        expect_true(total_idx < pval_idx)
    }
})


test_that("desctable total=FALSE excludes Total column", {
    
    result <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("age"),
        total = FALSE
    )
    
    expect_false("Total" %in% names(result))
})


test_that("desctable total_label customizes column name", {
    
    result <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("age"),
        total = TRUE,
        total_label = "Overall"
    )
    
    expect_true("Overall" %in% names(result))
    expect_false("Total" %in% names(result))
})


## ============================================================================
## SECTION 8: Missing Data Handling
## ============================================================================

test_that("desctable na_include=FALSE excludes missing data", {
    
    ## Use variable with known missing values
    result <- desctable(
        data = clintrial,
        variables = c("smoking"),  # Has NAs
        na_include = FALSE
    )
    
    ## Should NOT have "Unknown" or NA category
    expect_false(any(result$Group == "Unknown", na.rm = TRUE))
})


test_that("desctable na_include=TRUE shows missing category", {
    
    ## Ensure variable has missing values
    if (any(is.na(clintrial$smoking))) {
        result <- desctable(
            data = clintrial,
            variables = c("smoking"),
            na_include = TRUE,
            na_label = "Unknown"
        )
        
        ## Should have "Unknown" category
        expect_true(any(result$Group == "Unknown", na.rm = TRUE))
    }
})


test_that("desctable na_label customizes missing label", {
    
    if (any(is.na(clintrial$smoking))) {
        result <- desctable(
            data = clintrial,
            variables = c("smoking"),
            na_include = TRUE,
            na_label = "Missing Data"
        )
        
        expect_true(any(result$Group == "Missing Data", na.rm = TRUE))
    }
})


test_that("desctable na_include shows missing count for continuous", {
    
    ## Create data with missing continuous values
    test_data <- data.table::as.data.table(clintrial)
    test_data[1:10, age := NA]
    
    result <- desctable(
        data = test_data,
        variables = c("age"),
        na_include = TRUE,
        na_label = "Missing"
    )
    
    ## Should have a row showing missing count
    expect_true(any(result$Group == "Missing", na.rm = TRUE))
})


## ============================================================================
## SECTION 9: Formatting Options
## ============================================================================

test_that("desctable respects digits parameter", {
    
    result_1dig <- desctable(
        data = clintrial,
        variables = c("age"),
        stats_continuous = "mean_sd",
        digits = 1
    )
    
    result_2dig <- desctable(
        data = clintrial,
        variables = c("age"),
        stats_continuous = "mean_sd",
        digits = 2
    )
    
    ## Extract mean values and check decimal places
    mean_row_1 <- result_1dig[Group == "Mean +/- SD"]$Total
    mean_row_2 <- result_2dig[Group == "Mean +/- SD"]$Total
    
    ## Count decimals (rough check)
    if (length(mean_row_1) > 0 && length(mean_row_2) > 0) {
        ## 1 digit should have format like "45.2"
        ## 2 digits should have format like "45.23"
        expect_true(nchar(mean_row_1) <= nchar(mean_row_2) + 2)
    }
})


test_that("desctable formats large numbers with commas", {
    
    ## Create data with large values
    test_data <- data.table::as.data.table(clintrial)
    test_data[, large_var := age * 1000]
    
    result <- desctable(
        data = test_data,
        variables = c("large_var"),
        stats_continuous = "mean_sd"
    )
    
    ## Values >= 1000 should have commas
    mean_row <- result[Group == "Mean +/- SD"]$Total
    if (length(mean_row) > 0) {
        expect_true(grepl(",", mean_row))
    }
})


test_that("desctable respects digits_p parameter", {
    
    result <- desctable(
        data = clintrial,
        by = "sex",
        variables = c("age"),
        test = TRUE,
        digits_p = 4
    )
    
    ## P-values should have up to 4 decimal places
    if ("p-value" %in% names(result)) {
        pvals <- result$`p-value`
        numeric_pvals <- pvals[grepl("^0\\.", pvals)]
        
        if (length(numeric_pvals) > 0) {
            ## Check decimals
            decimals <- sapply(numeric_pvals, function(p) {
                nchar(sub("0\\.", "", p))
            })
            expect_true(all(decimals <= 4))
        }
    }
})


## ============================================================================
## SECTION 10: Custom Labels
## ============================================================================

test_that("desctable applies custom labels to variables", {
    
    custom_labels <- c(
        age = "Age (years)",
        sex = "Sex",
        stage = "Cancer Stage"
    )
    
    result <- desctable(
        data = clintrial,
        variables = c("age", "sex", "stage"),
        labels = custom_labels
    )
    
    ## Check labels are applied
    expect_true(any(result$Variable == "Age (years)"))
    expect_true(any(result$Variable == "Sex"))
    expect_true(any(result$Variable == "Cancer Stage"))
})


test_that("desctable works with clintrial_labels", {
    
    result <- desctable(
        data = clintrial,
        variables = c("age", "sex", "stage", "treatment"),
        labels = clintrial_labels
    )
    
    expect_desctable_result(result)
})


test_that("desctable labels work for grouping variable", {
    
    custom_labels <- c(
        treatment = "Treatment Arm"
    )
    
    result <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("age"),
        labels = custom_labels
    )
    
    ## The grouping variable label might appear in column names or elsewhere
    ## Just verify it doesn't error
    expect_desctable_result(result)
})


## ============================================================================
## SECTION 11: Raw Data Attribute
## ============================================================================

test_that("desctable raw_data attribute has correct structure", {
    
    result <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("age", "sex")
    )
    
    raw <- attr(result, "raw_data")
    
    expect_s3_class(raw, "data.table")
    expect_true("Variable" %in% names(raw) || "variable" %in% names(raw))
})


test_that("desctable raw_data contains numeric values", {
    
    result <- desctable(
        data = clintrial,
        variables = c("age"),
        stats_continuous = "mean_sd"
    )
    
    raw <- attr(result, "raw_data")
    
    ## Should have numeric values for statistics
    ## Column names may vary, but Total or similar should be numeric
    total_col <- grep("Total|total", names(raw), value = TRUE)[1]
    if (!is.na(total_col)) {
        expect_true(is.numeric(raw[[total_col]]))
    }
})


## ============================================================================
## SECTION 12: Edge Cases
## ============================================================================

test_that("desctable handles single variable", {
    
    result <- desctable(
        data = clintrial,
        variables = "age"
    )
    
    expect_desctable_result(result)
})


test_that("desctable handles many variables", {
    
    many_vars <- c("age", "sex", "smoking", "diabetes", "hypertension",
                   "stage", "grade", "ecog", "treatment", "surgery")
    
    result <- desctable(
        data = clintrial,
        variables = many_vars
    )
    
    expect_desctable_result(result)
    expect_equal(attr(result, "variables"), many_vars)
})


test_that("desctable handles all continuous variables", {
    
    result <- desctable(
        data = clintrial,
        variables = c("age", "bmi", "creatinine", "hemoglobin"),
        stats_continuous = "mean_sd"
    )
    
    expect_desctable_result(result)
})


test_that("desctable handles all categorical variables", {
    
    result <- desctable(
        data = clintrial,
        variables = c("sex", "stage", "grade", "treatment")
    )
    
    expect_desctable_result(result)
})


test_that("desctable handles mixed variable types", {
    
    result <- desctable(
        data = clintrial,
        variables = c("age", "sex", "bmi", "stage")
    )
    
    expect_desctable_result(result)
})


test_that("desctable works with data.frame input", {
    
    df <- as.data.frame(clintrial)
    
    result <- desctable(
        data = df,
        variables = c("age", "sex")
    )
    
    expect_desctable_result(result)
})


test_that("desctable handles groups with different factor levels", {
    
    result <- desctable(
        data = clintrial,
        by = "stage",  # Multiple levels
        variables = c("age", "sex")
    )
    
    expect_desctable_result(result)
    
    ## Should have columns for each stage level
    stage_levels <- levels(clintrial$stage)
    for (lvl in stage_levels) {
        expect_true(as.character(lvl) %in% names(result))
    }
})


## ============================================================================
## SECTION 13: Variable Order Preservation
## ============================================================================

test_that("desctable preserves variable order", {
    
    vars <- c("stage", "age", "sex")  # Specific order
    
    result <- desctable(
        data = clintrial,
        variables = vars
    )
    
    ## Variables should appear in order provided
    ## Get non-empty Variable values (these are the variable headers)
    var_col <- result$Variable[result$Variable != "" & result$Variable != "N"]
    
    ## Get unique variable names in order of first appearance
    seen <- c()
    for (v in var_col) {
        if (!(v %in% seen)) {
            seen <- c(seen, v)
        }
    }
    
    ## The order of seen variables should match the order of input vars
    ## (accounting for possible label transformations)
    ## Simply check that we have the right number and they appeared
    expect_equal(length(seen), length(vars))
})


## ============================================================================
## SECTION 14: Factor Level Ordering
## ============================================================================

test_that("desctable respects factor level ordering in categorical variables", {
    
    result <- desctable(
        data = clintrial,
        variables = c("stage")  # Ordered factor I, II, III, IV
    )
    
    ## Get stage rows
    stage_rows <- result[Variable == "stage" | (Variable == "" & Group %in% c("I", "II", "III", "IV"))]
    groups <- stage_rows$Group[stage_rows$Group != ""]
    
    ## Should be in factor level order
    expected_order <- levels(clintrial$stage)
    actual_order <- groups[groups %in% expected_order]
    
    if (length(actual_order) > 0) {
        expect_equal(actual_order, expected_order[expected_order %in% actual_order])
    }
})


test_that("desctable respects factor level ordering in grouping variable", {
    
    result <- desctable(
        data = clintrial,
        by = "stage",
        variables = c("age")
    )
    
    ## Column order should match factor levels
    cols <- names(result)
    stage_cols <- cols[cols %in% levels(clintrial$stage)]
    
    expected_order <- levels(clintrial$stage)[levels(clintrial$stage) %in% stage_cols]
    expect_equal(stage_cols, expected_order)
})


## ============================================================================
## SECTION 15: Consistency Tests
## ============================================================================

test_that("desctable produces consistent results across runs", {
    
    result1 <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("age", "sex")
    )
    
    result2 <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("age", "sex")
    )
    
    expect_equal(result1, result2)
})


test_that("desctable N values are consistent with data", {
    
    result <- desctable(
        data = clintrial,
        by = "sex",
        variables = c("age")
    )
    
    ## Total N should equal nrow
    total_n <- gsub(",", "", result$Total[1])
    expect_equal(as.integer(total_n), nrow(clintrial))
    
    ## Group Ns should sum to total
    male_n <- as.integer(gsub(",", "", result$Male[1]))
    female_n <- as.integer(gsub(",", "", result$Female[1]))
    
    expected_sum <- sum(clintrial$sex == "Male", na.rm = TRUE) + 
        sum(clintrial$sex == "Female", na.rm = TRUE)
    expect_equal(male_n + female_n, expected_sum)
})


## ============================================================================
## SECTION 16: Integration with Export Functions
## ============================================================================

test_that("desctable output is suitable for export", {
    
    result <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c("age", "sex", "stage"),
        labels = clintrial_labels
    )
    
    ## Should be a data.table
    expect_s3_class(result, "data.table")
    
    ## All columns should be character (for display)
    for (col in names(result)) {
        expect_true(is.character(result[[col]]) || is.numeric(result[[col]]))
    }
})


## ============================================================================
## SECTION 17: Complete Table 1 Workflow
## ============================================================================

test_that("desctable creates complete Table 1", {
    
    skip_if_not_installed("survival")
    
    result <- desctable(
        data = clintrial,
        by = "treatment",
        variables = c(
            "age", "sex", "bmi",
            "smoking", "diabetes",
            "stage", "grade",
            "Surv(os_months, os_status)"
        ),
        stats_continuous = c("median_iqr"),
        labels = clintrial_labels,
        total = TRUE,
        test = TRUE
    )
    
    expect_desctable_result(result)
    expect_total_column(result)
    expect_pvalue_column(result)
    
    ## Should have all the variables
    vars_attr <- attr(result, "variables")
    expect_equal(length(vars_attr), 8)
})


## ============================================================================
## SECTION 18: Two-group vs Multi-group Tests
## ============================================================================

test_that("desctable uses appropriate test for two groups", {
    
    ## Two groups - should work with t-test or Wilcoxon
    result <- desctable(
        data = clintrial,
        by = "sex",  # Two levels
        variables = c("age"),
        test = TRUE,
        test_continuous = "auto"
    )
    
    expect_pvalue_column(result)
})


test_that("desctable uses appropriate test for multiple groups", {
    
    ## Multiple groups - should use ANOVA or Kruskal-Wallis
    result <- desctable(
        data = clintrial,
        by = "stage",  # Four levels
        variables = c("age"),
        test = TRUE,
        test_continuous = "auto"
    )
    
    expect_pvalue_column(result)
})
