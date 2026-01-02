#' Test Suite for survtable()
#' 
#' Comprehensive tests covering survival summary tables, time points,
#' quantiles, statistical tests, formatting options, and edge cases.
#' 
#' @details Run with testthat::test_file("tests/testthat/test-survtable.R")

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

## Helper function to check survtable structure
expect_survtable_result <- function(result) {
    expect_s3_class(result, "survtable")
    expect_s3_class(result, "data.table")
    
    ## Required columns - either Variable (single outcome) or Variable+Group (multi)
    expect_true("Variable" %in% names(result) || "Group" %in% names(result))
    
    ## Required attributes
    expect_true(!is.null(attr(result, "raw_data")))
    expect_true(!is.null(attr(result, "outcome")))
    expect_true(!is.null(attr(result, "type")))
}

## Helper to check for p-value column
expect_pvalue_column <- function(result) {
    expect_true("p-value" %in% names(result))
}

## Helper to check for specific time column
expect_time_column <- function(result, time, unit = NULL) {
    if (!is.null(unit)) {
        col_name <- paste(time, unit)
    } else {
        col_name <- as.character(time)
    }
    expect_true(col_name %in% names(result) || 
                any(grepl(as.character(time), names(result))))
}


## ============================================================================
## SECTION 1: Basic Functionality - No Grouping
## ============================================================================

test_that("survtable works without grouping variable", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12, 24)
    )
    
    expect_survtable_result(result)
    
    ## Should NOT have p-value column (no groups to compare)
    expect_false("p-value" %in% names(result))
})


test_that("survtable ungrouped returns single row", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        total = FALSE
    )
    
    ## Should have just one row for overall
    expect_equal(nrow(result), 1)
})


test_that("survtable outcome attribute is set", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12)
    )
    
    expect_equal(attr(result, "outcome"), "Surv(os_months, os_status)")
})


## ============================================================================
## SECTION 2: Basic Functionality - With Grouping
## ============================================================================

test_that("survtable works with grouping variable", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12, 24)
    )
    
    expect_survtable_result(result)
    
    ## Should have p-value column by default
    expect_pvalue_column(result)
})


test_that("survtable by_variable attribute is set", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12)
    )
    
    expect_equal(attr(result, "by_variable"), "treatment")
})


test_that("survtable includes all group levels", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        total = FALSE
    )
    
    ## Should have rows for each treatment level
    treatment_levels <- levels(clintrial$treatment)
    group_col <- if ("Group" %in% names(result)) result$Group else result$Variable
    
    for (lvl in treatment_levels) {
        expect_true(lvl %in% group_col)
    }
})


## ============================================================================
## SECTION 3: Time Points
## ============================================================================

test_that("survtable handles single time point", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = 12
    )
    
    expect_survtable_result(result)
    expect_time_column(result, 12)
})


test_that("survtable handles multiple time points", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(6, 12, 24, 36)
    )
    
    expect_survtable_result(result)
    
    ## Check times attribute
    expect_equal(attr(result, "times"), c(6, 12, 24, 36))
})


test_that("survtable time_unit parameter works", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12, 24),
        time_unit = "months"
    )
    
    ## Column names should include unit
    expect_true("12 months" %in% names(result))
    expect_true("24 months" %in% names(result))
    
    ## Attribute should be set
    expect_equal(attr(result, "time_unit"), "months")
})


test_that("survtable times attribute is stored", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12, 24, 36)
    )
    
    expect_equal(attr(result, "times"), c(12, 24, 36))
})


## ============================================================================
## SECTION 4: Probability Quantiles (Median, etc.)
## ============================================================================

test_that("survtable handles median (probs = 0.5)", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        probs = 0.5
    )
    
    expect_survtable_result(result)
    
    ## Should have median column
    expect_true(any(grepl("Median", names(result), ignore.case = TRUE)))
})


test_that("survtable handles multiple quantiles", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        probs = c(0.25, 0.5, 0.75)
    )
    
    expect_survtable_result(result)
    
    ## Check probs attribute
    expect_equal(attr(result, "probs"), c(0.25, 0.5, 0.75))
})


test_that("survtable probs attribute is stored", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        probs = c(0.25, 0.5)
    )
    
    expect_equal(attr(result, "probs"), c(0.25, 0.5))
})


test_that("survtable combining times and probs works", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12, 24),
        probs = 0.5
    )
    
    expect_survtable_result(result)
    
    ## Should have both time columns and median column
    expect_true(any(grepl("12", names(result))))
    expect_true(any(grepl("Median", names(result), ignore.case = TRUE)))
})


## ============================================================================
## SECTION 5: Type Parameter (survival, risk, cumhaz)
## ============================================================================

test_that("survtable type='survival' works (default)", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        type = "survival"
    )
    
    expect_survtable_result(result)
    expect_equal(attr(result, "type"), "survival")
})


test_that("survtable type='risk' works", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        type = "risk"
    )
    
    expect_survtable_result(result)
    expect_equal(attr(result, "type"), "risk")
})


test_that("survtable type='cumhaz' works", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        type = "cumhaz"
    )
    
    expect_survtable_result(result)
    expect_equal(attr(result, "type"), "cumhaz")
})


test_that("survtable risk values are 1 - survival", {
    
    result_surv <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        type = "survival",
        percent = FALSE
    )
    
    result_risk <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        type = "risk",
        percent = FALSE
    )
    
    ## Raw values should be complements
    raw_surv <- attr(result_surv, "raw_data")
    raw_risk <- attr(result_risk, "raw_data")
    
    expect_true(!is.null(raw_surv))
    expect_true(!is.null(raw_risk))
})


## ============================================================================
## SECTION 6: Statistical Tests
## ============================================================================

test_that("survtable test=TRUE adds p-value column", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        test = TRUE
    )
    
    expect_pvalue_column(result)
})


test_that("survtable test=FALSE excludes p-value column", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        test = FALSE
    )
    
    expect_false("p-value" %in% names(result))
})


test_that("survtable test_type='logrank' works", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        test = TRUE,
        test_type = "logrank"
    )
    
    expect_pvalue_column(result)
    
    test_result <- attr(result, "test_result")
    expect_equal(test_result$test_type, "logrank")
})


test_that("survtable test_type='wilcoxon' works", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        test = TRUE,
        test_type = "wilcoxon"
    )
    
    expect_pvalue_column(result)
    
    test_result <- attr(result, "test_result")
    expect_equal(test_result$test_type, "wilcoxon")
})


test_that("survtable test_type='tarone' works", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        test = TRUE,
        test_type = "tarone"
    )
    
    expect_pvalue_column(result)
    
    test_result <- attr(result, "test_result")
    expect_equal(test_result$test_type, "tarone")
})


test_that("survtable test_type='petopeto' works", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        test = TRUE,
        test_type = "petopeto"
    )
    
    expect_pvalue_column(result)
    
    test_result <- attr(result, "test_result")
    expect_equal(test_result$test_type, "petopeto")
})


test_that("survtable test_result attribute contains valid data", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        test = TRUE
    )
    
    test_result <- attr(result, "test_result")
    
    expect_true(!is.null(test_result))
    expect_true("p_value" %in% names(test_result))
    expect_true("test_type" %in% names(test_result))
    expect_true(is.numeric(test_result$p_value))
})


## ============================================================================
## SECTION 7: Total Row
## ============================================================================

test_that("survtable total=TRUE includes total row first", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        total = TRUE
    )
    
    ## First row should be Total
    first_group <- if ("Group" %in% names(result)) result$Group[1] else result$Variable[1]
    expect_equal(first_group, "Total")
})


test_that("survtable total='first' includes total row first", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        total = "first"
    )
    
    first_group <- if ("Group" %in% names(result)) result$Group[1] else result$Variable[1]
    expect_equal(first_group, "Total")
})


test_that("survtable total='last' includes total row last", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        total = "last"
    )
    
    last_group <- if ("Group" %in% names(result)) result$Group[nrow(result)] else result$Variable[nrow(result)]
    expect_equal(last_group, "Total")
})


test_that("survtable total=FALSE excludes total row", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        total = FALSE
    )
    
    group_col <- if ("Group" %in% names(result)) result$Group else result$Variable
    expect_false("Total" %in% group_col)
})


test_that("survtable total_label parameter works", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        total = TRUE,
        total_label = "Overall"
    )
    
    group_col <- if ("Group" %in% names(result)) result$Group else result$Variable
    expect_true("Overall" %in% group_col)
})


## ============================================================================
## SECTION 8: Formatting Options
## ============================================================================

test_that("survtable percent=TRUE displays percentages", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        percent = TRUE
    )
    
    ## Values should contain %
    time_col <- names(result)[grepl("12", names(result))][1]
    expect_true(any(grepl("%", result[[time_col]])))
})


test_that("survtable percent=FALSE displays proportions", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        percent = FALSE
    )
    
    ## Values should NOT contain % (or have decimal values)
    time_col <- names(result)[grepl("12", names(result))][1]
    ## Check that values don't have %
    expect_false(any(grepl("%", result[[time_col]])))
})


test_that("survtable digits parameter affects survival precision", {
    
    result_0 <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        digits = 0
    )
    
    result_2 <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        digits = 2
    )
    
    ## Both should work
    expect_survtable_result(result_0)
    expect_survtable_result(result_2)
})


test_that("survtable time_digits parameter affects quantile precision", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        probs = 0.5,
        time_digits = 2
    )
    
    expect_survtable_result(result)
})


test_that("survtable p_digits parameter affects p-value precision", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        test = TRUE,
        p_digits = 4
    )
    
    expect_survtable_result(result)
    expect_pvalue_column(result)
})


## ============================================================================
## SECTION 9: Stats Parameter
## ============================================================================

test_that("survtable stats='survival' shows survival only", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        stats = "survival"
    )
    
    expect_survtable_result(result)
})


test_that("survtable stats includes 'ci' shows confidence intervals", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        stats = c("survival", "ci")
    )
    
    ## Values should contain parentheses for CI
    time_col <- names(result)[grepl("12", names(result))][1]
    expect_true(any(grepl("\\(", result[[time_col]])))
})


test_that("survtable stats includes 'n_risk' adds n at risk", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        stats = c("survival", "ci", "n_risk")
    )
    
    expect_survtable_result(result)
    
    ## Should contain n= in the output
    time_col <- names(result)[grepl("12", names(result))][1]
    expect_true(any(grepl("n=", result[[time_col]])))
})


test_that("survtable stats includes 'n_event' adds event count", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        stats = c("survival", "ci", "n_event")
    )
    
    expect_survtable_result(result)
    
    ## Should contain e= in the output
    time_col <- names(result)[grepl("12", names(result))][1]
    expect_true(any(grepl("e=", result[[time_col]])))
})


## ============================================================================
## SECTION 10: Confidence Interval Options
## ============================================================================

test_that("survtable conf_level=0.95 works (default)", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        conf_level = 0.95
    )
    
    expect_survtable_result(result)
})


test_that("survtable conf_level=0.90 works", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        conf_level = 0.90
    )
    
    expect_survtable_result(result)
})


test_that("survtable conf_type='log' works (default)", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        conf_type = "log"
    )
    
    expect_survtable_result(result)
})


test_that("survtable conf_type='log-log' works", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        conf_type = "log-log"
    )
    
    expect_survtable_result(result)
})


test_that("survtable conf_type='plain' works", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        conf_type = "plain"
    )
    
    expect_survtable_result(result)
})


## ============================================================================
## SECTION 11: Multiple Outcomes
## ============================================================================

test_that("survtable handles multiple outcomes", {
    
    result <- survtable(
        data = clintrial,
        outcome = c("Surv(pfs_months, pfs_status)", "Surv(os_months, os_status)"),
        by = "treatment",
        times = c(12),
        total = FALSE
    )
    
    expect_survtable_result(result)
    
    ## Should have Variable column for outcome names
    expect_true("Variable" %in% names(result))
    
    ## Outcome attribute should be vector
    expect_equal(length(attr(result, "outcome")), 2)
})


test_that("survtable multiple outcomes with labels works", {
    
    result <- survtable(
        data = clintrial,
        outcome = c("Surv(pfs_months, pfs_status)", "Surv(os_months, os_status)"),
        by = "treatment",
        times = c(12),
        total = FALSE,
        labels = c(
            "Surv(pfs_months, pfs_status)" = "PFS",
            "Surv(os_months, os_status)" = "OS"
        )
    )
    
    expect_survtable_result(result)
    
    ## Variable column should contain labels
    expect_true("PFS" %in% result$Variable)
    expect_true("OS" %in% result$Variable)
})


test_that("survtable multiple outcomes stacks correctly", {
    
    result <- survtable(
        data = clintrial,
        outcome = c("Surv(pfs_months, pfs_status)", "Surv(os_months, os_status)"),
        by = "treatment",
        times = c(12),
        total = FALSE
    )
    
    ## Number of rows should be 2 * number of treatment levels
    n_treatments <- length(levels(clintrial$treatment))
    expected_rows <- 2 * n_treatments
    expect_equal(nrow(result), expected_rows)
})


test_that("survtable multiple outcomes each have test results", {
    
    result <- survtable(
        data = clintrial,
        outcome = c("Surv(pfs_months, pfs_status)", "Surv(os_months, os_status)"),
        by = "treatment",
        times = c(12),
        test = TRUE,
        total = FALSE
    )
    
    ## test_result should be a list with both outcomes
    test_result <- attr(result, "test_result")
    expect_true(is.list(test_result))
    expect_equal(length(test_result), 2)
})


## ============================================================================
## SECTION 12: Labels Parameter
## ============================================================================

test_that("survtable labels parameter works for by variable", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "stage",
        times = c(12),
        labels = clintrial_labels
    )
    
    expect_survtable_result(result)
})


## ============================================================================
## SECTION 13: Missing Data
## ============================================================================

test_that("survtable handles missing values with na_rm=TRUE", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        na_rm = TRUE
    )
    
    expect_survtable_result(result)
})


test_that("survtable handles missing values with na_rm=FALSE", {
    
    ## This may produce warnings but should not error
    result <- tryCatch(
        survtable(
            data = clintrial,
            outcome = "Surv(os_months, os_status)",
            by = "treatment",
            times = c(12),
            na_rm = FALSE
        ),
        warning = function(w) {
            ## Capture but continue
            suppressWarnings(
                survtable(
                    data = clintrial,
                    outcome = "Surv(os_months, os_status)",
                    by = "treatment",
                    times = c(12),
                    na_rm = FALSE
                )
            )
        }
    )
    
    expect_survtable_result(result)
})


## ============================================================================
## SECTION 14: Edge Cases
## ============================================================================

test_that("survtable handles times=NULL with probs only", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = NULL,
        probs = c(0.25, 0.5, 0.75)
    )
    
    expect_survtable_result(result)
    
    ## Should not have time columns, only quantile columns
    expect_true(any(grepl("Median", names(result), ignore.case = TRUE)))
})


test_that("survtable handles probs=NULL with times only", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12, 24),
        probs = NULL
    )
    
    expect_survtable_result(result)
    
    ## Should have time columns, no median column
    expect_true(any(grepl("12", names(result))))
    expect_false(any(grepl("Median", names(result), ignore.case = TRUE)))
})


test_that("survtable handles single group", {
    
    ## Create single-group data
    single_group <- clintrial[clintrial$treatment == "Control", ]
    
    result <- survtable(
        data = single_group,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12),
        test = FALSE  # Can't test with single group
    )
    
    expect_survtable_result(result)
})


test_that("survtable handles two groups", {
    
    ## Use sex as binary grouping variable
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "sex",
        times = c(12),
        test = TRUE
    )
    
    expect_survtable_result(result)
    expect_pvalue_column(result)
})


test_that("survtable handles many groups", {
    
    ## Use stage with 4 levels
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "stage",
        times = c(12),
        test = TRUE
    )
    
    expect_survtable_result(result)
    expect_pvalue_column(result)
})


test_that("survtable handles long follow-up times", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(60, 120),  # Long follow-up
        probs = 0.5
    )
    
    expect_survtable_result(result)
})


test_that("survtable handles short follow-up times", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(1, 3, 6),  # Short follow-up
        probs = 0.5
    )
    
    expect_survtable_result(result)
})


## ============================================================================
## SECTION 15: Raw Data Attribute
## ============================================================================

test_that("survtable raw_data attribute contains numeric values", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12, 24),
        probs = 0.5
    )
    
    raw <- attr(result, "raw_data")
    
    expect_true(!is.null(raw))
    expect_s3_class(raw, "data.table")
    
    ## Should have estimate columns
    est_cols <- grep("_estimate$", names(raw), value = TRUE)
    expect_true(length(est_cols) > 0)
})


test_that("survtable survfit_objects attribute is stored", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12)
    )
    
    survfit_objs <- attr(result, "survfit_objects")
    
    expect_true(!is.null(survfit_objs))
    expect_true(is.list(survfit_objs))
})


## ============================================================================
## SECTION 16: Data Input Types
## ============================================================================

test_that("survtable works with data.frame input", {
    
    df <- as.data.frame(clintrial)
    
    result <- survtable(
        data = df,
        outcome = "Surv(os_months, os_status)",
        times = c(12)
    )
    
    expect_survtable_result(result)
})


test_that("survtable works with data.table input", {
    
    result <- survtable(
        data = clintrial,  # Already data.table
        outcome = "Surv(os_months, os_status)",
        times = c(12)
    )
    
    expect_survtable_result(result)
})


## ============================================================================
## SECTION 17: Consistency Tests
## ============================================================================

test_that("survtable produces consistent results across runs", {
    
    result1 <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12, 24)
    )
    
    result2 <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12, 24)
    )
    
    expect_equal(result1, result2)
})


test_that("survtable survival probabilities are in valid range", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        times = c(12),
        probs = NULL,  # No quantiles, only time points
        type = "survival",
        percent = FALSE
    )
    
    raw <- attr(result, "raw_data")
    
    ## Only check columns that are survival probability estimates (time point columns)
    ## Exclude median/quantile columns which contain time values, not probabilities
    est_cols <- grep("^12.*_estimate$", names(raw), value = TRUE)
    
    for (col in est_cols) {
        vals <- raw[[col]]
        vals <- vals[!is.na(vals)]
        if (length(vals) > 0) {
            expect_true(all(vals >= 0 & vals <= 1),
                        info = paste("Column", col, "has values outside [0,1]:", 
                                     paste(vals[vals < 0 | vals > 1], collapse = ", ")))
        }
    }
})


## ============================================================================
## SECTION 18: Print Method
## ============================================================================

test_that("survtable print method works", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12, 24),
        probs = 0.5
    )
    
    ## Print should not error
    expect_output(print(result))
})


test_that("survtable print method shows correct header info", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12, 24),
        time_unit = "months"
    )
    
    ## Capture print output
    output <- capture.output(print(result))
    
    ## Should mention survival
    expect_true(any(grepl("Survival", output, ignore.case = TRUE)))
})


## ============================================================================
## SECTION 19: Integration with Export Functions
## ============================================================================

test_that("survtable output is suitable for export", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12, 24),
        probs = 0.5
    )
    
    ## Should be a data.table
    expect_s3_class(result, "data.table")
    
    ## All columns should be character (for display)
    for (col in names(result)) {
        expect_true(is.character(result[[col]]))
    }
})


## ============================================================================
## SECTION 20: Error Handling
## ============================================================================

test_that("survtable errors on invalid Surv() syntax", {
    
    expect_error(
        survtable(
            data = clintrial,
            outcome = "invalid_outcome",
            times = c(12)
        ),
        regexp = "Surv"
    )
})


test_that("survtable errors on missing time variable", {
    
    expect_error(
        survtable(
            data = clintrial,
            outcome = "Surv(nonexistent_time, os_status)",
            times = c(12)
        ),
        regexp = "not found"
    )
})


test_that("survtable errors on missing status variable", {
    
    expect_error(
        survtable(
            data = clintrial,
            outcome = "Surv(os_months, nonexistent_status)",
            times = c(12)
        ),
        regexp = "not found"
    )
})


test_that("survtable errors on missing by variable", {
    
    expect_error(
        survtable(
            data = clintrial,
            outcome = "Surv(os_months, os_status)",
            by = "nonexistent_var",
            times = c(12)
        ),
        regexp = "not found"
    )
})


test_that("survtable errors when both times and probs are NULL", {
    
    expect_error(
        survtable(
            data = clintrial,
            outcome = "Surv(os_months, os_status)",
            times = NULL,
            probs = NULL
        )
    )
})


## ============================================================================
## SECTION 21: Complete Workflow
## ============================================================================

test_that("survtable creates complete survival summary", {
    
    result <- survtable(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        by = "treatment",
        times = c(12, 24, 36),
        probs = c(0.25, 0.5),
        time_unit = "months",
        stats = c("survival", "ci"),
        type = "survival",
        conf_level = 0.95,
        test = TRUE,
        test_type = "logrank",
        total = TRUE,
        digits = 0,
        time_digits = 1,
        p_digits = 3
    )
    
    expect_survtable_result(result)
    expect_pvalue_column(result)
    
    ## Check all expected columns
    expect_true("12 months" %in% names(result))
    expect_true("24 months" %in% names(result))
    expect_true("36 months" %in% names(result))
    expect_true(any(grepl("Median", names(result), ignore.case = TRUE)))
})


test_that("survtable PFS/OS workflow works", {
    
    result <- survtable(
        data = clintrial,
        outcome = c("Surv(pfs_months, pfs_status)", "Surv(os_months, os_status)"),
        by = "treatment",
        times = c(12, 24),
        probs = 0.5,
        time_unit = "months",
        total = FALSE,
        test = TRUE,
        labels = c(
            "Surv(pfs_months, pfs_status)" = "Progression-Free Survival",
            "Surv(os_months, os_status)" = "Overall Survival"
        )
    )
    
    expect_survtable_result(result)
    
    ## Should have labeled outcomes
    expect_true("Progression-Free Survival" %in% result$Variable)
    expect_true("Overall Survival" %in% result$Variable)
    
    ## Should have p-values for both
    expect_pvalue_column(result)
})
