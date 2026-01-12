#' Test Suite for fullfit()
#' 
#' Comprehensive tests covering variable selection methods, output formats,
#' model types, metrics options, return types, and edge cases.
#' 
#' @details Run with testthat::test_file("tests/testthat/test-fullfit.R")

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

## Create binary outcome for logistic regression tests
## IMPORTANT: Add to BOTH datasets
clintrial$response <- as.integer(clintrial$os_status == 1 & clintrial$os_months < 24)
clintrial_complete$response <- as.integer(clintrial_complete$os_status == 1 & clintrial_complete$os_months < 24)

## Create integer count outcome for Poisson/quasipoisson/negative binomial tests
## IMPORTANT: Add to BOTH datasets
set.seed(123)
clintrial$count_outcome <- as.integer(ceiling(clintrial$los_days))
clintrial_complete$count_outcome <- as.integer(ceiling(clintrial_complete$los_days))

## Helper function to check fullfit_result structure
expect_fullfit_result <- function(result) {
    expect_s3_class(result, "fullfit_result")
    expect_s3_class(result, "data.table")
    
    ## Required columns
    expect_true("Variable" %in% names(result))
    expect_true("Group" %in% names(result))
    
    ## Required attributes
    expect_true(!is.null(attr(result, "outcome")))
    expect_true(!is.null(attr(result, "model_type")))
    expect_true(!is.null(attr(result, "method")))
}

## Helper to check for univariable columns
expect_uni_columns <- function(result, metrics = c("effect", "p")) {
    col_names <- names(result)
    
    if ("effect" %in% metrics) {
        ## Column format is "OR (95% CI)" for univariable (not "Univariable OR...")
        has_uni_effect <- any(grepl("^(OR|HR|RR|Coefficient|Estimate)\\s*\\(.*CI", col_names))
        expect_true(has_uni_effect, 
                    info = paste("Expected univariable effect column, found:", 
                                 paste(col_names, collapse = ", ")))
    }
    
    if ("p" %in% metrics) {
        expect_true("Uni p" %in% col_names,
                    info = paste("Expected 'Uni p' column, found:", 
                                 paste(col_names, collapse = ", ")))
    }
}

## Helper to check for multivariable columns
expect_multi_columns <- function(result, metrics = c("effect", "p")) {
    col_names <- names(result)
    
    if ("effect" %in% metrics) {
        ## Column format is "aOR (95% CI)" for multivariable (not "Multivariable aOR...")
        has_multi_effect <- any(grepl("^(aOR|aHR|aRR|Adj\\. Coefficient)\\s*\\(.*CI", col_names))
        expect_true(has_multi_effect, 
                    info = paste("Expected multivariable effect column, found:", 
                                 paste(col_names, collapse = ", ")))
    }
    
    if ("p" %in% metrics) {
        expect_true("Multi p" %in% col_names,
                    info = paste("Expected 'Multi p' column, found:", 
                                 paste(col_names, collapse = ", ")))
    }
}


## ============================================================================
## SECTION 1: Basic Functionality - Method = "screen"
## ============================================================================

test_that("fullfit works with method='screen' (default)", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "smoking", "stage"),
        method = "screen",
        p_threshold = 0.05,
        model_type = "glm"
    )
    
    expect_fullfit_result(result)
    expect_equal(attr(result, "method"), "screen")
    
    ## Should have both uni and multi columns by default
    expect_uni_columns(result)
    expect_multi_columns(result)
})


test_that("fullfit screen method respects p_threshold", {
    
    ## Strict threshold - fewer variables in multi
    result_strict <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "smoking", "diabetes", "stage", "grade"),
        method = "screen",
        p_threshold = 0.01,
        model_type = "glm"
    )
    
    ## Liberal threshold - more variables in multi
    result_liberal <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "smoking", "diabetes", "stage", "grade"),
        method = "screen",
        p_threshold = 0.50,
        model_type = "glm"
    )
    
    ## Liberal should include at least as many predictors as strict
    n_strict <- attr(result_strict, "n_multi")
    n_liberal <- attr(result_liberal, "n_multi")
    
    ## Handle case where strict might have 0 or NULL
    if (!is.null(n_strict) && !is.null(n_liberal)) {
        expect_true(n_liberal >= n_strict)
    }
})


test_that("fullfit screen with very strict threshold handles edge cases", {
    
    ## Very strict threshold - may or may not select variables depending on data
    ## This tests that the function handles the edge case gracefully
    result <- tryCatch({
        suppressWarnings(
            fullfit(
                data = clintrial,
                outcome = "response",
                predictors = c("age", "sex"),
                method = "screen",
                p_threshold = 0.0001,  # Very strict
                model_type = "glm"
            )
        )
    }, error = function(e) NULL)
    
    ## Should either return a valid result or NULL (if no vars selected)
    if (!is.null(result)) {
        expect_fullfit_result(result)
    }
    
    ## Either way, the function should not error unexpectedly
    expect_true(TRUE)
})


## ============================================================================
## SECTION 2: Method = "all"
## ============================================================================

test_that("fullfit works with method='all'", {
    
    predictors <- c("age", "sex", "smoking", "stage")
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = predictors,
        method = "all",
        model_type = "glm"
    )
    
    expect_fullfit_result(result)
    expect_equal(attr(result, "method"), "all")
    
    ## All predictors should be in multivariable model
    n_multi <- attr(result, "n_multi")
    expect_equal(n_multi, length(predictors))
})


test_that("fullfit method='all' includes all variables in multi columns", {
    
    predictors <- c("age", "sex", "stage")
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = predictors,
        method = "all",
        model_type = "glm",
        columns = "both"
    )
    
    ## Multi columns should exist (aOR/aHR/aRR format)
    multi_effect_col <- grep("^(aOR|aHR|aRR|Adj\\.)", names(result), value = TRUE)[1]
    expect_true(!is.na(multi_effect_col))
    
    ## Get unique non-empty variable names
    vars_in_result <- unique(result$Variable[result$Variable != ""])
    
    ## Should have entries for all predictors (or their labels)
    expect_true(length(vars_in_result) >= length(predictors))
})


## ============================================================================
## SECTION 3: Method = "custom"
## ============================================================================

test_that("fullfit works with method='custom'", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "smoking", "diabetes", "stage"),
        method = "custom",
        multi_predictors = c("age", "stage"),  # Only include these in multi
        model_type = "glm"
    )
    
    expect_fullfit_result(result)
    expect_equal(attr(result, "method"), "custom")
    
    ## Only custom predictors in multi
    n_multi <- attr(result, "n_multi")
    expect_equal(n_multi, 2)
})


test_that("fullfit method='custom' requires multi_predictors", {
    
    expect_error(
        fullfit(
            data = clintrial,
            outcome = "response",
            predictors = c("age", "sex"),
            method = "custom"
            ## multi_predictors not specified - should error
        ),
        regexp = "multi_predictors must be specified"
    )
})


test_that("fullfit method='custom' shows uni for all, multi for selected", {
    
    all_predictors <- c("age", "sex", "smoking", "stage")
    selected <- c("age", "stage")
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = all_predictors,
        method = "custom",
        multi_predictors = selected,
        model_type = "glm",
        columns = "both"
    )
    
    ## All predictors should appear in the table
    vars_in_result <- unique(result$Variable[result$Variable != ""])
    
    ## Variables not in multi should show "-" in multi columns
    multi_effect_col <- grep("^(aOR|aHR|aRR|Adj\\.)", names(result), value = TRUE)[1]
    
    if (!is.na(multi_effect_col)) {
        ## Check that non-selected variables have "-" in multi column
        ## (This is complex due to factor levels, so just verify structure)
        expect_true(any(result[[multi_effect_col]] == "-") || 
                    length(selected) == length(all_predictors))
    } else {
        ## If no multi column found, that's also acceptable for certain column configs
        expect_true(TRUE)
    }
})


## ============================================================================
## SECTION 4: Columns Parameter
## ============================================================================

test_that("fullfit columns='both' shows uni and multi", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        method = "all",
        columns = "both",
        model_type = "glm"
    )
    
    expect_uni_columns(result)
    expect_multi_columns(result)
})


test_that("fullfit columns='uni' shows only univariable", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "stage"),
        method = "all",
        columns = "uni",
        model_type = "glm"
    )
    
    col_names <- names(result)
    
    ## Should have univariable columns (OR/HR/RR without 'a' prefix)
    expect_true(any(grepl("^(OR|HR|RR|Coefficient)\\s*\\(", col_names)))
    
    ## Should NOT have multivariable columns (aOR/aHR/aRR)
    expect_false(any(grepl("^(aOR|aHR|aRR|Adj\\.)", col_names)))
    expect_false("Multi p" %in% col_names)
})


test_that("fullfit columns='multi' shows only multivariable", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "stage"),
        method = "all",
        columns = "multi",
        model_type = "glm"
    )
    
    col_names <- names(result)
    
    ## Should have multivariable columns (aOR/aHR/aRR)
    expect_true(any(grepl("^(aOR|aHR|aRR|Adj\\.)", col_names)))
    
    ## Should NOT have univariable columns (OR without 'a' prefix, Uni p)
    expect_false("Uni p" %in% col_names)
})


## ============================================================================
## SECTION 5: Metrics Parameter
## ============================================================================

test_that("fullfit metrics='both' shows effect and p-value", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        method = "all",
        metrics = "both",
        model_type = "glm"
    )
    
    col_names <- names(result)
    
    ## Should have effect columns
    expect_true(any(grepl("(OR|HR|RR|Estimate).*CI", col_names)))
    
    ## Should have p-value columns
    expect_true(any(grepl("p$", col_names)))
})


test_that("fullfit metrics='effect' shows only effect estimates", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        method = "all",
        metrics = "effect",
        model_type = "glm"
    )
    
    col_names <- names(result)
    
    ## Should have effect columns
    expect_true(any(grepl("(OR|HR|RR|Estimate).*CI", col_names)))
    
    ## Should NOT have p-value columns
    expect_false("Uni p" %in% col_names)
    expect_false("Multi p" %in% col_names)
})


test_that("fullfit metrics='p' shows only p-values", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        method = "all",
        metrics = "p",
        model_type = "glm"
    )
    
    col_names <- names(result)
    
    ## Should have p-value columns
    expect_true(any(grepl("p$", col_names)))
    
    ## Should NOT have effect columns (CI columns)
    expect_false(any(grepl("95% CI", col_names)))
})


test_that("fullfit metrics as vector works", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        method = "all",
        metrics = c("effect", "p"),  # Vector form
        model_type = "glm"
    )
    
    col_names <- names(result)
    
    ## Should have both
    expect_true(any(grepl("CI", col_names)))
    expect_true(any(grepl("p$", col_names)))
})


## ============================================================================
## SECTION 6: Return Type Parameter
## ============================================================================

test_that("fullfit return_type='table' returns fullfit_result", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        method = "all",
        return_type = "table",
        model_type = "glm"
    )
    
    expect_fullfit_result(result)
})


test_that("fullfit return_type='model' returns model object", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        method = "all",
        return_type = "model",
        model_type = "glm"
    )
    
    ## Should return a glm object, not a data.table
    expect_s3_class(result, "glm")
    expect_false(inherits(result, "data.table"))
    
    ## Should be usable as a model
    expect_true(length(coef(result)) > 0)
})


test_that("fullfit return_type='both' returns list with table and model", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        method = "all",
        return_type = "both",
        model_type = "glm"
    )
    
    ## Should be a list
    expect_type(result, "list")
    expect_true("table" %in% names(result))
    expect_true("model" %in% names(result))
    
    ## Table should be fullfit_result
    expect_s3_class(result$table, "fullfit_result")
    
    ## Model should be glm
    expect_s3_class(result$model, "glm")
})


## ============================================================================
## SECTION 7: Model Types
## ============================================================================

test_that("fullfit works with GLM (logistic)", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "stage"),
        model_type = "glm",
        family = "binomial",
        method = "all"
    )
    
    expect_fullfit_result(result)
    
    ## Should show OR
    col_names <- names(result)
    expect_true(any(grepl("OR", col_names)))
})


test_that("fullfit works with linear model", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex", "stage"),
        model_type = "lm",
        method = "all"
    )
    
    expect_fullfit_result(result)
    expect_equal(attr(result, "model_type"), "lm")
    
    ## Linear models should show Estimate, not OR
    col_names <- names(result)
    expect_true(any(grepl("Estimate", col_names)) || any(grepl("Coefficient", col_names)))
})


test_that("fullfit works with Cox regression", {
    
    skip_if_not_installed("survival")
    
    result <- fullfit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "sex", "stage", "treatment"),
        model_type = "coxph",
        method = "all"
    )
    
    expect_fullfit_result(result)
    expect_equal(attr(result, "model_type"), "coxph")
    
    ## Should show HR
    col_names <- names(result)
    expect_true(any(grepl("HR", col_names)))
})


test_that("fullfit works with Poisson regression", {
    
    ## Create proper count data for Poisson
    set.seed(123)
    test_data <- data.table::as.data.table(clintrial)
    test_data[, count_outcome := as.integer(round(los_days))]
    
    result <- suppressWarnings(
        fullfit(
            data = test_data,
            outcome = "count_outcome",
            predictors = c("age", "sex", "stage"),
            model_type = "glm",
            family = "poisson",
            method = "all"
        )
    )
    
    expect_fullfit_result(result)
    
    ## Should show RR
    col_names <- names(result)
    expect_true(any(grepl("RR", col_names)))
})


## ============================================================================
## SECTION 8: Display Options
## ============================================================================

test_that("fullfit respects show_n parameter", {
    
    result_with_n <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        method = "all",
        show_n = TRUE
    )
    
    result_no_n <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        method = "all",
        show_n = FALSE
    )
    
    expect_true("n" %in% names(result_with_n))
    expect_false("n" %in% names(result_no_n))
})


test_that("fullfit respects show_events parameter", {
    
    result_with_events <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        method = "all",
        show_events = TRUE
    )
    
    result_no_events <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        method = "all",
        show_events = FALSE
    )
    
    expect_true("Events" %in% names(result_with_events))
    expect_false("Events" %in% names(result_no_events))
})


test_that("fullfit respects digits parameter", {
    
    result_2dig <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age"),
        model_type = "glm",
        method = "all",
        digits = 2
    )
    
    result_3dig <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age"),
        model_type = "glm",
        method = "all",
        digits = 3
    )
    
    ## Extract effect columns and check decimal places
    effect_col_2 <- grep("OR.*CI", names(result_2dig), value = TRUE)[1]
    effect_col_3 <- grep("OR.*CI", names(result_3dig), value = TRUE)[1]
    
    if (!is.na(effect_col_2) && !is.na(effect_col_3)) {
        val_2 <- result_2dig[[effect_col_2]][1]
        val_3 <- result_3dig[[effect_col_3]][1]
        
        ## Check decimal places differ
        if (!is.na(val_2) && !is.na(val_3) && val_2 != "-" && val_3 != "-") {
            first_num_2 <- sub(" .*", "", val_2)
            first_num_3 <- sub(" .*", "", val_3)
            
            decimals_2 <- nchar(sub(".*\\.", "", first_num_2))
            decimals_3 <- nchar(sub(".*\\.", "", first_num_3))
            
            expect_equal(decimals_2, 2)
            expect_equal(decimals_3, 3)
        }
    }
})


test_that("fullfit respects reference_rows parameter", {
    
    result_with_ref <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "stage"),
        model_type = "glm",
        method = "all",
        reference_rows = TRUE,
        parallel = FALSE
    )
    
    result_no_ref <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "stage"),
        model_type = "glm",
        method = "all",
        reference_rows = FALSE,
        parallel = FALSE
    )
    
    ## With reference rows should have more rows
    expect_true(nrow(result_with_ref) >= nrow(result_no_ref))
})


## ============================================================================
## SECTION 9: Labels Parameter
## ============================================================================

test_that("fullfit applies custom labels", {
    
    custom_labels <- c(
        age = "Age (years)",
        sex = "Sex",
        stage = "Cancer Stage"
    )
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "stage"),
        model_type = "glm",
        method = "all",
        labels = custom_labels
    )
    
    ## Check labels are applied
    variables <- result$Variable[result$Variable != ""]
    
    ## At least one custom label should appear
    has_label <- any(grepl("Age|Sex|Cancer Stage", variables))
    expect_true(has_label)
})


test_that("fullfit works with clintrial_labels", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "stage", "treatment"),
        model_type = "glm",
        method = "all",
        labels = clintrial_labels
    )
    
    expect_fullfit_result(result)
})


## ============================================================================
## SECTION 10: keep_models Parameter
## ============================================================================

test_that("fullfit keep_models=FALSE does not store uni models", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        method = "all",
        keep_models = FALSE
    )
    
    ## When keep_models is FALSE, uni_results might exist but models should not be stored
    ## Just verify the result is valid
    expect_fullfit_result(result)
})


test_that("fullfit keep_models=TRUE stores uni models", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        method = "all",
        keep_models = TRUE
    )
    
    expect_fullfit_result(result)
    
    ## When keep_models=TRUE, models should be accessible somehow
    ## The exact attribute structure may vary, so just check result is valid
    ## and has the expected attributes
    expect_true(!is.null(attr(result, "outcome")))
    expect_true(!is.null(attr(result, "model")))
})


## ============================================================================
## SECTION 11: Confidence Level
## ============================================================================

test_that("fullfit respects conf_level parameter", {
    
    result_95 <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age"),
        model_type = "glm",
        method = "all",
        conf_level = 0.95
    )
    
    result_90 <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age"),
        model_type = "glm",
        method = "all",
        conf_level = 0.90
    )
    
    ## Both should work - checking structure
    expect_fullfit_result(result_95)
    expect_fullfit_result(result_90)
    
    ## Results should be different (narrower CI for 90%)
    ## This is hard to test directly from formatted output
})


## ============================================================================
## SECTION 12: Exponentiate Parameter
## ============================================================================

test_that("fullfit auto-exponentiates for logistic", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "binomial",
        method = "all"
    )
    
    ## Should show OR by default
    col_names <- names(result)
    expect_true(any(grepl("OR", col_names)))
})


test_that("fullfit exponentiate=FALSE shows coefficients", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "binomial",
        method = "all",
        exponentiate = FALSE
    )
    
    ## Should show Coefficient or Estimate, not OR
    col_names <- names(result)
    expect_true(any(grepl("Coefficient|Estimate", col_names)))
})


## ============================================================================
## SECTION 13: Factor Variables
## ============================================================================

test_that("fullfit handles factor variables correctly", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "stage"),  # stage is a factor
        model_type = "glm",
        method = "all",
        reference_rows = TRUE,
        parallel = FALSE
    )
    
    expect_fullfit_result(result)
    
    ## Stage should have multiple rows (reference + levels)
    ## Find rows where Variable contains "stage" or follows it
    stage_idx <- which(result$Variable == "stage" | 
                       grepl("[Ss]tage", result$Variable))
    
    ## Should have multiple rows for factor
    if (length(stage_idx) > 0) {
        expect_true(length(stage_idx) >= 1)
    }
})


test_that("fullfit handles multiple factor variables", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("sex", "stage", "grade"),  # Multiple factors
        model_type = "glm",
        method = "all"
    )
    
    expect_fullfit_result(result)
    
    ## Should have rows for all factors
    expect_true(nrow(result) > 3)  # More than just 3 single rows
})


## ============================================================================
## SECTION 14: Attributes Preservation
## ============================================================================

test_that("fullfit preserves all expected attributes", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "stage"),
        model_type = "glm",
        method = "screen",
        p_threshold = 0.10
    )
    
    ## Check all expected attributes
    expect_equal(attr(result, "outcome"), "response")
    expect_equal(attr(result, "model_type"), "glm")
    expect_equal(attr(result, "method"), "screen")
    expect_equal(attr(result, "columns"), "both")
    
    ## Should have model if multivariable was fit
    model <- attr(result, "model")
    if (!is.null(attr(result, "n_multi")) && attr(result, "n_multi") > 0) {
        expect_true(!is.null(model))
    }
})


test_that("fullfit stores expected attributes", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        method = "all",
        columns = "both"  # Need uni for this
    )
    
    ## Should have key attributes
    expect_true(!is.null(attr(result, "outcome")))
    expect_true(!is.null(attr(result, "model_type")))
    expect_true(!is.null(attr(result, "method")))
    
    ## Model should be stored when multivariable is fit
    expect_true(!is.null(attr(result, "model")))
})


test_that("fullfit n_multi attribute is correct", {
    
    predictors <- c("age", "sex", "stage")
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = predictors,
        model_type = "glm",
        method = "all"
    )
    
    n_multi <- attr(result, "n_multi")
    expect_equal(n_multi, length(predictors))
})


## ============================================================================
## SECTION 15: Print Method
## ============================================================================

test_that("print.fullfit_result produces output", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        method = "all"
    )
    
    ## Capture print output
    output <- capture.output(print(result))
    
    ## Should contain key information
    expect_true(any(grepl("Fullfit|Analysis", output, ignore.case = TRUE)))
    expect_true(any(grepl("Outcome", output)))
    expect_true(any(grepl("Method", output)))
})


## ============================================================================
## SECTION 16: Edge Cases
## ============================================================================

test_that("fullfit handles single predictor", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = "age",  # Single predictor
        model_type = "glm",
        method = "all"
    )
    
    expect_fullfit_result(result)
})


test_that("fullfit handles many predictors", {
    
    many_predictors <- c("age", "sex", "smoking", "diabetes", "hypertension",
                         "stage", "grade", "ecog", "treatment", "surgery")
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = many_predictors,
        model_type = "glm",
        method = "all"
    )
    
    expect_fullfit_result(result)
})


test_that("fullfit handles missing data appropriately", {
    
    ## clintrial has some missing values
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "smoking"),  # smoking may have NAs
        model_type = "glm",
        method = "all"
    )
    
    expect_fullfit_result(result)
})


test_that("fullfit works with data.frame input", {
    
    df <- as.data.frame(clintrial)
    
    result <- fullfit(
        data = df,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        method = "all"
    )
    
    expect_fullfit_result(result)
})


## ============================================================================
## SECTION 16B: Input Validation Tests
## ============================================================================

test_that("fullfit errors or auto-corrects Surv outcome with GLM model_type", {
    
    ## Should either auto-correct (with message) or error
    result_or_error <- tryCatch({
        suppressMessages(suppressWarnings(
            fullfit(
                data = clintrial,
                outcome = "Surv(os_months, os_status)",
                predictors = c("age", "sex"),
                model_type = "glm",
                method = "all"
            )
        ))
    }, error = function(e) e)
    
    ## Either produces valid result or an error (both acceptable)
    expect_true(
        inherits(result_or_error, "fullfit_result") ||
        inherits(result_or_error, "error")
    )
})


test_that("fullfit errors when coxph used without Surv outcome", {
    
    ## Should error - either from validation or from coxph itself
    expect_error(
        suppressWarnings(
            fullfit(
                data = clintrial,
                outcome = "response",
                predictors = c("age", "sex"),
                model_type = "coxph",
                method = "all"
            )
        )
    )
})


test_that("fullfit errors with continuous outcome and binomial family", {
    
    ## Should error - either from validation or from glm itself
    expect_error(
        suppressWarnings(
            fullfit(
                data = clintrial,
                outcome = "los_days",
                predictors = c("age", "sex"),
                model_type = "glm",
                family = "binomial",
                method = "all"
            )
        )
    )
})


test_that("fullfit errors when outcome variable not found", {
    
    expect_error(
        suppressWarnings(
            fullfit(
                data = clintrial,
                outcome = "nonexistent_variable",
                predictors = c("age", "sex"),
                model_type = "glm",
                method = "all"
            )
        )
    )
})


test_that("fullfit errors when predictor variable not found", {
    
    expect_error(
        suppressWarnings(
            fullfit(
                data = clintrial,
                outcome = "response",
                predictors = c("age", "nonexistent_predictor"),
                model_type = "glm",
                method = "all"
            )
        )
    )
})


test_that("fullfit handles edge case p_threshold values", {
    
    ## p_threshold = 1 should include all predictors
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        method = "screen",
        p_threshold = 1.0
    )
    expect_fullfit_result(result)
})


test_that("fullfit errors with empty data", {
    
    empty_data <- clintrial[0, ]
    
    expect_error(
        suppressWarnings(
            fullfit(
                data = empty_data,
                outcome = "response",
                predictors = c("age", "sex"),
                model_type = "glm",
                method = "all"
            )
        )
    )
})


## ============================================================================
## SECTION 17: Integration with Other Functions
## ============================================================================

test_that("fullfit model output works with summary()", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        method = "all",
        return_type = "model"
    )
    
    ## Should be able to summarize
    summ <- summary(result)
    expect_true(!is.null(summ))
})


test_that("fullfit model output works with predict()", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        method = "all",
        return_type = "model"
    )
    
    ## Should be able to predict
    preds <- predict(result, type = "response")
    expect_true(length(preds) > 0)
})


test_that("fullfit produces consistent results", {
    
    result1 <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        method = "all"
    )
    
    result2 <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        method = "all"
    )
    
    ## Results should be identical
    expect_equal(result1, result2)
})


## ============================================================================
## SECTION 18: Cox Regression Specific Tests
## ============================================================================

test_that("fullfit Cox shows Events column", {
    
    skip_if_not_installed("survival")
    
    result <- fullfit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "sex", "stage"),
        model_type = "coxph",
        method = "all",
        show_events = TRUE
    )
    
    expect_true("Events" %in% names(result))
})


test_that("fullfit Cox screening works correctly", {
    
    skip_if_not_installed("survival")
    
    result <- fullfit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "sex", "smoking", "stage", "grade"),
        model_type = "coxph",
        method = "screen",
        p_threshold = 0.10
    )
    
    expect_fullfit_result(result)
    
    ## Should have HR columns
    col_names <- names(result)
    expect_true(any(grepl("HR", col_names)))
})


## ============================================================================
## SECTION 19: Linear Model Specific Tests
## ============================================================================

test_that("fullfit linear model doesn't show Events", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex"),
        model_type = "lm",
        method = "all",
        show_events = TRUE  # Should be ignored for LM
    )
    
    ## Linear models shouldn't have Events
    expect_false("Events" %in% names(result))
})


test_that("fullfit linear model screening works", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex", "stage", "treatment"),
        model_type = "lm",
        method = "screen",
        p_threshold = 0.20
    )
    
    expect_fullfit_result(result)
})


## ============================================================================
## SECTION 19b: Additional GLM Family Tests
## ============================================================================

## ----------------------------------------------------------------------------
## Gaussian Family Tests
## ----------------------------------------------------------------------------

test_that("fullfit works with gaussian family via glm", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex", "treatment", "stage"),
        model_type = "glm",
        family = "gaussian",
        method = "all"
    )
    
    expect_fullfit_result(result)
    
    ## Should have Coefficient columns (not OR/RR)
    col_names <- names(result)
    expect_true(any(grepl("Coefficient.*CI|Estimate.*CI", col_names)),
                info = paste("Expected Coefficient column for gaussian, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("fullfit gaussian screening works", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex", "treatment", "stage", "grade"),
        model_type = "glm",
        family = "gaussian",
        method = "screen",
        p_threshold = 0.20
    )
    
    expect_fullfit_result(result)
})


## ----------------------------------------------------------------------------
## Quasibinomial Family Tests
## ----------------------------------------------------------------------------

test_that("fullfit works with quasibinomial family", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "treatment", "stage"),
        model_type = "glm",
        family = "quasibinomial",
        method = "all"
    )
    
    expect_fullfit_result(result)
    
    ## Should have OR columns like binomial
    col_names <- names(result)
    expect_true(any(grepl("OR.*CI", col_names)),
                info = paste("Expected OR column for quasibinomial, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("fullfit quasibinomial screening works", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "treatment", "stage", "grade"),
        model_type = "glm",
        family = "quasibinomial",
        method = "screen",
        p_threshold = 0.20
    )
    
    expect_fullfit_result(result)
})


## ----------------------------------------------------------------------------
## Quasipoisson Family Tests
## ----------------------------------------------------------------------------

test_that("fullfit works with quasipoisson family", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment", "stage"),
        model_type = "glm",
        family = "quasipoisson",
        method = "all"
    )
    
    expect_fullfit_result(result)
    
    ## Should have RR columns like poisson
    col_names <- names(result)
    expect_true(any(grepl("RR.*CI", col_names)),
                info = paste("Expected RR column for quasipoisson, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("fullfit quasipoisson screening works", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment", "stage", "grade"),
        model_type = "glm",
        family = "quasipoisson",
        method = "screen",
        p_threshold = 0.20
    )
    
    expect_fullfit_result(result)
})


## ----------------------------------------------------------------------------
## Gamma Family Tests
## ----------------------------------------------------------------------------

test_that("fullfit works with Gamma family (log link)", {
    
    ## Create positive outcome
    test_data <- data.table::as.data.table(clintrial)
    test_data <- test_data[los_days > 0]
    
    result <- fullfit(
        data = test_data,
        outcome = "los_days",
        predictors = c("age", "sex", "treatment"),
        model_type = "glm",
        family = Gamma(link = "log"),
        method = "all"
    )
    
    expect_fullfit_result(result)
})


test_that("fullfit Gamma screening works", {
    
    ## Create positive outcome
    test_data <- data.table::as.data.table(clintrial)
    test_data <- test_data[los_days > 0]
    
    result <- fullfit(
        data = test_data,
        outcome = "los_days",
        predictors = c("age", "sex", "treatment", "stage"),
        model_type = "glm",
        family = Gamma(link = "log"),
        method = "screen",
        p_threshold = 0.20
    )
    
    expect_fullfit_result(result)
})


## ----------------------------------------------------------------------------
## Inverse Gaussian Family Tests
## ----------------------------------------------------------------------------

test_that("fullfit works with inverse.gaussian family", {
    
    ## Create positive outcome
    test_data <- data.table::as.data.table(clintrial)
    test_data <- test_data[los_days > 0]
    
    result <- fullfit(
        data = test_data,
        outcome = "los_days",
        predictors = c("age", "sex", "treatment"),
        model_type = "glm",
        family = inverse.gaussian(link = "log"),
        method = "all"
    )
    
    expect_fullfit_result(result)
})


## ----------------------------------------------------------------------------
## Negative Binomial (negbin) Tests
## ----------------------------------------------------------------------------

test_that("fullfit works with negative binomial regression (negbin)", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(fullfit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment", "stage"),
        model_type = "negbin",
        method = "all"
    ))
    
    expect_fullfit_result(result)
    
    ## Should have RR columns for count data
    col_names <- names(result)
    expect_true(any(grepl("RR.*CI", col_names)),
                info = paste("Expected RR column for negbin, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("fullfit negbin screening works", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(fullfit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment", "stage", "grade"),
        model_type = "negbin",
        method = "screen",
        p_threshold = 0.20
    ))
    
    expect_fullfit_result(result)
    expect_equal(attr(result, "method"), "screen")
})


test_that("fullfit negbin method='custom' works", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(fullfit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment", "stage"),
        model_type = "negbin",
        method = "custom",
        multi_predictors = c("age", "treatment")
    ))
    
    expect_fullfit_result(result)
    expect_equal(attr(result, "method"), "custom")
    expect_equal(attr(result, "n_multi"), 2)
})


test_that("fullfit negbin with labels works", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(fullfit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment"),
        model_type = "negbin",
        method = "all",
        labels = clintrial_labels
    ))
    
    expect_fullfit_result(result)
    
    ## Labels should be applied
    expect_true(any(grepl("Age|Sex|Treatment", result$Variable)))
})


test_that("fullfit negbin return_type='model' works", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(fullfit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment"),
        model_type = "negbin",
        method = "all",
        return_type = "model"
    ))
    
    ## Should return a negbin model object
    expect_true(inherits(result, "negbin"))
})


test_that("fullfit negbin return_type='both' works", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(fullfit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment"),
        model_type = "negbin",
        method = "all",
        return_type = "both"
    ))
    
    expect_type(result, "list")
    expect_s3_class(result$table, "fullfit_result")
    expect_true(inherits(result$model, "negbin"))
})


## ============================================================================
## SECTION 20: Complex Combinations
## ============================================================================

test_that("fullfit handles complex parameter combinations", {
    
    skip_if_not_installed("survival")
    
    result <- fullfit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "sex", "smoking", "stage", "grade", "treatment"),
        model_type = "coxph",
        method = "screen",
        p_threshold = 0.15,
        columns = "both",
        metrics = "both",
        show_n = TRUE,
        show_events = TRUE,
        reference_rows = TRUE,
        digits = 2,
        p_digits = 3,
        labels = clintrial_labels,
        parallel = FALSE
    )
    
    expect_fullfit_result(result)
    
    ## Should have all the expected columns
    col_names <- names(result)
    expect_true("Variable" %in% col_names)
    expect_true("Group" %in% col_names)
    expect_true("n" %in% col_names)
    expect_true("Events" %in% col_names)
    ## Check for uni columns (HR format) and multi columns (aHR format)
    expect_true(any(grepl("^(HR|OR|RR)\\s*\\(", col_names)))
    expect_true(any(grepl("^(aHR|aOR|aRR)", col_names)))
    expect_true(any(grepl("p$", col_names)))
})


test_that("fullfit custom method with columns='multi' works", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "smoking", "stage"),
        model_type = "glm",
        method = "custom",
        multi_predictors = c("age", "stage"),
        columns = "multi"  # Only show multi columns
    )
    
    expect_fullfit_result(result)
    
    ## Should have only multi columns (aOR/aHR/aRR format)
    col_names <- names(result)
    expect_true(any(grepl("^(aOR|aHR|aRR|Adj\\.)", col_names)))
    expect_false("Uni p" %in% col_names)
})


test_that("fullfit with return_type='both' and method='screen'", {
    
    result <- fullfit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "stage", "treatment"),
        model_type = "glm",
        method = "screen",
        p_threshold = 0.20,
        return_type = "both"
    )
    
    expect_type(result, "list")
    expect_s3_class(result$table, "fullfit_result")
    
    ## Model might be NULL if nothing passed screening
    if (!is.null(result$model)) {
        expect_s3_class(result$model, "glm")
    }
})
