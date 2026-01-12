#' Test Suite for multifit(), multiforest(), uniforest(), and uniscreen()
#' 
#' Comprehensive tests covering all model types, output formatting,
#' interactions, mixed effects, and visualization functions.
#' 
#' @details Run with testthat::test_file("tests/testthat/test_multivariate.R")

library(testthat)
library(data.table)
library(summata)

## ============================================================================
## Setup: Load test data and define helper functions
## ============================================================================

## Use package's built-in clinical trial data
data(clintrial)
data(clintrial_labels)

## Create a complete-case subset for tests requiring no missing data
clintrial_complete <- na.omit(clintrial)

## Create integer count outcome for Poisson/quasipoisson/negative binomial tests
set.seed(123)
clintrial$count_outcome <- as.integer(ceiling(clintrial$los_days))
clintrial_complete$count_outcome <- as.integer(ceiling(clintrial_complete$los_days))

## Define outcome sets for different model types
binary_outcomes <- c("surgery", "pfs_status", "os_status")
continuous_outcomes <- c("los_days", "biomarker_x")
surv_outcomes <- c("Surv(pfs_months, pfs_status)", "Surv(os_months, os_status)")

## Standard predictor sets
standard_predictors <- c("age", "sex", "treatment", "stage")
surgery_covariates <- c("age", "sex", "ecog", "smoking_status")


## ----------------------------------------------------------------------------
## Helper: Check uniscreen result structure
## ----------------------------------------------------------------------------
expect_uniscreen_result <- function(result) {
    expect_s3_class(result, "uniscreen_result")
    expect_s3_class(result, "data.table")
    
    ## Required columns
    expect_true("Variable" %in% names(result))
    expect_true("p-value" %in% names(result))
    
    ## Required attributes
    expect_true(!is.null(attr(result, "outcome")))
    expect_true(!is.null(attr(result, "model_type")))
    expect_true(!is.null(attr(result, "model_scope")))
    expect_true(!is.null(attr(result, "raw_data")))
}


## ----------------------------------------------------------------------------
## Helper: Check multifit result structure
## ----------------------------------------------------------------------------
expect_multifit_result <- function(result) {
    expect_s3_class(result, "multifit_result")
    expect_s3_class(result, "data.table")
    
    ## Required columns
    expect_true("Outcome" %in% names(result))
    
    ## Required attributes
    expect_true(!is.null(attr(result, "predictor")))
    expect_true(!is.null(attr(result, "outcomes")))
    expect_true(!is.null(attr(result, "model_type")))
    expect_true(!is.null(attr(result, "columns")))
    expect_true(!is.null(attr(result, "raw_data")))
}


## ----------------------------------------------------------------------------
## Helper: Check effect column exists with correct type
## ----------------------------------------------------------------------------
expect_effect_column <- function(result, effect_type = NULL) {
    col_names <- names(result)
    
    ## Find effect column (contains OR, HR, RR, Coefficient, or Estimate with CI)
    effect_cols <- grep("(OR|HR|RR|Coefficient|Estimate).*CI", col_names, value = TRUE)
    expect_true(length(effect_cols) >= 1, 
                info = paste("Expected effect column, found:", paste(col_names, collapse = ", ")))
    
    if (!is.null(effect_type)) {
        ## For linear models, both "Coefficient" and "Estimate" are valid
        if (effect_type == "Coefficient") {
            expect_true(any(grepl("Coefficient|Estimate", effect_cols)),
                        info = paste("Expected Coefficient or Estimate in columns, found:", 
                                     paste(effect_cols, collapse = ", ")))
        } else {
            expect_true(any(grepl(effect_type, effect_cols)),
                        info = paste("Expected", effect_type, "in columns, found:", 
                                     paste(effect_cols, collapse = ", ")))
        }
    }
}


## ----------------------------------------------------------------------------
## Helper: Check p-value column formatting
## ----------------------------------------------------------------------------
expect_pvalue_column <- function(result) {
    expect_true("p-value" %in% names(result))
    
    ## P-values should be formatted as numeric string, < threshold, or -
    pvals <- result$`p-value`
    valid_pvals <- grepl("^[0-9]\\.[0-9]+$|^< 0\\.0*1$|^-$|^$", pvals)
    expect_true(all(valid_pvals), 
                info = paste("Invalid p-values found:", 
                             paste(pvals[!valid_pvals], collapse = ", ")))
}


## ----------------------------------------------------------------------------
## Helper: Check raw_data attribute structure for uniscreen
## ----------------------------------------------------------------------------
expect_uniscreen_raw_data <- function(result) {
    raw <- attr(result, "raw_data")
    expect_s3_class(raw, "data.table")
    
    ## Should have effect estimate columns
    has_effect <- any(c("OR", "HR", "RR", "Estimate", "Coefficient") %in% names(raw))
    expect_true(has_effect, info = "raw_data should have effect estimate column")
    
    ## Should have CI columns (lowercase)
    expect_true("ci_lower" %in% names(raw))
    expect_true("ci_upper" %in% names(raw))
}


## ----------------------------------------------------------------------------
## Helper: Check raw_data attribute structure for multifit
## ----------------------------------------------------------------------------
expect_multifit_raw_data <- function(result, columns = "adjusted") {
    raw <- attr(result, "raw_data")
    expect_s3_class(raw, "data.table")
    
    ## Required columns depend on columns parameter
    if (columns == "both") {
        expect_true("exp_coef_unadj" %in% names(raw) || "exp_coef_adj" %in% names(raw))
    } else {
        expect_true("exp_coef" %in% names(raw) || 
                    "exp_coef_unadj" %in% names(raw) || 
                    "exp_coef_adj" %in% names(raw))
    }
}


## ============================================================================
## SECTION 1: uniscreen() Basic Functionality Tests
## ============================================================================

test_that("uniscreen works with basic logistic regression", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = standard_predictors,
        model_type = "glm",
        family = "binomial"
    )
    
    expect_uniscreen_result(result)
    expect_effect_column(result, "OR")
    expect_pvalue_column(result)
    expect_uniscreen_raw_data(result)
    
    ## All predictors should be represented
    expect_true(any(grepl("age|Age", result$Variable)))
    expect_true(any(grepl("sex|Sex", result$Variable)))
})


test_that("uniscreen works with linear regression", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "los_days",
        predictors = standard_predictors,
        model_type = "lm"
    )
    
    expect_uniscreen_result(result)
    expect_effect_column(result, "Coefficient")
    expect_equal(attr(result, "model_scope"), "Univariable")
})


test_that("uniscreen works with Cox regression", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = standard_predictors,
        model_type = "coxph"
    )
    
    expect_uniscreen_result(result)
    expect_effect_column(result, "HR")
    
    ## model_type attribute uses display name
    expect_equal(attr(result, "model_type"), "Cox PH")
})


test_that("uniscreen p_threshold filtering works", {
    
    ## All results
    result_all <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex", "treatment", "stage", "ecog"),
        model_type = "glm",
        p_threshold = 1
    )
    
    ## Filtered to p < 0.20
    result_filtered <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex", "treatment", "stage", "ecog"),
        model_type = "glm",
        p_threshold = 0.20
    )
    
    expect_lte(nrow(result_filtered), nrow(result_all))
})


test_that("uniscreen applies labels correctly", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex", "treatment"),
        model_type = "glm",
        labels = clintrial_labels
    )
    
    ## Labels should be applied to Variable column
    expect_true(any(grepl("Age|Sex|Treatment", result$Variable)))
})


test_that("uniscreen keep_models stores model objects", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex"),
        model_type = "glm",
        keep_models = TRUE
    )
    
    models <- attr(result, "models")
    expect_type(models, "list")
    expect_true(length(models) >= 2)
    expect_s3_class(models[[1]], "glm")
})


test_that("uniscreen exponentiate parameter works", {
    
    ## Force exponentiate = TRUE
    result_exp <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = "age",
        model_type = "glm",
        exponentiate = TRUE
    )
    
    ## Force exponentiate = FALSE
    result_coef <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = "age",
        model_type = "glm",
        exponentiate = FALSE
    )
    
    ## Column names should differ
    expect_true(any(grepl("OR", names(result_exp))))
    expect_true(any(grepl("Coefficient", names(result_coef))))
})


test_that("uniscreen show_n and show_events parameters work", {
    
    result_full <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "treatment"),
        model_type = "glm",
        show_n = TRUE,
        show_events = TRUE
    )
    
    expect_true("n" %in% names(result_full) || any(grepl("^n$", names(result_full))))
    
    result_minimal <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "treatment"),
        model_type = "glm",
        show_n = FALSE,
        show_events = FALSE
    )
    
    expect_false("n" %in% names(result_minimal))
})


test_that("uniscreen handles single predictor", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = "age",
        model_type = "glm"
    )
    
    expect_uniscreen_result(result)
    expect_true(nrow(result) >= 1)
})


test_that("uniscreen respects digits parameter", {
    
    result_2dig <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = "age",
        model_type = "glm",
        digits = 2
    )
    
    result_4dig <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = "age",
        model_type = "glm",
        digits = 4
    )
    
    ## 4-digit result should have more characters in effect column
    effect_col_2 <- grep("OR.*CI", names(result_2dig), value = TRUE)[1]
    effect_col_4 <- grep("OR.*CI", names(result_4dig), value = TRUE)[1]
    
    expect_true(nchar(result_4dig[[effect_col_4]][1]) > nchar(result_2dig[[effect_col_2]][1]))
})


## ============================================================================
## SECTION 1A: uniscreen() Additional GLM Family Tests
## ============================================================================

## ----------------------------------------------------------------------------
## Gaussian Family Tests
## ----------------------------------------------------------------------------

test_that("uniscreen works with gaussian family", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "los_days",
        predictors = standard_predictors,
        model_type = "glm",
        family = "gaussian"
    )
    
    expect_uniscreen_result(result)
    
    ## Should produce coefficients
    col_names <- names(result)
    expect_true(any(grepl("Coefficient.*CI|Estimate.*CI", col_names)),
                info = paste("Expected Coefficient column for gaussian, found:", 
                             paste(col_names, collapse = ", ")))
})


## ----------------------------------------------------------------------------
## Quasibinomial Family Tests
## ----------------------------------------------------------------------------

test_that("uniscreen works with quasibinomial family", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = standard_predictors,
        model_type = "glm",
        family = "quasibinomial"
    )
    
    expect_uniscreen_result(result)
    
    ## Should produce odds ratios
    col_names <- names(result)
    expect_true(any(grepl("OR.*CI", col_names)),
                info = paste("Expected OR column for quasibinomial, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("uniscreen quasibinomial point estimates match binomial", {
    
    result_quasi <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = "age",
        model_type = "glm",
        family = "quasibinomial",
        keep_models = TRUE
    )
    
    result_binom <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = "age",
        model_type = "glm",
        family = "binomial",
        keep_models = TRUE
    )
    
    ## Point estimates should match
    raw_quasi <- attr(result_quasi, "raw_data")
    raw_binom <- attr(result_binom, "raw_data")
    
    expect_equal(raw_quasi$OR[1], raw_binom$OR[1], tolerance = 1e-6)
})


## ----------------------------------------------------------------------------
## Quasipoisson Family Tests
## ----------------------------------------------------------------------------

test_that("uniscreen works with quasipoisson family", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "count_outcome",
        predictors = standard_predictors,
        model_type = "glm",
        family = "quasipoisson"
    )
    
    expect_uniscreen_result(result)
    
    ## Should produce rate ratios
    col_names <- names(result)
    expect_true(any(grepl("RR.*CI", col_names)),
                info = paste("Expected RR column for quasipoisson, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("uniscreen quasipoisson point estimates match poisson", {
    
    result_quasi <- uniscreen(
        data = clintrial,
        outcome = "count_outcome",
        predictors = "age",
        model_type = "glm",
        family = "quasipoisson",
        keep_models = TRUE
    )
    
    result_pois <- uniscreen(
        data = clintrial,
        outcome = "count_outcome",
        predictors = "age",
        model_type = "glm",
        family = "poisson",
        keep_models = TRUE
    )
    
    ## Point estimates should match
    raw_quasi <- attr(result_quasi, "raw_data")
    raw_pois <- attr(result_pois, "raw_data")
    
    expect_equal(raw_quasi$RR[1], raw_pois$RR[1], tolerance = 1e-6)
})


## ----------------------------------------------------------------------------
## Gamma Family Tests
## ----------------------------------------------------------------------------

test_that("uniscreen works with Gamma family", {
    
    ## Create positive outcome
    test_data <- data.table::as.data.table(clintrial)
    test_data <- test_data[los_days > 0]
    
    result <- uniscreen(
        data = test_data,
        outcome = "los_days",
        predictors = c("age", "sex", "treatment"),
        model_type = "glm",
        family = Gamma(link = "log")
    )
    
    expect_uniscreen_result(result)
})


## ----------------------------------------------------------------------------
## Inverse Gaussian Family Tests
## ----------------------------------------------------------------------------

test_that("uniscreen works with inverse.gaussian family", {
    
    ## Create positive outcome
    test_data <- data.table::as.data.table(clintrial)
    test_data <- test_data[los_days > 0]
    
    result <- uniscreen(
        data = test_data,
        outcome = "los_days",
        predictors = c("age", "sex", "treatment"),
        model_type = "glm",
        family = inverse.gaussian(link = "log")
    )
    
    expect_uniscreen_result(result)
})


## ----------------------------------------------------------------------------
## Negative Binomial (negbin) Tests
## ----------------------------------------------------------------------------

test_that("uniscreen works with negative binomial regression (negbin)", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(uniscreen(
        data = clintrial,
        outcome = "count_outcome",
        predictors = standard_predictors,
        model_type = "negbin",
        parallel = FALSE
    ))
    
    expect_uniscreen_result(result)
    
    ## Should produce rate ratios
    col_names <- names(result)
    expect_true(any(grepl("RR.*CI", col_names)),
                info = paste("Expected RR column for negbin, found:", 
                             paste(col_names, collapse = ", ")))
    
    ## All predictors should be represented
    expect_true(any(grepl("age|Age", result$Variable)))
    expect_true(any(grepl("sex|Sex", result$Variable)))
})


test_that("uniscreen negbin p_threshold filtering works", {
    
    skip_if_not_installed("MASS")
    
    ## All results
    result_all <- suppressWarnings(uniscreen(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment", "stage", "ecog"),
        model_type = "negbin",
        p_threshold = 1,
        parallel = FALSE
    ))
    
    ## Filtered to p < 0.20
    result_filtered <- suppressWarnings(uniscreen(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment", "stage", "ecog"),
        model_type = "negbin",
        p_threshold = 0.20,
        parallel = FALSE
    ))
    
    expect_lte(nrow(result_filtered), nrow(result_all))
})


test_that("uniscreen negbin with labels works", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(uniscreen(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment"),
        model_type = "negbin",
        labels = clintrial_labels,
        parallel = FALSE
    ))
    
    expect_uniscreen_result(result)
    expect_true(any(grepl("Age|Sex|Treatment", result$Variable)))
})


test_that("uniscreen negbin keep_models stores negbin objects", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(uniscreen(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex"),
        model_type = "negbin",
        keep_models = TRUE,
        parallel = FALSE
    ))
    
    models <- attr(result, "models")
    expect_type(models, "list")
    expect_true(length(models) >= 2)
    expect_true(inherits(models[[1]], "negbin"))
})


test_that("uniscreen negbin coefficients match direct model", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(uniscreen(
        data = clintrial_complete,
        outcome = "count_outcome",
        predictors = "age",
        model_type = "negbin",
        keep_models = TRUE,
        parallel = FALSE
    ))
    
    ## Fit directly
    direct_model <- suppressWarnings(MASS::glm.nb(count_outcome ~ age, data = clintrial_complete))
    
    ## Get raw data
    raw <- attr(result, "raw_data")
    
    ## RR should match exp(coefficient)
    expected_rr <- exp(coef(direct_model)["age"])
    expect_equal(raw$RR[1], unname(expected_rr), tolerance = 1e-6)
})


## ============================================================================
## SECTION 1B: uniscreen() Input Validation Tests
## ============================================================================

test_that("uniscreen errors or auto-corrects Surv outcome with GLM model_type", {
    
    ## Should either auto-correct (with message) or error
    result_or_error <- tryCatch({
        suppressMessages(suppressWarnings(
            uniscreen(
                data = clintrial,
                outcome = "Surv(os_months, os_status)",
                predictors = c("age", "sex"),
                model_type = "glm"
            )
        ))
    }, error = function(e) e)
    
    ## Either produces valid result or an error (both acceptable)
    expect_true(
        inherits(result_or_error, "uniscreen_result") ||
        inherits(result_or_error, "error")
    )
})


test_that("uniscreen errors when coxph used without Surv outcome", {
    
    ## Should error - either from validation or from coxph itself
    expect_error(
        suppressWarnings(
            uniscreen(
                data = clintrial,
                outcome = "os_status",
                predictors = c("age", "sex"),
                model_type = "coxph"
            )
        )
    )
})


test_that("uniscreen errors with continuous outcome and binomial family", {
    
    ## Should error - either from validation or from glm itself
    expect_error(
        suppressWarnings(
            uniscreen(
                data = clintrial,
                outcome = "los_days",
                predictors = c("age", "sex"),
                model_type = "glm",
                family = "binomial"
            )
        )
    )
})


test_that("uniscreen errors when outcome variable not found", {
    
    expect_error(
        suppressWarnings(
            uniscreen(
                data = clintrial,
                outcome = "nonexistent_variable",
                predictors = c("age", "sex"),
                model_type = "glm"
            )
        )
    )
})


test_that("uniscreen errors when predictor variable not found", {
    
    expect_error(
        suppressWarnings(
            uniscreen(
                data = clintrial,
                outcome = "os_status",
                predictors = c("age", "nonexistent_predictor"),
                model_type = "glm"
            )
        )
    )
})


test_that("uniscreen handles edge case p_threshold values", {
    
    ## p_threshold = 1 should include all predictors
    result <- uniscreen(
        data = clintrial,
        outcome = "os_status",
        predictors = c("age", "sex"),
        model_type = "glm",
        p_threshold = 1.0
    )
    expect_uniscreen_result(result)
})


test_that("uniscreen errors with empty data", {
    
    empty_data <- clintrial[0, ]
    
    expect_error(
        suppressWarnings(
            uniscreen(
                data = empty_data,
                outcome = "os_status",
                predictors = c("age", "sex"),
                model_type = "glm"
            )
        )
    )
})


test_that("uniscreen handles binary outcome with lm model_type", {
    
    ## Linear probability model - should work (may warn)
    result <- suppressWarnings(
        uniscreen(
            data = clintrial,
            outcome = "os_status",
            predictors = c("age", "sex"),
            model_type = "lm"
        )
    )
    
    expect_uniscreen_result(result)
})


## ============================================================================
## SECTION 2: multifit() Basic Functionality Tests
## ============================================================================

test_that("multifit works with basic logistic regression - unadjusted", {
    
    result <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    expect_equal(attr(result, "columns"), "unadjusted")
    expect_effect_column(result, "OR")
    expect_pvalue_column(result)
    
    ## Should have rows for each outcome
    expect_equal(length(unique(result$Outcome)), 3)
})


test_that("multifit works with basic logistic regression - adjusted", {
    
    result <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        covariates = c("age", "sex", "stage"),
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    expect_equal(attr(result, "columns"), "adjusted")
    expect_effect_column(result, "aOR")
    
    ## Should store covariates attribute
    expect_equal(attr(result, "covariates"), c("age", "sex", "stage"))
})


test_that("multifit works with columns = 'both'", {
    
    result <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        covariates = c("age", "sex"),
        columns = "both",
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    expect_equal(attr(result, "columns"), "both")
    
    ## Should have both p-value columns
    col_names <- names(result)
    expect_true("Uni p" %in% col_names)
    expect_true("Multi p" %in% col_names)
})


test_that("multifit handles continuous predictor", {
    
    result <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "age",
        covariates = c("sex", "stage"),
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    
    ## One row per outcome for continuous predictor
    expect_equal(nrow(result), 3)
})


test_that("multifit handles factor predictor with multiple levels", {
    
    result <- multifit(
        data = clintrial,
        outcomes = "surgery",
        predictor = "treatment",
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    
    ## Should have Predictor column with levels
    expect_true("Predictor" %in% names(result))
})


test_that("multifit show_n and show_events parameters work", {
    
    result_full <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        show_n = TRUE,
        show_events = TRUE,
        parallel = FALSE
    )
    
    expect_true("n" %in% names(result_full))
    expect_true("Events" %in% names(result_full))
    
    result_no_n <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        show_n = FALSE,
        show_events = FALSE,
        parallel = FALSE
    )
    
    expect_false("n" %in% names(result_no_n))
    expect_false("Events" %in% names(result_no_n))
})


test_that("multifit respects digits and p_digits parameters", {
    
    result <- multifit(
        data = clintrial,
        outcomes = "surgery",
        predictor = "treatment",
        digits = 3,
        p_digits = 4,
        parallel = FALSE
    )
    
    expect_multifit_result(result)
})


## ============================================================================
## SECTION 3: Cox Proportional Hazards Tests
## ============================================================================

test_that("uniscreen works with Cox model", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = standard_predictors,
        model_type = "coxph"
    )
    
    expect_uniscreen_result(result)
    expect_effect_column(result, "HR")
})


test_that("multifit works with Cox model", {
    
    result <- multifit(
        data = clintrial,
        outcomes = surv_outcomes,
        predictor = "treatment",
        covariates = c("age", "sex"),
        model_type = "coxph",
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    expect_effect_column(result, "HR")
})


test_that("multifit Cox with strata parameter", {
    
    result <- multifit(
        data = clintrial,
        outcomes = surv_outcomes[1],
        predictor = "treatment",
        covariates = c("age", "sex"),
        strata = "site",
        model_type = "coxph",
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    expect_equal(attr(result, "strata"), "site")
})


test_that("multifit Cox with cluster parameter", {
    
    result <- multifit(
        data = clintrial,
        outcomes = surv_outcomes[1],
        predictor = "treatment",
        covariates = c("age", "sex"),
        cluster = "site",
        model_type = "coxph",
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    expect_equal(attr(result, "cluster"), "site")
})


test_that("strata parameter validation for non-Cox models", {
    
    expect_error(
        multifit(clintrial, binary_outcomes, "treatment", 
                 strata = "site", model_type = "glm", parallel = FALSE),
        "'strata' is only supported for coxph"
    )
})


test_that("cluster parameter validation for non-Cox models", {
    
    expect_error(
        multifit(clintrial, binary_outcomes, "treatment", 
                 cluster = "site", model_type = "glm", parallel = FALSE),
        "'cluster' is only supported for coxph"
    )
})


## ============================================================================
## SECTION 4: Linear Model Tests
## ============================================================================

test_that("uniscreen works with linear model", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "los_days",
        predictors = standard_predictors,
        model_type = "lm"
    )
    
    expect_uniscreen_result(result)
    expect_effect_column(result, "Coefficient")
})


test_that("multifit works with linear model", {
    
    result <- multifit(
        data = clintrial,
        outcomes = "los_days",
        predictor = "treatment",
        covariates = c("age", "sex", "stage"),
        model_type = "lm",
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    expect_effect_column(result, "Coefficient")
    
    ## Events column should not be present for LM
    expect_false("Events" %in% names(result))
})


## ============================================================================
## SECTION 5: Poisson Regression Tests
## ============================================================================

test_that("multifit works with Poisson regression", {
    
    ## Create count variable
    clintrial$los_count <- ceiling(clintrial$los_days / 5)
    
    result <- multifit(
        data = clintrial,
        outcomes = "los_count",
        predictor = "treatment",
        covariates = c("age", "sex"),
        family = "poisson",
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    expect_effect_column(result, "RR")
})


## ----------------------------------------------------------------------------
## Gaussian Family Tests
## ----------------------------------------------------------------------------

test_that("multifit works with gaussian family", {
    
    result <- multifit(
        data = clintrial,
        outcomes = "los_days",
        predictor = "treatment",
        covariates = c("age", "sex"),
        family = "gaussian",
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    
    ## Should produce coefficients
    col_names <- names(result)
    expect_true(any(grepl("Coefficient.*CI|Estimate.*CI", col_names)),
                info = paste("Expected Coefficient column for gaussian, found:", 
                             paste(col_names, collapse = ", ")))
})


## ----------------------------------------------------------------------------
## Quasibinomial Family Tests
## ----------------------------------------------------------------------------

test_that("multifit works with quasibinomial family", {
    
    result <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        covariates = c("age", "sex"),
        family = "quasibinomial",
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    
    ## Should produce odds ratios
    col_names <- names(result)
    expect_true(any(grepl("OR.*CI", col_names)),
                info = paste("Expected OR column for quasibinomial, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("multifit quasibinomial with columns='both' works", {
    
    result <- multifit(
        data = clintrial,
        outcomes = "surgery",
        predictor = "treatment",
        covariates = c("age", "sex"),
        family = "quasibinomial",
        columns = "both",
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    
    ## Should have both OR and aOR columns
    col_names <- names(result)
    expect_true(any(grepl("^OR", col_names)) && any(grepl("^aOR", col_names)))
})


## ----------------------------------------------------------------------------
## Quasipoisson Family Tests
## ----------------------------------------------------------------------------

test_that("multifit works with quasipoisson family", {
    
    result <- multifit(
        data = clintrial,
        outcomes = "count_outcome",
        predictor = "treatment",
        covariates = c("age", "sex"),
        family = "quasipoisson",
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    
    ## Should produce rate ratios
    col_names <- names(result)
    expect_true(any(grepl("RR.*CI", col_names)),
                info = paste("Expected RR column for quasipoisson, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("multifit quasipoisson with columns='both' works", {
    
    result <- multifit(
        data = clintrial,
        outcomes = "count_outcome",
        predictor = "treatment",
        covariates = c("age", "sex"),
        family = "quasipoisson",
        columns = "both",
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    
    ## Should have both RR and aRR columns
    col_names <- names(result)
    expect_true(any(grepl("^RR", col_names)) && any(grepl("^aRR", col_names)))
})


## ----------------------------------------------------------------------------
## Gamma Family Tests
## ----------------------------------------------------------------------------

test_that("multifit works with Gamma family", {
    
    ## Create positive outcome
    test_data <- data.table::as.data.table(clintrial)
    test_data <- test_data[los_days > 0]
    
    result <- multifit(
        data = test_data,
        outcomes = "los_days",
        predictor = "treatment",
        covariates = c("age", "sex"),
        family = Gamma(link = "log"),
        parallel = FALSE
    )
    
    expect_multifit_result(result)
})


## ----------------------------------------------------------------------------
## Inverse Gaussian Family Tests
## ----------------------------------------------------------------------------

test_that("multifit works with inverse.gaussian family", {
    
    ## Create positive outcome
    test_data <- data.table::as.data.table(clintrial)
    test_data <- test_data[los_days > 0]
    
    result <- multifit(
        data = test_data,
        outcomes = "los_days",
        predictor = "treatment",
        covariates = c("age", "sex"),
        family = inverse.gaussian(link = "log"),
        parallel = FALSE
    )
    
    expect_multifit_result(result)
})


## ----------------------------------------------------------------------------
## Negative Binomial (negbin) Tests
## ----------------------------------------------------------------------------

test_that("multifit works with negative binomial regression (negbin)", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(multifit(
        data = clintrial,
        outcomes = "count_outcome",
        predictor = "treatment",
        covariates = c("age", "sex"),
        model_type = "negbin",
        parallel = FALSE
    ))
    
    expect_multifit_result(result)
    
    ## Should produce rate ratios
    col_names <- names(result)
    expect_true(any(grepl("RR.*CI", col_names)),
                info = paste("Expected RR column for negbin, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("multifit negbin works with multiple outcomes", {
    
    skip_if_not_installed("MASS")
    
    ## Create second count variable
    clintrial$count_outcome2 <- as.integer(ceiling(clintrial$los_days / 5))
    
    result <- suppressWarnings(multifit(
        data = clintrial,
        outcomes = c("count_outcome", "count_outcome2"),
        predictor = "treatment",
        covariates = c("age", "sex"),
        model_type = "negbin",
        parallel = FALSE
    ))
    
    expect_multifit_result(result)
    
    ## Should have rows for both outcomes
    expect_equal(length(unique(result$Outcome)), 2)
})


test_that("multifit negbin works with columns = 'both'", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(multifit(
        data = clintrial,
        outcomes = "count_outcome",
        predictor = "treatment",
        covariates = c("age", "sex"),
        model_type = "negbin",
        columns = "both",
        parallel = FALSE
    ))
    
    expect_multifit_result(result)
    
    ## Should have unadjusted (RR) and adjusted (aRR) columns
    col_names <- names(result)
    expect_true(any(grepl("^RR", col_names)) && any(grepl("^aRR", col_names)),
                info = paste("Expected both RR and aRR columns, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("multifit negbin with interactions works", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(multifit(
        data = clintrial,
        outcomes = "count_outcome",
        predictor = "treatment",
        covariates = c("age", "sex", "stage"),
        interactions = c("treatment:stage"),
        model_type = "negbin",
        parallel = FALSE
    ))
    
    expect_multifit_result(result)
    expect_equal(attr(result, "interactions"), "treatment:stage")
})


test_that("multifit negbin with labels works", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(multifit(
        data = clintrial,
        outcomes = "count_outcome",
        predictor = "treatment",
        covariates = c("age", "sex"),
        model_type = "negbin",
        labels = clintrial_labels,
        parallel = FALSE
    ))
    
    expect_multifit_result(result)
    
    ## Labels should be applied
    expect_true(any(grepl("Treatment", result$Predictor)))
})


test_that("multifit negbin keep_models stores negbin objects", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(multifit(
        data = clintrial,
        outcomes = "count_outcome",
        predictor = "treatment",
        covariates = c("age", "sex"),
        model_type = "negbin",
        keep_models = TRUE,
        parallel = FALSE
    ))
    
    expect_multifit_result(result)
    
    models <- attr(result, "models")
    expect_type(models, "list")
    
    ## Should have negbin model
    adj_model <- models$count_outcome$adjusted
    expect_true(inherits(adj_model, "negbin"))
})


## ============================================================================
## SECTION 6: Interaction Terms Tests
## ============================================================================

test_that("multifit works with factor * factor interaction", {
    
    result <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        covariates = c("age", "sex", "stage"),
        interactions = c("treatment:stage"),
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    expect_equal(attr(result, "interactions"), "treatment:stage")
    
    ## Interaction terms should contain * symbol
    expect_true(any(grepl("[:×*]", result$Predictor)))
})


test_that("multifit works with continuous * factor interaction", {
    
    result <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "age",
        covariates = c("sex", "stage"),
        interactions = c("age:sex"),
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    expect_true(any(grepl("[:×*]", result$Predictor)))
})


test_that("multifit works with multiple interactions", {
    
    result <- multifit(
        data = clintrial,
        outcomes = "surgery",
        predictor = "treatment",
        covariates = c("age", "sex", "stage"),
        interactions = c("treatment:sex", "treatment:stage"),
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    expect_equal(attr(result, "interactions"), c("treatment:sex", "treatment:stage"))
})


test_that("Cox model with interaction works", {
    
    result <- multifit(
        data = clintrial,
        outcomes = surv_outcomes[1],
        predictor = "treatment",
        covariates = c("age", "sex", "stage"),
        interactions = c("treatment:stage"),
        model_type = "coxph",
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    expect_true(any(grepl("[:×*]", result$Predictor)))
})


## ============================================================================
## SECTION 7: Mixed Effects Models Tests
## ============================================================================

test_that("multifit works with glmer", {
    skip_if_not_installed("lme4")
    
    result <- suppressWarnings(multifit(
        data = clintrial,
        outcomes = binary_outcomes[1],
        predictor = "treatment",
        random = "(1|site)",
        model_type = "glmer",
        parallel = FALSE
    ))
    
    expect_multifit_result(result)
    expect_equal(attr(result, "random"), "(1|site)")
})


test_that("multifit works with lmer", {
    skip_if_not_installed("lme4")
    
    result <- suppressWarnings(multifit(
        data = clintrial,
        outcomes = "los_days",
        predictor = "treatment",
        random = "(1|site)",
        model_type = "lmer",
        parallel = FALSE
    ))
    
    expect_multifit_result(result)
})


test_that("multifit works with nested random effects", {
    skip_if_not_installed("lme4")
    
    ## Create nested structure
    clintrial$patient_id <- 1:nrow(clintrial)
    
    result <- suppressWarnings(multifit(
        data = clintrial,
        outcomes = "surgery",
        predictor = "treatment",
        random = "(1|site)",
        model_type = "glmer",
        parallel = FALSE
    ))
    
    expect_multifit_result(result)
})


test_that("multifit works with coxme", {
    skip_if_not_installed("coxme")
    
    result <- suppressWarnings(multifit(
        data = clintrial,
        outcomes = surv_outcomes[1],
        predictor = "treatment",
        random = "(1|site)",
        model_type = "coxme",
        parallel = FALSE
    ))
    
    expect_multifit_result(result)
})


test_that("glmer with interaction and random effects works", {
    skip_if_not_installed("lme4")
    
    result <- suppressWarnings(multifit(
        data = clintrial,
        outcomes = "surgery",
        predictor = "treatment",
        covariates = c("age", "sex"),
        interactions = c("treatment:sex"),
        random = "(1|site)",
        model_type = "glmer",
        parallel = FALSE
    ))
    
    expect_multifit_result(result)
})


## ============================================================================
## SECTION 8: Labels Tests
## ============================================================================

test_that("uniscreen applies labels to Variable column", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex", "treatment"),
        model_type = "glm",
        labels = clintrial_labels
    )
    
    ## Labels should appear in Variable column
    has_labels <- any(sapply(clintrial_labels, function(l) any(grepl(l, result$Variable))))
    expect_true(has_labels)
})


test_that("multifit applies labels to Outcome column", {
    
    result <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        labels = c(
            surgery = "Surgical Resection",
            pfs_status = "Progression Event",
            os_status = "Death"
        ),
        parallel = FALSE
    )
    
    expect_true("Surgical Resection" %in% result$Outcome)
    expect_true("Progression Event" %in% result$Outcome)
    expect_true("Death" %in% result$Outcome)
})


test_that("multifit predictor_label parameter works", {
    
    result <- multifit(
        data = clintrial,
        outcomes = "surgery",
        predictor = "age",
        predictor_label = "Patient Age (years)",
        parallel = FALSE
    )
    
    ## predictor_label affects formatting, not necessarily stored as attribute
    ## Just verify the function runs without error
    expect_multifit_result(result)
})


## ============================================================================
## SECTION 9: P-value Filtering Tests
## ============================================================================

test_that("multifit p_threshold filtering works", {
    
    result_all <- multifit(
        data = clintrial,
        outcomes = c("surgery", "pfs_status", "os_status"),
        predictor = "age",
        covariates = c("sex"),
        p_threshold = 1,
        parallel = FALSE
    )
    
    result_filtered <- multifit(
        data = clintrial,
        outcomes = c("surgery", "pfs_status", "os_status"),
        predictor = "age",
        covariates = c("sex"),
        p_threshold = 0.05,
        parallel = FALSE
    )
    
    expect_lte(nrow(result_filtered), nrow(result_all))
})


## ============================================================================
## SECTION 10: Keep Models Tests
## ============================================================================

test_that("uniscreen keep_models stores all models", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex", "treatment"),
        model_type = "glm",
        keep_models = TRUE
    )
    
    models <- attr(result, "models")
    expect_type(models, "list")
    expect_true(length(models) >= 3)
})


test_that("multifit keep_models stores model objects", {
    
    result <- multifit(
        data = clintrial,
        outcomes = binary_outcomes[1:2],
        predictor = "treatment",
        covariates = c("age", "sex"),
        keep_models = TRUE,
        parallel = FALSE
    )
    
    models <- attr(result, "models")
    expect_type(models, "list")
    expect_equal(length(models), 2)
})


test_that("stored models can be used for predictions", {
    
    result <- multifit(
        data = clintrial_complete,
        outcomes = "surgery",
        predictor = "treatment",
        covariates = c("age", "sex"),
        keep_models = TRUE,
        parallel = FALSE
    )
    
    models <- attr(result, "models")
    model <- models$surgery$adjusted
    
    preds <- predict(model, type = "response")
    expect_true(length(preds) > 0)
    expect_true(all(preds >= 0 & preds <= 1))
})


## ============================================================================
## SECTION 11: Raw Data Attribute Tests
## ============================================================================

test_that("uniscreen raw_data contains expected columns", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "treatment"),
        model_type = "glm"
    )
    
    expect_uniscreen_raw_data(result)
    
    raw <- attr(result, "raw_data")
    expect_true("p_value" %in% names(raw))
})


test_that("multifit raw_data contains expected columns", {
    
    result <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "age",
        parallel = FALSE
    )
    
    expect_multifit_raw_data(result)
    
    raw <- attr(result, "raw_data")
    expect_true("outcome" %in% names(raw))
    expect_true("predictor" %in% names(raw))
})


test_that("multifit raw_data has correct structure for 'both' columns", {
    
    result <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        covariates = c("age", "sex"),
        columns = "both",
        parallel = FALSE
    )
    
    raw <- attr(result, "raw_data")
    
    ## Should have both adjusted and unadjusted columns
    expect_true("exp_coef_unadj" %in% names(raw))
    expect_true("exp_coef_adj" %in% names(raw))
})


## ============================================================================
## SECTION 12: uniforest() Tests
## ============================================================================

test_that("uniforest creates forest plot from uniscreen result", {
    skip_if_not_installed("ggplot2")
    
    table <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex", "treatment"),
        model_type = "glm"
    )
    
    p <- uniforest(table, title = "Univariable Screening")
    
    expect_s3_class(p, "ggplot")
    expect_true(!is.null(attr(p, "recommended_dims")))
})


test_that("uniforest works with Cox model results", {
    skip_if_not_installed("ggplot2")
    
    table <- uniscreen(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "sex", "treatment", "stage"),
        model_type = "coxph"
    )
    
    p <- uniforest(table, title = "Survival Analysis")
    
    expect_s3_class(p, "ggplot")
})


test_that("uniforest works with linear model results", {
    skip_if_not_installed("ggplot2")
    
    table <- uniscreen(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex", "treatment"),
        model_type = "lm"
    )
    
    p <- uniforest(table, title = "Linear Model")
    
    expect_s3_class(p, "ggplot")
})


test_that("uniforest show_n and show_events parameters work", {
    skip_if_not_installed("ggplot2")
    
    table <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex"),
        model_type = "glm"
    )
    
    p1 <- uniforest(table, show_n = TRUE, show_events = TRUE)
    p2 <- uniforest(table, show_n = FALSE, show_events = FALSE)
    
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
})


test_that("uniforest indent_groups parameter works", {
    skip_if_not_installed("ggplot2")
    
    table <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("treatment", "stage"),
        model_type = "glm"
    )
    
    p1 <- uniforest(table, indent_groups = FALSE)
    p2 <- uniforest(table, indent_groups = TRUE)
    
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
})


test_that("uniforest bold_variables parameter works", {
    skip_if_not_installed("ggplot2")
    
    table <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex"),
        model_type = "glm"
    )
    
    p1 <- uniforest(table, bold_variables = TRUE)
    p2 <- uniforest(table, bold_variables = FALSE)
    
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
})


test_that("uniforest zebra_stripes parameter works", {
    skip_if_not_installed("ggplot2")
    
    table <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex", "treatment"),
        model_type = "glm"
    )
    
    p1 <- uniforest(table, zebra_stripes = TRUE)
    p2 <- uniforest(table, zebra_stripes = FALSE)
    
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
})


test_that("uniforest custom color works", {
    skip_if_not_installed("ggplot2")
    
    table <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex"),
        model_type = "glm"
    )
    
    p <- uniforest(table, color = "#FF5733")
    
    expect_s3_class(p, "ggplot")
})


test_that("uniforest labels parameter works", {
    skip_if_not_installed("ggplot2")
    
    table <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex", "treatment"),
        model_type = "glm"
    )
    
    p <- uniforest(table, labels = clintrial_labels)
    
    expect_s3_class(p, "ggplot")
})


## ============================================================================
## SECTION 13: multiforest() Tests
## ============================================================================

test_that("multiforest creates forest plot from multifit result", {
    skip_if_not_installed("ggplot2")
    
    table <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        parallel = FALSE
    )
    
    p <- multiforest(table, title = "Multi-Outcome Analysis")
    
    expect_s3_class(p, "ggplot")
    expect_true(!is.null(attr(p, "recommended_dims")))
})


test_that("multiforest works with adjusted results", {
    skip_if_not_installed("ggplot2")
    
    table <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        covariates = c("age", "sex"),
        parallel = FALSE
    )
    
    p <- multiforest(table, column = "adjusted", covariates_footer = TRUE)
    
    expect_s3_class(p, "ggplot")
})


test_that("multiforest works with survival outcomes", {
    skip_if_not_installed("ggplot2")
    
    table <- multifit(
        data = clintrial,
        outcomes = surv_outcomes,
        predictor = "treatment",
        model_type = "coxph",
        parallel = FALSE
    )
    
    p <- multiforest(table)
    
    expect_s3_class(p, "ggplot")
})


test_that("multiforest works with linear model outcomes", {
    skip_if_not_installed("ggplot2")
    
    table <- multifit(
        data = clintrial,
        outcomes = c("los_days", "biomarker_x"),
        predictor = "treatment",
        model_type = "lm",
        parallel = FALSE
    )
    
    p <- multiforest(table)
    
    expect_s3_class(p, "ggplot")
})


test_that("multiforest bold_variables parameter works", {
    skip_if_not_installed("ggplot2")
    
    table <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        parallel = FALSE
    )
    
    p1 <- multiforest(table, bold_variables = TRUE)
    p2 <- multiforest(table, bold_variables = FALSE)
    
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
})


test_that("multiforest indent_predictor parameter works", {
    skip_if_not_installed("ggplot2")
    
    table <- multifit(
        data = clintrial,
        outcomes = binary_outcomes[1:2],
        predictor = "treatment",
        parallel = FALSE
    )
    
    p1 <- multiforest(table, indent_predictor = FALSE)
    p2 <- multiforest(table, indent_predictor = TRUE)
    
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
})


test_that("multiforest zebra_stripes parameter works", {
    skip_if_not_installed("ggplot2")
    
    table <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        parallel = FALSE
    )
    
    p1 <- multiforest(table, zebra_stripes = TRUE)
    p2 <- multiforest(table, zebra_stripes = FALSE)
    
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
})


test_that("multiforest custom color works", {
    skip_if_not_installed("ggplot2")
    
    table <- multifit(
        data = clintrial,
        outcomes = binary_outcomes[1:2],
        predictor = "treatment",
        parallel = FALSE
    )
    
    p <- multiforest(table, color = "#E74C3C")
    
    expect_s3_class(p, "ggplot")
})


test_that("multiforest with labels works", {
    skip_if_not_installed("ggplot2")
    
    table <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        labels = clintrial_labels,
        parallel = FALSE
    )
    
    p <- multiforest(table, labels = clintrial_labels)
    
    expect_s3_class(p, "ggplot")
})


test_that("multiforest with mixed-effects model results works", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("ggplot2")
    
    table <- suppressWarnings(multifit(
        data = clintrial,
        outcomes = binary_outcomes[1:2],
        predictor = "treatment",
        random = "(1|site)",
        model_type = "glmer",
        parallel = FALSE
    ))
    
    p <- multiforest(table)
    
    expect_s3_class(p, "ggplot")
})


## ============================================================================
## SECTION 14: Input Validation Tests
## ============================================================================

test_that("uniscreen validates required inputs", {
    
    expect_error(
        uniscreen(predictors = c("age"), model_type = "glm"),
        "data.*required|missing"
    )
})


test_that("uniscreen errors on non-existent predictor", {
    
    expect_error(
        uniscreen(
            data = clintrial,
            outcome = "surgery",
            predictors = c("nonexistent_var"),
            model_type = "glm"
        )
    )
})


test_that("multifit validates required inputs", {
    
    expect_error(
        multifit(outcomes = binary_outcomes, predictor = "treatment", parallel = FALSE),
        "'data' is required"
    )
    
    expect_error(
        multifit(data = clintrial, predictor = "treatment", parallel = FALSE),
        "'outcomes' is required"
    )
    
    expect_error(
        multifit(data = clintrial, outcomes = binary_outcomes, parallel = FALSE),
        "'predictor' is required"
    )
})


test_that("multifit validates predictor length", {
    
    expect_error(
        multifit(clintrial, binary_outcomes, c("age", "sex"), parallel = FALSE),
        "'predictor' must be a single variable name"
    )
})


test_that("multifit warns on missing random for glmer", {
    skip_if_not_installed("lme4")
    
    ## multifit uses tryCatch and warns instead of erroring
    expect_warning(
        multifit(
            data = clintrial,
            outcomes = c("surgery"),
            predictor = "treatment",
            model_type = "glmer",
            parallel = FALSE
        ),
        regexp = "random.*required|failed"
    )
})


test_that("multifit warns on non-existent outcome variable", {
    
    ## multifit uses tryCatch and warns instead of erroring
    expect_warning(
        multifit(
            data = clintrial,
            outcomes = c("nonexistent_var"),
            predictor = "treatment",
            model_type = "glm",
            parallel = FALSE
        ),
        regexp = "failed|not found"
    )
})


test_that("uniforest errors on non-uniscreen input", {
    skip_if_not_installed("ggplot2")
    
    expect_error(uniforest(data.frame(x = 1:5)))
})


test_that("multiforest errors on non-multifit input", {
    skip_if_not_installed("ggplot2")
    
    expect_error(multiforest(data.frame(x = 1:5)))
})


## ============================================================================
## SECTION 15: Edge Cases
## ============================================================================

test_that("uniscreen handles single predictor", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = "age",
        model_type = "glm"
    )
    
    expect_uniscreen_result(result)
})


test_that("multifit handles single outcome", {
    
    result <- multifit(
        data = clintrial,
        outcomes = "surgery",
        predictor = "treatment",
        parallel = FALSE
    )
    
    expect_multifit_result(result)
    expect_equal(length(unique(result$Outcome)), 1)
})


test_that("uniscreen handles missing data appropriately", {
    
    dt_missing <- data.table::as.data.table(clintrial)
    set.seed(123)
    dt_missing[sample(1:.N, 20), age := NA]
    
    expect_no_error({
        result <- uniscreen(
            data = dt_missing,
            outcome = "surgery",
            predictors = c("age", "sex"),
            model_type = "glm"
        )
    })
})


## ============================================================================
## SECTION 16: Integration Tests
## ============================================================================

test_that("full pipeline: uniscreen -> uniforest", {
    skip_if_not_installed("ggplot2")
    
    ## Screen predictors
    screen_results <- uniscreen(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "sex", "treatment", "stage", "ecog"),
        model_type = "coxph",
        labels = clintrial_labels
    )
    
    ## Create forest plot
    p <- uniforest(
        screen_results,
        title = "Univariable Screening: Overall Survival",
        show_n = TRUE,
        show_events = TRUE,
        zebra_stripes = TRUE,
        bold_variables = TRUE
    )
    
    expect_uniscreen_result(screen_results)
    expect_s3_class(p, "ggplot")
})


test_that("full pipeline: multifit -> multiforest with covariates", {
    skip_if_not_installed("ggplot2")
    
    ## Multi-outcome analysis
    multi_results <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        covariates = c("age", "sex"),
        labels = clintrial_labels,
        parallel = FALSE
    )
    
    ## Create forest plot
    p <- multiforest(
        multi_results,
        title = "Treatment Effects Across Outcomes",
        covariates_footer = TRUE,
        zebra_stripes = TRUE,
        bold_variables = TRUE
    )
    
    expect_multifit_result(multi_results)
    expect_s3_class(p, "ggplot")
})


test_that("full pipeline: multifit with interactions -> multiforest", {
    skip_if_not_installed("ggplot2")
    
    ## Multi-outcome analysis with interaction
    multi_results <- multifit(
        data = clintrial,
        outcomes = binary_outcomes[1:2],
        predictor = "treatment",
        covariates = c("age", "sex"),
        interactions = c("treatment:sex"),
        labels = clintrial_labels,
        parallel = FALSE
    )
    
    ## Create forest plot
    p <- multiforest(
        multi_results,
        title = "Treatment * Sex Interaction",
        zebra_stripes = TRUE
    )
    
    expect_multifit_result(multi_results)
    expect_s3_class(p, "ggplot")
})


test_that("full pipeline: multifit with mixed effects -> multiforest", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("ggplot2")
    
    ## Multi-outcome analysis with random effects
    multi_results <- suppressWarnings(multifit(
        data = clintrial,
        outcomes = binary_outcomes[1:2],
        predictor = "treatment",
        covariates = c("age", "sex"),
        random = "(1|site)",
        model_type = "glmer",
        labels = clintrial_labels,
        parallel = FALSE
    ))
    
    ## Create forest plot
    p <- multiforest(
        multi_results,
        title = "Treatment Effects with Random Effects",
        covariates_footer = TRUE,
        zebra_stripes = TRUE,
        bold_variables = TRUE
    )
    
    expect_multifit_result(multi_results)
    expect_s3_class(p, "ggplot")
})


## ============================================================================
## SECTION 17: Parallel Processing Tests
## ============================================================================

## test_that("parallel processing produces same results as sequential", {
##     skip_on_cran()
##     
##     set.seed(123)
##     result_seq <- multifit(
##         data = clintrial,
##         outcomes = binary_outcomes,
##         predictor = "treatment",
##         covariates = c("age", "sex"),
##         parallel = FALSE
##     )
##     
##     set.seed(123)
##     result_par <- multifit(
##         data = clintrial,
##         outcomes = binary_outcomes,
##         predictor = "treatment",
##         covariates = c("age", "sex"),
##         parallel = TRUE,
##         n_cores = 2
##     )
##     
##     ## Results should be identical (may differ in row order)
##     expect_equal(nrow(result_seq), nrow(result_par))
##     expect_equal(sort(result_seq$Outcome), sort(result_par$Outcome))
## })


## ============================================================================
## SECTION 18: Attribute Preservation Tests
## ============================================================================

test_that("uniscreen preserves all expected attributes", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex"),
        model_type = "glm"
    )
    
    expect_true(!is.null(attr(result, "raw_data")))
    expect_true(!is.null(attr(result, "outcome")))
    expect_true(!is.null(attr(result, "model_type")))
    expect_true(!is.null(attr(result, "model_scope")))
})


test_that("multifit preserves all expected attributes", {
    
    result <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        covariates = c("age", "sex"),
        interactions = c("treatment:sex"),
        parallel = FALSE
    )
    
    expect_equal(attr(result, "predictor"), "treatment")
    expect_equal(attr(result, "outcomes"), binary_outcomes)
    expect_equal(attr(result, "covariates"), c("age", "sex"))
    expect_equal(attr(result, "interactions"), c("treatment:sex"))
    expect_equal(attr(result, "model_type"), "glm")
})


test_that("uniforest recommended_dims attribute present", {
    skip_if_not_installed("ggplot2")
    
    table <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex"),
        model_type = "glm"
    )
    
    p <- uniforest(table)
    
    dims <- attr(p, "recommended_dims")
    expect_true(!is.null(dims))
    expect_true("width" %in% names(dims))
    expect_true("height" %in% names(dims))
})


test_that("multiforest recommended_dims attribute present", {
    skip_if_not_installed("ggplot2")
    
    table <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        parallel = FALSE
    )
    
    p <- multiforest(table)
    
    dims <- attr(p, "recommended_dims")
    expect_true(!is.null(dims))
    expect_true("width" %in% names(dims))
    expect_true("height" %in% names(dims))
})


## ============================================================================
## SECTION 19: Print Method Tests
## ============================================================================

test_that("uniscreen print method displays correctly", {
    
    result <- uniscreen(
        data = clintrial,
        outcome = "surgery",
        predictors = c("age", "sex", "treatment"),
        model_type = "glm"
    )
    
    output <- capture.output(print(result))
    
    expect_true(any(grepl("Univariable", output)))
})


test_that("multifit print method displays correctly", {
    
    result <- multifit(
        data = clintrial,
        outcomes = binary_outcomes,
        predictor = "treatment",
        covariates = c("age", "sex"),
        parallel = FALSE
    )
    
    output <- capture.output(print(result))
    
    expect_true(any(grepl("Multivariate|Multi", output)))
    expect_true(any(grepl("Predictor", output)))
})


test_that("multifit print method shows interactions when present", {
    
    result <- multifit(
        data = clintrial,
        outcomes = "surgery",
        predictor = "treatment",
        covariates = c("stage"),
        interactions = c("treatment:stage"),
        parallel = FALSE
    )
    
    output <- capture.output(print(result))
    expect_true(any(grepl("Interaction", output)))
})


## ============================================================================
## SECTION 20: Coefficient Accuracy Tests
## ============================================================================

test_that("multifit coefficients match direct model fitting", {
    
    ## Fit using multifit
    result <- multifit(
        data = clintrial_complete,
        outcomes = "surgery",
        predictor = "treatment",
        covariates = c("age", "sex"),
        keep_models = TRUE,
        parallel = FALSE
    )
    
    ## Fit directly
    direct_model <- glm(surgery ~ treatment + age + sex, 
                        data = clintrial_complete, family = binomial)
    
    ## Get multifit model coefficients
    mf_model <- attr(result, "models")$surgery$adjusted
    mf_coefs <- coef(mf_model)
    direct_coefs <- coef(direct_model)
    
    ## Coefficients should match
    expect_equal(mf_coefs, direct_coefs, tolerance = 1e-10)
})


test_that("uniscreen coefficients match direct model fitting", {
    
    ## Fit using uniscreen
    result <- uniscreen(
        data = clintrial_complete,
        outcome = "surgery",
        predictors = "age",
        model_type = "glm",
        keep_models = TRUE
    )
    
    ## Fit directly
    direct_model <- glm(surgery ~ age, data = clintrial_complete, family = binomial)
    
    ## Get raw data
    raw <- attr(result, "raw_data")
    
    ## OR should match exp(coefficient)
    expected_or <- exp(coef(direct_model)["age"])
    expect_equal(raw$OR[1], unname(expected_or), tolerance = 1e-6)
})
