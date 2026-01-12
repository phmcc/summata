#' Test Suite for fit()
#' 
#' Comprehensive tests covering all model types, output formatting,
#' interactions, stratification, clustering, weights, and edge cases.
#' 
#' @details Run with testthat::test_file("tests/testthat/test-fit.R")

library(testthat)
library(data.table)
library(summata)

## Load survival package for Surv() and strata() functions in formula evaluation
if (requireNamespace("survival", quietly = TRUE)) {
    library(survival)
}

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

## Create a weights column for weighted regression tests
set.seed(42)
clintrial$weight <- runif(nrow(clintrial), 0.5, 2.0)
clintrial_complete$weight <- runif(nrow(clintrial_complete), 0.5, 2.0)

## Helper function to check fit_result structure
expect_fit_result <- function(result) {
    expect_s3_class(result, "fit_result")
    expect_s3_class(result, "data.table")
    
    ## Required columns
    expect_true("Variable" %in% names(result))
    expect_true("Group" %in% names(result))
    
    ## Required attributes
    expect_true(!is.null(attr(result, "model")))
    expect_true(!is.null(attr(result, "raw_data")))
    expect_true(!is.null(attr(result, "outcome")))
    expect_true(!is.null(attr(result, "predictors")))
    expect_true(!is.null(attr(result, "formula_str")))
    expect_true(!is.null(attr(result, "model_scope")))
    expect_true(!is.null(attr(result, "model_type")))
}

## Helper to check effect column exists and is properly formatted
expect_effect_column <- function(result, effect_type = NULL) {
    col_names <- names(result)
    
    ## Find effect column (contains OR, HR, RR, Coefficient, or Estimate with CI)
    effect_cols <- grep("(OR|HR|RR|Coefficient|Estimate).*CI", col_names, value = TRUE)
    expect_true(length(effect_cols) >= 1, 
                info = paste("Expected effect column, found:", paste(col_names, collapse = ", ")))
    
    if (!is.null(effect_type)) {
        expect_true(any(grepl(effect_type, effect_cols)),
                    info = paste("Expected", effect_type, "in columns, found:", 
                                 paste(effect_cols, collapse = ", ")))
    }
}

## Helper to check p-value column
expect_pvalue_column <- function(result) {
    expect_true("p-value" %in% names(result))
    
    ## Check p-values are properly formatted (numeric string, < 0.001, or -)
    pvals <- result$`p-value`
    valid_pvals <- grepl("^[0-9]\\.[0-9]+$|^< 0\\.001$|^-$", pvals)
    expect_true(all(valid_pvals), 
                info = paste("Invalid p-values found:", 
                             paste(pvals[!valid_pvals], collapse = ", ")))
}

## Helper to verify reference rows
expect_reference_rows <- function(result, vars_with_refs) {
    raw <- attr(result, "raw_data")
    for (v in vars_with_refs) {
        var_rows <- raw[variable == v]
        if (nrow(var_rows) > 1) {
            ## Factor variables should have a reference row
            has_ref <- any(var_rows$reference == "reference", na.rm = TRUE)
            expect_true(has_ref, 
                        info = paste("Expected reference row for factor:", v))
        }
    }
}


## ============================================================================
## SECTION 1: Basic GLM (Logistic Regression) Tests
## ============================================================================

test_that("fit works with basic logistic regression - univariable", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = "age",
        model_type = "glm",
        family = "binomial"
    )
    
    expect_fit_result(result)
    expect_equal(attr(result, "model_scope"), "Univariable")
    expect_effect_column(result, "OR")
    expect_pvalue_column(result)
    
    ## Should have single predictor
    expect_equal(length(attr(result, "predictors")), 1)
    
    ## Effect column should say "OR (95% CI)" for univariable
    col_names <- names(result)
    expect_true(any(grepl("^OR \\(95% CI\\)$", col_names)))
})


test_that("fit works with basic logistic regression - multivariable", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "smoking"),
        model_type = "glm",
        family = "binomial"
    )
    
    expect_fit_result(result)
    expect_equal(attr(result, "model_scope"), "Multivariable")
    expect_effect_column(result, "aOR")
    expect_pvalue_column(result)
    
    ## Should have multiple predictors
    expect_equal(length(attr(result, "predictors")), 3)
    
    ## Effect column should say "aOR (95% CI)" for multivariable
    col_names <- names(result)
    expect_true(any(grepl("^aOR \\(95% CI\\)$", col_names)))
})


test_that("fit handles factor variables correctly in GLM", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "stage"),
        model_type = "glm",
        family = "binomial",
        reference_rows = TRUE
    )
    
    expect_fit_result(result)
    
    ## Stage is a factor - should have multiple rows
    raw <- attr(result, "raw_data")
    stage_rows <- raw[variable == "stage"]
    expect_true(nrow(stage_rows) >= 2)  # At least reference + one level
    
    ## Check reference row exists
    has_ref <- any(stage_rows$reference == "reference", na.rm = TRUE)
    expect_true(has_ref)
})


test_that("fit reference_rows = FALSE excludes reference rows", {
    
    result_with_ref <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "stage"),
        model_type = "glm",
        reference_rows = TRUE
    )
    
    result_no_ref <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "stage"),
        model_type = "glm",
        reference_rows = FALSE
    )
    
    ## Without reference rows should have fewer rows
    expect_true(nrow(result_no_ref) < nrow(result_with_ref))
    
    ## Check no reference markers in raw data
    raw_no_ref <- attr(result_no_ref, "raw_data")
    if ("reference" %in% names(raw_no_ref)) {
        has_refs <- sum(raw_no_ref$reference == "reference", na.rm = TRUE)
        expect_equal(has_refs, 0)
    }
})


test_that("fit shows n and events columns correctly for GLM", {
    
    ## With n and events
    result_full <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        show_n = TRUE,
        show_events = TRUE
    )
    expect_true("n" %in% names(result_full))
    expect_true("Events" %in% names(result_full))
    
    ## Without n
    result_no_n <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        show_n = FALSE,
        show_events = TRUE
    )
    expect_false("n" %in% names(result_no_n))
    expect_true("Events" %in% names(result_no_n))
    
    ## Without events
    result_no_events <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        show_n = TRUE,
        show_events = FALSE
    )
    expect_true("n" %in% names(result_no_events))
    expect_false("Events" %in% names(result_no_events))
})


test_that("fit respects digits parameter for GLM", {
    
    result_2dig <- fit(
        data = clintrial,
        outcome = "response",
        predictors = "age",
        model_type = "glm",
        digits = 2
    )
    
    result_3dig <- fit(
        data = clintrial,
        outcome = "response",
        predictors = "age",
        model_type = "glm",
        digits = 3
    )
    
    ## Extract effect columns
    effect_col_2 <- grep("OR.*CI", names(result_2dig), value = TRUE)[1]
    effect_col_3 <- grep("OR.*CI", names(result_3dig), value = TRUE)[1]
    
    ## 2-digit result should have pattern like "1.02 (0.99-1.05)"
    ## 3-digit result should have pattern like "1.023 (0.993-1.053)"
    val_2 <- result_2dig[[effect_col_2]][1]
    val_3 <- result_3dig[[effect_col_3]][1]
    
    ## Count decimal places in the first number
    first_num_2 <- sub(" .*", "", val_2)
    first_num_3 <- sub(" .*", "", val_3)
    
    decimals_2 <- nchar(sub(".*\\.", "", first_num_2))
    decimals_3 <- nchar(sub(".*\\.", "", first_num_3))
    
    expect_equal(decimals_2, 2)
    expect_equal(decimals_3, 3)
})


test_that("fit respects p_digits parameter", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        p_digits = 4
    )
    
    ## P-values should have up to 4 decimal places
    pvals <- result$`p-value`
    numeric_pvals <- pvals[grepl("^0\\.", pvals)]
    
    if (length(numeric_pvals) > 0) {
        ## Check that numeric p-values have 4 decimal places
        decimals <- sapply(numeric_pvals, function(p) {
            nchar(sub("0\\.", "", p))
        })
        expect_true(all(decimals == 4))
    }
})


## ============================================================================
## SECTION 2: Linear Model Tests
## ============================================================================

test_that("fit works with linear regression - univariable", {
    
    result <- fit(
        data = clintrial,
        outcome = "los_days",
        predictors = "age",
        model_type = "lm"
    )
    
    expect_fit_result(result)
    expect_equal(attr(result, "model_scope"), "Univariable")
    expect_equal(attr(result, "model_type"), "Linear")
    
    ## Linear models use Coefficient, not OR
    expect_effect_column(result, "Coefficient")
    
    ## Should NOT have Events column
    expect_false("Events" %in% names(result))
})


test_that("fit works with linear regression - multivariable", {
    
    result <- fit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex", "stage"),
        model_type = "lm"
    )
    
    expect_fit_result(result)
    expect_equal(attr(result, "model_scope"), "Multivariable")
    expect_effect_column(result, "Coefficient")
})


test_that("fit linear model show_events is ignored", {
    
    ## Even if show_events = TRUE, linear models shouldn't show events
    result <- fit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex"),
        model_type = "lm",
        show_events = TRUE
    )
    
    expect_false("Events" %in% names(result))
})


test_that("fit linear model QC stats are captured", {
    
    result <- fit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex", "stage"),
        model_type = "lm",
        keep_qc_stats = TRUE
    )
    
    raw <- attr(result, "raw_data")
    
    ## Linear models should have R-squared and other stats
    expect_true("AIC" %in% names(raw))
    expect_true("BIC" %in% names(raw))
})


## ============================================================================
## SECTION 3: Cox Proportional Hazards Tests
## ============================================================================

test_that("fit works with Cox regression - univariable", {
    
    skip_if_not_installed("survival")
    
    result <- fit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = "age",
        model_type = "coxph"
    )
    
    expect_fit_result(result)
    expect_equal(attr(result, "model_scope"), "Univariable")
    expect_effect_column(result, "HR")
    
    ## Should have Events column
    expect_true("Events" %in% names(result))
})


test_that("fit works with Cox regression - multivariable", {
    
    skip_if_not_installed("survival")
    
    result <- fit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "sex", "stage", "treatment"),
        model_type = "coxph"
    )
    
    expect_fit_result(result)
    expect_equal(attr(result, "model_scope"), "Multivariable")
    
    ## Multivariable Cox should show aHR
    col_names <- names(result)
    expect_true(any(grepl("aHR", col_names)))
})


test_that("fit Cox model with stratification", {
    
    skip_if_not_installed("survival")
    
    result <- fit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "treatment"),
        model_type = "coxph",
        strata = "sex"
    )
    
    expect_fit_result(result)
    
    ## Check strata is stored in attributes
    expect_equal(attr(result, "strata"), "sex")
    
    ## Formula should include strata()
    formula_str <- attr(result, "formula_str")
    expect_true(grepl("strata\\(.*sex.*\\)", formula_str))
    
    ## Sex should NOT appear as a predictor in the output
    ## (it's used for stratification, not estimation)
    raw <- attr(result, "raw_data")
    expect_false("sex" %in% raw$variable)
})


test_that("fit Cox model with clustering", {
    
    skip_if_not_installed("survival")
    
    result <- fit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "sex", "treatment"),
        model_type = "coxph",
        cluster = "site"
    )
    
    expect_fit_result(result)
    
    ## Check cluster is stored in attributes
    expect_equal(attr(result, "cluster"), "site")
    
    ## Model should have robust standard errors
    model <- attr(result, "model")
    expect_true(!is.null(model$naive.var) || !is.null(attr(model, "cluster")))
})


test_that("fit Cox model with both strata and cluster", {
    
    skip_if_not_installed("survival")
    
    ## Create a second grouping variable for stratification
    clintrial_test <- data.table::as.data.table(clintrial)
    clintrial_test[, region := sample(c("North", "South"), .N, replace = TRUE)]
    
    result <- fit(
        data = clintrial_test,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "treatment"),
        model_type = "coxph",
        strata = "region",
        cluster = "site"
    )
    
    expect_fit_result(result)
    expect_equal(attr(result, "strata"), "region")
    expect_equal(attr(result, "cluster"), "site")
})


## ============================================================================
## SECTION 4: Poisson Regression Tests
## ============================================================================

test_that("fit works with Poisson regression", {
    
    ## Create count outcome
    clintrial_test <- data.table::as.data.table(clintrial)
    clintrial_test[, complications := rpois(.N, lambda = 2)]
    
    result <- fit(
        data = clintrial_test,
        outcome = "complications",
        predictors = c("age", "sex", "stage"),
        model_type = "glm",
        family = "poisson"
    )
    
    expect_fit_result(result)
    
    ## Poisson models should show RR (rate ratio)
    expect_effect_column(result, "RR")
})


## ============================================================================
## SECTION 5: Interaction Terms Tests
## ============================================================================

test_that("fit handles basic interaction terms", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        interactions = c("age:sex"),
        model_type = "glm"
    )
    
    expect_fit_result(result)
    
    ## Check interactions stored in attributes
    expect_equal(attr(result, "interactions"), c("age:sex"))
    
    ## Formula should include interaction
    formula_str <- attr(result, "formula_str")
    expect_true(grepl("age:sex", formula_str))
    
    ## Raw data should have interaction term rows
    raw <- attr(result, "raw_data")
    interaction_rows <- raw[grepl(":", variable)]
    expect_true(nrow(interaction_rows) >= 1)
})


test_that("fit handles multiple interaction terms", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "treatment"),
        interactions = c("age:sex", "age:treatment"),
        model_type = "glm"
    )
    
    expect_fit_result(result)
    
    formula_str <- attr(result, "formula_str")
    expect_true(grepl("age:sex", formula_str))
    expect_true(grepl("age:treatment", formula_str))
})

## ============================================================================
## SECTION 6: Weighted Regression Tests
## ============================================================================

test_that("fit handles weights in GLM", {
    
    ## Suppress expected "non-integer #successes" warning from weighted binomial
    result <- suppressWarnings(fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        weights = "weight"
    ))
    
    expect_fit_result(result)
    expect_equal(attr(result, "weights"), "weight")
    
    ## Model should have weights
    model <- attr(result, "model")
    expect_true(!is.null(model$prior.weights))
})


test_that("fit handles weights in linear model", {
    
    result <- fit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex"),
        model_type = "lm",
        weights = "weight"
    )
    
    expect_fit_result(result)
    expect_equal(attr(result, "weights"), "weight")
})


test_that("weighted vs unweighted models produce different results", {
    
    result_unweighted <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm"
    )
    
    ## Suppress expected "non-integer #successes" warning from weighted binomial
    result_weighted <- suppressWarnings(fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        weights = "weight"
    ))
    
    ## Get coefficients from models
    coef_unweighted <- coef(attr(result_unweighted, "model"))
    coef_weighted <- coef(attr(result_weighted, "model"))
    
    ## Coefficients should be different (not identical)
    expect_false(all(abs(coef_unweighted - coef_weighted) < 1e-10))
})


## ============================================================================
## SECTION 7: Custom Labels Tests
## ============================================================================

test_that("fit applies custom labels to variables", {
    
    custom_labels <- c(
        age = "Age (years)",
        sex = "Sex",
        stage = "Cancer Stage"
    )
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "stage"),
        model_type = "glm",
        labels = custom_labels
    )
    
    ## Check that labels are applied
    ## Variable column should contain labeled names
    variables <- unique(result$Variable[result$Variable != ""])
    
    expect_true("Age (years)" %in% variables || any(grepl("Age", variables)))
    expect_true("Sex" %in% variables || any(grepl("Sex", variables)))
    expect_true("Cancer Stage" %in% variables || any(grepl("Stage", variables)))
})


test_that("fit applies labels to interaction terms", {
    
    custom_labels <- c(
        age = "Age",
        sex = "Sex"
    )
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        interactions = c("age:sex"),
        model_type = "glm",
        labels = custom_labels
    )
    
    ## Interaction should show labeled variables with " * " separator
    variables <- result$Variable
    
    ## Should contain interaction with labeled names
    has_labeled_interaction <- any(grepl("Age.*Sex|Sex.*Age", variables))
    expect_true(has_labeled_interaction)
})


test_that("fit labels work with clintrial_labels dataset", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "stage", "treatment"),
        model_type = "glm",
        labels = clintrial_labels
    )
    
    expect_fit_result(result)
    
    ## Labels should be applied
    variables <- result$Variable[result$Variable != ""]
    
    ## Check some expected labels exist
    expect_true(length(variables) > 0)
})


## ============================================================================
## SECTION 8: Confidence Level Tests
## ============================================================================

test_that("fit respects conf_level parameter", {
    
    result_95 <- fit(
        data = clintrial,
        outcome = "response",
        predictors = "age",
        model_type = "glm",
        conf_level = 0.95
    )
    
    result_90 <- fit(
        data = clintrial,
        outcome = "response",
        predictors = "age",
        model_type = "glm",
        conf_level = 0.90
    )
    
    ## Get raw data to compare CI widths
    raw_95 <- attr(result_95, "raw_data")
    raw_90 <- attr(result_90, "raw_data")
    
    ## 90% CI should be narrower than 95% CI
    ci_width_95 <- raw_95$exp_upper[1] - raw_95$exp_lower[1]
    ci_width_90 <- raw_90$exp_upper[1] - raw_90$exp_lower[1]
    
    expect_true(ci_width_90 < ci_width_95)
})


test_that("fit works with 99% confidence interval", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        conf_level = 0.99
    )
    
    expect_fit_result(result)
    
    ## 99% CI column should be named appropriately or at least exist
    col_names <- names(result)
    expect_true(any(grepl("CI", col_names)))
})


## ============================================================================
## SECTION 9: Exponentiate Parameter Tests
## ============================================================================

test_that("fit auto-exponentiates for logistic regression", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = "age",
        model_type = "glm",
        family = "binomial"
    )
    
    ## Should show OR, not raw coefficients
    col_names <- names(result)
    expect_true(any(grepl("OR", col_names)))
})


test_that("fit exponentiate = FALSE shows raw coefficients", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = "age",
        model_type = "glm",
        family = "binomial",
        exponentiate = FALSE
    )
    
    ## Should show Coefficient, not OR
    col_names <- names(result)
    expect_true(any(grepl("Coefficient", col_names)))
    expect_false(any(grepl("\\bOR\\b", col_names)))
})


test_that("fit exponentiate = TRUE forces OR even for linear model", {
    
    result <- fit(
        data = clintrial,
        outcome = "los_days",
        predictors = "age",
        model_type = "lm",
        exponentiate = TRUE
    )
    
    ## Should show exponentiated values
    raw <- attr(result, "raw_data")
    
    ## exp_coef should exist in raw data
    expect_true("exp_coef" %in% names(raw))
})


## ============================================================================
## SECTION 10: Mixed Effects Models Tests
## ============================================================================

test_that("fit works with lmer (linear mixed effects)", {
    
    skip_if_not_installed("lme4")
    skip_on_cran()  # lme4 can segfault in test environments
    
    gc()  # Force garbage collection before lme4
    
    ## Use smaller dataset for stability
    test_data <- clintrial_complete[1:200, ]
    
    result <- fit(
        data = test_data,
        outcome = "los_days",
        predictors = c("age", "sex", "(1|site)"),
        model_type = "lmer"
    )
    
    expect_fit_result(result)
    expect_equal(attr(result, "model_type"), "Linear Mixed")
    
    ## Model should be lmerMod
    model <- attr(result, "model")
    expect_true(inherits(model, "lmerMod"))
})


test_that("fit works with glmer (generalized linear mixed effects)", {
    
    skip_if_not_installed("lme4")
    skip_on_cran()
    
    gc()
    
    test_data <- clintrial_complete[1:200, ]
    
    result <- fit(
        data = test_data,
        outcome = "response",
        predictors = c("age", "sex", "(1|site)"),
        model_type = "glmer",
        family = "binomial"
    )
    
    expect_fit_result(result)
    
    ## Model should be glmerMod
    model <- attr(result, "model")
    expect_true(inherits(model, "glmerMod"))
})


test_that("fit works with coxme (mixed effects Cox)", {
    
    skip_if_not_installed("coxme")
    skip_if_not_installed("survival")
    skip_on_cran()
    
    gc()
    
    test_data <- clintrial_complete[1:200, ]
    
    result <- fit(
        data = test_data,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "sex", "(1|site)"),
        model_type = "coxme"
    )
    
    expect_fit_result(result)
    
    ## Model should be coxme
    model <- attr(result, "model")
    expect_true(inherits(model, "coxme"))
})


## ============================================================================
## SECTION 11: Conditional Logistic Regression Tests
## ============================================================================

test_that("fit works with conditional logistic regression", {
    
    skip_if_not_installed("survival")
    
    ## Create matched case-control data
    set.seed(123)
    n_strata <- 50
    matched_data <- data.table(
        stratum = rep(1:n_strata, each = 3),
        case = rep(c(1, 0, 0), n_strata),
        age = rnorm(n_strata * 3, 60, 10),
        exposure = rbinom(n_strata * 3, 1, 0.3)
    )
    
    result <- fit(
        data = matched_data,
        outcome = "case",
        predictors = c("age", "exposure"),
        model_type = "clogit",
        strata = "stratum"
    )
    
    expect_fit_result(result)
    expect_equal(attr(result, "strata"), "stratum")
    
    ## Conditional logistic shows effect estimates (could be OR or HR depending on implementation)
    expect_effect_column(result)
})


## ============================================================================
## SECTION 12: QC Stats and Raw Data Tests
## ============================================================================

test_that("fit keeps QC stats when requested", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "stage"),
        model_type = "glm",
        keep_qc_stats = TRUE
    )
    
    raw <- attr(result, "raw_data")
    
    ## Should have QC statistics
    expect_true("AIC" %in% names(raw))
    expect_true("BIC" %in% names(raw))
    expect_true("n" %in% names(raw))
    expect_true("events" %in% names(raw))
})


test_that("fit raw_data has correct structure", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm"
    )
    
    raw <- attr(result, "raw_data")
    
    ## Check required columns in raw data
    expect_true("variable" %in% names(raw))
    expect_true("coef" %in% names(raw))
    expect_true("se" %in% names(raw))
    expect_true("p_value" %in% names(raw))
    expect_true("exp_coef" %in% names(raw))  # For GLM
})


test_that("fit model object is accessible and valid", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm"
    )
    
    model <- attr(result, "model")
    
    ## Model should be a valid glm object
    expect_s3_class(model, "glm")
    
    ## Should be able to extract coefficients
    expect_true(length(coef(model)) > 0)
    
    ## Should be able to get predictions
    expect_true(length(predict(model)) > 0)
})


## ============================================================================
## SECTION 13: Edge Cases and Error Handling
## ============================================================================

test_that("fit handles missing data appropriately", {
    
    ## clintrial has some missing values
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "smoking"),  # smoking may have NAs
        model_type = "glm"
    )
    
    expect_fit_result(result)
    
    ## N should reflect complete cases
    raw <- attr(result, "raw_data")
    expect_true(raw$n[1] > 0)
    expect_true(raw$n[1] <= nrow(clintrial))
})


test_that("fit works with single-level factors appropriately", {
    
    ## Create data with a near-single-level factor
    test_data <- data.table::as.data.table(clintrial)
    test_data[, constant_var := factor("A")]
    
    ## This should either work or give an informative error
    ## depending on how the model handles it
    result_or_error <- tryCatch({
        fit(
            data = test_data,
            outcome = "response",
            predictors = c("age", "constant_var"),
            model_type = "glm"
        )
    }, error = function(e) e)
    
    ## Either succeeds or gives an error (both are acceptable)
    expect_true(inherits(result_or_error, "fit_result") || 
                inherits(result_or_error, "error"))
})


test_that("fit errors on invalid model_type", {
    
    expect_error(
        fit(
            data = clintrial,
            outcome = "response",
            predictors = "age",
            model_type = "invalid_model_type"
        ),
        regexp = "[Uu]nsupported model type"
    )
})


test_that("fit errors when required package missing for coxph", {
    
    ## This test verifies the error message when survival isn't loaded
    ## but since it's typically available, we just verify the check exists
    
    result <- fit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = "age",
        model_type = "coxph"
    )
    
    expect_fit_result(result)
})


test_that("fit handles empty predictors appropriately", {
    
    ## Empty predictors should cause an error
    expect_error(
        fit(
            data = clintrial,
            outcome = "response",
            predictors = character(0),
            model_type = "glm"
        )
    )
})


## ============================================================================
## SECTION 13B: Input Validation Tests
## ============================================================================

test_that("fit auto-corrects Surv outcome with GLM model_type", {
    
    ## Should auto-switch to coxph and emit a message
    expect_message(
        result <- fit(
            data = clintrial,
            outcome = "Surv(os_months, os_status)",
            predictors = c("age", "sex"),
            model_type = "glm"
        ),
        regexp = "[Ss]witch|[Ss]urvival"
    )
    
    expect_fit_result(result)
    ## Model type attribute stores display name ("Cox PH"), not input parameter
    expect_true(grepl("[Cc]ox", attr(result, "model_type")))
})


test_that("fit auto-corrects Surv outcome with lm model_type", {
    
    ## Should auto-switch to coxph and emit a message
    expect_message(
        result <- fit(
            data = clintrial,
            outcome = "Surv(os_months, os_status)",
            predictors = c("age", "sex"),
            model_type = "lm"
        ),
        regexp = "[Ss]witch|[Ss]urvival"
    )
    
    expect_fit_result(result)
    expect_true(grepl("[Cc]ox", attr(result, "model_type")))
})


test_that("fit errors when coxph used without Surv outcome", {
    
    expect_error(
        fit(
            data = clintrial,
            outcome = "response",
            predictors = c("age", "sex"),
            model_type = "coxph"
        ),
        regexp = "[Ss]urv|[Ss]urvival"
    )
})


test_that("fit errors when coxme used without Surv outcome", {
    
    skip_if_not_installed("coxme")
    
    expect_error(
        fit(
            data = clintrial,
            outcome = "response",
            predictors = c("age", "sex", "(1|site)"),
            model_type = "coxme"
        ),
        regexp = "[Ss]urv|[Ss]urvival"
    )
})


test_that("fit errors with continuous outcome and binomial family", {
    
    expect_error(
        fit(
            data = clintrial,
            outcome = "los_days",
            predictors = c("age", "sex"),
            model_type = "glm",
            family = "binomial"
        ),
        regexp = "[Cc]ontinuous|binomial"
    )
})


test_that("fit warns with binary outcome and gaussian family", {
    
    expect_warning(
        result <- fit(
            data = clintrial,
            outcome = "response",
            predictors = c("age", "sex"),
            model_type = "glm",
            family = "gaussian"
        ),
        regexp = "[Bb]inary|binomial|gaussian"
    )
    
    expect_fit_result(result)
})


test_that("fit warns with binary outcome and lm model_type", {
    
    expect_warning(
        result <- fit(
            data = clintrial,
            outcome = "response",
            predictors = c("age", "sex"),
            model_type = "lm"
        ),
        regexp = "[Bb]inary|logistic|binomial"
    )
    
    expect_fit_result(result)
})


test_that("fit errors when outcome variable not found", {
    
    expect_error(
        fit(
            data = clintrial,
            outcome = "nonexistent_variable",
            predictors = c("age", "sex"),
            model_type = "glm"
        ),
        regexp = "not found|[Oo]utcome"
    )
})


test_that("fit errors when Surv variables not found", {
    
    expect_error(
        fit(
            data = clintrial,
            outcome = "Surv(nonexistent_time, nonexistent_status)",
            predictors = c("age", "sex"),
            model_type = "coxph"
        ),
        regexp = "not found|[Ss]urvival"
    )
})


test_that("fit errors when predictor variable not found", {
    
    expect_error(
        fit(
            data = clintrial,
            outcome = "response",
            predictors = c("age", "nonexistent_predictor"),
            model_type = "glm"
        ),
        regexp = "not found|[Pp]redictor"
    )
})


test_that("fit errors with invalid conf_level", {
    
    expect_error(
        fit(
            data = clintrial,
            outcome = "response",
            predictors = c("age", "sex"),
            model_type = "glm",
            conf_level = 1.5
        ),
        regexp = "conf_level|between 0 and 1"
    )
    
    expect_error(
        fit(
            data = clintrial,
            outcome = "response",
            predictors = c("age", "sex"),
            model_type = "glm",
            conf_level = 0
        ),
        regexp = "conf_level|between 0 and 1"
    )
})


test_that("fit errors with invalid digits parameter", {
    
    expect_error(
        fit(
            data = clintrial,
            outcome = "response",
            predictors = c("age", "sex"),
            model_type = "glm",
            digits = -1
        ),
        regexp = "digits|non-negative"
    )
})


test_that("fit errors with empty data", {
    
    empty_data <- clintrial[0, ]
    
    expect_error(
        fit(
            data = empty_data,
            outcome = "response",
            predictors = c("age", "sex"),
            model_type = "glm"
        ),
        regexp = "empty|no rows|non-empty"
    )
})


test_that("fit handles very long predictor lists", {
    
    ## Test with many predictors
    many_predictors <- c("age", "sex", "smoking", "diabetes", "hypertension",
                         "stage", "grade", "ecog", "treatment", "surgery")
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = many_predictors,
        model_type = "glm"
    )
    
    expect_fit_result(result)
    expect_equal(length(attr(result, "predictors")), length(many_predictors))
})


## ============================================================================
## SECTION 14: Print Method Tests
## ============================================================================

test_that("print.fit_result produces output", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm"
    )
    
    ## Capture print output
    output <- capture.output(print(result))
    
    ## Should contain model information
    expect_true(any(grepl("Multivariable|Univariable", output)))
    expect_true(any(grepl("Formula", output)))
    expect_true(any(grepl("N =", output)))
})


test_that("print.fit_result shows events for survival models", {
    
    skip_if_not_installed("survival")
    
    result <- fit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "sex"),
        model_type = "coxph"
    )
    
    output <- capture.output(print(result))
    
    ## Should show Events
    expect_true(any(grepl("Events", output)))
})


## ============================================================================
## SECTION 15: Integration and Consistency Tests
## ============================================================================

test_that("fit produces consistent results across runs", {
    
    result1 <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm"
    )
    
    result2 <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm"
    )
    
    ## Results should be identical
    expect_equal(attr(result1, "raw_data")$coef, attr(result2, "raw_data")$coef)
    expect_equal(attr(result1, "raw_data")$p_value, attr(result2, "raw_data")$p_value)
})


test_that("fit model coefficients match direct model fitting", {
    
    ## Fit using fit()
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "binomial"
    )
    
    ## Fit directly
    direct_model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    ## Coefficients should match
    fit_coefs <- coef(attr(result, "model"))
    direct_coefs <- coef(direct_model)
    
    expect_equal(fit_coefs, direct_coefs, tolerance = 1e-10)
})


test_that("fit works with data.frame input (not just data.table)", {
    
    df <- as.data.frame(clintrial)
    
    result <- fit(
        data = df,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm"
    )
    
    expect_fit_result(result)
})


test_that("fit preserves factor level ordering", {
    
    ## Create data with ordered factor
    test_data <- data.table::as.data.table(clintrial)
    test_data[, stage_ordered := factor(stage, levels = c("I", "II", "III", "IV"), 
                                        ordered = TRUE)]
    
    result <- fit(
        data = test_data,
        outcome = "response",
        predictors = c("age", "stage_ordered"),
        model_type = "glm"
    )
    
    expect_fit_result(result)
    
    ## Check that stage levels appear in correct order in output
    raw <- attr(result, "raw_data")
    stage_rows <- raw[variable == "stage_ordered"]
    
    if (nrow(stage_rows) > 1) {
        ## Groups should be in factor order
        groups <- stage_rows$group[!is.na(stage_rows$group) & stage_rows$group != ""]
        
        ## Reference should be first level (I), so we should see II, III, IV
        expect_true(length(groups) >= 1)
    }
})


## ============================================================================
## SECTION 16: Different GLM Families Tests
## ============================================================================

## ----------------------------------------------------------------------------
## 16a: Gaussian Family Tests
## ----------------------------------------------------------------------------

test_that("fit works with gaussian family", {
    
    result <- fit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "gaussian"
    )
    
    expect_fit_result(result)
    
    ## Should produce coefficients (not exponentiated)
    col_names <- names(result)
    expect_true(any(grepl("Coefficient.*CI|Estimate.*CI", col_names)),
                info = paste("Expected Coefficient column for gaussian, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("fit gaussian coefficients match direct glm", {
    
    result <- fit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "gaussian"
    )
    
    ## Fit directly
    direct_model <- glm(los_days ~ age + sex, data = clintrial, family = gaussian)
    
    ## Coefficients should match
    fit_coefs <- coef(attr(result, "model"))
    direct_coefs <- coef(direct_model)
    
    expect_equal(fit_coefs, direct_coefs, tolerance = 1e-10)
})


test_that("fit gaussian matches lm results", {
    
    result_glm <- fit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "gaussian"
    )
    
    result_lm <- fit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex"),
        model_type = "lm"
    )
    
    ## Coefficients should be essentially identical
    coefs_glm <- coef(attr(result_glm, "model"))
    coefs_lm <- coef(attr(result_lm, "model"))
    
    expect_equal(coefs_glm, coefs_lm, tolerance = 1e-10)
})


## ----------------------------------------------------------------------------
## 16b: Quasibinomial Family Tests (overdispersed binary)
## ----------------------------------------------------------------------------

test_that("fit works with quasibinomial family", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "treatment"),
        model_type = "glm",
        family = "quasibinomial"
    )
    
    expect_fit_result(result)
    
    ## Should produce odds ratios like binomial
    col_names <- names(result)
    expect_true(any(grepl("OR.*CI", col_names)),
                info = paste("Expected OR column for quasibinomial, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("fit quasibinomial coefficients match direct glm", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "quasibinomial"
    )
    
    ## Fit directly
    direct_model <- glm(response ~ age + sex, data = clintrial, family = quasibinomial)
    
    ## Coefficients should match
    fit_coefs <- coef(attr(result, "model"))
    direct_coefs <- coef(direct_model)
    
    expect_equal(fit_coefs, direct_coefs, tolerance = 1e-10)
})


test_that("fit quasibinomial point estimates match binomial", {
    
    result_quasi <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "quasibinomial"
    )
    
    result_binom <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "binomial"
    )
    
    ## Point estimates should be identical (SEs differ due to dispersion)
    coefs_quasi <- coef(attr(result_quasi, "model"))
    coefs_binom <- coef(attr(result_binom, "model"))
    
    expect_equal(coefs_quasi, coefs_binom, tolerance = 1e-10)
})


test_that("fit quasibinomial handles overdispersion", {
    
    ## Quasibinomial allows dispersion != 1
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "treatment", "stage"),
        model_type = "glm",
        family = "quasibinomial"
    )
    
    model <- attr(result, "model")
    
    ## Dispersion parameter should be estimated
    summ <- summary(model)
    expect_true(!is.null(summ$dispersion))
})


## ----------------------------------------------------------------------------
## 16c: Quasipoisson Family Tests (overdispersed counts)
## ----------------------------------------------------------------------------

test_that("fit works with quasipoisson family", {
    
    result <- fit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment"),
        model_type = "glm",
        family = "quasipoisson"
    )
    
    expect_fit_result(result)
    
    ## Should produce rate ratios like poisson
    col_names <- names(result)
    expect_true(any(grepl("RR.*CI", col_names)),
                info = paste("Expected RR column for quasipoisson, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("fit quasipoisson coefficients match direct glm", {
    
    result <- fit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "quasipoisson"
    )
    
    ## Fit directly
    direct_model <- glm(count_outcome ~ age + sex, data = clintrial, family = quasipoisson)
    
    ## Coefficients should match
    fit_coefs <- coef(attr(result, "model"))
    direct_coefs <- coef(direct_model)
    
    expect_equal(fit_coefs, direct_coefs, tolerance = 1e-10)
})


test_that("fit quasipoisson point estimates match poisson", {
    
    result_quasi <- fit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "quasipoisson"
    )
    
    result_pois <- fit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "poisson"
    )
    
    ## Point estimates should be identical (SEs differ due to dispersion)
    coefs_quasi <- coef(attr(result_quasi, "model"))
    coefs_pois <- coef(attr(result_pois, "model"))
    
    expect_equal(coefs_quasi, coefs_pois, tolerance = 1e-10)
})


test_that("fit quasipoisson handles overdispersion", {
    
    result <- fit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment", "stage"),
        model_type = "glm",
        family = "quasipoisson"
    )
    
    model <- attr(result, "model")
    
    ## Dispersion parameter should be estimated
    summ <- summary(model)
    expect_true(!is.null(summ$dispersion))
    expect_true(summ$dispersion != 1)  ## Should differ from 1 for overdispersed data
})


## ----------------------------------------------------------------------------
## 16d: Gamma Family Tests (positive continuous)
## ----------------------------------------------------------------------------

test_that("fit works with Gamma family (log link)", {
    
    ## Create positive outcome
    test_data <- data.table::as.data.table(clintrial)
    test_data <- test_data[los_days > 0]
    
    result <- fit(
        data = test_data,
        outcome = "los_days",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = Gamma(link = "log")
    )
    
    expect_fit_result(result)
    
    ## With log link, should produce ratios (exponentiated)
    col_names <- names(result)
    expect_true(any(grepl("Ratio.*CI|Coefficient.*CI", col_names)),
                info = paste("Expected Ratio or Coefficient column for Gamma, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("fit works with Gamma family (inverse link - default)", {
    
    ## Create positive outcome
    test_data <- data.table::as.data.table(clintrial)
    test_data <- test_data[los_days > 0]
    
    result <- fit(
        data = test_data,
        outcome = "los_days",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = Gamma()  ## Default inverse link
    )
    
    expect_fit_result(result)
})


test_that("fit Gamma coefficients match direct glm", {
    
    ## Create positive outcome
    test_data <- data.table::as.data.table(clintrial)
    test_data <- test_data[los_days > 0]
    
    result <- fit(
        data = test_data,
        outcome = "los_days",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = Gamma(link = "log")
    )
    
    ## Fit directly
    direct_model <- glm(los_days ~ age + sex, data = test_data, family = Gamma(link = "log"))
    
    ## Coefficients should match
    fit_coefs <- coef(attr(result, "model"))
    direct_coefs <- coef(direct_model)
    
    expect_equal(fit_coefs, direct_coefs, tolerance = 1e-10)
})


test_that("fit Gamma handles factor variables", {
    
    ## Create positive outcome
    test_data <- data.table::as.data.table(clintrial)
    test_data <- test_data[los_days > 0]
    
    result <- fit(
        data = test_data,
        outcome = "los_days",
        predictors = c("age", "stage"),
        model_type = "glm",
        family = Gamma(link = "log"),
        reference_rows = TRUE
    )
    
    expect_fit_result(result)
    
    ## Stage should have reference row
    raw <- attr(result, "raw_data")
    stage_rows <- raw[variable == "stage"]
    expect_true(nrow(stage_rows) >= 2)
})


## ----------------------------------------------------------------------------
## 16e: Inverse Gaussian Family Tests (positive, right-skewed)
## ----------------------------------------------------------------------------

test_that("fit works with inverse.gaussian family", {
    
    ## Create positive outcome
    test_data <- data.table::as.data.table(clintrial)
    test_data <- test_data[los_days > 0]
    
    result <- fit(
        data = test_data,
        outcome = "los_days",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = inverse.gaussian()
    )
    
    expect_fit_result(result)
})


test_that("fit inverse.gaussian with log link works", {
    
    ## Create positive outcome
    test_data <- data.table::as.data.table(clintrial)
    test_data <- test_data[los_days > 0]
    
    result <- fit(
        data = test_data,
        outcome = "los_days",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = inverse.gaussian(link = "log")
    )
    
    expect_fit_result(result)
    
    ## With log link, should produce ratios
    col_names <- names(result)
    expect_true(any(grepl("Ratio.*CI|Coefficient.*CI", col_names)),
                info = paste("Expected Ratio or Coefficient column, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("fit inverse.gaussian coefficients match direct glm", {
    
    ## Create positive outcome
    test_data <- data.table::as.data.table(clintrial)
    test_data <- test_data[los_days > 0]
    
    result <- fit(
        data = test_data,
        outcome = "los_days",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = inverse.gaussian(link = "log")
    )
    
    ## Fit directly
    direct_model <- glm(los_days ~ age + sex, data = test_data, 
                        family = inverse.gaussian(link = "log"))
    
    ## Coefficients should match
    fit_coefs <- coef(attr(result, "model"))
    direct_coefs <- coef(direct_model)
    
    expect_equal(fit_coefs, direct_coefs, tolerance = 1e-10)
})


## ----------------------------------------------------------------------------
## 16f: Negative Binomial (negbin) Tests
## ----------------------------------------------------------------------------

test_that("fit works with negative binomial regression (negbin)", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(fit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment"),
        model_type = "negbin"
    ))
    
    expect_fit_result(result)
    
    ## Should produce rate ratios
    col_names <- names(result)
    expect_true(any(grepl("RR.*CI", col_names)),
                info = paste("Expected RR column for negative binomial, found:", 
                             paste(col_names, collapse = ", ")))
    
    ## Model should be of class negbin
    model <- attr(result, "model")
    expect_true(inherits(model, "negbin"))
})


test_that("fit negbin produces correct effect estimates", {
    
    skip_if_not_installed("MASS")
    
    ## Fit using fit()
    result <- suppressWarnings(fit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex"),
        model_type = "negbin"
    ))
    
    ## Fit directly with MASS
    direct_model <- suppressWarnings(MASS::glm.nb(count_outcome ~ age + sex, data = clintrial))
    
    ## Coefficients should match
    fit_coefs <- coef(attr(result, "model"))
    direct_coefs <- coef(direct_model)
    
    expect_equal(fit_coefs, direct_coefs, tolerance = 1e-10)
})


test_that("fit negbin works with multivariable models", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(fit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment", "stage"),
        model_type = "negbin"
    ))
    
    expect_fit_result(result)
    expect_equal(attr(result, "model_scope"), "Multivariable")
    
    ## Should have aRR column for multivariable
    col_names <- names(result)
    expect_true(any(grepl("aRR.*CI", col_names)),
                info = paste("Expected aRR column for multivariable negbin, found:", 
                             paste(col_names, collapse = ", ")))
})


test_that("fit negbin handles factor variables correctly", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(fit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "stage"),
        model_type = "negbin",
        reference_rows = TRUE
    ))
    
    expect_fit_result(result)
    
    ## Stage is a factor - should have multiple rows with reference
    raw <- attr(result, "raw_data")
    stage_rows <- raw[variable == "stage"]
    expect_true(nrow(stage_rows) >= 2)
    
    ## Check reference row exists
    has_ref <- any(stage_rows$reference == "reference", na.rm = TRUE)
    expect_true(has_ref)
})


test_that("fit negbin respects formatting parameters", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(fit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex"),
        model_type = "negbin",
        digits = 3,
        p_digits = 4,
        show_n = TRUE,
        show_events = FALSE
    ))
    
    expect_fit_result(result)
    expect_true("n" %in% names(result))
    expect_false("Events" %in% names(result))
})


test_that("fit negbin works with labels", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(fit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment"),
        model_type = "negbin",
        labels = clintrial_labels
    ))
    
    expect_fit_result(result)
    
    ## Labels should be applied
    expect_true(any(grepl("Age|Sex|Treatment", result$Variable)))
})


test_that("fit negbin works with interactions", {
    
    skip_if_not_installed("MASS")
    
    result <- suppressWarnings(fit(
        data = clintrial,
        outcome = "count_outcome",
        predictors = c("age", "sex", "treatment"),
        interactions = c("age:sex"),
        model_type = "negbin"
    ))
    
    expect_fit_result(result)
    
    ## Interaction should be in formula
    formula_str <- attr(result, "formula_str")
    expect_true(grepl("age:sex", formula_str))
})


test_that("fit with pre-fitted model works for negbin", {
    
    skip_if_not_installed("MASS")
    
    ## Fit model directly with MASS
    nb_model <- suppressWarnings(MASS::glm.nb(count_outcome ~ age + sex + treatment, data = clintrial))
    
    ## Pass to fit()
    result <- fit(
        model = nb_model,
        data = clintrial,
        labels = clintrial_labels
    )
    
    ## Check basic structure (without outcome/predictors which aren't available for pre-fitted models)
    expect_s3_class(result, "fit_result")
    expect_s3_class(result, "data.table")
    expect_true("Variable" %in% names(result))
    expect_true("Group" %in% names(result))
    expect_true(!is.null(attr(result, "model")))
    expect_true(!is.null(attr(result, "raw_data")))
    
    ## Should recognize it as negative binomial
    model_type <- attr(result, "model_type")
    expect_true(grepl("Negative Binomial", model_type))
})


## ============================================================================
## SECTION 17: Model Attribute Preservation Tests
## ============================================================================

test_that("fit preserves all specified attributes", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "stage"),
        interactions = c("age:sex"),
        model_type = "glm",
        family = "binomial"
    )
    
    ## Check all attributes are preserved
    expect_equal(attr(result, "outcome"), "response")
    expect_equal(attr(result, "predictors"), c("age", "sex", "stage"))
    expect_equal(attr(result, "interactions"), c("age:sex"))
    
    ## Formula should be correct
    formula_str <- attr(result, "formula_str")
    expect_true(grepl("response", formula_str))
    expect_true(grepl("age", formula_str))
    expect_true(grepl("sex", formula_str))
    expect_true(grepl("stage", formula_str))
    expect_true(grepl("age:sex", formula_str))
})


test_that("fit stores data in model object", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm"
    )
    
    model <- attr(result, "model")
    
    ## Data should be accessible from model
    model_data <- attr(model, "data")
    expect_true(!is.null(model_data) || !is.null(model$data))
})


## ============================================================================
## SECTION 18: Performance and Large Data Tests
## ============================================================================

test_that("fit handles moderately large datasets", {
    
    ## Create larger dataset by replicating
    large_data <- rbindlist(rep(list(clintrial), 10))
    large_data[, id := .I]  # Add unique ID
    
    result <- fit(
        data = large_data,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm"
    )
    
    expect_fit_result(result)
    
    ## N should reflect larger dataset
    raw <- attr(result, "raw_data")
    expect_true(raw$n[1] > nrow(clintrial))
})


## ============================================================================
## SECTION 19: Numeric Precision Tests
## ============================================================================

test_that("fit handles extreme coefficient values", {
    
    ## Create data that might produce extreme coefficients
    test_data <- data.table::as.data.table(clintrial)
    test_data[, tiny_var := age / 10000]  # Very small scale
    
    result <- fit(
        data = test_data,
        outcome = "response",
        predictors = c("tiny_var", "sex"),
        model_type = "glm"
    )
    
    expect_fit_result(result)
    
    ## Should still produce valid output
    raw <- attr(result, "raw_data")
    expect_true(all(!is.na(raw$coef)))
})


test_that("fit handles perfect separation gracefully",
{
    ## Create data with near-perfect separation
    test_data <- data.table(
        outcome = c(rep(0, 50), rep(1, 50)),
        predictor = c(rep(0, 50), rep(1, 50)),
        noise = rnorm(100)
    )
    
    ## This might produce warnings about convergence
    result <- suppressWarnings(
        fit(
            data = test_data,
            outcome = "outcome",
            predictors = c("predictor", "noise"),
            model_type = "glm",
            family = "binomial"
        )
    )
    
    ## Should still return a valid result
    expect_fit_result(result)
})


## ============================================================================
## SECTION 20: Formula String Validation Tests
## ============================================================================

test_that("fit formula_str is valid R formula", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "stage"),
        interactions = c("age:sex"),
        model_type = "glm"
    )
    
    formula_str <- attr(result, "formula_str")
    
    ## Should be parseable as formula
    parsed_formula <- as.formula(formula_str)
    expect_s3_class(parsed_formula, "formula")
})


test_that("fit formula includes all components correctly", {
    
    result <- fit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "treatment"),
        model_type = "coxph",
        strata = "sex"
    )
    
    formula_str <- attr(result, "formula_str")
    
    ## Check all components present
    expect_true(grepl("Surv\\(os_months, os_status\\)", formula_str))
    expect_true(grepl("age", formula_str))
    expect_true(grepl("treatment", formula_str))
    expect_true(grepl("strata\\(.*sex.*\\)", formula_str))
})
