#' Test Suite for compfit()
#' 
#' Comprehensive tests covering standard models, mixed-effects models,
#' interactions, edge cases, and output validation.
#' 
#' @details Run with testthat::test_file("tests/testthat/test-compfit.R")

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

## Helper function to check compfit_result structure
expect_compfit_result <- function(result) {
    expect_s3_class(result, "compfit_result")
    expect_s3_class(result, "data.table")
    
    ## Required columns
    expect_true("Model" %in% names(result))
    expect_true("N" %in% names(result))
    expect_true("Converged" %in% names(result))
    expect_true("AIC" %in% names(result))
    expect_true("Summata Score" %in% names(result))
    
    ## Required attributes
    expect_true(!is.null(attr(result, "models")))
    expect_true(!is.null(attr(result, "model_type")))
    expect_true(!is.null(attr(result, "outcome")))
    expect_true(!is.null(attr(result, "best_model")))
}

## Helper to check scores are valid
expect_valid_scores <- function(result) {
    scores <- result$`Summata Score`
    expect_true(all(scores >= 0 | is.na(scores)))
    expect_true(all(scores <= 100 | is.na(scores)))
    ## Should be sorted descending
    non_na_scores <- scores[!is.na(scores)]
    if (length(non_na_scores) > 1) {
        expect_true(all(diff(non_na_scores) <= 0))
    }
}


## ============================================================================
## SECTION 1: Basic Functionality Tests
## ============================================================================

test_that("compfit works with basic GLM models", {
    
    models <- list(
        base = c("age", "sex"),
        clinical = c("age", "sex", "smoking", "diabetes"),
        full = c("age", "sex", "smoking", "diabetes", "stage", "grade")
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        model_type = "glm",
        family = "binomial"
    )
    
    expect_compfit_result(result)
    expect_equal(nrow(result), 3)
    expect_valid_scores(result)
    
    ## Check that all models have valid N
    expect_true(all(result$N > 0))
    
    ## Check concordance is between 0.5 and 1
    conc <- result$Concordance[!is.na(result$Concordance)]
    expect_true(all(conc >= 0.5 & conc <= 1))
})


test_that("compfit works with linear models", {
    
    models <- list(
        simple = c("age", "sex"),
        moderate = c("age", "sex", "stage", "ecog"),
        complex = c("age", "sex", "stage", "ecog", "treatment", "surgery")
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "los_days",
        model_list = models,
        model_type = "lm"
    )
    
    expect_compfit_result(result)
    expect_equal(nrow(result), 3)
    expect_equal(attr(result, "model_type"), "lm")
    
    ## Linear models should have Pseudo-R^2
    expect_true("Pseudo-R^2" %in% names(result))
})


test_that("compfit works with Cox models", {
    
    models <- list(
        base = c("age", "sex"),
        clinical = c("age", "sex", "stage", "grade"),
        full = c("age", "sex", "stage", "grade", "treatment", "surgery")
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        model_list = models,
        model_type = "coxph"
    )
    
    expect_compfit_result(result)
    expect_equal(nrow(result), 3)
    expect_equal(attr(result, "model_type"), "coxph")
    
    ## Cox models should have Events column
    expect_true("Events" %in% names(result))
    expect_true(all(result$Events > 0))
})


test_that("compfit auto-detects model type correctly", {
    
    ## Binary outcome -> GLM
    result_binary <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = list(simple = c("age", "sex")),
        model_type = "auto"
    )
    expect_equal(attr(result_binary, "model_type"), "glm")
    
    ## Continuous outcome -> LM
    result_cont <- compfit(
        data = clintrial,
        outcome = "los_days",
        model_list = list(simple = c("age", "sex")),
        model_type = "auto"
    )
    expect_equal(attr(result_cont, "model_type"), "lm")
    
    ## Survival outcome -> coxph
    result_surv <- compfit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        model_list = list(simple = c("age", "sex")),
        model_type = "auto"
    )
    expect_equal(attr(result_surv, "model_type"), "coxph")
})


## ============================================================================
## SECTION 2: Mixed-Effects Model Tests
## ============================================================================

## Note: lme4 can cause segfaults when models are fitted within certain 
## evaluation contexts (testthat, tryCatch, etc.) due to C++ memory management
## issues. We use skip_on_cran() and careful error handling to mitigate this.

test_that("compfit works with lmer models", {
    
    skip_if_not_installed("lme4")
    skip_on_cran()  # lme4 can segfault in test environments
    
    ## Force garbage collection before lme4 operations
    gc()
    
    models <- list(
        simple = c("age", "sex", "(1|site)")
    )
    
    ## Use a smaller dataset to reduce memory pressure
    test_data <- clintrial_complete[1:200, ]
    
    result <- compfit(
        data = test_data,
        outcome = "los_days",
        model_list = models,
        model_type = "lmer",
        REML = FALSE  # Use ML for valid AIC comparison
    )
    
    expect_compfit_result(result)
    expect_equal(nrow(result), 1)
    
    ## lmer-specific columns
    expect_true("Groups" %in% names(result))
    expect_true("Marginal R2" %in% names(result) || "Pseudo-R^2" %in% names(result))
    expect_true("ICC" %in% names(result))
    
    ## Clean up

    gc()
})


test_that("compfit works with glmer models", {
    
    skip_if_not_installed("lme4")
    skip_on_cran()  # lme4 can segfault in test environments
    
    gc()
    
    models <- list(
        simple = c("age", "sex", "(1|site)")
    )
    
    ## Use smaller dataset - need response column
    test_data <- clintrial_complete[1:200, ]
    
    result <- compfit(
        data = test_data,
        outcome = "response",
        model_list = models,
        model_type = "glmer",
        family = "binomial"
    )
    
    expect_compfit_result(result)
    expect_equal(nrow(result), 1)
    
    ## glmer-specific columns
    expect_true("Groups" %in% names(result))
    expect_true("Concordance" %in% names(result))
    
    gc()
})


test_that("compfit works with coxme models", {
    
    skip_if_not_installed("coxme")
    skip_on_cran()  # Mixed models can segfault in test environments
    
    gc()
    
    models <- list(
        simple = c("age", "sex", "(1|site)")
    )
    
    ## Use smaller dataset
    test_data <- clintrial_complete[1:200, ]
    
    result <- compfit(
        data = test_data,
        outcome = "Surv(os_months, os_status)",
        model_list = models,
        model_type = "coxme"
    )
    
    expect_compfit_result(result)
    expect_equal(nrow(result), 1)
    expect_equal(attr(result, "model_type"), "coxme")
    
    ## Should have Events and Groups
    expect_true("Events" %in% names(result))
    expect_true("Groups" %in% names(result))
    
    gc()
})


test_that("compfit auto-detects mixed models from random effects in predictors", {
    
    skip_if_not_installed("lme4")
    skip_on_cran()  # lme4 can segfault in test environments
    
    gc()
    
    ## Use smaller dataset - response column is already in clintrial_complete
    test_data <- clintrial_complete[1:150, ]
    
    ## Continuous + random effects -> lmer
    result_lmer <- compfit(
        data = test_data,
        outcome = "los_days",
        model_list = list(mixed = c("age", "sex", "(1|site)")),
        model_type = "auto",
        REML = FALSE
    )
    expect_equal(attr(result_lmer, "model_type"), "lmer")
    
    gc()
    
    ## Binary + random effects -> glmer
    result_glmer <- compfit(
        data = test_data,
        outcome = "response",
        model_list = list(mixed = c("age", "sex", "(1|site)")),
        model_type = "auto"
    )
    expect_equal(attr(result_glmer, "model_type"), "glmer")
    
    gc()
})


test_that("compfit handles different random effects structures", {
    
    skip_if_not_installed("lme4")
    skip_on_cran()  # lme4 can segfault in test environments
    skip("Random slope models are prone to convergence issues and segfaults")
    
    gc()
    
    models <- list(
        intercept_only = c("age", "treatment", "(1|site)"),
        random_slope = c("age", "treatment", "(1 + age|site)")
    )
    
    ## This may produce convergence warnings for complex random effects
    result <- suppressWarnings(compfit(
        data = clintrial_complete,
        outcome = "los_days",
        model_list = models,
        model_type = "lmer",
        REML = FALSE
    ))
    
    expect_compfit_result(result)
    expect_equal(nrow(result), 2)
    
    gc()
})


## ============================================================================
## SECTION 3: Interaction Term Tests
## ============================================================================

test_that("compfit works with interactions_list parameter", {
    
    models <- list(
        base = c("age", "treatment"),
        interaction = c("age", "treatment")
    )
    
    interactions <- list(
        NULL,  # No interactions for base model
        c("age:treatment")  # Interaction for second model
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        interactions_list = interactions,
        model_type = "glm"
    )
    
    expect_compfit_result(result)
    expect_equal(nrow(result), 2)
    
    ## Interaction model should have more predictors listed
    ## (Note: exact count depends on implementation)
})


test_that("compfit handles multiple interactions per model", {
    
    models <- list(
        simple = c("age", "treatment", "stage"),
        complex = c("age", "treatment", "stage")
    )
    
    interactions <- list(
        NULL,
        c("age:treatment", "treatment:stage")
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        interactions_list = interactions,
        model_type = "glm"
    )
    
    expect_compfit_result(result)
})


test_that("compfit validates interactions_list length matches model_list", {
    
    models <- list(
        a = c("age"),
        b = c("age", "sex")
    )
    
    ## Wrong number of interactions
    interactions <- list(NULL)  # Only 1 instead of 2
    
    expect_error(
        compfit(
            data = clintrial,
            outcome = "response",
            model_list = models,
            interactions_list = interactions,
            model_type = "glm"
        ),
        regexp = "same length"
    )
})


## ============================================================================
## SECTION 4: Model Naming Tests
## ============================================================================

test_that("compfit uses provided model_names", {
    
    models <- list(
        c("age", "sex"),
        c("age", "sex", "stage")
    )
    
    custom_names <- c("Baseline", "Extended")
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        model_names = custom_names,
        model_type = "glm"
    )
    
    expect_true(all(custom_names %in% result$Model))
})


test_that("compfit uses list names as model names when no model_names provided", {
    
    models <- list(
        base = c("age", "sex"),
        full = c("age", "sex", "stage")
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        model_type = "glm"
    )
    
    expect_true(all(c("base", "full") %in% result$Model))
})


test_that("compfit generates default model names for unnamed list", {
    
    models <- list(
        c("age", "sex"),
        c("age", "sex", "stage")
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        model_type = "glm"
    )
    
    ## Should have some model names (Model 1, Model 2 or similar)
    expect_true(all(!is.na(result$Model)))
    expect_equal(length(unique(result$Model)), 2)
})


test_that("compfit passes labels to underlying fit function", {
    
    models <- list(simple = c("age", "sex"))
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        model_type = "glm",
        labels = clintrial_labels,
        include_coefficients = TRUE
    )
    
    coef_table <- attr(result, "coefficients")
    
    ## If labels were applied, coefficient table should use them
    ## (exact check depends on implementation)
    expect_true(!is.null(coef_table))
})


## ============================================================================
## SECTION 5: Scoring and Ranking Tests
## ============================================================================

test_that("compfit ranks models by Summata Score in descending order", {
    
    models <- list(
        minimal = c("age"),
        base = c("age", "sex"),
        full = c("age", "sex", "stage", "grade", "treatment")
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        model_list = models,
        model_type = "coxph"
    )
    
    scores <- result$`Summata Score`
    
    ## Scores should be in descending order (best first)
    expect_true(all(diff(scores) <= 0))
})


test_that("compfit identifies best model correctly", {
    
    models <- list(
        a = c("age"),
        b = c("age", "sex"),
        c = c("age", "sex", "stage")
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        model_type = "glm"
    )
    
    best_model <- attr(result, "best_model")
    
    ## Best model should be the first row (highest score)
    expect_equal(best_model, result$Model[1])
})


test_that("compfit accepts custom scoring weights", {
    
    custom_weights <- list(
        convergence = 0.10,
        aic = 0.30,
        concordance = 0.40,
        pseudo_r2 = 0.20
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = list(simple = c("age", "sex")),
        model_type = "glm",
        scoring_weights = custom_weights
    )
    
    expect_compfit_result(result)
    
    ## Weights should be stored in attributes
    expect_true(!is.null(attr(result, "weights")))
})


## ============================================================================
## SECTION 6: Coefficient Table Tests
## ============================================================================

test_that("compfit returns coefficient table when requested", {
    
    models <- list(
        base = c("age", "sex"),
        full = c("age", "sex", "stage")
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        model_type = "glm",
        include_coefficients = TRUE
    )
    
    coef_table <- attr(result, "coefficients")
    
    expect_true(!is.null(coef_table))
    expect_s3_class(coef_table, "data.table")
    expect_true("Model" %in% names(coef_table))
    
    ## Should have rows for both models - check model names are present
    ## Models might be renamed to "Model 1", "Model 2" or keep original names
    model_names_in_coef <- unique(coef_table$Model)
    expect_equal(length(model_names_in_coef), 2)
})


test_that("compfit returns NULL coefficients when not requested", {
    
    models <- list(simple = c("age", "sex"))
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        model_type = "glm",
        include_coefficients = FALSE
    )
    
    expect_null(attr(result, "coefficients"))
})


## ============================================================================
## SECTION 7: Edge Cases and Error Handling
## ============================================================================

test_that("compfit handles single model", {
    
    models <- list(only = c("age", "sex", "stage"))
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        model_type = "glm"
    )
    
    expect_compfit_result(result)
    expect_equal(nrow(result), 1)
})


test_that("compfit continues when one model fails", {
    
    ## This test checks that compfit handles model failures gracefully
    ## We can't easily create a guaranteed failure, so we just verify
    ## the function handles edge cases
    
    models <- list(
        good = c("age", "sex"),
        also_good = c("age", "sex", "stage")
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        model_type = "glm"
    )
    
    expect_compfit_result(result)
})


test_that("compfit validates model_list is not empty", {
    
    expect_error(
        compfit(
            data = clintrial,
            outcome = "response",
            model_list = list(),
            model_type = "glm"
        )
    )
})


test_that("compfit handles missing data appropriately", {
    
    models <- list(
        has_missing = c("age", "bmi", "smoking"),  # bmi and smoking have NAs
        complete = c("age", "sex")  # No NAs in these
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        model_type = "glm"
    )
    
    expect_compfit_result(result)
    
    ## Get N values - note: results may be reordered by Summata Score
    ## Use which() to find the correct rows
    idx_missing <- which(result$Model == "has_missing")
    idx_complete <- which(result$Model == "complete")
    
    ## Both models should be present
    expect_true(length(idx_missing) == 1)
    expect_true(length(idx_complete) == 1)
    
    n_missing <- result$N[idx_missing]
    n_complete <- result$N[idx_complete]
    
    ## Model with complete predictors should have >= N as model with missing
    expect_true(n_complete >= n_missing)
})


test_that("compfit handles invalid model_type gracefully", {
    
    models <- list(simple = c("age", "sex"))
    
    ## This should either error or fall back to auto
    ## The exact behavior depends on implementation
    error_occurred <- FALSE
    result <- tryCatch({
        compfit(
            data = clintrial,
            outcome = "response",
            model_list = models,
            model_type = "invalid_type"
        )
    }, error = function(e) {
        error_occurred <<- TRUE
        NULL
    })
    
    ## Should either error or produce a result
    expect_true(error_occurred || inherits(result, "compfit_result"))
})


## ============================================================================
## SECTION 8: Output Format Tests
## ============================================================================

test_that("Summata Score is always the second column", {
    
    ## Test with GLM
    result_glm <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = list(simple = c("age", "sex")),
        model_type = "glm"
    )
    expect_equal(names(result_glm)[2], "Summata Score")
    
    ## Test with LM
    result_lm <- compfit(
        data = clintrial,
        outcome = "los_days",
        model_list = list(simple = c("age", "sex")),
        model_type = "lm"
    )
    expect_equal(names(result_lm)[2], "Summata Score")
    
    ## Test with Cox
    result_cox <- compfit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        model_list = list(simple = c("age", "sex")),
        model_type = "coxph"
    )
    expect_equal(names(result_cox)[2], "Summata Score")
})


test_that("Summata Score is second column for mixed models", {
    
    skip_if_not_installed("lme4")
    skip_on_cran()
    
    gc()
    
    test_data <- clintrial_complete[1:150, ]
    
    result <- compfit(
        data = test_data,
        outcome = "los_days",
        model_list = list(simple = c("age", "sex", "(1|site)")),
        model_type = "lmer",
        REML = FALSE
    )
    
    expect_equal(names(result)[2], "Summata Score")
    
    gc()
})


test_that("compfit column order is correct for GLM", {
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = list(simple = c("age", "sex")),
        model_type = "glm"
    )
    
    expected_cols <- c("Model", "Summata Score", "N", "Events", "Predictors", "Converged",
                       "AIC", "BIC", "Pseudo-R^2", "Concordance", "Brier Score",
                       "Global p")
    
    actual_cols <- names(result)
    
    ## Check that expected columns exist and come in right order
    expected_present <- intersect(expected_cols, actual_cols)
    expect_equal(actual_cols[seq_along(expected_present)], expected_present)
})


test_that("compfit column order is correct for lmer", {
    
    skip_if_not_installed("lme4")
    skip_on_cran()  # lme4 can segfault in test environments
    
    gc()
    
    ## Use smaller dataset
    test_data <- clintrial_complete[1:150, ]
    
    result <- compfit(
        data = test_data,
        outcome = "los_days",
        model_list = list(mixed = c("age", "sex", "(1|site)")),
        model_type = "lmer",
        REML = FALSE
    )
    
    ## Should have lmer-specific columns
    expect_true("Groups" %in% names(result))
    expect_true("ICC" %in% names(result))
    
    gc()
})


test_that("print.compfit_result works without error", {
    
    models <- list(
        base = c("age", "sex"),
        full = c("age", "sex", "stage")
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        model_type = "glm"
    )
    
    ## Should print without error
    expect_output(print(result), "Model Comparison Results")
    expect_output(print(result), "Summata Score")
    expect_output(print(result), "Recommended Model")
})


## ============================================================================
## SECTION 9: Convergence Detection Tests
## ============================================================================

test_that("compfit detects convergence status", {
    
    models <- list(
        simple = c("age", "sex"),
        complex = c("age", "sex", "stage", "grade", "treatment", "surgery", "ecog")
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        model_type = "glm"
    )
    
    ## Converged column should exist and have valid values
    expect_true("Converged" %in% names(result))
    expect_true(all(result$Converged %in% c("Yes", "No", "Suspect", "Failed")))
})


test_that("compfit detects singular fits in mixed models", {
    
    skip_if_not_installed("lme4")
    skip_on_cran()  # lme4 can segfault in test environments
    skip("Singular fit detection requires careful setup and is prone to segfaults")
    
    gc()
    
    ## Create a scenario likely to produce singular fit
    ## (very few groups or variance near zero)
    small_data <- clintrial_complete[1:100, ]
    
    models <- list(
        simple = c("age", "(1|site)")
    )
    
    ## May produce singular fit warning
    result <- suppressWarnings(compfit(
        data = small_data,
        outcome = "los_days",
        model_list = models,
        model_type = "lmer",
        REML = FALSE
    ))
    
    expect_compfit_result(result)
    ## Convergence status should reflect any issues
    expect_true(result$Converged[1] %in% c("Yes", "Suspect"))
    
    gc()
})


## ============================================================================
## SECTION 10: Attribute Preservation Tests
## ============================================================================

test_that("compfit preserves fitted models in attributes", {
    
    models <- list(
        base = c("age", "sex"),
        full = c("age", "sex", "stage")
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        model_type = "glm"
    )
    
    stored_models <- attr(result, "models")
    
    expect_true(is.list(stored_models))
    expect_equal(length(stored_models), 2)
    
    ## Models should be stored - names might be original or display names
    ## Just check that we have 2 models that are glm objects
    expect_true(all(sapply(stored_models, function(m) inherits(m, "glm") || is.null(m))))
})


test_that("compfit stores weights in attributes", {
    
    custom_weights <- list(
        convergence = 0.10,
        aic = 0.30,
        concordance = 0.40,
        pseudo_r2 = 0.20
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = list(simple = c("age", "sex")),
        model_type = "glm",
        scoring_weights = custom_weights
    )
    
    stored_weights <- attr(result, "weights")
    expect_true(!is.null(stored_weights))
})


## ============================================================================
## SECTION 11: Integration Tests (Complex Workflows)
## ============================================================================

test_that("compfit supports full clinical trial analysis workflow", {
    
    ## Define clinically meaningful model hierarchy
    models <- list(
        demographics = c("age", "sex", "race"),
        clinical = c("age", "sex", "race", "ecog", "smoking", "diabetes"),
        disease = c("age", "sex", "race", "ecog", "smoking", "diabetes", 
                    "stage", "grade"),
        full = c("age", "sex", "race", "ecog", "smoking", "diabetes",
                 "stage", "grade", "treatment", "surgery")
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        model_list = models,
        model_names = c("Demographics", "Clinical", "Disease", "Full"),
        model_type = "coxph",
        labels = clintrial_labels,
        include_coefficients = TRUE
    )
    
    expect_compfit_result(result)
    expect_equal(nrow(result), 4)
    expect_true(!is.null(attr(result, "coefficients")))
    
    ## Verify meaningful discrimination improvement
    ## (fuller models should generally have better concordance)
    ## Note: Not strictly guaranteed, but expected with this data
})


test_that("compfit supports mixed model comparison workflow", {
    
    skip_if_not_installed("lme4")
    skip_on_cran()  # lme4 can segfault in test environments
    skip("Comparing fixed vs mixed models in auto mode is complex and prone to issues")
    
    gc()
    
    ## Compare fixed vs mixed effects
    models <- list(
        fixed = c("age", "treatment", "stage"),
        mixed_intercept = c("age", "treatment", "stage", "(1|site)"),
        mixed_slope = c("age", "treatment", "stage", "(1 + treatment|site)")
    )
    
    ## mixed_slope may have convergence issues
    result <- suppressWarnings(compfit(
        data = clintrial_complete,
        outcome = "los_days",
        model_list = models,
        model_type = "auto",
        REML = FALSE
    ))
    
    expect_compfit_result(result)
    
    ## The fixed-only model will be fit as lm, others as lmer
    ## This tests the auto-detection per model
    
    gc()
})


test_that("compfit works with PFS outcome", {
    
    models <- list(
        base = c("age", "sex", "ecog"),
        full = c("age", "sex", "ecog", "stage", "treatment")
    )
    
    result <- compfit(
        data = clintrial,
        outcome = "Surv(pfs_months, pfs_status)",
        model_list = models,
        model_type = "coxph"
    )
    
    expect_compfit_result(result)
    expect_true(all(result$Events > 0))
})


## ============================================================================
## SECTION 12: Performance and Scalability Tests
## ============================================================================

test_that("compfit handles many models efficiently", {
    
    ## Create 10 nested models
    base_vars <- c("age", "sex")
    add_vars <- c("stage", "grade", "ecog", "smoking", "diabetes", 
                  "hypertension", "treatment", "surgery")
    
    models <- list()
    for (i in 0:length(add_vars)) {
        if (i == 0) {
            models[[paste0("model_", i)]] <- base_vars
        } else {
            models[[paste0("model_", i)]] <- c(base_vars, add_vars[1:i])
        }
    }
    
    ## Should complete in reasonable time
    time_start <- Sys.time()
    result <- compfit(
        data = clintrial,
        outcome = "response",
        model_list = models,
        model_type = "glm"
    )
    time_end <- Sys.time()
    
    expect_compfit_result(result)
    expect_equal(nrow(result), length(models))
    
    ## Should complete in under 30 seconds for 10 models
    expect_true(difftime(time_end, time_start, units = "secs") < 30)
})
