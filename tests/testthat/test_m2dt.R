#' Test Suite for m2dt()
#' 
#' Comprehensive tests covering model-to-data.table conversion for
#' GLM, LM, Cox, mixed effects models, with QC stats, reference rows,
#' and various output options.
#' 
#' @details Run with testthat::test_file("tests/testthat/test-m2dt.R")

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
clintrial$response <- as.integer(clintrial$os_status == 1 & clintrial$os_months < 24)
clintrial_complete$response <- as.integer(clintrial_complete$os_status == 1 & clintrial_complete$os_months < 24)

## Helper function to check m2dt result structure
expect_m2dt_result <- function(result) {
    expect_s3_class(result, "data.table")
    
    ## Required columns
    expect_true("model_scope" %in% names(result))
    expect_true("model_type" %in% names(result))
    expect_true("term" %in% names(result))
    expect_true("coefficient" %in% names(result))
    expect_true("se" %in% names(result))
    expect_true("p_value" %in% names(result))
    
    ## Required attributes
    expect_true(!is.null(attr(result, "model_class")))
    expect_true(!is.null(attr(result, "formula_str")))
}

## Helper to check for effect column
expect_effect_column <- function(result, effect_type) {
    expect_true(effect_type %in% names(result),
                info = paste("Expected", effect_type, "column, found:", 
                             paste(names(result), collapse = ", ")))
}

## Helper to check for QC stats
expect_qc_stats <- function(result, stats) {
    for (stat in stats) {
        expect_true(stat %in% names(result),
                    info = paste("Expected QC stat", stat, "not found"))
    }
}


## ============================================================================
## SECTION 1: Basic GLM - Logistic Regression
## ============================================================================

test_that("m2dt works with logistic regression", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_m2dt_result(result)
    expect_effect_column(result, "OR")
    expect_equal(attr(result, "model_class"), "glm")
})


test_that("m2dt logistic shows correct model type", {
    
    model <- glm(response ~ age, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_equal(unique(result$model_type), "Logistic")
})


test_that("m2dt logistic detects univariable vs multivariable", {
    
    ## Univariable
    model_uni <- glm(response ~ age, data = clintrial, family = binomial)
    result_uni <- m2dt(data = clintrial, model = model_uni)
    expect_equal(unique(result_uni$model_scope), "Univariable")
    
    ## Multivariable
    model_multi <- glm(response ~ age + sex + stage, data = clintrial, family = binomial)
    result_multi <- m2dt(data = clintrial, model = model_multi)
    expect_equal(unique(result_multi$model_scope), "Multivariable")
})


test_that("m2dt logistic includes events count", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_true("events" %in% names(result))
    expect_true(all(!is.na(result$events)))
})


## ============================================================================
## SECTION 2: Linear Model
## ============================================================================

test_that("m2dt works with linear model", {
    
    model <- lm(los_days ~ age + sex, data = clintrial)
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_m2dt_result(result)
    expect_effect_column(result, "Coefficient")
    expect_equal(attr(result, "model_class"), "lm")
})


test_that("m2dt linear shows correct model type", {
    
    model <- lm(los_days ~ age, data = clintrial)
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_equal(unique(result$model_type), "Linear")
})


test_that("m2dt linear does not exponentiate", {
    
    model <- lm(los_days ~ age, data = clintrial)
    
    result <- m2dt(data = clintrial, model = model)
    
    ## Coefficient should equal coefficient (not exponentiated)
    expect_equal(result$Coefficient, result$coefficient)
})


## ============================================================================
## SECTION 3: Poisson Regression
## ============================================================================

test_that("m2dt works with Poisson regression", {
    
    ## Create integer count data
    test_data <- data.table::as.data.table(clintrial)
    test_data[, count_outcome := as.integer(round(los_days))]
    
    model <- suppressWarnings(
        glm(count_outcome ~ age + sex, data = test_data, family = poisson)
    )
    
    result <- m2dt(data = test_data, model = model)
    
    expect_m2dt_result(result)
    expect_effect_column(result, "RR")
})


test_that("m2dt Poisson shows correct model type", {
    
    test_data <- data.table::as.data.table(clintrial)
    test_data[, count_outcome := as.integer(round(los_days))]
    
    model <- suppressWarnings(
        glm(count_outcome ~ age, data = test_data, family = poisson)
    )
    
    result <- m2dt(data = test_data, model = model)
    
    expect_equal(unique(result$model_type), "Poisson")
})


## ============================================================================
## SECTION 4: Cox Proportional Hazards
## ============================================================================

test_that("m2dt works with Cox regression", {
    
    skip_if_not_installed("survival")
    
    model <- survival::coxph(
                           survival::Surv(os_months, os_status) ~ age + sex + stage,
                           data = clintrial
                       )
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_m2dt_result(result)
    expect_effect_column(result, "HR")
    expect_equal(attr(result, "model_class"), "coxph")
})


test_that("m2dt Cox shows correct model type", {
    
    skip_if_not_installed("survival")
    
    model <- survival::coxph(
                           survival::Surv(os_months, os_status) ~ age,
                           data = clintrial
                       )
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_equal(unique(result$model_type), "Cox PH")
})


test_that("m2dt Cox includes events count", {
    
    skip_if_not_installed("survival")
    
    model <- survival::coxph(
                           survival::Surv(os_months, os_status) ~ age + sex,
                           data = clintrial
                       )
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_true("events" %in% names(result))
    expect_true(all(!is.na(result$events)))
})


## ============================================================================
## SECTION 5: Conditional Logistic Regression
## ============================================================================

test_that("m2dt works with conditional logistic regression", {
    
    skip_if_not_installed("survival")
    
    ## Create strata variable
    test_data <- data.table::as.data.table(clintrial)
    test_data[, strata_var := rep(1:100, length.out = .N)]
    
    model <- survival::clogit(
                           response ~ age + sex + strata(strata_var),
                           data = test_data
                       )
    
    result <- m2dt(data = test_data, model = model)
    
    expect_m2dt_result(result)
    expect_effect_column(result, "HR")  # clogit uses HR
    expect_equal(attr(result, "model_class"), "clogit")
})


## ============================================================================
## SECTION 6: Mixed Effects Models - lmer
## ============================================================================

test_that("m2dt works with lmer", {
    
    skip_if_not_installed("lme4")
    
    model <- lme4::lmer(los_days ~ age + sex + (1|site), data = clintrial)
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_m2dt_result(result)
    expect_effect_column(result, "Coefficient")
    expect_true(attr(result, "model_class") %in% c("lmerMod", "lmer"))
})


test_that("m2dt lmer shows correct model type", {
    
    skip_if_not_installed("lme4")
    
    model <- lme4::lmer(los_days ~ age + (1|site), data = clintrial)
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_true(grepl("Mixed|lmer", unique(result$model_type), ignore.case = TRUE))
})


## ============================================================================
## SECTION 7: Mixed Effects Models - glmer
## ============================================================================

test_that("m2dt works with glmer (binomial)", {
    
    skip_if_not_installed("lme4")
    
    model <- suppressWarnings(
        lme4::glmer(response ~ age + sex + (1|site), 
                    data = clintrial, family = binomial)
    )
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_m2dt_result(result)
    expect_effect_column(result, "OR")
    expect_true(attr(result, "model_class") %in% c("glmerMod", "glmer"))
})


## ============================================================================
## SECTION 8: Mixed Effects Cox - coxme
## ============================================================================

test_that("m2dt works with coxme", {
    
    skip_if_not_installed("coxme")
    skip_if_not_installed("survival")
    
    model <- coxme::coxme(
                        survival::Surv(os_months, os_status) ~ age + sex + (1|site),
                        data = clintrial
                    )
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_m2dt_result(result)
    expect_effect_column(result, "HR")
    expect_equal(attr(result, "model_class"), "coxme")
})


## ============================================================================
## SECTION 9: QC Statistics - keep_qc_stats Parameter
## ============================================================================

test_that("m2dt keep_qc_stats=TRUE includes QC statistics for GLM", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, keep_qc_stats = TRUE)
    
    expect_qc_stats(result, c("AIC", "BIC"))
})


test_that("m2dt keep_qc_stats=FALSE excludes QC statistics", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, keep_qc_stats = FALSE)
    
    expect_false("AIC" %in% names(result))
    expect_false("BIC" %in% names(result))
})


test_that("m2dt includes QC stats for linear model", {
    
    model <- lm(los_days ~ age + sex, data = clintrial)
    
    result <- m2dt(data = clintrial, model = model, keep_qc_stats = TRUE)
    
    expect_qc_stats(result, c("R2", "adj_R2", "AIC", "BIC"))
})


test_that("m2dt includes QC stats for Cox model", {
    
    skip_if_not_installed("survival")
    
    model <- survival::coxph(
                           survival::Surv(os_months, os_status) ~ age + sex,
                           data = clintrial
                       )
    
    result <- m2dt(data = clintrial, model = model, keep_qc_stats = TRUE)
    
    expect_qc_stats(result, c("concordance", "likelihood_ratio_test"))
})


## ============================================================================
## SECTION 10: Intercept Handling
## ============================================================================

test_that("m2dt include_intercept=TRUE includes intercept", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, include_intercept = TRUE)
    
    expect_true("(Intercept)" %in% result$term)
})


test_that("m2dt include_intercept=FALSE excludes intercept", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, include_intercept = FALSE)
    
    expect_false("(Intercept)" %in% result$term)
})


## ============================================================================
## SECTION 11: Reference Rows for Factors
## ============================================================================

test_that("m2dt reference_rows=TRUE adds reference rows", {
    
    model <- glm(response ~ age + stage, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, reference_rows = TRUE)
    
    ## Should have a reference row
    expect_true("reference" %in% names(result))
    expect_true(any(result$reference == "reference"))
})


test_that("m2dt reference_rows=FALSE excludes reference rows", {
    
    model <- glm(response ~ age + stage, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, reference_rows = FALSE)
    
    ## Result should be valid
    expect_m2dt_result(result)
    
    ## Should NOT have reference rows with "reference" label
    if ("reference" %in% names(result)) {
        expect_true(all(result$reference == "" | is.na(result$reference)))
    } else {
        ## If no reference column at all, that's also correct
        expect_true(TRUE)
    }
})


test_that("m2dt reference rows have correct values", {
    
    model <- glm(response ~ stage, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, reference_rows = TRUE)
    
    ## Reference row should have OR = 1, coefficient = 0
    ref_rows <- result[reference == "reference"]
    
    if (nrow(ref_rows) > 0) {
        expect_true(all(ref_rows$OR == 1))
        expect_true(all(ref_rows$coefficient == 0))
    }
})


test_that("m2dt reference_label customizes reference label", {
    
    model <- glm(response ~ stage, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, 
                   reference_rows = TRUE,
                   reference_label = "ref_level")
    
    expect_true(any(result$reference == "ref_level"))
})


## ============================================================================
## SECTION 12: Term Exclusion
## ============================================================================

test_that("m2dt terms_to_exclude removes specified terms", {
    
    model <- glm(response ~ age + sex + stage, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, 
                   terms_to_exclude = c("sexMale"))
    
    expect_false("sexMale" %in% result$term)
})


test_that("m2dt terms_to_exclude with multiple terms", {
    
    model <- glm(response ~ age + sex + stage, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, 
                   terms_to_exclude = c("age", "sexMale"))
    
    expect_false("age" %in% result$term)
    expect_false("sexMale" %in% result$term)
})


## ============================================================================
## SECTION 13: Confidence Level
## ============================================================================

test_that("m2dt respects conf_level parameter", {
    
    model <- glm(response ~ age, data = clintrial, family = binomial)
    
    result_95 <- m2dt(data = clintrial, model = model, conf_level = 0.95)
    result_90 <- m2dt(data = clintrial, model = model, conf_level = 0.90)
    
    ## 90% CI should be narrower than 95% CI
    ## Column names are ci_lower/ci_upper (lowercase)
    ci_width_95 <- result_95$ci_upper[1] - result_95$ci_lower[1]
    ci_width_90 <- result_90$ci_upper[1] - result_90$ci_lower[1]
    
    expect_true(ci_width_90 < ci_width_95)
})


## ============================================================================
## SECTION 14: Significance Indicators
## ============================================================================

test_that("m2dt adds significance indicators", {
    
    model <- glm(response ~ age + stage, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_true("sig" %in% names(result))
    expect_true("sig_binary" %in% names(result))
})


test_that("m2dt sig_binary is correct", {
    
    model <- glm(response ~ age + stage, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model)
    
    ## sig_binary should be TRUE where p < 0.05
    result_check <- result[!is.na(p_value)]
    expect_equal(result_check$sig_binary, result_check$p_value < 0.05)
})


test_that("m2dt sig markers are correct", {
    
    model <- glm(response ~ age + stage, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model)
    
    ## Check significance markers
    for (i in seq_len(nrow(result))) {
        p <- result$p_value[i]
        sig <- result$sig[i]
        
        if (is.na(p)) {
            expect_equal(sig, "")
        } else if (p < 0.001) {
            expect_equal(sig, "***")
        } else if (p < 0.01) {
            expect_equal(sig, "**")
        } else if (p < 0.05) {
            expect_equal(sig, "*")
        } else if (p < 0.1) {
            expect_equal(sig, ".")
        } else {
            expect_equal(sig, "")
        }
    }
})


## ============================================================================
## SECTION 15: Variable and Group Parsing
## ============================================================================

test_that("m2dt parses factor variables correctly", {
    
    model <- glm(response ~ stage, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, 
                   include_intercept = FALSE,
                   reference_rows = TRUE)
    
    ## All stage terms should have variable = "stage"
    stage_terms <- result[variable == "stage"]
    expect_true(nrow(stage_terms) > 1)  # Multiple levels
    
    ## Group should contain level names
    expect_true(all(stage_terms$group %in% c("I", "II", "III", "IV")))
})


test_that("m2dt parses continuous variables correctly", {
    
    model <- glm(response ~ age + bmi, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, include_intercept = FALSE)
    
    ## Continuous variables should have variable = term name
    age_row <- result[variable == "age"]
    expect_equal(nrow(age_row), 1)
    expect_equal(age_row$group, "")
})


## ============================================================================
## SECTION 16: n and Events Columns
## ============================================================================

test_that("m2dt includes n column", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_true("n" %in% names(result))
    expect_true(all(!is.na(result$n)))
})


test_that("m2dt includes n_group for factor variables", {
    
    model <- glm(response ~ stage, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, reference_rows = TRUE)
    
    expect_true("n_group" %in% names(result))
})


test_that("m2dt events_group for factor variables in survival model", {
    
    skip_if_not_installed("survival")
    
    model <- survival::coxph(
                           survival::Surv(os_months, os_status) ~ stage,
                           data = clintrial
                       )
    
    result <- m2dt(data = clintrial, model = model, reference_rows = TRUE)
    
    expect_true("events_group" %in% names(result))
})


## ============================================================================
## SECTION 17: Edge Cases
## ============================================================================

test_that("m2dt handles single predictor", {
    
    model <- glm(response ~ age, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_m2dt_result(result)
})


test_that("m2dt handles many predictors", {
    
    model <- glm(response ~ age + sex + smoking + diabetes + hypertension + stage,
                 data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_m2dt_result(result)
})


test_that("m2dt handles interaction terms", {
    
    model <- glm(response ~ age + sex + age:sex, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_m2dt_result(result)
    
    ## Should have interaction term
    expect_true(any(grepl(":", result$term)))
})


test_that("m2dt handles factor-factor interaction", {
    
    model <- glm(response ~ sex + stage + sex:stage, data = clintrial, family = binomial)
    
    ## Suppress Hoslem test warning about bins
    result <- suppressWarnings(m2dt(data = clintrial, model = model))
    
    expect_m2dt_result(result)
})


test_that("m2dt requires data parameter", {
    
    model <- glm(response ~ age, data = clintrial, family = binomial)
    
    expect_error(
        m2dt(model = model),
        regexp = "data parameter is required"
    )
})


test_that("m2dt works with data.frame input", {
    
    df <- as.data.frame(clintrial)
    model <- glm(response ~ age + sex, data = df, family = binomial)
    
    result <- m2dt(data = df, model = model)
    
    expect_m2dt_result(result)
})


## ============================================================================
## SECTION 18: Attributes
## ============================================================================

test_that("m2dt sets model_class attribute", {
    
    model <- glm(response ~ age, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_equal(attr(result, "model_class"), "glm")
})


test_that("m2dt sets formula_str attribute", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model)
    
    formula_str <- attr(result, "formula_str")
    expect_true(grepl("response", formula_str))
    expect_true(grepl("age", formula_str))
})


test_that("m2dt sets family attributes for GLM", {
    
    model <- glm(response ~ age, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_equal(attr(result, "model_family"), "binomial")
    expect_equal(attr(result, "model_link"), "logit")
})


## ============================================================================
## SECTION 19: Coefficient and CI Columns
## ============================================================================

test_that("m2dt includes both raw and exponentiated coefficients", {
    
    model <- glm(response ~ age, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model)
    
    ## Should have raw coefficients
    expect_true("coefficient" %in% names(result))
    expect_true("coef" %in% names(result))
    expect_true("coef_lower" %in% names(result))
    expect_true("coef_upper" %in% names(result))
    
    ## Should have exponentiated coefficients
    expect_true("exp_coef" %in% names(result))
    expect_true("exp_lower" %in% names(result))
    expect_true("exp_upper" %in% names(result))
})


test_that("m2dt OR equals exp(coefficient) for logistic", {
    
    model <- glm(response ~ age, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, include_intercept = FALSE)
    
    expect_equal(result$OR, exp(result$coefficient), tolerance = 1e-10)
})


test_that("m2dt HR equals exp(coefficient) for Cox", {
    
    skip_if_not_installed("survival")
    
    model <- survival::coxph(
                           survival::Surv(os_months, os_status) ~ age,
                           data = clintrial
                       )
    
    result <- m2dt(data = clintrial, model = model)
    
    expect_equal(result$HR, exp(result$coefficient), tolerance = 1e-10)
})


## ============================================================================
## SECTION 20: Consistency Tests
## ============================================================================

test_that("m2dt produces consistent results across runs", {
    
    model <- glm(response ~ age + sex + stage, data = clintrial, family = binomial)
    
    result1 <- m2dt(data = clintrial, model = model)
    result2 <- m2dt(data = clintrial, model = model)
    
    expect_equal(result1$coefficient, result2$coefficient)
    expect_equal(result1$p_value, result2$p_value)
})


test_that("m2dt coefficients match summary(model)", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, include_intercept = TRUE)
    
    model_coefs <- coef(model)
    
    for (term_name in names(model_coefs)) {
        result_coef <- result[term == term_name]$coefficient
        expected_coef <- unname(model_coefs[term_name])
        expect_equal(result_coef, expected_coef, 
                     tolerance = 1e-10,
                     info = paste("Coefficient mismatch for", term_name))
    }
})


test_that("m2dt p-values match summary(model)", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    result <- m2dt(data = clintrial, model = model, include_intercept = TRUE)
    
    model_summary <- summary(model)$coefficients
    
    for (term_name in rownames(model_summary)) {
        result_p <- result[term == term_name]$p_value
        summary_p <- unname(model_summary[term_name, "Pr(>|z|)"])
        expect_equal(result_p, summary_p, 
                     tolerance = 1e-10,
                     info = paste("P-value mismatch for", term_name))
    }
})
