#' Test Suite for Forest Plot Functions
#' 
#' Comprehensive tests covering glmforest, lmforest, coxforest, autoforest,
#' uniforest, and multiforest. Tests output structure, parameters, visual
#' elements, mixed-effects models, and interaction terms.
#' 
#' @details Run with testthat::test_file("tests/testthat/test_forest.R")

library(testthat)
library(data.table)
library(ggplot2)
library(summata)

## Load survival for Surv() function
if (requireNamespace("survival", quietly = TRUE)) {
    library(survival)
}

## ============================================================================
## Setup: Create test data and helper functions
## ============================================================================

## Use package's built-in clinical trial data
data(clintrial)
data(clintrial_labels)

## Create complete-case subset for mixed-effects tests
clintrial_complete <- na.omit(clintrial)

## Create binary response for logistic tests
clintrial$response <- as.integer(clintrial$os_status == 1 & clintrial$os_months < 24)
clintrial_complete$response <- as.integer(clintrial_complete$os_status == 1 & 
                                          clintrial_complete$os_months < 24)

## ============================================================================
## Helper Functions
## ============================================================================

## Helper to check forest plot output structure
expect_forest_plot <- function(p, has_dims = TRUE) {
    expect_true(inherits(p, "ggplot"), 
                info = "Forest plot should be a ggplot object")
    
    if (has_dims) {
        dims <- attr(p, "recommended_dims")
        expect_true(!is.null(dims), 
                    info = "Forest plot should have recommended_dims attribute")
        expect_true("width" %in% names(dims), 
                    info = "recommended_dims should have width")
        expect_true("height" %in% names(dims), 
                    info = "recommended_dims should have height")
        expect_true(dims$width > 0, info = "Width should be positive")
        expect_true(dims$height > 0, info = "Height should be positive")
    }
}

## Helper to check plot layers
expect_has_layers <- function(p, min_layers = 5) {
    n_layers <- length(p$layers)
    expect_true(n_layers >= min_layers,
                info = paste("Expected at least", min_layers, "layers, got", n_layers))
}

## ============================================================================
## SECTION 1: glmforest() Basic Tests
## ============================================================================

test_that("glmforest works with basic GLM", {
    
    model <- glm(response ~ age + sex + smoking, 
                 data = clintrial, family = binomial)
    
    p <- suppressMessages(glmforest(model, data = clintrial, 
                                    title = "Test GLM Forest"))
    
    expect_forest_plot(p)
    expect_has_layers(p)
})


test_that("glmforest works with fit_result input", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "smoking"),
        model_type = "glm"
    )
    
    p <- suppressMessages(glmforest(result, title = "GLM from fit()"))
    
    expect_forest_plot(p)
})


test_that("glmforest respects digits parameter", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    p1 <- suppressMessages(glmforest(model, data = clintrial, digits = 1))
    p2 <- suppressMessages(glmforest(model, data = clintrial, digits = 4))
    
    expect_forest_plot(p1)
    expect_forest_plot(p2)
})


test_that("glmforest respects p_digits parameter", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    p1 <- suppressMessages(glmforest(model, data = clintrial, p_digits = 2))
    p2 <- suppressMessages(glmforest(model, data = clintrial, p_digits = 4))
    
    expect_forest_plot(p1)
    expect_forest_plot(p2)
})


test_that("glmforest respects conf_level parameter", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    ## Test various confidence levels
    p_90 <- suppressMessages(glmforest(model, data = clintrial, conf_level = 0.90))
    p_95 <- suppressMessages(glmforest(model, data = clintrial, conf_level = 0.95))
    p_99 <- suppressMessages(glmforest(model, data = clintrial, conf_level = 0.99))
    
    expect_forest_plot(p_90)
    expect_forest_plot(p_95)
    expect_forest_plot(p_99)
})


test_that("glmforest respects show_n parameter", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    p_with_n <- suppressMessages(glmforest(model, data = clintrial, show_n = TRUE))
    p_without_n <- suppressMessages(glmforest(model, data = clintrial, show_n = FALSE))
    
    expect_forest_plot(p_with_n)
    expect_forest_plot(p_without_n)
})


test_that("glmforest respects show_events parameter", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    p_with_events <- suppressMessages(glmforest(model, data = clintrial, show_events = TRUE))
    p_without_events <- suppressMessages(glmforest(model, data = clintrial, show_events = FALSE))
    
    expect_forest_plot(p_with_events)
    expect_forest_plot(p_without_events)
})


test_that("glmforest respects indent_groups parameter", {
    
    model <- glm(response ~ age + stage, data = clintrial, family = binomial)
    
    p_flat <- suppressMessages(glmforest(model, data = clintrial, indent_groups = FALSE))
    p_indented <- suppressMessages(glmforest(model, data = clintrial, indent_groups = TRUE))
    
    expect_forest_plot(p_flat)
    expect_forest_plot(p_indented)
})


test_that("glmforest respects condense_table parameter", {
    
    model <- glm(response ~ age + sex + smoking, data = clintrial, family = binomial)
    
    p_full <- suppressMessages(glmforest(model, data = clintrial, condense_table = FALSE))
    p_condensed <- suppressMessages(glmforest(model, data = clintrial, condense_table = TRUE))
    
    expect_forest_plot(p_full)
    expect_forest_plot(p_condensed)
})


test_that("glmforest respects zebra_stripes parameter", {
    
    model <- glm(response ~ age + sex + smoking, data = clintrial, family = binomial)
    
    p_with_stripes <- suppressMessages(glmforest(model, data = clintrial, zebra_stripes = TRUE))
    p_without_stripes <- suppressMessages(glmforest(model, data = clintrial, zebra_stripes = FALSE))
    
    expect_forest_plot(p_with_stripes)
    expect_forest_plot(p_without_stripes)
})


test_that("glmforest respects qc_footer parameter", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    p_with_footer <- suppressMessages(glmforest(model, data = clintrial, qc_footer = TRUE))
    p_without_footer <- suppressMessages(glmforest(model, data = clintrial, qc_footer = FALSE))
    
    expect_forest_plot(p_with_footer)
    expect_forest_plot(p_without_footer)
})


test_that("glmforest respects labels parameter", {
    
    model <- glm(response ~ age + sex + smoking, data = clintrial, family = binomial)
    
    p <- suppressMessages(glmforest(model, data = clintrial, 
                                    labels = clintrial_labels))
    
    expect_forest_plot(p)
})


test_that("glmforest respects color parameter", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    p_default <- suppressMessages(glmforest(model, data = clintrial))
    p_red <- suppressMessages(glmforest(model, data = clintrial, color = "#E74C3C"))
    p_blue <- suppressMessages(glmforest(model, data = clintrial, color = "blue"))
    
    expect_forest_plot(p_default)
    expect_forest_plot(p_red)
    expect_forest_plot(p_blue)
})


test_that("glmforest respects bold_variables parameter", {
    
    model <- glm(response ~ age + sex + stage, data = clintrial, family = binomial)
    
    p_plain <- suppressMessages(glmforest(model, data = clintrial, bold_variables = FALSE))
    p_bold <- suppressMessages(glmforest(model, data = clintrial, bold_variables = TRUE))
    
    expect_forest_plot(p_plain)
    expect_forest_plot(p_bold)
})


test_that("glmforest handles Poisson regression", {
    
    ## Create count outcome
    clintrial$event_count <- rpois(nrow(clintrial), lambda = 3)
    
    model <- glm(event_count ~ age + sex + smoking, 
                 data = clintrial, family = poisson)
    
    p <- suppressMessages(glmforest(model, data = clintrial, 
                                    title = "Poisson Model"))
    
    expect_forest_plot(p)
})


test_that("glmforest handles interaction terms", {
    
    model <- glm(response ~ age + sex + treatment + sex:treatment,
                 data = clintrial, family = binomial)
    
    p <- suppressMessages(glmforest(model, data = clintrial,
                                    title = "GLM with Interactions",
                                    labels = clintrial_labels))
    
    expect_forest_plot(p)
})


## ============================================================================
## SECTION 2: lmforest() Basic Tests
## ============================================================================

test_that("lmforest works with basic LM", {
    
    model <- lm(los_days ~ age + sex + smoking, data = clintrial)
    
    p <- suppressMessages(lmforest(model, data = clintrial,
                                   title = "Test LM Forest"))
    
    expect_forest_plot(p)
    expect_has_layers(p)
})


test_that("lmforest works with fit_result input", {
    
    result <- fit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex", "smoking"),
        model_type = "lm"
    )
    
    p <- suppressMessages(lmforest(result, title = "LM from fit()"))
    
    expect_forest_plot(p)
})


test_that("lmforest respects key parameters", {
    
    model <- lm(los_days ~ age + sex + stage, data = clintrial)
    
    ## Test multiple parameters at once
    p <- suppressMessages(lmforest(model, data = clintrial,
                                   digits = 3,
                                   p_digits = 4,
                                   conf_level = 0.90,
                                   show_n = TRUE,
                                   indent_groups = TRUE,
                                   zebra_stripes = TRUE,
                                   bold_variables = TRUE,
                                   labels = clintrial_labels))
    
    expect_forest_plot(p)
})


test_that("lmforest handles interaction terms", {
    
    model <- lm(los_days ~ age + sex + treatment + age:treatment,
                data = clintrial)
    
    p <- suppressMessages(lmforest(model, data = clintrial,
                                   title = "LM with Interactions"))
    
    expect_forest_plot(p)
})


## ============================================================================
## SECTION 3: coxforest() Basic Tests
## ============================================================================

test_that("coxforest works with basic coxph", {
    
    skip_if_not_installed("survival")
    
    model <- survival::coxph(Surv(os_months, os_status) ~ age + sex + smoking,
                             data = clintrial)
    
    p <- suppressMessages(coxforest(model, data = clintrial,
                                    title = "Test Cox Forest"))
    
    expect_forest_plot(p)
    expect_has_layers(p)
})


test_that("coxforest works with fit_result input", {
    
    skip_if_not_installed("survival")
    
    result <- fit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "sex", "smoking"),
        model_type = "coxph"
    )
    
    p <- suppressMessages(coxforest(result, title = "Cox from fit()"))
    
    expect_forest_plot(p)
})


test_that("coxforest respects key parameters", {
    
    skip_if_not_installed("survival")
    
    model <- survival::coxph(Surv(os_months, os_status) ~ age + sex + stage,
                             data = clintrial)
    
    p <- suppressMessages(coxforest(model, data = clintrial,
                                    digits = 3,
                                    p_digits = 4,
                                    conf_level = 0.90,
                                    show_n = TRUE,
                                    show_events = TRUE,
                                    indent_groups = TRUE,
                                    labels = clintrial_labels))
    
    expect_forest_plot(p)
})


test_that("coxforest handles interaction terms", {
    
    skip_if_not_installed("survival")
    
    model <- survival::coxph(Surv(os_months, os_status) ~ age + sex + treatment + 
                                 treatment:stage,
                             data = clintrial)
    
    p <- suppressMessages(coxforest(model, data = clintrial,
                                    title = "Cox with Interactions"))
    
    expect_forest_plot(p)
})


## ============================================================================
## SECTION 4: autoforest() Tests
## ============================================================================

test_that("autoforest correctly routes GLM", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    p <- suppressMessages(autoforest(model, data = clintrial))
    
    expect_forest_plot(p)
})


test_that("autoforest correctly routes LM", {
    
    model <- lm(los_days ~ age + sex, data = clintrial)
    
    p <- suppressMessages(autoforest(model, data = clintrial))
    
    expect_forest_plot(p)
})


test_that("autoforest correctly routes coxph", {
    
    skip_if_not_installed("survival")
    
    model <- survival::coxph(Surv(os_months, os_status) ~ age + sex,
                             data = clintrial)
    
    p <- suppressMessages(autoforest(model, data = clintrial))
    
    expect_forest_plot(p)
})


test_that("autoforest correctly routes fit_result model", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex"),
        model_type = "glm"
    )
    
    ## Extract model from fit_result and pass to autoforest
    model <- attr(result, "model")
    p <- suppressMessages(autoforest(model, data = clintrial))
    
    expect_forest_plot(p)
})


test_that("autoforest passes through parameters", {
    
    model <- glm(response ~ age + sex + stage, data = clintrial, family = binomial)
    
    p <- suppressMessages(autoforest(model, data = clintrial,
                                     title = "Custom Title",
                                     digits = 3,
                                     p_digits = 4,
                                     conf_level = 0.90,
                                     labels = clintrial_labels))
    
    expect_forest_plot(p)
})


## ============================================================================
## SECTION 5: uniforest() Tests
## ============================================================================

test_that("uniforest works with uniscreen result", {
    
    uni_result <- uniscreen(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "smoking", "stage"),
        model_type = "glm",
        parallel = FALSE
    )
    
    p <- suppressMessages(uniforest(uni_result, 
                                    title = "Univariable Screening"))
    
    expect_forest_plot(p)
})


test_that("uniforest respects key parameters", {
    
    uni_result <- uniscreen(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "smoking"),
        model_type = "glm",
        parallel = FALSE
    )
    
    p <- suppressMessages(uniforest(uni_result,
                                    digits = 3,
                                    p_digits = 4,
                                    conf_level = 0.90,
                                    show_n = TRUE,
                                    show_events = TRUE,
                                    indent_groups = TRUE,
                                    zebra_stripes = TRUE,
                                    labels = clintrial_labels))
    
    expect_forest_plot(p)
})


test_that("uniforest works with Cox uniscreen", {
    
    skip_if_not_installed("survival")
    
    uni_result <- uniscreen(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "sex", "smoking"),
        model_type = "coxph",
        parallel = FALSE
    )
    
    p <- suppressMessages(uniforest(uni_result,
                                    title = "Univariable Cox"))
    
    expect_forest_plot(p)
})


test_that("uniforest works with LM uniscreen", {
    
    uni_result <- uniscreen(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex", "smoking"),
        model_type = "lm",
        parallel = FALSE
    )
    
    p <- suppressMessages(uniforest(uni_result,
                                    title = "Univariable Linear"))
    
    expect_forest_plot(p)
})


## ============================================================================
## SECTION 6: multiforest() Tests
## ============================================================================

test_that("multiforest works with multifit result", {
    
    multi_result <- multifit(
        data = clintrial,
        outcomes = c("response", "os_status"),
        predictor = "treatment",
        covariates = c("age", "sex"),
        model_type = "glm",
        family = "binomial"
    )
    
    p <- suppressMessages(multiforest(multi_result,
                                      title = "Multi-Outcome Analysis"))
    
    expect_forest_plot(p)
})


test_that("multiforest respects key parameters", {
    
    multi_result <- multifit(
        data = clintrial,
        outcomes = c("response", "os_status"),
        predictor = "treatment",
        covariates = c("age", "sex"),
        model_type = "glm",
        family = "binomial"
    )
    
    p <- suppressMessages(multiforest(multi_result,
                                      digits = 3,
                                      p_digits = 4,
                                      conf_level = 0.90,
                                      show_n = TRUE,
                                      show_events = TRUE,
                                      zebra_stripes = TRUE,
                                      labels = clintrial_labels))
    
    expect_forest_plot(p)
})


## ============================================================================
## SECTION 7: Mixed-Effects Model Tests
## ============================================================================

test_that("glmforest works with glmer models", {
    
    skip_if_not_installed("lme4")
    skip_on_cran()
    
    gc()
    
    test_data <- clintrial_complete[1:200, ]
    
    model <- suppressWarnings(
        lme4::glmer(os_status ~ age + sex + smoking + (1|site),
                    data = test_data, family = binomial)
    )
    
    ## Suppress nobars deprecation warning from lme4
    p <- suppressMessages(suppressWarnings(
        glmforest(model, data = test_data, title = "Mixed-Effects Logistic")
    ))
    
    expect_forest_plot(p)
    
    gc()
})


test_that("lmforest works with lmer models", {
    
    skip_if_not_installed("lme4")
    skip_on_cran()
    
    gc()
    
    test_data <- clintrial_complete[1:200, ]
    
    model <- lme4::lmer(los_days ~ age + sex + smoking + (1|site),
                        data = test_data)
    
    p <- suppressMessages(lmforest(model, data = test_data,
                                   title = "Mixed-Effects Linear"))
    
    expect_forest_plot(p)
    
    gc()
})


test_that("coxforest works with coxme models", {
    
    skip_if_not_installed("coxme")
    skip_if_not_installed("survival")
    skip_on_cran()
    
    gc()
    
    test_data <- clintrial_complete[1:200, ]
    
    model <- coxme::coxme(Surv(os_months, os_status) ~ age + sex + (1|site),
                          data = test_data)
    
    p <- suppressMessages(coxforest(model, data = test_data,
                                    title = "Mixed-Effects Cox"))
    
    expect_forest_plot(p)
    
    gc()
})


## ============================================================================
## SECTION 8: Edge Cases and Error Handling
## ============================================================================

test_that("glmforest handles single predictor", {
    
    model <- glm(response ~ age, data = clintrial, family = binomial)
    
    p <- suppressMessages(glmforest(model, data = clintrial))
    
    expect_forest_plot(p)
})


test_that("glmforest handles many predictors", {
    
    model <- glm(response ~ age + sex + smoking + stage + treatment + 
                     diabetes + hypertension + surgery,
                 data = clintrial, family = binomial)
    
    p <- suppressMessages(glmforest(model, data = clintrial))
    
    expect_forest_plot(p)
})


test_that("glmforest handles all continuous predictors", {
    
    model <- glm(response ~ age + bmi + biomarker_x,
                 data = clintrial, family = binomial)
    
    p <- suppressMessages(glmforest(model, data = clintrial))
    
    expect_forest_plot(p)
})


test_that("glmforest handles all categorical predictors", {
    
    model <- glm(response ~ sex + smoking + stage + treatment,
                 data = clintrial, family = binomial)
    
    p <- suppressMessages(glmforest(model, data = clintrial))
    
    expect_forest_plot(p)
})


test_that("forest functions handle extreme confidence levels", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    ## Very narrow CI
    p_80 <- suppressMessages(glmforest(model, data = clintrial, conf_level = 0.80))
    expect_forest_plot(p_80)
    
    ## Very wide CI
    p_99 <- suppressMessages(glmforest(model, data = clintrial, conf_level = 0.99))
    expect_forest_plot(p_99)
})


test_that("forest functions handle custom ref_label", {
    
    model <- glm(response ~ age + sex + stage, data = clintrial, family = binomial)
    
    p <- suppressMessages(glmforest(model, data = clintrial,
                                    ref_label = "ref"))
    
    expect_forest_plot(p)
})


test_that("forest functions handle unit conversions", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    p_in <- suppressMessages(glmforest(model, data = clintrial, units = "in"))
    p_cm <- suppressMessages(glmforest(model, data = clintrial, units = "cm"))
    
    expect_forest_plot(p_in)
    expect_forest_plot(p_cm)
    
    ## Check that dimensions are appropriately scaled
    dims_in <- attr(p_in, "recommended_dims")
    dims_cm <- attr(p_cm, "recommended_dims")
    
    ## cm should be larger numbers (2.54x inches)
    expect_true(dims_cm$width > dims_in$width)
})


## ============================================================================
## SECTION 9: Plot Dimension Tests
## ============================================================================

test_that("forest plots provide reasonable dimension recommendations", {
    
    ## Small model
    model_small <- glm(response ~ age + sex, data = clintrial, family = binomial)
    p_small <- suppressMessages(glmforest(model_small, data = clintrial))
    dims_small <- attr(p_small, "recommended_dims")
    
    ## Larger model
    model_large <- glm(response ~ age + sex + smoking + stage + treatment + 
                           diabetes + hypertension,
                       data = clintrial, family = binomial)
    p_large <- suppressMessages(glmforest(model_large, data = clintrial))
    dims_large <- attr(p_large, "recommended_dims")
    
    ## Larger model should have greater height
    expect_true(dims_large$height > dims_small$height,
                info = "Larger model should have taller recommended height")
    
    ## Both should have reasonable widths (8-16 inches is typical)
    expect_true(dims_small$width >= 7 && dims_small$width <= 20)
    expect_true(dims_large$width >= 7 && dims_large$width <= 20)
})


test_that("forest plots respect explicit dimensions", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    p <- suppressMessages(glmforest(model, data = clintrial,
                                    plot_width = 10,
                                    plot_height = 6))
    
    expect_forest_plot(p)
})


## ============================================================================
## SECTION 10: table_width Parameter Tests
## ============================================================================

test_that("forest plots respect table_width parameter", {
    
    model <- glm(response ~ age + sex + stage, data = clintrial, family = binomial)
    
    ## Default (0.6)
    p_default <- suppressMessages(glmforest(model, data = clintrial))
    
    ## Narrow table (0.4)
    p_narrow <- suppressMessages(glmforest(model, data = clintrial, table_width = 0.4))
    
    ## Wide table (0.7)
    p_wide <- suppressMessages(glmforest(model, data = clintrial, table_width = 0.7))
    
    expect_forest_plot(p_default)
    expect_forest_plot(p_narrow)
    expect_forest_plot(p_wide)
})


## ============================================================================
## SECTION 11: Integration with fit() workflow
## ============================================================================

test_that("complete fit() to forest workflow works for GLM", {
    
    result <- fit(
        data = clintrial,
        outcome = "response",
        predictors = c("age", "sex", "smoking", "stage"),
        model_type = "glm",
        labels = clintrial_labels
    )
    
    p <- suppressMessages(glmforest(result,
                                    title = "Predictors of Response",
                                    conf_level = 0.95,
                                    indent_groups = TRUE,
                                    qc_footer = TRUE))
    
    expect_forest_plot(p)
})


test_that("complete fit() to forest workflow works for Cox", {
    
    skip_if_not_installed("survival")
    
    result <- fit(
        data = clintrial,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "sex", "treatment", "stage"),
        model_type = "coxph",
        labels = clintrial_labels
    )
    
    p <- suppressMessages(coxforest(result,
                                    title = "Survival Analysis",
                                    conf_level = 0.95,
                                    indent_groups = TRUE))
    
    expect_forest_plot(p)
})


test_that("complete fit() to forest workflow works for LM", {
    
    result <- fit(
        data = clintrial,
        outcome = "los_days",
        predictors = c("age", "sex", "surgery", "diabetes"),
        model_type = "lm",
        labels = clintrial_labels
    )
    
    p <- suppressMessages(lmforest(result,
                                   title = "Length of Stay Predictors",
                                   conf_level = 0.95))
    
    expect_forest_plot(p)
})


## ============================================================================
## SECTION 12: Font Size Parameter Tests
## ============================================================================

test_that("forest plots respect font_size parameter", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    p_small <- suppressMessages(glmforest(model, data = clintrial, font_size = 0.8))
    p_normal <- suppressMessages(glmforest(model, data = clintrial, font_size = 1.0))
    p_large <- suppressMessages(glmforest(model, data = clintrial, font_size = 1.2))
    
    expect_forest_plot(p_small)
    expect_forest_plot(p_normal)
    expect_forest_plot(p_large)
})


test_that("forest plots respect annot_size and header_size", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    p <- suppressMessages(glmforest(model, data = clintrial,
                                    annot_size = 4.0,
                                    header_size = 6.0,
                                    title_size = 24))
    
    expect_forest_plot(p)
})


## ============================================================================
## SECTION 13: exponentiate Parameter Tests
## ============================================================================

test_that("glmforest respects exponentiate = FALSE", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    p_exp <- suppressMessages(glmforest(model, data = clintrial, exponentiate = TRUE))
    p_raw <- suppressMessages(glmforest(model, data = clintrial, exponentiate = FALSE))
    
    expect_forest_plot(p_exp)
    expect_forest_plot(p_raw)
})


## ============================================================================
## SECTION 14: Consistency Tests
## ============================================================================

test_that("autoforest produces same result as direct function call", {
    
    model <- glm(response ~ age + sex, data = clintrial, family = binomial)
    
    p_direct <- suppressMessages(glmforest(model, data = clintrial, 
                                           title = "Test",
                                           conf_level = 0.90))
    p_auto <- suppressMessages(autoforest(model, data = clintrial,
                                          title = "Test",
                                          conf_level = 0.90))
    
    ## Both should be valid forest plots
    expect_forest_plot(p_direct)
    expect_forest_plot(p_auto)
    
    ## Both should have similar dimensions
    dims_direct <- attr(p_direct, "recommended_dims")
    dims_auto <- attr(p_auto, "recommended_dims")
    
    expect_equal(dims_direct$width, dims_auto$width, tolerance = 0.1)
    expect_equal(dims_direct$height, dims_auto$height, tolerance = 0.1)
})
