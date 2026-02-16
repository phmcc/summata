#' Test Suite for International Number Formatting
#'
#' Tests for locale-aware number formatting across all summata functions:
#' desctable, survtable, fit, fullfit, uniscreen, compfit, multifit,
#' and forest plot functions. Tests exercise the public API; a small
#' section uses \code{:::} for edge cases that are difficult to trigger
#' through exported functions alone.
#'
#' @details Run with testthat::test_file("tests/testthat/test_number_format.R")

library(testthat)
library(data.table)
library(summata)

## ============================================================================
## Setup: Create test data
## ============================================================================

data(clintrial)
data(clintrial_labels)

## Small dataset for fast regression tests
set.seed(42)
test_dt <- data.table(
    age = rnorm(200, 55, 15),
    bmi = rnorm(200, 27, 5),
    sbp = rnorm(200, 130, 20),
    treatment = factor(sample(c("Drug", "Placebo"), 200, replace = TRUE)),
    sex = factor(sample(c("Male", "Female"), 200, replace = TRUE)),
    outcome = rbinom(200, 1, 0.3),
    surv_time = rexp(200, 0.02),
    surv_status = rbinom(200, 1, 0.4),
    site = factor(sample(paste0("Site", 1:5), 200, replace = TRUE))
)


## ============================================================================
## Section 1: Validation — invalid number_format rejected everywhere
## ============================================================================

test_that("desctable rejects invalid number_format", {
    expect_error(desctable(clintrial, variables = "age",
                           number_format = "french"))
    expect_error(desctable(clintrial, variables = "age",
                           number_format = 123))
})

test_that("fit rejects invalid number_format", {
    expect_error(fit(test_dt, outcome = "outcome", predictors = "age",
                     number_format = "invalid"))
})

test_that("survtable rejects invalid number_format", {
    skip_if_not_installed("survival")
    expect_error(survtable(clintrial,
                           outcome = "Surv(surv_time, surv_status)",
                           times = 12, number_format = "invalid"))
})

test_that("compfit rejects invalid number_format", {
    expect_error(compfit(test_dt, outcome = "outcome",
                         model_list = list("M1" = "age"),
                         number_format = c(",", ",")))
})

test_that("multifit rejects invalid number_format", {
    expect_error(multifit(test_dt, outcomes = "outcome",
                          predictor = "treatment",
                          number_format = c("x")))
})


## ============================================================================
## Section 2: desctable() with number_format
## ============================================================================

test_that("desctable number_format=NULL runs without error", {
    result <- desctable(clintrial, variables = "age", number_format = NULL)
    expect_s3_class(result, "data.table")
})

test_that("desctable eu continuous stats use comma decimal", {
    result <- desctable(clintrial, variables = c("age", "bmi"),
                        number_format = "eu")
    expect_s3_class(result, "data.table")

    ## Stat rows should contain comma (eu decimal mark)
    all_cells <- unlist(result[, -c("Variable", "Group"), with = FALSE])
    non_empty <- all_cells[nchar(all_cells) > 0]
    expect_true(any(grepl(",", non_empty)))
})

test_that("desctable eu categorical percentages use comma decimal", {
    result <- desctable(clintrial, variables = "sex",
                        number_format = "eu")

    all_cells <- unlist(result[, -c("Variable", "Group"), with = FALSE])
    pct_cells <- all_cells[grepl("%", all_cells)]
    if (length(pct_cells) > 0) {
        expect_true(any(grepl(",", pct_cells)))
    }
})

test_that("desctable eu p-values use comma decimal", {
    result <- desctable(clintrial, by = "treatment",
                        variables = "age", test = TRUE,
                        number_format = "eu")

    if ("p-value" %in% names(result)) {
        p_vals <- result[`p-value` != "" & `p-value` != "-", `p-value`]
        if (length(p_vals) > 0) {
            ## No bare "0.xxx" — should be "0,xxx" or "< 0,xxx"
            expect_false(any(grepl("^0\\.", p_vals)))
        }
    }
})

test_that("desctable eu N row uses dot thousands separator", {
    big_dt <- rbindlist(replicate(10, clintrial, simplify = FALSE))
    result <- desctable(big_dt, by = "treatment", variables = "age",
                        number_format = "eu")

    n_row <- result[Variable == "N"]
    all_n <- unlist(n_row[, -c("Variable", "Group"), with = FALSE])
    numeric_n <- all_n[all_n != ""]
    ## If N >= 1000, thousands separator should be dot (not comma)
    if (any(nchar(gsub("[^0-9]", "", numeric_n)) >= 4)) {
        expect_true(any(grepl("\\.", numeric_n)))
    }
})

test_that("desctable space format runs without error", {
    result <- desctable(clintrial, variables = "age",
                        number_format = "space")
    expect_s3_class(result, "data.table")
})

test_that("desctable custom format runs without error", {
    result <- desctable(clintrial, variables = "age",
                        number_format = c("'", "."))
    expect_s3_class(result, "data.table")
})

test_that("desctable global option summata.number_format works", {
    old <- getOption("summata.number_format")
    options(summata.number_format = "eu")
    on.exit(options(summata.number_format = old), add = TRUE)

    result <- desctable(clintrial, variables = "age")
    all_cells <- unlist(result[, -c("Variable", "Group"), with = FALSE])
    expect_true(any(grepl(",", all_cells)))
})

test_that("desctable US and eu produce same structure", {
    result_us <- desctable(clintrial, by = "treatment",
                           variables = c("age", "sex"),
                           number_format = "us")
    result_eu <- desctable(clintrial, by = "treatment",
                             variables = c("age", "sex"),
                             number_format = "eu")

    expect_equal(nrow(result_us), nrow(result_eu))
    expect_equal(ncol(result_us), ncol(result_eu))
    expect_equal(names(result_us), names(result_eu))
    expect_equal(result_us$Variable, result_eu$Variable)
    expect_equal(result_us$Group, result_eu$Group)
})


## ============================================================================
## Section 3: desctable() — negative values and range separators
## ============================================================================

test_that("desctable negative range uses 'to' separator", {
    dt_neg <- data.table(
        x = c(-10, -5, 0, 5, 10, 15),
        grp = factor(c("A", "A", "A", "B", "B", "B"))
    )

    result <- desctable(dt_neg, variables = "x",
                        stats_continuous = "range")

    all_cells <- unlist(result[, -c("Variable", "Group"), with = FALSE])
    range_cells <- all_cells[nchar(all_cells) > 3]
    if (length(range_cells) > 0) {
        expect_true(any(grepl(" to ", range_cells)))
    }
})

test_that("desctable eu positive range uses en-dash", {
    dt_pos <- data.table(x = c(1, 5, 10, 15, 20))

    result <- desctable(dt_pos, variables = "x",
                        stats_continuous = "range",
                        number_format = "eu")

    all_cells <- unlist(result[, -c("Variable", "Group"), with = FALSE])
    range_cells <- all_cells[nchar(all_cells) > 3]
    if (length(range_cells) > 0) {
        expect_true(any(grepl("\u2013", range_cells)))
    }
})

test_that("desctable median_iqr with negatives avoids ambiguous output", {
    dt_neg <- data.table(x = rnorm(100, -5, 3))

    result <- desctable(dt_neg, variables = "x",
                        stats_continuous = "median_iqr")

    all_cells <- unlist(result[, -c("Variable", "Group"), with = FALSE])
    iqr_cells <- all_cells[grepl("\\[", all_cells)]
    if (length(iqr_cells) > 0) {
        ## Should never have double-dash like "[-8.2--3.1]"
        expect_false(any(grepl("--", iqr_cells, fixed = TRUE)))
    }
})


## ============================================================================
## Section 4: survtable() with number_format
## ============================================================================

test_that("survtable works with eu number_format", {
    skip_if_not_installed("survival")
    result <- survtable(clintrial,
                        outcome = "Surv(os_months, os_status)",
                        times = c(12, 24), number_format = "eu")
    expect_s3_class(result, "data.table")
})

test_that("survtable works with space number_format", {
    skip_if_not_installed("survival")
    result <- survtable(clintrial,
                        outcome = "Surv(os_months, os_status)",
                        times = c(12, 24), number_format = "space")
    expect_s3_class(result, "data.table")
})


## ============================================================================
## Section 5: fit() with number_format
## ============================================================================

test_that("fit eu effect column uses comma decimal", {
    result <- fit(test_dt, outcome = "outcome",
                  predictors = c("age", "sex"),
                  number_format = "eu")
    expect_s3_class(result, "data.table")

    effect_col <- grep("OR|HR|RR|Coeff", names(result), value = TRUE)
    if (length(effect_col) > 0) {
        vals <- result[[effect_col[1]]]
        non_ref <- vals[vals != "" & vals != "reference"]
        if (length(non_ref) > 0) {
            expect_true(any(grepl(",", non_ref)))
        }
    }
})

test_that("fit eu p-values use comma decimal", {
    result <- fit(test_dt, outcome = "outcome",
                  predictors = "age", number_format = "eu")

    if ("p-value" %in% names(result)) {
        p_vals <- result[`p-value` != "" & `p-value` != "-", `p-value`]
        if (length(p_vals) > 0) {
            expect_false(any(grepl("^0\\.", p_vals)))
        }
    }
})

test_that("fit linear model never produces double-dash CIs", {
    result <- fit(test_dt, outcome = "sbp",
                  predictors = c("age", "treatment"),
                  family = "gaussian",
                  number_format = "eu")

    effect_col <- grep("Coeff", names(result), value = TRUE)
    if (length(effect_col) > 0) {
        vals <- result[[effect_col[1]]]
        non_ref <- vals[vals != "" & vals != "reference"]
        expect_false(any(grepl("--", non_ref, fixed = TRUE)))
    }
})

test_that("fit space format runs without error", {
    result <- fit(test_dt, outcome = "outcome",
                  predictors = c("age", "sex"),
                  number_format = "space")
    expect_s3_class(result, "data.table")
})

test_that("fit with pre-fitted model accepts number_format", {
    mod <- glm(outcome ~ age + sex, data = test_dt, family = binomial)
    result <- fit(test_dt, model = mod, number_format = "eu")
    expect_s3_class(result, "data.table")
})

test_that("fit US and eu produce same structure", {
    result_us <- fit(test_dt, outcome = "outcome",
                     predictors = c("age", "sex"),
                     number_format = "us")
    result_eu <- fit(test_dt, outcome = "outcome",
                       predictors = c("age", "sex"),
                       number_format = "eu")

    expect_equal(nrow(result_us), nrow(result_eu))
    expect_equal(ncol(result_us), ncol(result_eu))
})


## ============================================================================
## Section 6: uniscreen() with number_format
## ============================================================================

test_that("uniscreen eu runs without error", {
    result <- uniscreen(test_dt, outcome = "outcome",
                        predictors = c("age", "sex"),
                        number_format = "eu")
    expect_s3_class(result, "data.table")
})

test_that("uniscreen eu effect column uses comma", {
    result <- uniscreen(test_dt, outcome = "outcome",
                        predictors = c("age", "sex"),
                        number_format = "eu")

    effect_col <- grep("OR|HR|RR|Coeff", names(result), value = TRUE)
    if (length(effect_col) > 0) {
        vals <- result[[effect_col[1]]]
        non_ref <- vals[vals != "" & vals != "reference"]
        if (length(non_ref) > 0) {
            expect_true(any(grepl(",", non_ref)))
        }
    }
})


## ============================================================================
## Section 7: fullfit() with number_format
## ============================================================================

test_that("fullfit eu runs without error", {
    result <- fullfit(test_dt, outcome = "outcome",
                      predictors = c("age", "sex"),
                      method = "all", number_format = "eu")
    expect_s3_class(result, "data.table")
})

test_that("fullfit eu output contains comma decimals", {
    result <- fullfit(test_dt, outcome = "outcome",
                      predictors = c("age", "sex"),
                      method = "all", columns = "both",
                      number_format = "eu")

    all_cells <- unlist(result[, -c("Variable", "Group"), with = FALSE])
    non_empty <- all_cells[all_cells != "" & all_cells != "reference" &
                           all_cells != "-"]
    if (length(non_empty) > 0) {
        expect_true(any(grepl(",", non_empty)))
    }
})


## ============================================================================
## Section 8: compfit() with number_format
## ============================================================================

test_that("compfit eu runs without error", {
    result <- compfit(test_dt, outcome = "outcome",
                      model_list = list("M1" = "age",
                                        "M2" = c("age", "sex")),
                      number_format = "eu")
    expect_s3_class(result, "data.table")
})

test_that("compfit eu Global p uses comma decimal", {
    result <- compfit(test_dt, outcome = "outcome",
                      model_list = list("M1" = "age",
                                        "M2" = c("age", "sex")),
                      number_format = "eu")

    if ("Global p" %in% names(result)) {
        gp <- result$`Global p`
        gp_valid <- gp[!is.na(gp) & gp != "-"]
        if (length(gp_valid) > 0) {
            expect_false(any(grepl("^0\\.", gp_valid)))
        }
    }
})


## ============================================================================
## Section 9: multifit() with number_format
## ============================================================================

test_that("multifit eu runs without error", {
    result <- multifit(test_dt, outcomes = "outcome",
                       predictor = "treatment",
                       covariates = "age", columns = "adjusted",
                       number_format = "eu")
    expect_s3_class(result, "data.table")
})

test_that("multifit eu effect column uses comma", {
    result <- multifit(test_dt, outcomes = "outcome",
                       predictor = "treatment",
                       covariates = "age", columns = "adjusted",
                       number_format = "eu")

    effect_col <- grep("OR|HR|RR|Coeff", names(result), value = TRUE)
    if (length(effect_col) > 0) {
        vals <- result[[effect_col[1]]]
        non_empty <- vals[vals != "" & vals != "-"]
        if (length(non_empty) > 0) {
            expect_true(any(grepl(",", non_empty)))
        }
    }
})


## ============================================================================
## Section 10: Forest plots with number_format
## ============================================================================

test_that("glmforest eu runs without error", {
    skip_if_not_installed("ggplot2")
    mod <- glm(outcome ~ age + sex, data = test_dt, family = binomial)
    p <- glmforest(mod, data = test_dt, number_format = "eu")
    expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("lmforest eu runs without error", {
    skip_if_not_installed("ggplot2")
    mod <- lm(sbp ~ age + sex, data = test_dt)
    p <- lmforest(mod, data = test_dt, number_format = "eu")
    expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("coxforest eu runs without error", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("survival")
    mod <- survival::coxph(
        survival::Surv(surv_time, surv_status) ~ age + sex,
        data = test_dt)
    p <- coxforest(mod, data = test_dt, number_format = "eu")
    expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("autoforest passes number_format through", {
    skip_if_not_installed("ggplot2")
    mod <- glm(outcome ~ age + sex, data = test_dt, family = binomial)
    p <- autoforest(mod, data = test_dt, number_format = "eu")
    expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("uniforest eu runs without error", {
    skip_if_not_installed("ggplot2")
    uni_result <- uniscreen(test_dt, outcome = "outcome",
                            predictors = c("age", "sex"))
    p <- uniforest(uni_result, number_format = "eu")
    expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("multiforest eu runs without error", {
    skip_if_not_installed("ggplot2")
    mf_result <- multifit(test_dt, outcomes = "outcome",
                          predictor = "treatment",
                          covariates = "age", columns = "adjusted")
    p <- multiforest(mf_result, number_format = "eu")
    expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})


## ============================================================================
## Section 11: Edge cases via ::: (hard to trigger through public API)
## ============================================================================

## These internal functions have tricky logic that is difficult to exercise
## reliably through the public API with controlled inputs.

test_that("negative zero is fixed for all presets", {
    apply_dm <- summata:::apply_decimal_mark
    for (preset in c("us", "eu", "space", "none")) {
        m <- summata:::resolve_number_marks(preset)
        result <- apply_dm("-0.00", m)
        expect_false(startsWith(result, "-"),
                     info = paste("negative zero not fixed for:", preset))
    }
})

test_that("resolve_separator: negatives always produce ' to '", {
    resolve_sep <- summata:::resolve_separator
    for (preset in c("us", "eu", "space", "none")) {
        m <- summata:::resolve_number_marks(preset)
        expect_equal(resolve_sep(-5, 3, m), " to ",
                     info = paste("failed for:", preset))
        expect_equal(resolve_sep(-5, -3, m), " to ",
                     info = paste("both-negative failed for:", preset))
    }
})

test_that("resolve_separator: eu positive uses en-dash", {
    resolve_sep <- summata:::resolve_separator
    m <- summata:::resolve_number_marks("eu")
    expect_equal(resolve_sep(1, 3, m), "\u2013")
})

test_that("resolve_separator: US positive uses hyphen", {
    resolve_sep <- summata:::resolve_separator
    m <- summata:::resolve_number_marks("us")
    expect_equal(resolve_sep(1, 3, m), "-")
})

test_that("validate_number_format rejects identical custom separators", {
    validate <- summata:::validate_number_format
    expect_error(validate(c(",", ",")), "cannot be the same")
})

test_that("validate_number_format rejects wrong-length vector", {
    validate <- summata:::validate_number_format
    expect_error(validate(c(",", ".", ";")))
    expect_error(validate(c(",")))
})

test_that("format_pvalues_fit: eu threshold string uses comma", {
    fmt_p <- summata:::format_pvalues_fit
    m <- summata:::resolve_number_marks("eu")
    result <- fmt_p(c(0.042, 0.0001, NA), 3, m)
    expect_equal(result[1], "0,042")
    expect_equal(result[2], "< 0,001")
    expect_equal(result[3], "-")
})

test_that("format_pvalues_fit: NULL marks falls back to US", {
    fmt_p <- summata:::format_pvalues_fit
    result <- fmt_p(c(0.042, 0.0001, NA), 3)
    expect_equal(result[1], "0.042")
    expect_equal(result[2], "< 0.001")
    expect_equal(result[3], "-")
})

test_that("forest_ci_separator: eu positive uses en-dash", {
    fcs <- summata:::forest_ci_separator
    m <- summata:::resolve_number_marks("eu")
    expect_equal(fcs(FALSE, m), "\u2013")
})

test_that("forest_ci_separator: negatives always use comma-space", {
    fcs <- summata:::forest_ci_separator
    for (preset in c("us", "eu", "space")) {
        m <- summata:::resolve_number_marks(preset)
        expect_equal(fcs(TRUE, m), ", ",
                     info = paste("failed for:", preset))
    }
})

test_that("forest_ci_separator: NULL marks defaults to US", {
    fcs <- summata:::forest_ci_separator
    expect_equal(fcs(FALSE, NULL), "-")
    expect_equal(fcs(TRUE, NULL), ", ")
})

test_that("global option fallback works when number_format is NULL", {
    resolve <- summata:::resolve_number_marks

    old <- getOption("summata.number_format")
    options(summata.number_format = "eu")
    on.exit(options(summata.number_format = old), add = TRUE)

    m <- resolve(NULL)
    expect_equal(m$decimal.mark, ",")

    options(summata.number_format = NULL)
    m <- resolve(NULL)
    expect_equal(m$decimal.mark, ".")
})
