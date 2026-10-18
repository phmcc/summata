## When glmforest() and lmforest() receive a model rather than a fit() result,
## they compute its intervals as m2dt() does for the regression tables. Each
## interval must belong to the coefficient in its row, also when a coefficient
## is aliased, and the method must follow the summata.conf_method option, use
## Wald intervals for quasi-likelihood families, and fall back to Wald
## intervals when profiling fails.

estimated_rows <- function(p) {
    tbl <- attr(p, "table_data")
    tbl[!is.na(tbl$conf_low)]
}


test_that("forest intervals stay with their coefficients when one is aliased", {

    withr::local_options(summata.conf_method = NULL)
    df <- fresh_clintrial()
    df$age_months <- df$age * 12

    models <- list(
        lm = stats::lm(bmi ~ age + age_months + sex, data = df),
        glm = stats::glm(os_status ~ age + age_months + sex, data = df,
                         family = stats::binomial)
    )
    plotters <- list(lm = lmforest, glm = glmforest)

    for (type in names(models)) {
        model <- models[[type]]
        expect_true(is.na(stats::coef(model)[["age_months"]]))

        p <- expect_no_warning(suppressMessages(plotters[[type]](model, data = df)))
        tbl <- estimated_rows(p)
        ci <- suppressMessages(stats::confint(model, parm = c("age", "sexMale")))

        expect_equal(tbl$var, c("age", "sex"), info = type)
        expect_equal(tbl$conf_low, unname(ci[, 1]), info = type)
        expect_equal(tbl$conf_high, unname(ci[, 2]), info = type)
    }
})


test_that("forest intervals follow summata.conf_method", {

    df <- fresh_clintrial()
    model <- stats::glm(os_status ~ age + sex, data = df, family = stats::binomial)
    wald <- stats::confint.default(model, parm = c("age", "sexMale"))

    tbl <- withr::with_options(list(summata.conf_method = "wald"),
                               estimated_rows(suppressMessages(glmforest(model, data = df))))
    expect_equal(tbl$conf_low, unname(wald[, 1]))
    expect_equal(tbl$conf_high, unname(wald[, 2]))

    tbl <- withr::with_options(list(summata.conf_method = "profile"),
                               estimated_rows(suppressMessages(glmforest(model, data = df))))
    profile <- suppressMessages(stats::confint(model, parm = c("age", "sexMale")))
    expect_equal(tbl$conf_low, unname(profile[, 1]))
    expect_equal(tbl$conf_high, unname(profile[, 2]))
})


test_that("forest intervals of quasi-likelihood models are Wald intervals, as in tables", {

    df <- fresh_clintrial()
    model <- stats::glm(os_status ~ age + sex, data = df, family = stats::quasibinomial)
    tbl <- estimated_rows(suppressMessages(glmforest(model, data = df)))

    wald <- stats::confint.default(model, parm = c("age", "sexMale"))
    expect_equal(tbl$conf_low, unname(wald[, 1]))
    expect_equal(tbl$conf_high, unname(wald[, 2]))

    raw <- m2dt(df, model, include_intercept = FALSE)
    raw <- raw[!is.na(raw$ci_lower)]
    expect_equal(exp(tbl$conf_low), raw$ci_lower)
    expect_equal(exp(tbl$conf_high), raw$ci_upper)
})


test_that("forest intervals fall back to Wald intervals when profiling fails", {

    ## x separates the outcome completely, so profiling fails
    d <- data.frame(y = rep(c(0, 1), each = 20), x = c(1:20, 31:50),
                    g = factor(rep(c("a", "b"), 20)))
    model <- suppressWarnings(stats::glm(y ~ x + g, data = d, family = stats::binomial))
    profile <- tryCatch(suppressWarnings(suppressMessages(stats::confint(model))),
                        error = function(e) e)
    skip_if_not(inherits(profile, "error"), "profiling succeeds for this model")

    p <- suppressWarnings(suppressMessages(glmforest(model, data = d)))
    expect_s3_class(p, "ggplot")
    tbl <- estimated_rows(p)
    wald <- stats::confint.default(model, parm = c("x", "gb"))
    expect_equal(tbl$conf_low, unname(wald[, 1]))
    expect_equal(tbl$conf_high, unname(wald[, 2]))
})
