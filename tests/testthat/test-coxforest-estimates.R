## coxforest() draws the fixed effects of a coxph() or coxme() model. Its
## intervals and p-values must come from the variances of those fixed
## effects, and its estimates must be matched to their rows when variable
## names are not syntactic.

test_that("coxforest uses the fixed-effect variances of a coxme model", {

    skip_if_not_installed("coxme")
    df <- fresh_clintrial()

    model <- coxme::coxme(survival::Surv(os_months, os_status) ~ age + sex + (1 | site),
                          data = df)
    p <- suppressMessages(suppressWarnings(coxforest(model, data = df)))
    tbl <- attr(p, "table_data")
    tbl <- tbl[!is.na(tbl$p_value)]

    beta <- coxme::fixef(model)
    se <- sqrt(diag(as.matrix(stats::vcov(model))))
    z <- stats::qnorm(0.975)

    expect_equal(tbl$estimate, unname(beta))
    expect_equal(tbl$conf_low, unname(beta - z * se))
    expect_equal(tbl$conf_high, unname(beta + z * se))
    expect_equal(tbl$p_value, unname(2 * stats::pnorm(abs(beta / se), lower.tail = FALSE)))

    ## The same intervals as the regression table for the same model
    raw <- attr(fit(df, "Surv(os_months, os_status)", c("age", "sex", "(1|site)"),
                    model_type = "coxme"), "raw_data")
    raw <- raw[!is.na(raw$p_value)]
    expect_equal(exp(tbl$conf_low), raw$ci_lower, tolerance = 1e-6)
    expect_equal(exp(tbl$conf_high), raw$ci_upper, tolerance = 1e-6)
})


test_that("coxforest matches coefficients of non-syntactic variable names", {

    df <- fresh_clintrial()
    df$`age years` <- df$age
    df$`trial arm` <- df$treatment

    model <- survival::coxph(survival::Surv(os_months, os_status) ~ `age years` + `trial arm`,
                             data = df)
    p <- suppressMessages(suppressWarnings(coxforest(model, data = df)))
    tbl <- attr(p, "table_data")

    expect_equal(tbl$estimate[tbl$var == "age years"], unname(stats::coef(model)[1]))
    expect_equal(tbl$estimate[tbl$var == "trial arm" & tbl$level == "Drug B"],
                 unname(stats::coef(model)[3]))
})
