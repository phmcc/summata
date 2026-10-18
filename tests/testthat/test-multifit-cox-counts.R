## For Cox models, stats::nobs() returns the number of events, and the first
## element of the n component of a coxme fit is also the number of events.
## multifit() must report observations in its n column and events in its
## Events column.

test_that("multifit reports observations and events for Cox models", {

    df <- fresh_clintrial()
    used <- stats::complete.cases(df[, c("os_months", "os_status", "age", "sex")])

    result <- multifit(df, outcomes = "Surv(os_months, os_status)", predictor = "age",
                       covariates = "sex", model_type = "coxph", parallel = FALSE)
    raw <- attr(result, "raw_data")
    expect_equal(raw$n, sum(used))
    expect_equal(raw$events, sum(df$os_status[used]))

    skip_if_not_installed("coxme")
    used <- stats::complete.cases(df[, c("os_months", "os_status", "age", "sex", "site")])
    result <- multifit(df, outcomes = "Surv(os_months, os_status)", predictor = "age",
                       covariates = "sex", random = "(1|site)", model_type = "coxme",
                       parallel = FALSE)
    raw <- attr(result, "raw_data")
    expect_equal(raw$n, sum(used))
    expect_equal(raw$events, sum(df$os_status[used]))
})
