## The n and events columns that m2dt() returns hold group counts for the
## levels of a factor. Model-level counts are recorded separately, and
## compfit() and the fit() print method must use those rather than the
## first row of the table, which describes a single level whenever the first
## predictor is a factor.

test_that("m2dt records the model-level sample size and events", {

    df <- fresh_clintrial()

    model <- stats::glm(readmission_30d ~ sex + age, data = df, family = stats::binomial)
    raw <- m2dt(df, model, include_intercept = FALSE)
    expect_equal(attr(raw, "model_n"), stats::nobs(model))
    expect_equal(attr(raw, "model_events"), sum(model$y))
    ## The rows of sex hold its group counts, which partition the model n
    expect_equal(sum(raw$n[raw$variable == "sex"]), attr(raw, "model_n"))
    expect_true(all(raw$n[raw$variable == "sex"] < attr(raw, "model_n")))

    model <- survival::coxph(survival::Surv(os_months, os_status) ~ sex + age, data = df)
    raw <- m2dt(df, model)
    expect_equal(attr(raw, "model_n"), model$n)
    expect_equal(attr(raw, "model_events"), model$nevent)
})


test_that("compfit reports model-level counts whatever the predictor order", {

    df <- fresh_clintrial()
    df <- df[stats::complete.cases(df[, c("readmission_30d", "age", "sex", "stage")]), ]

    result <- expect_no_warning(
        suppressMessages(compfit(df, "readmission_30d",
                                 list(a = c("sex", "age"), b = c("stage", "age")))),
        message = "different numbers of observations"
    )
    events <- sum(df$readmission_30d == levels(df$readmission_30d)[2])

    expect_equal(result$N, c(nrow(df), nrow(df)))
    expect_equal(result$Events, c(events, events))
})


test_that("the fit print method reports the model-level count when counts are unavailable", {

    skip_if_not_installed("coxme")
    df <- fresh_clintrial()

    result <- suppressMessages(fit(df, "Surv(os_months, os_status)",
                                   c("sex", "age", "(1|site)"), model_type = "coxme"))
    attr(result, "analysis_counts") <- NULL
    model <- attr(result, "model")

    printed <- utils::capture.output(print(result))
    expect_true(any(grepl(paste0("Observations analyzed: ", model$n[2]), printed, fixed = TRUE)))
    expect_true(any(grepl(paste0("Events analyzed: ", model$n[1]), printed, fixed = TRUE)))
})


test_that("coxme models fitted by fit() count the observations used in fitting", {

    ## The call of a coxme model fitted by fit() names the variable "formula",
    ## from which no model variable could be read: group counts included
    ## observations with a missing covariate, and the print method stated no
    ## denominator for the analyzed sample
    skip_if_not_installed("coxme")

    df <- fresh_clintrial()
    result <- fit(data = df, outcome = "Surv(os_months, os_status)",
                  predictors = c("age", "hypertension", "diabetes"),
                  model_type = "coxme", random = "(1|site)")
    raw <- attr(result, "raw_data")
    counts <- attr(result, "analysis_counts")
    used <- stats::complete.cases(df[, c("os_months", "os_status", "age",
                                         "hypertension", "diabetes", "site")])

    expect_equal(counts$n_supplied, nrow(df))
    expect_equal(counts$n_analyzed, sum(used))
    for (v in c("hypertension", "diabetes")) {
        expect_equal(sum(raw$n[raw$variable == v]), sum(used), info = v)
    }
    expect_output(print(result),
                  sprintf("Observations analyzed: %d of %d", sum(used), nrow(df)),
                  fixed = TRUE)
})
