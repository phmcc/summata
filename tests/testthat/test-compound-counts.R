## Compound tables place unadjusted and adjusted estimates side by side under a
## single sample size column. That column must describe the observations behind
## the adjusted estimate, since the adjusted model is fitted to complete cases
## across every covariate and is therefore the smaller sample. Reporting the
## unadjusted sample alongside an adjusted estimate overstates it.

test_that("multifit both-column tables report the adjusted sample size", {

    ## Missingness confined to a covariate, so the adjusted model is fitted to
    ## fewer observations than the unadjusted one.
    df <- fresh_clintrial()
    df$bmi[seq_len(120)] <- NA
    dt <- data.table::as.data.table(df)

    result <- multifit(
        data = dt,
        outcomes = "readmission_30d",
        predictor = "sex",
        covariates = c("age", "bmi"),
        columns = "both",
        model_type = "glm",
        family = "binomial"
    )

    raw <- attr(result, "raw_data")

    ## multifit reports one row per non-reference level, so the count for the
    ## single row is the size of that level within the adjusted sample.
    adj_model <- stats::glm(readmission_30d ~ sex + age + bmi,
                            data = dt, family = stats::binomial)
    adj_frame <- stats::model.frame(adj_model)
    expected <- sum(adj_frame$sex == "Male")

    expect_equal(as.numeric(raw[["n"]][1]), as.numeric(expected))

    ## The unadjusted sample is strictly larger here, so the two are
    ## distinguishable and the assertion above is meaningful.
    unadj_model <- stats::glm(readmission_30d ~ sex,
                              data = dt, family = stats::binomial)
    expect_gt(stats::nobs(unadj_model), stats::nobs(adj_model))
})
