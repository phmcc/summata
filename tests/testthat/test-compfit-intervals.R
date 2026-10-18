## compfit() displays confidence intervals only in the coefficient table
## requested with include_coefficients = TRUE. Without it, the comparison must
## not depend on the interval method, and the returned models must not carry
## intervals computed for speed rather than for display.

compfit_data <- function() {
    df <- fresh_clintrial()
    df[stats::complete.cases(df[, c("readmission_30d", "age", "sex", "stage")]), ]
}

compfit_models <- list(a = c("age", "sex"), b = c("age", "sex", "stage"))


test_that("the comparison does not depend on whether coefficients are requested", {

    df <- compfit_data()

    without <- suppressMessages(compfit(df, "readmission_30d", compfit_models))
    with <- suppressMessages(compfit(df, "readmission_30d", compfit_models,
                                     include_coefficients = TRUE))

    plain <- function(x) {
        x <- as.data.frame(x)
        attributes(x) <- attributes(x)[c("names", "row.names", "class")]
        x
    }
    expect_identical(plain(without), plain(with))
})


test_that("requested coefficients carry profile likelihood intervals", {

    withr::local_options(summata.conf_method = NULL)
    df <- compfit_data()

    result <- suppressMessages(compfit(df, "readmission_30d", compfit_models,
                                       include_coefficients = TRUE))
    model <- attr(result, "models")$b
    expect_false(is.null(attr(model, "cached_confint")))

    ## Intervals are computed for the predictors, not the intercept
    profile <- suppressMessages(stats::confint(model))[-1, ]
    expect_equal(attr(model, "cached_confint")[rownames(profile), ], profile)
})


test_that("models returned without coefficients carry no cached intervals", {

    withr::local_options(summata.conf_method = NULL)
    df <- compfit_data()

    result <- suppressMessages(compfit(df, "readmission_30d", compfit_models))
    for (model in attr(result, "models")) {
        expect_null(attr(model, "cached_confint"))
    }

    ## The session option is restored
    expect_null(getOption("summata.conf_method"))
})
