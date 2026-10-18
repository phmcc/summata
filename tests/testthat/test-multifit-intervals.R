## multifit() reports profile likelihood intervals for the predictor's terms.
## Only those terms are profiled; the intervals must equal those of confint()
## on the full model, including a predictor with a single coefficient, for
## which confint() returns a vector rather than a matrix.

test_that("multifit intervals match confint() for the predictor terms", {

    withr::local_options(summata.conf_method = NULL)
    df <- fresh_clintrial()

    for (predictor in c("smoking", "age")) {
        result <- multifit(df, outcomes = "readmission_30d", predictor = predictor,
                           covariates = c("sex", "stage"), keep_models = TRUE,
                           parallel = FALSE)
        raw <- attr(result, "raw_data")
        model <- attr(result, "models")[["readmission_30d"]]$adjusted
        terms <- grep(paste0("^", predictor), names(stats::coef(model)), value = TRUE)
        profile <- suppressMessages(stats::confint(model))[terms, , drop = FALSE]

        expect_equal(raw$ci_lower, unname(exp(profile[, 1])))
        expect_equal(raw$ci_upper, unname(exp(profile[, 2])))
    }
})
