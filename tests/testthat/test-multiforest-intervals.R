## multiforest() draws the estimates and intervals reported by multifit():
## the profile likelihood intervals of generalized linear models and the
## exact t intervals of linear models, in either column of a two-column
## result. A different confidence level than the one used by multifit()
## yields Wald intervals at the requested level.

forest_rows <- function(p) {
    tbl <- attr(p, "table_data")
    tbl[!is.na(tbl$estimate) & tbl$level != "-"]
}


test_that("multiforest draws the intervals of a logistic multifit table", {

    df <- fresh_clintrial()
    outcomes <- c("readmission_30d", "any_complication")

    for (columns in c("adjusted", "both")) {
        result <- multifit(df, outcomes = outcomes, predictor = "treatment",
                           covariates = c("age", "sex"), columns = columns,
                           parallel = FALSE)
        raw <- attr(result, "raw_data")
        lower <- if (columns == "both") raw$ci_lower_adj else raw$ci_lower
        upper <- if (columns == "both") raw$ci_upper_adj else raw$ci_upper

        rows <- forest_rows(suppressMessages(multiforest(result)))
        expect_equal(exp(rows$conf_low), lower)
        expect_equal(exp(rows$conf_high), upper)
    }
})


test_that("multiforest draws the coefficients of a linear multifit table", {

    df <- fresh_clintrial()

    result <- multifit(df, outcomes = c("los_days", "bmi"), predictor = "treatment",
                       covariates = c("age", "sex"), model_type = "lm",
                       columns = "both", parallel = FALSE)
    raw <- attr(result, "raw_data")

    for (column in c("adjusted", "unadjusted")) {
        rows <- forest_rows(suppressMessages(multiforest(result, column = column)))
        suffix <- if (column == "adjusted") "_adj" else "_unadj"
        expect_equal(rows$estimate, raw[[paste0("exp_coef", suffix)]])
        expect_equal(rows$conf_low, raw[[paste0("ci_lower", suffix)]])
        expect_equal(rows$conf_high, raw[[paste0("ci_upper", suffix)]])
    }
})


test_that("a different confidence level gives Wald intervals at that level", {

    df <- fresh_clintrial()

    result <- multifit(df, outcomes = "readmission_30d", predictor = "treatment",
                       covariates = c("age", "sex"), parallel = FALSE)
    raw <- attr(result, "raw_data")

    rows <- forest_rows(suppressMessages(multiforest(result, conf_level = 0.9)))
    z <- stats::qnorm(0.95)
    expect_equal(rows$conf_low, raw$coefficient - z * raw$se)
    expect_equal(rows$conf_high, raw$coefficient + z * raw$se)
})


test_that("a different confidence level gives t intervals for linear models", {

    ## The table holds the exact t intervals of each linear model, so that
    ## intervals at another level use the t quantile on its residual degrees
    ## of freedom, in either column of a two-column result
    df <- fresh_clintrial()
    result <- multifit(df, outcomes = c("los_days", "bmi"), predictor = "treatment",
                       covariates = c("age", "sex"), model_type = "lm",
                       columns = "both", parallel = FALSE)

    for (column in c("adjusted", "unadjusted")) {
        rows <- forest_rows(suppressMessages(
            multiforest(result, column = column, conf_level = 0.9)))
        for (outcome in c("los_days", "bmi")) {
            terms <- if (column == "adjusted") c("treatment", "age", "sex") else "treatment"
            model <- stats::lm(stats::reformulate(terms, outcome), data = df)
            expected <- stats::confint(model, level = 0.9)[
                c("treatmentDrug A", "treatmentDrug B"), , drop = FALSE]
            drawn <- rows[rows$var == outcome]
            expect_equal(drawn$conf_low, unname(expected[, 1]),
                         info = paste(column, outcome))
            expect_equal(drawn$conf_high, unname(expected[, 2]),
                         info = paste(column, outcome))
        }
    }
})
