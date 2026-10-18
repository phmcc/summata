## multifit() lists outcomes in the order given and the levels of a factor
## predictor in their defined order, whichever columns are displayed. The
## combined unadjusted and adjusted table previously came out of merge()
## sorted alphabetically by outcome and level.

test_that("both columns keep the order of the outcomes and of the levels", {

    df <- fresh_clintrial()
    outcomes <- c("readmission_30d", "any_complication", "icu_admission")

    both <- multifit(df, outcomes = outcomes, predictor = "smoking",
                     covariates = c("age", "sex"), columns = "both", parallel = FALSE)
    adjusted <- multifit(df, outcomes = outcomes, predictor = "smoking",
                         covariates = c("age", "sex"), parallel = FALSE)
    unadjusted <- multifit(df, outcomes = outcomes, predictor = "smoking",
                           covariates = c("age", "sex"), columns = "unadjusted",
                           parallel = FALSE)

    expect_identical(both$Outcome, adjusted$Outcome)
    expect_identical(both$Predictor, adjusted$Predictor)
    expect_identical(both$Predictor, unadjusted$Predictor)
    expect_identical(unique(both$Outcome), outcomes)

    raw <- attr(both, "raw_data")
    expect_identical(raw$group, rep(c("Former", "Current"), 3))
    expect_null(data.table::key(raw))
    expect_false(anyNA(raw$effect_type))

    ## Estimates are those of the single-column results
    expect_equal(raw$exp_coef_adj, attr(adjusted, "raw_data")$exp_coef)
    expect_equal(raw$exp_coef_unadj, attr(unadjusted, "raw_data")$exp_coef)
})
