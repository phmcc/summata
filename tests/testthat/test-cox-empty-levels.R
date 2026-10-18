## Cox models keep factor levels without observations, which lm() and glm()
## drop. With an empty first level, the other levels were estimated against
## whichever of them was aliased while the tables named the empty level as the
## reference, and coxme() failed to fit. fit(), uniscreen() and multifit() now
## drop empty levels before fitting a Cox model, which then takes the first
## level with observations as its reference, as the other models do.

## Stage I is removed, while the factor keeps the level
without_stage_i <- function() {
    subset(fresh_clintrial(), stage != "I")
}


test_that("Cox models take the first level with observations as the reference", {

    df <- without_stage_i()
    expect_identical(levels(df$stage)[1], "I")

    cox <- attr(fit(data = df, outcome = "Surv(os_months, os_status)",
                    predictors = c("age", "stage"), model_type = "coxph"),
                "raw_data")
    logistic <- attr(fit(data = df, outcome = "os_status",
                         predictors = c("age", "stage"), model_type = "glm"),
                     "raw_data")

    expect_identical(cox$group[cox$variable == "stage"], c("II", "III", "IV"))
    expect_identical(cox$reference[cox$variable == "stage"],
                     c("reference", "", ""))
    expect_identical(cox$group[cox$variable == "stage"],
                     logistic$group[logistic$variable == "stage"])

    ## The estimates are those of the model fitted without the empty level
    reference <- survival::coxph(survival::Surv(os_months, os_status) ~ age + stage,
                                 data = droplevels(df))
    expect_equal(cox$coefficient[cox$group %in% c("III", "IV")],
                 unname(stats::coef(reference)[c("stageIII", "stageIV")]))

    ## The supplied data keep their levels
    expect_identical(levels(df$stage), c("I", "II", "III", "IV"))
})


test_that("uniscreen and multifit fit Cox models without the empty levels", {

    df <- without_stage_i()

    screened <- attr(uniscreen(df, "Surv(os_months, os_status)", "stage",
                               model_type = "coxph", parallel = FALSE),
                     "raw_data")
    expect_identical(screened$group, c("II", "III", "IV"))

    reference <- survival::coxph(survival::Surv(os_months, os_status) ~ stage + age,
                                 data = droplevels(df))
    multi <- attr(multifit(df, outcomes = "Surv(os_months, os_status)",
                           predictor = "stage", covariates = "age",
                           model_type = "coxph", parallel = FALSE),
                  "raw_data")
    expect_identical(multi$group, c("III", "IV"))
    expect_equal(multi$coefficient,
                 unname(stats::coef(reference)[c("stageIII", "stageIV")]))
})


test_that("coxme models are fitted when a factor has an empty level", {

    skip_if_not_installed("coxme")

    result <- fit(data = without_stage_i(), outcome = "Surv(os_months, os_status)",
                  predictors = c("age", "stage"), model_type = "coxme",
                  random = "(1|site)")
    raw <- attr(result, "raw_data")

    expect_identical(raw$group[raw$variable == "stage"], c("II", "III", "IV"))
})


test_that("levels are empty when no complete case of the model takes them", {

    ## Progression-free survival is missing for every stage I patient, so that
    ## the model of that outcome uses no stage I observation although the data
    ## have some
    df <- fresh_clintrial()
    df$pfs_months[df$stage %in% "I"] <- NA

    result <- multifit(df, outcomes = c("Surv(os_months, os_status)",
                                        "Surv(pfs_months, pfs_status)"),
                       predictor = "stage", covariates = "age",
                       model_type = "coxph", parallel = FALSE)
    raw <- attr(result, "raw_data")

    used <- df[stats::complete.cases(df[, c("pfs_months", "pfs_status", "age",
                                            "stage")]), ]
    reference <- survival::coxph(survival::Surv(pfs_months, pfs_status) ~ stage + age,
                                 data = droplevels(used))
    pfs <- raw[raw$outcome == "Surv(pfs_months, pfs_status)"]
    expect_identical(pfs$group, c("III", "IV"))
    expect_equal(pfs$coefficient,
                 unname(stats::coef(reference)[c("stageIII", "stageIV")]))

    ## The model of overall survival keeps stage I as its reference
    os <- raw[raw$outcome == "Surv(os_months, os_status)"]
    expect_identical(os$group, c("II", "III", "IV"))
})


test_that("the supplied data table keeps its levels", {

    dt <- data.table::as.data.table(without_stage_i())
    levels_before <- levels(dt$stage)

    fit(data = dt, outcome = "Surv(os_months, os_status)",
        predictors = c("age", "stage"), model_type = "coxph")

    expect_identical(levels(dt$stage), levels_before)
})


test_that("coxme models are fitted when a level occurs only among incomplete cases", {

    skip_if_not_installed("coxme")

    df <- fresh_clintrial()
    df$age[df$stage %in% "I"] <- NA

    result <- fit(data = df, outcome = "Surv(os_months, os_status)",
                  predictors = c("age", "stage"), model_type = "coxme",
                  random = "(1|site)")
    raw <- attr(result, "raw_data")

    expect_identical(raw$group[raw$variable == "stage"], c("II", "III", "IV"))
})
