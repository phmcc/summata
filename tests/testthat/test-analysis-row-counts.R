## Group-level sample sizes must be counted over the observations that entered
## the model. Fitted objects that retain no model frame (survival::coxph() by
## default, coxme, and the memory-conserving model = FALSE fits used internally)
## previously fell back to the unrestricted data, so categorical predictors
## reported totals over every supplied row while continuous predictors reported
## complete cases only.

## Sample sizes for a single variable, including its reference level.
##
## Indexing is performed on extracted vectors rather than through the data.table
## '[' method: an argument named for a column is otherwise resolved against that
## column rather than against the calling environment.
group_sizes <- function(raw, var_name) {
    raw[["n"]][raw[["variable"]] == var_name]
}

## Event counts for a single variable, including its reference level.
group_events <- function(raw, var_name) {
    raw[["events"]][raw[["variable"]] == var_name]
}

## Count events using the same convention as m2dt(): for factor outcomes the
## second level is the event.
count_events <- function(x) {
    if (is.factor(x)) {
        return(sum(as.integer(x) == 2L, na.rm = TRUE))
    }
    return(sum(x, na.rm = TRUE))
}


test_that("uniscreen group sizes exclude observations with a missing response", {

    dt <- with_missing_outcome("readmission_30d", 100L)
    n_complete <- sum(!is.na(dt$readmission_30d))

    result <- uniscreen(
        data = dt,
        outcome = "readmission_30d",
        predictors = c("age", "sex", "race", "ethnicity"),
        model_type = "glm",
        family = "binomial",
        parallel = FALSE
    )

    raw <- attr(result, "raw_data")

    ## Continuous predictor: n is the fitted sample size
    expect_equal(unique(group_sizes(raw, "age")), n_complete)

    ## Categorical predictors: group sizes partition the fitted sample size
    for (pred in c("sex", "race", "ethnicity")) {
        expect_equal(sum(group_sizes(raw, pred)), n_complete)
    }
})


test_that("uniscreen event counts remain consistent with group sizes", {

    dt <- with_missing_outcome("readmission_30d", 100L)
    total_events <- count_events(dt$readmission_30d)

    result <- uniscreen(
        data = dt,
        outcome = "readmission_30d",
        predictors = c("sex", "race"),
        model_type = "glm",
        family = "binomial",
        parallel = FALSE
    )

    raw <- attr(result, "raw_data")

    for (pred in c("sex", "race")) {
        expect_equal(sum(group_events(raw, pred)), total_events)
        expect_true(all(group_events(raw, pred) <= group_sizes(raw, pred)))
    }
})


test_that("multivariable fit group sizes match the fitted sample size", {

    dt <- with_missing_outcome("readmission_30d", 100L)

    result <- fit(
        data = dt,
        outcome = "readmission_30d",
        predictors = c("age", "sex", "race"),
        model_type = "glm",
        family = "binomial"
    )

    raw <- attr(result, "raw_data")
    n_model <- stats::nobs(attr(result, "model"))

    for (pred in c("sex", "race")) {
        expect_equal(sum(group_sizes(raw, pred)), n_model)
    }
    expect_equal(unique(group_sizes(raw, "age")), n_model)
})


test_that("Cox models count only the observations used in fitting", {

    dt <- with_missing_outcome("os_status", 100L)

    result <- fit(
        data = dt,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "sex", "race"),
        model_type = "coxph"
    )

    raw <- attr(result, "raw_data")
    model <- attr(result, "model")

    for (pred in c("sex", "race")) {
        expect_equal(sum(group_sizes(raw, pred)), model$n)
        expect_equal(sum(group_events(raw, pred)), model$nevent)
    }
})


test_that("group counts do not depend on whether the model frame was retained", {

    dt <- with_missing_outcome("readmission_30d", 100L)

    with_frame <- stats::glm(readmission_30d ~ age + sex + race,
                             data = dt, family = stats::binomial,
                             model = TRUE)
    without_frame <- stats::glm(readmission_30d ~ age + sex + race,
                                data = dt, family = stats::binomial,
                                model = FALSE, y = TRUE)

    kept <- m2dt(dt, with_frame, keep_qc_stats = FALSE)
    dropped <- m2dt(dt, without_frame, keep_qc_stats = FALSE)

    expect_equal(kept[["n"]], dropped[["n"]])
    expect_equal(kept[["events"]], dropped[["events"]])
})


test_that("factor levels with no analysis rows report zero, not the model total", {

    ## Every observation in one race category has a missing event indicator, so
    ## that level contributes no analysis rows while remaining in the factor.
    ## Unlike stats::glm(), survival::coxph() retains the singular coefficient,
    ## so the emptied level still appears as a row in the output.
    df <- fresh_clintrial()
    df$os_status[df$race == "Other"] <- NA
    dt <- data.table::as.data.table(df)

    result <- fit(
        data = dt,
        outcome = "Surv(os_months, os_status)",
        predictors = c("age", "race"),
        model_type = "coxph"
    )

    raw <- attr(result, "raw_data")
    model <- attr(result, "model")

    ## Group sizes still partition the fitted sample size
    expect_equal(sum(group_sizes(raw, "race")), model$n)
    expect_equal(sum(group_events(raw, "race")), model$nevent)

    ## The emptied level reports zero rather than inheriting the total
    empty_row <- raw[["variable"]] == "race" & raw[["group"]] == "Other"

    if (any(empty_row)) {
        expect_equal(raw[["n"]][empty_row], 0)
        expect_equal(raw[["events"]][empty_row], 0)
    }
})


test_that("complete data is unaffected by the analysis-row restriction", {

    dt <- data.table::as.data.table(fresh_clintrial())
    n_total <- nrow(dt)

    result <- uniscreen(
        data = dt,
        outcome = "readmission_30d",
        predictors = c("sex", "race"),
        model_type = "glm",
        family = "binomial",
        parallel = FALSE
    )

    raw <- attr(result, "raw_data")

    for (pred in c("sex", "race")) {
        expect_equal(sum(group_sizes(raw, pred)), n_total)
    }
})
