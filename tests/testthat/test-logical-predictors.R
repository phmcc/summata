## R codes a logical predictor as a factor with the levels FALSE and TRUE,
## naming its coefficient by TRUE ("olderTRUE") with FALSE as the reference.
## Tables and forest plots show it as they show a factor: by level, with a
## reference row for FALSE and the counts of each level.


### * Fixtures

logical_data <- function() {
    dt <- data.table::as.data.table(fresh_clintrial())
    dt[, older := age > 65]
    dt[, smoker := smoking == "Current"]
    dt[]
}


### * Tables

test_that("fit shows a logical predictor by level with FALSE as reference", {

    dt <- logical_data()
    expect_no_warning(
        result <- fit(dt, "readmission_30d", c("older", "sex", "smoker"),
                      model_type = "glm")
    )
    raw <- attr(result, "raw_data")
    model <- stats::glm(readmission_30d ~ older + sex + smoker, data = dt,
                        family = stats::binomial)
    counts <- counts_by_level(dt, "older", "readmission_30d",
                              c("readmission_30d", "older", "sex", "smoker"))

    ## Variables in the order given, each logical by level
    expect_identical(result$Variable[nzchar(result$Variable)],
                     c("older", "sex", "smoker"))
    older <- raw[raw$variable == "older"]
    expect_identical(older$group, c("FALSE", "TRUE"))
    expect_identical(older$reference, c("reference", ""))
    expect_equal(older$n, counts$n)
    expect_equal(older$events, counts$events)
    expect_equal(older$coefficient[2], unname(stats::coef(model)["olderTRUE"]))
})


test_that("m2dt parses a logical coefficient into its variable and level", {

    dt <- logical_data()
    model <- stats::glm(readmission_30d ~ older + sex, data = dt,
                        family = stats::binomial)

    with_ref <- m2dt(dt, model)
    without_ref <- m2dt(dt, model, reference_rows = FALSE)

    expect_identical(unique(with_ref$variable), c("(Intercept)", "older", "sex"))
    expect_identical(without_ref$group[without_ref$variable == "older"], "TRUE")
})


test_that("uniscreen and Cox tables show logical predictors by level", {

    dt <- logical_data()

    screened <- uniscreen(dt, "readmission_30d", c("older", "sex"),
                          model_type = "glm", parallel = FALSE)
    expect_identical(screened$Group[1:2], c("FALSE", "TRUE"))

    cox <- fit(dt, "Surv(os_months, os_status)", c("older", "sex"),
               model_type = "coxph")
    raw <- attr(cox, "raw_data")
    counts <- counts_by_level(dt, "older", "os_status",
                              c("os_months", "os_status", "older", "sex"))

    older <- raw[raw$variable == "older"]
    expect_identical(older$group, c("FALSE", "TRUE"))
    expect_equal(older$events, counts$events)
})


test_that("mixed models show logical predictors by level", {

    skip_if_not_installed("lme4")

    dt <- logical_data()
    model <- suppressMessages(suppressWarnings(
        lme4::glmer(readmission_30d ~ older + sex + (1 | site), data = dt,
                    family = stats::binomial)
    ))
    result <- m2dt(dt, model)

    expect_identical(result$group[result$variable == "older"], c("FALSE", "TRUE"))
    expect_identical(result$reference[result$variable == "older"],
                     c("reference", ""))
})


test_that("coxme models show logical predictors by level", {

    skip_if_not_installed("coxme")

    dt <- logical_data()
    model <- coxme::coxme(
        survival::Surv(os_months, os_status) ~ older + sex + (1 | site),
        data = dt
    )
    result <- m2dt(dt, model)
    plot <- suppressMessages(suppressWarnings(coxforest(model, data = dt)))
    rows <- attr(plot, "table_data")

    expect_identical(result$group[result$variable == "older"], c("FALSE", "TRUE"))
    expect_identical(rows$level[rows$var == "older"], c("FALSE", "TRUE"))
})


### * Forest plots

test_that("forest plots draw a logical predictor by level", {

    dt <- logical_data()

    glm_model <- stats::glm(readmission_30d ~ older + sex, data = dt,
                            family = stats::binomial)
    cox_model <- survival::coxph(
        survival::Surv(os_months, os_status) ~ older + sex, data = dt
    )
    lm_model <- stats::lm(los_days ~ older + sex, data = dt)

    plots <- list(
        glm = suppressMessages(suppressWarnings(glmforest(glm_model, data = dt))),
        cox = suppressMessages(suppressWarnings(coxforest(cox_model, data = dt))),
        lm = suppressMessages(suppressWarnings(lmforest(lm_model, data = dt)))
    )
    estimates <- list(glm = stats::coef(glm_model)["olderTRUE"],
                      cox = stats::coef(cox_model)["olderTRUE"],
                      lm = stats::coef(lm_model)["olderTRUE"])
    outcomes <- list(glm = "readmission_30d", cox = c("os_months", "os_status"),
                     lm = "los_days")

    for (name in names(plots)) {
        rows <- attr(plots[[name]], "table_data")
        older <- rows[rows$var == "older"]
        counts <- counts_by_level(dt, "older",
                                  outcomes[[name]][length(outcomes[[name]])],
                                  c(outcomes[[name]], "older", "sex"))
        expect_identical(older$level, c("FALSE", "TRUE"), info = name)
        expect_equal(older$estimate, c(0, unname(estimates[[name]])),
                     info = name)
        expect_equal(older$N, counts$n, info = name)
    }

    screened <- uniscreen(dt, "readmission_30d", c("older", "sex"),
                          model_type = "glm", parallel = FALSE)
    uni_rows <- attr(suppressMessages(suppressWarnings(uniforest(screened))),
                     "table_data")
    expect_identical(uni_rows$level[uni_rows$predictor == "older"],
                     c("FALSE", "TRUE"))
})


test_that("a condensed forest plot draws a logical predictor as one row", {

    dt <- logical_data()
    model <- stats::glm(readmission_30d ~ older + sex, data = dt,
                        family = stats::binomial)

    rows <- attr(suppressMessages(suppressWarnings(
        glmforest(model, data = dt, condense_table = TRUE))), "table_data")

    expect_identical(rows$var[1], "older")
    expect_equal(rows$estimate[1], unname(stats::coef(model)["olderTRUE"]))
})


test_that("interactions with logical predictors keep their rows", {

    dt <- logical_data()
    model <- stats::glm(readmission_30d ~ older * sex, data = dt,
                        family = stats::binomial)

    result <- m2dt(dt, model)
    rows <- attr(suppressMessages(suppressWarnings(glmforest(model, data = dt))),
                 "table_data")

    expect_true("olderTRUE:sexMale" %in% result$variable)
    expect_identical(result$group[result$variable == "older"], c("FALSE", "TRUE"))
    expect_true("older (TRUE) × sex (Male)" %in% rows$var)
})


test_that("a logical expression is drawn and counted by its coefficient", {

    ## I(age > 65) is a logical term but not a variable of the data, so that
    ## its levels cannot be counted there; it keeps its single row
    dt <- logical_data()
    cox <- survival::coxph(survival::Surv(os_months, os_status) ~ I(age > 65) +
                               sex, data = dt)
    glm_model <- stats::glm(readmission_30d ~ I(age > 65) + smoker, data = dt,
                            family = stats::binomial, model = FALSE)

    rows <- attr(suppressMessages(suppressWarnings(coxforest(cox, data = dt))),
                 "table_data")
    expect_identical(rows$var[1], "I(age > 65)TRUE")

    raw <- m2dt(dt, cox)
    expect_identical(raw$variable[1], "I(age > 65)TRUE")
    expect_equal(raw$n[1], nrow(dt))

    glm_rows <- attr(suppressMessages(suppressWarnings(
        glmforest(glm_model, data = dt))), "table_data")
    expect_identical(glm_rows$var, c("I(age > 65)TRUE", "smoker", "smoker"))
    expect_identical(glm_rows$level[2:3], c("FALSE", "TRUE"))
})
