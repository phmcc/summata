## R names the coefficient of a factor level by the variable and the level
## ("shift12:00") and joins the components of an interaction with colons
## ("shift12:00:sexMale"). A level that contains a colon, such as a time of
## day or a ratio, is read as a level, not as an interaction, wherever its
## levels are known.


### * Fixtures

colon_data <- function() {
    dt <- data.table::as.data.table(fresh_clintrial())
    set.seed(20261016)
    dt[, shift := factor(sample(c("08:00", "12:00", "18:00"), .N, replace = TRUE),
                         levels = c("08:00", "12:00", "18:00"))]
    dt[, ratio := factor(sample(c("1:1", "2:1"), .N, replace = TRUE))]
    dt[]
}


### * Tables

test_that("fit and m2dt read a level that contains a colon as a level", {

    dt <- colon_data()
    result <- fit(dt, "readmission_30d", c("shift", "sex", "ratio"),
                  model_type = "glm")
    raw <- attr(result, "raw_data")
    model <- stats::glm(readmission_30d ~ shift + sex + ratio, data = dt,
                        family = stats::binomial)
    counts <- counts_by_level(dt, "shift", "readmission_30d")

    expect_identical(result$Variable[nzchar(result$Variable)],
                     c("shift", "sex", "ratio"))
    shift <- raw[raw$variable == "shift"]
    expect_identical(shift$group, levels(dt$shift))
    expect_identical(shift$reference, c("reference", "", ""))
    expect_equal(shift$n, counts$n)
    expect_equal(shift$events, counts$events)
    expect_equal(shift$coefficient[-1],
                 unname(stats::coef(model)[c("shift12:00", "shift18:00")]))
    expect_identical(raw$group[raw$variable == "ratio"], c("1:1", "2:1"))

    without_ref <- m2dt(dt, model, reference_rows = FALSE)
    expect_identical(without_ref$group[without_ref$variable == "shift"],
                     c("12:00", "18:00"))
})


test_that("a model of one factor with colon levels is univariable", {

    dt <- colon_data()
    single <- stats::glm(readmission_30d ~ shift, data = dt,
                         family = stats::binomial)
    two <- stats::glm(readmission_30d ~ shift + ratio, data = dt,
                      family = stats::binomial)

    expect_identical(summata:::detect_model_type(single), "Univariable")
    expect_identical(summata:::detect_model_type(two), "Multivariable")
})


test_that("uniscreen and multifit read levels that contain a colon", {

    dt <- colon_data()
    counts <- counts_by_level(dt, "shift", "readmission_30d")

    screened <- uniscreen(dt, "readmission_30d", c("shift", "sex"),
                          model_type = "glm", parallel = FALSE)
    expect_identical(screened$Variable[1:3], c("shift", "", ""))
    expect_identical(screened$Group[1:3], levels(dt$shift))

    multi <- multifit(dt, outcomes = "readmission_30d", predictor = "shift",
                      covariates = "sex", parallel = FALSE)
    expect_identical(multi$Predictor, c("shift (12:00)", "shift (18:00)"))
    expect_equal(attr(multi, "raw_data")$n, counts$n[-1])
})


### * Interactions

test_that("an interaction with a factor of colon levels is split by its levels", {

    dt <- colon_data()
    components <- summata:::split_interaction_term(
        "shift12:00:sexMale", c("shift", "sex"),
        list(shift = levels(dt$shift), sex = levels(dt$sex))
    )

    expect_identical(vapply(components, `[[`, "", "var"), c("shift", "sex"))
    expect_identical(vapply(components, `[[`, "", "level"), c("12:00", "Male"))

    model <- stats::glm(readmission_30d ~ shift * sex, data = dt,
                        family = stats::binomial)
    raw <- m2dt(dt, model)
    both <- dt$shift == "12:00" & dt$sex == "Male"
    events <- as.integer(dt$readmission_30d) == 2L

    row <- raw[raw$term == "shift12:00:sexMale"]
    expect_equal(row$n, sum(both))
    expect_equal(row$events, sum(events[both]))
})


test_that("a level that contains a colon is read in a later component", {

    ## The parts of a level were joined again in the first component alone,
    ## so that "sexMale:shift12:00" was read as three components
    dt <- colon_data()
    components <- summata:::split_interaction_term(
        "sexMale:shift12:00", c("sex", "shift"),
        list(sex = levels(dt$sex), shift = levels(dt$shift))
    )

    expect_identical(vapply(components, `[[`, "", "var"), c("sex", "shift"))
    expect_identical(vapply(components, `[[`, "", "level"), c("Male", "12:00"))

    model <- stats::glm(readmission_30d ~ sex * shift, data = dt,
                        family = stats::binomial)
    rows <- attr(suppressMessages(suppressWarnings(glmforest(model, data = dt))),
                 "table_data")
    both <- dt$shift == "12:00" & dt$sex == "Male"

    label <- "sex (Male) × shift (12:00)"
    expect_true(label %in% rows$var)
    expect_equal(rows$N[rows$var == label], sum(both))
})


test_that("multifit reports the interactions of factors with colon levels", {

    dt <- colon_data()
    model <- stats::glm(readmission_30d ~ shift + sex + shift:sex, data = dt,
                        family = stats::binomial)

    ## With the predictor's levels split at their colons, its interactions
    ## were not found among the model's coefficients
    as_predictor <- multifit(dt, outcomes = "readmission_30d", predictor = "shift",
                             covariates = "sex", interactions = "shift:sex",
                             columns = "adjusted", parallel = FALSE)
    expect_identical(as_predictor$Predictor,
                     c("shift (12:00)", "shift (18:00)",
                       "shift (12:00) × sex (Male)",
                       "shift (18:00) × sex (Male)"))
    expect_equal(attr(as_predictor, "raw_data")$coefficient,
                 unname(stats::coef(model)[c("shift12:00", "shift18:00",
                                             "shift12:00:sexMale",
                                             "shift18:00:sexMale")]))

    ## A covariate's levels are read in the labels of the interactions too
    as_covariate <- multifit(dt, outcomes = "readmission_30d", predictor = "sex",
                             covariates = "shift", interactions = "sex:shift",
                             columns = "adjusted", parallel = FALSE)
    expect_identical(as_covariate$Predictor,
                     c("sex (Male)", "sex (Male) × shift (12:00)",
                       "sex (Male) × shift (18:00)"))
})


### * Forest plots

test_that("forest plots draw each level that contains a colon once", {

    dt <- colon_data()
    glm_model <- stats::glm(readmission_30d ~ shift + sex, data = dt,
                            family = stats::binomial)
    cox_model <- survival::coxph(
        survival::Surv(os_months, os_status) ~ shift + sex, data = dt
    )
    lm_model <- stats::lm(los_days ~ shift + sex, data = dt)

    plots <- list(
        glm = suppressMessages(suppressWarnings(glmforest(glm_model, data = dt))),
        cox = suppressMessages(suppressWarnings(coxforest(cox_model, data = dt))),
        lm = suppressMessages(suppressWarnings(lmforest(lm_model, data = dt)))
    )

    for (name in names(plots)) {
        rows <- attr(plots[[name]], "table_data")
        expect_identical(rows$var, rep(c("shift", "sex"), c(3L, 2L)), info = name)
        expect_identical(rows$level[rows$var == "shift"], levels(dt$shift),
                         info = name)
    }

    int_model <- stats::glm(readmission_30d ~ shift * sex, data = dt,
                            family = stats::binomial)
    int_rows <- attr(suppressMessages(suppressWarnings(
        glmforest(int_model, data = dt))), "table_data")
    both <- dt$shift == "12:00" & dt$sex == "Male"

    label <- "shift (12:00) × sex (Male)"
    expect_true(label %in% int_rows$var)
    expect_equal(int_rows$N[int_rows$var == label], sum(both))
})
