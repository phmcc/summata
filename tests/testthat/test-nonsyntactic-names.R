## The fitting functions take their variables as column names, from which they
## write the model formula. A column name that is not a syntactic name, such as
## "tumor stage", must be quoted with backticks to be read in a formula as one
## variable; the functions quote such names themselves and accept them with or
## without backticks.


### * Fixtures

## clintrial with non-syntactic names for a numeric and a factor predictor,
## the time of the survival outcome, a binary outcome, and a stratum
nonsyntactic_data <- function() {
    dt <- data.table::as.data.table(fresh_clintrial())
    data.table::setnames(dt, c("age", "stage", "os_months"),
                         c("age years", "tumor stage", "os months"))
    dt[, `readmitted 30d` := readmission_30d]
    dt[, `site name` := site]
    dt[]
}


### * Formula names

test_that("formula_names quotes the column names that require it", {

    vars <- c("age", "tumor stage", "sex", "os months", "os_status", "if")
    terms <- c("age", "tumor stage", "`tumor stage`", "tumor stage:sex",
               "(1|site)", "I(age^2)", "Surv(os months, os_status)",
               "Surv(os_months,os_status)", "if", "age:sex", "unknown")

    expect_identical(
        summata:::formula_names(terms, vars),
        c("age", "`tumor stage`", "`tumor stage`", "`tumor stage`:sex",
          "(1|site)", "I(age^2)", "Surv(`os months`, os_status)",
          "Surv(os_months,os_status)", "`if`", "age:sex", "unknown")
    )

    ## Without quoting, backticks are removed so that the names match columns
    expect_identical(
        summata:::formula_names(c("`tumor stage`", "`tumor stage`:sex", "age"),
                                vars, quote = FALSE),
        c("tumor stage", "tumor stage:sex", "age")
    )
    expect_null(summata:::formula_names(NULL, vars))
})


### * fit()

test_that("fit accepts non-syntactic predictor names", {

    dt <- nonsyntactic_data()
    result <- fit(dt, "readmission_30d", c("age years", "tumor stage", "sex"),
                  model_type = "glm")
    raw <- attr(result, "raw_data")

    model <- stats::glm(readmission_30d ~ `age years` + `tumor stage` + sex,
                        data = dt, family = stats::binomial)
    estimable <- raw[!is.na(raw$se)]

    expect_identical(result$Variable[nzchar(result$Variable)],
                     c("age years", "tumor stage", "sex"))
    expect_equal(estimable$coefficient, unname(stats::coef(model)[-1]))
})


test_that("fit accepts names given with backticks", {

    dt <- nonsyntactic_data()
    plain <- fit(dt, "readmission_30d", c("age years", "tumor stage"),
                 model_type = "glm")
    quoted <- fit(dt, "readmission_30d", c("`age years`", "`tumor stage`"),
                  model_type = "glm")

    expect_identical(quoted$Variable, plain$Variable)
    expect_identical(quoted[["aOR (95% CI)"]], plain[["aOR (95% CI)"]])
})


test_that("fit accepts a non-syntactic outcome name", {

    dt <- nonsyntactic_data()
    result <- fit(dt, "readmitted 30d", c("age years", "sex"), model_type = "glm")
    reference <- fit(dt, "readmission_30d", c("age years", "sex"),
                     model_type = "glm")

    expect_identical(result[["aOR (95% CI)"]], reference[["aOR (95% CI)"]])
})


test_that("fit accepts non-syntactic names in Surv(), strata, and interactions", {

    dt <- nonsyntactic_data()

    cox <- fit(dt, "Surv(os months, os_status)", c("age years", "tumor stage"),
               model_type = "coxph")
    cox_quoted <- fit(dt, "Surv(`os months`, os_status)",
                      c("age years", "tumor stage"), model_type = "coxph")
    model <- survival::coxph(
        survival::Surv(`os months`, os_status) ~ `age years` + `tumor stage`,
        data = dt
    )
    raw <- attr(cox, "raw_data")

    expect_equal(raw$coefficient[!is.na(raw$se)], unname(stats::coef(model)))
    expect_identical(cox_quoted[["aHR (95% CI)"]], cox[["aHR (95% CI)"]])

    stratified <- fit(dt, "Surv(os months, os_status)", c("age years", "sex"),
                      strata = "site name", model_type = "coxph")
    stratified_model <- survival::coxph(
        survival::Surv(`os months`, os_status) ~ `age years` + sex +
            strata(`site name`),
        data = dt
    )
    raw_strat <- attr(stratified, "raw_data")
    expect_equal(raw_strat$coefficient[!is.na(raw_strat$se)],
                 unname(stats::coef(stratified_model)))

    interacting <- fit(dt, "readmission_30d", c("tumor stage", "sex"),
                       interactions = "tumor stage:sex", model_type = "glm")
    int_model <- stats::glm(readmission_30d ~ `tumor stage` * sex, data = dt,
                            family = stats::binomial)
    raw_int <- attr(interacting, "raw_data")
    expect_equal(sort(raw_int$coefficient[!is.na(raw_int$se)]),
                 sort(unname(stats::coef(int_model)[-1])))
})


### * uniscreen() and fullfit()

test_that("uniscreen screens non-syntactic predictor names", {

    dt <- nonsyntactic_data()
    result <- uniscreen(dt, "Surv(os months, os_status)",
                        c("age years", "tumor stage"), model_type = "coxph",
                        parallel = FALSE)
    raw <- attr(result, "raw_data")

    for (pred in c("age years", "tumor stage")) {
        model <- survival::coxph(
            stats::as.formula(paste0("survival::Surv(`os months`, os_status) ~ `",
                                     pred, "`")),
            data = dt
        )
        rows <- raw$predictor == pred & !is.na(raw$se)
        expect_equal(raw$coefficient[rows], unname(stats::coef(model)),
                     info = pred)
    }

    ## The analyzed sample of each model is counted from its variables
    complete <- vapply(c("age years", "tumor stage"), function(pred) {
        vars <- c("os months", "os_status", pred)
        sum(stats::complete.cases(dt[, vars, with = FALSE]))
    }, numeric(1))
    expect_equal(range(attr(result, "analysis_counts")$n_analyzed),
                 range(complete))
})


test_that("fullfit screens and fits non-syntactic predictor names", {

    dt <- nonsyntactic_data()
    result <- suppressMessages(
        fullfit(dt, "readmission_30d", c("age years", "tumor stage", "sex"),
                model_type = "glm", parallel = FALSE)
    )

    expect_identical(result$Variable[nzchar(result$Variable)],
                     c("age years", "tumor stage", "sex"))
    expect_false(anyNA(attr(result, "multi_raw")$coefficient))
})


### * multifit()

test_that("multifit accepts non-syntactic names and labels their interactions", {

    dt <- nonsyntactic_data()
    result <- multifit(dt, outcomes = c("readmitted 30d", "icu_admission"),
                       predictor = "tumor stage",
                       covariates = c("age years", "sex"),
                       interactions = "tumor stage:sex", parallel = FALSE)
    raw <- attr(result, "raw_data")

    model <- stats::glm(`readmitted 30d` ~ `tumor stage` + `age years` + sex +
                            `tumor stage`:sex,
                        data = dt, family = stats::binomial)
    own <- raw$outcome == "readmitted 30d" & !grepl(":", raw$group, fixed = TRUE)

    expect_identical(raw$group[own], c("II", "III", "IV"))
    expect_equal(raw$coefficient[own],
                 unname(stats::coef(model)[paste0("`tumor stage`", c("II", "III", "IV"))]))
    expect_true("tumor stage (II) × sex (Male)" %in% result$Predictor)
})


test_that("interaction labels split each component by its variable name", {

    ## Variable names outside lowercase letters, digits, and underscores
    expect_identical(
        summata:::format_interaction_term("SmokingYes:sexMale",
                                          var_names = c("Smoking", "sex")),
        "Smoking (Yes) × sex (Male)"
    )
    expect_identical(
        summata:::format_interaction_term("tumor stageII:age years",
                                          labels = c("age years" = "Age (years)"),
                                          var_names = c("tumor stage", "age years")),
        "tumor stage (II) × Age (years)"
    )
})


### * Forest plots

test_that("forest plots draw the results of non-syntactic names", {

    dt <- nonsyntactic_data()

    fitted <- fit(dt, "readmission_30d", c("age years", "tumor stage"),
                  model_type = "glm")
    p_fit <- suppressMessages(suppressWarnings(glmforest(fitted)))
    expect_identical(unique(attr(p_fit, "table_data")$var),
                     c("age years", "tumor stage"))

    screened <- uniscreen(dt, "readmission_30d", c("age years", "tumor stage"),
                          model_type = "glm", parallel = FALSE)
    expect_s3_class(suppressMessages(suppressWarnings(uniforest(screened))),
                    "ggplot")

    multi <- multifit(dt, outcomes = c("readmitted 30d", "icu_admission"),
                      predictor = "tumor stage", covariates = "age years",
                      parallel = FALSE)
    expect_s3_class(suppressMessages(suppressWarnings(multiforest(multi))),
                    "ggplot")
})


### * survtable()

test_that("survtable accepts non-syntactic time and grouping names", {

    dt <- nonsyntactic_data()
    dt[, `treatment arm` := treatment]

    result <- survtable(dt, "Surv(os months, os_status)", by = "treatment arm",
                        times = 12)
    quoted <- survtable(dt, "Surv(`os months`, os_status)", by = "`treatment arm`",
                        times = 12)
    reference <- survtable(fresh_clintrial(), "Surv(os_months, os_status)",
                           by = "treatment", times = 12)

    expect_identical(result[["12"]], reference[["12"]])
    expect_identical(quoted[["12"]], reference[["12"]])
    expect_identical(result$Variable, reference$Variable)
    expect_equal(attr(result, "test_result")$p_value,
                 attr(reference, "test_result")$p_value)
})


### * Names that contain a colon

## A colon joins the components of an interaction term, but a column name
## may contain one too; such a variable is a main effect like any other
colon_name_data <- function() {
    dt <- data.table::as.data.table(fresh_clintrial())
    dt[, `age:group` := factor(ifelse(age > 60, "Older", "Younger"),
                               levels = c("Younger", "Older"))]
    dt[, `dose:mg` := round(bmi * 10)]
    dt[]
}


test_that("a variable whose name contains a colon is not taken for an interaction", {

    ## fit() split such a name into the names of an interaction's components,
    ## which it did not find among the columns
    dt <- colon_name_data()
    vars <- c("readmission_30d", "age:group", "dose:mg", "sex")
    result <- fit(dt, "readmission_30d", c("age:group", "dose:mg", "sex"),
                  model_type = "glm", labels = c(sex = "Sex"))
    raw <- attr(result, "raw_data")
    model <- stats::glm(readmission_30d ~ `age:group` + `dose:mg` + sex,
                        data = dt, family = stats::binomial)
    counts <- counts_by_level(dt, "age:group", "readmission_30d", vars)

    expect_identical(result$Variable[nzchar(result$Variable)],
                     c("age:group", "dose:mg", "Sex"))
    group <- raw[raw$variable == "age:group"]
    expect_identical(group$group, c("Younger", "Older"))
    expect_equal(group$n, counts$n)
    expect_equal(raw$coefficient[raw$variable == "dose:mg"],
                 unname(stats::coef(model)["`dose:mg`"]))

    single <- stats::glm(readmission_30d ~ `age:group`, data = dt,
                         family = stats::binomial)
    expect_identical(summata:::detect_model_type(single), "Univariable")

    screened <- uniscreen(dt, "readmission_30d", c("age:group", "sex"),
                          model_type = "glm", labels = c(sex = "Sex"),
                          parallel = FALSE)
    expect_identical(screened$Variable[nzchar(screened$Variable)],
                     c("age:group", "Sex"))
})


test_that("forest plots draw a variable whose name contains a colon", {

    ## The forest plots dropped such a variable from their terms
    dt <- colon_name_data()
    models <- list(
        glm = stats::glm(readmission_30d ~ `age:group` + sex, data = dt,
                         family = stats::binomial),
        lm = stats::lm(los_days ~ `age:group` + sex, data = dt),
        cox = survival::coxph(survival::Surv(os_months, os_status) ~
                                  `age:group` + sex, data = dt)
    )
    forests <- list(glm = glmforest, lm = lmforest, cox = coxforest)

    for (name in names(models)) {
        rows <- attr(suppressMessages(suppressWarnings(
            forests[[name]](models[[name]], data = dt))), "table_data")
        expect_identical(rows$level[rows$var == "age:group"],
                         c("Younger", "Older"), info = name)
    }

    ## The variable's levels are read in the components of its interactions
    interacting <- stats::glm(readmission_30d ~ `age:group` * sex, data = dt,
                              family = stats::binomial)
    rows <- attr(suppressMessages(suppressWarnings(
        glmforest(interacting, data = dt))), "table_data")
    complete <- stats::complete.cases(dt[, c("readmission_30d", "age:group", "sex")])
    both <- complete & dt$`age:group` == "Older" & dt$sex == "Male"

    label <- "age:group (Older) × sex (Male)"
    expect_true(label %in% rows$var)
    expect_equal(rows$N[rows$var == label], sum(both))
})
