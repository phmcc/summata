## The rows of a regression table follow the model: variables in the order of
## the model's terms, which is the order of the formula with interactions after
## the main effects, and the levels of a factor in the order of its levels.
## Tables of a model passed in carry no list of predictors, so their order
## comes from the model itself rather than from sorting the variable names.


### * Fixtures

## Formula order that is neither alphabetical nor reverse alphabetical
order_glm <- function(data = fresh_clintrial()) {
    stats::glm(readmission_30d ~ stage + age + treatment + sex,
               data = data, family = stats::binomial)
}

order_terms <- c("stage", "age", "treatment", "sex")


## Variable names of the table rows that start a variable, in order
table_variables <- function(result) {
    vars <- result$Variable
    vars[nzchar(vars)]
}


### * m2dt()

test_that("m2dt lists variables in the order of the model's terms", {

    df <- fresh_clintrial()
    model <- order_glm(df)

    with_ref <- m2dt(df, model)
    without_ref <- m2dt(df, model, reference_rows = FALSE)

    expect_identical(unique(with_ref$variable), c("(Intercept)", order_terms))
    expect_identical(unique(without_ref$variable), unique(with_ref$variable))
})


test_that("reference rows precede the levels of their variable", {

    df <- fresh_clintrial()
    result <- m2dt(df, order_glm(df))

    for (var in c("stage", "treatment", "sex")) {
        rows <- result[result$variable == var]
        expect_identical(rows$group, levels(df[[var]]), info = var)
        expect_identical(rows$reference[1], "reference", info = var)
    }
})


test_that("interaction terms follow the main effects of a Cox model", {

    df <- fresh_clintrial()
    model <- survival::coxph(
        survival::Surv(os_months, os_status) ~ treatment + age + stage * sex,
        data = df
    )

    result <- m2dt(df, model)
    vars <- unique(result$variable)
    interactions <- grep(":", vars, fixed = TRUE, value = TRUE)

    expect_identical(setdiff(vars, interactions),
                     c("treatment", "age", "stage", "sex"))
    ## The interactions come last, in the order of their coefficients
    expect_identical(vars[seq_along(interactions) + 4L], interactions)
    expect_identical(interactions,
                     grep(":", names(stats::coef(model)), fixed = TRUE,
                          value = TRUE))
})


test_that("mixed models list their fixed effects in the order of the formula", {

    skip_if_not_installed("lme4")

    df <- fresh_clintrial()
    model <- suppressMessages(suppressWarnings(
        lme4::glmer(readmission_30d ~ stage + age + sex + (1 | site),
                    data = df, family = stats::binomial)
    ))

    result <- m2dt(df, model)

    expect_identical(unique(result$variable),
                     c("(Intercept)", "stage", "age", "sex"))
})


### * fit()

test_that("fit with a model passed in follows the order of its formula", {

    df <- fresh_clintrial()
    result <- fit(data = df, model = order_glm(df))

    expect_identical(table_variables(result), order_terms)
})


test_that("fit with a model passed in matches the table it would fit itself", {

    df <- fresh_clintrial()

    passed <- fit(data = df, model = order_glm(df))
    fitted <- fit(data = df, outcome = "readmission_30d",
                  predictors = order_terms, model_type = "glm")

    expect_identical(passed$Variable, fitted$Variable)
    expect_identical(passed$Group, fitted$Group)
})
