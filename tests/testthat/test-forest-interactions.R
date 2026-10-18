## An interaction between two factors is drawn as one row labeled, for
## example, "hypertension (Yes) x diabetes (Yes)". With indented or condensed
## layouts, that row must still show the interaction estimate and the number
## of observations in both levels, not those of a main effect whose variable
## and level happen to match the label.

interaction_row <- function(p) {
    tbl <- attr(p, "table_data")
    tbl[grepl("×", tbl$var, fixed = TRUE)]
}


test_that("interaction rows keep their own estimate in every layout", {

    df <- fresh_clintrial()
    term <- "hypertensionYes:diabetesYes"

    models <- list(
        glm = stats::glm(os_status ~ age + hypertension * diabetes, data = df,
                         family = stats::binomial),
        lm = stats::lm(bmi ~ age + hypertension * diabetes, data = df),
        cox = survival::coxph(survival::Surv(os_months, os_status) ~ age +
                                  hypertension * diabetes, data = df)
    )
    plotters <- list(glm = glmforest, lm = lmforest, cox = coxforest)

    for (type in names(models)) {
        model <- models[[type]]
        used <- df[stats::complete.cases(df[, all.vars(stats::formula(model))]), ]
        n_both <- sum(used$hypertension == "Yes" & used$diabetes == "Yes")

        for (layout in list(c(FALSE, FALSE), c(TRUE, FALSE), c(TRUE, TRUE))) {
            p <- suppressMessages(suppressWarnings(
                plotters[[type]](model, data = df, indent_groups = layout[1],
                                 condense_table = layout[2])
            ))
            row <- interaction_row(p)
            expect_equal(nrow(row), 1L)
            expect_equal(row$estimate, unname(stats::coef(model)[term]),
                         info = paste(type, "indent", layout[1], "condense", layout[2]))
            expect_equal(row$N, n_both,
                         info = paste(type, "indent", layout[1], "condense", layout[2]))
        }
    }
})


test_that("interactions of non-syntactic names keep their estimate and counts", {

    ## Coefficients of names that are not syntactic carry backticks
    ## ("`high bp`Yes:`has diabetes`Yes"), which the variable names do not, and
    ## the interaction row was drawn as a reference category with no count
    df <- fresh_clintrial()
    df$`high bp` <- df$hypertension
    df$`has diabetes` <- df$diabetes
    term <- "`high bp`Yes:`has diabetes`Yes"

    models <- list(
        glm = stats::glm(os_status ~ age + `high bp` * `has diabetes`, data = df,
                         family = stats::binomial),
        lm = stats::lm(bmi ~ age + `high bp` * `has diabetes`, data = df),
        cox = survival::coxph(survival::Surv(os_months, os_status) ~ age +
                                  `high bp` * `has diabetes`, data = df)
    )
    plotters <- list(glm = glmforest, lm = lmforest, cox = coxforest)

    for (type in names(models)) {
        model <- models[[type]]
        used <- df[stats::complete.cases(df[, all.vars(stats::formula(model))]), ]
        both <- used$`high bp` == "Yes" & used$`has diabetes` == "Yes"

        for (layout in list(c(FALSE, FALSE), c(TRUE, FALSE), c(TRUE, TRUE))) {
            info <- paste(type, "indent", layout[1], "condense", layout[2])
            p <- suppressMessages(suppressWarnings(
                plotters[[type]](model, data = df, indent_groups = layout[1],
                                 condense_table = layout[2])
            ))
            row <- interaction_row(p)
            expect_equal(nrow(row), 1L, info = info)
            expect_identical(row$var, "high bp (Yes) × has diabetes (Yes)",
                             info = info)
            expect_equal(row$estimate, unname(stats::coef(model)[term]), info = info)
            expect_equal(row$N, sum(both), info = info)
            if (type != "lm") {
                expect_equal(row$Events, sum(used$os_status[both]), info = info)
            }
        }
    }
})


test_that("interaction components are assigned to the longest variable name", {

    ## "sex" precedes "sexuality" in the formula, the order in which the
    ## component "sexualityOther" was assigned to "sex" with level
    ## "ualityOther", leaving no observations in the interaction
    set.seed(20240103)
    df <- data.frame(
        y = stats::rbinom(300, 1, 0.5),
        sex = factor(sample(c("Female", "Male"), 300, replace = TRUE)),
        sexuality = factor(sample(c("Hetero", "Other"), 300, replace = TRUE))
    )
    both <- df$sex == "Male" & df$sexuality == "Other"

    model <- stats::glm(y ~ sex * sexuality, data = df, family = stats::binomial)
    row <- interaction_row(suppressMessages(suppressWarnings(
        glmforest(model, data = df))))

    expect_identical(row$var, "sex (Male) × sexuality (Other)")
    expect_equal(row$estimate, unname(stats::coef(model)["sexMale:sexualityOther"]))
    expect_equal(row$N, sum(both))
    expect_equal(row$Events, sum(df$y[both]))
})


test_that("coxme forest plots draw the main effects of an interaction written with *", {

    ## The predictors of a coxme model were read by splitting the formula at
    ## "+", which kept "hypertension * diabetes" as one continuous term drawn
    ## with the reference label, and left out both main effects
    skip_if_not_installed("coxme")

    df <- fresh_clintrial()
    model <- coxme::coxme(survival::Surv(os_months, os_status) ~ age +
                              hypertension * diabetes + (1 | site), data = df)
    estimates <- coxme::fixef(model)

    p <- suppressMessages(suppressWarnings(coxforest(model, data = df)))
    tbl <- attr(p, "table_data")

    expect_identical(unique(tbl$var[!grepl("×", tbl$var, fixed = TRUE)]),
                     c("age", "hypertension", "diabetes"))
    expect_equal(tbl$estimate[tbl$var == "hypertension" & tbl$level == "Yes"],
                 unname(estimates["hypertensionYes"]))

    used <- df[stats::complete.cases(df[, c("os_months", "os_status", "age",
                                            "hypertension", "diabetes", "site")]), ]
    both <- used$hypertension == "Yes" & used$diabetes == "Yes"
    row <- interaction_row(p)

    expect_equal(row$estimate, unname(estimates["hypertensionYes:diabetesYes"]))
    expect_equal(row$N, sum(both))
    expect_equal(row$Events, sum(used$os_status[both]))
})


test_that("interactions with a logical variable are drawn once and counted", {

    ## The coefficients of a logical variable were found by a pattern that
    ## also matched its interactions, which were drawn a second time, and the
    ## interaction counted the observations of the factor level alone
    df <- fresh_clintrial()
    df$older <- df$age > 65
    used <- df[stats::complete.cases(df[, c("os_status", "age", "sex")]), ]
    both <- used$older & used$sex == "Male"

    model <- stats::glm(os_status ~ older * sex, data = df,
                        family = stats::binomial)
    p <- suppressMessages(suppressWarnings(glmforest(model, data = df)))
    tbl <- attr(p, "table_data")

    expect_false(any(grepl(":", tbl$var, fixed = TRUE)))
    row <- interaction_row(p)
    expect_equal(nrow(row), 1L)
    expect_equal(row$estimate, unname(stats::coef(model)["olderTRUE:sexMale"]))
    expect_equal(row$N, sum(both))
    expect_equal(row$Events, sum(used$os_status[both]))
})


test_that("a variable named after a logical variable is drawn once", {

    ## The coefficients of the logical smoker were those whose names begin with
    ## "smoker", which included the continuous smoker_years
    df <- fresh_clintrial()
    df$smoker <- df$age > 60
    df$smoker_years <- seq_len(nrow(df)) %% 30

    model <- stats::glm(os_status ~ smoker + smoker_years, data = df,
                        family = stats::binomial)
    tbl <- attr(suppressMessages(suppressWarnings(glmforest(model, data = df))),
                "table_data")

    ## The logical is drawn by level, and smoker_years once
    expect_identical(tbl$var, c("smoker", "smoker", "smoker_years"))
    expect_identical(tbl$level, c("FALSE", "TRUE", "-"))
})
