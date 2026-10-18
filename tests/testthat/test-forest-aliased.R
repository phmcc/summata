## A coefficient that a model cannot estimate because its column of the model
## matrix is a linear combination of others is aliased: lm(), glm() and coxph()
## report it as NA, and lme4 drops it before fitting. The forest plots drew each
## such row as a reference category, and the Cox regression tables showed it as
## an empty row beside a reference row. Both now omit it, as the tables of the
## other models did, and a factor none of whose levels was estimated is omitted
## with its reference row.

## age2 is aliased with age, every level of sex2 with sex, and level "c" of
## grp3 with sex == "Male", so that grp3b is estimated and grp3c is not
aliased_data <- function() {
    df <- fresh_clintrial()
    df$age2 <- 2 * df$age
    df$sex2 <- df$sex
    df$grp3 <- factor(ifelse(df$sex == "Male", "c",
                             rep(c("a", "b"), length.out = nrow(df))))
    df
}

aliased_predictors <- c("age", "age2", "sex", "sex2", "grp3")

forest_rows <- function(p) {
    tbl <- attr(p, "table_data")
    tbl[, c("var", "level", "estimate", "conf_low")]
}


test_that("forest plots omit aliased coefficients", {

    df <- aliased_data()
    rhs <- paste(aliased_predictors, collapse = " + ")

    models <- list(
        glm = stats::glm(stats::as.formula(paste("os_status ~", rhs)), data = df,
                         family = stats::binomial),
        lm = stats::lm(stats::as.formula(paste("bmi ~", rhs)), data = df),
        cox = survival::coxph(stats::as.formula(
            paste("survival::Surv(os_months, os_status) ~", rhs)), data = df)
    )
    plotters <- list(glm = glmforest, lm = lmforest, cox = coxforest)

    for (type in names(models)) {
        model <- models[[type]]
        expect_identical(names(which(is.na(stats::coef(model)))),
                         c("age2", "sex2Male", "grp3c"), info = type)

        rows <- forest_rows(suppressWarnings(plotters[[type]](model, data = df)))

        ## Neither the aliased variables nor the aliased level remain
        expect_identical(unique(rows$var), c("age", "sex", "grp3"), info = type)
        expect_identical(rows$level[rows$var == "grp3"], c("a", "b"), info = type)

        ## Every remaining row is either a reference category or an estimate
        coefs <- stats::coef(model)
        expect_equal(rows$estimate[rows$var == "grp3"],
                     c(0, unname(coefs["grp3b"])), info = type)
        expect_equal(sum(is.na(rows$conf_low)), 2L, info = type)
    }
})


test_that("condensed and indented layouts omit aliased coefficients", {

    df <- aliased_data()
    model <- stats::glm(os_status ~ age + age2 + sex + sex2 + grp3, data = df,
                        family = stats::binomial)

    for (layout in list(c(TRUE, FALSE), c(TRUE, TRUE))) {
        p <- suppressWarnings(glmforest(model, data = df, indent_groups = layout[1],
                                        condense_table = layout[2]))
        labels <- attr(p, "table_data")$var
        info <- paste("indent", layout[1], "condense", layout[2])

        expect_false(any(grepl("age2|sex2", labels)), info = info)
        expect_false(any(trimws(labels) == "c"), info = info)
    }
})


test_that("mixed-model forest plots omit the coefficients lme4 drops", {

    skip_if_not_installed("lme4")

    df <- aliased_data()

    glmer_model <- suppressMessages(lme4::glmer(
        os_status ~ age + age2 + sex + sex2 + grp3 + (1 | site), data = df,
        family = stats::binomial))
    lmer_model <- suppressMessages(lme4::lmer(
        bmi ~ age + age2 + sex + sex2 + grp3 + (1 | site), data = df))

    for (model in list(glmer_model, lmer_model)) {
        dropped <- names(attr(lme4::getME(model, "X"), "col.dropped"))
        expect_setequal(dropped, c("age2", "sex2Male", "grp3c"))

        plotter <- if (inherits(model, "glmerMod")) glmforest else lmforest
        rows <- forest_rows(suppressWarnings(plotter(model, data = df)))

        expect_identical(unique(rows$var), c("age", "sex", "grp3"))
        expect_identical(rows$level[rows$var == "grp3"], c("a", "b"))
    }
})


test_that("Cox regression tables omit aliased coefficients as other tables do", {

    df <- aliased_data()

    tables <- list(
        glm = fit(data = df, outcome = "os_status", predictors = aliased_predictors,
                  model_type = "glm"),
        cox = fit(data = df, outcome = "Surv(os_months, os_status)",
                  predictors = aliased_predictors, model_type = "coxph")
    )

    for (type in names(tables)) {
        result <- tables[[type]]
        raw <- attr(result, "raw_data")

        expect_false(anyNA(raw$coefficient), info = type)
        expect_false(any(c("age2", "sex2") %in% raw$variable), info = type)
        expect_identical(raw$group[raw$variable == "grp3"], c("a", "b"), info = type)
    }

    ## Both tables show the same rows
    expect_identical(tables$cox$Variable, tables$glm$Variable)
    expect_identical(tables$cox$Group, tables$glm$Group)
})


test_that("a factor whose first level is empty is shown against its effective reference", {

    ## coxph() keeps the empty level I, so that the indicators of II, III and
    ## IV sum to one: stageIV is aliased and II and III are estimated against
    ## IV, which the table and plot showed as an empty row or omitted while
    ## labeling the empty level I as the reference
    df <- subset(fresh_clintrial(), stage != "I")
    model <- survival::coxph(survival::Surv(os_months, os_status) ~ age + stage,
                             data = df)
    expect_true(is.na(stats::coef(model)["stageIV"]))

    raw <- attr(fit(data = df, model = model), "raw_data")
    stage <- raw[raw$variable == "stage"]
    expect_identical(stage$group, c("IV", "II", "III"))
    expect_identical(stage$reference, c("reference", "", ""))
    expect_equal(stage$n[1], sum(df$stage == "IV"))

    for (indent in c(FALSE, TRUE)) {
        tbl <- attr(suppressWarnings(coxforest(model, data = df,
                                               indent_groups = indent)),
                    "table_data")
        levels_drawn <- if (indent) trimws(tbl$var[-(1:2)]) else tbl$level[-1]
        expect_identical(levels_drawn, c("IV", "II", "III"), info = indent)
        expect_equal(tbl$estimate[!is.na(tbl$conf_low) & tbl$var != "age"],
                     unname(stats::coef(model)[c("stageII", "stageIII")]),
                     info = indent)
    }
})


## Levels II, III and IV observed and I empty; a covariate that indicates one
## level makes that level aliased as well
several_aliased_data <- function(indicated) {
    df <- subset(fresh_clintrial(), stage != "I")
    df$indicator <- as.numeric(df$stage == indicated)
    df
}


test_that("a factor with several aliased levels is shown against its effective reference", {

    ## With the indicator of IV before stage, stageIII is aliased by the
    ## constant and stageIV by the indicator. II is estimated against III,
    ## whose rows the indicator leaves at zero; the table and plot labeled
    ## the empty level I as the reference
    df <- several_aliased_data("IV")
    model <- survival::coxph(survival::Surv(os_months, os_status) ~ age +
                                 indicator + stage, data = df)
    expect_true(all(is.na(stats::coef(model)[c("stageIII", "stageIV")])))

    df_ref <- df
    df_ref$stage <- stats::relevel(droplevels(df_ref$stage), "III")
    refit <- survival::coxph(survival::Surv(os_months, os_status) ~ age +
                                 indicator + stage, data = df_ref)

    raw <- attr(fit(data = df, model = model), "raw_data")
    stage <- raw[raw$variable == "stage"]
    expect_identical(stage$group, c("III", "II"))
    expect_identical(stage$reference, c("reference", ""))
    expect_equal(stage$n[1], sum(df$stage == "III"))
    expect_equal(stage$coefficient[2], unname(stats::coef(refit)["stageII"]))

    tbl <- attr(suppressWarnings(coxforest(model, data = df)), "table_data")
    expect_identical(tbl$level[tbl$var == "stage"], c("III", "II"))
})


test_that("the effective reference may be the last of the aliased levels", {

    ## With the indicator of II, stageII is aliased by the indicator and
    ## stageIV by the constant: III is estimated against IV
    df <- several_aliased_data("II")
    model <- survival::coxph(survival::Surv(os_months, os_status) ~ indicator +
                                 stage, data = df)
    expect_true(all(is.na(stats::coef(model)[c("stageII", "stageIV")])))

    raw <- m2dt(df, model)
    stage <- raw[raw$variable == "stage"]
    expect_identical(stage$group, c("IV", "III"))
    expect_identical(stage$reference, c("reference", ""))

    ## The strata of a stratified model take the place of the constant
    stratified <- survival::coxph(survival::Surv(os_months, os_status) ~
                                      indicator + stage + survival::strata(site),
                                  data = df)
    raw_strat <- m2dt(df, stratified)
    expect_identical(raw_strat$group[raw_strat$variable == "stage"], c("IV", "III"))
})


test_that("a factor estimated against no single level is shown without a reference", {

    ## pair has the empty first level A and the levels k1, r1, k2 and r2. The
    ## indicator of "k1 or r1" leaves r1 aliased by the indicator and r2 by
    ## the constant, so that k1 is estimated against r1 and k2 against r2:
    ## no level is the reference of every estimate
    df <- subset(fresh_clintrial(), stage != "I")
    df$pair <- factor(c(II = "k1", III = "r1", IV = "k2")[as.character(df$stage)],
                      levels = c("A", "k1", "r1", "k2", "r2"))
    df$pair[df$stage == "IV" & seq_len(nrow(df)) %% 2 == 0] <- "r2"
    df$indicator <- as.numeric(df$pair %in% c("k1", "r1"))

    model <- survival::coxph(survival::Surv(os_months, os_status) ~ indicator +
                                 pair, data = df)
    expect_true(all(is.na(stats::coef(model)[c("pairr1", "pairr2")])))

    expect_warning(raw <- m2dt(df, model), "not estimated against a single")
    expect_identical(raw$group[raw$variable == "pair"], c("k1", "k2"))
    expect_false(any(nzchar(raw$reference[raw$variable == "pair"])))

    expect_warning(tbl <- attr(suppressMessages(coxforest(model, data = df)),
                               "table_data"),
                   "not estimated against a single")
    expect_identical(tbl$level[tbl$var == "pair"], c("k1", "k2"))
})


test_that("the effective reference is found for a factor in an interaction", {

    ## The columns of the interactions of stage are built from its levels
    ## and are set aside with its own, so that III is the reference as it is
    ## without the interaction; no reference was found for such a factor
    df <- several_aliased_data("IV")
    model <- survival::coxph(survival::Surv(os_months, os_status) ~ indicator +
                                 stage * sex, data = df)
    expect_true(all(is.na(stats::coef(model)[c("stageIII", "stageIV")])))

    expect_no_warning(raw <- m2dt(df, model))
    stage <- raw[raw$variable == "stage"]
    expect_identical(stage$group, c("III", "II"))
    expect_identical(stage$reference, c("reference", ""))

    tbl <- attr(suppressMessages(coxforest(model, data = df)), "table_data")
    expect_identical(tbl$level[tbl$var == "stage"], c("III", "II"))
})


test_that("no reference is shown when the model matrix cannot be rebuilt", {

    ## The model matrix is rebuilt from the data of the model's call, which
    ## no longer hold the observations analyzed
    df <- several_aliased_data("IV")
    model_data <- df
    model <- survival::coxph(survival::Surv(os_months, os_status) ~ indicator +
                                 stage, data = model_data)
    model_data <- model_data[1:100, ]

    expect_null(summata:::determined_levels(model, "stage", levels(df$stage)))
    expect_warning(raw <- m2dt(df, model), "could not be determined")
    expect_identical(raw$group[raw$variable == "stage"], "II")
    expect_false(any(nzchar(raw$reference[raw$variable == "stage"])))
})


test_that("a model with no estimable predictor stops with the cause", {

    df <- fresh_clintrial()
    df$const <- 1

    expect_error(fit(data = df, outcome = "os_status", predictors = "const",
                     model_type = "glm"),
                 "None of the predictors of the model could be estimated")
    expect_error(fit(data = df, outcome = "Surv(os_months, os_status)",
                     predictors = "const", model_type = "coxph"),
                 "None of the predictors of the model could be estimated")
})


test_that("forest plots of a model with no estimable predictor stop with the cause", {

    df <- fresh_clintrial()
    df$const <- 1
    message <- "None of the predictors of the model could be estimated"

    expect_error(suppressWarnings(glmforest(
        stats::glm(os_status ~ const, data = df, family = stats::binomial),
        data = df)), message)
    expect_error(suppressWarnings(lmforest(
        stats::lm(bmi ~ const, data = df), data = df)), message)
    expect_error(suppressWarnings(coxforest(
        survival::coxph(survival::Surv(os_months, os_status) ~ const, data = df),
        data = df)), message)
})
