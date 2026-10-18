## Coefficient names are matched against variable names by prefix. Matching in
## formula order allowed a short name to claim a coefficient belonging to a
## longer name that starts with the same characters, and matching by regular
## expression mishandled names containing metacharacters.

test_that("coefficients are assigned to the longest matching variable name", {

    ## Declared with the longer name first, which is the order in which the
    ## shorter name would previously have overwritten the correct assignment.
    xlevels <- list(
        sexuality = c("Hetero", "Other"),
        sex = c("Female", "Male")
    )

    parsed <- summata:::parse_term(c("sexualityOther", "sexMale"), xlevels)

    expect_identical(parsed$variable, c("sexuality", "sex"))
    expect_identical(parsed$group, c("Other", "Male"))
})


test_that("variable names containing regular expression metacharacters parse", {

    xlevels <- list(`factor(stage)` = c("I", "II", "III"))

    parsed <- summata:::parse_term(c("factor(stage)II", "factor(stage)III"),
                                   xlevels)

    expect_identical(parsed$variable, rep("factor(stage)", 2L))
    expect_identical(parsed$group, c("II", "III"))
})


test_that("interaction terms are left unparsed regardless of name overlap", {

    xlevels <- list(
        sex = c("Female", "Male"),
        sexuality = c("Hetero", "Other")
    )

    parsed <- summata:::parse_term("sexMale:sexualityOther", xlevels)

    expect_identical(parsed$variable, "sexMale:sexualityOther")
    expect_identical(parsed$group, "")
})


test_that("variable counting distinguishes names that share a prefix", {

    ## The shorter name is written first, which is the order in which it would
    ## previously have absorbed the coefficients of the longer name and
    ## understated the number of variables in the model.
    set.seed(20240101)
    df <- data.frame(
        y = stats::rbinom(200, 1, 0.5),
        sex = factor(sample(c("Female", "Male"), 200, replace = TRUE)),
        sexuality = factor(sample(c("Hetero", "Other"), 200, replace = TRUE))
    )

    model <- stats::glm(y ~ sex + sexuality, data = df,
                        family = stats::binomial)

    expect_identical(summata:::detect_model_type(model), "Multivariable")
})


test_that("prefix overlap does not disturb group sizes", {

    set.seed(20240102)
    df <- data.frame(
        y = stats::rbinom(200, 1, 0.5),
        sex = factor(sample(c("Female", "Male"), 200, replace = TRUE)),
        sexuality = factor(sample(c("Hetero", "Other"), 200, replace = TRUE))
    )

    result <- fit(
        data = df,
        outcome = "y",
        predictors = c("sex", "sexuality"),
        model_type = "glm",
        family = "binomial"
    )

    raw <- attr(result, "raw_data")

    for (pred in c("sex", "sexuality")) {
        sizes <- raw[["n"]][raw[["variable"]] == pred]
        expect_length(sizes, 2L)
        expect_equal(sum(sizes), nrow(df))
    }
})


test_that("interaction components are matched to the variable they name", {

    ## The shorter name is listed first, the order in which it previously
    ## claimed the component of the longer name
    components <- summata:::split_interaction_term(
        "sexMale:sexualityOther:age", c("sex", "sexuality", "age"),
        list(sex = c("Female", "Male"), sexuality = c("Hetero", "Other")))

    expect_identical(vapply(components, `[[`, "", "var"),
                     c("sex", "sexuality", "age"))
    expect_identical(vapply(components, `[[`, "", "level"),
                     c("Male", "Other", NA))
})


test_that("interaction components of non-syntactic names are matched", {

    components <- summata:::split_interaction_term(
        "`high bp`Yes:`has diabetes`Yes", c("high bp", "has diabetes"),
        list(`high bp` = c("No", "Yes"), `has diabetes` = c("No", "Yes")))

    expect_identical(vapply(components, `[[`, "", "var"),
                     c("high bp", "has diabetes"))
    expect_identical(vapply(components, `[[`, "", "level"), c("Yes", "Yes"))
})


test_that("declared levels decide between names of which one begins the other", {

    ## "abX" is level "bX" of a, which the longest name alone would assign to
    ## ab as an undeclared level "X"
    components <- summata:::split_interaction_term(
        "abX:abY", c("a", "ab"), list(a = c("bX", "c"), ab = c("Y", "Z")))

    expect_identical(vapply(components, `[[`, "", "var"), c("a", "ab"))
    expect_identical(vapply(components, `[[`, "", "level"), c("bX", "Y"))
})


test_that("coefficients of non-syntactic names are parsed into variable and level", {

    ## The backticks of "`high bp`Yes" prevented the match with the factor
    ## name, so the term was shown whole, without a reference row, and with
    ## the model sample size as its count
    df <- fresh_clintrial()
    df$`high bp` <- df$hypertension
    df$`has diabetes` <- df$diabetes

    model <- stats::glm(os_status ~ age + `high bp` * `has diabetes`, data = df,
                        family = stats::binomial)
    raw <- attr(fit(data = df, model = model), "raw_data")
    used <- df[stats::complete.cases(
        df[, c("os_status", "age", "high bp", "has diabetes")]), ]

    expect_identical(raw$group[raw$variable == "high bp"], c("No", "Yes"))
    expect_equal(raw$n[raw$variable == "high bp"],
                 as.numeric(table(used$`high bp`)))

    ## The interaction counts the observations in both levels
    interaction <- raw[grepl(":", raw$variable, fixed = TRUE)]
    expect_equal(interaction$n,
                 sum(used$`high bp` == "Yes" & used$`has diabetes` == "Yes"))
})


test_that("prefix overlap does not disturb interaction counts", {

    set.seed(20240103)
    df <- data.frame(
        y = stats::rbinom(300, 1, 0.5),
        sex = factor(sample(c("Female", "Male"), 300, replace = TRUE)),
        sexuality = factor(sample(c("Hetero", "Other"), 300, replace = TRUE))
    )
    both <- df$sex == "Male" & df$sexuality == "Other"

    result <- fit(data = df, outcome = "y", predictors = c("sex", "sexuality"),
                  interactions = "sex:sexuality", model_type = "glm")
    interaction <- attr(result, "raw_data")[term == "sexMale:sexualityOther"]

    expect_equal(interaction$n, sum(both))
    expect_equal(interaction$events, sum(df$y[both]))
})


test_that("a model of one non-syntactic factor is univariable", {

    df <- fresh_clintrial()
    df$`disease stage` <- df$stage
    model <- stats::glm(os_status ~ `disease stage`, data = df,
                        family = stats::binomial)

    expect_identical(summata:::detect_model_type(model), "Univariable")
})


test_that("a column named like a coefficient does not disturb interaction counts", {

    ## The column sexMale is not in the model, and its name must not claim the
    ## component sexMale of the interaction, as every column of the counting
    ## data once could
    df <- fresh_clintrial()
    df$sexMale <- as.numeric(df$sex == "Male")
    model <- survival::coxph(survival::Surv(os_months, os_status) ~ age +
                                 sex * treatment, data = df)
    raw <- attr(fit(data = df, model = model), "raw_data")
    used <- df[stats::complete.cases(df[, c("os_months", "os_status", "age",
                                            "sex", "treatment")]), ]
    both <- used$sex == "Male" & used$treatment == "Drug A"

    interaction <- raw[raw$term == "sexMale:treatmentDrug A"]
    expect_equal(interaction$n, sum(both))
    expect_equal(interaction$events, sum(used$os_status[both]))
})


test_that("a logical variable named after a factor is not parsed as its level", {

    ## The coefficient "graderTRUE" of the logical grader begins with the name
    ## of the factor grade, and was shown as its level "rTRUE"
    df <- fresh_clintrial()
    df$grader <- df$age > 65
    model <- stats::glm(os_status ~ grade + grader, data = df,
                        family = stats::binomial)
    raw <- attr(fit(data = df, model = model), "raw_data")

    expect_identical(raw$group[raw$variable == "grade"], levels(df$grade))
    expect_identical(raw$group[raw$variable == "grader"], c("FALSE", "TRUE"))
    expect_identical(summata:::detect_model_type(model), "Multivariable")
})


test_that("a continuous variable named after a factor is not parsed as its level", {

    ## "smoking_years" begins with the factor name "smoking"; the table showed
    ## it as level "_years" of smoking, and its interactions with no
    ## observations
    df <- fresh_clintrial()
    df$smoking_years <- seq_len(nrow(df)) %% 41

    model <- stats::glm(os_status ~ smoking + treatment * smoking_years,
                        data = df, family = stats::binomial)
    raw <- attr(fit(data = df, model = model), "raw_data")
    used <- df[stats::complete.cases(
        df[, c("os_status", "smoking", "treatment", "smoking_years")]), ]

    expect_identical(raw$group[raw$variable == "smoking"],
                     c("Never", "Former", "Current"))
    expect_true("smoking_years" %in% raw$variable)

    ## The interaction with a continuous variable counts the observations in
    ## the level of the factor, as the forest plots do
    interaction <- raw[raw$term == "treatmentDrug A:smoking_years"]
    expect_equal(interaction$n, sum(used$treatment == "Drug A"))
    expect_equal(interaction$events,
                 sum(used$os_status[used$treatment == "Drug A"]))

    expect_identical(summata:::detect_model_type(
        stats::glm(os_status ~ smoking + smoking_years, data = df,
                   family = stats::binomial)), "Multivariable")
})


test_that("reference terms are named as R names coefficients", {

    df <- fresh_clintrial()
    df$`high bp` <- df$hypertension
    model <- stats::glm(os_status ~ age + `high bp`, data = df,
                        family = stats::binomial)
    raw <- attr(fit(data = df, model = model), "raw_data")

    expect_identical(raw$term[raw$variable == "high bp"],
                     c("`high bp`No", "`high bp`Yes"))
})


test_that("mixed models count the variables of their fixed-effect terms", {

    ## The formula text was split at "+", which kept "hypertension * diabetes"
    ## as one term, and its random effects were removed by a pattern that
    ## reduced "log(age) + sex + (1 | site)" to "log"
    skip_if_not_installed("lme4")

    df <- fresh_clintrial()
    scope <- function(formula) {
        model <- suppressMessages(lme4::lmer(formula, data = df))
        summata:::detect_model_type(model)
    }

    expect_identical(scope(bmi ~ hypertension * diabetes + (1 | site)), "Multivariable")
    expect_identical(scope(bmi ~ log(age) + sex + (1 | site)), "Multivariable")
    expect_identical(scope(bmi ~ sex + (1 | site)), "Univariable")
    expect_identical(scope(bmi ~ stage + (1 | site) + (1 | treatment)), "Univariable")
})


test_that("coxme models fitted by fit() count their variables from their terms", {

    ## The formula of the call of such a model is only a name, from which no
    ## factor was found, so that each level of stage counted as a variable
    skip_if_not_installed("coxme")

    result <- fit(data = fresh_clintrial(), outcome = "Surv(os_months, os_status)",
                  predictors = "stage", model_type = "coxme", random = "(1|site)")

    expect_identical(summata:::detect_model_type(attr(result, "model")),
                     "Univariable")
    expect_true("HR (95% CI)" %in% names(result))
})
