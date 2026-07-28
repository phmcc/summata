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
