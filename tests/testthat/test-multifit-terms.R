## multifit() reports the coefficients of its predictor: the predictor's own
## term and the interactions of which the predictor is a component. Other
## terms of the adjusted model are not reported, however their names relate
## to the predictor's.


### * Fixtures

terms_data <- function() {
    dt <- data.table::as.data.table(fresh_clintrial())
    dt[, age_group := cut(age, c(0, 50, 65, Inf),
                          labels = c("Young", "Middle", "Old"))]
    dt[, smoker := smoking == "Current"]
    dt[]
}


### * Covariates and interactions of other variables

test_that("a covariate whose name begins with the predictor's is not reported", {

    dt <- terms_data()
    result <- multifit(dt, outcomes = "readmission_30d", predictor = "age",
                       covariates = c("age_group", "sex"), parallel = FALSE)
    raw <- attr(result, "raw_data")

    model <- stats::glm(readmission_30d ~ age + age_group + sex, data = dt,
                        family = stats::binomial)

    expect_identical(raw$group, "-")
    expect_equal(raw$coefficient, unname(stats::coef(model)["age"]))
})


test_that("an interaction of other variables is not reported under the predictor", {

    dt <- terms_data()

    ## "stage" contains "age", but the interaction does not involve age
    result <- multifit(dt, outcomes = "readmission_30d", predictor = "age",
                       covariates = c("stage", "sex"), interactions = "stage:sex",
                       parallel = FALSE)
    raw <- attr(result, "raw_data")

    expect_identical(raw$group, "-")
})


### * The predictor's own terms

test_that("each level of a factor predictor is reported once", {

    dt <- terms_data()
    dt[, treatment_cycles := ifelse(age > 60, 4, 6)]

    result <- multifit(dt, outcomes = "readmission_30d", predictor = "treatment",
                       covariates = c("treatment_cycles", "age"), parallel = FALSE)
    raw <- attr(result, "raw_data")

    model <- stats::glm(readmission_30d ~ treatment + treatment_cycles + age,
                        data = dt, family = stats::binomial)
    estimates <- stats::coef(model)[c("treatmentDrug A", "treatmentDrug B")]

    expect_identical(raw$group, c("Drug A", "Drug B"))
    expect_equal(raw$coefficient, unname(estimates))
})


test_that("a logical predictor is reported by its TRUE coefficient", {

    dt <- terms_data()
    result <- multifit(dt, outcomes = "readmission_30d", predictor = "smoker",
                       covariates = c("age", "sex"), parallel = FALSE)
    raw <- attr(result, "raw_data")

    model <- stats::glm(readmission_30d ~ smoker + age + sex, data = dt,
                        family = stats::binomial)

    expect_identical(raw$group, "TRUE")
    expect_equal(raw$coefficient, unname(stats::coef(model)["smokerTRUE"]))
})


test_that("interactions that involve the predictor are reported in either order", {

    dt <- terms_data()
    model <- stats::glm(readmission_30d ~ treatment + age + sex + treatment:sex,
                        data = dt, family = stats::binomial)
    interaction_terms <- grep(":", names(stats::coef(model)), fixed = TRUE,
                              value = TRUE)

    for (int in c("treatment:sex", "sex:treatment")) {
        result <- multifit(dt, outcomes = "readmission_30d",
                           predictor = "treatment", covariates = c("age", "sex"),
                           interactions = int, parallel = FALSE)
        raw <- attr(result, "raw_data")
        rows <- grepl(":", raw$group, fixed = TRUE)

        expect_identical(raw$group[!rows], c("Drug A", "Drug B"), info = int)
        expect_length(raw$group[rows], length(interaction_terms))
        expect_equal(sort(raw$coefficient[rows]),
                     sort(unname(stats::coef(model)[interaction_terms])),
                     info = int)
    }
})
