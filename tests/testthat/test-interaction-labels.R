## The labels of a regression table name the variables of its interaction
## terms as well: each component of an interaction is matched to the variable
## it names, and only the components of labeled variables are relabeled.


### * Fixtures

label_data <- function() {
    dt <- data.table::as.data.table(fresh_clintrial())
    set.seed(20261016)
    dt[, sexuality := factor(sample(c("Hetero", "Other"), .N, replace = TRUE))]
    dt[, shift := factor(sample(c("08:00", "12:00", "18:00"), .N, replace = TRUE),
                         levels = c("08:00", "12:00", "18:00"))]
    dt[]
}


## The Variable entries of the interaction rows of a fit() table
interaction_labels <- function(result) {
    raw <- attr(result, "raw_data")
    result$Variable[grepl(":", raw$variable, fixed = TRUE)]
}


### * Components

test_that("a component is labeled by the variable it names", {

    ## "sexualityOther" begins with the labeled "sex", but names the
    ## unlabeled sexuality, whose component keeps its name
    dt <- label_data()
    result <- fit(dt, "readmission_30d", c("treatment", "sexuality", "sex"),
                  interactions = "treatment:sexuality", model_type = "glm",
                  labels = c(sex = "Sex", treatment = "Treatment"))

    expect_identical(interaction_labels(result),
                     c("Treatment (Drug A) × sexualityOther",
                       "Treatment (Drug B) × sexualityOther"))
})


test_that("a level that contains a colon stays in its component", {

    dt <- label_data()
    labels <- c(shift = "Shift start", sex = "Sex")
    first <- fit(dt, "readmission_30d", c("shift", "sex"),
                 interactions = "shift:sex", model_type = "glm", labels = labels)
    second <- fit(dt, "readmission_30d", c("sex", "shift"),
                  interactions = "sex:shift", model_type = "glm", labels = labels)

    expect_identical(interaction_labels(first),
                     c("Shift start (12:00) × Sex (Male)",
                       "Shift start (18:00) × Sex (Male)"))
    expect_identical(interaction_labels(second),
                     c("Sex (Male) × Shift start (12:00)",
                       "Sex (Male) × Shift start (18:00)"))
})


test_that("a label that contains a colon is not taken for an interaction", {

    dt <- label_data()
    result <- fit(dt, "readmission_30d", c("shift", "sex"), model_type = "glm",
                  labels = c(shift = "Shift (hh:mm)", sex = "Sex"))

    expect_identical(result$Variable[nzchar(result$Variable)],
                     c("Shift (hh:mm)", "Sex"))
})


test_that("interactions of labeled variables are labeled as before", {

    dt <- label_data()
    result <- fit(dt, "readmission_30d", c("treatment", "sex"),
                  interactions = "treatment:sex", model_type = "glm",
                  labels = c(sex = "Sex", treatment = "Treatment"))

    expect_identical(interaction_labels(result),
                     c("Treatment (Drug A) × Sex (Male)",
                       "Treatment (Drug B) × Sex (Male)"))
})
