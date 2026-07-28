## Labels supplied to a fitting function travel with its result, so that the
## forest plot functions apply the same labels without being given them again.

label_fixture <- function() {
    env <- new.env(parent = emptyenv())
    utils::data("clintrial_labels", package = "summata", envir = env)
    env$clintrial_labels
}


test_that("fitting functions attach the labels they were given", {

    dt <- data.table::as.data.table(fresh_clintrial())
    labs <- label_fixture()

    fitted <- fit(dt, "readmission_30d", c("age", "sex"),
                  model_type = "glm", family = "binomial", labels = labs)
    screened <- uniscreen(dt, "readmission_30d", c("age", "sex"),
                          model_type = "glm", family = "binomial",
                          labels = labs, parallel = FALSE)

    expect_identical(attr(fitted, "labels"), labs)
    expect_identical(attr(screened, "labels"), labs)

    ## A result built without labels carries none
    plain <- fit(dt, "readmission_30d", c("age", "sex"),
                 model_type = "glm", family = "binomial")
    expect_null(attr(plain, "labels"))
})


test_that("forest plots inherit labels from the result they are given", {

    dt <- data.table::as.data.table(fresh_clintrial())
    labs <- label_fixture()

    screened <- uniscreen(dt, "readmission_30d", c("age", "sex"),
                          model_type = "glm", family = "binomial",
                          labels = labs, parallel = FALSE)

    inherited <- suppressMessages(suppressWarnings(uniforest(screened)))
    explicit <- suppressMessages(suppressWarnings(
        uniforest(screened, labels = labs)))

    ## The two are equivalent, which is the point of the inheritance
    expect_identical(attr(inherited, "table_data"),
                     attr(explicit, "table_data"))
})


test_that("an explicit labels argument overrides the inherited one", {

    dt <- data.table::as.data.table(fresh_clintrial())
    labs <- label_fixture()
    override <- c(age = "AGE OVERRIDE", sex = "SEX OVERRIDE")

    screened <- uniscreen(dt, "readmission_30d", c("age", "sex"),
                          model_type = "glm", family = "binomial",
                          labels = labs, parallel = FALSE)

    p <- suppressMessages(suppressWarnings(
        uniforest(screened, labels = override)))

    tbl <- attr(p, "table_data")
    skip_if(is.null(tbl), "uniforest does not expose table_data")

    expect_true(any(grepl("AGE OVERRIDE", unlist(tbl), fixed = TRUE)))
})
