## tablesave() and forestsave() are counterparts and validate, report, and
## return on the same terms.

tablesave_fixture <- function() {
    dt <- data.table::as.data.table(fresh_clintrial())
    desctable(dt, by = "treatment", variables = c("age", "sex"))
}


test_that("tablesave writes delimited formats without further dependencies", {

    tbl <- tablesave_fixture()

    csv <- tempfile(fileext = ".csv")
    expect_message(tablesave(tbl, csv), "Table saved to")
    expect_true(file.exists(csv))
    expect_gt(nrow(data.table::fread(csv)), 0)

    tsv <- tempfile(fileext = ".tsv")
    tablesave(tbl, tsv, quiet = TRUE)
    expect_true(file.exists(tsv))

    ## The tab separator is applied, not merely the extension
    expect_true(any(grepl("\t", readLines(tsv, n = 1), fixed = TRUE)))
})


test_that("tablesave returns its path invisibly and can be silenced", {

    tbl <- tablesave_fixture()
    out <- tempfile(fileext = ".csv")

    expect_no_message(result <- tablesave(tbl, out, quiet = TRUE))
    expect_identical(result, out)
})


test_that("tablesave validates its arguments as forestsave does", {

    tbl <- tablesave_fixture()

    expect_error(tablesave("not a table", tempfile(fileext = ".csv")),
                 "must be a data frame")
    expect_error(tablesave(tbl, character(0)), "non-empty character")
    expect_error(tablesave(tbl, ""), "non-empty character")
    expect_error(tablesave(tbl, tempfile()), "no extension")
    expect_error(tablesave(tbl, tempfile(fileext = ".xyz")),
                 "Unsupported file extension")
})


test_that("the two save functions agree on their contract", {

    tbl <- tablesave_fixture()
    dt <- data.table::as.data.table(fresh_clintrial())
    model <- stats::glm(readmission_30d ~ age + sex,
                        data = dt, family = stats::binomial)
    p <- suppressMessages(suppressWarnings(glmforest(model, data = dt)))

    ## Both are quiet on request, and both return the path they were given
    t_out <- tempfile(fileext = ".csv")
    f_out <- tempfile(fileext = ".pdf")

    expect_no_message(t_result <- tablesave(tbl, t_out, quiet = TRUE))
    expect_no_message(f_result <- forestsave(p, f_out, quiet = TRUE))

    expect_identical(t_result, t_out)
    expect_identical(f_result, f_out)

    ## Both name the file in their message
    expect_message(tablesave(tbl, tempfile(fileext = ".csv")), "saved to")
    expect_message(forestsave(p, tempfile(fileext = ".pdf")), "saved to")
})
