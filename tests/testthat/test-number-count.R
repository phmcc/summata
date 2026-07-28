## format_count() is the package's single count formatter. It is vectorized,
## applies the grouping mark only from four digits, and must not fall back to
## exponential notation for large counts.

test_that("format_count is vectorized and formats each element independently", {

    marks <- summata:::resolve_number_marks("us")

    result <- summata:::format_count(c(5, 999, 1000, 20000), marks)

    expect_length(result, 4L)
    ## No padding to a common width
    expect_identical(result, c("5", "999", "1,000", "20,000"))
})


test_that("format_count does not fall back to exponential notation", {

    marks <- summata:::resolve_number_marks("us")

    expect_identical(summata:::format_count(100000, marks), "100,000")
    expect_identical(summata:::format_count(1000000, marks), "1,000,000")
})


test_that("format_count follows the locale and resolves it when absent", {

    expect_identical(
        summata:::format_count(20000, summata:::resolve_number_marks("eu")),
        "20.000")

    withr::local_options(summata.number_format = "eu")
    expect_identical(summata:::format_count(20000), "20.000")
})


test_that("format_count preserves missing values", {

    marks <- summata:::resolve_number_marks("us")

    result <- summata:::format_count(c(1000, NA, 5), marks)

    expect_true(is.na(result[2]))
    expect_identical(result[c(1, 3)], c("1,000", "5"))
})


test_that("format_count emits no marks warning under any locale", {

    for (fmt in c("us", "eu")) {
        marks <- summata:::resolve_number_marks(fmt)
        expect_no_warning(summata:::format_count(c(999, 1000, 250000), marks))
    }
})


test_that("format_count substitutes missing values per its caller's convention", {

    marks <- summata:::resolve_number_marks("us")
    counts <- c(1000, NA, 5)

    ## Tables show a dash, plot annotations are left blank, and the default
    ## preserves the missing value
    expect_identical(summata:::format_count(counts, marks, na = "-"),
                     c("1,000", "-", "5"))
    expect_identical(summata:::format_count(counts, marks, na = ""),
                     c("1,000", "", "5"))
    expect_true(is.na(summata:::format_count(counts, marks)[2]))
})
