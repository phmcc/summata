## format_num() formats the summary statistics of desctable(). Values of 1000
## or more receive a thousands mark, and must keep the requested number of
## decimal places and the decimal mark of the locale.

test_that("format_num keeps the decimal places of values of 1000 or more", {

    us <- summata:::resolve_number_marks("us")

    expect_identical(summata:::format_num(1234.5, "%.2f", us), "1,234.50")
    expect_identical(summata:::format_num(1234.4, "%.0f", us), "1,234")
    expect_identical(summata:::format_num(-1234.567, "%.1f", us), "-1,234.6")
    expect_identical(summata:::format_num(1000000.25, "%.2f", us), "1,000,000.25")
    expect_identical(summata:::format_num(999.04, "%.1f", us), "999.0")

    ## A value that rounds to 1000 is marked like one that is 1000
    expect_identical(summata:::format_num(999.96, "%.1f", us), "1,000.0")
})


test_that("format_num applies both marks of every locale", {

    x <- 1234.5

    expect_identical(summata:::format_num(x, "%.1f", summata:::resolve_number_marks("eu")),
                     "1.234,5")
    expect_identical(summata:::format_num(x, "%.1f", summata:::resolve_number_marks("space")),
                     "1 234.5")
    expect_identical(summata:::format_num(x, "%.1f", summata:::resolve_number_marks("none")),
                     "1234.5")
    expect_identical(summata:::format_num(x, "%.2f", summata:::resolve_number_marks(c("'", "."))),
                     "1'234.50")
})


test_that("desctable reports large values with the requested digits", {

    df <- fresh_clintrial()
    df$cost <- df$age * 1000.37

    result <- desctable(df, by = "treatment", variables = "cost",
                        stats_continuous = "mean_sd", digits = 2, number_format = "eu")
    expected_mean <- formatC(mean(df$cost), format = "f", digits = 2,
                             big.mark = ".", decimal.mark = ",")
    expect_true(startsWith(result$Total[2], expected_mean))
})
