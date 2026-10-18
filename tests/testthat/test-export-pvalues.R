## The export functions bold significant p-values by reading them back from
## the formatted table. They must recognize p-values in every number format
## and at every p_digits setting, not only "< 0.001" and decimal points.

test_that("is_significant_pvalue reads formatted p-values in any locale", {

    x <- c("0.023", "0,023", "0.051", "0,051", "< 0.001", "< 0,001", "< 0.0001",
           "< 0.01", "1.000", "-", "", NA, "reference")
    expect_identical(summata:::is_significant_pvalue(x, 0.05),
                     c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE,
                       TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))

    ## A value reported only as below a threshold is significant when the
    ## threshold does not exceed the significance level
    expect_identical(summata:::is_significant_pvalue(c("< 0.01", "< 0.001"), 0.005),
                     c(FALSE, TRUE))
})


test_that("HTML and LaTeX exports bold significant p-values in every format", {

    skip_if_not_installed("xtable")
    df <- fresh_clintrial()
    out <- tempfile(fileext = ".html")
    tex <- tempfile(fileext = ".tex")

    for (args in list(list(number_format = "eu"), list(p_digits = 4), list(p_digits = 2))) {
        result <- do.call(fit, c(list(data = df, outcome = "os_status",
                                      predictors = c("age", "stage")), args))
        n_significant <- sum(summata:::is_significant_pvalue(result$`p-value`, 0.05))
        expect_gt(n_significant, 0)

        table2html(result, out, quiet = TRUE)
        expect_equal(sum(grepl("<b>", readLines(out), fixed = TRUE)), n_significant)

        table2tex(result, tex, quiet = TRUE)
        expect_equal(sum(grepl("\\textbf{", readLines(tex), fixed = TRUE)), n_significant)
    }
})
