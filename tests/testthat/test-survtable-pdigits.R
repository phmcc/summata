## survtable() documents p_digits as the number of decimal places of its
## p-values; it must reach the p-value column for one outcome and several,
## and the test line of the printed summary.

test_that("survtable formats p-values with p_digits", {

    df <- fresh_clintrial()

    expected <- survival::survdiff(survival::Surv(os_months, os_status) ~ sex,
                                   data = df)$pvalue
    for (p_digits in c(2, 3, 5)) {
        result <- survtable(df, "Surv(os_months, os_status)", by = "sex",
                            times = 12, p_digits = p_digits)
        p <- result$`p-value`[result$`p-value` != ""]
        if (expected < 10^-p_digits) {
            expect_identical(p, paste0("< 0.", strrep("0", p_digits - 1), "1"))
        } else {
            expect_identical(p, sprintf(paste0("%.", p_digits, "f"), expected))
        }
    }

    result <- survtable(df, c("Surv(os_months, os_status)", "Surv(pfs_months, pfs_status)"),
                        by = "stage", times = 12, p_digits = 5)
    p <- result$`p-value`[result$`p-value` != ""]
    expect_true(all(grepl("^< 0\\.00001$|^0\\.[0-9]{5}$", p)))
})


test_that("the printed test line formats the p-value with p_digits", {

    df <- fresh_clintrial()

    for (p_digits in c(2, 5)) {
        result <- survtable(df, "Surv(os_months, os_status)", by = "sex",
                            times = 12, p_digits = p_digits)
        p <- result$`p-value`[result$`p-value` != ""]
        line <- if (startsWith(p, "<")) paste0("(p ", p, ")") else paste0("(p = ", p, ")")

        printed <- utils::capture.output(print(result))
        expect_true(any(grepl(paste0("Test: Log-rank ", line), printed, fixed = TRUE)),
                    info = paste("p_digits =", p_digits))
    }
})
