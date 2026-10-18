## desctable() reports the median survival time with a confidence interval at
## the requested level. survfit() names the interval columns after the level
## (for example "0.9LCL"), and the interval must be found at every level, not
## only those written with two decimals.

test_that("the median survival interval is reported at any confidence level", {

    df <- fresh_clintrial()

    for (level in c(0.8, 0.9, 0.95, 0.99)) {
        result <- desctable(df, by = "treatment", variables = "Surv(os_months, os_status)",
                            conf_level = level)
        raw <- attr(result, "raw_data")

        fit <- survival::survfit(survival::Surv(os_months, os_status) ~ 1, data = df,
                                 conf.int = level)
        tab <- summary(fit)$table
        expect_equal(unname(raw$Total_ci_lower), unname(tab[paste0(level, "LCL")]))
        expect_equal(unname(raw$Total_ci_upper), unname(tab[paste0(level, "UCL")]))

        ## No group interval is missing
        expect_false(any(grepl("NA", unlist(result[2, -1]))))
        expect_match(result$Group[2], paste0(round(level * 100), "% CI"), fixed = TRUE)
    }
})
