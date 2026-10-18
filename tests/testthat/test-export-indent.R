## Indented exports place the levels of a factor under the variable name.
## The p-values of the levels of a regression table belong to those levels
## and must survive the indentation, whatever the effect measure and
## confidence level in the column name. Tables whose interval columns name a
## statistic rather than an effect measure, as in survtable(), are not
## regression tables.

level_pvalues <- function(table, ...) {
    out <- summata:::format_indented_groups(data.table::as.data.table(table), ...)
    out[["p-value"]][grepl("^\\s+", out$Variable)]
}


test_that("indented regression tables keep the p-values of factor levels", {

    df <- fresh_clintrial()

    for (conf_level in c(0.95, 0.9)) {
        result <- fit(df, "readmission_30d", c("age", "sex", "stage"),
                      conf_level = conf_level)
        ## Multivariable effect column, e.g. "aOR (95% CI)" or "aOR (90% CI)"
        expect_true(any(grepl("^aOR \\(", names(result))))

        expected <- result$`p-value`[result$Group %in% c("Female", "Male", "I", "II", "III", "IV")]
        expect_identical(level_pvalues(result), expected)
    }
})


test_that("indented descriptive tables report the p-value on the variable row", {

    df <- fresh_clintrial()
    result <- desctable(df, by = "treatment", variables = c("sex", "stage"))
    out <- summata:::format_indented_groups(data.table::as.data.table(result))

    expect_identical(out$`p-value`[out$Variable %in% c("sex", "stage")],
                     result$`p-value`[result$Variable %in% c("sex", "stage")])
    expect_true(all(out$`p-value`[grepl("^\\s+", out$Variable)] == ""))
})


test_that("effect columns are recognized by their measure", {

    col_names <- c("Variable", "n", "OR (95% CI)", "aOR (90% CI)",
                   "Adj. Coefficient (95% CI)", "Exp(Coef) (99% CI)",
                   "Median (95% CI)", "50% Survival Time (95% CI)", "p-value")
    expect_identical(summata:::effect_ci_columns(col_names),
                     c("OR (95% CI)", "aOR (90% CI)", "Adj. Coefficient (95% CI)",
                       "Exp(Coef) (99% CI)"))
})


test_that("indented survival tables keep one p-value per outcome", {

    df <- fresh_clintrial()
    result <- survtable(df, c("Surv(os_months, os_status)", "Surv(pfs_months, pfs_status)"),
                        by = "treatment", times = 12)
    out <- summata:::format_indented_groups(data.table::as.data.table(result))

    outcome_rows <- out$Variable %in% result$Variable[result$Variable != ""]
    expect_identical(out$`p-value`[outcome_rows],
                     result$`p-value`[result$Variable != ""])
    expect_true(all(out$`p-value`[!outcome_rows] == ""))
})
