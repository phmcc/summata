## condense_table = TRUE shows a binary variable on one row. In a descriptive
## table, that row must report the second level, which regression models
## estimate against the first, unless only the second is a negative category
## such as "No": "Yes" for a Yes/No variable in either order, labeled by the
## variable name alone, and the level named otherwise ("sex (Male)"). It must
## carry the p-value of that row in a regression table, and the labels of
## condensed continuous and survival rows must state their statistic.

condensed <- function(table) {
    summata:::condense_table_rows(data.table::as.data.table(table))
}


test_that("a condensed descriptive table reports the second level", {

    df <- fresh_clintrial()
    result <- desctable(df, by = "treatment", variables = c("hypertension", "sex"))
    out <- condensed(result)

    yes <- result[result$Group == "Yes", ]
    male <- result[result$Group == "Male", ]
    expect_identical(out$Total[out$Variable == "hypertension"], yes$Total)
    expect_identical(out$Control[out$Variable == "hypertension"], yes$Control)
    expect_identical(out$Total[out$Variable == "sex (Male)"], male$Total)

    ## The variable's p-value is kept
    expect_identical(out$`p-value`[out$Variable == "hypertension"],
                     result$`p-value`[result$Variable == "hypertension"])
})


test_that("a condensed regression table keeps the p-value of the displayed level", {

    df <- fresh_clintrial()
    result <- fit(df, "readmission_30d", c("age", "sex", "surgery"))
    out <- condensed(result)

    expect_identical(out$`p-value`[out$Variable == "sex (Male)"],
                     result$`p-value`[result$Group == "Male"])
    expect_identical(out$`p-value`[out$Variable == "surgery"],
                     result$`p-value`[result$Group == "Yes"])
})


test_that("condensed continuous and survival rows state their statistic", {

    df <- fresh_clintrial()
    result <- desctable(df, by = "treatment",
                        variables = c("age", "Surv(os_months, os_status)"),
                        stats_continuous = "mean_sd", conf_level = 0.9)

    for (out in list(condensed(result),
                     summata:::condense_quantitative_rows(data.table::as.data.table(result)))) {
        expect_true("age, mean ± SD" %in% out$Variable)
        expect_true("Surv(os_months, os_status), median (90% CI)" %in% out$Variable)
    }
})


test_that("a condensed Yes/No variable reports Yes in either order of occurrence", {

    df <- fresh_clintrial()
    df$htn <- as.character(df$hypertension)

    ## The levels of a character variable appear in order of occurrence
    for (first in c("Yes", "No")) {
        d <- df[order(df$htn != first), ]
        result <- desctable(d, by = "treatment", variables = "htn")
        expect_identical(result$Group[result$Variable == "htn"], first)

        out <- condensed(result)
        expect_identical(out$Total[out$Variable == "htn"],
                         result$Total[result$Group == "Yes"], info = first)
    }
})


test_that("survival tables are condensed as descriptive tables", {

    df <- fresh_clintrial()
    result <- survtable(df, "Surv(os_months, os_status)", by = "treatment", times = 12)

    expect_no_warning(summata:::condense_quantitative_rows(data.table::as.data.table(result)))
})
