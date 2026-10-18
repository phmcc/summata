## compfit() ranks candidate models by a composite score. Each model must be
## assessed on its own merits: a rank-deficient candidate is still assessed,
## every weighted component contributes to the score, and the score does not
## depend on how numbers are formatted for display.

test_that("compfit assesses a model with an aliased coefficient", {

    df <- fresh_clintrial()
    df <- df[stats::complete.cases(df[, c("readmission_30d", "age", "sex")]), ]
    df$age_doubled <- 2 * df$age

    result <- suppressMessages(compfit(df, "readmission_30d",
                                       list(a = c("age", "sex"),
                                            b = c("age", "age_doubled", "sex"))))

    expect_setequal(result$Model, c("a", "b"))
    expect_false(any(is.na(result$CMS)))
    expect_true(all(result$Converged == "Yes"))
})


test_that("every weighted component contributes to the logistic score", {

    df <- fresh_clintrial()
    df <- df[stats::complete.cases(df[, c("readmission_30d", "age", "sex", "stage")]), ]

    result <- suppressMessages(compfit(df, "readmission_30d",
                                       list(a = "age", b = c("age", "sex", "stage"))))
    scores <- attr(result, "detailed_scores")
    weights <- attr(result, "weights")

    expect_equal(sum(unlist(weights)), 1)
    expected <- scores$conv_score * weights$convergence +
        scores$aic_score * weights$aic +
        scores$concordance_score * weights$concordance +
        scores$pseudo_r2_score * weights$pseudo_r2 +
        scores$brier_score * weights$brier
    ## The detailed scores keep the input order; the table is sorted by score
    expect_equal(sort(result$CMS), sort(round(expected, 1)))
})


test_that("the score does not depend on the number format", {

    df <- fresh_clintrial()
    df <- df[stats::complete.cases(df[, c("os_months", "os_status", "age", "sex", "stage")]), ]
    models <- list(a = "age", b = c("age", "sex", "stage"))

    us <- suppressMessages(compfit(df, "Surv(os_months, os_status)", models,
                                   model_type = "coxph"))
    eu <- suppressMessages(compfit(df, "Surv(os_months, os_status)", models,
                                   model_type = "coxph", number_format = "eu"))
    fine <- suppressMessages(compfit(df, "Surv(os_months, os_status)", models,
                                     model_type = "coxph", p_digits = 6))

    expect_equal(eu$CMS, us$CMS)
    expect_equal(fine$CMS, us$CMS)
})


test_that("linear mixed models are compared on maximum likelihood fits", {

    skip_if_not_installed("lme4")
    df <- fresh_clintrial()
    df <- df[stats::complete.cases(df[, c("bmi", "age", "sex", "smoking", "site")]), ]

    result <- suppressMessages(compfit(df, "bmi",
                                       list(a = c("age", "(1|site)"),
                                            b = c("age", "sex", "smoking", "(1|site)")),
                                       model_type = "lmer"))
    result <- result[order(result$Model)]

    full <- suppressMessages(lme4::lmer(bmi ~ age + sex + smoking + (1 | site), data = df,
                                        REML = FALSE))
    null <- suppressMessages(lme4::lmer(bmi ~ 1 + (1 | site), data = df, REML = FALSE))
    lr <- suppressMessages(stats::anova(null, full))

    expect_equal(result$AIC[2], round(stats::AIC(full), 1))
    expect_equal(result$BIC[2], round(stats::BIC(full), 1))
    expect_identical(result$`Global p`[2],
                     summata:::format_pvalues_fit(lr$`Pr(>Chisq)`[2], 3))
})
