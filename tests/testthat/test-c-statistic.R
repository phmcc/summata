## The c-statistic of logistic models is computed internally rather than with
## the suggested pROC package, so that concordance, and the composite score of
## compfit() that weights it, do not depend on which packages are installed.
## It must equal the area under the empirical ROC curve, oriented as
## pROC::roc(direction = "auto") orients it. The reference below counts
## concordant case-control pairs directly.

pairwise_auc <- function(y, p) {
    cases <- p[y == 1]
    controls <- p[y == 0]
    concordant <- outer(cases, controls, ">") + 0.5 * outer(cases, controls, "==")
    auc <- mean(concordant)
    if (stats::median(controls) > stats::median(cases)) 1 - auc else auc
}


test_that("c_statistic equals the pairwise concordance", {

    set.seed(20261009)

    for (i in 1:200) {
        n <- sample(c(10, 60, 500), 1)
        y <- stats::rbinom(n, 1, 0.4)
        if (length(unique(y)) < 2) next
        p <- switch(i %% 3 + 1,
                    stats::plogis(stats::rnorm(n) + y * stats::runif(1, -1, 2)),
                    round(stats::plogis(stats::rnorm(n) + y), 1),
                    stats::plogis(stats::rnorm(n)))
        expect_equal(summata:::c_statistic(y, p), pairwise_auc(y, p), tolerance = 1e-12)
    }
})


test_that("c_statistic handles ties, missing values, and a single outcome", {

    expect_equal(summata:::c_statistic(c(0, 0, 1, 1), c(0.1, 0.5, 0.5, 0.9)), 0.875)
    expect_equal(summata:::c_statistic(c(0, 1, 0, 1), rep(0.5, 4)), 0.5)
    expect_equal(summata:::c_statistic(c(0, 1, NA, 1), c(0.2, 0.8, 0.5, NA)), 1)
    expect_true(is.na(summata:::c_statistic(c(1, 1, 1), c(0.2, 0.4, 0.6))))
})


test_that("logistic model summaries report the c-statistic", {

    df <- fresh_clintrial()
    model <- stats::glm(readmission_30d ~ age + sex, data = df, family = stats::binomial)

    raw <- m2dt(df, model)
    expect_equal(unique(raw$c_statistic),
                 summata:::c_statistic(model$y, stats::fitted(model)))

    result <- suppressMessages(compfit(df[stats::complete.cases(df[, c("readmission_30d", "age", "sex")]), ],
                                       "readmission_30d", list(a = "age", b = c("age", "sex"))))
    expect_false(any(is.na(result$Concordance)))
})
