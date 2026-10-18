## survtable() compares survival curves with the log-rank test or one of three
## weighted log-rank tests, each of which must compute the test it names:
## Gehan-Breslow ("wilcoxon"), Tarone-Ware ("tarone") and Peto-Peto
## ("petopeto"). "wilcoxon" and "petopeto" previously both computed
## survival::survdiff(rho = 1), and "tarone" survdiff(rho = 0.5), none of
## which is the named test. The published values below come from Klein and
## Moeschberger (2003), Table 7.3, and the SAS PROC LIFETEST documentation.

## Kidney dialysis data of Klein and Moeschberger (2003), Section 1.4: months to
## exit-site infection by catheter placement (1 = surgical, 2 = percutaneous)
kidney_data <- function() {
    data.frame(
        time = c(
            ## Surgical placement: 15 infections, then 28 censored
            1.5, 3.5, 4.5, 4.5, 5.5, 8.5, 8.5, 9.5, 10.5, 11.5, 15.5, 16.5,
            18.5, 23.5, 26.5,
            2.5, 2.5, 3.5, 3.5, 3.5, 4.5, 5.5, 6.5, 6.5, 7.5, 7.5, 7.5, 7.5,
            8.5, 9.5, 10.5, 11.5, 12.5, 12.5, 13.5, 14.5, 14.5, 21.5, 21.5,
            22.5, 22.5, 25.5, 27.5,
            ## Percutaneous placement: 11 infections, then 65 censored
            0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 2.5, 2.5, 3.5, 6.5, 15.5,
            0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.5, 1.5, 1.5,
            1.5, 2.5, 2.5, 2.5, 2.5, 2.5, 3.5, 3.5, 3.5, 3.5, 3.5, 4.5, 4.5,
            4.5, 5.5, 5.5, 5.5, 5.5, 5.5, 6.5, 7.5, 7.5, 7.5, 8.5, 8.5, 8.5,
            9.5, 9.5, 10.5, 10.5, 10.5, 11.5, 11.5, 12.5, 12.5, 12.5, 12.5,
            14.5, 14.5, 16.5, 16.5, 18.5, 19.5, 19.5, 19.5, 20.5, 22.5, 24.5,
            25.5, 26.5, 26.5, 28.5),
        status = rep(c(1, 0, 1, 0), c(15, 28, 11, 65)),
        type = rep(1:2, c(43, 76))
    )
}

## Survival times in days of the two treatment groups of the SAS PROC LIFETEST
## "Getting Started" example, coded here 1 and 0; the codes do not enter the
## tests
exposure_data <- function() {
    data.frame(
        days = c(179, 256, 262, 256, 255, 224, 225, 287, 319, 264,
                 378, 355, 319, 256, 171, 325, 325, 217, 255, 256,
                 237, 156, 270, 257, 242, 157, 249, 180, 226, 268,
                 291, 323, 253, 206, 206, 237, 211, 229, 234, 209),
        status = c(1, 1, 1, 1, 1, 0, 1, 1, 1, 1,
                   0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   0, 1, 1, 1, 1, 1, 1, 1, 1, 0,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        treatment = rep(c(1, 0), each = 20)
    )
}


test_that("the survival tests reproduce Klein and Moeschberger, Table 7.3", {

    d <- kidney_data()
    surv <- survival::Surv(d$time, d$status)

    ## Z_1(tau), its variance, chi-square and p-value as published
    published <- list(logrank = c(3.96, 6.21, 2.53, 0.112),
                      wilcoxon = c(-9, 38862, 0.002, 0.964),
                      tarone = c(13.20, 432.83, 0.40, 0.526),
                      petopeto = c(2.47, 4.36, 1.40, 0.237))
    digits <- list(logrank = c(2, 2, 2, 3), wilcoxon = c(0, 0, 3, 3),
                   tarone = c(2, 2, 2, 3), petopeto = c(2, 2, 2, 3))

    for (type in names(published)) {
        result <- summata:::perform_survival_test(surv, d$type, type)
        test <- result$test_object
        computed <- c(test$obs[1] - test$exp[1], test$var[1, 1],
                      result$statistic, result$p_value)
        expect_equal(round(computed, digits[[type]]), published[[type]], info = type)
        expect_identical(result$df, 1L, info = type)
    }
})


test_that("the Gehan-Breslow test reproduces the SAS Wilcoxon test", {

    d <- exposure_data()
    surv <- survival::Surv(d$days, d$status)

    ## Test of Equality over Strata: Log-Rank 5.6485 (p = 0.0175),
    ## Wilcoxon 5.0312 (p = 0.0249)
    logrank <- summata:::perform_survival_test(surv, d$treatment, "logrank")
    gehan <- summata:::perform_survival_test(surv, d$treatment, "wilcoxon")
    expect_equal(round(logrank$statistic, 4), 5.6485)
    expect_equal(round(logrank$p_value, 4), 0.0175)
    expect_equal(round(gehan$statistic, 4), 5.0312)
    expect_equal(round(gehan$p_value, 4), 0.0249)
})


test_that("the weighted statistic equals survdiff for the G-rho weights", {

    df <- fresh_clintrial()
    surv <- survival::Surv(df$os_months, df$os_status)
    km_before <- function(rho) {
        function(n_risk, n_event) {
            c(1, cumprod(1 - n_event / n_risk))[seq_along(n_risk)]^rho
        }
    }

    for (group in c("sex", "treatment", "stage")) {
        for (rho in c(0, 0.5, 1)) {
            weighted <- summata:::weighted_logrank_test(surv, df[[group]], km_before(rho))
            reference <- survival::survdiff(surv ~ df[[group]], rho = rho)
            expect_equal(weighted$chisq, reference$chisq, tolerance = 1e-10,
                         info = paste(group, rho))
            expect_equal(weighted$pvalue, reference$pvalue, tolerance = 1e-10,
                         info = paste(group, rho))
        }
    }
})


test_that("survtable reports the test it names", {

    df <- fresh_clintrial()
    surv <- survival::Surv(df$os_months, df$os_status)
    p_values <- numeric()

    for (type in c("logrank", "wilcoxon", "tarone", "petopeto")) {
        result <- survtable(df, "Surv(os_months, os_status)", by = "treatment",
                            times = 12, test_type = type)
        p_values[type] <- attr(result, "test_result")$p_value
        expect_equal(p_values[[type]],
                     summata:::perform_survival_test(surv, df$treatment, type)$p_value,
                     info = type)
    }

    ## Four different tests, where "wilcoxon" and "petopeto" were the same
    expect_length(unique(p_values), 4L)

    printed <- utils::capture.output(print(
        survtable(df, "Surv(os_months, os_status)", by = "treatment", times = 12,
                  test_type = "wilcoxon")))
    expect_true(any(grepl("Test: Gehan-Breslow", printed, fixed = TRUE)))
})


test_that("groups without subjects at risk at any event time are left out", {

    ## Group c is censored before the first event, so that survdiff() tests
    ## two groups on one degree of freedom; the log-rank result reported two
    d <- data.frame(time = c(1, 2, 3, 4, 0.1, 0.2), status = c(1, 1, 0, 1, 0, 0),
                    group = c("a", "a", "b", "b", "c", "c"))
    surv <- survival::Surv(d$time, d$status)
    reference <- survival::survdiff(surv ~ d$group)

    for (type in c("logrank", "wilcoxon", "tarone", "petopeto")) {
        result <- summata:::perform_survival_test(surv, d$group, type)
        expect_identical(result$df, 1L, info = type)
        expect_false(is.na(result$p_value), info = type)
    }

    logrank <- summata:::perform_survival_test(surv, d$group, "logrank")
    expect_equal(logrank$p_value,
                 stats::pchisq(reference$chisq, logrank$df, lower.tail = FALSE))
})


test_that("no test is reported when fewer than two groups expect events", {

    ## Without events no group expects any, and the statistic, its degrees of
    ## freedom and the p-value are all missing rather than undefined
    d <- data.frame(time = c(1, 2, 3, 4), status = c(0, 0, 0, 0),
                    group = c("a", "a", "b", "b"))
    surv <- survival::Surv(d$time, d$status)

    for (type in c("logrank", "wilcoxon", "tarone", "petopeto")) {
        result <- suppressWarnings(
            summata:::perform_survival_test(surv, d$group, type))
        expect_true(is.na(result$statistic), info = type)
        expect_identical(result$df, NA_integer_, info = type)
        expect_true(is.na(result$p_value), info = type)
    }
})
