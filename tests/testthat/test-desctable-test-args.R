## desctable() documents '...' as arguments for the underlying statistical
## tests. A single call can summarize continuous, categorical, and
## time-to-event variables, so each argument must reach the tests that accept
## it and no others, and an argument that no test accepts must not vanish
## silently. The p-values below are checked against direct calls to the test
## functions.

test_that("var.equal reaches the t-test", {

    df <- fresh_clintrial()

    welch <- desctable(df, by = "sex", variables = "age", test_continuous = "t")
    student <- desctable(df, by = "sex", variables = "age", test_continuous = "t",
                         var.equal = TRUE)

    expect_equal(raw_p(welch), stats::t.test(age ~ sex, data = df)$p.value)
    expect_equal(raw_p(student),
                 stats::t.test(age ~ sex, data = df, var.equal = TRUE)$p.value)
    expect_false(isTRUE(all.equal(raw_p(welch), raw_p(student))))
})


test_that("arguments for different tests can be combined in one call", {

    df <- fresh_clintrial()

    set.seed(42)
    result <- expect_no_warning(
        desctable(df, by = "treatment", variables = c("age", "sex"),
                  test_continuous = "aov", var.equal = TRUE,
                  simulate.p.value = TRUE, B = 500)
    )
    p <- raw_p(result)

    ## ANOVA takes neither argument and still returns its p-value
    expect_equal(p[1], summary(stats::aov(age ~ treatment, data = df))[[1]][["Pr(>F)"]][1])

    ## The categorical test receives the simulation arguments
    set.seed(42)
    simulated <- stats::chisq.test(table(df$sex, df$treatment),
                                   simulate.p.value = TRUE, B = 500)$p.value
    expect_equal(p[2], simulated)
})


test_that("rho reaches the survival curve comparison", {

    df <- fresh_clintrial()

    result <- desctable(df, by = "treatment",
                        variables = "Surv(os_months, os_status)", rho = 1)

    expected <- survival::survdiff(survival::Surv(os_months, os_status) ~ treatment,
                                   data = df, rho = 1)$pvalue
    expect_equal(raw_p(result), expected)
})


test_that("arguments that no test accepts are reported", {

    df <- fresh_clintrial()

    expect_warning(
        desctable(df, by = "sex", variables = "age", var.eqaul = TRUE),
        "var.eqaul"
    )

    ## An argument accepted by some test is not reported even when no variable
    ## in the call uses that test
    expect_no_warning(
        desctable(df, by = "sex", variables = "age", simulate.p.value = TRUE)
    )
})


test_that("filter_test_args keeps only the arguments a test accepts", {

    args <- list(var.equal = TRUE, simulate.p.value = TRUE, B = 100, p_digits = 3)

    expect_identical(summata:::filter_test_args(args, "t"), list(var.equal = TRUE))
    expect_identical(summata:::filter_test_args(args, "fisher"),
                     list(simulate.p.value = TRUE, B = 100))
    expect_identical(summata:::filter_test_args(args, "fisher.test"),
                     summata:::filter_test_args(args, "fisher"))
    expect_identical(summata:::filter_test_args(args, "kwt"), setNames(list(), character(0)))
    expect_identical(summata:::filter_test_args(list(), "t"), list())
    expect_identical(summata:::filter_test_args(list(TRUE), "t"), list())
})
