## uniscreen() and multifit() fit one model per predictor or outcome, either
## sequentially or in parallel. A model that cannot be fitted must be reported
## by name, in the same way on both paths; in a forked worker the failure
## previously surfaced only as "subscript out of bounds" when the results were
## combined.

test_that("uniscreen names a predictor that cannot be fitted", {

    df <- fresh_clintrial()
    df$one_level <- factor("only")

    for (use_parallel in c(FALSE, TRUE)) {
        expect_error(
            uniscreen(df, "os_status", c("age", "one_level", "sex"),
                      parallel = use_parallel, n_cores = 2),
            "Model fitting failed for 1 of 3 predictors:\n  one_level: contrasts"
        )
    }
})


test_that("stop_on_failed_fits reports every failure with its reason", {

    ok <- list(raw = 1)
    captured <- summata:::capture_fit_error(function(x) stop("cannot fit ", x))("b")
    forked <- structure("Error : worker failed\n", class = "try-error",
                        condition = simpleError("worker failed"))

    expect_identical(summata:::stop_on_failed_fits(list(ok, ok), c("a", "b"), "outcome"),
                     list(ok, ok))
    expect_error(
        summata:::stop_on_failed_fits(list(ok, captured, forked, NULL), c("a", "b", "c", "d"),
                                      "outcome"),
        paste0("Model fitting failed for 3 of 4 outcomes:\n",
               "  b: cannot fit b\n",
               "  c: worker failed\n",
               "  d: the worker process returned no result"),
        fixed = TRUE
    )
})


test_that("resolve_core_count respects the check limit and a requested number", {

    withr::local_envvar(c("_R_CHECK_LIMIT_CORES_" = "true"))
    expect_identical(summata:::resolve_core_count(8), 2L)
    expect_identical(summata:::resolve_core_count(1), 1L)
    expect_lte(summata:::resolve_core_count(NULL), 2L)

    withr::local_envvar(c("_R_CHECK_LIMIT_CORES_" = ""))
    expect_identical(summata:::resolve_core_count(3), 3L)
    expect_gte(summata:::resolve_core_count(NULL), 1L)
})


test_that("resolve_core_count uses one core when the core count is unknown", {

    local_mocked_bindings(detectCores = function(...) NA_integer_, .package = "parallel")
    withr::local_envvar(c("_R_CHECK_LIMIT_CORES_" = ""))
    expect_identical(summata:::resolve_core_count(NULL), 1L)
})
