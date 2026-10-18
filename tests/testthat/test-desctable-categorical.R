## The categorical test in desctable() compares the observed levels of a
## variable across the defined groups. Observations with a missing group, and
## levels or groups without observations, must not enter the test, and an
## exact test that cannot be completed must not leave a silently blank
## p-value.

test_that("observations with a missing group do not enter the test", {

    df <- fresh_clintrial()
    df$treatment[1:40] <- NA

    with_na_row <- desctable(df, by = "treatment", variables = "stage", na_include = TRUE)
    without_na_row <- desctable(df, by = "treatment", variables = "stage")

    observed <- table(df$stage, df$treatment)
    expect_equal(raw_p(with_na_row), stats::chisq.test(observed)$p.value)
    expect_equal(raw_p(with_na_row), raw_p(without_na_row))
})


test_that("levels without observations do not change the test", {

    df <- fresh_clintrial()
    df$stage_unused <- factor(as.character(df$stage),
                              levels = c(levels(df$stage), "Unused"))

    result <- desctable(df, by = "treatment", variables = c("stage", "stage_unused"))
    p <- raw_p(result)

    expect_length(p, 2L)
    expect_equal(p[2], p[1])
    expect_equal(p[1], stats::chisq.test(table(df$stage, df$treatment))$p.value)
})


test_that("groups without observations do not change the test", {

    df <- fresh_clintrial()
    df$treatment_unused <- factor(as.character(df$treatment),
                                  levels = c(levels(df$treatment), "Placebo"))

    used <- desctable(df, by = "treatment", variables = c("sex", "stage"))
    unused <- desctable(df, by = "treatment_unused", variables = c("sex", "stage"))

    expect_equal(raw_p(unused), raw_p(used))
})


test_that("an exact test that cannot be completed falls back to simulation", {

    df <- fresh_clintrial()
    set.seed(1)
    df$sparse <- factor(sample(paste0("L", 1:12), nrow(df), replace = TRUE,
                               prob = rev(1:12)^2), levels = paste0("L", 1:12))
    tab <- table(df$sparse, df$treatment)
    tab <- tab[rowSums(tab) > 0, ]

    ## The network algorithm exhausts its default workspace on this table
    exact <- tryCatch(stats::fisher.test(tab), error = function(e) e)
    skip_if_not(inherits(exact, "error"), "fisher.test() completes this table")
    expect_match(conditionMessage(exact), "FEXACT")

    for (test_type in c("auto", "fisher")) {
        set.seed(2)
        result <- expect_warning(
            desctable(df, by = "treatment", variables = "sparse",
                      test_categorical = test_type),
            "Monte Carlo"
        )
        set.seed(2)
        expected <- stats::fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value
        expect_equal(raw_p(result), expected)
    }

    ## A requested simulation is used directly, with the requested size
    set.seed(3)
    result <- expect_no_warning(
        desctable(df, by = "treatment", variables = "sparse",
                  simulate.p.value = TRUE, B = 2000)
    )
    set.seed(3)
    expect_equal(raw_p(result),
                 stats::fisher.test(tab, simulate.p.value = TRUE, B = 2000)$p.value)
})


test_that("the automatic choice follows the expected counts", {

    df <- fresh_clintrial()

    ## All expected counts are at least 5: chi-square test
    result <- desctable(df, by = "treatment", variables = "sex")
    expect_equal(raw_p(result), stats::chisq.test(table(df$sex, df$treatment))$p.value)

    ## Some expected counts are below 5: Fisher's exact test
    small <- df[df$site %in% c("Site Alpha", "Site Beta"), ]
    small$ecog <- droplevels(small$ecog)
    tab <- table(small$ecog, small$treatment)
    expect_true(any(suppressWarnings(stats::chisq.test(tab))$expected < 5))
    result <- desctable(small, by = "treatment", variables = "ecog")
    expect_equal(raw_p(result), stats::fisher.test(tab)$p.value)
})


test_that("the total column counts observations with a missing group", {

    df <- fresh_clintrial()
    df$treatment[1:40] <- NA

    result <- desctable(df, by = "treatment", variables = c("sex", "stage"))
    raw <- attr(result, "raw_data")

    ## The total column describes the same observations as the N row
    expect_identical(result$Total[1], "850")
    expect_equal(raw$Total[raw$Variable == "sex" | raw$Group == "Male"],
                 as.vector(table(df$sex)))
    expect_equal(sum(raw$Total[raw$stat_type == "category"]),
                 sum(!is.na(df$sex)) + sum(!is.na(df$stage)))

    ## The counts do not depend on whether missing values are displayed
    shown <- desctable(df, by = "treatment", variables = c("sex", "stage"),
                       na_include = TRUE)
    shown_raw <- attr(shown, "raw_data")
    keep <- shown_raw$Group != "Unknown"
    expect_equal(shown_raw$Total[keep], raw$Total)
})


test_that("a group without observations shows counts without a percentage", {

    df <- fresh_clintrial()
    df$treatment <- factor(as.character(df$treatment),
                           levels = c(levels(df$treatment), "Placebo"))

    result <- desctable(df, by = "treatment", variables = "sex")
    expect_identical(result$Placebo, c("0", "0", "0"))

    result <- desctable(df, by = "treatment", variables = "sex",
                        stats_categorical = "percent")
    expect_identical(result$Placebo[-1], c("", ""))
})
