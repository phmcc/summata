## survtable() reports each group from its own stratum of the stratified
## survfit() fit, including groups that survfit() leaves out or does not name,
## and lists the groups in the order of the data.


### * Fixtures

## A grouping variable whose groups first appear in neither sorted nor reverse
## sorted order
groups_data <- function() {
    set.seed(20261011)
    n <- 90
    data.frame(
        time = round(stats::rexp(n, 0.08), 1),
        status = stats::rbinom(n, 1, 0.75),
        arm = sample(rep(c("Drug B", "Zcontrol", "Drug A"), each = 30)),
        stringsAsFactors = FALSE
    )
}


## Estimates of one group fitted on its own
own_fit <- function(df, group) {
    survival::survfit(survival::Surv(time, status) ~ 1,
                      data = df[df$arm == group, ])
}


### * Estimates of each group

test_that("each group's estimates come from its own stratum when a level is empty", {

    df <- groups_data()
    df$arm <- factor(df$arm, levels = c("Drug A", "Placebo", "Drug B", "Zcontrol"))

    result <- survtable(df, "Surv(time, status)", by = "arm", times = 12,
                        probs = c(0.25, 0.5), total = FALSE, test = FALSE)
    raw <- attr(result, "raw_data")

    for (g in c("Drug A", "Drug B", "Zcontrol")) {
        fit <- own_fit(df, g)
        row <- raw$Variable == g
        expect_equal(raw[["12_estimate"]][row],
                     summary(fit, times = 12, extend = TRUE)$surv, info = g)
        expect_equal(raw[["Median (95% CI)_estimate"]][row],
                     unname(stats::quantile(fit, 0.5)$quantile), info = g)
        expect_equal(raw[["75% Survival Time (95% CI)_estimate"]][row],
                     unname(stats::quantile(fit, 0.25)$quantile), info = g)
    }

    ## The empty level is listed, as desctable() lists it, without estimates
    expect_identical(result$Variable, c("Drug A", "Placebo", "Drug B", "Zcontrol"))
    empty <- result$Variable == "Placebo"
    expect_identical(result[["12"]][empty], "-")
    expect_identical(result[["Median (95% CI)"]][empty], "-")
    expect_identical(result[["75% Survival Time (95% CI)"]][empty], "-")
    expect_true(is.na(raw[["Median (95% CI)_estimate"]][empty]))
})


test_that("a single group with observations is reported from its own fit", {

    df <- groups_data()
    df$arm <- factor(df$arm)
    single <- df[df$arm == "Drug B", ]
    fit <- own_fit(single, "Drug B")

    ## survfit() returns no strata when one group has observations; the
    ## factor keeps its other levels, which have no estimates
    result <- survtable(single, "Surv(time, status)", by = "arm", times = 12,
                        probs = c(0.25, 0.5), test = FALSE)
    raw <- attr(result, "raw_data")
    row <- raw$Variable == "Drug B"

    expect_equal(raw[["12_estimate"]][row],
                 summary(fit, times = 12, extend = TRUE)$surv)
    expect_equal(raw[["Median (95% CI)_estimate"]][row],
                 unname(stats::quantile(fit, 0.5)$quantile))
    expect_equal(raw[["75% Survival Time (95% CI)_estimate"]][row],
                 unname(stats::quantile(fit, 0.25)$quantile))

    others <- result$Variable %in% c("Drug A", "Zcontrol")
    expect_true(all(result[["12"]][others] == "-"))
    expect_true(all(result[["Median (95% CI)"]][others] == "-"))

    ## The same holds for a character variable with a single value
    single$arm <- as.character(single$arm)
    result_chr <- survtable(single, "Surv(time, status)", by = "arm", times = 12,
                            probs = c(0.25, 0.5), total = FALSE, test = FALSE)
    raw_chr <- attr(result_chr, "raw_data")

    expect_identical(result_chr$Variable, "Drug B")
    expect_equal(raw_chr[["12_estimate"]],
                 summary(fit, times = 12, extend = TRUE)$surv)
    expect_equal(raw_chr[["Median (95% CI)_estimate"]],
                 unname(stats::quantile(fit, 0.5)$quantile))
})


### * Order of the groups

## As in desctable(): the levels of a factor as given, and the values of other
## grouping variables in order of first appearance

test_that("character groups appear in order of first appearance", {

    df <- groups_data()
    appearance <- unique(df$arm)

    ## The fixture is ordered neither alphabetically nor in reverse
    expect_false(identical(appearance, sort(appearance)))
    expect_false(identical(appearance, rev(sort(appearance))))

    first <- survtable(df, "Surv(time, status)", by = "arm", times = 12)
    expect_identical(first$Variable, c("Total", appearance))

    last <- survtable(df, "Surv(time, status)", by = "arm", probs = 0.5,
                      total = "last")
    expect_identical(last$Variable, c(appearance, "Total"))

    ## Times and quantiles together take the merged path
    both <- survtable(df, "Surv(time, status)", by = "arm", times = 12,
                      probs = c(0.25, 0.5), total = FALSE)
    expect_identical(both$Variable, appearance)
})


test_that("each reordered group keeps its own estimates", {

    df <- groups_data()
    result <- survtable(df, "Surv(time, status)", by = "arm", times = c(6, 12),
                        probs = c(0.25, 0.5), total = FALSE, test = FALSE)
    raw <- attr(result, "raw_data")

    for (g in unique(df$arm)) {
        fit <- own_fit(df, g)
        row <- raw$Variable == g
        expect_equal(c(raw[["6_estimate"]][row], raw[["12_estimate"]][row]),
                     summary(fit, times = c(6, 12), extend = TRUE)$surv,
                     info = g)
        expect_equal(raw[["Median (95% CI)_estimate"]][row],
                     unname(stats::quantile(fit, 0.5)$quantile), info = g)
        expect_equal(raw[["75% Survival Time (95% CI)_estimate"]][row],
                     unname(stats::quantile(fit, 0.25)$quantile), info = g)
    }
})


test_that("logical and numeric groups appear in order of first appearance", {

    df <- groups_data()
    ## TRUE first and doses unsorted, so that neither order is the sorted one
    df$flag <- df$arm == df$arm[1]
    df$dose <- unname(stats::setNames(c(20, 0, 10), unique(df$arm))[df$arm])
    expect_true(df$flag[1])
    expect_false(identical(unique(df$dose), sort(unique(df$dose))))

    for (by in c("flag", "dose")) {
        result <- survtable(df, "Surv(time, status)", by = by, times = 12,
                            total = FALSE)
        expect_identical(result$Variable, as.character(unique(df[[by]])),
                         info = by)
    }
})


test_that("factor groups follow their levels", {

    df <- groups_data()
    df$arm <- factor(df$arm, levels = c("Zcontrol", "Drug A", "Drug B"))

    result <- survtable(df, "Surv(time, status)", by = "arm", times = 12,
                        total = FALSE)

    expect_identical(result$Variable, c("Zcontrol", "Drug A", "Drug B"))
})


test_that("the stored survfit() strata follow the order of the table", {

    df <- groups_data()
    result <- survtable(df, "Surv(time, status)", by = "arm", times = 12)
    fit <- attr(result, "survfit_objects")$stratified

    expect_identical(names(fit$strata), paste0("arm=", unique(df$arm)))
})


test_that("labels and missing groups follow the order of appearance", {

    df <- groups_data()
    df$arm[c(1, 5)] <- NA
    appearance <- unique(stats::na.omit(df$arm))
    labels <- c("Drug A" = "Active A", "Drug B" = "Active B",
                "Zcontrol" = "Control")

    result <- survtable(df, "Surv(time, status)", by = "arm", times = 12,
                        total = FALSE, labels = labels)

    expect_identical(result$Variable, unname(labels[appearance]))
})


test_that("the order of the groups leaves the test unchanged", {

    df <- groups_data()
    result <- survtable(df, "Surv(time, status)", by = "arm", times = 12)

    expected <- survival::survdiff(survival::Surv(time, status) ~ arm, data = df)
    expect_equal(attr(result, "test_result")$p_value, expected$pvalue)
})


test_that("the grouping variable of the caller's data is not modified", {

    dt <- data.table::as.data.table(groups_data())

    for (na_rm in c(TRUE, FALSE)) {
        result <- survtable(dt, "Surv(time, status)", by = "arm", times = 12,
                            na_rm = na_rm)
        expect_type(dt$arm, "character")
        expect_identical(result$Variable, c("Total", unique(dt$arm)),
                         info = paste("na_rm =", na_rm))
    }
})


test_that("each outcome of a stacked table lists its groups in order of appearance", {

    df <- groups_data()
    df$time2 <- df$time * 1.5

    result <- survtable(df, c("Surv(time, status)", "Surv(time2, status)"),
                        by = "arm", times = 12, total = FALSE)

    expect_identical(result$Group, rep(unique(df$arm), 2))
})
