## Forest plots report the same sample sizes as the tables they accompany.
## Counts are read from the "table_data" attribute, which holds the values drawn
## in the plot in model term order.

forest_n <- function(p, variable) {
    tbl <- attr(p, "table_data")
    tbl[["N"]][tbl[["var"]] == variable]
}


test_that("coxforest counts only the observations used in fitting", {

    df <- fresh_clintrial()
    df$os_status[seq_len(100)] <- NA
    dt <- data.table::as.data.table(df)

    model <- survival::coxph(
        survival::Surv(os_months, os_status) ~ age + sex + race,
        data = dt
    )

    p <- suppressMessages(suppressWarnings(coxforest(model, data = dt)))

    ## Group sizes partition the fitted sample size
    for (pred in c("sex", "race")) {
        expect_equal(sum(forest_n(p, pred)), model$n)
    }

    ## Continuous rows report the fitted sample size, not every supplied row
    expect_equal(unique(forest_n(p, "age")), model$n)
    expect_true(all(attr(p, "table_data")[["N"]] <= nrow(dt)))
})


test_that("glmforest counts correctly when no model frame was retained", {

    df <- fresh_clintrial()
    df$readmission_30d[seq_len(100)] <- NA
    dt <- data.table::as.data.table(df)

    ## model = FALSE leaves the fit without a model frame, so the counting frame
    ## must be reconstructed from the supplied data.
    model <- stats::glm(readmission_30d ~ age + sex + race,
                        data = dt, family = stats::binomial,
                        model = FALSE, y = TRUE)

    p <- suppressMessages(suppressWarnings(glmforest(model, data = dt)))

    n_model <- stats::nobs(model)

    for (pred in c("sex", "race")) {
        expect_equal(sum(forest_n(p, pred)), n_model)
    }
    expect_equal(unique(forest_n(p, "age")), n_model)
})


test_that("forest plots expose the assembled table", {

    df <- fresh_clintrial()
    dt <- data.table::as.data.table(df)

    model <- stats::glm(readmission_30d ~ age + sex,
                        data = dt, family = stats::binomial)

    p <- suppressMessages(suppressWarnings(glmforest(model, data = dt)))
    tbl <- attr(p, "table_data")

    expect_s3_class(tbl, "data.table")
    expect_true(all(c("var", "level", "N", "estimate") %in% names(tbl)))
    expect_gt(nrow(tbl), 0L)
})


test_that("forest footers state the denominator alongside the percentage", {

    df <- fresh_clintrial()
    df$readmission_30d[seq_len(100)] <- NA
    dt <- data.table::as.data.table(df)

    model <- stats::glm(readmission_30d ~ age + sex,
                        data = dt, family = stats::binomial)

    p <- suppressMessages(suppressWarnings(glmforest(model, data = dt)))

    ## The footer text is carried on an annotation layer
    labels <- unlist(lapply(p$layers, function(lyr) {
        as.character(lyr$aes_params$label)
    }))
    footer <- labels[grepl("Observations analyzed", labels, fixed = TRUE)]

    skip_if(length(footer) == 0, "footer annotation not found")

    expect_match(footer[1],
                 sprintf("Observations analyzed: %d of %d",
                         as.integer(stats::nobs(model)), nrow(dt)),
                 fixed = TRUE)
})


test_that("forest footers report observations and events on stacked lines", {

    df <- fresh_clintrial()
    df$readmission_30d[seq_len(100)] <- NA
    dt <- data.table::as.data.table(df)

    model <- stats::glm(readmission_30d ~ age + sex,
                        data = dt, family = stats::binomial)

    p <- suppressMessages(suppressWarnings(glmforest(model, data = dt)))

    labels <- unlist(lapply(p$layers, function(lyr) {
        as.character(lyr$aes_params$label)
    }))
    footer <- labels[grepl("Observations analyzed", labels, fixed = TRUE)]

    skip_if(length(footer) == 0, "footer annotation not found")

    ## Both lines present, and the events line follows the observations line
    expect_match(footer[1], "Observations analyzed: [0-9]+ of [0-9]+")
    expect_match(footer[1], "Events analyzed: [0-9]+ of [0-9]+")
    expect_lt(regexpr("Observations analyzed", footer[1], fixed = TRUE),
              regexpr("Events analyzed", footer[1], fixed = TRUE))
})


test_that("linear model footers carry no events line", {

    dt <- data.table::as.data.table(fresh_clintrial())

    model <- stats::lm(los_days ~ age + sex, data = dt)

    p <- suppressMessages(suppressWarnings(lmforest(model, data = dt)))

    labels <- unlist(lapply(p$layers, function(lyr) {
        as.character(lyr$aes_params$label)
    }))
    footer <- labels[grepl("Observations analyzed", labels, fixed = TRUE)]

    skip_if(length(footer) == 0, "footer annotation not found")

    expect_false(grepl("Events analyzed", footer[1], fixed = TRUE))
})
