## The print methods disclose how much of the supplied data reached the model.
## The line is unconditional: a complete sample is stated rather than left to be
## distinguished from an absent disclosure.

test_that("fit() records the analyzed sample and the reason for exclusions", {

    dt <- with_missing_outcome("readmission_30d", 100L)

    result <- fit(
        data = dt,
        outcome = "readmission_30d",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "binomial"
    )

    counts <- attr(result, "analysis_counts")

    expect_equal(counts$n_supplied, nrow(dt))
    expect_equal(counts$n_analyzed, as.numeric(stats::nobs(attr(result, "model"))))
    expect_equal(counts$n_missing_outcome, 100)
    expect_equal(counts$n_analyzed + counts$n_missing_outcome +
                 counts$n_missing_predictor, counts$n_supplied)
})


test_that("print methods report the analyzed sample when observations are excluded", {

    dt <- with_missing_outcome("readmission_30d", 100L)

    result <- fit(
        data = dt,
        outcome = "readmission_30d",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "binomial"
    )

    counts <- attr(result, "analysis_counts")
    expected <- sprintf("Observations analyzed: %d of %d",
                        as.integer(counts$n_analyzed),
                        as.integer(counts$n_supplied))

    expect_output(print(result), expected, fixed = TRUE)
})


test_that("a complete sample is reported rather than left unstated", {

    ## Reporting only on loss would leave the reader unable to distinguish a
    ## complete sample from an absent disclosure, so the line is unconditional.
    df <- fresh_clintrial()
    keep <- c("readmission_30d", "age", "sex")
    dt <- data.table::as.data.table(df[stats::complete.cases(df[, keep]), keep])

    result <- fit(
        data = dt,
        outcome = "readmission_30d",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "binomial"
    )

    counts <- attr(result, "analysis_counts")
    expect_equal(counts$n_analyzed, counts$n_supplied)

    expect_output(print(result),
                  sprintf("Observations analyzed: %d of %d (100.0%%)",
                          as.integer(counts$n_analyzed),
                          as.integer(counts$n_supplied)),
                  fixed = TRUE)
})


test_that("uniscreen reports a range when predictors differ in missingness", {

    df <- fresh_clintrial()
    df$readmission_30d[seq_len(100)] <- NA
    df$bmi[seq_len(150)] <- NA
    dt <- data.table::as.data.table(df)

    result <- uniscreen(
        data = dt,
        outcome = "readmission_30d",
        predictors = c("age", "bmi"),
        model_type = "glm",
        family = "binomial",
        parallel = FALSE
    )

    counts <- attr(result, "analysis_counts")

    ## One analyzed count per predictor, and bmi contributes the smaller one
    expect_length(counts$n_analyzed, 2L)
    expect_gt(max(counts$n_analyzed), min(counts$n_analyzed))

    expect_output(print(result), "Observations analyzed: [0-9]+-[0-9]+ of ")
})


test_that("a blank line separates the header from the descriptors", {

    dt <- data.table::as.data.table(fresh_clintrial())

    screened <- uniscreen(
        data = dt,
        outcome = "readmission_30d",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "binomial",
        parallel = FALSE
    )

    out <- utils::capture.output(print(screened))

    expect_equal(out[2], "Univariable Screening Results")
    expect_equal(out[3], "")
    expect_match(out[4], "^Outcome: ")
})


test_that("fit() print reports observations and events on the same footing", {

    dt <- with_missing_outcome("readmission_30d", 100L)

    result <- fit(
        data = dt,
        outcome = "readmission_30d",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "binomial"
    )

    counts <- attr(result, "analysis_counts")
    raw <- attr(result, "raw_data")

    ## The events denominator counts events in the supplied data using the same
    ## coding as the model, so it is never smaller than the events analyzed
    expect_true(!is.null(counts$events_supplied))
    expect_gte(counts$events_supplied, counts$events_analyzed)
    expect_equal(counts$events_analyzed, as.numeric(raw$events[1]))

    out <- utils::capture.output(print(result))

    expect_true(any(grepl(sprintf("Observations analyzed: %d of %d",
                                  as.integer(counts$n_analyzed),
                                  as.integer(counts$n_supplied)),
                          out, fixed = TRUE)))
    expect_true(any(grepl(sprintf("Events analyzed: %d of %d",
                                  as.integer(counts$events_analyzed),
                                  as.integer(counts$events_supplied)),
                          out, fixed = TRUE)))

    ## The bare "n = " line is gone, its figure now carried by the first line
    expect_false(any(grepl("^n = ", out)))
})


test_that("linear models report no events line", {

    dt <- data.table::as.data.table(fresh_clintrial())

    result <- fit(
        data = dt,
        outcome = "los_days",
        predictors = c("age", "sex"),
        model_type = "lm"
    )

    out <- utils::capture.output(print(result))

    expect_true(any(grepl("Observations analyzed:", out, fixed = TRUE)))
    expect_false(any(grepl("Events analyzed:", out, fixed = TRUE)))
})


test_that("screening workflows report events alongside observations", {

    df <- fresh_clintrial()
    df$readmission_30d[seq_len(100)] <- NA
    df$bmi[seq_len(150)] <- NA
    dt <- data.table::as.data.table(df)

    screened <- uniscreen(
        data = dt,
        outcome = "readmission_30d",
        predictors = c("age", "bmi"),
        model_type = "glm",
        family = "binomial",
        parallel = FALSE
    )

    counts <- attr(screened, "analysis_counts")

    ## One value per screened predictor, with a shared denominator
    expect_length(counts$events_analyzed, 2L)
    expect_length(counts$events_supplied, 1L)
    expect_true(all(counts$events_analyzed <= counts$events_supplied))

    out <- utils::capture.output(print(screened))
    expect_true(any(grepl("Observations analyzed:", out, fixed = TRUE)))
    expect_true(any(grepl("Events analyzed:", out, fixed = TRUE)))
})


test_that("multifit reports observations only", {

    df <- fresh_clintrial()
    df$bmi[seq_len(100)] <- NA
    dt <- data.table::as.data.table(df)

    result <- multifit(
        data = dt,
        outcomes = c("readmission_30d", "os_status"),
        predictor = "sex",
        covariates = c("age", "bmi"),
        model_type = "glm",
        family = "binomial"
    )

    counts <- attr(result, "analysis_counts")

    ## Each outcome has its own event total, so no single figure could describe
    ## them all and events are deliberately omitted
    expect_null(counts$events_analyzed)

    out <- utils::capture.output(print(result))
    expect_true(any(grepl("Observations analyzed:", out, fixed = TRUE)))
    expect_false(any(grepl("Events analyzed:", out, fixed = TRUE)))
})


test_that("the event resolver accepts a family in any of its forms", {

    outcome_vars <- "readmission_30d"

    ## A name, a resolved family object, and the generator function itself
    expect_equal(summata:::get_event_variable_for_counts(
        outcome_vars, "glm", "binomial"), "readmission_30d")
    expect_equal(summata:::get_event_variable_for_counts(
        outcome_vars, "glm", stats::binomial()), "readmission_30d")
    expect_equal(summata:::get_event_variable_for_counts(
        outcome_vars, "glm", stats::binomial), "readmission_30d")

    ## Families that carry no event count, including one already resolved to an
    ## object by the caller
    expect_null(summata:::get_event_variable_for_counts(
        outcome_vars, "glm", "gaussian"))
    expect_null(summata:::get_event_variable_for_counts(
        outcome_vars, "glm", stats::Gamma(link = "log")))
    expect_null(summata:::get_event_variable_for_counts(
        outcome_vars, "lm", NULL))
})


test_that("Gamma and inverse Gaussian screens report no events line", {

    dt <- data.table::as.data.table(fresh_clintrial())

    result <- suppressWarnings(uniscreen(
        data = dt,
        outcome = "los_days",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "Gamma",
        parallel = FALSE
    ))

    counts <- attr(result, "analysis_counts")
    expect_null(counts$events_analyzed)

    out <- utils::capture.output(print(result))
    expect_true(any(grepl("Observations analyzed:", out, fixed = TRUE)))
    expect_false(any(grepl("Events analyzed:", out, fixed = TRUE)))
})


test_that("printed counts follow the number_format setting", {

    ## A dataset large enough for a thousands separator to appear
    df <- fresh_clintrial()
    df <- do.call(rbind, replicate(3, df, simplify = FALSE))
    df$readmission_30d[seq_len(200)] <- NA
    dt <- data.table::as.data.table(df)

    us <- fit(
        data = dt,
        outcome = "readmission_30d",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "binomial",
        number_format = "us"
    )

    eu <- fit(
        data = dt,
        outcome = "readmission_30d",
        predictors = c("age", "sex"),
        model_type = "glm",
        family = "binomial",
        number_format = "eu"
    )

    ## format() warns when big.mark matches the decimal mark it would default
    ## to, so both marks must be supplied together
    expect_no_warning(us_out <- utils::capture.output(print(us)))
    expect_no_warning(eu_out <- utils::capture.output(print(eu)))

    us_line <- us_out[grepl("Observations analyzed", us_out, fixed = TRUE)]
    eu_line <- eu_out[grepl("Observations analyzed", eu_out, fixed = TRUE)]

    ## Counts carry the thousands separator and percentages the decimal mark
    expect_match(us_line, "[0-9],[0-9]{3}")
    expect_match(us_line, "\\.[0-9]%")
    expect_match(eu_line, "[0-9]\\.[0-9]{3}")
    expect_match(eu_line, ",[0-9]%")
})
