## The events of a Cox model are those of the event argument of its Surv()
## response, coded as Surv() codes them: the second argument of
## Surv(time, event), the third of Surv(start, stop, event), or the argument
## named event, whether the status is coded 0 and 1 or 1 and 2 or is given as
## an expression. clogit() writes its response as Surv(rep(1, N), case), whose
## event is case. Tables, forest plots, and the counts of analyzed events all
## count them so.


### * Fixtures

## survival::lung codes its status 1 (censored) and 2 (dead)
lung_data <- function() {
    dt <- data.table::as.data.table(survival::lung)
    dt[, sex_f := factor(sex, labels = c("Male", "Female"))]
    dt[, dead := status == 2]
    dt[]
}

## A matched case-control design with three subjects per set
case_control_data <- function() {
    set.seed(20261016)
    n_sets <- 60
    dt <- data.table::data.table(set = rep(seq_len(n_sets), each = 3),
                                 case = rep(c(1, 0, 0), n_sets))
    dt[, exposure := factor(ifelse(stats::runif(.N) < 0.3 + 0.3 * case,
                                   "Yes", "No"), levels = c("No", "Yes"))]
    dt[, dose := round(stats::rnorm(.N, 10 + 3 * case, 2), 1)]
    dt[]
}


### * Status coding

test_that("a status coded 1 and 2 is counted as Surv() codes it", {

    dt <- lung_data()
    dead <- dt$status == 2

    result <- fit(dt, "Surv(time, status)", c("age", "sex_f"),
                  model_type = "coxph")
    raw <- attr(result, "raw_data")
    counts <- attr(result, "analysis_counts")

    expect_equal(raw$events[raw$variable == "sex_f"],
                 counts_by_level(dt, "sex_f", "dead")$events)
    expect_equal(counts$events_supplied, sum(dead))
    expect_equal(counts$events_analyzed, sum(dead))

    plot <- suppressMessages(suppressWarnings(
        coxforest(attr(result, "model"), data = dt)))
    rows <- attr(plot, "table_data")
    expect_equal(rows$Events[rows$var == "sex_f"],
                 counts_by_level(dt, "sex_f", "dead")$events)
    expect_equal(rows$Events[rows$var == "age"], sum(dead))
})


test_that("an event given as an expression or by name is counted", {

    dt <- lung_data()
    dead <- dt$status == 2
    expected <- counts_by_level(dt, "sex_f", "dead")$events

    expression_model <- survival::coxph(survival::Surv(time, status == 2) ~ sex_f,
                                        data = dt)
    named_model <- survival::coxph(survival::Surv(time = time, event = status) ~ sex_f,
                                   data = dt)

    for (model in list(expression_model, named_model)) {
        raw <- m2dt(dt, model)
        expect_equal(raw$events, expected)
    }
})


test_that("observations whose status Surv() codes as missing are not counted", {

    ## Surv() codes a status of 9 among statuses of 0 and 1 as missing, with
    ## a warning, and the model omits those observations; they were counted
    ## in the observations of each level and of the analysis
    dt <- lung_data()
    dt[, status01 := as.integer(status == 2)]
    dt[seq(1L, .N, by = 10L), status01 := 9L]
    kept <- dt[status01 != 9L]
    expected <- counts_by_level(kept, "sex_f", "status01")

    model <- suppressWarnings(
        survival::coxph(survival::Surv(time, status01) ~ age + sex_f, data = dt)
    )
    expect_equal(model$n, nrow(kept))

    raw <- m2dt(dt, model)
    expect_equal(raw$n[raw$variable == "sex_f"], expected$n)
    expect_equal(raw$events[raw$variable == "sex_f"], expected$events)

    rows <- attr(suppressMessages(suppressWarnings(coxforest(model, data = dt))),
                 "table_data")
    expect_equal(rows$N[rows$var == "sex_f"], expected$n)
    expect_equal(rows$N[rows$var == "age"], nrow(kept))

    result <- suppressWarnings(fit(dt, "Surv(time, status01)", c("age", "sex_f"),
                                   model_type = "coxph"))
    counts <- attr(result, "analysis_counts")
    expect_equal(counts$n_analyzed, nrow(kept))
    expect_equal(counts$events_analyzed, sum(kept$status01))

    screened <- suppressWarnings(uniscreen(dt, "Surv(time, status01)",
                                           c("age", "sex_f"), model_type = "coxph",
                                           parallel = FALSE))
    expect_true(all(attr(screened, "analysis_counts")$n_analyzed == nrow(kept)))
})


test_that("the screening and full-fit counts code the status as Surv() does", {

    dt <- lung_data()
    dead <- sum(dt$status == 2)

    screened <- uniscreen(dt, "Surv(time, status)", c("age", "sex_f"),
                          model_type = "coxph", parallel = FALSE)
    full <- suppressMessages(
        fullfit(dt, "Surv(time, status)", c("age", "sex_f"),
                model_type = "coxph", parallel = FALSE)
    )

    for (counts in list(attr(screened, "analysis_counts"),
                        attr(full, "analysis_counts"))) {
        expect_equal(counts$events_supplied, dead)
        expect_true(all(counts$events_analyzed == dead))
    }
})


### * Forms of the response

test_that("the events of a counting-process response are its event argument", {

    ## The second variable of Surv(start, stop, event) is the stop time
    dt <- data.table::as.data.table(survival::heart)
    dt[, surgery_f := factor(surgery, labels = c("No", "Yes"))]

    model <- survival::coxph(survival::Surv(start, stop, event) ~ age + surgery_f,
                             data = dt)
    raw <- m2dt(dt, model)
    expect_equal(raw$events[raw$variable == "surgery_f"],
                 counts_by_level(dt, "surgery_f", "event")$events)

    rows <- attr(suppressMessages(suppressWarnings(coxforest(model, data = dt))),
                 "table_data")
    expect_equal(rows$Events[rows$var == "surgery_f"],
                 counts_by_level(dt, "surgery_f", "event")$events)

    result <- fit(dt, "Surv(start, stop, event)", c("age", "surgery_f"),
                  model_type = "coxph")
    expect_equal(attr(result, "analysis_counts")$events_supplied, sum(dt$event))
})


test_that("the events of a clogit model are its cases", {

    ## clogit() writes its response as Surv(rep(1, N), case); the second
    ## variable of its formula is the first predictor
    dt <- case_control_data()
    model <- survival::clogit(case ~ dose + exposure + survival::strata(set),
                              data = dt)

    raw <- m2dt(dt, model)
    expect_equal(raw$events[raw$variable == "exposure"],
                 counts_by_level(dt, "exposure", "case")$events)

    rows <- attr(suppressMessages(suppressWarnings(coxforest(model, data = dt))),
                 "table_data")
    expect_equal(rows$Events[rows$var == "dose"], sum(dt$case))
    expect_equal(rows$Events[rows$var == "exposure"],
                 counts_by_level(dt, "exposure", "case")$events)
})


test_that("the events of a coxme model are coded as Surv() codes them", {

    skip_if_not_installed("coxme")

    dt <- lung_data()[!is.na(inst)]
    dead <- dt$status == 2
    model <- coxme::coxme(survival::Surv(time, status) ~ sex_f + (1 | inst),
                          data = dt)

    raw <- m2dt(dt, model)
    rows <- attr(suppressMessages(suppressWarnings(coxforest(model, data = dt))),
                 "table_data")

    expect_equal(raw$events, counts_by_level(dt, "sex_f", "dead")$events)
    expect_equal(rows$Events, counts_by_level(dt, "sex_f", "dead")$events)
})


test_that("the event argument is read from the Surv() call", {

    read <- function(text) {
        expr <- summata:::surv_event_expression(text)
        if (is.null(expr)) NULL else deparse(expr)
    }

    expect_identical(read("Surv(time, status)"), "status")
    expect_identical(read("Surv(start, stop, event)"), "event")
    expect_identical(read("Surv(rep(1, 150L), case)"), "case")
    expect_identical(read("Surv(time = t, event = s)"), "s")
    expect_identical(read("survival::Surv(time, status == 2)"), "status == 2")
    expect_null(read("outcome"))
})
