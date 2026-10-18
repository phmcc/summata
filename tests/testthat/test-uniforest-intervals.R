## uniforest() draws the intervals of a uniscreen() table, which are computed
## at the table's confidence level. Asked for another level, it showed that
## level in the header but drew the table's intervals; it now computes Wald
## intervals at the requested level, as multiforest() does.

screened_fixture <- function() {
    uniscreen(fresh_clintrial(), "readmission_30d", c("age", "sex", "stage"),
              model_type = "glm", parallel = FALSE)
}

drawn_rows <- function(p) {
    tbl <- attr(p, "table_data")
    tbl[!is.na(tbl$conf_low)]
}

plot_labels <- function(p) {
    unlist(lapply(p$layers, function(layer) as.character(layer$aes_params$label)))
}


test_that("uniforest draws the table's intervals at the table's level", {

    screened <- screened_fixture()
    expect_identical(attr(screened, "conf_level"), 0.95)

    raw <- attr(screened, "raw_data")
    raw <- raw[!is.na(raw$se)]
    tbl <- drawn_rows(suppressMessages(uniforest(screened)))

    expect_equal(tbl$conf_low, raw$coef_lower)
    expect_equal(tbl$conf_high, raw$coef_upper)
})


test_that("uniforest computes intervals at a requested level", {

    screened <- screened_fixture()
    raw <- attr(screened, "raw_data")
    raw <- raw[!is.na(raw$se)]

    p <- suppressMessages(uniforest(screened, conf_level = 0.9))
    tbl <- drawn_rows(p)
    z <- stats::qnorm(0.95)

    expect_equal(tbl$estimate, raw$coefficient)
    expect_equal(tbl$conf_low, raw$coefficient - z * raw$se)
    expect_equal(tbl$conf_high, raw$coefficient + z * raw$se)

    ## The header states the level of the drawn intervals
    labels <- plot_labels(p)
    expect_gt(length(labels), 0L)
    expect_true(any(grepl("90% CI", labels, fixed = TRUE)))
    expect_false(any(grepl("95% CI", labels, fixed = TRUE)))
})


test_that("uniforest reads the level of a result without the attribute from its header", {

    ## Results of earlier versions carry no "conf_level" attribute
    screened <- screened_fixture()
    data.table::setattr(screened, "conf_level", NULL)
    raw <- attr(screened, "raw_data")
    raw <- raw[!is.na(raw$se)]

    tbl <- drawn_rows(suppressMessages(uniforest(screened)))
    expect_equal(tbl$conf_low, raw$coef_lower)
    expect_equal(tbl$conf_high, raw$coef_upper)

    tbl <- drawn_rows(suppressMessages(uniforest(screened, conf_level = 0.9)))
    expect_equal(tbl$conf_low, raw$coefficient - stats::qnorm(0.95) * raw$se)
})


test_that("uniforest computes t intervals for linear models at a requested level", {

    ## The table holds the exact t intervals of each linear model, so that
    ## intervals at another level use the t quantile on its residual degrees
    ## of freedom
    df <- fresh_clintrial()
    screened <- uniscreen(df, "los_days", c("age", "stage"), model_type = "lm",
                          parallel = FALSE)
    raw <- attr(screened, "raw_data")
    raw <- raw[!is.na(raw$se)]

    tbl <- drawn_rows(suppressMessages(uniforest(screened, conf_level = 0.9)))

    for (pred in c("age", "stage")) {
        model <- stats::lm(stats::reformulate(pred, "los_days"), data = df)
        expected <- stats::confint(model, level = 0.9)[-1, , drop = FALSE]
        rows <- raw$predictor == pred
        expect_equal(tbl$conf_low[rows], unname(expected[, 1]), info = pred)
        expect_equal(tbl$conf_high[rows], unname(expected[, 2]), info = pred)
    }
})


test_that("uniforest keeps Wald intervals for linear models screened with them", {

    df <- fresh_clintrial()
    screened <- uniscreen(df, "los_days", "age", model_type = "lm",
                          conf_method = "wald", parallel = FALSE)
    raw <- attr(screened, "raw_data")

    tbl <- drawn_rows(suppressMessages(uniforest(screened, conf_level = 0.9)))

    expect_equal(tbl$conf_low, raw$coefficient - stats::qnorm(0.95) * raw$se)
})
