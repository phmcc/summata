## Factor levels are drawn in the order the model established, not in
## alphabetical order. The distinction is visible only where the two differ,
## so the fixtures below use levels that are deliberately not alphabetical.

order_fixture <- function() {
    dt <- data.table::as.data.table(fresh_clintrial())
    ## "Well" < "Moderately" < "Poorly" clinically, but alphabetically the
    ## order is Moderately, Poorly, Well
    dt[]
}

plot_levels <- function(p, variable) {
    tbl <- attr(p, "table_data")
    rows <- trimws(tbl[["display_text"]])
    ## Levels follow their variable's header row
    start <- which(rows == variable)
    if (length(start) == 0) return(character(0))
    after <- rows[(start[1] + 1):length(rows)]
    headers <- trimws(tbl[["display_text"]])[tbl[["row_type"]] == "header"]
    keep <- character(0)
    for (r in after) {
        if (r %in% headers) break
        keep <- c(keep, r)
    }
    return(keep)
}


test_that("uniforest preserves factor level order when indenting", {

    dt <- order_fixture()
    skip_if(!"grade" %in% names(dt), "grade not present")

    screened <- uniscreen(
        data = dt,
        outcome = "readmission_30d",
        predictors = c("age", "grade"),
        model_type = "glm",
        family = "binomial",
        parallel = FALSE
    )

    p <- suppressMessages(suppressWarnings(
        uniforest(screened, indent_groups = TRUE)))

    drawn <- plot_levels(p, "grade")
    skip_if(length(drawn) < 2, "grade levels not found in table_data")

    ## The drawn order matches the factor levels, not sort()
    expect_identical(drawn, levels(dt$grade))
    expect_false(identical(drawn, sort(levels(dt$grade))))
})


test_that("indenting does not change which levels are drawn", {

    dt <- order_fixture()

    screened <- uniscreen(
        data = dt,
        outcome = "readmission_30d",
        predictors = c("age", "grade"),
        model_type = "glm",
        family = "binomial",
        parallel = FALSE
    )

    flat <- suppressMessages(suppressWarnings(
        uniforest(screened, indent_groups = FALSE)))
    indented <- suppressMessages(suppressWarnings(
        uniforest(screened, indent_groups = TRUE)))

    flat_n <- attr(flat, "table_data")[["n_formatted"]]
    ind_tbl <- attr(indented, "table_data")
    ind_n <- ind_tbl[["n_formatted"]][ind_tbl[["row_type"]] != "header"]

    ## Header rows aside, the same rows appear in the same order
    expect_identical(flat_n, ind_n)
})
