## desctable() lists rows and columns in the order of the data rather than in
## sorted order: variables and statistics in the order requested, the
## categories of a factor in the order of its levels, and the categories of a
## character or logical variable, like the groups of such a grouping variable,
## in order of first appearance, with the missing category last.

order_data <- function() {
    data.frame(
        grp = c("y", "x", "y", "z", "x", "y", NA, "z", "x", "y", "z", "x"),
        num_grp = c(2, 1, 2, 3, 1, 2, 3, 3, 1, 2, 3, 1),
        chr = c("b", NA, "c", "b", "a", "a", "c", "b", "a", "c", "b", "a"),
        lgl = c(TRUE, FALSE, TRUE, NA, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE,
                FALSE, TRUE),
        fct = factor(rep(c("lo", "hi"), 6), levels = c("lo", "mid", "hi")),
        num = c(5, 3, 4, 8, 1, 2, 6, 7, 9, 10, 11, 12),
        stringsAsFactors = FALSE
    )
}

## The Group entries of one variable's rows
variable_groups <- function(result, variable) {
    first <- which(result$Variable == variable)
    following <- which(result$Variable != "" & seq_len(nrow(result)) > first)
    last <- if (length(following) > 0L) following[1] - 1L else nrow(result)
    result$Group[first:last]
}


test_that("categories of character and logical variables appear in order of appearance", {

    result <- desctable(order_data(), variables = c("chr", "lgl"))

    expect_identical(variable_groups(result, "chr"), c("b", "c", "a"))
    expect_identical(variable_groups(result, "lgl"), c("TRUE", "FALSE"))
})


test_that("categories of a factor follow its levels, including empty levels", {

    result <- desctable(order_data(), variables = "fct")

    expect_identical(variable_groups(result, "fct"), c("lo", "mid", "hi"))
})


test_that("the missing category appears last", {

    ## The first missing value of chr precedes the first "c" and "a"
    result <- desctable(order_data(), variables = "chr", na_include = TRUE,
                        na_label = "Missing")

    expect_identical(variable_groups(result, "chr"), c("b", "c", "a", "Missing"))
})


test_that("group columns appear in order of appearance", {

    df <- order_data()
    for (by in c("grp", "num_grp", "lgl")) {
        result <- desctable(df, by = by, variables = "num")
        groups <- setdiff(names(result), c("Variable", "Group", "Total", "p-value"))
        expect_identical(groups, as.character(unique(stats::na.omit(df[[by]]))),
                         info = by)
    }
})


test_that("group columns of a factor follow its levels, including empty levels", {

    df <- order_data()
    df$fct_grp <- df$fct
    result <- desctable(df, by = "fct_grp", variables = "num")

    groups <- setdiff(names(result), c("Variable", "Group", "Total", "p-value"))
    expect_identical(groups, c("lo", "mid", "hi"))
    expect_identical(result$mid[result$Variable == "N"], "0")
})


test_that("variables and statistics appear in the order requested", {

    result <- desctable(order_data(), variables = c("num", "chr"),
                        stats_continuous = c("range", "mean_sd"))

    expect_identical(result$Variable[result$Variable != ""],
                     c("N", "num", "chr"))
    expect_identical(variable_groups(result, "num"),
                     c("Range", "Mean ± SD"))
})


test_that("groups are derived in display order and located once", {

    ## Distinct values in order of first appearance, without missing values
    expect_identical(summata:::derive_group_levels(c("b", NA, "a", "b")),
                     c("b", "a"))
    expect_identical(summata:::derive_group_levels(c(2, 1, NA, 2)), c(2, 1))

    ## The levels of a factor, including levels without observations
    grp <- factor(c("b", "a", NA, "b"), levels = c("c", "b", "a"))
    expect_identical(summata:::derive_group_levels(grp), c("c", "b", "a"))

    ## Row positions per group, in the order of the groups; a missing group
    ## belongs to none
    expect_identical(summata:::derive_group_index(grp, c("c", "b", "a")),
                     list(integer(0), c(1L, 4L), 2L))
})
