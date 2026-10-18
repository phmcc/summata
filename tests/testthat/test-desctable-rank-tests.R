## desctable() computes the Kruskal-Wallis and large-sample Wilcoxon rank sum
## p-values from a single sort rather than through rank() and table(), which
## dominate the running time of kruskal.test() and wilcox.test() on large
## samples. The p-values must agree with the stats implementations, including
## tie corrections, unequal group sizes, unused factor levels, and degenerate
## data.

rank_test_designs <- function() {
    set.seed(20261009)
    values <- list(
        continuous = function(n) rnorm(n, 50, 10),
        rounded = function(n) round(rnorm(n, 50, 10), 1),
        integer = function(n) rpois(n, 3),
        heavy_ties = function(n) sample(c(1, 2, 2.5, 7), n, replace = TRUE),
        negative = function(n) -abs(rnorm(n))
    )
    designs <- list()
    for (kind in names(values)) {
        for (n in c(30, 99, 100, 101, 2500)) {
            for (k in 2:4) {
                g <- factor(sample(letters[seq_len(k)], n, replace = TRUE, prob = seq_len(k)),
                            levels = letters[seq_len(k + 1)])
                designs[[length(designs) + 1L]] <- list(x = values[[kind]](n), g = g)
            }
        }
    }
    designs
}


test_that("the Kruskal-Wallis p-value matches kruskal.test()", {

    for (d in rank_test_designs()) {
        expect_equal(summata:::kruskal_wallis_pvalue(d$x, d$g),
                     stats::kruskal.test(d$x, d$g)$p.value,
                     tolerance = 1e-12)
    }
})


test_that("the Wilcoxon rank sum p-value matches wilcox.test()", {

    for (d in rank_test_designs()) {
        g <- droplevels(d$g)
        if (nlevels(g) != 2L) next
        expect_equal(suppressWarnings(summata:::rank_sum_pvalue(d$x, g)),
                     suppressWarnings(stats::wilcox.test(d$x ~ g)$p.value),
                     tolerance = 1e-12)
    }
})


test_that("degenerate data give the same result as the stats functions", {

    tied <- rep(5, 200)
    g2 <- rep(c("a", "b"), 100)
    g3 <- rep(c("a", "b", "c"), length.out = 200)

    expect_identical(summata:::kruskal_wallis_pvalue(tied, g3),
                     stats::kruskal.test(tied, factor(g3))$p.value)
    expect_identical(summata:::rank_sum_pvalue(tied, g2),
                     suppressWarnings(stats::wilcox.test(tied ~ factor(g2))$p.value))
    expect_error(summata:::kruskal_wallis_pvalue(1:10, rep("a", 10)), "same group")
    expect_error(summata:::rank_sum_pvalue(1:30, g3[1:30]), "exactly 2 levels")
})


test_that("average_ranks matches rank() and counts ties", {

    x <- c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5)
    ranks <- summata:::average_ranks(x)

    expect_identical(as.vector(ranks), rank(x))
    expect_identical(attr(ranks, "ties"), as.vector(table(x)))
})


test_that("desctable reports the rank test p-values", {

    df <- fresh_clintrial()

    three <- desctable(df, by = "treatment", variables = "age")
    expect_equal(attr(three, "raw_data")$p_value[1],
                 stats::kruskal.test(age ~ treatment, data = df)$p.value,
                 tolerance = 1e-12)

    two <- desctable(df, by = "sex", variables = "age")
    expect_equal(attr(two, "raw_data")$p_value[1],
                 stats::wilcox.test(age ~ sex, data = df)$p.value,
                 tolerance = 1e-12)

    ## Additional arguments reach wilcox.test() itself
    shifted <- desctable(df, by = "sex", variables = "age", test_continuous = "wrs",
                         correct = FALSE)
    expect_equal(attr(shifted, "raw_data")$p_value[1],
                 stats::wilcox.test(age ~ sex, data = df, correct = FALSE)$p.value)
})
