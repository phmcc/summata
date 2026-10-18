### * Variable sorting

#' Process variable wrapper
#' 
#' Routes variable processing to appropriate handler based on variable type
#' (continuous, categorical, or survival). Returns both formatted display
#' strings and raw numeric values.
#' 
#' @param data Data.table containing the variable.
#' @param var Character string naming the variable to process.
#' @param group_var Optional character string naming the grouping variable.
#' @param stats_continuous Character vector of statistics for continuous variables.
#' @param stats_categorical Character vector of statistics for categorical variables.
#' @param digits Integer number of decimal places for continuous statistics.
#' @param conf_level Numeric confidence level for survival confidence intervals.
#' @param na_include Logical whether to include missing values as a category.
#' @param na_label Character string label for missing values.
#' @param test Logical whether to perform statistical tests.
#' @param test_continuous Character string specifying test type for continuous variables.
#' @param test_categorical Character string specifying test type for categorical variables.
#' @param total Logical or character controlling total column display.
#' @param total_label Character string label for total column.
#' @param labels Named character vector of variable labels.
#' @param na_percent Logical whether to include NA in percentage denominators.
#' @param p_per_stat Logical whether to show separate \emph{p}-values per statistic for
#'   continuous variables. Default \code{FALSE} for better performance.
#' @param marks List with \code{big.mark} and \code{decimal.mark} as returned
#'   by \code{\link{resolve_number_marks}}.
#' @param groups Optional vector of groups, as returned by
#'   \code{derive_group_levels()}. Derived from the data when \code{NULL}.
#' @param group_index Optional list of row positions for each group, as
#'   returned by \code{derive_group_index()}. Derived from the data when
#'   \code{NULL}.
#' @param ... Additional arguments passed to test functions.
#' @return List with 'formatted' and 'raw' data.table components.
#' @keywords internal
process_variable <- function(data, var, group_var = NULL, 
                             stats_continuous, stats_categorical,
                             digits, conf_level = 0.95, na_include, na_label,
                             test, test_continuous, test_categorical,
                             total, total_label, labels, na_percent, 
                             p_per_stat = FALSE, marks = NULL,
                             groups = NULL, group_index = NULL, ...) {
    
    ## Get variable label
    var_label <- if (!is.null(labels) && var %chin% names(labels)) {
                     labels[var]
                 } else {
                     var
                 }
    
    ## Determine variable type and process accordingly
    if (grepl("^Surv\\(", var)) {
        return(process_survival(
            data = data,
            var = var,
            var_label = var_label,
            group_var = group_var,
            digits = digits,
            conf_level = conf_level,
            na_include = na_include,
            na_label = na_label,
            test = test,
            total = total,
            total_label = total_label,
            marks = marks,
            groups = groups,
            group_index = group_index,
            ...
        ))
    } else if (is.numeric(data[[var]])) {
        return(process_continuous(
            data = data,
            var = var,
            var_label = var_label,
            group_var = group_var,
            stats = stats_continuous,
            digits = digits,
            na_include = na_include,
            na_label = na_label,
            test = test,
            test_type = test_continuous,
            total = total,
            total_label = total_label,
            p_per_stat = p_per_stat,
            marks = marks,
            groups = groups,
            group_index = group_index,
            ...
        ))
    } else {
        return(process_categorical(
            data = data,
            var = var,
            var_label = var_label,
            group_var = group_var,
            stats = stats_categorical,
            na_include = na_include,
            na_label = na_label,
            test = test,
            test_type = test_categorical,
            total = total,
            total_label = total_label,
            na_percent = na_percent,
            marks = marks,
            groups = groups,
            ...
        ))
    }
}


### * Group structure

#' Derive the groups of a grouping variable
#'
#' Returns the groups in display order: the levels of a factor, including
#' levels without observations, or otherwise the distinct non-missing values
#' in order of first appearance.
#'
#' @param grp_vec Grouping vector.
#' @return Vector of groups.
#' @keywords internal
derive_group_levels <- function(grp_vec) {
    if (is.factor(grp_vec)) {
        levels(grp_vec)
    } else {
        groups <- unique(grp_vec)
        groups[!is.na(groups)]
    }
}


#' Locate the observations of each group
#'
#' Finds the row positions of each group once, so that the summaries of
#' every variable can subset them instead of comparing the full grouping
#' vector with each group again.
#'
#' @param grp_vec Grouping vector.
#' @param groups Vector of groups, as returned by
#'   \code{derive_group_levels()}.
#' @return List with one integer vector of row positions per group, in the
#'   order of \code{groups}. Observations with a missing group appear in none.
#' @keywords internal
derive_group_index <- function(grp_vec, groups) {
    lapply(groups, function(g) which(grp_vec == g))
}


### * Test arguments

#' Name the arguments that each descriptive test accepts
#'
#' Returns the names of the arguments that \code{desctable()} forwards
#' from \code{...} to a given statistical test. A single call can summarize
#' continuous, categorical, and time-to-event variables, so an argument
#' intended for one test (\emph{e.g.,} \code{var.equal} for the
#' \emph{t}-test) must not reach a test that does not accept it. The lists are
#' explicit because the test functions are S3 generics whose formal arguments
#' do not include the arguments of their methods.
#'
#' @param test_type Character string naming the test: \code{"t"},
#'   \code{"wrs"}, \code{"aov"}, \code{"kwt"}, \code{"chisq"},
#'   \code{"fisher"}, or \code{"survival"}. When \code{NULL}, the names
#'   accepted by any of the tests are returned.
#' @return Character vector of argument names.
#' @keywords internal
test_arg_names <- function(test_type = NULL) {
    accepted <- list(
        t = c("alternative", "mu", "var.equal", "conf.level"),
        wrs = c("alternative", "mu", "exact", "correct", "conf.int",
                "conf.level", "tol.root", "digits.rank"),
        aov = character(0),
        kwt = character(0),
        chisq = c("correct", "simulate.p.value", "B"),
        fisher = c("workspace", "hybrid", "hybridPars", "control", "or",
                   "alternative", "conf.int", "conf.level",
                   "simulate.p.value", "B"),
        survival = c("rho", "timefix")
    )
    if (is.null(test_type)) {
        return(unique(unlist(accepted, use.names = FALSE)))
    }
    test_type <- switch(test_type,
                        "fisher.test" = "fisher",
                        "chisq.test" = "chisq",
                        test_type)
    accepted[[test_type]] %||% character(0)
}


#' Select the test arguments that a given test accepts
#'
#' Subsets the arguments collected from \code{...} in \code{desctable()}
#' to those accepted by one statistical test, as listed by
#' \code{test_arg_names()}.
#'
#' @param test_args Named list of arguments collected from \code{...}.
#' @param test_type Character string naming the test, as in
#'   \code{test_arg_names()}.
#' @return Named list of the arguments accepted by the test (possibly empty).
#' @keywords internal
filter_test_args <- function(test_args, test_type) {
    if (length(test_args) == 0L || is.null(names(test_args))) {
        return(list())
    }
    test_args[names(test_args) %chin% test_arg_names(test_type)]
}


### * Continuous variables

#' Process continuous variable
#' 
#' Calculates descriptive statistics for continuous numeric variables, with
#' optional grouping and statistical testing. Supports multiple summary
#' statistics (mean \eqn{\pm} SD, median [IQR], range) and various hypothesis tests.
#' 
#' @param data Data.table containing the variable.
#' @param var Character string naming the variable to process.
#' @param var_label Character string label for display.
#' @param group_var Optional character string naming the grouping variable.
#' @param stats Character vector of statistics to calculate.
#' @param digits Integer number of decimal places.
#' @param na_include Logical whether to include missing values.
#' @param na_label Character string label for missing values.
#' @param test Logical whether to perform statistical tests.
#' @param test_type Character string specifying test type.
#' @param total Logical or character controlling total column display.
#' @param total_label Character string label for total column.
#' @param p_per_stat Logical. If TRUE, calculate separate \emph{p}-values for each 
#'   statistic type (\emph{e.g.,} t-test for means, Wilcoxon for medians). If \code{FALSE}
#'   (default), calculate a single \emph{p}-value based on the first statistic type.
#' @param groups Optional vector of groups (see \code{process_variable()}).
#' @param group_index Optional list of row positions for each group (see
#'   \code{process_variable()}).
#' @param ... Additional arguments passed to test functions.
#' @return List with 'formatted' and 'raw' data.table components.
#' @keywords internal
process_continuous <- function(data, var, var_label, group_var, stats, digits,
                               na_include, na_label, test, test_type,
                               total, total_label, p_per_stat = FALSE,
                               marks = NULL, groups = NULL, group_index = NULL,
                               ...) {
    
    ## Pre-extract the variable column
    var_vec <- data[[var]]
    not_na <- !is.na(var_vec)
    total_vals <- var_vec[not_na]
    total_n <- length(total_vals)
    
    ## Calculate total statistics
    if (total_n > 0) {
        total_fivenum <- stats::fivenum(total_vals)
        total_stats <- list(
            mean = mean(total_vals),
            sd = stats::sd(total_vals),
            median = total_fivenum[3],
            q1 = total_fivenum[2],
            q3 = total_fivenum[4],
            min = total_fivenum[1],
            max = total_fivenum[5],
            n = total_n
        )
    } else {
        total_stats <- list(mean = NA_real_, sd = NA_real_, median = NA_real_,
                            q1 = NA_real_, q3 = NA_real_, min = NA_real_, 
                            max = NA_real_, n = 0L)
    }
    
    ## Calculate group statistics (if grouped)
    group_stats_list <- NULL
    
    if (!is.null(group_var)) {
        grp_vec <- data[[group_var]]
        if (is.null(groups)) groups <- derive_group_levels(grp_vec)
        if (is.null(group_index)) group_index <- derive_group_index(grp_vec, groups)
        
        ## Compute statistics for each group
        group_stats_list <- lapply(group_index, function(idx) {
            idx <- idx[not_na[idx]]
            if (length(idx) == 0) {
                return(list(mean = NA_real_, sd = NA_real_, median = NA_real_,
                           q1 = NA_real_, q3 = NA_real_, min = NA_real_, 
                           max = NA_real_, n = 0L))
            }
            vals <- var_vec[idx]
            fn <- stats::fivenum(vals)
            list(
                mean = mean(vals),
                sd = stats::sd(vals),
                median = fn[3],
                q1 = fn[2],
                q3 = fn[4],
                min = fn[1],
                max = fn[5],
                n = length(vals)
            )
        })
        names(group_stats_list) <- as.character(groups)
    }
    
    ## Calculate p-values
    p_values <- list()
    p_value_single <- NULL
    if (test && !is.null(group_var) && length(groups) >= 2) {
        grp_vec <- data[[group_var]]
        test_args <- list(...)
        
        if (p_per_stat) {
            ## Calculate separate p-value for each statistic type
            for (stat_type in stats) {
                if (stat_type == "range") {
                    p_values[[stat_type]] <- NULL
                } else {
                    p_values[[stat_type]] <- perform_continuous_test(
                        var_vec, grp_vec, test_type, stat_type, test_args)
                }
            }
        } else {
            ## Calculate single p-value based on first non-range statistic
            first_stat <- stats[stats != "range"][1]
            if (is.na(first_stat)) first_stat <- "median_iqr"
            p_value_single <- perform_continuous_test(
                var_vec, grp_vec, test_type, first_stat, test_args)
        }
    }
    
    ## Build output
    n_stats <- length(stats)
    has_missing <- na_include && any(!not_na)
    n_rows <- n_stats + if (has_missing) 1L else 0L
    
    ## Compute format string
    fmt_str <- paste0("%.", digits, "f")
    
    formatted_list <- vector("list", n_rows)
    raw_list <- vector("list", n_rows)
    
    for (i in seq_along(stats)) {
        stat_type <- stats[i]
        
        formatted_row <- list(
            variable = if (i == 1L) var_label else "",
            level = get_stat_label(stat_type)
        )
        
        raw_row <- list(
            variable = if (i == 1L) var_label else "",
            level = stat_type,
            stat_type = stat_type
        )
        
        ## Add total column
        if (!isFALSE(total)) {
            formatted_row[[total_label]] <- format_continuous_stat(
                total_stats, stat_type, fmt_str, marks)
            raw_row <- add_raw_stats(raw_row, total_label, total_stats, stat_type)
        }
        
        ## Add group columns
        if (!is.null(group_var)) {
            for (g in groups) {
                grp_col <- as.character(g)
                grp_stats <- group_stats_list[[grp_col]]
                formatted_row[[grp_col]] <- format_continuous_stat(
                    grp_stats, stat_type, fmt_str, marks)
                raw_row <- add_raw_stats(raw_row, grp_col, grp_stats, stat_type)
            }
            
            ## Add p-value
            if (test && stat_type != "range") {
                if (p_per_stat && !is.null(p_values[[stat_type]])) {
                    ## Per-statistic p-values
                    formatted_row[["p_value"]] <- p_values[[stat_type]]
                    raw_row[["p_value"]] <- p_values[[stat_type]]
                } else if (!p_per_stat && i == 1L && !is.null(p_value_single)) {
                    ## Single p-value on first row only
                    formatted_row[["p_value"]] <- p_value_single
                    raw_row[["p_value"]] <- p_value_single
                }
            }
        }
        
        formatted_list[[i]] <- data.table::as.data.table(formatted_row)
        raw_list[[i]] <- data.table::as.data.table(raw_row)
    }
    
    ## Missing row
    if (has_missing) {
        n_miss_total <- sum(!not_na)
        miss_formatted <- list(variable = "", level = na_label)
        miss_raw <- list(variable = "", level = na_label, stat_type = "missing")
        
        if (!isFALSE(total)) {
            miss_formatted[[total_label]] <- format_count(n_miss_total, marks)
            miss_raw[[total_label]] <- n_miss_total
        }
        
        if (!is.null(group_var)) {
            for (j in seq_along(groups)) {
                grp_col <- as.character(groups[j])
                n_miss <- sum(!not_na[group_index[[j]]])
                miss_formatted[[grp_col]] <- format_count(n_miss, marks)
                miss_raw[[grp_col]] <- n_miss
            }
        }
        
        formatted_list[[n_rows]] <- data.table::as.data.table(miss_formatted)
        raw_list[[n_rows]] <- data.table::as.data.table(miss_raw)
    }
    
    list(
        formatted = data.table::rbindlist(formatted_list, fill = TRUE),
        raw = data.table::rbindlist(raw_list, fill = TRUE)
    )
}


#' Add raw statistics to row
#' 
#' Appends raw numeric statistics to a data.table row for downstream processing.
#' Used to preserve underlying values alongside formatted display strings.
#' 
#' @param row Data.table row to modify.
#' @param col Character string column name for the statistics.
#' @param stats Named list of numeric statistics (mean, sd, median, \emph{etc.}).
#' @param stat_type Character string indicating which statistic is displayed.
#' @return Modified row (by reference).
#' @keywords internal
add_raw_stats <- function(row, col, stats, stat_type) {
    if (stat_type == "mean_sd") {
        row[[col]] <- stats$mean
        row[[paste0(col, "_sd")]] <- stats$sd
        row[[paste0(col, "_n")]] <- stats$n
    } else if (stat_type == "median_iqr") {
        row[[col]] <- stats$median
        row[[paste0(col, "_q1")]] <- stats$q1
        row[[paste0(col, "_q3")]] <- stats$q3
        row[[paste0(col, "_n")]] <- stats$n
    } else if (stat_type == "median_range") {
        row[[col]] <- stats$median
        row[[paste0(col, "_min")]] <- stats$min
        row[[paste0(col, "_max")]] <- stats$max
        row[[paste0(col, "_n")]] <- stats$n
    } else if (stat_type == "range") {
        row[[col]] <- stats$min
        row[[paste0(col, "_max")]] <- stats$max
        row[[paste0(col, "_n")]] <- stats$n
    }
    row  # Return the modified row
}


#' Perform statistical tests for continuous variables
#' 
#' Conducts hypothesis tests comparing continuous variables across groups.
#' Supports t-tests, Wilcoxon tests, ANOVA, and Kruskal-Wallis tests with
#' automatic selection based on number of groups.
#' 
#' @param var_vec Numeric vector of the continuous variable.
#' @param grp_vec Factor or character vector defining groups.
#' @param test_type Character string: "parametric", "nonparametric", or "auto".
#' @param stat_type Character string indicating primary statistic being tested.
#' @param test_args Named list of additional test arguments from \code{...};
#'   each test receives only the arguments it accepts (see
#'   \code{filter_test_args()}).
#' @return Numeric \emph{p}-value from the hypothesis test.
#' @keywords internal
perform_continuous_test <- function(var_vec, grp_vec, test_type, stat_type,
                                    test_args = list()) {
    ## Remove NAs
    valid <- !is.na(var_vec) & !is.na(grp_vec)
    if (sum(valid) < 3) return(NA_real_)
    
    x <- var_vec[valid]
    g <- grp_vec[valid]
    
    groups <- unique(g)
    n_groups <- length(groups)
    if (n_groups < 2) return(NA_real_)
    
    ## Auto-select test
    if (test_type == "auto") {
        if (grepl("mean", stat_type)) {
            test_type <- if (n_groups == 2) "t" else "aov"
        } else {
            test_type <- if (n_groups == 2) "wrs" else "kwt"
        }
    }
    
    args <- filter_test_args(test_args, test_type)

    tryCatch({
        switch(test_type,
            "t" = do.call(stats::t.test, c(list(x ~ g), args))$p.value,
            "wrs" = rank_sum_pvalue(x, g, args),
            "aov" = {
                fit <- stats::aov(x ~ g)
                summary(fit)[[1]][["Pr(>F)"]][1]
            },
            "kwt" = kruskal_wallis_pvalue(x, g),
            NA_real_
        )
    }, error = function(e) NA_real_)
}


#' Average ranks with tie counts
#'
#' Ranks a numeric vector with ties assigned their average rank, as
#' \code{rank()} does by default, from a single sort with
#' \code{data.table::frank()}. The \eqn{j}th smallest distinct value, shared
#' by \eqn{t_j} observations, occupies sorted positions
#' \eqn{c_{j-1} + 1, \dots, c_j}{c[j-1] + 1, ..., c[j]}, where \eqn{c} is the
#' cumulative sum of the tie counts, so its average rank is
#' \eqn{c_j - (t_j - 1)/2}{c[j] - (t[j] - 1)/2}. Every quantity is an integer
#' or half-integer, so the ranks are exact.
#'
#' @param x Numeric vector without missing values.
#' @return Numeric vector of ranks, with attribute \code{"ties"} holding the
#'   number of observations at each distinct value in increasing order.
#' @keywords internal
average_ranks <- function(x) {
    dense <- data.table::frank(x, ties.method = "dense")
    ties <- tabulate(dense)
    ranks <- (cumsum(ties) - (ties - 1) / 2)[dense]
    attr(ranks, "ties") <- ties
    ranks
}


#' Kruskal-Wallis rank sum test p-value
#'
#' Computes the \emph{p}-value of the Kruskal-Wallis rank sum test with the
#' statistic, tie correction, and degrees of freedom of
#' \code{stats::kruskal.test()}, using ranks and tie counts from
#' \code{average_ranks()}. \code{kruskal.test()} ranks with \code{rank()}
#' and counts ties with \code{table()}, which converts every value to a
#' character string; on large samples these two steps account for nearly all
#' of its running time.
#'
#' @param x Numeric vector without missing values.
#' @param g Grouping vector of the same length without missing values.
#' @return Numeric \emph{p}-value.
#' @keywords internal
kruskal_wallis_pvalue <- function(x, g) {
    ## Integer group codes; groups without observations are dropped, as
    ## kruskal.test() does through factor()
    g <- if (is.factor(g)) as.integer(g) else as.integer(factor(g))
    n_i <- tabulate(g)
    n_i <- n_i[n_i > 0L]
    k <- length(n_i)
    n <- length(x)
    if (k < 2L) stop("all observations are in the same group")
    if (n < 2L) stop("not enough observations")

    ranks <- average_ranks(x)
    ties <- attr(ranks, "ties")

    statistic <- sum(rowsum(as.vector(ranks), g)^2 / n_i)
    statistic <- (12 * statistic / (n * (n + 1)) - 3 * (n + 1)) /
        (1 - sum(ties^3 - ties) / (n^3 - n))
    stats::pchisq(statistic, k - 1L, lower.tail = FALSE)
}


#' Wilcoxon rank sum test p-value
#'
#' Computes the two-sided \emph{p}-value of the Wilcoxon rank sum test for
#' two groups. When either group has 50 or more observations,
#' \code{stats::wilcox.test()} uses the normal approximation with
#' continuity correction and tie correction by default; that computation is
#' reproduced here with ranks and tie counts from \code{average_ranks()},
#' because \code{wilcox.test()} ranks with \code{rank()} and counts ties with
#' \code{table()}, which converts every rank to a character string; on large
#' samples these two steps dominate its running time. Smaller samples, for
#' which \code{wilcox.test()} computes the exact distribution by default, and
#' calls with additional test arguments are passed to \code{wilcox.test()}
#' itself, as are data without variation, whose \emph{p}-value differs
#' between versions of \R.
#'
#' \code{wilcox.test()} ranks the values rounded to the significant digits of
#' its argument \code{digits.rank}, which defaults to 7 in recent versions of
#' \R and is absent from earlier ones, which do not round. The values are
#' rounded as the installed \code{wilcox.test()} rounds them by default, so
#' that the ties, and the \emph{p}-value, are those it finds.
#'
#' @param x Numeric vector without missing values.
#' @param g Grouping vector of the same length without missing values, with
#'   exactly two groups.
#' @param test_args Named list of additional arguments for
#'   \code{wilcox.test()}.
#' @return Numeric \emph{p}-value.
#' @keywords internal
rank_sum_pvalue <- function(x, g, test_args = list()) {
    g <- factor(g)
    if (nlevels(g) != 2L) stop("grouping factor must have exactly 2 levels")
    first <- g == levels(g)[1L]
    n_x <- as.double(sum(first))
    n_y <- as.double(length(x)) - n_x

    if (length(test_args) > 0L || (n_x < 50 && n_y < 50)) {
        return(do.call(stats::wilcox.test,
                       c(list(x = x[first], y = x[!first]), test_args))$p.value)
    }

    ## Rounded as wilcox.test() rounds before ranking
    digits_rank <- wilcox_digits_rank()
    ranks <- average_ranks(if (is.finite(digits_rank)) signif(x, digits_rank) else x)
    ties <- attr(ranks, "ties")

    z <- sum(ranks[first]) - n_x * (n_x + 1) / 2 - n_x * n_y / 2
    sigma <- sqrt((n_x * n_y / 12) *
                  ((n_x + n_y + 1) -
                   sum(ties^3 - ties) / ((n_x + n_y) * (n_x + n_y - 1))))
    if (!is.finite(sigma) || sigma == 0) {
        return(stats::wilcox.test(x = x[first], y = x[!first])$p.value)
    }
    z <- (z - sign(z) * 0.5) / sigma
    2 * min(stats::pnorm(z), stats::pnorm(z, lower.tail = FALSE))
}


#' Significant digits to which wilcox.test() rounds before ranking
#'
#' @return The default of the \code{digits.rank} argument of
#'   \code{stats::wilcox.test()}, or \code{Inf} where it has none.
#' @keywords internal
wilcox_digits_rank <- function() {
    default <- formals(utils::getS3method("wilcox.test", "default"))$digits.rank
    value <- tryCatch(eval(default, baseenv()), error = function(e) NULL)
    if (is.numeric(value) && length(value) == 1L && !is.na(value)) value else Inf
}


#' Format continuous statistic for display
#'
#' Converts numeric summary statistics into formatted display strings following
#' standard conventions (mean \eqn{\pm} SD, median [IQR], range, \emph{etc.}).
#' 
#' @param stats Named list of numeric statistics (mean, sd, median, q1, q3, min, max).
#' @param stat_type Character string: "mean_sd", "median_iqr", "median_range", or "range".
#' @param fmt_str Character string format specification for sprintf.
#' @param marks List with \code{big.mark} and \code{decimal.mark} as returned
#'   by \code{\link{resolve_number_marks}}.
#' @return Character string with formatted statistic.
#' @keywords internal
format_continuous_stat <- function(stats, stat_type, fmt_str, marks) {
    if (stats$n == 0) return("")
    
    switch(stat_type,
        "mean_sd" = paste0(
            format_num(stats$mean, fmt_str, marks), " \u00B1 ",
            format_num(stats$sd, fmt_str, marks)
        ),
        "median_iqr" = {
            sep <- resolve_separator(stats$q1, stats$q3, marks)
            paste0(
                format_num(stats$median, fmt_str, marks), " [",
                format_num(stats$q1, fmt_str, marks), sep,
                format_num(stats$q3, fmt_str, marks), "]"
            )
        },
        "median_range" = {
            sep <- resolve_separator(stats$min, stats$max, marks)
            paste0(
                format_num(stats$median, fmt_str, marks), " (",
                format_num(stats$min, fmt_str, marks), sep,
                format_num(stats$max, fmt_str, marks), ")"
            )
        },
        "range" = {
            sep <- resolve_separator(stats$min, stats$max, marks)
            paste0(format_num(stats$min, fmt_str, marks), sep, 
                   format_num(stats$max, fmt_str, marks))
        },
        ""
    )
}


#' Get display label for statistic type
#' 
#' Converts internal statistic type codes to formatted display labels for
#' table column headers.
#' 
#' @param stat_type Character string: "mean_sd", "median_iqr", "median_range",
#'   "range", "n_miss", or custom type.
#' @return Character string with formatted label (\emph{e.g.,} "Mean \eqn{\pm} SD").
#' @keywords internal
get_stat_label <- function(stat_type) {
    switch(stat_type,
           "mean_sd" = "Mean \u00B1 SD",
           "median_iqr" = "Median [IQR]",
           "median_range" = "Median (Range)",
           "range" = "Range",
           "n_miss" = "Missing",
           stat_type
    )
}


### * Categorical variables

#' Process categorical variable
#' 
#' Calculates frequency and percentage statistics for categorical variables,
#' with optional grouping and chi-square/Fisher's exact testing. Handles
#' factor levels, missing values, and custom labeling.
#' 
#' @param data Data.table containing the variable.
#' @param var Character string naming the variable to process.
#' @param var_label Character string label for display.
#' @param group_var Optional character string naming the grouping variable.
#' @param stats Character vector of statistics to calculate.
#' @param na_include Logical whether to include missing values as a category.
#' @param na_label Character string label for missing values.
#' @param test Logical whether to perform statistical tests.
#' @param test_type Character string specifying test type.
#' @param total Logical or character controlling total column display.
#' @param total_label Character string label for total column.
#' @param na_percent Logical whether to include NA in percentage denominators.
#' @param groups Optional vector of groups (see \code{process_variable()}).
#' @param ... Additional arguments passed to test functions.
#' @return List with 'formatted' and 'raw' data.table components.
#' @keywords internal
process_categorical <- function(data, var, var_label, group_var, stats,
                                na_include, na_label, test, test_type,
                                total, total_label, na_percent,
                                marks = NULL, groups = NULL, ...) {
    
    ## Pre-extract vectors
    var_vec <- data[[var]]
    
    ## Get levels
    if (is.factor(var_vec)) {
        levels_to_show <- levels(var_vec)
    } else {
        levels_to_show <- unique(var_vec)
        levels_to_show <- levels_to_show[!is.na(levels_to_show)]
    }
    
    ## Add NA level if requested
    has_na <- anyNA(var_vec)
    if (na_include && has_na) {
        levels_to_show <- c(levels_to_show, NA)
    }
    
    n_levels <- length(levels_to_show)
    
    ## Pre-compute all counts with single table() call
    use_na <- if (na_include) "ifany" else "no"
    
    if (!is.null(group_var)) {
        grp_vec <- data[[group_var]]
        if (is.null(groups)) groups <- derive_group_levels(grp_vec)
        
        ## Single table call for cross-tabulation
        tab <- table(var_vec, grp_vec, useNA = use_na)
        grp_totals <- colSums(tab)
        
        ## The total column counts every observation, as the N row and the
        ## continuous and survival summaries do, including observations whose
        ## group is missing (which the cross-tabulation omits unless
        ## na_include = TRUE)
        total_tab <- table(var_vec, useNA = use_na)
        total_na <- is.na(names(total_tab))

        ## Calculate denominators
        if (na_percent) {
            total_denom <- sum(total_tab)
            grp_denoms <- grp_totals
        } else {
            ## Exclude NA row from denominators
            non_na_rows <- !is.na(rownames(tab))
            total_denom <- sum(total_tab[!total_na])
            grp_denoms <- colSums(tab[non_na_rows, , drop = FALSE])
        }
        
        ## Calculate p-value once
        p_value <- if (test) {
            perform_categorical_test(tab, test_type, list(...), var_name = var)
                   } else {
            NULL
        }
        
        ## Build output
        formatted_list <- vector("list", n_levels)
        raw_list <- vector("list", n_levels)
        
        for (i in seq_along(levels_to_show)) {
            lvl <- levels_to_show[i]
            lvl_label <- if (is.na(lvl)) na_label else as.character(lvl)
            lvl_char <- if (is.na(lvl)) NA_character_ else as.character(lvl)
            
            level_formatted <- list(
                variable = if (i == 1L) var_label else "",
                level = lvl_label
            )
            
            level_raw <- list(
                variable = if (i == 1L) var_label else "",
                level = lvl_label,
                stat_type = "category"
            )
            
            ## Get row from table
            if (is.na(lvl)) {
                n_total <- sum(total_tab[total_na])
            } else {
                n_total <- total_tab[lvl_char]
                if (is.na(n_total)) n_total <- 0L
            }
            
            ## Add total column
            if (!isFALSE(total)) {
                ## NA row: use total denom if na_percent, otherwise count only
                if (is.na(lvl)) {
                    denom <- if (na_percent) total_denom else sum(total_tab)
                    na_stat <- if (na_percent) stats else "n"
                } else {
                    denom <- total_denom
                    na_stat <- stats
                }
                level_formatted[[total_label]] <- format_categorical_stat(
                    as.integer(n_total), as.integer(denom), na_stat, marks)
                level_raw[[total_label]] <- as.integer(n_total)
                level_raw[[paste0(total_label, "_total")]] <- as.integer(denom)
            }
            
            ## Add group columns
            for (g in groups) {
                grp_col <- as.character(g)
                
                if (is.na(lvl)) {
                    n <- sum(tab[is.na(rownames(tab)), grp_col])
                    denom <- if (na_percent) grp_denoms[grp_col] else grp_totals[grp_col]
                    na_stat <- if (na_percent) stats else "n"
                } else {
                    n <- tab[lvl_char, grp_col]
                    if (is.na(n)) n <- 0L
                    denom <- grp_denoms[grp_col]
                    na_stat <- stats
                }
                
                level_formatted[[grp_col]] <- format_categorical_stat(
                    as.integer(n), as.integer(denom), na_stat, marks)
                level_raw[[grp_col]] <- as.integer(n)
                level_raw[[paste0(grp_col, "_total")]] <- as.integer(denom)
            }
            
            ## Add p-value to first level only
            if (i == 1L && !is.null(p_value)) {
                level_formatted[["p_value"]] <- p_value
                level_raw[["p_value"]] <- p_value
            }
            
            formatted_list[[i]] <- data.table::as.data.table(level_formatted)
            raw_list[[i]] <- data.table::as.data.table(level_raw)
        }
        
    } else {
        ## Ungrouped version
        total_tab <- table(var_vec, useNA = use_na)
        total_denom <- if (na_percent) sum(total_tab) else sum(total_tab[!is.na(names(total_tab))])
        
        formatted_list <- vector("list", n_levels)
        raw_list <- vector("list", n_levels)
        
        for (i in seq_along(levels_to_show)) {
            lvl <- levels_to_show[i]
            lvl_label <- if (is.na(lvl)) na_label else as.character(lvl)
            
            level_formatted <- list(
                variable = if (i == 1L) var_label else "",
                level = lvl_label
            )
            
            level_raw <- list(
                variable = if (i == 1L) var_label else "",
                level = lvl_label,
                stat_type = "category"
            )
            
            if (!isFALSE(total)) {
                if (is.na(lvl)) {
                    n <- sum(is.na(var_vec))
                    na_stat <- if (na_percent) stats else "n"
                } else {
                    n <- total_tab[as.character(lvl)]
                    if (is.na(n)) n <- 0L
                    na_stat <- stats
                }
                
                level_formatted[[total_label]] <- format_categorical_stat(
                    as.integer(n), as.integer(total_denom), na_stat, marks)
                level_raw[[total_label]] <- as.integer(n)
                level_raw[[paste0(total_label, "_total")]] <- as.integer(total_denom)
            }
            
            formatted_list[[i]] <- data.table::as.data.table(level_formatted)
            raw_list[[i]] <- data.table::as.data.table(level_raw)
        }
    }
    
    list(
        formatted = data.table::rbindlist(formatted_list, fill = TRUE),
        raw = data.table::rbindlist(raw_list, fill = TRUE)
    )
}


#' Perform statistical tests for categorical variables
#' 
#' Conducts chi-square or Fisher's exact tests for categorical variables
#' across groups. Automatically selects Fisher's exact test for small
#' expected frequencies. The test compares the observed levels of the variable
#' across the defined groups. The missing-value row and column that
#' \code{table(useNA = "ifany")} adds are removed, as are levels and groups
#' ithout observations: an empty row or column contributes nothing to either
#' test, but its zero expected counts would otherwise select the exact test
#' under \code{"auto"} and make the chi-square statistic undefined.
#'
#' @param tab Contingency table (matrix or table object).
#' @param test_type Character string: "chisq" for chi-square, "fisher" for
#'   Fisher's exact, or "auto" for automatic selection.
#' @param test_args Named list of additional test arguments from \code{...};
#'   each test receives only the arguments it accepts (see
#'   \code{filter_test_args()}).
#' @param var_name Character string naming the variable, used in warnings.
#' @return Numeric \emph{p}-value from the hypothesis test.
#' @keywords internal
perform_categorical_test <- function(tab, test_type, test_args = list(),
                                     var_name = NULL) {
    ## Keep the observed levels of the variable and the defined groups only
    tab_test <- tab[!is.na(rownames(tab)), !is.na(colnames(tab)), drop = FALSE]
    tab_test <- tab_test[rowSums(tab_test) > 0, colSums(tab_test) > 0, drop = FALSE]
    
    if (any(dim(tab_test) < 2)) return(NA_real_)
    
    ## Auto-select test. When the chi-square test is selected, the fit used for
    ## the selection also supplies its p-value; no warning is suppressed then,
    ## because every expected count is at least 5.
    chisq_fit <- NULL
    if (test_type == "auto") {
        chisq_fit <- suppressWarnings(stats::chisq.test(tab_test))
        test_type <- if (any(chisq_fit$expected < 5)) "fisher" else "chisq"
    }
    
    args <- filter_test_args(test_args, test_type)

    tryCatch({
        switch(test_type,
            "fisher" = ,
            "fisher.test" = fisher_test_pvalue(tab_test, args, var_name),
            "chisq" = ,
            "chisq.test" = if (!is.null(chisq_fit) && length(args) == 0L) {
                               chisq_fit$p.value
                           } else {
                               do.call(stats::chisq.test, c(list(tab_test), args))$p.value
                           },
            NA_real_
        )
    }, error = function(e) NA_real_)
}


#' Fisher's exact test with a Monte Carlo fallback
#'
#' Computes the \emph{p}-value of Fisher's exact test. For tables larger than
#' 2 \eqn{\times}{x} 2, \code{stats::fisher.test()} uses a network
#' algorithm with a bounded workspace (\code{workspace = 200000} by default),
#' which is exhausted by many tables with a few hundred observations or more.
#' A larger workspace is not a practical remedy, because computation time
#' grows steeply and the algorithm still fails for many such tables. When the
#' exact computation fails, the \emph{p}-value is estimated by Monte Carlo
#' simulation conditional on the margins, as \code{fisher.test()} itself
#' recommends, and a warning reports the substitution.
#'
#' @param tab Contingency table without missing-value rows or columns.
#' @param args Named list of arguments for \code{fisher.test()}. With
#'   \code{simulate.p.value = TRUE}, the simulation is used directly. A
#'   supplied \code{B} sets the number of simulated tables for the fallback,
#'   which otherwise uses 10,000.
#' @param var_name Character string naming the variable, used in the warning.
#' @return Numeric \emph{p}-value.
#' @keywords internal
fisher_test_pvalue <- function(tab, args = list(), var_name = NULL) {
    args <- utils::modifyList(list(workspace = 2e5), args)
    tryCatch(
        do.call(stats::fisher.test, c(list(tab), args))$p.value,
        error = function(e) {
            if (isTRUE(args$simulate.p.value)) stop(e)
            reason <- sub("\\..*$", "",
                          strsplit(conditionMessage(e), "\n", fixed = TRUE)[[1L]][1L])
            args$simulate.p.value <- TRUE
            args$B <- args$B %||% 10000
            warning("Fisher's exact test",
                    if (!is.null(var_name)) paste0(" for '", var_name, "'"),
                    " could not be computed exactly (", reason, "); the p-value ",
                    "is a Monte Carlo estimate from ",
                    format(args$B, big.mark = ",", scientific = FALSE),
                    " simulated tables. Set a seed for a reproducible value, or ",
                    "pass a larger 'workspace' through '...' to attempt the ",
                    "exact computation.", call. = FALSE)
            do.call(stats::fisher.test, c(list(tab), args))$p.value
        })
}


### * Time-to-event variables

#' Process time-to-event variable
#' 
#' Calculates survival statistics including median survival times with
#' confidence intervals, with optional grouping and log-rank testing. Parses
#' Surv() expressions and uses \pkg{survival} package functions.
#' 
#' @param data Data.table containing the survival variables.
#' @param var Character string with Surv() expression (\emph{e.g.,} "Surv(time, status)").
#' @param var_label Character string label for display.
#' @param group_var Optional character string naming the grouping variable.
#' @param digits Integer number of decimal places.
#' @param conf_level Numeric confidence level for confidence intervals.
#' @param na_include Logical whether to include missing values.
#' @param na_label Character string label for missing values.
#' @param test Logical whether to perform log-rank test.
#' @param total Logical or character controlling total column display.
#' @param total_label Character string label for total column.
#' @param groups Optional vector of groups (see \code{process_variable()}).
#' @param group_index Optional list of row positions for each group (see
#'   \code{process_variable()}).
#' @param ... Additional arguments; \code{rho} and \code{timefix} are passed
#'   to \code{survival::survdiff()}.
#' @return List with 'formatted' and 'raw' data.table components.
#' @keywords internal
process_survival <- function(data, var, var_label, group_var, digits,
                             conf_level = 0.95, na_include, na_label, 
                             test, total, total_label, marks = NULL,
                             groups = NULL, group_index = NULL, ...) {
    
    ## Parse Surv() expression
    surv_match <- regexec("Surv\\(([^,]+),\\s*([^)]+)\\)", var)
    surv_parts <- regmatches(var, surv_match)[[1]]
    
    if (length(surv_parts) < 3) {
        stop("Invalid Surv() syntax: ", var)
    }
    
    time_var <- trimws(surv_parts[2])
    status_var <- trimws(surv_parts[3])
    
    if (!requireNamespace("survival", quietly = TRUE)) {
        stop("Package 'survival' required for survival analysis")
    }
    
    ## Build CI label based on conf_level
    ci_pct <- round(conf_level * 100)
    ci_label <- paste0("Median (", ci_pct, "% CI)")
    
    ## Build CI column name suffixes for extracting from survfit summary.
    ## survfit names the CI columns paste(conf.int, "LCL", sep = ""), e.g.
    ## "0.95LCL" or "0.9LCL", so the same coercion is used here; formatting
    ## with nsmall = 2 gave "0.90LCL", which matches no column.
    ci_lower_name <- paste0(conf_level, "LCL")
    ci_upper_name <- paste0(conf_level, "UCL")
    
    fmt_str <- paste0("%.", digits, "f")
    
    if (!is.null(group_var)) {
        grp_vec <- data[[group_var]]
        if (is.null(groups)) groups <- derive_group_levels(grp_vec)
        if (is.null(group_index)) group_index <- derive_group_index(grp_vec, groups)
        
        ## Calculate p-value (log-rank test, or another member of the
        ## G-rho family when rho is supplied through ...)
        p_value <- if (test) {
            args <- filter_test_args(list(...), "survival")
            tryCatch({
                surv_obj <- survival::Surv(data[[time_var]], data[[status_var]])
                do.call(survival::survdiff, c(list(surv_obj ~ grp_vec), args))$pvalue
            }, error = function(e) NA_real_)
        } else {
            NULL
        }
        
        formatted_row <- list(
            variable = var_label,
            level = ci_label
        )
        
        raw_row <- list(
            variable = var_label,
            level = "median",
            stat_type = "survival"
        )
        
        ## Add total column
        if (!isFALSE(total)) {
            surv_obj <- survival::Surv(data[[time_var]], data[[status_var]])
            fit <- survival::survfit(surv_obj ~ 1, conf.int = conf_level)
            table <- summary(fit)$table
            
            formatted_row[[total_label]] <- format_survival_ci(
                table["median"], table[ci_lower_name], table[ci_upper_name],
                fmt_str, marks
            )
            
            raw_row[[total_label]] <- table["median"]
            raw_row[[paste0(total_label, "_ci_lower")]] <- table[ci_lower_name]
            raw_row[[paste0(total_label, "_ci_upper")]] <- table[ci_upper_name]
        }
        
        ## Add group columns
        for (j in seq_along(groups)) {
            grp_col <- as.character(groups[j])
            grp_idx <- group_index[[j]]
            grp_time <- data[[time_var]][grp_idx]
            grp_status <- data[[status_var]][grp_idx]
            
            surv_obj <- survival::Surv(grp_time, grp_status)
            fit <- survival::survfit(surv_obj ~ 1, conf.int = conf_level)
            table <- summary(fit)$table
            
            formatted_row[[grp_col]] <- format_survival_ci(
                table["median"], table[ci_lower_name], table[ci_upper_name],
                fmt_str, marks
            )
            
            raw_row[[grp_col]] <- table["median"]
            raw_row[[paste0(grp_col, "_ci_lower")]] <- table[ci_lower_name]
            raw_row[[paste0(grp_col, "_ci_upper")]] <- table[ci_upper_name]
        }
        
        if (test && !is.null(p_value)) {
            formatted_row[["p_value"]] <- p_value
            raw_row[["p_value"]] <- p_value
        }
        
        formatted_result <- data.table::as.data.table(formatted_row)
        raw_result <- data.table::as.data.table(raw_row)
        
    } else {
        formatted_row <- list(
            variable = var_label,
            level = ci_label
        )
        
        raw_row <- list(
            variable = var_label,
            level = "median",
            stat_type = "survival"
        )
        
        if (!isFALSE(total)) {
            surv_obj <- survival::Surv(data[[time_var]], data[[status_var]])
            fit <- survival::survfit(surv_obj ~ 1, conf.int = conf_level)
            table <- summary(fit)$table
            
            formatted_row[[total_label]] <- format_survival_ci(
                table["median"], table[ci_lower_name], table[ci_upper_name],
                fmt_str, marks
            )
            
            raw_row[[total_label]] <- table["median"]
            raw_row[[paste0(total_label, "_ci_lower")]] <- table[ci_lower_name]
            raw_row[[paste0(total_label, "_ci_upper")]] <- table[ci_upper_name]
        }
        
        formatted_result <- data.table::as.data.table(formatted_row)
        raw_result <- data.table::as.data.table(raw_row)
    }
    
    list(formatted = formatted_result, raw = raw_result)
}


### * Table formatting

#' Format \emph{p}-values for descriptive tables
#' 
#' Converts numeric p-values to formatted strings with appropriate precision.
#' Handles very small p-values with threshold notation (\emph{e.g.,} "< 0.001").
#' 
#' @param result Data.table with 'p_value' column to format.
#' @param p_digits Integer number of decimal places for p-values.
#' @param marks List with \code{big.mark} and \code{decimal.mark} as returned
#'   by \code{\link{resolve_number_marks}}.
#' @return Modified data.table with 'p_value' column (formatted strings).
#' @keywords internal
format_pvalues_desctable <- function(result, p_digits, marks) {
    if ("p_value" %chin% names(result)) {
        threshold <- 10^(-p_digits)
        threshold_str <- paste0("< 0", marks$decimal.mark,
                                strrep("0", p_digits - 1), "1")
        fmt_str <- paste0("%.", p_digits, "f")
        
        ## Format p-values vectorized
        p_str <- sprintf(fmt_str, result$p_value)
        if (marks$decimal.mark != ".") {
            p_str <- sub(".", marks$decimal.mark, p_str, fixed = TRUE)
        }
        
        result[, p_formatted := data.table::fifelse(
            is.na(p_value) | Variable == "N", 
            "",
            data.table::fifelse(
                p_value < threshold,
                threshold_str,
                p_str
            )
        )]
        
        ## Remove old column, rename new
        result[, p_value := NULL]
        data.table::setnames(result, "p_formatted", "p-value")
    }
    result
}


#' Reorder columns to position total column
#' 
#' Rearranges data.table columns to place the total column in the specified
#' position (first, last, or default). Ensures proper ordering of Variable,
#' Group, total, group columns, and \emph{p}-value.
#' 
#' @param result Data.table with columns to reorder.
#' @param total Logical or character: \code{TRUE}/"first" (total first), "last" (total last).
#' @param total_label Character string name of the total column.
#' @return Modified data.table with reordered columns.
#' @keywords internal
reorder_total_column <- function(result, total, total_label) {
    if (total_label %chin% names(result)) {
        cols <- names(result)
        
        base_cols <- c("Variable", "Group")
        group_cols <- cols[!cols %chin% c(base_cols, total_label, "p-value")]
        p_col <- if ("p-value" %chin% cols) "p-value" else NULL
        
        if (isTRUE(total) || identical(total, "first")) {
            new_order <- c(base_cols, total_label, group_cols, p_col)
        } else if (identical(total, "last")) {
            new_order <- c(base_cols, group_cols, total_label, p_col)
        } else {
            new_order <- c(base_cols, total_label, group_cols, p_col)
        }
        
        new_order <- new_order[new_order %chin% cols]
        data.table::setcolorder(result, new_order)
    }
    result
}
