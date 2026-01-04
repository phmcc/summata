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
#' @param na_include Logical whether to include missing values as a category.
#' @param na_label Character string label for missing values.
#' @param test Logical whether to perform statistical tests.
#' @param test_continuous Character string specifying test type for continuous variables.
#' @param test_categorical Character string specifying test type for categorical variables.
#' @param total Logical or character controlling total column display.
#' @param total_label Character string label for total column.
#' @param labels Named character vector of variable labels.
#' @param na_percent Logical whether to include NA in percentage denominators.
#' @param p_per_stat Logical whether to show separate p-values per statistic for
#'   continuous variables. Default FALSE for better performance.
#' @param ... Additional arguments passed to test functions.
#' @return List with 'formatted' and 'raw' data.table components.
#' @keywords internal
process_variable <- function(data, var, group_var = NULL, 
                             stats_continuous, stats_categorical,
                             digits, na_include, na_label,
                             test, test_continuous, test_categorical,
                             total, total_label, labels, na_percent, 
                             p_per_stat = FALSE, ...) {
    
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
            na_include = na_include,
            na_label = na_label,
            test = test,
            total = total,
            total_label = total_label,
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
            ...
        ))
    }
}

#' Process continuous variable
#' 
#' @param p_per_stat Logical. If TRUE, calculate separate p-values for each 
#'   statistic type (e.g., t-test for means, Wilcoxon for medians). If FALSE
#'   (default), calculate a single p-value based on the first statistic type.
#' @keywords internal
process_continuous <- function(data, var, var_label, group_var, stats, digits,
                               na_include, na_label, test, test_type,
                               total, total_label, p_per_stat = FALSE, ...) {
    
    ## Pre-extract the variable column
    var_vec <- data[[var]]
    not_na <- !is.na(var_vec)
    total_vals <- var_vec[not_na]
    total_n <- length(total_vals)
    
    ## Calculate total statistics
    if (total_n > 0) {
        total_fivenum <- fivenum(total_vals)
        total_stats <- list(
            mean = mean(total_vals),
            sd = sd(total_vals),
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
    groups <- NULL
    
    if (!is.null(group_var)) {
        grp_vec <- data[[group_var]]
        if (is.factor(grp_vec)) {
            groups <- levels(grp_vec)
        } else {
            groups <- unique(grp_vec)
            groups <- groups[!is.na(groups)]
        }
        
        ## Compute statistics for each group
        group_stats_list <- lapply(groups, function(g) {
            idx <- which(grp_vec == g & not_na)
            if (length(idx) == 0) {
                return(list(mean = NA_real_, sd = NA_real_, median = NA_real_,
                           q1 = NA_real_, q3 = NA_real_, min = NA_real_, 
                           max = NA_real_, n = 0L))
            }
            vals <- var_vec[idx]
            fn <- fivenum(vals)
            list(
                mean = mean(vals),
                sd = sd(vals),
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
        
        if (p_per_stat) {
            ## Calculate separate p-value for each statistic type
            for (stat_type in stats) {
                if (stat_type == "range") {
                    p_values[[stat_type]] <- NULL
                } else {
                    p_values[[stat_type]] <- perform_continuous_test(
                        var_vec, grp_vec, test_type, stat_type)
                }
            }
        } else {
            ## Calculate single p-value based on first non-range statistic
            first_stat <- stats[stats != "range"][1]
            if (is.na(first_stat)) first_stat <- "median_iqr"
            p_value_single <- perform_continuous_test(
                var_vec, grp_vec, test_type, first_stat)
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
                total_stats, stat_type, fmt_str)
            raw_row <- add_raw_stats(raw_row, total_label, total_stats, stat_type)
        }
        
        ## Add group columns
        if (!is.null(group_var)) {
            for (g in groups) {
                grp_col <- as.character(g)
                grp_stats <- group_stats_list[[grp_col]]
                formatted_row[[grp_col]] <- format_continuous_stat(
                    grp_stats, stat_type, fmt_str)
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
            miss_formatted[[total_label]] <- format_count(n_miss_total)
            miss_raw[[total_label]] <- n_miss_total
        }
        
        if (!is.null(group_var)) {
            grp_vec <- data[[group_var]]
            for (g in groups) {
                grp_col <- as.character(g)
                n_miss <- sum(!not_na & grp_vec == g, na.rm = TRUE)
                miss_formatted[[grp_col]] <- format_count(n_miss)
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

#' Perform continuous tests
#' @keywords internal
perform_continuous_test <- function(var_vec, grp_vec, test_type, stat_type) {
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
    
    tryCatch({
        switch(test_type,
            "t" = t.test(x ~ g)$p.value,
            "wrs" = wilcox.test(x ~ g)$p.value,
            "aov" = {
                fit <- aov(x ~ g)
                summary(fit)[[1]][["Pr(>F)"]][1]
            },
            "kwt" = kruskal.test(x ~ g)$p.value,
            NA_real_
        )
    }, error = function(e) NA_real_)
}


#' Format continuous statistic
#' @keywords internal
format_continuous_stat <- function(stats, stat_type, fmt_str) {
    if (stats$n == 0) return("")
    
    switch(stat_type,
        "mean_sd" = paste0(
            format_num(stats$mean, fmt_str), " +/- ", 
            format_num(stats$sd, fmt_str)
        ),
        "median_iqr" = paste0(
            format_num(stats$median, fmt_str), " [",
            format_num(stats$q1, fmt_str), "-",
            format_num(stats$q3, fmt_str), "]"
        ),
        "median_range" = {
            sep <- if (stats$min < 0 || stats$max < 0) ", " else "-"
            paste0(
                format_num(stats$median, fmt_str), " (",
                format_num(stats$min, fmt_str), sep,
                format_num(stats$max, fmt_str), ")"
            )
        },
        "range" = {
            sep <- if (stats$min < 0 || stats$max < 0) " to " else "-"
            paste0(format_num(stats$min, fmt_str), sep, 
                   format_num(stats$max, fmt_str))
        },
        ""
    )
}


#' Format numbers with commas for large values
#' @keywords internal
format_num <- function(x, fmt_str) {
    if (is.na(x)) return("")
    if (abs(x) >= 1000) {
        format(round(x, 1), big.mark = ",")
    } else {
        result <- sprintf(fmt_str, x)
        ## Fix negative zero: "-0.0", "-0.00", etc. -> "0.0", "0.00"
        gsub("(?<![0-9])-0(\\.0+)(?![0-9])", "0\\1", result, perl = TRUE)
    }
}


#' Format count
#' @keywords internal
format_count <- function(n) {
    if (n >= 1000) format(n, big.mark = ",") else as.character(n)
}


#' Get statistic label
#' @keywords internal
get_stat_label <- function(stat_type) {
    switch(stat_type,
           "mean_sd" = "Mean +/- SD",
           "median_iqr" = "Median [IQR]",
           "median_range" = "Median (Range)",
           "range" = "Range",
           "n_miss" = "Missing",
           stat_type
    )
}


#' Process categorical variable
#' @keywords internal
process_categorical <- function(data, var, var_label, group_var, stats,
                                na_include, na_label, test, test_type,
                                total, total_label, na_percent, ...) {
    
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
        
        if (is.factor(grp_vec)) {
            groups <- levels(grp_vec)
        } else {
            groups <- unique(grp_vec)
            groups <- groups[!is.na(groups)]
        }
        
        ## Single table call for cross-tabulation
        tab <- table(var_vec, grp_vec, useNA = use_na)
        total_tab <- rowSums(tab)
        grp_totals <- colSums(tab)
        
        ## Calculate denominators
        if (na_percent) {
            total_denom <- sum(tab)
            grp_denoms <- grp_totals
        } else {
            ## Exclude NA row from denominators
            non_na_rows <- !is.na(rownames(tab))
            total_denom <- sum(tab[non_na_rows, , drop = FALSE])
            grp_denoms <- colSums(tab[non_na_rows, , drop = FALSE])
        }
        
        ## Calculate p-value once
        p_value <- if (test) {
            perform_categorical_test(tab, test_type)
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
                tab_row <- tab[is.na(rownames(tab)), , drop = FALSE]
                n_total <- sum(tab_row)
            } else {
                n_total <- total_tab[lvl_char]
                if (is.na(n_total)) n_total <- 0L
            }
            
            ## Add total column
            if (!isFALSE(total)) {
                denom <- if (is.na(lvl)) sum(tab) else total_denom
                level_formatted[[total_label]] <- format_categorical_stat(
                    as.integer(n_total), as.integer(denom), stats)
                level_raw[[total_label]] <- as.integer(n_total)
                level_raw[[paste0(total_label, "_total")]] <- as.integer(denom)
            }
            
            ## Add group columns
            for (g in groups) {
                grp_col <- as.character(g)
                
                if (is.na(lvl)) {
                    n <- sum(tab[is.na(rownames(tab)), grp_col])
                    denom <- grp_totals[grp_col]
                } else {
                    n <- tab[lvl_char, grp_col]
                    if (is.na(n)) n <- 0L
                    denom <- grp_denoms[grp_col]
                }
                
                level_formatted[[grp_col]] <- format_categorical_stat(
                    as.integer(n), as.integer(denom), stats)
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
                } else {
                    n <- total_tab[as.character(lvl)]
                    if (is.na(n)) n <- 0L
                }
                
                level_formatted[[total_label]] <- format_categorical_stat(
                    as.integer(n), as.integer(total_denom), stats)
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


#' Perform categorical test
#' @keywords internal
perform_categorical_test <- function(tab, test_type) {
    ## Remove NA rows for testing
    non_na_rows <- !is.na(rownames(tab))
    tab_test <- tab[non_na_rows, , drop = FALSE]
    
    if (any(dim(tab_test) < 2)) return(NA_real_)
    
    ## Auto-select test
    if (test_type == "auto") {
        expected <- suppressWarnings(chisq.test(tab_test)$expected)
        test_type <- if (any(expected < 5)) "fisher" else "chisq"
    }
    
    tryCatch({
        switch(test_type,
            "fisher" = ,
            "fisher.test" = fisher.test(tab_test, workspace = 2e5)$p.value,
            "chisq" = ,
            "chisq.test" = chisq.test(tab_test)$p.value,
            NA_real_
        )
    }, error = function(e) NA_real_)
}


#' Format categorical statistic
#' @keywords internal
format_categorical_stat <- function(n, total, stat_type) {
    if (length(stat_type) > 1) stat_type <- stat_type[1]
    
    n_formatted <- if (n >= 1000) format(n, big.mark = ",") else as.character(n)
    
    switch(stat_type,
           "n" = n_formatted,
           "percent" = sprintf("%.1f%%", 100 * n / total),
           "n_percent" = sprintf("%s (%.1f%%)", n_formatted, 100 * n / total),
           n_formatted
    )
}


#' Process survival variable
#' @keywords internal
process_survival <- function(data, var, var_label, group_var, digits,
                             na_include, na_label, test, total, total_label, ...) {
    
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
    
    fmt_str <- paste0("%.", digits, "f")
    
    if (!is.null(group_var)) {
        grp_vec <- data[[group_var]]
        if (is.factor(grp_vec)) {
            groups <- levels(grp_vec)
        } else {
            groups <- unique(grp_vec)
            groups <- groups[!is.na(groups)]
        }
        
        ## Calculate p-value (log-rank test)
        p_value <- if (test) {
            tryCatch({
                surv_obj <- survival::Surv(data[[time_var]], data[[status_var]])
                survival::survdiff(surv_obj ~ grp_vec)$pvalue
            }, error = function(e) NA_real_)
        } else {
            NULL
        }
        
        formatted_row <- list(
            variable = var_label,
            level = "Median (95% CI)"
        )
        
        raw_row <- list(
            variable = var_label,
            level = "median",
            stat_type = "survival"
        )
        
        ## Add total column
        if (!isFALSE(total)) {
            surv_obj <- survival::Surv(data[[time_var]], data[[status_var]])
            fit <- survival::survfit(surv_obj ~ 1)
            table <- summary(fit)$table
            
            formatted_row[[total_label]] <- sprintf(
                paste0(fmt_str, " (", fmt_str, "-", fmt_str, ")"),
                table["median"], table["0.95LCL"], table["0.95UCL"]
            )
            
            raw_row[[total_label]] <- table["median"]
            raw_row[[paste0(total_label, "_ci_lower")]] <- table["0.95LCL"]
            raw_row[[paste0(total_label, "_ci_upper")]] <- table["0.95UCL"]
        }
        
        ## Add group columns
        for (g in groups) {
            grp_col <- as.character(g)
            grp_idx <- which(grp_vec == g)
            grp_time <- data[[time_var]][grp_idx]
            grp_status <- data[[status_var]][grp_idx]
            
            surv_obj <- survival::Surv(grp_time, grp_status)
            fit <- survival::survfit(surv_obj ~ 1)
            table <- summary(fit)$table
            
            formatted_row[[grp_col]] <- sprintf(
                paste0(fmt_str, " (", fmt_str, "-", fmt_str, ")"),
                table["median"], table["0.95LCL"], table["0.95UCL"]
            )
            
            raw_row[[grp_col]] <- table["median"]
            raw_row[[paste0(grp_col, "_ci_lower")]] <- table["0.95LCL"]
            raw_row[[paste0(grp_col, "_ci_upper")]] <- table["0.95UCL"]
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
            level = "Median (95% CI)"
        )
        
        raw_row <- list(
            variable = var_label,
            level = "median",
            stat_type = "survival"
        )
        
        if (!isFALSE(total)) {
            surv_obj <- survival::Surv(data[[time_var]], data[[status_var]])
            fit <- survival::survfit(surv_obj ~ 1)
            table <- summary(fit)$table
            
            formatted_row[[total_label]] <- sprintf(
                paste0(fmt_str, " (", fmt_str, "-", fmt_str, ")"),
                table["median"], table["0.95LCL"], table["0.95UCL"]
            )
            
            raw_row[[total_label]] <- table["median"]
            raw_row[[paste0(total_label, "_ci_lower")]] <- table["0.95LCL"]
            raw_row[[paste0(total_label, "_ci_upper")]] <- table["0.95UCL"]
        }
        
        formatted_result <- data.table::as.data.table(formatted_row)
        raw_result <- data.table::as.data.table(raw_row)
    }
    
    list(formatted = formatted_result, raw = raw_result)
}


#' Format p-values in result table
#' @keywords internal
format_pvalues_desctable <- function(result, p_digits) {
    if ("p_value" %chin% names(result)) {
        threshold <- 10^(-p_digits)
        threshold_str <- paste0("< 0.", strrep("0", p_digits - 1), "1")
        fmt_str <- paste0("%.", p_digits, "f")
        
        result[, p_formatted := data.table::fifelse(
            is.na(p_value) | Variable == "N", 
            "",
            data.table::fifelse(
                p_value < threshold,
                threshold_str,
                sprintf(fmt_str, p_value)
            )
        )]
        
        ## Remove old column, rename new
        result[, p_value := NULL]
        data.table::setnames(result, "p_formatted", "p-value")
    }
    result
}


#' Reorder columns to position total column
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
