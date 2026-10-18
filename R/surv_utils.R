### * Data extraction

#' Process survival probability quantiles
#'
#' Extracts survival time quantiles from survfit objects.
#'
#' @param survfit_objects List of survfit objects.
#' @param probs Numeric vector of probabilities.
#' @param groups Character vector of group names.
#' @param group_labels Character vector of group display labels.
#' @param time_digits Integer decimal places for time values.
#' @param total Logical or character controlling total column.
#' @param total_label Character label for total column.
#' @param median_label Character label for median row.
#' @param by Character name of stratifying variable.
#' @param data Data.table with the source data.
#' @return List with formatted and raw data.tables.
#' @keywords internal
process_survival_probs <- function(survfit_objects,
                                   probs,
                                   groups,
                                   group_labels,
                                   time_digits,
                                   total,
                                   total_label,
                                   median_label,
                                   by,
                                   data,
                                   conf_level = 0.95,
                                   marks = NULL) {

    ## Build column headers for quantiles
    ci_pct <- round(conf_level * 100)
    quantile_cols <- ifelse(
        probs == 0.5,
        median_label,
        sprintf("%d%% Survival Time (%d%% CI)", round((1 - probs) * 100), ci_pct)
    )

    n_probs <- length(probs)
    fmt_str <- paste0("%.", time_digits, "f")

    ## Extract all quantile data at once
    if (!is.null(by) && !is.null(groups) && "stratified" %in% names(survfit_objects)) {
        fit <- survfit_objects$stratified
        quant <- stats::quantile(fit, probs = probs)

        n_rows <- length(groups)
        row_labels <- group_labels

        ## quant$quantile is a matrix of strata by probabilities. Its rows are
        ## matched to the groups by stratum name, as in
        ## process_survival_times(): survfit() orders the strata itself and
        ## omits groups without observations, which keep missing values
        strata_rows <- match(paste0(by, "=", groups), names(fit$strata))
        est_mat <- matrix(quant$quantile, ncol = n_probs)[strata_rows, , drop = FALSE]
        lower_mat <- matrix(quant$lower, ncol = n_probs)[strata_rows, , drop = FALSE]
        upper_mat <- matrix(quant$upper, ncol = n_probs)[strata_rows, , drop = FALSE]
        no_stratum <- is.na(strata_rows)
    } else {
        ## Overall only
        fit <- survfit_objects$overall
        quant <- stats::quantile(fit, probs = probs)

        n_rows <- 1L
        row_labels <- total_label

        est_mat <- matrix(quant$quantile, nrow = 1)
        lower_mat <- matrix(quant$lower, nrow = 1)
        upper_mat <- matrix(quant$upper, nrow = 1)
        no_stratum <- FALSE
    }

    ## Build formatted output
    formatted_list <- vector("list", n_probs + 1L)
    raw_list <- vector("list", n_probs * 3L + 1L)

    formatted_list[[1L]] <- row_labels
    raw_list[[1L]] <- row_labels
    names(formatted_list)[1L] <- "Group"
    names(raw_list)[1L] <- "Group"

    ## Formatting for each quantile column
    for (j in seq_len(n_probs)) {
        col_name <- quantile_cols[j]

        est_vec <- est_mat[, j]
        lower_vec <- lower_mat[, j]
        upper_vec <- upper_mat[, j]

        ## Cell formatting. A missing quantile is one the survival curve
        ## did not reach, except in groups without observations, which have
        ## no curve
        cells <- format_quantile_cells(
            est_vec, lower_vec, upper_vec, fmt_str, marks
        )
        cells[no_stratum] <- "-"
        formatted_list[[j + 1L]] <- cells
        names(formatted_list)[j + 1L] <- col_name

        ## Raw values
        raw_idx <- (j - 1L) * 3L + 2L
        raw_list[[raw_idx]] <- est_vec
        raw_list[[raw_idx + 1L]] <- lower_vec
        raw_list[[raw_idx + 2L]] <- upper_vec
        names(raw_list)[raw_idx:(raw_idx + 2L)] <- paste0(col_name, c("_estimate", "_lower", "_upper"))
    }

    formatted <- data.table::as.data.table(formatted_list)
    raw <- data.table::as.data.table(raw_list)

    ## Add overall row if stratified and total requested
    if (!is.null(by) && !isFALSE(total) && "overall" %in% names(survfit_objects)) {
        fit_overall <- survfit_objects$overall
        quant_overall <- stats::quantile(fit_overall, probs = probs)

        est_overall <- as.numeric(quant_overall$quantile)
        lower_overall <- as.numeric(quant_overall$lower)
        upper_overall <- as.numeric(quant_overall$upper)

        ## Build overall row
        overall_formatted <- list(Group = total_label)
        overall_raw <- list(Group = total_label)

        for (j in seq_len(n_probs)) {
            col_name <- quantile_cols[j]
            overall_formatted[[col_name]] <- format_quantile_cells(
                est_overall[j], lower_overall[j], upper_overall[j], fmt_str, marks
            )

            overall_raw[[paste0(col_name, "_estimate")]] <- est_overall[j]
            overall_raw[[paste0(col_name, "_lower")]] <- lower_overall[j]
            overall_raw[[paste0(col_name, "_upper")]] <- upper_overall[j]
        }

        overall_formatted <- data.table::as.data.table(overall_formatted)
        overall_raw <- data.table::as.data.table(overall_raw)

        ## Position total row
        if (isTRUE(total) || identical(total, "first")) {
            formatted <- rbind(overall_formatted, formatted, fill = TRUE)
            raw <- rbind(overall_raw, raw, fill = TRUE)
        } else if (identical(total, "last")) {
            formatted <- rbind(formatted, overall_formatted, fill = TRUE)
            raw <- rbind(raw, overall_raw, fill = TRUE)
        }
    }

    list(formatted = formatted, raw = raw)
}


### * Statistical calculations

#' Perform survival comparison test
#'
#' Performs statistical test comparing survival curves across groups. The
#' log-rank test is computed by \code{survival::survdiff()}. The Gehan-Breslow
#' (\code{"wilcoxon"}), Tarone-Ware (\code{"tarone"}), and Peto-Peto
#' (\code{"petopeto"}) tests weight each event time by the number at risk, its
#' square root, and the Peto-Peto estimate of the pooled survival function,
#' respectively (see \code{weighted_logrank_test()}). \code{survdiff()} offers
#' only the weights \eqn{\hat{S}(t-)^\rho}{S(t-)^rho} of the Fleming-Harrington
#' tests, of which \eqn{\rho = 1}{rho = 1} approximates the Peto-Peto test and
#' none gives the other two. Groups without expected events, such as a group
#' whose observations are all censored before the first event, are left out of
#' the test, as \code{survdiff()} leaves them out, and the degrees of freedom
#' count the groups that remain.
#'
#' @param surv_obj Survival object created by Surv().
#' @param group_var Vector of group assignments.
#' @param test_type Character string specifying test type.
#' @return List with test statistic, degrees of freedom, \emph{p}-value, and
#'   test type. The statistic, degrees of freedom, and \emph{p}-value are
#'   \code{NA} when fewer than two groups have expected events.
#' @keywords internal
perform_survival_test <- function(surv_obj, group_var, test_type) {
    tryCatch({
        if (identical(test_type, "logrank")) {
            test_obj <- survival::survdiff(surv_obj ~ group_var)
            statistic <- test_obj$chisq
            ## Groups without expected events are left out of the test, so
            ## the degrees of freedom are counted as survdiff() counts them
            df <- sum(test_obj$exp > 0) - 1L
            p_value <- test_obj$pvalue
        } else {
            weights <- switch(test_type,
                              "wilcoxon" = function(n_risk, n_event) n_risk,
                              "tarone" = function(n_risk, n_event) sqrt(n_risk),
                              "petopeto" = function(n_risk, n_event) {
                                  cumprod(1 - n_event / (n_risk + 1))
                              },
                              stop("unknown test type: ", test_type))
            test_obj <- weighted_logrank_test(surv_obj, group_var, weights)
            statistic <- test_obj$chisq
            df <- test_obj$df
            p_value <- test_obj$pvalue
        }

        ## Fewer than two groups with expected events leave nothing to test
        if (df < 1L) {
            statistic <- NA_real_
            df <- NA_integer_
            p_value <- NA_real_
        }
        list(
            statistic = statistic,
            df = df,
            p_value = p_value,
            test_type = test_type,
            test_object = test_obj
        )
    }, error = function(e) {
        list(
            statistic = NA_real_,
            df = NA_integer_,
            p_value = NA_real_,
            test_type = test_type,
            test_object = NULL,
            error = e$message
        )
    })
}


#' Weighted log-rank test
#'
#' Compares survival curves across groups with a weighted log-rank statistic.
#' At each distinct event time \eqn{t_j}, with \eqn{n_j} subjects at risk and
#' \eqn{d_j} events in the pooled sample, the observed and expected numbers of
#' events in each group are compared with weight, with hypergeometric
#' covariance. The statistic \eqn{v' V^{-1} v}, computed after removing one
#' group, has a \eqn{\chi^2} distribution with one degree of freedom fewer than
#' the number of groups. These are the definitions of the SAS LIFETEST and Stata
#' \code{sts test} procedures; with unit weights the statistic equals that of
#' \code{survival::survdiff()}, and times are aligned with
#' \code{survival::aeqSurv()} as \code{survdiff()} aligns them. Groups with no
#' subjects at risk at any event time are excluded, as in \code{survdiff()};
#' when fewer than two groups remain, the statistic and \emph{p}-value are
#' \code{NA}.
#'
#' @param surv_obj Survival object created by \code{survival::Surv()}, with
#'   right-censored data.
#' @param group_var Vector of group assignments of the same length.
#' @param weights Function of the numbers at risk and of events at the
#'   ordered distinct event times of the pooled sample, returning the weight
#'   of each event time.
#' @return List with the components of a \code{survdiff} object: the number of
#'   subjects in each group (\code{n}), the weighted observed and expected
#'   numbers of events in each group (\code{obs}, \code{exp}), their
#'   covariance matrix (\code{var}), the test statistic (\code{chisq}), and
#'   the \emph{p}-value (\code{pvalue}), with its degrees of freedom
#'   (\code{df}).
#' @keywords internal
weighted_logrank_test <- function(surv_obj, group_var, weights) {
    surv_obj <- survival::aeqSurv(surv_obj)
    time <- surv_obj[, "time"]
    status <- surv_obj[, "status"]
    keep <- !is.na(time) & !is.na(status) & !is.na(group_var)
    time <- time[keep]
    status <- status[keep]
    group <- droplevels(as.factor(group_var[keep]))
    n_groups <- nlevels(group)
    if (n_groups < 2L) {
        stop("There is only 1 group")
    }
    group_sizes <- table(group, dnn = NULL)
    group <- as.integer(group)

    ## Numbers at risk and events in each group at the distinct event times
    event_times <- sort(unique(time[status == 1]))
    n_times <- length(event_times)
    at_risk <- matrix(0, n_times, n_groups)
    events <- matrix(0, n_times, n_groups)
    for (k in seq_len(n_groups)) {
        group_times <- sort(time[group == k])
        at_risk[, k] <- length(group_times) -
            findInterval(event_times, group_times, left.open = TRUE)
        events[, k] <- tabulate(match(time[group == k & status == 1], event_times),
                                nbins = n_times)
    }
    n_risk <- rowSums(at_risk)
    n_event <- rowSums(events)
    w <- weights(n_risk, n_event)

    share <- at_risk / n_risk
    observed <- colSums(w * events)
    expected <- colSums(w * n_event * share)
    ## With a single subject at risk the hypergeometric variance is zero
    var_weight <- ifelse(n_risk > 1, w^2 * n_event * (n_risk - n_event) / (n_risk - 1), 0)
    v <- diag(colSums(var_weight * share), n_groups) - crossprod(share, var_weight * share)

    tested <- expected > 0
    df <- sum(tested) - 1L
    if (df < 1L) {
        chisq <- NA_real_
        p_value <- NA_real_
    } else {
        u <- (observed - expected)[tested][-1L]
        v_tested <- v[tested, tested, drop = FALSE][-1L, -1L, drop = FALSE]
        chisq <- sum(solve(v_tested, u) * u)
        p_value <- stats::pchisq(chisq, df, lower.tail = FALSE)
    }

    list(n = group_sizes,
         obs = observed,
         exp = expected,
         var = v,
         chisq = chisq,
         pvalue = p_value,
         df = df)
}


#' Process survival at specified time points
#'
#' Extracts survival probabilities at specified time points from survfit objects.
#'
#' @param survfit_objects List of survfit objects.
#' @param times Numeric vector of time points.
#' @param groups Character vector of group names.
#' @param group_labels Character vector of group display labels.
#' @param stats Character vector of statistics to include.
#' @param type Character string specifying probability type.
#' @param digits Integer decimal places for percentages.
#' @param percent Logical whether to display as percentages.
#' @param total Logical or character controlling total column.
#' @param total_label Character label for total column.
#' @param time_label Character template for time column headers.
#' @param time_unit Character time unit for column headers.
#' @param by Character name of stratifying variable.
#' @param data Data.table with the source data.
#' @return List with formatted and raw data.tables.
#' @keywords internal
process_survival_times <- function(survfit_objects,
                                   times,
                                   groups,
                                   group_labels,
                                   stats,
                                   type,
                                   digits,
                                   percent,
                                   total,
                                   total_label,
                                   time_label,
                                   time_unit,
                                   by,
                                   data,
                                   marks = NULL) {

    ## Build column headers for times
    time_cols <- vapply(times, function(t) {
        label <- gsub("\\{time\\}", t, time_label)
        if (!is.null(time_unit)) {
            label <- gsub("\\{unit\\}", time_unit, label)
        } else {
            label <- gsub("\\s*\\{unit\\}", "", label)
        }
        label
    }, character(1))

    n_times <- length(times)

    ## Pre-compute format strings
    if (percent) {
        fmt_est <- paste0("%.", digits, "f%%")
        fmt_ci_lower <- paste0("%.", digits, "f")
        fmt_ci_upper <- paste0("%.", digits, "f%%")
    } else {
        fmt_est <- paste0("%.", digits + 2L, "f")
        fmt_ci_lower <- paste0("%.", digits + 2L, "f")
        fmt_ci_upper <- paste0("%.", digits + 2L, "f")
    }

    ## Extract all survival data at once from stratified fit
    if (!is.null(by) && !is.null(groups) && "stratified" %in% names(survfit_objects)) {
        fit <- survfit_objects$stratified
        summ <- summary(fit, times = times, extend = TRUE)

        ## Build lookup for strata indices
        strata_vec <- as.character(summ$strata)
        strata_lookup <- split(seq_along(strata_vec), strata_vec)

        n_rows <- length(groups)
        row_labels <- group_labels

        ## Pre-allocate result matrices
        est_mat <- matrix(NA_real_, nrow = n_rows, ncol = n_times)
        lower_mat <- matrix(NA_real_, nrow = n_rows, ncol = n_times)
        upper_mat <- matrix(NA_real_, nrow = n_rows, ncol = n_times)
        n_risk_mat <- matrix(NA_integer_, nrow = n_rows, ncol = n_times)
        n_event_mat <- matrix(NA_integer_, nrow = n_rows, ncol = n_times)

        ## Extract data for each stratum
        for (i in seq_along(groups)) {
            strata_name <- paste0(by, "=", groups[i])
            idx <- strata_lookup[[strata_name]]

            if (!is.null(idx) && length(idx) > 0) {
                est_mat[i, ] <- summ$surv[idx]
                lower_mat[i, ] <- summ$lower[idx]
                upper_mat[i, ] <- summ$upper[idx]
                n_risk_mat[i, ] <- summ$n.risk[idx]
                n_event_mat[i, ] <- summ$n.event[idx]
            }
        }

        ## Apply type transformation
        if (type == "risk") {
            est_mat <- 1 - est_mat
            temp <- 1 - lower_mat
            lower_mat <- 1 - upper_mat
            upper_mat <- temp
        } else if (type == "cumhaz") {
            est_mat <- -log(est_mat)
            temp <- -log(lower_mat)
            lower_mat <- -log(upper_mat)
            upper_mat <- temp
        }

    } else {
        ## Overall only
        fit <- survfit_objects$overall
        summ <- summary(fit, times = times, extend = TRUE)

        n_rows <- 1L
        row_labels <- total_label

        est_mat <- matrix(summ$surv, nrow = 1)
        lower_mat <- matrix(summ$lower, nrow = 1)
        upper_mat <- matrix(summ$upper, nrow = 1)
        n_risk_mat <- matrix(summ$n.risk, nrow = 1)
        n_event_mat <- matrix(summ$n.event, nrow = 1)

        ## Apply type transformation
        if (type == "risk") {
            est_mat <- 1 - est_mat
            temp <- 1 - lower_mat
            lower_mat <- 1 - upper_mat
            upper_mat <- temp
        } else if (type == "cumhaz") {
            est_mat <- -log(est_mat)
            temp <- -log(lower_mat)
            lower_mat <- -log(upper_mat)
            upper_mat <- temp
        }
    }

    ## Build formatted output
    formatted_list <- vector("list", n_times + 1L)
    raw_list <- vector("list", n_times * 5L + 1L)

    formatted_list[[1L]] <- row_labels
    raw_list[[1L]] <- row_labels
    names(formatted_list)[1L] <- "Group"
    names(raw_list)[1L] <- "Group"

    ## Formatting for each time column
    for (j in seq_len(n_times)) {
        col_name <- time_cols[j]

        est_vec <- est_mat[, j]
        lower_vec <- lower_mat[, j]
        upper_vec <- upper_mat[, j]
        n_risk_vec <- n_risk_mat[, j]
        n_event_vec <- n_event_mat[, j]

        ## Cell formatting
        formatted_list[[j + 1L]] <- format_survival_cells(
            est_vec, lower_vec, upper_vec, n_risk_vec, n_event_vec,
            stats, fmt_est, fmt_ci_lower, fmt_ci_upper, percent, marks
        )
        names(formatted_list)[j + 1L] <- col_name

        ## Raw values
        raw_idx <- (j - 1L) * 5L + 2L
        raw_list[[raw_idx]] <- est_vec
        raw_list[[raw_idx + 1L]] <- lower_vec
        raw_list[[raw_idx + 2L]] <- upper_vec
        raw_list[[raw_idx + 3L]] <- n_risk_vec
        raw_list[[raw_idx + 4L]] <- n_event_vec
        names(raw_list)[raw_idx:(raw_idx + 4L)] <- paste0(col_name, c("_estimate", "_lower", "_upper", "_n_risk", "_n_event"))
    }

    formatted <- data.table::as.data.table(formatted_list)
    raw <- data.table::as.data.table(raw_list)

    ## Add overall row if stratified and total requested
    if (!is.null(by) && !isFALSE(total) && "overall" %in% names(survfit_objects)) {
        fit_overall <- survfit_objects$overall
        summ_overall <- summary(fit_overall, times = times, extend = TRUE)

        est_overall <- summ_overall$surv
        lower_overall <- summ_overall$lower
        upper_overall <- summ_overall$upper
        n_risk_overall <- summ_overall$n.risk
        n_event_overall <- summ_overall$n.event

        ## Apply type transformation
        if (type == "risk") {
            est_overall <- 1 - est_overall
            temp <- 1 - lower_overall
            lower_overall <- 1 - upper_overall
            upper_overall <- temp
        } else if (type == "cumhaz") {
            est_overall <- -log(est_overall)
            temp <- -log(lower_overall)
            lower_overall <- -log(upper_overall)
            upper_overall <- temp
        }

        ## Build overall row
        overall_formatted <- list(Group = total_label)
        overall_raw <- list(Group = total_label)

        for (j in seq_len(n_times)) {
            col_name <- time_cols[j]
            overall_formatted[[col_name]] <- format_survival_cells(
                est_overall[j], lower_overall[j], upper_overall[j],
                n_risk_overall[j], n_event_overall[j],
                stats, fmt_est, fmt_ci_lower, fmt_ci_upper, percent, marks
            )

            overall_raw[[paste0(col_name, "_estimate")]] <- est_overall[j]
            overall_raw[[paste0(col_name, "_lower")]] <- lower_overall[j]
            overall_raw[[paste0(col_name, "_upper")]] <- upper_overall[j]
            overall_raw[[paste0(col_name, "_n_risk")]] <- n_risk_overall[j]
            overall_raw[[paste0(col_name, "_n_event")]] <- n_event_overall[j]
        }

        overall_formatted <- data.table::as.data.table(overall_formatted)
        overall_raw <- data.table::as.data.table(overall_raw)

        ## Position total row
        if (isTRUE(total) || identical(total, "first")) {
            formatted <- rbind(overall_formatted, formatted, fill = TRUE)
            raw <- rbind(overall_raw, raw, fill = TRUE)
        } else if (identical(total, "last")) {
            formatted <- rbind(formatted, overall_formatted, fill = TRUE)
            raw <- rbind(raw, overall_raw, fill = TRUE)
        }
    }

    list(formatted = formatted, raw = raw)
}


### * Table formatting

#' Survival cell formatting
#'
#' Formats survival probability cells for multiple rows at once. Uses
#' locale-aware decimal marks and safe CI separators that avoid ambiguity
#' with negative values or decimal commas.
#'
#' @param est Numeric vector of estimates.
#' @param lower Numeric vector of lower CI bounds.
#' @param upper Numeric vector of upper CI bounds.
#' @param n_risk Integer vector of numbers at risk.
#' @param n_event Integer vector of event counts.
#' @param stats Character vector of statistics to include.
#' @param fmt_est Format string for estimate.
#' @param fmt_ci_lower Format string for lower CI bound.
#' @param fmt_ci_upper Format string for upper CI bound.
#' @param percent Logical whether percentages.
#' @param marks List with \code{big.mark} and \code{decimal.mark} as returned
#'   by \code{\link{resolve_number_marks}}.
#' @return Character vector of formatted cells.
#' @keywords internal
format_survival_cells <- function(est, lower, upper, n_risk, n_event,
                                             stats, fmt_est, fmt_ci_lower,
                                             fmt_ci_upper, percent,
                                             marks = NULL) {
    n <- length(est)
    result <- character(n)

    ## Handle NA values
    na_mask <- is.na(est)

    if (all(na_mask)) {
        return(rep("-", n))
    }

    ## Format non-NA values
    valid_idx <- which(!na_mask)

    if (percent) {
        est_vals <- est[valid_idx] * 100
        lower_vals <- lower[valid_idx] * 100
        upper_vals <- upper[valid_idx] * 100
    } else {
        est_vals <- est[valid_idx]
        lower_vals <- lower[valid_idx]
        upper_vals <- upper[valid_idx]
    }

    ## Format individual components
    est_fmt <- sprintf(fmt_est, est_vals)
    lower_fmt <- sprintf(fmt_ci_lower, lower_vals)
    upper_fmt <- sprintf(fmt_ci_upper, upper_vals)

    ## Apply locale decimal mark and fix negative zeros
    if (!is.null(marks)) {
        est_fmt <- apply_decimal_mark(est_fmt, marks)
        lower_fmt <- apply_decimal_mark(lower_fmt, marks)
        upper_fmt <- apply_decimal_mark(upper_fmt, marks)
    }

    ## Determine CI separator (consistent within column)
    any_negative <- any(lower_vals < 0, upper_vals < 0, na.rm = TRUE)
    if (any_negative) {
        sep <- " to "
    } else if (!is.null(marks) && marks$decimal.mark == ",") {
        sep <- "\u2013"
    } else {
        sep <- "-"
    }

    ## Build formatted strings
    if ("survival" %in% stats && "ci" %in% stats) {
        formatted <- paste0(est_fmt, " (", lower_fmt, sep, upper_fmt, ")")
    } else if ("survival" %in% stats) {
        formatted <- est_fmt
    } else if ("ci" %in% stats) {
        formatted <- paste0("(", lower_fmt, sep, upper_fmt, ")")
    } else {
        formatted <- est_fmt
    }

    ## Add n_risk if requested
    if ("n_risk" %in% stats) {
        n_risk_fmt <- if (!is.null(marks)) {
            vapply(n_risk[valid_idx], format_count, character(1), marks = marks)
        } else {
            as.character(n_risk[valid_idx])
        }
        formatted <- paste0(formatted, " [n = ", n_risk_fmt, "]")
    }

    ## Add n_event if requested
    if ("n_event" %in% stats) {
        n_event_fmt <- if (!is.null(marks)) {
            vapply(n_event[valid_idx], format_count, character(1), marks = marks)
        } else {
            as.character(n_event[valid_idx])
        }
        formatted <- paste0(formatted, " [e = ", n_event_fmt, "]")
    }

    result[valid_idx] <- formatted
    result[na_mask] <- "-"

    result
}


#' Quantile cell formatting
#'
#' Formats survival quantile cells for multiple rows at once.
#'
#' @param est Numeric vector of estimates.
#' @param lower Numeric vector of lower CI bounds.
#' @param upper Numeric vector of upper CI bounds.
#' @param fmt_str Format string for numeric values.
#' @param marks List with \code{big.mark} and \code{decimal.mark} as returned
#'   by \code{\link{resolve_number_marks}}.
#' @return Character vector of formatted cells.
#' @keywords internal
format_quantile_cells <- function(est, lower, upper, fmt_str,
                                             marks = NULL) {
    n <- length(est)
    result <- character(n)

    ## Handle NA (not reached) values
    na_mask <- is.na(est)

    if (all(na_mask)) {
        return(rep("NR", n))
    }

    ## Format non-NA values
    valid_idx <- which(!na_mask)

    est_fmt <- sprintf(fmt_str, est[valid_idx])
    lower_fmt <- ifelse(is.na(lower[valid_idx]), "NR", sprintf(fmt_str, lower[valid_idx]))
    upper_fmt <- ifelse(is.na(upper[valid_idx]), "NR", sprintf(fmt_str, upper[valid_idx]))

    if (!is.null(marks)) {
        ## Apply locale decimal mark and fix negative zeros
        est_fmt <- apply_decimal_mark(est_fmt, marks)
        lower_fmt <- ifelse(lower_fmt == "NR", "NR",
                            apply_decimal_mark(lower_fmt, marks))
        upper_fmt <- ifelse(upper_fmt == "NR", "NR",
                            apply_decimal_mark(upper_fmt, marks))

        ## Determine CI separator
        ## Use worst-case across all valid elements for consistency within a column
        any_negative <- any(lower[valid_idx] < 0, upper[valid_idx] < 0, na.rm = TRUE)
        if (any_negative) {
            sep <- " to "
        } else if (marks$decimal.mark == ",") {
            sep <- "\u2013"
        } else {
            sep <- "-"
        }
    } else {
        ## Fallback: check for negatives only
        any_negative <- any(lower[valid_idx] < 0, upper[valid_idx] < 0, na.rm = TRUE)
        sep <- if (any_negative) " to " else "-"
    }

    result[valid_idx] <- paste0(est_fmt, " (", lower_fmt, sep, upper_fmt, ")")
    result[na_mask] <- "NR"

    result
}


#' Add \emph{p}-value column to result table
#'
#' Adds formatted \emph{p}-value column to the survtable result.
#'
#' @param result Data.table result.
#' @param p_value Numeric \emph{p}-value.
#' @param p_digits Integer decimal places for \emph{p}-value.
#' @param marks List with \code{big.mark} and \code{decimal.mark} as returned
#'   by \code{\link{resolve_number_marks}}.
#' @return Data.table with \emph{p}-value column added.
#' @keywords internal
add_pvalue_column <- function(result, p_value, p_digits, marks = NULL) {
    n_rows <- nrow(result)
    p_col <- c(format_pvalue_survtable(p_value, p_digits, marks), rep("", n_rows - 1L))
    result[, `p-value` := p_col]
    result
}


#' Format \emph{p}-value for survtable
#'
#' Provides \emph{p}-value formatting to the survtable result.
#'
#' @param p Numeric \emph{p}-value.
#' @param digits Integer decimal places.
#' @param marks List with \code{big.mark} and \code{decimal.mark} as returned
#'   by \code{\link{resolve_number_marks}}.
#' @return Character-formatted \emph{p}-value.
#' @keywords internal
format_pvalue_survtable <- function(p, digits, marks = NULL) {
    if (is.na(p)) {
        return("-")
    }

    ## Use marks-aware formatter if available
    if (!is.null(marks)) {
        return(format_pvalue(p, digits, marks))
    }

    ## Fallback (e.g., from print method without marks)
    threshold <- 10^(-digits)
    if (p < threshold) {
        return(paste0("< ", format(threshold, scientific = FALSE)))
    }

    sprintf(paste0("%.", digits, "f"), p)
}
