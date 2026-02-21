#' Create Publication-Ready Survival Summary Tables
#'
#' Generates comprehensive survival summary tables with survival probabilities
#' at specified time points, median survival times, and optional group
#' comparisons with statistical testing. Designed for creating survival
#' summaries commonly used in clinical and epidemiological research
#' publications.
#'
#' @param data Data frame or data.table containing the survival dataset.
#'   Automatically converted to a data.table for efficient processing.
#'
#' @param outcome Character string or character vector specifying one or more
#'   survival outcomes using \code{Surv()} syntax (\emph{e.g.,}
#'   \code{"Surv(os_months, os_status)"}). When multiple outcomes are
#'   provided, results are stacked into a single table with outcome labels
#'   as row headers.
#'
#' @param by Character string specifying the column name of the stratifying
#'   variable for group comparisons (\emph{e.g.,} treatment arm, risk group). When
#'   \code{NULL} (default), produces overall survival summaries only.
#'
#' @param times Numeric vector of time points at which to estimate survival
#'   probabilities. For example, \code{c(12, 24, 36)} for 1-, 2-, and 3-year
#'   survival when time is measured in months. Default is \code{NULL}.
#'
#' @param probs Numeric vector of survival probabilities for which to estimate
#'   corresponding survival times (quantiles). Values must be between 0 and 1.
#'   For example, \code{c(0.5)} returns median survival time,
#'   \code{c(0.25, 0.5, 0.75)} returns quartiles. Default is \code{0.5}
#'   (median only).
#'
#' @param stats Character vector specifying which statistics to display:
#'   \itemize{
#'     \item \code{"survival"} - Survival probability at specified times
#'     \item \code{"ci"} - Confidence interval for survival probability
#'     \item \code{"n_risk"} - Number at risk at each time point
#'     \item \code{"n_event"} - Cumulative number of events by each time
#'       point
#'   }
#'   Default is \code{c("survival", "ci")}.
#'
#' @param type Character string specifying the type of probability to report:
#'   \itemize{
#'     \item \code{"survival"} - Survival probability S(t) [default]
#'     \item \code{"risk"} - Cumulative incidence/risk 1 - S(t)
#'     \item \code{"cumhaz"} - Cumulative hazard -log(S(t))
#'   }
#'
#' @param conf_level Numeric confidence level for confidence intervals. Must
#'   be between 0 and 1. Default is 0.95 (95\% confidence intervals).
#'
#' @param conf_type Character string specifying the confidence interval type
#'   for survival estimates:
#'   \itemize{
#'     \item \code{"log"} - Log transformation (default, recommended)
#'     \item \code{"log-log"} - Log-log transformation
#'     \item \code{"plain"} - Linear/identity (can produce CIs outside
#'       [0, 1])
#'     \item \code{"logit"} - Logit transformation
#'     \item \code{"arcsin"} - Arcsin square root transformation
#'   }
#'
#' @param digits Integer specifying the number of decimal places for survival
#'   probabilities (as percentages). Default is 0 (whole percentages).
#'
#' @param time_digits Integer specifying the number of decimal places for
#'   survival time estimates (median, quantiles). Default is 1.
#'
#' @param p_digits Integer specifying the number of decimal places for
#'   \emph{p}-values. Values smaller than \code{10^(-p_digits)} are displayed
#'   as \code{"< 0.001"} (for \code{p_digits = 3}), \code{"< 0.0001"} (for
#'   \code{p_digits = 4}), etc. Default is 3.
#'
#' @param percent Logical. If \code{TRUE} (default), displays survival
#'   probabilities as percentages (\emph{e.g.,} \code{"85\%"}). If \code{FALSE},
#'   displays as proportions (\emph{e.g.,} \code{"0.85"}).
#'
#' @param test Logical. If \code{TRUE} (default), performs a survival curve
#'   comparison test and adds a \emph{p}-value column. Requires \code{by}
#'   to be specified.
#'
#' @param test_type Character string specifying the statistical test for
#'   comparing survival curves:
#'   \itemize{
#'     \item \code{"logrank"} - Log-rank test (default)
#'     \item \code{"wilcoxon"} - Wilcoxon (Breslow) test
#'     \item \code{"tarone"} - Tarone-Ware test
#'     \item \code{"petopeto"} - Peto-Peto test
#'   }
#'
#' @param total Logical or character string controlling the total/overall
#'   column:
#'   \itemize{
#'     \item \code{TRUE} or \code{"first"} - Include total column first
#'       [default]
#'     \item \code{"last"} - Include total column last (before
#'       \emph{p}-value)
#'     \item \code{FALSE} - Exclude total column
#'   }
#'
#' @param total_label Character string for the total/overall row label.
#'   Default is \code{"Total"}.
#'
#' @param time_unit Character string specifying the time unit for display
#'   in column headers and labels (\emph{e.g.,} \code{"months"}, \code{"days"},
#'   \code{"years"}). When specified, time column headers become
#'   "\{time\} \{time_unit\}" (\emph{e.g.,} "12 months"). Default is \code{NULL}
#'   (no unit shown).
#'
#' @param time_label Character string template for time column headers when
#'   \code{times} is specified. Use \code{"\{time\}"} as placeholder for
#'   the time value and \code{"\{unit\}"} for the time unit. Default is
#'   \code{"\{time\} \{unit\}"} when \code{time_unit} is specified,
#'   otherwise just \code{"\{time\}"}.
#'
#' @param median_label Character string for the median survival row label.
#'   Default is \code{NULL}, which auto-constructs from \code{conf_level}
#'   (\emph{e.g.,} \code{"Median (95\% CI)"} for \code{conf_level = 0.95}).
#'
#' @param labels Named character vector or list providing custom display
#'   labels. For stratified analyses, names should match levels of the
#'   \code{by} variable. For multiple outcomes, names should match the
#'   \code{Surv()} expressions. Default is \code{NULL}.
#'
#' @param by_label Character string providing a custom label for the
#'   stratifying variable (used in output attributes and headers).
#'   Default is \code{NULL} (uses variable name).
#'
#' @param na_rm Logical. If \code{TRUE} (default), observations with missing
#'   values in time, status, or the stratifying variable are excluded.
#'
#' @param number_format Character string or two-element character vector
#'   controlling thousand and decimal separators in formatted output. Named
#'   presets:
#'   \itemize{
#'     \item \code{"us"} - Comma thousands, period decimal: \code{1,234.56} [default]
#'     \item \code{"eu"} - Period thousands, comma decimal: \code{1.234,56}
#'     \item \code{"space"} - Thin-space thousands, period decimal: \code{1 234.56}
#'       (SI/ISO 31-0)
#'     \item \code{"none"} - No thousands separator: \code{1234.56}
#'   }
#'   Or provide a custom two-element vector \code{c(big.mark, decimal.mark)},
#'   \emph{e.g.}, \code{c("'", ".")} for Swiss-style: \verb{1'234.56}.
#'
#'   When \code{NULL} (default), uses
#'   \code{getOption("summata.number_format", "us")}. Set the global option
#'   once per session to avoid passing this argument repeatedly:
#'   \preformatted{
#'     options(summata.number_format = "eu")
#'   }
#'
#' @param ... Additional arguments passed to
#'   \code{\link[survival]{survfit}}.
#'
#' @return A data.table with S3 class \code{"survtable"} containing formatted
#'   survival statistics. The table structure depends on parameters:
#'
#'   \strong{When \code{times} is specified (survival at time points):}
#'   \describe{
#'     \item{Variable/Group}{Row identifier -- stratifying variable levels}
#'     \item{Time columns}{Survival statistics at each requested time point}
#'     \item{\emph{p}-value}{Test \emph{p}-value (if \code{test = TRUE} and
#'       \code{by} specified)}
#'   }
#'
#'   \strong{When only \code{probs} is specified (survival quantiles):}
#'   \describe{
#'     \item{Variable/Group}{Row identifier -- stratifying variable levels}
#'     \item{Quantile columns}{Time to reach each survival probability}
#'     \item{\emph{p}-value}{Test \emph{p}-value (if \code{test = TRUE} and
#'       \code{by} specified)}
#'   }
#'
#'   All numeric output (probabilities, times, counts, \emph{p}-values)
#'   respects the \code{number_format} setting for locale-appropriate
#'   formatting.
#'
#'   The returned object includes the following attributes:
#'   \describe{
#'     \item{raw_data}{Data.table with unformatted numeric values}
#'     \item{survfit_objects}{List of survfit objects for each stratum}
#'     \item{by_variable}{The stratifying variable name}
#'     \item{times}{The time points requested}
#'     \item{probs}{The probability quantiles requested}
#'     \item{test_result}{Full test result object (if test performed)}
#'   }
#'
#' @details
#' \strong{Survival Probability Estimation:}
#'
#' Survival probabilities are estimated using the Kaplan-Meier method via
#' \code{\link[survival]{survfit}}. At each specified time point, the function
#' reports the estimated probability of surviving beyond that time.
#'
#' \strong{Confidence Intervals:}
#'
#' The default \code{"log"} transformation for confidence intervals is
#' recommended as it ensures intervals remain within [0, 1] and has good
#' statistical properties. The \code{"log-log"} transformation is also
#' commonly used and may perform better in the tails.
#'
#' \strong{Statistical Testing:}
#'
#' The log-rank test (default) tests the null hypothesis that survival curves
#' are identical across groups. Alternative tests weight different parts of
#' the survival curve:
#' \itemize{
#'   \item Log-rank: Equal weights (best for proportional hazards)
#'   \item Wilcoxon: Weights by number at risk (sensitive to early
#'     differences)
#'   \item Tarone-Ware: Weights by square root of number at risk
#'   \item Peto-Peto: Modified Wilcoxon weights
#' }
#'
#' \strong{Formatting:}
#'
#' All numeric output respects the \code{number_format} parameter.
#' Separators within confidence intervals adapt automatically to avoid
#' ambiguity:
#' \itemize{
#'   \item Survival probabilities: \code{"85\% (80\%-89\%)"} (US) or
#'     \code{"85\% (80\%–89\%)"} (EU, en-dash separator)
#'   \item Median survival: \code{"24.5 (21.2-28.9)"} (US) or
#'     \code{"24,5 (21,2–28,9)"} (EU)
#'   \item Counts ≥ 1000: \code{"1,234"} (US) or
#'     \code{"1.234"} (EU)
#'   \item \emph{p}-values: \code{"< 0.001"} (US) or
#'     \code{"< 0,001"} (EU)
#' }
#'
#' @seealso
#' \code{\link{desctable}} for baseline characteristics tables,
#' \code{\link{fit}} for regression analysis,
#' \code{\link{table2pdf}} for PDF export,
#' \code{\link{table2docx}} for Word export,
#' \code{\link[survival]{survfit}} for underlying survival estimation,
#' \code{\link[survival]{survdiff}} for survival curve comparison tests
#'
#' @examples
#' # Load example data
#' data(clintrial)
#'
#' # Example 1: Survival at specific time points by treatment
#' survtable(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     by = "treatment",
#'     times = c(12, 24, 36),
#'     time_unit = "months"
#' )
#'
#' \donttest{
#' # Example 2: Median survival only
#' survtable(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     by = "treatment",
#'     times = NULL,
#'     probs = 0.5
#' )
#'
#' # Example 3: Multiple quantiles (quartiles)
#' survtable(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     by = "stage",
#'     times = NULL,
#'     probs = c(0.25, 0.5, 0.75)
#' )
#'
#' # Example 4: Both time points and median
#' survtable(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     by = "treatment",
#'     times = c(12, 24),
#'     probs = 0.5,
#'     time_unit = "months"
#' )
#'
#' # Example 5: Cumulative incidence (1 - survival)
#' survtable(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     by = "treatment",
#'     times = c(12, 24),
#'     type = "risk"
#' )
#'
#' # Example 6: Include number at risk
#' survtable(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     by = "treatment",
#'     times = c(12, 24),
#'     stats = c("survival", "ci", "n_risk")
#' )
#'
#' # Example 7: Overall survival without stratification
#' survtable(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     times = c(12, 24, 36, 48)
#' )
#'
#' # Example 8: Without total row
#' survtable(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     by = "treatment",
#'     times = c(12, 24),
#'     total = FALSE
#' )
#'
#' # Example 9: Custom labels
#' survtable(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     by = "treatment",
#'     times = c(12, 24),
#'     labels = c("Drug A" = "Treatment A", "Drug B" = "Treatment B"),
#'     time_unit = "months"
#' )
#'
#' # Example 10: Different confidence interval type
#' survtable(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     by = "treatment",
#'     times = c(12, 24),
#'     conf_type = "log-log"
#' )
#'
#' # Example 11: Wilcoxon test instead of log-rank
#' survtable(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     by = "treatment",
#'     times = c(12, 24),
#'     test_type = "wilcoxon"
#' )
#'
#' # Example 12: Access raw data for custom analysis
#' result <- survtable(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     by = "treatment",
#'     times = c(12, 24)
#' )
#' raw <- attr(result, "raw_data")
#' print(raw)
#'
#' # Example 13: Access survfit objects for plotting
#' fits <- attr(result, "survfit_objects")
#' # plot(fits$overall)  # Plot overall survival curve
#'
#' # Example 14: Multiple survival outcomes stacked
#' survtable(
#'     data = clintrial,
#'     outcome = c("Surv(pfs_months, pfs_status)", "Surv(os_months, os_status)"),
#'     by = "treatment",
#'     times = c(12, 24),
#'     probs = 0.5,
#'     time_unit = "months",
#'     total = FALSE,
#'     labels = c(
#'         "Surv(pfs_months, pfs_status)" = "Progression-Free Survival",
#'         "Surv(os_months, os_status)" = "Overall Survival"
#'     )
#' )
#'
#' # Example 15: European number formatting
#' survtable(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     by = "treatment",
#'     times = c(12, 24),
#'     number_format = "eu"
#' )
#'
#' }
#'
#' @family descriptive functions
#' @export
survtable <- function(data,
                      outcome,
                      by = NULL,
                      times = NULL,
                      probs = 0.5,
                      stats = c("survival", "ci"),
                      type = "survival",
                      conf_level = 0.95,
                      conf_type = "log",
                      digits = 0,
                      time_digits = 1,
                      p_digits = 3,
                      percent = TRUE,
                      test = TRUE,
                      test_type = "logrank",
                      total = TRUE,
                      total_label = "Total",
                      time_unit = NULL,
                      time_label = NULL,
                      median_label = NULL,
                      labels = NULL,
                      by_label = NULL,
                      na_rm = TRUE,
                      number_format = NULL,
                      ...) {

    ## Validate inputs
    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }

    ## Resolve number formatting marks
    validate_number_format(number_format)
    marks <- resolve_number_marks(number_format)

    if (is.null(times) && is.null(probs)) {
        stop("At least one of 'times' or 'probs' must be specified")
    }

    if (!is.null(probs)) {
        if (any(probs <= 0 | probs >= 1)) {
            stop("'probs' must contain values between 0 and 1 (exclusive)")
        }
    }

    ## Build time_label if not provided
    if (is.null(time_label)) {
        if (!is.null(time_unit)) {
            time_label <- "{time} {unit}"
        } else {
            time_label <- "{time}"
        }
    }

    type <- match.arg(type, c("survival", "risk", "cumhaz"))
    conf_type <- match.arg(conf_type, c("log", "log-log", "plain", "logit", "arcsin"))
    test_type <- match.arg(test_type, c("logrank", "wilcoxon", "tarone", "petopeto"))
    stats <- match.arg(stats, c("survival", "ci", "n_risk", "n_event"), several.ok = TRUE)

    if (!requireNamespace("survival", quietly = TRUE)) {
        stop("Package 'survival' required for survival analysis")
    }

    ## Auto-construct median_label from conf_level if not provided
    if (is.null(median_label)) {
        ci_pct <- round(conf_level * 100)
        median_label <- paste0("Median (", ci_pct, "% CI)")
    }

    ## Handle multiple outcomes
    if (length(outcome) > 1) {
        ## Process each outcome and stack results
        all_results <- vector("list", length(outcome))
        all_raw <- vector("list", length(outcome))
        all_survfit <- list()
        all_tests <- list()

        for (i in seq_along(outcome)) {
            ## Get label for this outcome
            outcome_label <- if (!is.null(labels) && outcome[i] %in% names(labels)) {
                                 labels[[outcome[i]]]
                             } else {
                                 outcome[i]
                             }

            ## Process single outcome
            single_result <- process_single_outcome(
                data = data,
                outcome = outcome[i],
                outcome_label = outcome_label,
                by = by,
                times = times,
                probs = probs,
                stats = stats,
                type = type,
                conf_level = conf_level,
                conf_type = conf_type,
                digits = digits,
                time_digits = time_digits,
                percent = percent,
                test = test,
                test_type = test_type,
                total = total,
                total_label = total_label,
                time_unit = time_unit,
                time_label = time_label,
                median_label = median_label,
                labels = labels,
                na_rm = na_rm,
                marks = marks,
                ...
            )

            all_results[[i]] <- single_result$formatted
            all_raw[[i]] <- single_result$raw
            all_survfit[[outcome[i]]] <- single_result$survfit_objects
            all_tests[[outcome[i]]] <- single_result$test_result
        }

        ## Stack all results
        result <- data.table::rbindlist(all_results, fill = TRUE)
        raw_result <- data.table::rbindlist(all_raw, fill = TRUE)

        ## Set attributes for multiple outcomes
        data.table::setattr(result, "raw_data", raw_result)
        data.table::setattr(result, "survfit_objects", all_survfit)
        data.table::setattr(result, "outcome", outcome)
        data.table::setattr(result, "by_variable", by)
        data.table::setattr(result, "by_label", by_label)
        data.table::setattr(result, "times", times)
        data.table::setattr(result, "probs", probs)
        data.table::setattr(result, "time_unit", time_unit)
        data.table::setattr(result, "test_result", all_tests)
        data.table::setattr(result, "type", type)

        class(result) <- c("survtable", class(result))
        return(result)
    }

    ## Single outcome - process directly
    single_result <- process_single_outcome(
        data = data,
        outcome = outcome,
        outcome_label = NULL,  
        by = by,
        times = times,
        probs = probs,
        stats = stats,
        type = type,
        conf_level = conf_level,
        conf_type = conf_type,
        digits = digits,
        time_digits = time_digits,
        percent = percent,
        test = test,
        test_type = test_type,
        total = total,
        total_label = total_label,
        time_unit = time_unit,
        time_label = time_label,
        median_label = median_label,
        labels = labels,
        na_rm = na_rm,
        marks = marks,
        ...
    )

    result <- single_result$formatted
    raw_result <- single_result$raw
    
    ## For single outcome, rename Group to Variable for consistency with desctable
    ## (Group column contains stratum labels which become the row identifiers)
    if ("Group" %in% names(result) && !"Variable" %in% names(result)) {
        data.table::setnames(result, "Group", "Variable")
        data.table::setnames(raw_result, "Group", "Variable")
    }

    ## Set attributes
    data.table::setattr(result, "raw_data", raw_result)
    data.table::setattr(result, "survfit_objects", single_result$survfit_objects)
    data.table::setattr(result, "outcome", outcome)
    data.table::setattr(result, "by_variable", by)
    data.table::setattr(result, "by_label", by_label)
    data.table::setattr(result, "times", times)
    data.table::setattr(result, "probs", probs)
    data.table::setattr(result, "time_unit", time_unit)
    data.table::setattr(result, "test_result", single_result$test_result)
    data.table::setattr(result, "type", type)

    class(result) <- c("survtable", class(result))

    return(result)
}


#' Process a single survival outcome
#' @keywords internal
process_single_outcome <- function(data,
                                   outcome,
                                   outcome_label,
                                   by,
                                   times,
                                   probs,
                                   stats,
                                   type,
                                   conf_level,
                                   conf_type,
                                   digits,
                                   time_digits,
                                   percent,
                                   test,
                                   test_type,
                                   total,
                                   total_label,
                                   time_unit,
                                   time_label,
                                   median_label,
                                   labels,
                                   na_rm,
                                   marks = NULL,
                                   ...) {

    ## Parse Surv() expression
    surv_match <- regexec("Surv\\(([^,]+),\\s*([^)]+)\\)", outcome)
    surv_parts <- regmatches(outcome, surv_match)[[1]]

    if (length(surv_parts) < 3) {
        stop("Invalid Surv() syntax: ", outcome, 
             "\nExpected format: Surv(time_var, status_var)")
    }

    time <- trimws(surv_parts[2])
    status <- trimws(surv_parts[3])

    if (!time %in% names(data)) {
        stop("Time variable '", time, "' not found in data")
    }
    if (!status %in% names(data)) {
        stop("Status variable '", status, "' not found in data")
    }
    if (!is.null(by) && !by %in% names(data)) {
        stop("Stratifying variable '", by, "' not found in data")
    }

    ## Handle missing data
    if (na_rm) {
        complete_vars <- c(time, status)
        if (!is.null(by)) complete_vars <- c(complete_vars, by)
        complete_idx <- complete.cases(data[, ..complete_vars])
        data <- data[complete_idx]
    }

    ## Create survival object
    surv_obj <- survival::Surv(data[[time]], data[[status]])

    ## Store survfit objects
    survfit_objects <- list()

    ## Get group information
    groups <- NULL
    group_labels <- NULL
    if (!is.null(by)) {
        grp_vec <- data[[by]]
        if (is.factor(grp_vec)) {
            groups <- levels(grp_vec)
        } else {
            groups <- sort(unique(grp_vec))
            groups <- groups[!is.na(groups)]
        }

        ## Apply labels to groups
        group_labels <- as.character(groups)
        if (!is.null(labels)) {
            for (i in seq_along(groups)) {
                grp_char <- as.character(groups[i])
                if (grp_char %in% names(labels)) {
                    group_labels[i] <- labels[[grp_char]]
                }
            }
        }
    }

    ## Fit survival models
    if (!is.null(by)) {
        ## Stratified analysis
        fit_formula <- stats::as.formula(paste("surv_obj ~", by))
        fit_stratified <- survival::survfit(fit_formula, data = data,
                                            conf.int = conf_level,
                                            conf.type = conf_type, ...)
        survfit_objects$stratified <- fit_stratified

        ## Also fit overall
        if (!isFALSE(total)) {
            fit_overall <- survival::survfit(surv_obj ~ 1, data = data,
                                             conf.int = conf_level,
                                             conf.type = conf_type, ...)
            survfit_objects$overall <- fit_overall
        }
    } else {
        ## Overall only
        fit_overall <- survival::survfit(surv_obj ~ 1, data = data,
                                         conf.int = conf_level,
                                         conf.type = conf_type, ...)
        survfit_objects$overall <- fit_overall
    }

    ## Perform statistical test
    test_result <- NULL
    p_value <- NULL
    if (test && !is.null(by) && length(groups) >= 2) {
        test_result <- perform_survival_test(surv_obj, data[[by]], test_type)
        p_value <- test_result$p_value
    }

    ## Build results table
    result_list <- list()
    raw_list <- list()

    ## Process time points (survival probabilities)
    if (!is.null(times)) {
        time_results <- process_survival_times(
            survfit_objects = survfit_objects,
            times = times,
            groups = groups,
            group_labels = group_labels,
            stats = stats,
            type = type,
            digits = digits,
            percent = percent,
            total = total,
            total_label = total_label,
            time_label = time_label,
            time_unit = time_unit,
            by = by,
            data = data,
            marks = marks
        )
        result_list$times <- time_results$formatted
        raw_list$times <- time_results$raw
    }

    ## Process probability quantiles (median, quartiles, etc.)
    if (!is.null(probs)) {
        prob_results <- process_survival_probs(
            survfit_objects = survfit_objects,
            probs = probs,
            groups = groups,
            group_labels = group_labels,
            time_digits = time_digits,
            total = total,
            total_label = total_label,
            median_label = median_label,
            by = by,
            data = data,
            conf_level = conf_level,
            marks = marks
        )
        result_list$probs <- prob_results$formatted
        raw_list$probs <- prob_results$raw
    }

    ## Combine results - merge times and probs by Group
    if (length(result_list) == 1) {
        result <- result_list[[1]]
        raw_result <- raw_list[[1]]
    } else {
        ## Both times and probs - merge by Group column
        result <- merge(result_list$times, result_list$probs, 
                        by = "Group", all = TRUE, sort = FALSE)
        raw_result <- merge(raw_list$times, raw_list$probs,
                            by = "Group", all = TRUE, sort = FALSE)
        
        ## Restore original row order (Total first if present, then groups)
        if (!isFALSE(total) && total_label %in% result$Group) {
            total_row <- result[Group == total_label]
            other_rows <- result[Group != total_label]
            if (isTRUE(total) || identical(total, "first")) {
                result <- rbind(total_row, other_rows)
            } else {
                result <- rbind(other_rows, total_row)
            }
            
            ## Same for raw
            total_row_raw <- raw_result[Group == total_label]
            other_rows_raw <- raw_result[Group != total_label]
            if (isTRUE(total) || identical(total, "first")) {
                raw_result <- rbind(total_row_raw, other_rows_raw)
            } else {
                raw_result <- rbind(other_rows_raw, total_row_raw)
            }
        }
    }

    ## Add p-value column if test performed
    if (!is.null(p_value)) {
        result <- add_pvalue_column(result, p_value, 3, marks)  # p_digits handled in main
    }

    ## Rename Group to match desctable pattern when stacking multiple outcomes
    ## For multi-outcome tables: Variable = outcome label, Group = stratum label
    if (!is.null(outcome_label)) {
        ## Create Variable column with outcome label only in first row
        n_rows <- nrow(result)
        result[, Variable := c(outcome_label, rep("", n_rows - 1))]
        raw_result[, Variable := c(outcome_label, rep("", n_rows - 1))]
        
        ## Rename existing Group column to preserve stratum labels
        ## (Group column already exists from process_survival_times/probs)
        
        ## Reorder columns: Variable first, then Group, then rest
        other_cols <- setdiff(names(result), c("Variable", "Group"))
        data.table::setcolorder(result, c("Variable", "Group", other_cols))
        
        other_cols_raw <- setdiff(names(raw_result), c("Variable", "Group"))
        data.table::setcolorder(raw_result, c("Variable", "Group", other_cols_raw))
    }

    list(
        formatted = result,
        raw = raw_result,
        survfit_objects = survfit_objects,
        test_result = test_result
    )
}


#' Print method for survtable
#' @family descriptive functions
#' @export
#' @keywords internal
print.survtable <- function(x, ...) {
    outcome <- attr(x, "outcome")
    by_var <- attr(x, "by_variable")
    by_label <- attr(x, "by_label")
    times <- attr(x, "times")
    probs <- attr(x, "probs")
    time_unit <- attr(x, "time_unit")
    type <- attr(x, "type")
    test_result <- attr(x, "test_result")

    cat("\nSurvival Summary Table\n")
    
    if (!is.null(outcome)) {
        if (length(outcome) > 1) {
            cat("Outcomes: ", length(outcome), "\n", sep = "")
        } else {
            cat("Outcome: ", outcome, "\n", sep = "")
        }
    }

    if (!is.null(by_var)) {
        display_by <- if (!is.null(by_label)) by_label else by_var
        cat("Stratified by: ", display_by, "\n", sep = "")
    }

    if (!is.null(times)) {
        times_str <- paste(times, collapse = ", ")
        if (!is.null(time_unit)) {
            times_str <- paste(times_str, time_unit)
        }
        cat("Time points: ", times_str, "\n", sep = "")
    }
    if (!is.null(probs)) {
        cat("Quantiles: ", paste(probs * 100, "%", sep = "", collapse = ", "), "\n", sep = "")
    }

    type_label <- switch(type,
                         "survival" = "Survival probability",
                         "risk" = "Cumulative incidence",
                         "cumhaz" = "Cumulative hazard")
    cat("Statistic: ", type_label, "\n", sep = "")

    ## Handle test result - could be single or list for multiple outcomes
    if (!is.null(test_result)) {
        if (is.list(test_result) && !is.null(test_result$test_type)) {
            ## Single test result
            test_name <- switch(test_result$test_type,
                                "logrank" = "Log-rank",
                                "wilcoxon" = "Wilcoxon",
                                "tarone" = "Tarone-Ware",
                                "petopeto" = "Peto-Peto")
            cat("Test: ", test_name, " (p = ", format_pvalue_survtable(test_result$p_value, 3), ")\n", sep = "")
        } else if (is.list(test_result) && length(test_result) > 0) {
            ## Multiple test results - show first one's type
            first_test <- test_result[[1]]
            if (!is.null(first_test)) {
                test_name <- switch(first_test$test_type,
                                    "logrank" = "Log-rank",
                                    "wilcoxon" = "Wilcoxon",
                                    "tarone" = "Tarone-Ware",
                                    "petopeto" = "Peto-Peto")
                cat("Test: ", test_name, "\n", sep = "")
            }
        }
    }

    cat("\n")
    NextMethod("print", x)
    invisible(x)
}
