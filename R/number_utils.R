#' Number Formatting Utilities
#'
#' Internal utilities for locale-aware number formatting across all summata
#' output functions. Supports preset locales (US, European, SI/ISO, plain) and
#' fully custom separator definitions.
#'
#' @section Global Option:
#' The default number format can be set once per session:
#' \preformatted{
#'   options(summata.number_format = "eu")
#' }
#' This avoids passing \code{number_format} to every function call.
#'
#' @name number_format
#' @keywords internal
NULL


#' Resolve number format marks
#'
#' Converts a \code{number_format} specification into a list of \code{big.mark}
#' and \code{decimal.mark} values used by all downstream formatting functions.
#' Supports named presets, custom two-element vectors, and the global
#' \code{summata.number_format} option.
#'
#' @param number_format Character string specifying a named preset, a
#'   two-element character vector \code{c(big.mark, decimal.mark)}, or
#'   \code{NULL} to use the global option (falling back to \code{"us"}).
#'
#'   Named presets:
#'   \describe{
#'     \item{\code{"us"}}{Comma thousands, period decimal: 1,234.56}
#'     \item{\code{"eu"}}{Period thousands, comma decimal: 1.234,56}
#'     \item{\code{"space"}}{Thin-space thousands, period decimal: 1 234.56
#'       (SI/ISO 31-0 standard)}
#'     \item{\code{"none"}}{No thousands separator, period decimal: 1234.56}
#'   }
#'
#'   Custom vector: \code{c(",", ".")} or \code{c(".", ",")} \emph{etc.} The first
#'   element is \code{big.mark}, the second is \code{decimal.mark}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{big.mark}}{Character string for thousands separator.}
#'     \item{\code{decimal.mark}}{Character string for decimal separator.}
#'   }
#' 
#' @keywords internal
resolve_number_marks <- function(number_format = NULL) {

    ## Fall back to global option, then to "us"
    if (is.null(number_format)) {
        number_format <- getOption("summata.number_format", "us")
    }

    ## Custom vector: c(big.mark, decimal.mark)
    if (is.character(number_format) && length(number_format) == 2L) {
        return(list(big.mark = number_format[1L], decimal.mark = number_format[2L]))
    }

    ## Named presets (use thin space U+202F for SI style)
    if (is.character(number_format) && length(number_format) == 1L) {
        return(switch(number_format,
            "us"    = list(big.mark = ",",      decimal.mark = "."),
            "eu"  = list(big.mark = ".",      decimal.mark = ","),
            "space" = list(big.mark = "\u202F", decimal.mark = "."),
            "none"  = list(big.mark = "",       decimal.mark = "."),
            stop("Unknown number_format preset: '", number_format,
                 "'. Use 'us', 'eu', 'space', 'none', or a ",
                 "two-element character vector c(big.mark, decimal.mark).",
                 call. = FALSE)
        ))
    }

    stop("'number_format' must be a character string preset or a ",
         "two-element character vector.", call. = FALSE)
}


#' Format a numeric value with locale-aware separators
#'
#' General-purpose number formatter used by all display functions. For values
#' ≥ 1000 (in absolute value), inserts the appropriate
#' thousands separator. Fixes negative-zero display artefacts.
#'
#' @param x Numeric value to format.
#' @param fmt_str Character string \code{sprintf} format specification
#'   (\emph{e.g.,} \code{"\%.1f"}). Used only for |x| < 1000.
#' @param marks List with \code{big.mark} and \code{decimal.mark} as returned
#'   by \code{\link{resolve_number_marks}}.
#' @return Character string with the formatted number.
#' @keywords internal
format_num <- function(x, fmt_str, marks) {
    if (is.na(x)) return("")
    if (abs(x) >= 1000) {
        formatted <- format(round(x, 1), big.mark = marks$big.mark,
                            scientific = FALSE)
        ## Replace decimal mark if non-standard
        if (marks$decimal.mark != ".") {
            formatted <- sub(".", marks$decimal.mark, formatted, fixed = TRUE)
        }
        trimws(formatted)
    } else {
        apply_decimal_mark(sprintf(fmt_str, x), marks)
    }
}


#' Format an integer count with locale-aware thousands separator
#'
#' Formats integer values with a thousands separator for display in tables.
#' Values below 1000 are returned as plain character strings.
#'
#' @param n Integer count value.
#' @param marks List with \code{big.mark} and \code{decimal.mark} as returned
#'   by \code{\link{resolve_number_marks}}.
#' @return Character string with the formatted count.
#' @keywords internal
format_count <- function(n, marks) {
    if (n >= 1000) {
        trimws(format(n, big.mark = marks$big.mark,
                      decimal.mark = marks$decimal.mark))
    } else {
        as.character(n)
    }
}


#' Format a categorical statistic for display
#'
#' Converts frequency counts into formatted display strings following
#' standard conventions (n, n (\%), \% only) with locale-aware decimal marks.
#'
#' @param n Integer count for the category.
#' @param total Integer total count for percentage calculation.
#' @param stat_type Character string: \code{"n"}, \code{"n_percent"}, or
#'   \code{"percent"}.
#' @param marks List with \code{big.mark} and \code{decimal.mark} as returned
#'   by \code{\link{resolve_number_marks}}.
#' @return Character string with the formatted statistic.
#' @keywords internal
format_categorical_stat <- function(n, total, stat_type, marks) {
    if (length(stat_type) > 1) stat_type <- stat_type[1]

    n_formatted <- format_count(n, marks)
    pct <- 100 * n / total

    switch(stat_type,
           "n" = n_formatted,
           "percent" = {
               pct_str <- sprintf("%.1f", pct)
               if (marks$decimal.mark != ".") {
                   pct_str <- sub(".", marks$decimal.mark, pct_str, fixed = TRUE)
               }
               paste0(pct_str, "%")
           },
           "n_percent" = {
               pct_str <- sprintf("%.1f", pct)
               if (marks$decimal.mark != ".") {
                   pct_str <- sub(".", marks$decimal.mark, pct_str, fixed = TRUE)
               }
               paste0(n_formatted, " (", pct_str, "%)")
           },
           n_formatted
    )
}


#' Format a \emph{p}-value with locale-aware decimal mark
#'
#' Converts a numeric \emph{p}-value to a display string with the correct
#' decimal separator and threshold notation (\emph{e.g.,} "< 0.001" or "< 0,001").
#'
#' @param p Numeric \emph{p}-value.
#' @param digits Integer number of decimal places.
#' @param marks List with \code{big.mark} and \code{decimal.mark} as returned
#'   by \code{\link{resolve_number_marks}}.
#' @return Character string with the formatted \emph{p}-value.
#' @keywords internal
format_pvalue <- function(p, digits, marks) {
    if (is.na(p)) return("")

    threshold <- 10^(-digits)
    if (p < threshold) {
        threshold_str <- paste0("< 0", marks$decimal.mark,
                                strrep("0", digits - 1), "1")
        return(threshold_str)
    }

    fmt_str <- paste0("%.", digits, "f")
    result <- sprintf(fmt_str, p)
    if (marks$decimal.mark != ".") {
        result <- sub(".", marks$decimal.mark, result, fixed = TRUE)
    }
    result
}


#' Validate number_format parameter
#'
#' Checks that a \code{number_format} value is valid before use. Called early
#' in top-level functions to fail fast with a clear error message.
#'
#' @param number_format Value to validate.
#' @return Invisibly returns \code{TRUE} if valid.
#' @keywords internal
validate_number_format <- function(number_format) {
    if (is.null(number_format)) return(invisible(TRUE))

    if (!is.character(number_format)) {
        stop("'number_format' must be a character string or character vector.",
             call. = FALSE)
    }

    if (length(number_format) == 1L) {
        valid_presets <- c("us", "eu", "space", "none")
        if (!number_format %in% valid_presets) {
            stop("Unknown number_format preset: '", number_format,
                 "'. Valid presets are: ",
                 paste(paste0("'", valid_presets, "'"), collapse = ", "),
                 ". Or use a two-element vector c(big.mark, decimal.mark).",
                 call. = FALSE)
        }
    } else if (length(number_format) == 2L) {
        if (number_format[1] == number_format[2] && nchar(number_format[1]) > 0) {
            stop("big.mark and decimal.mark cannot be the same non-empty character.",
                 call. = FALSE)
        }
    } else {
        stop("Custom 'number_format' must be a two-element character vector ",
             "c(big.mark, decimal.mark).", call. = FALSE)
    }

    invisible(TRUE)
}


#' Resolve the CI or range separator
#'
#' Determines the appropriate separator character between two numeric bounds
#' (\emph{e.g.,} CI lower-upper, range min-max) based on whether either bound is
#' negative and the current locale's decimal mark. This avoids ambiguous
#' output like \code{"(-5--3)"} or \code{"(1,2-3,4)"} with European commas.
#'
#' Rules:
#' \enumerate{
#'   \item If either bound is negative, use \code{" to "} to avoid
#'     double-hyphen ambiguity (\emph{e.g.,} \code{"-5 to -3"} not \code{"-5--3"}).
#'   \item If the decimal mark is a comma (EU locale), use \code{"\u2013"}
#'     (en-dash) to avoid confusion between decimal commas and separating
#'     commas (\emph{e.g.,} \code{"1,2\u20133,4"} not \code{"1,2-3,4"}).
#'   \item Otherwise, use a plain hyphen \code{"-"}.
#' }
#'
#' @param lower Numeric value of the lower bound (or minimum).
#' @param upper Numeric value of the upper bound (or maximum).
#' @param marks List with \code{big.mark} and \code{decimal.mark} as returned
#'   by \code{\link{resolve_number_marks}}.
#' @return Character string separator.
#' @keywords internal
resolve_separator <- function(lower, upper, marks) {
    has_negative <- (!is.na(lower) && lower < 0) ||
                    (!is.na(upper) && upper < 0)
    if (has_negative) return(" to ")
    if (marks$decimal.mark == ",") return("\u2013")
    "-"
}


#' Apply locale decimal mark to a sprintf-formatted string
#'
#' Replaces the period decimal mark in a \code{sprintf}-formatted string with
#' the locale-appropriate decimal mark, and fixes negative-zero artefacts.
#'
#' @param x Character string (already formatted with \code{sprintf}).
#' @param marks List with \code{big.mark} and \code{decimal.mark} as returned
#'   by \code{\link{resolve_number_marks}}.
#' @return Character string with corrected decimal marks and no negative zeros.
#' @keywords internal
apply_decimal_mark <- function(x, marks) {
    if (marks$decimal.mark != ".") {
        x <- gsub(".", marks$decimal.mark, x, fixed = TRUE)
    }
    ## Fix negative zero
    dec_esc <- gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", marks$decimal.mark)
    pattern <- paste0("(?<![0-9])-0(", dec_esc, "0+)(?![0-9])")
    gsub(pattern, "0\\1", x, perl = TRUE)
}


#' Format survival median with CI for display
#'
#' Formats a survival median estimate with confidence interval using
#' locale-aware decimal marks and safe CI separators. Used by
#' \code{process_survival} in descriptive tables.
#'
#' @param median Numeric median survival time.
#' @param lower Numeric lower CI bound.
#' @param upper Numeric upper CI bound.
#' @param fmt_str Character string \code{sprintf} format specification.
#' @param marks List with \code{big.mark} and \code{decimal.mark} as returned
#'   by \code{\link{resolve_number_marks}}.
#' @return Character string with formatted "median (lower-upper)".
#' @keywords internal
format_survival_ci <- function(median, lower, upper, fmt_str, marks) {
    sep <- resolve_separator(lower, upper, marks)
    fmt_one <- function(x) {
        result <- sprintf(fmt_str, x)
        apply_decimal_mark(result, marks)
    }
    paste0(fmt_one(median), " (", fmt_one(lower), sep, fmt_one(upper), ")")
}
