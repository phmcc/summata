### * Table formatting

#' Calculate table layout for forest plots
#' 
#' Computes column widths and positions for the table portion of a forest plot.
#' Determines spacing based on content width, font size, and desired table/forest
#' proportion. Returns positions in log-scale units for plot coordinate system.
#' 
#' @param to_show_exp_clean Data.table with formatted display data for the plot.
#' @param show_n Logical whether to include sample size column.
#' @param show_events Logical whether to include events column.
#' @param indent_groups Logical whether groups are indented (affects level column).
#' @param condense_table Logical whether table is condensed (affects level column).
#' @param effect_label Character string describing effect measure type.
#' @param ref_label Character string label for reference categories.
#' @param font_size Numeric font size for width calculations.
#' @param table_width Numeric proportion of total width for table (0-1).
#' @param rangeb Numeric vector of length 2 with plot x-axis range.
#' @param center_padding Numeric additional padding for effect column.
#' @return List with table_width, forest_width, positions, rangeplot_start,
#'   total_width, and effect_abbrev components.
#' @keywords internal
calculate_forest_layout <- function(to_show_exp_clean, show_n, show_events, 
                                    indent_groups, condense_table, 
                                    effect_label, ref_label, font_size, 
                                    table_width = 0.6, rangeb, center_padding) {
    
    char_to_inch <- 0.08 * font_size
    margin <- 0.3  # inches between columns
    
    ## Build list of active columns with their widths
    columns <- list()
    
    ## Always have variable column
    var_width_chars <- max(nchar(to_show_exp_clean$var_display), nchar("Variable"), na.rm = TRUE)
    columns$var <- var_width_chars * char_to_inch
    
    ## Conditionally add level column
    if (!(indent_groups || condense_table)) {
        level_width_chars <- max(nchar(to_show_exp_clean$level), nchar("Group"), na.rm = TRUE)
        columns$level <- level_width_chars * char_to_inch
    }
    
    ## Conditionally add n column
    if (show_n) {
        n_width_chars <- max(nchar(to_show_exp_clean$n_formatted), nchar("____"), na.rm = TRUE)
        columns$n <- n_width_chars * char_to_inch
    }
    
    ## Conditionally add events column
    if (show_events) {
        events_width_chars <- max(nchar(to_show_exp_clean$events_formatted), nchar("___"), na.rm = TRUE)
        columns$events <- events_width_chars * char_to_inch
    }
    
    ## Always have effect column
    ## Create abbreviation from effect_label - handle common cases and custom labels
    ## For multivariable models, use adjusted abbreviations
    effect_abbrev <- if(effect_label == "Odds Ratio") "aOR" 
                     else if(effect_label == "Risk Ratio") "aRR" 
                     else if(effect_label == "Rate Ratio") "aRR"
                     else if(effect_label == "Hazard Ratio") "aHR"
                     else if(effect_label == "Coefficient") "Adj. Coefficient"
                     else if(effect_label == "Exp(Coefficient)") "Exp(Adj. Coef)"
                     else effect_label  ## Use custom label as-is for the header
    
    effect_header <- paste0(effect_abbrev, " (95% CI); p-value")
    
    ## Calculate max effect string length (without expression parsing)
    effect_lengths <- nchar(data.table::fifelse(
                                            is.na(to_show_exp_clean$estimate),
                                            ref_label,
                                            paste0(to_show_exp_clean$effect_formatted, " (",
                                                   to_show_exp_clean$conf_low_formatted, "-",
                                                   to_show_exp_clean$conf_high_formatted, "); p = ",
                                                   to_show_exp_clean$p_formatted)
                                        )) + center_padding
    
    effect_width_chars <- max(effect_lengths, nchar(effect_header), na.rm = TRUE)
    columns$effect <- effect_width_chars * char_to_inch
    
    ## Calculate total table width
    calc_table_width <- sum(unlist(columns)) + length(columns) * margin * 2
    
    ## Calculate forest width based on table_width proportion
    calc_forest_width <- calc_table_width * ((1 - table_width) / table_width)
    
    ## Convert table width to log scale units
    calc_range_width <- diff(rangeb)
    calc_table_width_log_units <- (calc_table_width / calc_forest_width) * calc_range_width
    
    ## Calculate positions in log scale
    ## Start from the extended left edge
    rangeplot_start <- rangeb[1] - calc_table_width_log_units
    
    ## Convert inch measurements to log scale units
    inch_to_log <- calc_range_width / calc_forest_width
    
    positions <- list()
    current_pos <- rangeplot_start
    
    for (name in names(columns)) {
        current_pos <- current_pos + margin * inch_to_log
        positions[[name]] <- current_pos
        current_pos <- current_pos + columns[[name]] * inch_to_log + margin * inch_to_log
    }
    
    return(list(
        calc_table_width = calc_table_width,
        calc_forest_width = calc_forest_width,
        positions = positions,
        rangeplot_start = rangeplot_start,
        total_width = calc_table_width + calc_forest_width,
        effect_abbrev = effect_abbrev
    ))
}


#' Find non-reference row for binary variable condensing
#' 
#' Identifies the non-reference row in a binary categorical variable by checking
#' for NA estimates (reference rows have NA). This is more robust than assuming
#' row position or matching specific strings like "Yes" or "Positive".
#' 
#' @param var_rows Data.table containing rows for a single variable.
#' @param estimate_col Character string naming the estimate column 
#'   (\emph{e.g.,} "estimate", "coef"). Default is "estimate".
#' @return Integer index of the non-reference row within var_rows, or \code{NULL} if
#'   cannot be determined (\emph{e.g.,} no NA estimates found, or multiple non-NA rows).
#' @keywords internal
find_non_reference_row <- function(var_rows, estimate_col = "estimate") {
    if (!estimate_col %in% names(var_rows)) {
        ## Fallback: try common column names
        possible_cols <- c("estimate", "coef", "coefficient", "hr", "or", "rr")
        estimate_col <- intersect(tolower(names(var_rows)), possible_cols)[1]
        if (is.na(estimate_col)) return(NULL)
        ## Match actual case
        estimate_col <- names(var_rows)[tolower(names(var_rows)) == estimate_col][1]
    }
    
    estimates <- var_rows[[estimate_col]]
    
    ## Reference row has NA estimate; non-reference has actual value
    non_ref_idx <- which(!is.na(estimates))
    
    if (length(non_ref_idx) == 1) {
        return(non_ref_idx)
    } else if (length(non_ref_idx) == 0) {
        ## All NA - shouldn't happen for binary, but fallback to row 2
        return(2L)
    } else {
        ## Multiple non-NA values - not a simple binary reference pattern
        ## Return NULL to signal condensing should not occur
        return(NULL)
    }
}


#' Check if category name is a standard reference/negative value
#' 
#' Determines whether a category name represents a standard reference or negative
#' value that indicates absence. Used to suppress redundant category names when
#' condensing binary variables.
#' 
#' @param category Character string with the category name.
#' @param label Optional character string with the variable label. If provided,
#'   checks if category is "No [label]" or similar patterns.
#' @param norm_category Optional pre-normalized category (lowercase, trimmed).
#'   If provided, skips normalization for performance.
#' @param norm_label Optional pre-normalized label (lowercase, trimmed).
#'   If provided, skips normalization for performance.
#' @return Logical indicating whether category is a reference/negative value.
#' @keywords internal
is_reference_category <- function(category, label = NULL, 
                                   norm_category = NULL, norm_label = NULL) {
    if (is.null(category) || is.na(category) || category == "") {
        return(FALSE)
    }
    
    ## Use pre-normalized value if provided, otherwise normalize
    if (is.null(norm_category)) {
        norm_category <- tolower(trimws(category))
    }
    
    ## Fast check against standard reference values using switch for common cases
    if (norm_category == "no" || norm_category == "none" || 
        norm_category == "0" || norm_category == "false" ||
        norm_category == "absent" || norm_category == "negative" ||
        norm_category == "normal" || norm_category == "-") {
        return(TRUE)
    }
    
    ## Check prefix patterns with startsWith (faster than regex for simple prefixes)
    if (startsWith(norm_category, "no ") || 
        startsWith(norm_category, "non-") || 
        startsWith(norm_category, "non ") ||
        startsWith(norm_category, "without ")) {
        return(TRUE)
    }
    
    ## Check suffix patterns with endsWith (faster than regex)
    if (endsWith(norm_category, " absent") || endsWith(norm_category, " negative")) {
        return(TRUE)
    }
    
    ## If label provided, check if category is "No [label]" pattern
    if (!is.null(label) && !is.na(label) && label != "") {
        if (is.null(norm_label)) {
            norm_label <- tolower(trimws(label))
        }
        
        ## "No [label]" exact pattern
        if (norm_category == paste0("no ", norm_label)) {
            return(TRUE)
        }
        
        ## "No [partial label]" pattern - only check if starts with "no "
        if (startsWith(norm_category, "no ")) {
            category_suffix <- substring(norm_category, 4L)  # Faster than sub()
            if (grepl(category_suffix, norm_label, fixed = TRUE) ||
                grepl(norm_label, category_suffix, fixed = TRUE)) {
                return(TRUE)
            }
        }
    }
    
    return(FALSE)
}


#' Check if category name should be suppressed in condensed label
#' 
#' Determines whether a category name should be suppressed when condensing
#' binary variables. Returns \code{TRUE} for standard affirmative values (\emph{e.g.,} "Yes", 
#' "1", "Positive"), standard reference values (\emph{e.g.,} "No", "Absent", "None"),
#' or when the category name essentially matches the variable label 
#' (case-insensitive comparison).
#' 
#' @param category Character string with the category name.
#' @param label Optional character string with the variable label. If provided,
#'   returns \code{TRUE} when category is a case-insensitive match or substring.
#' @param norm_category Optional pre-normalized category (lowercase, trimmed).
#'   If provided, skips normalization for performance.
#' @param norm_label Optional pre-normalized label (lowercase, trimmed).
#'   If provided, skips normalization for performance.
#' @return Logical indicating whether category should be suppressed.
#' @keywords internal
is_affirmative_category <- function(category, label = NULL,
                                     norm_category = NULL, norm_label = NULL) {
    if (is.null(category) || is.na(category) || category == "") {
        return(FALSE)
    }
    
    ## Use pre-normalized value if provided, otherwise normalize once
    if (is.null(norm_category)) {
        norm_category <- tolower(trimws(category))
    }
    
    ## Fast check against standard affirmative values
    if (norm_category == "yes" || norm_category == "1" || 
        norm_category == "true" || norm_category == "present" ||
        norm_category == "positive" || norm_category == "+") {
        return(TRUE)
    }
    
    ## Pre-normalize label if needed (do this before calling is_reference_category)
    if (!is.null(label) && !is.na(label) && label != "" && is.null(norm_label)) {
        norm_label <- tolower(trimws(label))
    }
    
    ## Check if it's a reference category (pass pre-normalized values)
    if (is_reference_category(category, label, norm_category, norm_label)) {
        return(TRUE)
    }
    
    ## Check if category matches or is contained in the label (case-insensitive)
    if (!is.null(norm_label)) {
        ## Check for exact match
        if (norm_category == norm_label) {
            return(TRUE)
        }
        
        ## Check if category is a substantial substring of label
        nc_len <- nchar(norm_category)
        if (nc_len >= 3L && grepl(norm_category, norm_label, fixed = TRUE)) {
            return(TRUE)
        }
        
        ## Check reverse: label is substring of category
        nl_len <- nchar(norm_label)
        if (nl_len >= 3L && grepl(norm_label, norm_category, fixed = TRUE)) {
            return(TRUE)
        }
    }
    
    return(FALSE)
}


#' Check if a binary variable should be condensed without category suffix
#' 
#' Uses a greedy/liberal approach to determine if a binary variable's condensed
#' display should omit the category name. Returns \code{TRUE} if EITHER level of the
#' binary variable is a standard reference/affirmative value, OR if either level
#' matches/contains the variable label.
#' 
#' This function is designed for binary (2-level) categorical variables where
#' one level is a reference and one is the "event" or "condition" level.
#' 
#' @param ref_category Character string with the reference category name
#'   (the level with NA estimate).
#' @param non_ref_category Character string with the non-reference category name
#'   (the level with the actual estimate).
#' @param label Optional character string with the variable label. Used for
#'   intelligent matching (\emph{e.g.,} "30-Day Readmission" label with 
#'   "30-day readmission" / "No 30-day readmission" levels).
#' @return Logical indicating whether the binary variable should be condensed
#'   without appending the category name.
#' @keywords internal
should_condense_binary <- function(ref_category, non_ref_category, label = NULL) {
    ## Pre-normalize all strings once to avoid redundant normalization
    norm_ref <- if (!is.null(ref_category) && !is.na(ref_category) && ref_category != "") {
                    tolower(trimws(ref_category))
                } else {
                    NULL
                }
    
    norm_non_ref <- if (!is.null(non_ref_category) && !is.na(non_ref_category) && non_ref_category != "") {
                        tolower(trimws(non_ref_category))
                    } else {
                        NULL
                    }
    
    norm_label <- if (!is.null(label) && !is.na(label) && label != "") {
                      tolower(trimws(label))
                  } else {
                      NULL
                  }
    
    ## Check reference category first (more likely to be "No", "0", etc.)
    if (!is.null(norm_ref)) {
        if (is_reference_category(ref_category, label, norm_ref, norm_label)) {
            return(TRUE)
        }
    }
    
    ## Check non-reference category
    if (!is.null(norm_non_ref)) {
        if (is_affirmative_category(non_ref_category, label, norm_non_ref, norm_label)) {
            return(TRUE)
        }
    }
    
    return(FALSE)
}


### * Text formatting

#' Format numeric value with fixed decimal places
#' 
#' Formats a numeric value to a specified number of decimal places, removing
#' leading/trailing whitespace and fixing negative zero display (\emph{e.g.,} "-0.00"
#' becomes "0.00"). When \code{marks} is supplied, applies locale-appropriate
#' decimal mark substitution.
#' 
#' @param x Numeric value to format.
#' @param digits Integer number of decimal places.
#' @param marks Optional list with \code{big.mark} and \code{decimal.mark} as
#'   returned by \code{\link{resolve_number_marks}}.
#' @return Character string with formatted value.
#' @keywords internal
format_number <- function(x, digits, marks = NULL) {
    result <- trimws(format(round(x, digits), nsmall = digits))
    if (!is.null(marks)) {
        return(apply_decimal_mark(result, marks))
    }
    ## Fallback: fix negative zero only (US locale)
    gsub("(?<![0-9])-0(\\.0+)(?![0-9])", "0\\1", result, perl = TRUE)
}


#' Determine CI separator for forest plot text annotations
#'
#' Returns the appropriate separator string between CI lower and upper bounds
#' in forest plot annotations. Considers whether values may be negative and
#' the current locale's decimal mark.
#'
#' @param has_negatives Logical indicating whether any CI bounds are negative.
#' @param marks Optional list with \code{big.mark} and \code{decimal.mark} as
#'   returned by \code{\link{resolve_number_marks}}.
#' @return Character string separator.
#' @keywords internal
forest_ci_separator <- function(has_negatives, marks = NULL) {
    if (has_negatives) return(", ")
    if (!is.null(marks) && marks$decimal.mark == ",") return("\u2013")
    "-"
}


#' Format a p-value for forest plot annotations
#'
#' Returns a formatted p-value string suitable for forest plot display,
#' using locale-aware decimal marks when \code{marks} is provided.
#'
#' @param p Numeric p-value.
#' @param p_digits Integer decimal places.
#' @param marks Optional list with \code{big.mark} and \code{decimal.mark} as
#'   returned by \code{\link{resolve_number_marks}}.
#' @return Character string.
#' @keywords internal
format_p_forest <- function(p, p_digits, marks = NULL) {
    threshold <- 10^(-p_digits)
    if (p < threshold) {
        if (!is.null(marks)) {
            return(paste0("< 0", marks$decimal.mark,
                          strrep("0", p_digits - 1), "1"))
        }
        return(paste0("< ", format(threshold, scientific = FALSE)))
    }
    format_number(p, p_digits, marks)
}


### * CI caching

#' Retrieve confidence intervals with cache support
#'
#' Returns confidence intervals for a fitted model, using a cached result
#' when available from upstream \code{fit()} or \code{m2dt()} calls to avoid
#' redundant computation. This is particularly beneficial for GLM and negative
#' binomial models where \code{confint()} performs profile likelihood profiling,
#' which can be expensive for models with many parameters.  Otherwise the
#' intervals are computed by \code{coef_confint()}, as in \code{m2dt()}, by the
#' method that \code{getOption("summata.conf_method")} selects. They are
#' computed for the rows of the coefficient table, which omits aliased
#' coefficients, except the intercept, which the forest plots do not show; its
#' row holds \code{NA}.
#'
#' @param model Fitted model object. If the model carries a
#'   \code{"cached_confint"} attribute (set by \code{fit()} during table
#'   generation) and the cached confidence level matches \code{conf_level},
#'   the cached result is returned directly.
#' @param conf_level Numeric confidence level. Must match the cached level
#'   for the cache to be used.
#' @return A matrix with one row per row of \code{summary(model)$coefficients}
#'   and two columns (lower and upper bounds).
#' @keywords internal
get_cached_confint <- function(model, conf_level = 0.95) {

    ## Check for cached result from upstream fit()/m2dt()
    cached <- attr(model, "cached_confint")
    if (!is.null(cached) &&
        isTRUE(attr(model, "cached_confint_level") == conf_level)) {
        return(cached)
    }

    ## The callers pair these rows with those of the coefficient table by
    ## position, so the rows must be the same
    terms <- rownames(summary(model)$coefficients)
    parm <- setdiff(terms, "(Intercept)")
    conf_int <- matrix(NA_real_, nrow = length(terms), ncol = 2L,
                       dimnames = list(terms, NULL))
    if (length(parm) > 0L) {
        computed <- coef_confint(model, parm, conf_level,
                                 getOption("summata.conf_method", "profile"))
        conf_int[parm, ] <- computed
        colnames(conf_int) <- colnames(computed)
    }
    conf_int
}


#' Quantile for confidence intervals recomputed at another level
#'
#' \code{uniforest()} and \code{multiforest()} draw the intervals of their
#' tables, and recompute them from the coefficients and their standard
#' errors when another confidence level is requested. The tables hold the
#' exact \emph{t} intervals of linear models (unless Wald intervals were
#' requested) and record their residual degrees of freedom; the other
#' intervals are recomputed as Wald intervals.
#'
#' @param conf_level Numeric confidence level.
#' @param df_residual Numeric vector of residual degrees of freedom, one per
#'   row, missing where the table's intervals are not \emph{t} intervals, or
#'   \code{NULL} when the table records none.
#' @return Numeric vector of quantiles, one per row, or a single normal
#'   quantile when \code{df_residual} is \code{NULL}.
#' @keywords internal
interval_quantile <- function(conf_level, df_residual = NULL) {
    p <- (1 + conf_level) / 2
    if (is.null(df_residual)) {
        return(stats::qnorm(p))
    }
    ifelse(is.na(df_residual), stats::qnorm(p), stats::qt(p, df_residual))
}


### * Model terms

#' Classes of the fixed-effect variables of a mixed model
#'
#' The counterpart, for \pkg{lme4} models, of the \code{"dataClasses"}
#' attribute of the terms of other models: the class of each variable of
#' the fixed-effect terms, as the model frame records it. A variable is an
#' expression of the formula, such as \code{log(age)}, rather than a name
#' within one, so that the variables are named as the coefficients of
#' their terms begin.
#'
#' @param model Fitted \pkg{lme4} model.
#' @return Named character vector of the classes of the fixed-effect
#'   variables, without the response.
#' @keywords internal
fixed_effect_classes <- function(model) {
    fixed_terms <- stats::terms(model, fixed.only = TRUE)
    variables <- as.list(attr(fixed_terms, "variables"))[-1L]
    response <- attr(fixed_terms, "response")
    if (length(response) == 1L && response > 0L) {
        variables <- variables[-response]
    }

    ## Named as model.frame() names the columns of the model frame
    var_names <- vapply(variables, function(v) {
        paste(deparse(v, width.cutoff = 500L,
                      backtick = !is.symbol(v) && is.language(v)),
              collapse = " ")
    }, character(1))

    frame <- model@frame
    frame_classes <- attr(attr(frame, "terms"), "dataClasses")
    classes <- vapply(var_names, function(v) {
        if (v %in% names(frame_classes)) {
            frame_classes[[v]]
        } else if (v %in% names(frame)) {
            stats::.MFclass(frame[[v]])
        } else {
            "numeric"
        }
    }, character(1))
    stats::setNames(classes, var_names)
}


### * Aliased coefficients

#' Names of the coefficients that a model could not estimate
#'
#' Returns the names of the aliased coefficients of a model, those whose
#' columns of the model matrix are linear combinations of other columns.
#' \code{stats::coef()} reports them as \code{NA} for \code{lm}, \code{glm},
#' negative binomial and \code{coxph} models; \pkg{lme4} drops their columns
#' from the model matrix before fitting and records them in its
#' \code{"col.dropped"} attribute. The backticks that enclose non-syntactic
#' variable names in coefficient names are removed.
#'
#' @param model Fitted model object.
#' @return Character vector of coefficient names, empty when every
#'   coefficient was estimated.
#' @keywords internal
aliased_coefficients <- function(model) {
    aliased <- if (inherits(model, "merMod")) {
                   names(attr(lme4::getME(model, "X"), "col.dropped"))
               } else {
                   estimates <- tryCatch(stats::coef(model),
                                         error = function(e) NULL)
                   names(estimates)[is.na(estimates)]
               }
    gsub("`", "", as.character(aliased), fixed = TRUE)
}


#' Level that serves as the reference of a factor whose first level is empty
#'
#' Cox models keep factor levels without observations, where \code{lm()}
#' and \code{glm()} drop them. When the first level, the reference of the
#' treatment contrasts, has no observations, the indicators of the other
#' levels sum to one, which the baseline hazard absorbs. One of them is then
#' aliased, its coefficient fixed at zero, and the others are estimated
#' against it: that level, not the empty first level, is the reference of
#' the estimates. A level may also be aliased because other terms of the
#' model determine it, as when a covariate is the indicator of that level.
#' Its effect is then carried by those terms rather than fixed at zero, so
#' that it is not the reference. Where several levels are aliased, the
#' levels that other terms determine (see \code{determined_levels()}) are
#' set aside.
#'
#' @param levels Character vector of the levels of a factor, in order.
#' @param counts Numeric vector of the number of observations at each level.
#' @param aliased Logical vector indicating, for each level, whether its
#'   coefficient is aliased.
#' @param determined Optional logical vector indicating, for each level,
#'   whether other terms of the model determine its indicator.
#' @return The level that serves as the reference, or \code{NA} when the
#'   first level has observations or no single aliased level with
#'   observations can serve.
#' @keywords internal
effective_reference_level <- function(levels, counts, aliased, determined = NULL) {
    if (length(levels) < 2L || is.na(counts[1L]) || counts[1L] > 0) {
        return(NA_character_)
    }
    candidate <- aliased & !is.na(counts) & counts > 0
    if (!is.null(determined)) {
        candidate <- candidate & !determined
    }
    candidates <- levels[-1L][candidate[-1L]]
    if (length(candidates) == 1L) as.character(candidates) else NA_character_
}


#' Levels of a factor that the other terms of a model determine
#'
#' Tests, for each level of a factor, whether its indicator lies in the
#' span of the columns of the model's other terms and of the constant, or
#' of the strata of a stratified model, whose baseline hazards absorb
#' constants within strata. The columns are centered within strata (or
#' about their means), which takes the place of the constant or the strata,
#' and the span is read from the model matrix, whose rows are the
#' observations analyzed. The columns of the interactions that the factor
#' enters are set aside with its own, as they are built from its levels.
#'
#' @param model Fitted \code{coxph} or \code{clogit} model.
#' @param var Character string naming the factor.
#' @param levels Character vector of its levels.
#' @return Logical vector, one per level, \code{TRUE} where the other terms
#'   determine the level, or \code{NULL} when the model matrix cannot be
#'   rebuilt for the observations analyzed.
#' @keywords internal
determined_levels <- function(model, var, levels) {
    model_terms <- tryCatch(stats::terms(model), error = function(e) NULL)
    X <- tryCatch(stats::model.matrix(model), error = function(e) NULL)
    if (is.null(model_terms) || is.null(X)) {
        return(NULL)
    }

    ## The model matrix is rebuilt from the data of the model's call, which
    ## must still hold the observations analyzed
    n_used <- get_model_nobs(model, class(model)[1L])
    if (!is.na(n_used) && nrow(X) != n_used) {
        return(NULL)
    }

    ## The terms that the factor enters: its main effect and its interactions
    factors <- attr(model_terms, "factors")
    var_row <- match(var, gsub("`", "", rownames(factors), fixed = TRUE))
    if (is.na(var_row)) {
        return(NULL)
    }
    own <- attr(X, "assign") %in% which(factors[var_row, ] > 0)

    ## Centered within strata, or about the mean of a model without strata
    strata_index <- attr(model_terms, "specials")$strata
    groups <- if (length(strata_index) > 0L) {
                  frame <- tryCatch(stats::model.frame(model), error = function(e) NULL)
                  if (is.null(frame) || nrow(frame) != nrow(X)) {
                      return(NULL)
                  }
                  interaction(frame[strata_index], drop = TRUE)
              } else {
                  rep(1L, nrow(X))
              }
    center <- function(M) {
        M - apply(M, 2L, function(column) stats::ave(column, groups))
    }

    base <- center(X[, !own, drop = FALSE])
    base_rank <- if (ncol(base) > 0L) qr(base)$rank else 0L

    columns <- match(paste0(var, levels), gsub("`", "", colnames(X), fixed = TRUE))
    vapply(columns, function(column) {
        if (is.na(column)) {
            return(FALSE)
        }
        indicator <- center(X[, column, drop = FALSE])
        qr(cbind(base, indicator))$rank == base_rank
    }, logical(1))
}


#' Find the reference level of a factor whose first level is empty
#'
#' Applies \code{effective_reference_level()}, setting aside the levels that
#' other terms determine where several aliased levels could serve. Some
#' aliasing leaves no single reference: the levels may each be estimated
#' against a different aliased level. No reference is then named, and a
#' warning says so, as when the levels that other terms determine cannot be
#' found because the model matrix cannot be rebuilt.
#'
#' @param model Fitted model, or \code{NULL}.
#' @param var Character string naming the factor.
#' @param levels,counts,aliased As for \code{effective_reference_level()}.
#' @return List with \code{reference}, the reference level or \code{NA},
#'   and \code{ambiguous}, \code{TRUE} when the first level is empty but no
#'   single level serves as the reference.
#' @keywords internal
factor_reference <- function(model, var, levels, counts, aliased) {
    reference <- effective_reference_level(levels, counts, aliased)
    first_empty <- length(levels) >= 2L && !is.na(counts[1L]) && counts[1L] == 0
    n_candidates <- sum(aliased[-1L] & !is.na(counts[-1L]) & counts[-1L] > 0)
    if (!is.na(reference) || !first_empty || n_candidates < 2L) {
        return(list(reference = reference, ambiguous = FALSE))
    }

    determined <- if (!is.null(model)) determined_levels(model, var, levels)
    if (is.null(determined)) {
        warning("The reference level of '", var, "' could not be determined: ",
                "its first level has no observations and several of its levels ",
                "are aliased. No reference level is shown.", call. = FALSE)
        return(list(reference = NA_character_, ambiguous = TRUE))
    }
    reference <- effective_reference_level(levels, counts, aliased, determined)
    if (is.na(reference)) {
        warning("The levels of '", var, "' are not estimated against a single ",
                "reference level: its first level has no observations and ",
                "several of its levels are aliased with other terms of the ",
                "model. No reference level is shown.", call. = FALSE)
    }
    list(reference = reference, ambiguous = is.na(reference))
}


#' Omit the rows of aliased coefficients from a forest plot
#'
#' Removes the rows that display aliased coefficients, which have no
#' estimate and would otherwise be drawn as reference categories, and the
#' rows of any factor none of whose levels remains estimated, which would
#' otherwise be shown as a reference category compared with nothing. The
#' regression tables omit both. A factor whose first level has no
#' observations is estimated against an aliased level instead (see
#' \code{effective_reference_level()}). That level is kept, first among
#' the rows of the factor, as its reference category, and the empty first
#' level is removed, as in the regression tables.
#'
#' @param terms_df Data table with one row per displayed term, with columns
#'   \code{var}, \code{level}, \code{N} and \code{pos} and, for interaction
#'   rows, \code{.int_term}.
#' @param aliased Character vector of aliased coefficient names, as returned
#'   by \code{aliased_coefficients()}.
#' @param coef_terms Character vector of the names of the model's
#'   coefficients.
#' @param model Optional fitted model, from which the reference of a factor
#'   with several aliased levels is found (see \code{factor_reference()}).
#' @return \code{terms_df} without the omitted rows.
#' @keywords internal
drop_aliased_rows <- function(terms_df, aliased, coef_terms, model = NULL) {
    if (length(aliased) == 0L || nrow(terms_df) == 0L) {
        return(terms_df)
    }

    ## The coefficient that each row displays, keyed as the forest plots key
    ## the merge with the coefficient table: the stored term of an interaction
    ## row, and otherwise the variable name followed by its level
    key <- data.table::fifelse(terms_df$level == "-", terms_df$var,
                               paste0(terms_df$var, terms_df$level))
    if (".int_term" %in% names(terms_df)) {
        key <- data.table::fifelse(is.na(terms_df$.int_term), key,
                                   terms_df$.int_term)
    }
    key <- gsub("`", "", key, fixed = TRUE)
    aliased_row <- key %in% aliased

    ## Factors whose first level is empty are estimated against an aliased
    ## level, which is kept as their reference in place of the empty level.
    ## Where no single level serves, the empty level is omitted all the same.
    affected <- unique(terms_df$var[aliased_row & terms_df$level != "-"])
    reference_row <- rep(FALSE, nrow(terms_df))
    empty_row <- rep(FALSE, nrow(terms_df))
    for (v in affected) {
        rows <- which(terms_df$var == v & terms_df$level != "-")
        rows <- rows[order(terms_df$pos[rows])]
        found <- factor_reference(model, v, terms_df$level[rows],
                                  terms_df$N[rows], aliased_row[rows])
        if (!is.na(found$reference)) {
            reference_row[rows[terms_df$level[rows] == found$reference]] <- TRUE
        }
        ## The empty first level is not the reference of the estimates
        if (!is.na(found$reference) || found$ambiguous) {
            empty_row[rows[1L]] <- TRUE
        }
    }

    ## Only factors that lost a level are examined, so that a factor whose
    ## coefficients are named otherwise (by non-treatment contrasts, for
    ## example) is left as it is
    estimated <- setdiff(gsub("`", "", coef_terms, fixed = TRUE), aliased)
    unestimated <- affected[vapply(affected, function(v) {
        !any(key[terms_df$var == v & !aliased_row] %in% estimated)
    }, logical(1))]

    omit <- (aliased_row & !reference_row) | empty_row |
        terms_df$var %in% unestimated
    result <- terms_df[!omit]

    ## The reference category is shown first, as for other factors
    moved <- which(reference_row[!omit])
    if (length(moved) > 0L) {
        order_rows <- seq_len(nrow(result))
        for (i in moved) {
            rows <- which(result$var == result$var[i])
            order_rows[rows] <- c(i, setdiff(rows, i))
        }
        result <- result[order_rows]
        for (v in unique(result$var[moved])) {
            result[var == v, pos := seq_len(.N)]
        }
    }

    result
}


### * Font selection

#' Detect available sans-serif font for plots
#' 
#' Checks for commonly available sans-serif fonts in order of preference
#' (Helvetica, Arial, Helvetica Neue) and returns the first available one.
#' Falls back to "sans" if none are found or if systemfonts is unavailable.
#' When ragg is being used as the graphics device (detected via options or
#' knitr settings), font detection works in non-interactive sessions since
#' ragg handles font rendering independently of the R graphics system.
#' 
#' @return Character string with the font family name to use.
#' @keywords internal
detect_plot_font <- function() {
    
    ## Check if ragg is being used - if so, can detect fonts even in
    ## non-interactive sessions since ragg handles fonts via systemfonts.
    ## knitr is only suggested; its chunk options matter only while it is
    ## running, when its namespace is loaded.
    knitr_dev <- if (isNamespaceLoaded("knitr")) knitr::opts_chunk$get("dev") else NULL
    ragg_in_use <- isTRUE(getOption("summata.use_ragg")) ||
        identical(knitr_dev, "ragg_png") ||
        identical(knitr_dev, "agg_png")
    
    ## Only fall back to "sans" in non-interactive sessions when ragg is not in use
    if (!interactive() && !ragg_in_use) {
        return("sans")
    }
    
    if (requireNamespace("systemfonts", quietly = TRUE)) {
        available_fonts <- systemfonts::system_fonts()$family
        ## Check fonts in order of preference
        if ("Helvetica" %in% available_fonts) return("Helvetica")
        if ("Nimbus Sans L" %in% available_fonts) return("Nimbus Sans L")
        if ("Nimbus Sans" %in% available_fonts) return("Nimbus Sans")
        if ("Liberation Sans" %in% available_fonts) return("Liberation Sans")
        if ("Arial" %in% available_fonts) return("Arial")
        if ("Helvetica Neue" %in% available_fonts) return("Helvetica Neue")
    }
    ## Fallback to generic sans-serif
    return("sans")
}


### * Unit conversion

#' Convert between units
#' 
#' Converts measurements between different unit systems commonly used in
#' graphics (inches, centimeters, millimeters, pixels, points).
#' 
#' @param value Numeric value to convert.
#' @param from Character string specifying source unit ("in", "cm", "mm", "px", "pt").
#' @param to Character string specifying target unit ("in", "cm", "mm", "px", "pt").
#' @param dpi Integer dots per inch for pixel conversions (default 96).
#' @return Numeric value in target units.
#' @keywords internal
convert_units <- function(value, from = "in", to = "in", dpi = 96) {
    to_inches <- c(
        "in" = 1,
        "cm" = 1/2.54,    # 1 cm = 1/2.54 inches
        "mm" = 1/25.4,    # 1 mm = 1/25.4 inches
        "px" = 1/dpi,     # 1 px = 1/dpi inches
        "pt" = 1/72       # 1 pt = 1/72 inches
    )
    
    value_in_inches <- value * to_inches[[from]]
    value_in_target <- value_in_inches / to_inches[[to]]
    
    return(value_in_target)
}
