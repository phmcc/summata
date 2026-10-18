### * Table formatting

#' Format model results for publication-ready display
#' 
#' Transforms raw model coefficient data into a formatted table suitable for
#' publication. Handles effect measure formatting (OR, HR, RR, Estimate),
#' confidence intervals, \emph{p}-values, sample sizes, and variable labels.
#' Supports interaction terms and mixed-effects models.
#' 
#' @param data Data.table containing raw model results with coefficient columns.
#' @param effect_col Optional character string specifying the effect column name.
#'   If \code{NULL}, auto-detects from OR, HR, RR, or Estimate columns.
#' @param digits Integer number of decimal places for effect estimates.
#' @param p_digits Integer number of decimal places for \emph{p}-values.
#' @param labels Optional named character vector mapping variable names to display labels.
#'   Supports automatic labeling of interaction terms.
#' @param show_n Logical whether to include sample size column.
#' @param show_events Logical whether to include events column (ignored for linear models).
#' @param reference_label Character string to display for reference categories.
#' @param exponentiate Optional logical to force exponentiated (TRUE) or raw (FALSE)
#'   coefficient display. If \code{NULL}, uses existing columns.
#' @return Formatted data.table with publication-ready columns.
#' @keywords internal
format_model_table <- function(data, 
                               effect_col = NULL,
                               digits = 2, 
                               p_digits = 3, 
                               labels = NULL,
                               show_n = TRUE,
                               show_events = TRUE,
                               reference_label = "reference",
                               exponentiate = NULL,
                               conf_level = 0.95,
                               marks = NULL) {
    
    ## Determine the required columns, to avoid copying the whole table
    ## Start with columns that will be in output
    needed_cols <- c("variable", "group", "n", "n_group", "events", "events_group",
                     "p_value", "ci_lower", "ci_upper", "reference", "model_type",
                     "model_scope", "OR", "HR", "RR", "Coefficient",
                     "coef", "coef_lower", "coef_upper", "exp_coef", "exp_lower", "exp_upper")
    
    ## Only copy columns that exist and are needed
    available_cols <- intersect(needed_cols, names(data))
    result <- data[, ..available_cols]
    
    ## Make a shallow copy to avoid modifying original
    result <- data.table::copy(result)

    ## Disallow "Events" column for models where events don't make sense
    ## This includes linear models and GLMs with continuous distributions
    if ("model_type" %in% names(result)) {
        model_type <- result$model_type[1]
        ## Linear, Gamma, Gaussian, and Inverse Gaussian models don't have meaningful events
        continuous_patterns <- c("Linear", "Gamma", "Gaussian", "Inverse")
        if (any(vapply(continuous_patterns, function(p) grepl(p, model_type, ignore.case = TRUE), logical(1)))) {
            show_events <- FALSE
        }
    }
    
    ## Also check if all events are NA - if so, don't show the column
    if (show_events && "events" %in% names(result)) {
        if (all(is.na(result$events))) {
            show_events <- FALSE
        }
    }
    
    ## Standardize column names
    if ("variable" %in% names(result)) {
        data.table::setnames(result, "variable", "Variable")
    }
    if ("group" %in% names(result)) {
        data.table::setnames(result, "group", "Group")
    }

    ## Handle the exponentiate parameter to choose which columns to use
    if (!is.null(exponentiate)) {
        if (exponentiate && "exp_coef" %in% names(result)) {
            ## Check model type
            if ("OR" %in% names(result)) {
                effect_col <- "OR"
            } else if ("HR" %in% names(result)) {
                effect_col <- "HR"
            } else if ("RR" %in% names(result)) {
                effect_col <- "RR"
            } else {
                ## Generic model - create OR/RR based on model type
                model_type <- result$model_type[1]
                if (grepl("Logistic", model_type)) {
                    result[, `:=`(
                        OR = exp_coef,
                        ci_lower = exp_lower,
                        ci_upper = exp_upper
                    )]
                    effect_col <- "OR"
                } else if (grepl("Poisson", model_type)) {
                    result[, `:=`(
                        RR = exp_coef,
                        ci_lower = exp_lower,
                        ci_upper = exp_upper
                    )]
                    effect_col <- "RR"
                } else {
                    result[, `:=`(
                        Coefficient = exp_coef,
                        ci_lower = exp_lower,
                        ci_upper = exp_upper
                    )]
                    effect_col <- "Coefficient"
                }
            }
        } else if (!exponentiate && "coef" %in% names(result)) {
            result[, `:=`(
                Coefficient = coef,
                ci_lower = coef_lower,
                ci_upper = coef_upper
            )]
            effect_col <- "Coefficient"
        }
    }
    
    ## Auto-detect effect column if not specified
    if (is.null(effect_col)) {
        effect_col <- intersect(c("OR", "HR", "RR", "Coefficient"), names(result))[1]
        if (length(effect_col) == 0) {
            stop("No effect measure column found (OR, HR, RR, or Coefficient)")
        }
    }
    
    ## Apply variable labels if provided
    if (!is.null(labels) && "Variable" %in% names(result) && length(labels) > 0) {
        
        ## Interaction terms join their components with ":". They are
        ## identified before the main effects are labeled, so that a label
        ## containing a colon is not taken for an interaction, and apart
        ## from the main-effect terms that m2dt() records, so that neither is
        ## the name of a variable that contains one.
        interaction_mask <- grepl(":", result$Variable, fixed = TRUE) &
            !result$Variable %in% attr(data, "main_terms")
        interaction_vars <- result$Variable[interaction_mask]

        ## The variables of the table and the levels of its factors, read from
        ## the rows of the main effects, by which the interaction terms are
        ## split into their components
        table_vars <- unique(result$Variable[!interaction_mask])
        table_levels <- lapply(split(result$Group[!interaction_mask],
                                     result$Variable[!interaction_mask]),
                               function(groups) setdiff(groups, c("", "-")))

        ## Create lookup table for main effects
        label_dt <- data.table::data.table(
                                    var_orig = names(labels),
                                    var_new = unname(unlist(labels))
                                )
        
        ## Update main effect variable names using merge
        result[label_dt, Variable := i.var_new, on = .(Variable = var_orig)]
        
        if (any(interaction_mask)) {
            ## Vectorized interaction labeling
            labeled_interactions <- vapply(interaction_vars, function(original_var) {
                ## Check for direct custom label first
                if (original_var %in% names(labels)) {
                    return(labels[[original_var]])
                }
                
                ## Match each component to the variable it names, the longest
                ## name first and preferring a variable of which the remainder
                ## is a level (see split_interaction_term()). Matching the
                ## labeled variables alone took the component of an unlabeled
                ## sexuality for a level of a labeled sex, and splitting at
                ## every colon divided a level such as "12:00".
                components <- split_interaction_term(original_var,
                                                     union(table_vars, names(labels)),
                                                     table_levels)
                
                labeled_parts <- vapply(components, function(component) {
                    ## Components of unlabeled variables keep their names
                    if (is.na(component$var) || !component$var %in% names(labels)) {
                        return(component$part)
                    }
                    if (is.na(component$level)) {
                        labels[[component$var]]
                    } else {
                        paste0(labels[[component$var]], " (", component$level, ")")
                    }
                }, character(1))
                
                paste(labeled_parts, collapse = " \u00d7 ")
            }, character(1))
            
            result$Variable[interaction_mask] <- labeled_interactions
        }
    }
    
    ## Clean up Group display (handle empty groups for continuous vars)
    if ("Group" %in% names(result)) {
        result[Group == "", Group := "-"]
    }

    ## Handle n column
    if ("n" %in% names(result)) {
        has_n_group <- "n_group" %in% names(result)
        
        if (has_n_group) {
            ## Use n_group where available, fallback to n
            n_vals <- data.table::fifelse(!is.na(result$n_group), result$n_group, result$n)
        } else {
            n_vals <- result$n
        }
        
        ## Vectorized formatting
        result[, n := as.character(n_vals)]
        if (!is.null(marks)) {
            result[n_vals >= 1000, n := vapply(n_vals[n_vals >= 1000], format_count, character(1), marks = marks)]
        } else {
            result[n_vals >= 1000, n := format(n_vals[n_vals >= 1000], big.mark = ",")]
        }
        result[is.na(n_vals), n := NA_character_]
    }

    ## Handle events column
    if ("events" %in% names(result)) {
        has_events_group <- "events_group" %in% names(result)
        
        if (has_events_group) {
            ## Use events_group where available, fallback to events
            ## Ensure both are numeric for comparison
            events_grp <- as.numeric(result$events_group)
            events_main <- as.numeric(result$events)
            event_vals <- data.table::fifelse(!is.na(events_grp), events_grp, events_main)
        } else {
            event_vals <- as.numeric(result$events)
        }
        
        ## Vectorized formatting
        evt_str <- as.character(round(event_vals))
        big_idx <- which(!is.na(event_vals) & event_vals >= 1000)
        if (length(big_idx) > 0) {
            if (!is.null(marks)) {
                evt_str[big_idx] <- vapply(round(event_vals[big_idx]), format_count, character(1), marks = marks)
            } else {
                evt_str[big_idx] <- format(round(event_vals[big_idx]), big.mark = ",", scientific = FALSE)
            }
        }
        evt_str[is.na(event_vals)] <- NA_character_
        result[, events := evt_str]
    }
    
    ## Eliminate repeated variable names
    if ("Variable" %in% names(result) && nrow(result) > 1) {
        vars <- result$Variable
        result[, Variable := data.table::fifelse(
                                             c(TRUE, vars[-length(vars)] != vars[-1]),
                                             Variable,
                                             ""
                                         )]
    }
    
    ## Create effect column label based on model scope
    model_scope <- if ("model_scope" %in% names(result)) result$model_scope[1] else "Effect"
    
    ## Format confidence level for display (e.g., 0.95 -> "95%", 0.90 -> "90%")
    ci_pct <- round(conf_level * 100)
    ci_label <- paste0(ci_pct, "% CI")
    
    ## Create appropriate label
    if (model_scope == "Univariable") {
        ## Univariable: OR, HR, RR, Coefficient
        effect_label <- paste0(effect_col, " (", ci_label, ")")
    } else if (model_scope == "Multivariable") {
        ## Multivariable: aOR, aHR, aRR, Adj. Coefficient
        adjusted_col <- switch(effect_col,
                               "OR" = "aOR",
                               "HR" = "aHR", 
                               "RR" = "aRR",
                               "Coefficient" = "Adj. Coefficient",
                               paste0("Adj. ", effect_col)
                               )
        effect_label <- paste0(adjusted_col, " (", ci_label, ")")
    } else {
        effect_label <- paste0(effect_col, " (", ci_label, ")")
    }
    
    ## Format effect sizes with CI
    if ("ci_lower" %in% names(result) && "ci_upper" %in% names(result)) {
        is_reference <- if ("reference" %in% names(result)) {
                            !is.na(result$reference) & result$reference == reference_label
                        } else {
                            rep(FALSE, nrow(result))
                        }
        
        effect_vals <- result[[effect_col]]
        ci_lower_vals <- result$ci_lower
        ci_upper_vals <- result$ci_upper
        
        ## Determine CI separator based on locale and values
        if (!is.null(marks)) {
            any_negative <- any(ci_lower_vals < 0 | ci_upper_vals < 0, na.rm = TRUE)
            ci_sep <- if (any_negative) " to " else if (marks$decimal.mark == ",") "\u2013" else "-"
        } else {
            ci_sep <- if (effect_col %in% c("OR", "HR", "RR")) "-" else ", "
        }
        
        ## Pre-compute format strings for individual components
        fmt_single <- paste0("%.", digits, "f")
        
        ## Vectorized sprintf for components
        eff_fmt <- sprintf(fmt_single, effect_vals)
        low_fmt <- sprintf(fmt_single, ci_lower_vals)
        hi_fmt  <- sprintf(fmt_single, ci_upper_vals)
        
        ## Apply locale decimal mark and fix negative zeros
        if (!is.null(marks)) {
            eff_fmt <- apply_decimal_mark(eff_fmt, marks)
            low_fmt <- apply_decimal_mark(low_fmt, marks)
            hi_fmt  <- apply_decimal_mark(hi_fmt, marks)
        } else {
            eff_fmt <- fix_negative_zero(eff_fmt)
            low_fmt <- fix_negative_zero(low_fmt)
            hi_fmt  <- fix_negative_zero(hi_fmt)
        }
        
        ## Assemble: "effect (lower-upper)" or "effect (lower, upper)"
        ## For non-exponentiated models without marks, use comma+space as separator
        if (is.null(marks) && !effect_col %in% c("OR", "HR", "RR")) {
            formatted_effects <- paste0(eff_fmt, " (", low_fmt, ", ", hi_fmt, ")")
        } else {
            formatted_effects <- paste0(eff_fmt, " (", low_fmt, ci_sep, hi_fmt, ")")
        }
        formatted_effects[is.na(effect_vals)] <- ""
        formatted_effects[is_reference] <- reference_label
        
        result[, (effect_label) := formatted_effects]
    }
    
    ## Format p-values
    if ("p_value" %in% names(result)) {
        result[, `p-value` := format_pvalues_fit(p_value, p_digits, marks)]
        
        if ("reference" %in% names(result)) {
            result[!is.na(reference) & reference == reference_label, `p-value` := "-"]
        }
    }

    ## Select columns for final output
    display_cols <- c(
        if ("Variable" %in% names(result)) "Variable",
        if ("Group" %in% names(result)) "Group",
        if (show_n && "n" %in% names(result)) "n",
        if (show_events && "events" %in% names(result)) "events",
        if (effect_label %in% names(result)) effect_label,
        if ("p-value" %in% names(result)) "p-value"
    )
    
    formatted <- result[, ..display_cols]

    if ("events" %in% names(formatted)) {
        data.table::setnames(formatted, "events", "Events")
    }

    return(formatted)
}


### * Ancillary formatters

#' Fix negative zero in formatted strings
#' 
#' Corrects floating-point rounding artifacts that produce "-0.00" or similar
#' negative zero strings. Works on character vectors, replacing patterns like
#' "-0.00", "-0.000", \emph{etc.} with their positive equivalents, even when embedded
#' within larger strings (\emph{e.g.,} "(-0.00, 1.23)" becomes "(0.00, 1.23)").
#' 
#' @param x Character vector of formatted numbers.
#' @param marks Optional list with \code{big.mark} and \code{decimal.mark} as
#'   returned by \code{\link{resolve_number_marks}}. When \code{NULL}, uses
#'   the default US period decimal.
#' @return Character vector with negative zeros corrected.
#' @keywords internal
fix_negative_zero <- function(x, marks = NULL) {
    if (!is.null(marks)) {
        return(apply_decimal_mark(x, marks))
    }
    ## Fallback: US locale (period decimal)
    gsub("(?<![0-9])-0(\\.0+)(?![0-9])", "0\\1", x, perl = TRUE)
}


#' Format \emph{p}-values for display
#' 
#' Converts numeric \emph{p}-values to formatted character strings using vectorized
#' operations. Values below the threshold (determined by digits parameter) 
#' display as "< 0.001" (for digits=3), "< 0.0001" (for digits=4), \emph{etc.}
#' NA values display as "-".
#' 
#' @param p Numeric vector of \emph{p}-values.
#' @param digits Integer number of decimal places. Also determines the threshold
#'   for "less than" display: threshold = 10^(-digits). Default is 3.
#' @param marks Optional list with \code{big.mark} and \code{decimal.mark} as
#'   returned by \code{\link{resolve_number_marks}}.
#' @return Character vector of formatted \emph{p}-values.
#' @keywords internal
format_pvalues_fit <- function(p, digits = 3, marks = NULL) {
    ## Calculate threshold based on digits
    threshold <- 10^(-digits)
    
    if (!is.null(marks)) {
        less_than_string <- paste0("< 0", marks$decimal.mark,
                                   strrep("0", digits - 1), "1")
    } else {
        less_than_string <- paste0("< ", format(threshold, scientific = FALSE))
    }
    
    ## Pre-compute format string
    fmt_str <- paste0("%.", digits, "f")
    
    ## Vectorized formatting
    result <- sprintf(fmt_str, p)
    result <- fix_negative_zero(result, marks)
    result[is.na(p)] <- "-"
    result[!is.na(p) & p < threshold] <- less_than_string
    
    result
}


### * Fit validation

#' Check if outcome is a Surv() expression
#' 
#' Tests whether an outcome specification string represents a survival outcome
#' by checking for the Surv() function pattern. Used to route model fitting
#' to Cox proportional hazards methods.
#' 
#' @param outcome Character string of the outcome specification.
#' @return Logical \code{TRUE} if outcome starts with "Surv(", \code{FALSE} otherwise.
#' @keywords internal
is_surv_outcome <- function(outcome) {
    grepl("^Surv\\s*\\(", outcome)
}

#' Detect outcome type from data
#' 
#' Automatically determines whether an outcome variable is binary, continuous,
#' or count-based by examining the data values. Used for automatic model type
#' selection and validation. Binary outcomes have 2 unique values, continuous
#' have many values or non-integers, counts have integers >= 0.
#' 
#' @param data Data frame or data.table containing the outcome variable.
#' @param outcome Character string naming the outcome variable.
#' @return Character string: "binary", "continuous", "count", or "unknown".
#' @keywords internal
detect_outcome_type <- function(data, outcome) {
    if (!outcome %in% names(data)) return("unknown")
    
    y <- data[[outcome]]
    
    if (is.factor(y) && length(levels(y)) == 2) return("binary")
    
    if (is.numeric(y)) {
        unique_vals <- unique(y[!is.na(y)])
        if (length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))) {
            return("binary")
        }
        if (all(y[!is.na(y)] >= 0) && 
            all(y[!is.na(y)] == floor(y[!is.na(y)])) &&
            max(y, na.rm = TRUE) > 1) {
            return("count")
        }
        return("continuous")
    }
    "unknown"
}

#' Validate model type matches outcome specification
#' 
#' Ensures consistency between the specified model type, outcome variable type,
#' and GLM family (if applicable). Detects common mismatches like using survival
#' outcomes with non-survival models or binary outcomes with linear models.
#' Can auto-correct fixable issues or raise informative errors.
#'
#' @param outcome Character string outcome specification (may include Surv()).
#' @param model_type Character string specified model type.
#' @param family GLM family object, function, or string if applicable.
#' @param data Data frame or data.table for outcome type detection.
#' @param auto_correct Logical whether to auto-correct fixable mismatches.
#' @return List with model_type, family, messages, auto_corrected flag.
#' @keywords internal
validate_model_outcome <- function(outcome, model_type, family = NULL, 
                                   data = NULL, auto_correct = TRUE) {
    
    corrected_type <- model_type
    is_survival <- is_surv_outcome(outcome)
    
    ## coxph and coxme require Surv() syntax
    ## clogit can work with binary outcome + strata, so don't require Surv()
    surv_required_models <- c("coxph", "coxme")
    non_survival_models <- c("glm", "lm", "lmer", "glmer")
    
    ## Surv() outcome with non-survival model
    if (is_survival && model_type %in% non_survival_models) {
        if (auto_correct) {
            corrected_type <- "coxph"
            message(sprintf(
                "Survival outcome detected but model_type='%s' specified. Switching to 'coxph'.",
                model_type))
        } else {
            stop(sprintf(
                "Survival outcome detected but model_type='%s' specified. Use 'coxph', 'clogit', or 'coxme'.",
                model_type), call. = FALSE)
        }
    }
    
    ## Non-Surv outcome with model that requires Surv()
    if (!is_survival && model_type %in% surv_required_models) {
        stop(sprintf(
            "model_type='%s' requires Surv() syntax.\nExample: outcome = \"Surv(time, status)\"\nGot: \"%s\"",
            model_type, outcome), call. = FALSE)
    }
    
    ## GLM family checks - only when family is a character string
    if (!is.null(data) && model_type == "glm" && !is.null(family) && !is_survival) {
        ## Convert family to string name if it's a family object
        family_name <- if (is.character(family)) {
            family
        } else if (is.function(family)) {
            ## family is a function like binomial, gaussian, etc.
            family()$family
        } else if (is.list(family) && "family" %in% names(family)) {
            ## family is already evaluated (e.g., Gamma(link="log"))
            family$family
        } else {
            NULL
        }
        
        if (!is.null(family_name)) {
            outcome_type <- detect_outcome_type(data, outcome)
            
            if (outcome_type == "categorical" && family_name == "binomial") {
                warning(sprintf(
                    "Categorical outcome '%s' has more than 2 levels. Binomial GLM will coerce this to binary (first level vs all others), which is likely not intended. Consider: (1) recoding to a true binary variable, (2) using multinomial regression (nnet::multinom), (3) using ordinal regression (MASS::polr or ordinal::clm) if levels are ordered, or (4) using a different outcome.",
                    outcome), call. = FALSE)
            }
            if (outcome_type == "continuous" && family_name == "binomial") {
                stop(sprintf(
                    "Continuous outcome '%s' with family='binomial'. Use model_type='lm' or family='gaussian'.",
                    outcome), call. = FALSE)
            }
            if (outcome_type == "binary" && family_name == "gaussian") {
                warning(sprintf(
                    "Binary outcome '%s' with family='gaussian'. Consider family='binomial'.",
                    outcome), call. = FALSE)
            }
        }
    }
    
    ## glmer family checks
    if (!is.null(data) && model_type == "glmer" && !is.null(family) && !is_survival) {
        family_name <- if (is.character(family)) {
            family
        } else if (is.function(family)) {
            family()$family
        } else if (is.list(family) && "family" %in% names(family)) {
            family$family
        } else {
            NULL
        }
        
        if (!is.null(family_name) && family_name == "binomial") {
            outcome_type <- detect_outcome_type(data, outcome)
            if (outcome_type == "categorical") {
                warning(sprintf(
                    "Categorical outcome '%s' has more than 2 levels. Binomial GLMER will coerce this to binary (first level vs all others), which is likely not intended. Consider recoding to a true binary variable or using a different outcome.",
                    outcome), call. = FALSE)
            }
        }
    }
    
    ## lm with binary outcome
    if (!is.null(data) && model_type == "lm" && !is_survival) {
        if (detect_outcome_type(data, outcome) == "binary") {
            warning(sprintf(
                "Binary outcome '%s' with model_type='lm'. Consider model_type='glm' with family='binomial'.",
                outcome), call. = FALSE)
        }
    }
    
    list(model_type = corrected_type, 
         family = family,
         auto_corrected = corrected_type != model_type)
}

#' Write variable names as a model formula reads them
#'
#' The variables of \code{fit()}, \code{uniscreen()}, \code{multifit()}, and
#' \code{fullfit()} are given as column names, from which the model formula
#' is written as text. A column name that is not a syntactic name, such as
#' \code{"tumor stage"}, is read in a formula as code unless it is quoted
#' with backticks. A term that is a column name, each component of an
#' interaction term (\code{"a:b"}) whose components are all column names,
#' and each argument of a \code{Surv()} outcome that is a column name are
#' quoted where the name requires it. A column name may also be given with
#' backticks. Other terms, such as function calls and random effects, are
#' kept as given, as is a term in which no name changes.
#'
#' @param terms Character vector of terms, or \code{NULL}.
#' @param var_names Character vector of column names.
#' @param quote Logical. If \code{TRUE} (default), names are written as a
#'   formula reads them. If \code{FALSE}, the backticks of column names are
#'   removed instead, so that the terms match the columns of the data;
#'   \code{Surv()} outcomes are then kept as given.
#' @return Character vector of terms, or \code{NULL} when \code{terms} is
#'   \code{NULL}.
#' @keywords internal
formula_names <- function(terms, var_names, quote = TRUE) {
    if (is.null(terms)) {
        return(NULL)
    }

    ## The column a name refers to, given with or without backticks, or NA
    column_of <- function(x) {
        x <- trimws(x)
        if (nchar(x) > 2L && startsWith(x, "`") && endsWith(x, "`")) {
            x <- substr(x, 2L, nchar(x) - 1L)
        }
        if (x %in% var_names) x else NA_character_
    }
    write_name <- function(name) {
        if (quote && make.names(name) != name) paste0("`", name, "`") else name
    }
    ## The parts with their column names written, or NULL when none changes
    rewrite <- function(parts) {
        trimmed <- trimws(parts)
        written <- trimmed
        for (k in seq_along(parts)) {
            column <- column_of(parts[k])
            if (!is.na(column)) {
                written[k] <- write_name(column)
            }
        }
        if (identical(written, trimmed)) NULL else written
    }

    vapply(terms, function(term) {
        column <- column_of(term)
        if (!is.na(column)) {
            return(write_name(column))
        }

        ## Arguments of a Surv() outcome
        surv <- regmatches(term, regexec("^\\s*Surv\\s*\\((.*)\\)\\s*$", term))[[1L]]
        if (length(surv) == 2L) {
            written <- if (quote) rewrite(strsplit(surv[2L], ",", fixed = TRUE)[[1L]])
            if (is.null(written)) {
                return(term)
            }
            return(paste0("Surv(", paste(written, collapse = ", "), ")"))
        }

        ## Components of an interaction term, all of which must be columns
        parts <- strsplit(term, ":", fixed = TRUE)[[1L]]
        if (length(parts) > 1L &&
            !anyNA(vapply(parts, column_of, character(1), USE.NAMES = FALSE))) {
            written <- rewrite(parts)
            if (!is.null(written)) {
                return(paste(written, collapse = ":"))
            }
        }
        term
    }, character(1), USE.NAMES = FALSE)
}


#' Validate outcome exists in data
#'
#' Checks that the specified outcome variable (or survival variables within
#' Surv() expression) exists in the dataset. Raises informative error if
#' variables are missing. Handles both simple outcomes and Surv() expressions.
#' 
#' @param data Data frame or data.table to check.
#' @param outcome Character string outcome specification (may include Surv()).
#' @return Invisible \code{TRUE} if validation passes, otherwise stops with error.
#' @keywords internal
validate_outcome_exists <- function(data, outcome) {
    if (is_surv_outcome(outcome)) {
        surv_content <- gsub("^Surv\\s*\\((.*)\\)$", "\\1", outcome)
        surv_vars <- trimws(unlist(strsplit(surv_content, ",")))
        ## A non-syntactic name may be quoted with backticks
        surv_vars <- gsub("^`(.*)`$", "\\1", surv_vars)
        missing <- surv_vars[!surv_vars %in% names(data)]
        if (length(missing) > 0) {
            stop("Survival variable(s) not found: ", paste(missing, collapse = ", "),
                 call. = FALSE)
        }
    } else if (!outcome %in% names(data)) {
        stop("Outcome '", outcome, "' not found in data.", call. = FALSE)
    }
    invisible(TRUE)
}

#' Validate predictors exist in data
#' 
#' Checks that all specified predictor variables exist in the dataset. Handles
#' interaction terms (splits on ":"), mixed-effects random effects (ignores
#' "|" syntax), and raises informative errors for missing variables. A name
#' that contains ":" but is the name of a column is a variable, not an
#' interaction.
#' 
#' @param data Data frame or data.table to check.
#' @param predictors Character vector of predictor variable names.
#' @return Invisible \code{TRUE} if validation passes, otherwise stops with error.
#' @keywords internal
validate_predictors_exist <- function(data, predictors) {
    ## Remove random effects and extract interaction components
    fixed <- predictors[!grepl("\\|", predictors)]
    is_interaction <- grepl(":", fixed, fixed = TRUE) & !fixed %in% names(data)
    if (any(is_interaction)) {
        interaction_vars <- unlist(strsplit(fixed[is_interaction], ":", fixed = TRUE))
        fixed <- c(fixed[!is_interaction], interaction_vars)
    }
    
    missing <- fixed[!fixed %in% names(data)]
    if (length(missing) > 0) {
        stop("Predictor(s) not found: ", paste(missing, collapse = ", "),
             call. = FALSE)
    }
    invisible(TRUE)
}

#' Complete input validation for fit functions
#'
#' Master validation function called by fit(), uniscreen(), fullfit(). Performs
#' comprehensive checks on data structure, variable existence, numeric parameter
#' ranges, and model-outcome consistency. Returns validated parameters with
#' auto-corrections applied when appropriate.
#'
#' @param data Data frame or data.table containing all variables.
#' @param outcome Character string outcome specification (may include Surv()).
#' @param predictors Character vector of predictor variable names.
#' @param model_type Character string model type to validate.
#' @param family GLM family object, function, or string if applicable.
#' @param conf_level Numeric confidence level (must be between 0 and 1).
#' @param digits Integer number of decimal places for effect estimates.
#' @param p_digits Integer number of decimal places for \emph{p}-values.
#' @param p_threshold Numeric \emph{p}-value threshold for significance highlighting.
#' @param auto_correct_model Logical whether to auto-correct model type mismatches.
#' @return List with validated model_type, family, auto_corrected flag.
#' @keywords internal
validate_fit_inputs <- function(data, outcome, predictors, model_type,
                                family = NULL, conf_level = 0.95,
                                digits = 2, p_digits = 3, p_threshold = NULL,
                                auto_correct_model = TRUE) {
    
    ## Basic data check
    if (is.null(data) || !is.data.frame(data) || nrow(data) == 0) {
        stop("'data' must be a non-empty data.frame.", call. = FALSE)
    }
    
    ## Check variables exist
    validate_outcome_exists(data, outcome)
    validate_predictors_exist(data, predictors)
    
    ## Numeric parameter checks
    if (conf_level <= 0 || conf_level >= 1) {
        stop("'conf_level' must be between 0 and 1.", call. = FALSE)
    }
    if (digits < 0 || digits != floor(digits)) {
        stop("'digits' must be a non-negative integer.", call. = FALSE)
    }
    if (p_digits < 0 || p_digits != floor(p_digits)) {
        stop("'p_digits' must be a non-negative integer.", call. = FALSE
)
    }
    if (!is.null(p_threshold) && (p_threshold < 0 || p_threshold > 1)) {
        stop("'p_threshold' must be between 0 and 1.", call. = FALSE)
    }
    
    ## Model-outcome validation
    validation <- validate_model_outcome(
        outcome = outcome,
        model_type = model_type,
        family = family,
        data = data,
        auto_correct = auto_correct_model
    )
    
    validation
}


### * Miscellaneous functions

#' Fit a model with selective warning suppression
#'
#' Wraps model fitting expressions to suppress routine warnings from mixed-effects
#' and GLM fitting (\emph{e.g.,} singular fits, convergence messages, separation
#' warnings) while allowing unexpected warnings through. When \code{verbose = TRUE},
#' all warnings are displayed.
#'
#' @param expr An unevaluated expression (model fitting call) to execute.
#' @param verbose Logical. If \code{TRUE}, all warnings are shown. If \code{FALSE},
#'   routine fitting warnings are suppressed. Default is \code{FALSE}.
#' @return The result of evaluating \code{expr}.
#' @keywords internal
quiet_fit <- function(expr, verbose = FALSE) {
    if (verbose) {
        return(eval.parent(expr))
    }

    ## Diagnostics that are expected during screening and are reported through
    ## the fitted result rather than the console
    fit_noise <- paste0("singular|convergence|step size|maxfun|",
                        "fitted probabilities numerically 0 or 1|",
                        "fitted rates numerically 0|",
                        "algorithm did not converge|",
                        "number of observations.*less than")

    ## Both condition classes are handled: lme4 raises a singular fit as a
    ## message rather than a warning, so a warning handler alone leaves
    ## "boundary (singular) fit" on the console.
    withCallingHandlers(
        eval.parent(expr),
        warning = function(w) {
            if (grepl(fit_noise, conditionMessage(w), ignore.case = TRUE)) {
                invokeRestart("muffleWarning")
            }
        },
        message = function(m) {
            if (grepl(fit_noise, conditionMessage(m), ignore.case = TRUE)) {
                invokeRestart("muffleMessage")
            }
        }
    )
}


#' Drop factor levels without observations before fitting a Cox model
#'
#' \code{survival::coxph()}, \code{survival::clogit()} and
#' \code{coxme::coxme()} keep the levels of a factor that have no
#' observations among the complete cases of the model, where \code{lm()},
#' \code{glm()} and \pkg{lme4} drop them from the model frame. An empty level
#' is then aliased; when it is the first level, the other levels are
#' estimated against whichever of them is aliased in turn (see
#' \code{effective_reference_level()}), and \code{coxme()} fails to fit. The
#' empty levels of the factors that a model uses are therefore dropped before
#' a Cox model is fitted, so that it takes the first level with observations
#' as its reference, as the other models do.
#'
#' A level counts as empty when no complete case of the model variables
#' takes it, as for the model frame of \code{glm()}, so that a level that
#' occurs only in observations excluded for missing values is dropped as
#' well. Those observations, which the model does not use, are then missing
#' for the factor; the attributes of the factor, such as a variable label,
#' are kept.
#'
#' @param data Data table.
#' @param terms Character vector of the terms of the model, such as the
#'   outcome, predictors, random effects, strata and cluster variable, from
#'   which the variables are read.
#' @return \code{data} when none of its factors among the variables has an
#'   empty level, and otherwise a copy in which they have none.
#' @keywords internal
drop_empty_levels <- function(data, terms) {
    vars <- unique(unlist(lapply(terms, function(term) {
        tryCatch(all.vars(str2lang(term)), error = function(e) term)
    }), use.names = FALSE))
    vars <- intersect(vars, names(data))
    if (length(vars) == 0L) {
        return(data)
    }

    complete <- do.call(stats::complete.cases,
                        unname(lapply(vars, function(v) data[[v]])))
    level_counts <- lapply(vars, function(v) {
        x <- data[[v]]
        if (is.factor(x)) tabulate(x[complete], nlevels(x)) else NULL
    })
    has_empty <- vapply(level_counts, function(counts) {
        !is.null(counts) && any(counts == 0L)
    }, logical(1))

    if (!any(has_empty)) {
        return(data)
    }

    data <- data.table::copy(data)
    for (k in which(has_empty)) {
        x <- data[[vars[k]]]
        kept <- levels(x)[level_counts[[k]] > 0L]
        y <- factor(as.character(x), levels = kept, ordered = is.ordered(x))
        for (a in setdiff(names(attributes(x)), c("levels", "class", "names"))) {
            attr(y, a) <- attr(x, a)
        }
        data.table::set(data, j = vars[k], value = y)
    }
    data
}


#' Stop when a model has no coefficient to report
#'
#' A model none of whose coefficients other than the intercept could be
#' estimated, such as a model whose only predictor is constant, leaves
#' \code{m2dt()} no row to report once aliased coefficients are omitted.
#' The table cannot be formatted, so the cause is reported instead.
#'
#' @param raw_data Data table returned by \code{m2dt()} without the
#'   intercept.
#' @return \code{NULL}, invisibly, when \code{raw_data} has rows.
#' @keywords internal
stop_if_no_estimates <- function(raw_data) {
    if (nrow(raw_data) == 0L) {
        stop("None of the predictors of the model could be estimated: each ",
             "is constant or a linear combination of other terms of the ",
             "model.", call. = FALSE)
    }
    invisible(NULL)
}


### * Parallel fitting

#' Resolve the number of worker processes
#'
#' Uses all but one of the detected cores when no number is requested.
#' \code{parallel::detectCores()} returns \code{NA} when the number of cores
#' cannot be determined, which would otherwise propagate into
#' \code{mc.cores}; one core is used then. Under \code{R CMD check}, the limit
#' of two cores signaled by \code{_R_CHECK_LIMIT_CORES_} applies to any
#' number, requested or detected.
#'
#' @param n_cores Requested number of cores, or \code{NULL}.
#' @return Positive integer.
#' @keywords internal
resolve_core_count <- function(n_cores = NULL) {
    if (is.null(n_cores)) {
        detected <- parallel::detectCores()
        n_cores <- if (is.na(detected)) 1L else max(1L, detected - 1L)
    }
    chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
    if (nzchar(chk) && chk == "true") {
        n_cores <- min(n_cores, 2L)
    }
    as.integer(n_cores)
}


#' Return model-fitting errors as values
#'
#' Wraps a function of one argument so that an error is returned as an
#' object of class \code{"summata_fit_error"} holding the condition message
#' instead of being signaled. Sequential and parallel runs then deliver every
#' failure in the same form to \code{stop_on_failed_fits()}.
#'
#' @param fun Function of one argument.
#' @return Function of one argument.
#' @keywords internal
capture_fit_error <- function(fun) {
    force(fun)
    function(x) {
        tryCatch(fun(x), error = function(e) {
            structure(list(message = conditionMessage(e)), class = "summata_fit_error")
        })
    }
}


#' Stop when any model could not be fitted
#'
#' Checks the results of a run over predictors or outcomes and stops with an
#' error that names every element whose fit failed, with its message. An
#' element has failed when its fit raised an error (captured by
#' \code{capture_fit_error()}), when a forked worker failed outside the
#' fitting function (\code{parallel::mclapply()} returns a
#' \code{"try-error"} object), or when a worker returned no result.
#'
#' @param results List of results, one per element of \code{inputs}.
#' @param inputs Character vector naming the elements.
#' @param what Character string naming one element in the message, such as
#'   \code{"predictor"} or \code{"outcome"}.
#' @return \code{results}, invisibly, when every element was fitted.
#' @keywords internal
stop_on_failed_fits <- function(results, inputs, what) {
    failed <- vapply(results, function(r) {
        is.null(r) || inherits(r, c("summata_fit_error", "try-error"))
    }, logical(1))
    if (!any(failed)) {
        return(invisible(results))
    }
    reasons <- vapply(results[failed], function(r) {
        if (is.null(r)) {
            "the worker process returned no result"
        } else if (inherits(r, "try-error")) {
            trimws(conditionMessage(attr(r, "condition")))
        } else {
            r$message
        }
    }, character(1))
    stop(sprintf("Model fitting failed for %d of %d %ss:\n%s",
                 sum(failed), length(results), what,
                 paste0("  ", inputs[failed], ": ", reasons, collapse = "\n")),
         call. = FALSE)
}
