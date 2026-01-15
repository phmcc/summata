#' Fix negative zero in formatted strings
#' 
#' Corrects floating-point rounding artifacts that produce "-0.00" or similar
#' negative zero strings. Works on character vectors, replacing patterns like
#' "-0.00", "-0.000", etc. with their positive equivalents, even when embedded
#' within larger strings (e.g., "(-0.00, 1.23)" becomes "(0.00, 1.23)").
#' 
#' @param x Character vector of formatted numbers.
#' @return Character vector with negative zeros corrected.
#' @keywords internal
fix_negative_zero <- function(x) {
    ## Match -0 followed by decimal point and zeros, with word boundaries
    ## Handles: "-0.0", "-0.00", "-0.000", etc. anywhere in string
    gsub("(?<![0-9])-0(\\.0+)(?![0-9])", "0\\1", x, perl = TRUE)
}

#' Format model results for publication-ready display
#' 
#' Transforms raw model coefficient data into a formatted table suitable for
#' publication. Handles effect measure formatting (OR, HR, RR, Estimate),
#' confidence intervals, p-values, sample sizes, and variable labels.
#' Supports interaction terms and mixed-effects models.
#' 
#' @param data Data.table containing raw model results with coefficient columns.
#' @param effect_col Optional character string specifying the effect column name.
#'   If NULL, auto-detects from OR, HR, RR, or Estimate columns.
#' @param digits Integer number of decimal places for effect estimates.
#' @param p_digits Integer number of decimal places for p-values.
#' @param labels Optional named character vector mapping variable names to display labels.
#'   Supports automatic labeling of interaction terms.
#' @param show_n Logical whether to include sample size column.
#' @param show_events Logical whether to include events column (ignored for linear models).
#' @param reference_label Character string to display for reference categories.
#' @param exponentiate Optional logical to force exponentiated (TRUE) or raw (FALSE)
#'   coefficient display. If NULL, uses existing columns.
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
                               exponentiate = NULL) {
    
    ## Determine which columns we actually need to avoid copying everything
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

    ## Disallow "Events" column if linear model 
    if ("model_type" %in% names(result)) {
        model_type <- result$model_type[1]
        if (grepl("Linear", model_type, ignore.case = TRUE)) {
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
        
        ## Create lookup table for main effects
        label_dt <- data.table::data.table(
                                    var_orig = names(labels),
                                    var_new = unname(unlist(labels))
                                )
        
        ## Update main effect variable names using merge
        result[label_dt, Variable := i.var_new, on = .(Variable = var_orig)]
        
        ## Handle interaction terms (contain ":")
        interaction_mask <- grepl(":", result$Variable, fixed = TRUE)
        
        if (any(interaction_mask)) {
            ## Process all interaction terms at once
            interaction_vars <- result$Variable[interaction_mask]
            
            ## Pre-sort labels by length (longest first) once
            sorted_label_names <- names(labels)[order(-nchar(names(labels)))]
            
            ## Vectorized interaction labeling
            labeled_interactions <- vapply(interaction_vars, function(original_var) {
                ## Check for direct custom label first
                if (original_var %in% names(labels)) {
                    return(labels[[original_var]])
                }
                
                ## Split and process components
                components <- strsplit(original_var, ":", fixed = TRUE)[[1]]
                
                labeled_parts <- vapply(components, function(comp) {
                    ## Try to match against base variable names
                    for (base_var in sorted_label_names) {
                        if (startsWith(comp, base_var)) {
                            suffix <- substring(comp, nchar(base_var) + 1)
                            if (nchar(suffix) == 0) {
                                return(labels[[base_var]])
                            } else {
                                return(paste0(labels[[base_var]], " (", suffix, ")"))
                            }
                        }
                    }
                    return(comp)  # No match, keep original
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
        result[n_vals >= 1000, n := format(n_vals[n_vals >= 1000], big.mark = ",")]
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
        
        ## Vectorized formatting using fifelse to avoid subsetting issues
        result[, events := data.table::fcase(
                                           is.na(event_vals), NA_character_,
                                           event_vals >= 1000, format(round(event_vals), big.mark = ",", scientific = FALSE),
                                           default = as.character(round(event_vals))
                                       )]
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
    
    ## Create appropriate label
    if (model_scope == "Univariable") {
        ## Univariable: OR, HR, RR, Coefficient
        effect_label <- paste0(effect_col, " (95% CI)")
    } else if (model_scope == "Multivariable") {
        ## Multivariable: aOR, aHR, aRR, Adj. Coefficient
        adjusted_col <- switch(effect_col,
                               "OR" = "aOR",
                               "HR" = "aHR", 
                               "RR" = "aRR",
                               "Coefficient" = "Adj. Coefficient",
                               paste0("Adj. ", effect_col)
                               )
        effect_label <- paste0(adjusted_col, " (95% CI)")
    } else {
        effect_label <- paste0(effect_col, " (95% CI)")
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
        
        ## Pre-compute format string once
        if (effect_col %in% c("OR", "HR", "RR")) {
            fmt_str <- paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)")
        } else {
            fmt_str <- paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)")
        }
        
        ## Vectorized sprintf is faster than per-row fcase with sprintf
        formatted_effects <- sprintf(fmt_str, effect_vals, ci_lower_vals, ci_upper_vals)
        formatted_effects <- fix_negative_zero(formatted_effects)
        formatted_effects[is.na(effect_vals)] <- ""
        formatted_effects[is_reference] <- reference_label
        
        result[, (effect_label) := formatted_effects]
    }
    
    ## Format p-values
    if ("p_value" %in% names(result)) {
        result[, `p-value` := format_pvalues_fit(p_value, p_digits)]
        
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

#' Format p-values for display
#' 
#' Converts numeric p-values to formatted character strings using vectorized
#' operations. Values below the threshold (determined by digits parameter) 
#' display as "< 0.001" (for digits=3), "< 0.0001" (for digits=4), etc.
#' NA values display as "-".
#' 
#' @param p Numeric vector of p-values.
#' @param digits Integer number of decimal places. Also determines the threshold
#'   for "less than" display: threshold = 10^(-digits). Default is 3.
#' @return Character vector of formatted p-values.
#' @keywords internal
format_pvalues_fit <- function(p, digits = 3) {
    ## Calculate threshold based on digits
    threshold <- 10^(-digits)
    less_than_string <- paste0("< ", format(threshold, scientific = FALSE))
    
    ## Pre-compute format string
    fmt_str <- paste0("%.", digits, "f")
    
    ## Vectorized formatting (faster than fcase for simple conditions)
    result <- sprintf(fmt_str, p)
    result <- fix_negative_zero(result)
    result[is.na(p)] <- "-"
    result[!is.na(p) & p < threshold] <- less_than_string
    
    result
}

#' Check if outcome is a Surv() expression
#' 
#' Tests whether an outcome specification string represents a survival outcome
#' by checking for the Surv() function pattern. Used to route model fitting
#' to Cox proportional hazards methods.
#' 
#' @param outcome Character string of the outcome specification.
#' @return Logical TRUE if outcome starts with "Surv(", FALSE otherwise.
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
#' Checks for mismatches and auto-corrects or errors as appropriate.
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

#' Validate outcome exists in data
#' 
#' Checks that the specified outcome variable (or survival variables within
#' Surv() expression) exists in the dataset. Raises informative error if
#' variables are missing. Handles both simple outcomes and Surv() expressions.
#' 
#' @param data Data frame or data.table to check.
#' @param outcome Character string outcome specification (may include Surv()).
#' @return Invisible TRUE if validation passes, otherwise stops with error.
#' @keywords internal
validate_outcome_exists <- function(data, outcome) {
    if (is_surv_outcome(outcome)) {
        surv_content <- gsub("^Surv\\s*\\((.*)\\)$", "\\1", outcome)
        surv_vars <- trimws(unlist(strsplit(surv_content, ",")))
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
#' "|" syntax), and raises informative errors for missing variables.
#' 
#' @param data Data frame or data.table to check.
#' @param predictors Character vector of predictor variable names.
#' @return Invisible TRUE if validation passes, otherwise stops with error.
#' @keywords internal
validate_predictors_exist <- function(data, predictors) {
    ## Remove random effects and extract interaction components
    fixed <- predictors[!grepl("\\|", predictors)]
    if (any(grepl(":", fixed))) {
        interaction_vars <- unlist(strsplit(fixed[grepl(":", fixed)], ":"))
        fixed <- c(fixed[!grepl(":", fixed)], interaction_vars)
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
#' @param p_digits Integer number of decimal places for p-values.
#' @param p_threshold Numeric p-value threshold for significance highlighting.
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
