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
    
    result <- data.table::copy(data)

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
                
                paste(labeled_parts, collapse = " * ")
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
            result[, n := data.table::fcase(
                                          !is.na(n_group) & n_group >= 1000, format(n_group, big.mark = ","),
                                          !is.na(n_group), as.character(n_group),
                                          !is.na(n) & n >= 1000, format(n, big.mark = ","),
                                          !is.na(n), as.character(n),
                                          default = NA_character_
                                      )]
        } else {
            result[, n := data.table::fcase(
                                          is.na(n), NA_character_,
                                          n >= 1000, format(n, big.mark = ","),
                                          default = as.character(n)
                                      )]
        }
    }

    ## Handle events column
    if ("events" %in% names(result)) {
        has_events_group <- "events_group" %in% names(result)
        
        if (has_events_group) {
            result[, events := data.table::fcase(
                                               !is.na(events_group) & events_group >= 1000, format(events_group, big.mark = ","),
                                               !is.na(events_group), as.character(events_group),
                                               !is.na(events) & events >= 1000, format(events, big.mark = ","),
                                               !is.na(events), as.character(events),
                                               default = NA_character_
                                           )]
        } else {
            result[, events := data.table::fcase(
                                               is.na(events), NA_character_,
                                               events >= 1000, format(events, big.mark = ","),
                                               default = as.character(events)
                                           )]
        }
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
        ci_lower <- result$ci_lower
        ci_upper <- result$ci_upper
        
        ## Use different CI format based on effect type
        if (effect_col %in% c("OR", "HR", "RR")) {
            result[, (effect_label) := data.table::fcase(
                                                       is_reference, reference_label,
                                                       !is.na(effect_vals), sprintf("%.*f (%.*f-%.*f)", 
                                                                                    digits, effect_vals, digits, ci_lower, digits, ci_upper),
                                                       default = ""
                                                   )]
        } else {
            result[, (effect_label) := data.table::fcase(
                                                       is_reference, reference_label,
                                                       !is.na(effect_vals), sprintf("%.*f (%.*f, %.*f)", 
                                                                                    digits, effect_vals, digits, ci_lower, digits, ci_upper),
                                                       default = ""
                                                   )]
        }
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
    ## Calculate threshold based on digits (e.g., digits=3 -> 0.001, digits=4 -> 0.0001)
    threshold <- 10^(-digits)
    less_than_string <- paste0("< ", format(threshold, scientific = FALSE))
    
    data.table::fcase(
                    is.na(p), "-",
                    p < threshold, less_than_string,
                    default = sprintf(paste0("%.", digits, "f"), p)
                )
}
