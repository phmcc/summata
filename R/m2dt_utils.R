### * Extraction functions

#' Get factor levels from model (works with S3 and S4)
#' 
#' Extracts factor level information from fitted model objects. Handles both
#' S3 models (glm, lm, coxph) via xlevels slot and S4 models (lme4) via
#' the model frame.
#' 
#' @param model Fitted model object (S3 or S4).
#' @return Named list of factor levels, or \code{NULL} if no factors present.
#' @keywords internal
get_model_xlevels <- function(model) {
    
    ## For regular S3 models (glm, lm, coxph, clogit)
    if (!isS4(model)) {
        ## Special handling for coxme - it doesn't store xlevels
        if (inherits(model, "coxme")) {
            return(NULL)  # Will be handled separately
        }
        return(model$xlevels)
    }
    
    ## For S4 lme4 models (glmer, lmer)
    if (inherits(model, "merMod")) {
        ## Extract from the model frame
        frame_data <- model@frame
        
        xlevels <- list()
        for (var_name in names(frame_data)) {
            if (is.factor(frame_data[[var_name]])) {
                xlevels[[var_name]] <- levels(frame_data[[var_name]])
            }
        }
        
        return(if (length(xlevels) > 0) xlevels else NULL)
    }
    
    ## Default: return NULL
    return(NULL)
}


#' Get data from model object (works with S3 and S4)
#' 
#' Retrieves the original data used to fit a model. Checks multiple locations
#' including model attributes, $data slot, $model slot, and @frame for S4.
#' 
#' @param model Fitted model object (S3 or S4).
#' @return Data frame or data.table used to fit the model, or \code{NULL} if unavailable.
#' @keywords internal
get_model_data <- function(model) {
    
    ## Try attribute first
    data_attr <- attr(model, "data")
    if (!is.null(data_attr)) {
        return(data_attr)
    }
    
    ## For S3 objects, try model$data
    if (!isS4(model)) {
        if (!is.null(model$data)) {
            return(model$data)
        }
        if (!is.null(model$model)) {
            return(model$model)
        }
    }
    
    ## For S4 lme4 objects, use model frame
    if (inherits(model, "merMod")) {
        return(model@frame)
    }
    
    ## Default
    return(NULL)
}


#' Identify the variables required by a fitted model
#'
#' Returns the names of every variable referenced by a model formula, including
#' the response, stratification factors, offsets and random-effect grouping
#' factors. Terms are expanded where possible, so a model specified with a dot
#' on the right-hand side (\emph{e.g.,} \code{y ~ .}) returns the complete set
#' of variables rather than the literal dot.
#'
#' @param model Fitted model object.
#' @param model_class Character string of the model class.
#' @return Character vector of variable names, or \code{character(0)} when the
#'   formula cannot be recovered.
#' @keywords internal
get_model_variables <- function(model, model_class) {

    ## Mixed-effects and coxme formulas carry the grouping factors, which must
    ## also be non-missing for an observation to enter the fit. The formula is
    ## therefore preferred over the terms object for these classes, because the
    ## terms object describes the fixed effects only.
    if (identical(model_class, "coxme") || inherits(model, "merMod")) {
        return(unique(tryCatch(all.vars(stats::formula(model)),
                               error = function(e) character(0))))
    }

    ## The terms object holds the expanded formula, which matters when the
    ## model was specified with a dot on the right-hand side.
    vars <- tryCatch(all.vars(stats::terms(model)),
                     error = function(e) NULL)

    if (is.null(vars)) {
        vars <- tryCatch(all.vars(stats::formula(model)),
                         error = function(e) character(0))
    }

    return(unique(vars))
}


#' Number of observations used to fit a model
#'
#' Provides a single accessor for the fitted sample size across the model
#' classes supported by \pkg{summata}. Survival models record this quantity in
#' class-specific components rather than through \code{stats::nobs()}.
#'
#' @param model Fitted model object.
#' @param model_class Character string of the model class.
#' @return Numeric scalar giving the number of observations used in fitting, or
#'   \code{NA_real_} when it cannot be determined.
#' @keywords internal
get_model_nobs <- function(model, model_class) {

    n_used <- NA_real_

    if (identical(model_class, "coxme")) {
        ## coxme stores c(events, observations)
        if (!is.null(model$n) && length(model$n) >= 2L) {
            n_used <- as.numeric(model$n[2])
        }
    } else if (model_class %in% c("coxph", "clogit")) {
        if (!is.null(model$n)) {
            n_used <- as.numeric(model$n[1])
        }
    }

    if (is.na(n_used)) {
        n_used <- tryCatch(as.numeric(stats::nobs(model)),
                           error = function(e) NA_real_)
    }

    return(n_used)
}


#' Restrict data to the observations used in model fitting
#'
#' Group-level sample sizes must be counted over the observations that entered
#' the model, not over every row supplied by the user. Several fitted objects
#' carry no model frame from which those observations can be read directly:
#' \code{survival::coxph()} does not store one by default, \pkg{coxme} models
#' never store one, and the memory-conserving fits performed internally by
#' \code{uniscreen()} and \code{multifit()} pass \code{model = FALSE}. In these
#' cases the supplied data still contains the observations dropped during
#' fitting, and counting them inflates group sizes whenever the response or a
#' covariate is missing.
#'
#' @param model Fitted model object.
#' @param model_class Character string of the model class.
#' @param data Data frame or data.table supplied to the model call.
#' @return Data.table restricted to the observations used in fitting. The input
#'   is returned unchanged when the analysis rows cannot be determined.
#' @keywords internal
get_analysis_data <- function(model, model_class, data) {

    if (is.null(data) || !inherits(data, "data.frame")) {
        return(data)
    }

    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }

    n_used <- get_model_nobs(model, model_class)

    ## No observations were dropped, so no restriction is required
    if (!is.na(n_used) && nrow(data) == n_used) {
        return(data)
    }

    model_vars <- get_model_variables(model, model_class)
    model_vars <- model_vars[model_vars %in% names(data)]

    if (length(model_vars) == 0L) {
        return(data)
    }

    keep <- stats::complete.cases(data[, model_vars, with = FALSE])

    return(data[keep])
}


#' Detect if model is univariable or multivariable
#' 
#' Determines whether a model contains one predictor (univariable) or multiple
#' predictors (multivariable) by analyzing coefficient names and factor structure.
#' Handles interactions and random effects appropriately.
#' 
#' @param model Fitted model object.
#' @return Character string: "Univariable" or "Multivariable".
#' @keywords internal
detect_model_type <- function(model) {
    ## For glmer/glmerMod/lmer/lmerMod models, handle random effects specially
    if (inherits(model, c("glmerMod", "lmerMod", "merMod"))) {
        ## Get the fixed effects formula (excluding random effects)
        if (!is.null(model@call$formula)) {
            formula_str <- as.character(model@call$formula)[3]  # RHS of formula
            ## Remove random effects terms (anything with |)
            formula_str <- gsub("\\([^|]*\\|[^)]*\\)", "", formula_str)
            ## Clean up extra spaces and plus signs
            formula_str <- gsub("\\s*\\+\\s*\\+", " +", formula_str)
            formula_str <- gsub("^\\s*\\+\\s*", "", formula_str)
            formula_str <- gsub("\\s*\\+\\s*$", "", formula_str)
            
            ## Split by + to get terms
            terms <- trimws(strsplit(formula_str, "\\+")[[1]])
            terms <- terms[terms != ""]
            
            ## Count unique base variables (handle interactions)
            unique_vars <- character()
            for (term in terms) {
                if (grepl(":", term)) {
                    ## Interaction term - split and add components
                    components <- trimws(strsplit(term, ":")[[1]])
                    unique_vars <- c(unique_vars, components)
                } else {
                    unique_vars <- c(unique_vars, term)
                }
            }
            unique_vars <- unique(unique_vars)
            
            return(data.table::fifelse(length(unique_vars) == 1, "Univariable", "Multivariable"))
        }
    }
    
    ## Original code for other model types
    ## Get coefficient names once
    coef_names <- names(stats::coef(model))
    
    ## Remove intercept from count
    if ("(Intercept)" %in% coef_names) {
        term_names <- coef_names[coef_names != "(Intercept)"]
    } else {
        term_names <- coef_names
    }
    
    ## Identify interaction terms (contain ":")
    is_interaction <- grepl(":", term_names, fixed = TRUE)
    interaction_terms <- term_names[is_interaction]
    main_terms <- term_names[!is_interaction]
    
    ## For models without factor variables, count unique base variables from main effects
    xlevels <- get_model_xlevels(model)
    
    ## Special handling for coxme - reconstruct factor info
    if (inherits(model, "coxme") && is.null(xlevels)) {
        data_source <- get_model_data(model)
        if (!is.null(data_source)) {
            xlevels <- list()
            formula_vars <- all.vars(stats::formula(model))
            if (length(formula_vars) >= 2) {
                predictor_vars <- formula_vars[-(1:2)]  # Remove Surv components
                predictor_vars <- predictor_vars[!grepl("\\|", predictor_vars)]  # Remove random effects
                
                for (var in predictor_vars) {
                    if (var %in% names(data_source) && is.factor(data_source[[var]])) {
                        xlevels[[var]] <- levels(data_source[[var]])
                    }
                }
            }
            if (length(xlevels) == 0) xlevels <- NULL
        }
    }
    
    if (is.null(xlevels) || length(xlevels) == 0) {
        ## Count main effect terms
        n_terms <- length(main_terms)
        
        ## Where interactions are present, extract their base variables as well
        if (length(interaction_terms) > 0) {
            ## Extract variables from interactions (split on ":")
            interaction_vars <- unique(unlist(strsplit(interaction_terms, ":", fixed = TRUE)))
            ## Count any interaction variables not already in main effects
            new_vars <- setdiff(interaction_vars, main_terms)
            n_terms <- n_terms + length(new_vars)
        }
        
        return(data.table::fifelse(n_terms == 1, "Univariable", "Multivariable"))
    }
    
    ## For models with factor variables, count unique base variables
    ##
    ## Each coefficient is assigned to the longest factor name it begins with.
    ## Returning the first name in the list would collapse "sexualityHetero"
    ## onto "sex" whenever both variables appear in the model, understating the
    ## variable count and misreporting a multivariable model as univariable.
    ## Names are compared literally, so a variable carrying an inline
    ## transformation is matched on its own characters.
    xlevel_names <- names(xlevels)
    xlevel_names <- xlevel_names[order(nchar(xlevel_names), decreasing = TRUE)]
    
    matched_var <- vapply(main_terms, function(term) {
        hits <- xlevel_names[startsWith(term, xlevel_names)]
        if (length(hits) == 0L) NA_character_ else hits[1]
    }, character(1), USE.NAMES = FALSE)
    
    factor_vars_present <- unique(matched_var[!is.na(matched_var)])
    n_vars <- length(factor_vars_present)
    
    ## Count unique continuous variables IN MAIN EFFECTS
    continuous_terms <- main_terms[is.na(matched_var)]
    n_continuous <- length(unique(continuous_terms))
    
    ## Total unique variables (from main effects only)
    ## Interactions are excluded from variable count
    n_terms <- n_vars + n_continuous
    
    return(data.table::fifelse(n_terms == 1, "Univariable", "Multivariable"))
}


#' Get readable model type name
#' 
#' Converts model class names to human-readable descriptions. For GLMs,
#' uses the family to provide specific names (\emph{e.g.,} "Logistic", "Poisson").
#' 
#' @param model Fitted model object.
#' @return Character string with readable model type name.
#' @keywords internal
get_model_type_name <- function(model) {
    model_class <- class(model)[1]
    
    ## Remove wrapper classes
    if (model_class == "mmodel") {
        model_class <- class(model)[2]
    }
    
    ## For GLM, be more specific based on family
    if (model_class == "glm") {
        family <- model$family$family
        
        ## Use switch instead of multiple ifelse()
        return(switch(family,
                      binomial = "Logistic",
                      poisson = "Poisson",
                      gaussian = "Linear (GLM)",
                      Gamma = "Gamma",
                      quasibinomial = "Quasi-Binomial",
                      quasipoisson = "Quasi-Poisson",
                      paste0(stringr::str_to_title(family), " GLM")
                      ))
    }
    
    ## Map to readable names for other model types
    return(switch(model_class,
                  lm = "Linear",
                  coxph = "Cox PH",
                  clogit = "Conditional Logistic",
                  coxme = "Mixed Effects Cox",
                  negbin = "Negative Binomial",
                  glmer = "glmerMod",  # Keep full name for clarity
                  glmerMod = "glmerMod",
                  lmer = "Linear Mixed",
                  lmerMod = "Linear Mixed",
                  model_class  # default
                  ))
}


#' Parse term into variable and group
#' 
#' Splits coefficient term names into base variable names and factor levels.
#' For example, "sexMale" becomes variable="sex" and group="Male".
#' Handles interaction terms and continuous variables appropriately.
#' 
#' @param terms Character vector of coefficient term names.
#' @param xlevels Named list of factor levels from the model.
#' @param model Optional model object for extracting factor info from coxme models.
#' @return Data.table with 'variable' and 'group' columns.
#' @keywords internal
parse_term <- function(terms, xlevels = NULL, model = NULL) {
    n_terms <- length(terms)
    
    ## Initialize result vectors
    variable <- character(n_terms)
    group <- character(n_terms)
    
    ## Interactions should not be parsed as factor variables
    is_interaction <- grepl(":", terms, fixed = TRUE)

    ## Special handling for coxme models - reconstruct xlevels if needed
    if (!is.null(model) && inherits(model, "coxme") && is.null(xlevels)) {
        data_source <- attr(model, "data")
        if (is.null(data_source)) {
            data_source <- get_model_data(model)
        }
        
        if (!is.null(data_source)) {
            xlevels <- list()
            
            ## Extract factor structure from coefficient names
            coef_names <- names(coxme::fixef(model))
            
            ## Check each column in the data
            for (col_name in names(data_source)) {
                ## Check if any coefficient starts with this column name
                if (any(grepl(paste0("^", col_name), coef_names))) {
                    ## This is a factor variable - get its levels
                    if (is.factor(data_source[[col_name]])) {
                        xlevels[[col_name]] <- levels(data_source[[col_name]])
                    } else if (is.character(data_source[[col_name]])) {
                        xlevels[[col_name]] <- sort(unique(data_source[[col_name]]))
                    }
                }
            }
            
            if (length(xlevels) == 0) xlevels <- NULL
        }
    }
    
    if (!is.null(xlevels) && length(xlevels) > 0) {
        ## Match each coefficient to the factor name it begins with.
        ##
        ## Names are compared literally rather than as regular expressions, so
        ## that variables carrying an inline transformation (\emph{e.g.,}
        ## \code{factor(stage)}) are matched on their own characters rather than
        ## on the grouping and alternation those characters would denote.
        ##
        ## Names are processed from shortest to longest so that the longest
        ## match is the one retained. Processing in formula order would assign
        ## the coefficient "sexualityHetero" to "sex" whenever both variables
        ## appear in the model and "sex" is written second.
        xlevel_names <- names(xlevels)
        xlevel_names <- xlevel_names[order(nchar(xlevel_names))]
        
        ## For each factor variable, find all matching terms at once, skip interaction terms
        for (var in xlevel_names) {
            ## Find which terms start with this variable name
            matches <- startsWith(terms, var) & !is_interaction  # Skip interactions
            
            if (any(matches)) {
                ## Extract levels for all matching terms at once
                variable[matches] <- var
                group[matches] <- substring(terms[matches], nchar(var) + 1L)
            }
        }
        
        ## Any remaining unmatched terms (including interactions) are treated as-is
        unmatched <- variable == ""
        if (any(unmatched)) {
            variable[unmatched] <- terms[unmatched]
            ## group already initialized to "" for these
        }
    } else {
        ## No factor variables - all continuous (including interactions)
        variable <- terms
        ## group already initialized to ""
    }
    
    return(data.table::data.table(variable = variable, group = group))
}


#' Extract event variable from survival model
#' 
#' Parses the Surv() expression in survival model formulas to extract
#' the event/status variable name. Works with coxph, clogit, and coxme models.
#' 
#' @param model Fitted survival model object.
#' @param model_class Character string of the model class.
#' @return Character string naming the event variable, or \code{NULL} if not found.
#' @keywords internal
get_event_variable <- function(model, model_class) {
    event_var <- NULL
    
    if (model_class %in% c("coxph", "clogit", "coxme")) {
        
        ## Get the formula string differently for each model type
        outcome_str <- NULL
        
        if (model_class == "coxme") {
            ## For coxme, use formulaList$fixed
            if (!is.null(model$formulaList$fixed)) {
                outcome_str <- tryCatch({
                    as.character(model$formulaList$fixed)[2]
                }, error = function(e) NULL)
            }
        } else {
            ## For coxph and clogit, use standard formula
            outcome_str <- tryCatch({
                as.character(stats::formula(model))[2]
            }, error = function(e) NULL)
        }
        
        ## Extract event variable from Surv()
        if (!is.null(outcome_str) && !is.na(outcome_str)) {
            if (grepl("Surv\\(", outcome_str)) {
                ## Remove "Surv(" from beginning and ")" from end
                surv_expr <- gsub("^Surv\\(", "", outcome_str)
                surv_expr <- gsub("\\)$", "", surv_expr)
                
                ## Split by comma to get time and event
                surv_parts <- trimws(strsplit(surv_expr, ",")[[1]])
                if (length(surv_parts) >= 2) {
                    event_var <- surv_parts[2]
                }
            }
        }
    }
    
    return(event_var)
}


#' Identify the response variables of a fitted model
#'
#' Returns the variables appearing on the left-hand side of a model formula. A
#' survival response contributes both the time and status variables.
#'
#' @param model Fitted model object.
#' @return Character vector of response variable names, or \code{character(0)}
#'   when the formula cannot be recovered.
#' @keywords internal
get_outcome_variables <- function(model) {

    fml <- tryCatch(stats::formula(model), error = function(e) NULL)

    if (is.null(fml) || length(fml) < 3L) {
        return(character(0))
    }

    return(tryCatch(all.vars(fml[[2]]), error = function(e) character(0)))
}


#' Count the observations available to an analysis
#'
#' Summarizes how many of the supplied observations could enter a model, and
#' separates those excluded for a missing response from those excluded for a
#' missing covariate. The distinction matters: a missing response and a missing
#' covariate raise different questions about whether a complete-case analysis is
#' appropriate.
#'
#' Several models may be described at once by passing a list of predictor sets,
#' as when a univariable screen fits one model per predictor. The analyzed count
#' is then returned as one element per set.
#'
#' @param data Data frame or data.table supplied to the analysis.
#' @param outcome_vars Character vector of response variable names.
#' @param predictor_vars Character vector of predictor variable names, or a list
#'   of such vectors when several models are described.
#' @param event_var Character name of the event indicator, or \code{NULL} for
#'   models that carry no event count. Factors are coded as the modeling
#'   functions code them, with any level beyond the first counting as an event.
#' @return Named list with \code{n_supplied}, \code{n_analyzed},
#'   \code{n_missing_outcome} and \code{n_missing_predictor}, and, where an
#'   event indicator was supplied, \code{events_supplied} and
#'   \code{events_analyzed}. \code{NULL} when the counts cannot be determined.
#' @keywords internal
get_analysis_counts <- function(data, outcome_vars, predictor_vars = NULL,
                                event_var = NULL) {

    if (is.null(data) || !inherits(data, "data.frame") || nrow(data) == 0) {
        return(NULL)
    }

    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }

    outcome_vars <- unique(outcome_vars[outcome_vars %in% names(data)])

    if (length(outcome_vars) == 0) {
        return(NULL)
    }

    if (!is.list(predictor_vars)) {
        predictor_vars <- list(predictor_vars)
    }

    n_supplied <- nrow(data)
    outcome_ok <- stats::complete.cases(data[, outcome_vars, with = FALSE])

    ## Event indicator, coded as the modeling functions code it. A count
    ## response contributes its counts rather than a binary flag, which
    ## mirrors the total that m2dt() reports for Poisson and negative
    ## binomial models.
    event_indicator <- NULL
    if (!is.null(event_var) && length(event_var) == 1 &&
        event_var %in% names(data)) {
        y <- data[[event_var]]
        event_indicator <- if (is.factor(y)) {
                               as.integer(y) > 1L
                           } else if (is.numeric(y) || is.logical(y)) {
                               y
                           } else {
                               NULL
                           }
    }

    counts <- lapply(predictor_vars, function(vars) {
        vars <- unique(vars[vars %in% names(data) & !vars %in% outcome_vars])

        predictor_ok <- if (length(vars) > 0) {
                            stats::complete.cases(data[, vars, with = FALSE])
                        } else {
                            rep(TRUE, n_supplied)
                        }

        analyzed <- outcome_ok & predictor_ok

        c(n_analyzed = sum(analyzed),
          n_missing_predictor = sum(outcome_ok & !predictor_ok),
          events_analyzed = if (is.null(event_indicator)) {
                                NA_real_
                            } else {
                                sum(event_indicator[analyzed], na.rm = TRUE)
                            })
    })

    result <- list(
        n_supplied = n_supplied,
        n_analyzed = vapply(counts, function(x) x[["n_analyzed"]], numeric(1)),
        n_missing_outcome = sum(!outcome_ok),
        n_missing_predictor = vapply(counts, function(x) x[["n_missing_predictor"]], numeric(1))
    )

    if (!is.null(event_indicator)) {
        result$events_supplied <- sum(event_indicator, na.rm = TRUE)
        result$events_analyzed <- vapply(counts, function(x) x[["events_analyzed"]],
                                         numeric(1))
    }

    return(result)
}


#' Count the observations available to a fitted model
#'
#' Convenience wrapper around \code{get_analysis_counts()} that derives the
#' response and predictor variables from the model itself.
#'
#' @param model Fitted model object.
#' @param model_class Character string of the model class.
#' @param data Data frame or data.table supplied to the model call.
#' @return See \code{get_analysis_counts()}.
#' @keywords internal
get_model_analysis_counts <- function(model, model_class, data) {

    outcome_vars <- get_outcome_variables(model)

    if (length(outcome_vars) == 0) {
        return(NULL)
    }

    predictor_vars <- setdiff(get_model_variables(model, model_class), outcome_vars)

    ## stats::family() resolves for glm, lm and merMod objects and errors for
    ## survival fits, which are identified by their class instead
    family_name <- tryCatch(stats::family(model)$family, error = function(e) NULL)
    event_var <- get_event_variable_for_counts(outcome_vars, model_class, family_name)

    return(get_analysis_counts(data, outcome_vars, predictor_vars, event_var))
}


#' Identify the event indicator among a model's response variables
#'
#' Determines which response variable carries the event count, and whether the
#' model has one at all. Linear and Gaussian models do not; the families that
#' count events are those for which \code{m2dt()} reports an events figure, so
#' that a denominator is only offered where a numerator exists.
#'
#' @param outcome_vars Character vector of response variable names.
#' @param model_type Character string of the model type or class.
#' @param family Model family, where applicable. Accepted as a name, a family
#'   object, or a generator function, since callers resolve the family at
#'   different points: \code{fit()} passes the name it was given, while
#'   \code{uniscreen()} may already have resolved \code{"Gamma"} to
#'   \code{Gamma(link = "log")}.
#' @return Character name of the event indicator, or \code{NULL} when the model
#'   carries no event count.
#' @keywords internal
get_event_variable_for_counts <- function(outcome_vars, model_type = NULL,
                                          family = NULL) {

    if (length(outcome_vars) == 0) {
        return(NULL)
    }

    ## The family may arrive as a name, a family object, or the generator
    ## function itself. Comparing an unresolved family object against a
    ## character vector yields one value per list element rather than one
    ## logical, so it is reduced to its name first.
    family_name <- if (inherits(family, "family")) {
                       family$family
                   } else if (is.function(family)) {
                       tryCatch(family()$family, error = function(e) NULL)
                   } else if (is.character(family)) {
                       family[1]
                   } else {
                       NULL
                   }

    survival_types <- c("coxph", "clogit", "coxme")
    count_types <- c("negbin", "glm.nb")
    glm_types <- c("glm", "glmer", "glmerMod")
    event_families <- c("binomial", "quasibinomial", "poisson", "quasipoisson")

    survival_response <- if (length(outcome_vars) >= 2) {
                             outcome_vars[2]
                         } else {
                             outcome_vars[1]
                         }

    ## Without a usable model type, a two-variable response is survival
    if (length(model_type) != 1 || is.na(model_type)) {
        return(survival_response)
    }

    ## A survival response carries the time variable first and the status
    ## variable second; the event indicator is the status variable
    if (model_type %in% survival_types) {
        return(survival_response)
    }

    if (model_type %in% count_types) {
        return(outcome_vars[1])
    }

    if (model_type %in% glm_types &&
        length(family_name) == 1 && family_name %in% event_families) {
        return(outcome_vars[1])
    }

    return(NULL)
}


### * Formatting functions

#' Describe the analyzed sample for console output
#'
#' Renders the counts produced by \code{get_analysis_counts()} as a single line
#' for the \code{print()} methods. The line is produced whenever the counts are
#' available, including when every supplied observation entered the analysis:
#' a complete sample is itself worth stating, and reporting it only on loss
#' would leave the reader to infer the difference between no exclusions and no
#' disclosure.
#'
#' @param counts List as returned by \code{get_analysis_counts()}.
#' @param label Character label preceding the counts.
#' @param marks List of number marks as returned by
#'   \code{resolve_number_marks()}. Counts and percentages are separated as
#'   the accompanying table separates them, so that a locale setting applies
#'   to the whole of the output rather than to the table alone. Resolved from
#'   the global option when not supplied.
#' @return Character string, or \code{NULL} when the counts are unavailable.
#' @keywords internal
format_analysis_counts <- function(counts, label = "Observations analyzed",
                                   marks = NULL) {

    if (is.null(counts) || length(counts$n_analyzed) == 0) {
        return(NULL)
    }

    n_supplied <- counts$n_supplied
    lo <- min(counts$n_analyzed)
    hi <- max(counts$n_analyzed)

    ## The denominator may be absent, as when a total is unavailable for the
    ## quantity being described. Nothing is reported rather than a bare
    ## numerator dressed as a proportion.
    if (length(n_supplied) == 0 || !is.finite(n_supplied) || n_supplied == 0) {
        return(NULL)
    }

    if (!is.finite(lo) || !is.finite(hi)) {
        return(NULL)
    }

    if (is.null(marks)) {
        marks <- tryCatch(resolve_number_marks(NULL),
                          error = function(e) list(big.mark = "",
                                                   decimal.mark = "."))
    }

    ## Counts and percentages are rendered by the same helpers the tables use,
    ## so that a figure in this line and the same figure in the n column below
    ## cannot be formatted differently
    count_str <- function(x) format_count(x, marks)

    pct_str <- function(x) apply_decimal_mark(sprintf("%.1f", x), marks)

    if (lo == hi) {
        return(sprintf("%s: %s of %s (%s%%)",
                       label, count_str(lo), count_str(n_supplied),
                       pct_str(100 * lo / n_supplied)))
    }

    return(sprintf("%s: %s-%s of %s (%s-%s%%)",
                   label, count_str(lo), count_str(hi), count_str(n_supplied),
                   pct_str(100 * lo / n_supplied),
                   pct_str(100 * hi / n_supplied)))
}


#' Describe the analyzed events for console output
#'
#' Renders the event counts held in \code{get_analysis_counts()} using the same
#' formatter as the observation counts, so that the two lines cannot drift apart
#' in rounding, separators or range handling.
#'
#' @param counts List as returned by \code{get_analysis_counts()}.
#' @param label Character label preceding the counts.
#' @param marks List of number marks as returned by
#'   \code{resolve_number_marks()}.
#' @return Character string, or \code{NULL} when the model carries no events.
#' @keywords internal
format_event_counts <- function(counts, label = "Events analyzed",
                                marks = NULL) {

    if (is.null(counts) || is.null(counts$events_analyzed)) {
        return(NULL)
    }

    return(format_analysis_counts(
        list(n_supplied = counts$events_supplied,
             n_analyzed = counts$events_analyzed),
        label = label,
        marks = marks
    ))
}
