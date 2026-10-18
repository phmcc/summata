### * Extraction functions

#' Get factor levels from model (works with S3 and S4)
#' 
#' Extracts factor level information from fitted model objects. Handles both
#' S3 models (glm, lm, coxph) via xlevels slot and S4 models (lme4) via
#' the model frame. Logical predictors are included with the levels
#' \code{"FALSE"} and \code{"TRUE"} (see \code{add_logical_levels()}).
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
        return(add_logical_levels(model$xlevels, predictor_classes(model$terms)))
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
        
        ## The first column of the frame is the response; expressions such
        ## as I(age > 65) are not variables of the data
        frame_classes <- vapply(frame_data[-1L], function(x) class(x)[1L],
                                character(1))
        frame_vars <- tryCatch(all.vars(stats::formula(model)),
                               error = function(e) character(0))
        xlevels <- add_logical_levels(
            xlevels, frame_classes[names(frame_classes) %in% frame_vars]
        )

        return(if (length(xlevels) > 0) xlevels else NULL)
    }
    
    ## Default: return NULL
    return(NULL)
}


#' Add the levels of logical predictors to the factor levels of a model
#'
#' R codes a logical predictor as a factor with the levels \code{FALSE} and
#' \code{TRUE}, naming its coefficient by the second (\code{olderTRUE}) with
#' \code{FALSE} as the reference, but records levels only for factors and
#' character variables (\code{xlevels}). With the two levels added, tables
#' and forest plots show a logical predictor by level, with a reference row
#' for \code{FALSE}, as they show a factor, rather than as a single term named
#' by its coefficient.
#'
#' @param xlevels Named list of factor levels, or \code{NULL}.
#' @param classes Named character vector of the classes of the predictors,
#'   as in the \code{"dataClasses"} attribute of a terms object.
#' @return \code{xlevels} with the levels \code{"FALSE"} and \code{"TRUE"}
#'   for each logical predictor it does not already contain.
#' @keywords internal
add_logical_levels <- function(xlevels, classes) {
    logical_vars <- setdiff(names(classes)[classes %in% "logical"], names(xlevels))
    if (length(logical_vars) == 0L) {
        return(xlevels)
    }
    if (is.null(xlevels)) {
        xlevels <- list()
    }
    for (var in logical_vars) {
        xlevels[[var]] <- c("FALSE", "TRUE")
    }
    xlevels
}


#' Read the classes of the predictors of a terms object
#'
#' Returns the classes of the variables of the model other than the
#' response, leaving out expressions such as \code{I(age > 65)}, which are
#' not variables of the data and whose levels cannot be counted there.
#'
#' @param model_terms Terms object of a fitted model, or \code{NULL}.
#' @return Named character vector of the \code{"dataClasses"} of those
#'   variables, or \code{NULL}.
#' @keywords internal
predictor_classes <- function(model_terms) {
    classes <- attr(model_terms, "dataClasses")
    if (is.null(classes)) {
        return(NULL)
    }
    response <- attr(model_terms, "response")
    if (isTRUE(response > 0L)) {
        classes <- classes[-response]
    }
    classes[names(classes) %in% all.vars(model_terms)]
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


#' Term labels of the fixed effects of a model
#'
#' Returns the labels of the fixed-effect terms of a model as R writes them
#' in coefficient names, with non-syntactic variable names enclosed in
#' backticks (\emph{e.g.,} "`high bp`") and calls written out (\emph{e.g.,}
#' "factor(stage)").
#'
#' @param model Fitted model object.
#' @return Character vector of term labels, or \code{character(0)} when the
#'   terms cannot be recovered.
#' @keywords internal
model_term_labels <- function(model) {
    model_terms <- if (inherits(model, "coxme")) {
                       tryCatch(stats::terms(model$formulaList$fixed),
                                error = function(e) NULL)
                   } else {
                       tryCatch(stats::terms(model), error = function(e) NULL)
                   }
    as.character(attr(model_terms, "term.labels"))
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
    ##
    ## A coxme model records its call, whose formula is the name of a variable
    ## when the model is fitted inside a function (fit() fits it from a
    ## formula named "formula"). Its fixed and random parts are read from the
    ## formula list that the model keeps instead.
    if (identical(model_class, "coxme") && !is.null(model$formulaList)) {
        vars <- tryCatch(c(all.vars(model$formulaList$fixed),
                           unlist(lapply(model$formulaList$random, all.vars))),
                         error = function(e) character(0))
        return(unique(vars))
    }
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


#' C-statistic of predictions for a binary outcome
#'
#' Computes the area under the empirical ROC curve as the Mann-Whitney
#' statistic, with tied predictions counted as one half, from a single sort
#' (see \code{average_ranks()}). The orientation matches that of
#' \code{pROC::roc()} with \code{direction = "auto"}, which this function
#' replaces: higher predictions indicate events unless the median prediction
#' of the events is below that of the non-events. The result therefore no
#' longer depends on whether the suggested package is installed.
#'
#' @param y Binary outcome coded 0 and 1, or logical.
#' @param p Numeric vector of predictions of the same length.
#' @return Numeric scalar; \code{NA} when either outcome is absent.
#' @keywords internal
c_statistic <- function(y, p) {
    ok <- !is.na(y) & !is.na(p)
    event <- as.numeric(y[ok]) == 1
    p <- p[ok]
    n_1 <- sum(event)
    n_0 <- length(event) - n_1
    if (n_1 == 0 || n_0 == 0) {
        return(NA_real_)
    }
    ranks <- average_ranks(p)
    auc <- (sum(ranks[event]) - n_1 * (n_1 + 1) / 2) / (as.double(n_1) * n_0)
    if (stats::median(p[!event]) > stats::median(p[event])) 1 - auc else auc
}


#' Confidence intervals for model coefficients
#'
#' Computes confidence intervals for the named coefficients by the method that
#' \code{conf_method} selects. With \code{"profile"}, \code{stats::confint()}
#' gives profile likelihood intervals for GLM and negative binomial models,
#' exact \emph{t} intervals for linear models, and the model's own intervals
#' for other classes (Wald intervals for Cox models); if it fails, as
#' profiling can under separation, Wald intervals are returned instead. With
#' \code{"wald"}, and for quasi-likelihood families, which have no likelihood
#' to profile, Wald intervals are returned. \code{m2dt()} and the forest
#' plot functions share this function so that a table and a plot of the same
#' model show the same intervals.
#'
#' @param model Fitted model object.
#' @param parm Character vector of coefficient names.
#' @param conf_level Numeric confidence level.
#' @param conf_method Character string, \code{"profile"} or \code{"wald"}.
#' @return A matrix with one row per element of \code{parm} and two columns
#'   (lower and upper bounds).
#' @keywords internal
coef_confint <- function(model, parm, conf_level, conf_method) {
    is_quasi <- inherits(model, "glm") &&
        grepl("^quasi", model$family$family)

    conf_int <- if (identical(conf_method, "wald") || is_quasi) {
        stats::confint.default(model, parm = parm, level = conf_level)
    } else {
        tryCatch(
            suppressMessages(suppressWarnings(
                stats::confint(model, parm = parm, level = conf_level)
            )),
            error = function(e) {
                stats::confint.default(model, parm = parm, level = conf_level)
            }
        )
    }

    ## confint() returns a named vector when parm names a single coefficient
    if (is.null(dim(conf_int))) {
        conf_int <- matrix(conf_int, nrow = 1,
                           dimnames = list(parm, names(conf_int)))
    }
    conf_int
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
    
    ## For glmer/glmerMod/lmer/lmerMod models, the variables of the fixed-effect
    ## terms are counted, the components of an interaction among them. The
    ## terms were formerly read by splitting the formula text at "+", which
    ## kept "hypertension * diabetes" as a single term, and after removing the
    ## random effects with a pattern that spanned parentheses, which reduced
    ## "log(age) + sex + (1 | site)" to "log"; both models were described as
    ## univariable. The terms of coxme models are read in the same way, from
    ## the formula list they keep: the formula of their call is only a name
    ## when they are fitted inside a function, so that the levels of each
    ## factor were counted as variables.
    if (inherits(model, c("glmerMod", "lmerMod", "merMod", "coxme"))) {
        term_labels <- model_term_labels(model)
        if (length(term_labels) > 0L) {
            unique_vars <- unique(unlist(strsplit(term_labels, ":", fixed = TRUE)))
            return(data.table::fifelse(length(unique_vars) == 1, "Univariable", "Multivariable"))
        }
    }
    
    ## Original code for other model types
    ## Get coefficient names once, without the backticks of non-syntactic
    ## variable names, which the names in xlevels do not carry
    coef_names <- gsub("`", "", names(stats::coef(model)), fixed = TRUE)
    
    ## Remove intercept from count
    if ("(Intercept)" %in% coef_names) {
        term_names <- coef_names[coef_names != "(Intercept)"]
    } else {
        term_names <- coef_names
    }
    
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
                    } else if (var %in% names(data_source) &&
                               is.logical(data_source[[var]])) {
                        xlevels[[var]] <- c("FALSE", "TRUE")
                    }
                }
            }
            if (length(xlevels) == 0) xlevels <- NULL
        }
    }
    
    ## Identify interaction terms, those containing ":" other than the
    ## coefficients of levels and the names of variables that contain one
    is_interaction <- is_interaction_term(term_names, xlevels,
                                          main_effect_labels(model))
    interaction_terms <- term_names[is_interaction]
    main_terms <- term_names[!is_interaction]

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
    
    ## A coefficient of another variable whose name begins with a factor's
    ## name (the continuous "smoking_years", or "smokerTRUE" of a logical
    ## "smoker", beside the factor "smoking") is not a level of that factor
    other_labels <- setdiff(gsub("`", "", model_term_labels(model), fixed = TRUE),
                            xlevel_names)

    matched_var <- vapply(main_terms, function(term) {
        hits <- xlevel_names[startsWith(term, xlevel_names)]
        hits <- hits[!vapply(hits, function(hit) {
            claimed_by_longer_term(term, hit, xlevels[[hit]], other_labels)
        }, logical(1))]
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
                      inverse.gaussian = "Inverse Gaussian",
                      quasibinomial = "Quasi-Binomial",
                      quasipoisson = "Quasi-Poisson",
                      paste0(toupper(substring(family, 1, 1)), substring(family, 2), " GLM")
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


#' Coefficients that belong to a variable whose name extends a factor's name
#'
#' Identifies the coefficients whose names begin with the name of a factor
#' but belong to another variable of the model whose name begins with it as
#' well, such as the continuous "smoking_years" or the coefficient
#' "smokerTRUE" of a logical variable "smoker" beside a factor "smoking". A
#' coefficient whose remainder after the factor's name is one of the
#' factor's levels is kept for the factor.
#'
#' @param terms Character vector of coefficient names, without backticks.
#' @param var Character string. Name of the factor.
#' @param var_levels Character vector of the levels of the factor.
#' @param labels Character vector of the term labels of the other
#'   variables of the model, without backticks.
#' @return Logical vector, \code{TRUE} for the coefficients that belong to a
#'   variable with a longer name.
#' @keywords internal
claimed_by_longer_term <- function(terms, var, var_levels, labels) {
    longer <- labels[nchar(labels) > nchar(var) & startsWith(labels, var)]
    if (length(longer) == 0L) {
        return(rep(FALSE, length(terms)))
    }
    claimed <- Reduce(`|`, lapply(longer, function(label) startsWith(terms, label)))
    claimed & !(substring(terms, nchar(var) + 1L) %in% var_levels)
}


#' Parse term into variable and group
#' 
#' Splits coefficient term names into base variable names and factor levels.
#' For example, "sexMale" becomes variable="sex" and group="Male".
#' Handles interaction terms and continuous variables appropriately.
#' 
#' A term that belongs to another variable of the model whose name begins
#' with a factor's name, such as a continuous variable "smoking_years" or the
#' coefficient "smokerTRUE" of a logical variable "smoker" beside a factor
#' "smoking", is not taken for a level of that factor (see
#' \code{claimed_by_longer_term()}).
#'
#' @param terms Character vector of coefficient term names.
#' @param xlevels Named list of factor levels from the model.
#' @param model Optional model object, from which the terms of the model are
#'   read and the factor levels of coxme models are reconstructed.
#' @return Data.table with 'variable' and 'group' columns.
#' @keywords internal
parse_term <- function(terms, xlevels = NULL, model = NULL) {
    n_terms <- length(terms)
    
    ## Initialize result vectors
    variable <- character(n_terms)
    group <- character(n_terms)
    
    ## Special handling for coxme models - reconstruct xlevels if needed
    if (!is.null(model) && inherits(model, "coxme") && is.null(xlevels)) {
        data_source <- attr(model, "data")
        if (is.null(data_source)) {
            data_source <- get_model_data(model)
        }
        
        if (!is.null(data_source)) {
            xlevels <- list()
            
            ## Extract factor structure from coefficient names, without the
            ## backticks of non-syntactic variable names
            coef_names <- gsub("`", "", names(coxme::fixef(model)), fixed = TRUE)
            
            ## Check each column in the data
            for (col_name in names(data_source)) {
                ## Check if any coefficient starts with this column name,
                ## compared literally rather than as a regular expression
                if (any(startsWith(coef_names, col_name))) {
                    ## This is a factor variable - get its levels
                    if (is.factor(data_source[[col_name]])) {
                        xlevels[[col_name]] <- levels(data_source[[col_name]])
                    } else if (is.character(data_source[[col_name]])) {
                        xlevels[[col_name]] <- sort(unique(data_source[[col_name]]))
                    } else if (is.logical(data_source[[col_name]])) {
                        xlevels[[col_name]] <- c("FALSE", "TRUE")
                    }
                }
            }
            
            if (length(xlevels) == 0) xlevels <- NULL
        }
    }
    
    ## Interactions should not be parsed as factor variables. A level that
    ## contains a colon ("shift08:00"), or the name of a variable that does,
    ## is not an interaction.
    is_interaction <- is_interaction_term(
        terms, xlevels, if (!is.null(model)) main_effect_labels(model))

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
        
        ## A coefficient of another variable whose name begins with a
        ## factor's name (the continuous "smoking_years", or "smokerTRUE" of a
        ## logical "smoker", beside the factor "smoking") is not a level of
        ## that factor
        other_labels <- if (is.null(model)) {
                            character(0)
                        } else {
                            setdiff(gsub("`", "", model_term_labels(model), fixed = TRUE),
                                    xlevel_names)
                        }

        ## For each factor variable, find all matching terms at once, skip interaction terms
        for (var in xlevel_names) {
            ## Find which terms start with this variable name
            matches <- startsWith(terms, var) & !is_interaction &  # Skip interactions
                !claimed_by_longer_term(terms, var, xlevels[[var]], other_labels)
            
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


#' Name the coefficients of the levels of a model's factors
#'
#' The coefficient of a factor level is named by the variable followed by
#' the level (\emph{e.g.,} "sexMale").
#'
#' @param xlevels Named list of factor levels, or \code{NULL}.
#' @return Character vector of the names that the coefficients of the
#'   levels in \code{xlevels} take.
#' @keywords internal
level_terms <- function(xlevels) {
    if (length(xlevels) == 0L) {
        return(character(0))
    }
    unlist(lapply(names(xlevels), function(var) paste0(var, xlevels[[var]])),
           use.names = FALSE)
}


#' Identify the coefficients of interactions
#'
#' R joins the components of an interaction coefficient with colons
#' (\emph{e.g.,} "treatmentDrug A:stageII"), but a factor level may contain
#' a colon as well (\emph{e.g.,} the level "08:00" of a factor \code{shift},
#' whose coefficient is "shift08:00"), as may the name of a variable. A name
#' with a colon is therefore taken for an interaction unless it is the
#' coefficient of a level in \code{xlevels} or the name of a main-effect
#' term of the model.
#'
#' @param terms Character vector of coefficient names, without the
#'   backticks of non-syntactic variable names.
#' @param xlevels Named list of factor levels, or \code{NULL}.
#' @param main_terms Character vector of the main-effect terms of the model,
#'   as returned by \code{main_effect_labels()}, or \code{NULL}.
#' @return Logical vector, \code{TRUE} for the coefficients of interactions.
#' @keywords internal
is_interaction_term <- function(terms, xlevels = NULL, main_terms = NULL) {
    grepl(":", terms, fixed = TRUE) &
        !terms %in% c(level_terms(xlevels), main_terms)
}


#' Name the main-effect terms of a model
#'
#' @param model Fitted model object.
#' @return Character vector of the labels of the terms of order one (the
#'   main effects), without the backticks of non-syntactic variable names.
#' @keywords internal
main_effect_labels <- function(model) {
    model_terms <- if (inherits(model, "coxme")) {
                       tryCatch(stats::terms(model$formulaList$fixed),
                                error = function(e) NULL)
                   } else {
                       tryCatch(stats::terms(model), error = function(e) NULL)
                   }
    labels <- attr(model_terms, "term.labels")
    term_order <- attr(model_terms, "order")
    if (is.null(labels) || length(term_order) != length(labels)) {
        return(character(0))
    }
    gsub("`", "", labels[term_order == 1L], fixed = TRUE)
}


#' Split an interaction coefficient into variables and levels
#'
#' Splits the name of an interaction coefficient (\emph{e.g.,}
#' "treatmentDrug A:stageII") at each colon and matches each component to
#' the variable it names. A component of a factor is the variable name
#' followed by one of its levels; a component of a continuous variable is
#' the variable name alone.
#'
#' The backticks that enclose non-syntactic variable names in coefficient
#' names (\emph{e.g.,} "`high bp`Yes") are removed first. Names are compared
#' literally, longest first, and a variable is preferred whose levels
#' include the remainder of the component, so that "sexualityHetero" is
#' assigned to \code{sexuality} rather than to \code{sex}. A level that
#' contains a colon (\emph{e.g.,} "08:00") is kept in one component when
#' it is among the levels in \code{xlevels}.
#'
#' @param term Character string. Name of an interaction coefficient.
#' @param var_names Character vector. Names of the variables to match.
#' @param xlevels Named list of factor levels, or \code{NULL}.
#' @return List with one element per component, each a list with elements
#'   \code{part} (the component), \code{var} (the matched variable, or
#'   \code{NA}) and \code{level} (the level, or \code{NA} for a continuous
#'   variable or an unmatched component).
#' @keywords internal
split_interaction_term <- function(term, var_names, xlevels = NULL) {
    parts <- strsplit(gsub("`", "", term, fixed = TRUE), ":", fixed = TRUE)[[1]]
    var_names <- as.character(var_names)
    var_names <- var_names[order(nchar(var_names), decreasing = TRUE)]

    ## Parts that a colon within a declared level split apart are joined
    ## again, the longest run that names a variable and one of its levels
    ## first
    declared_terms <- level_terms(xlevels)
    if (any(grepl(":", declared_terms, fixed = TRUE)) && length(parts) > 1L) {
        joined <- character(0)
        i <- 1L
        while (i <= length(parts)) {
            j <- i
            if (i < length(parts)) {
                for (end in seq.int(length(parts), i + 1L)) {
                    if (paste(parts[i:end], collapse = ":") %in% declared_terms) {
                        j <- end
                        break
                    }
                }
            }
            joined <- c(joined, paste(parts[i:j], collapse = ":"))
            i <- j + 1L
        }
        parts <- joined
    }

    lapply(parts, function(part) {
        if (part %in% var_names) {
            return(list(part = part, var = part, level = NA_character_))
        }

        candidates <- var_names[startsWith(part, var_names) &
                                nchar(part) > nchar(var_names)]
        if (length(candidates) == 0L) {
            return(list(part = part, var = NA_character_, level = NA_character_))
        }

        levels <- substring(part, nchar(candidates) + 1L)
        declared <- vapply(seq_along(candidates), function(k) {
            levels[k] %in% xlevels[[candidates[k]]]
        }, logical(1))
        k <- if (any(declared)) which(declared)[1L] else 1L

        list(part = part, var = candidates[k], level = levels[k])
    })
}


#' Read the event argument of a \code{Surv()} response
#'
#' The event of a \code{Surv()} response is its event argument: the second
#' argument of \code{Surv(time, event)}, the third of
#' \code{Surv(time, time2, event)}, or the argument named \code{event}.
#' \code{survival::clogit()} writes its response as
#' \code{Surv(rep(1, N), case)}, whose event is \code{case}. The argument
#' may be an expression, such as \code{status == 2}.
#'
#' @param response Left-hand side of a model formula, as a call, or a
#'   character string such as \code{"Surv(os_months, os_status)"}.
#' @return The expression of the event argument, or \code{NULL} when the
#'   response is not a call to \code{Surv()} with an event.
#' @keywords internal
surv_event_expression <- function(response) {
    if (is.character(response)) {
        response <- tryCatch(str2lang(response), error = function(e) NULL)
    }
    if (!is.call(response)) {
        return(NULL)
    }
    fn <- response[[1L]]
    is_surv <- identical(fn, as.name("Surv")) ||
        (is.call(fn) && identical(fn[[1L]], as.name("::")) &&
         identical(fn[[3L]], as.name("Surv")))
    if (!is_surv) {
        return(NULL)
    }
    args <- tryCatch(as.list(match.call(survival::Surv, response))[-1L],
                     error = function(e) NULL)
    if (!is.null(args$event)) args$event else args$time2
}


#' Code the events of a survival model
#'
#' Evaluates the event argument of the model's \code{Surv()} response (see
#' \code{surv_event_expression()}) on the rows of \code{data} and codes it
#' as \code{Surv()} codes it: 1 for an event and 0 for a censored
#' observation, whether the status is coded 0 and 1, \code{FALSE} and
#' \code{TRUE}, or 1 and 2, or is given as an expression such as
#' \code{status == 2}. A factor is coded by its second level. Summing the
#' variable named in the response instead counted a status coded 1 and 2
#' twice for each event and once for each censored observation, and
#' reading the second variable of the formula took the stop time of a
#' counting-process response and the first predictor of a
#' \code{clogit()} model.
#'
#' @param model Fitted \code{coxph}, \code{clogit}, or \code{coxme} model.
#' @param data Data frame in which the event argument is evaluated, such as
#'   the observations analyzed or the observations supplied.
#' @return Numeric vector of 0 and 1 (missing where the event is), one per
#'   row of \code{data}, or \code{NULL} when the events cannot be read.
#' @keywords internal
get_event_indicator <- function(model, data) {
    model_formula <- if (inherits(model, "coxme")) {
                         model$formulaList$fixed
                     } else {
                         tryCatch(stats::formula(model), error = function(e) NULL)
                     }
    if (is.null(model_formula) || length(model_formula) < 3L) {
        return(NULL)
    }
    surv_event_indicator(model_formula[[2L]], data, environment(model_formula))
}


#' Add the coded events of a survival model to its analysis rows
#'
#' The group counts of \code{m2dt()} sum an event column by group. The
#' events are coded from the model's \code{Surv()} response (see
#' \code{get_event_indicator()}) into a column of a copy of the model's
#' variables, so that the data supplied are not modified. Observations
#' whose event \code{Surv()} codes as missing, such as a status of 9 among
#' statuses of 0 and 1, did not enter the model and are not counted.
#'
#' @param model Fitted survival model.
#' @param model_class Character string of the model class.
#' @param data Data table of the analysis rows.
#' @return List with \code{data}, the data with an event column where the
#'   events could be coded, and \code{event_var}, the name of that column or
#'   \code{NULL}.
#' @keywords internal
add_event_column <- function(model, model_class, data) {
    event_var <- ".summata_event"
    if (is.null(data) || event_var %in% names(data)) {
        return(list(data = data,
                    event_var = if (!is.null(data)) event_var))
    }
    events <- get_event_indicator(model, data)
    if (is.null(events)) {
        return(list(data = data, event_var = NULL))
    }
    cols <- intersect(get_model_variables(model, model_class), names(data))
    keep <- !is.na(events)
    data <- data[keep, cols, with = FALSE]
    data.table::set(data, j = event_var, value = events[keep])
    list(data = data, event_var = event_var)
}


#' Code the events of a \code{Surv()} response
#'
#' @param response Left-hand side of a model formula, as a call, or a
#'   character string such as \code{"Surv(os_months, os_status)"}.
#' @param data Data frame in which the event argument is evaluated.
#' @param env Environment in which variables not in \code{data} are found.
#' @return See \code{get_event_indicator()}.
#' @keywords internal
surv_event_indicator <- function(response, data, env = parent.frame()) {
    event_expr <- surv_event_expression(response)
    if (is.null(event_expr) || is.null(data) || !inherits(data, "data.frame")) {
        return(NULL)
    }
    if (is.null(env)) {
        env <- parent.frame()
    }

    event <- tryCatch(eval(event_expr, data, env), error = function(e) NULL)
    if (is.null(event) || length(event) != nrow(data)) {
        return(NULL)
    }

    if (is.factor(event)) {
        return(as.numeric(as.integer(event) == 2L))
    }
    if (!is.numeric(event) && !is.logical(event)) {
        return(NULL)
    }

    ## Coded by Surv() itself, which converts values it does not accept to
    ## missing, as the model did when it was fitted
    status <- tryCatch(
        suppressWarnings(unclass(survival::Surv(rep(1, length(event)), event))[, "status"]),
        error = function(e) NULL
    )
    if (is.null(status)) {
        return(NULL)
    }
    as.numeric(status)
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

    ## The call of a coxme model fitted inside a function names a variable
    ## rather than the formula, which its formula list keeps
    fml <- if (inherits(model, "coxme") && !is.null(model$formulaList$fixed)) {
               model$formulaList$fixed
           } else {
               tryCatch(stats::formula(model), error = function(e) NULL)
           }

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
#' @param event_indicator Numeric vector of the events, one per row of
#'   \code{data}, coded 1 for an event and 0 otherwise, which is used in place
#'   of \code{event_var} when given. Survival models pass the events coded
#'   from their \code{Surv()} response (see \code{get_event_indicator()}).
#' @return Named list with \code{n_supplied}, \code{n_analyzed},
#'   \code{n_missing_outcome} and \code{n_missing_predictor}, and, where an
#'   event indicator was supplied, \code{events_supplied} and
#'   \code{events_analyzed}. \code{NULL} when the counts cannot be determined.
#' @keywords internal
get_analysis_counts <- function(data, outcome_vars, predictor_vars = NULL,
                                event_var = NULL, event_indicator = NULL) {

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
    if (!is.null(event_indicator) && length(event_indicator) != n_supplied) {
        event_indicator <- NULL
    }
    ## An event that Surv() codes as missing leaves the response missing
    if (!is.null(event_indicator)) {
        outcome_ok <- outcome_ok & !is.na(event_indicator)
    }
    if (is.null(event_indicator) && !is.null(event_var) && length(event_var) == 1 &&
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

    ## The events of a survival model are coded from its Surv() response
    event_indicator <- if (model_class %in% c("coxph", "clogit", "coxme")) {
                           get_event_indicator(model, data)
                       }

    return(get_analysis_counts(data, outcome_vars, predictor_vars, event_var,
                               event_indicator))
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
