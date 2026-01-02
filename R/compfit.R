#' Compare Multiple Regression Models
#'
#' Fits multiple regression models and provides a comprehensive comparison table
#' with model quality metrics, convergence diagnostics, and selection guidance.
#' Computes a composite score combining multiple quality metrics to facilitate 
#' rapid model comparison and selection.
#'
#' @param data A data.frame or data.table containing the dataset.
#' @param outcome Character string specifying the outcome variable. For survival
#'   analysis, use Surv() syntax (e.g., "Surv(time, status)").
#' @param model_list List of character vectors, each containing predictor names
#'   for one model. Can also be a single character vector to auto-generate nested models.
#' @param model_names Character vector of names for each model. If NULL, uses
#'   "Model 1", "Model 2", etc. Default is NULL.
#' @param interactions_list List of character vectors specifying interaction
#'   terms for each model. Each element corresponds to one model in model_list.
#'   Use NULL for models without interactions. Use colon notation for interactions
#'   (e.g., c("age:treatment")). If NULL, no interactions are added to any model.
#'   Default is NULL.
#' @param model_type Character string specifying model type. If "auto", detects
#'   based on outcome. Options: "auto", "glm", "lm", "coxph", "clogit".
#'   Default is "auto".
#' @param family For GLM models, the error distribution family. Default is "binomial".
#' @param conf_level Numeric confidence level for intervals. Default is 0.95.
#' @param p_digits Integer specifying the number of decimal places for p-values.
#'   P-values smaller than \code{10^(-p_digits)} are displayed as "< 0.001"
#'   (for \code{p_digits = 3}), "< 0.0001" (for \code{p_digits = 4}), etc.
#'   Default is 3.
#' @param include_coefficients Logical. If TRUE, includes a second table with
#'   coefficient estimates. Default is FALSE.
#' @param scoring_weights Named list of scoring weights. Each weight should be
#'   between 0 and 1, and they should sum to 1. Available metrics depend on model
#'   type. If NULL, uses sensible defaults. See Details for available metrics.
#' @param labels Named character vector for custom variable labels. Default is NULL.
#' @param ... Additional arguments passed to model fitting functions.
#'
#' @return A data.table with class "compfit_result" containing:
#'   \describe{
#'     \item{Model}{Model name/identifier}
#'     \item{Summata Score}{Composite metric for model selection (higher is better)}
#'     \item{N}{Sample size}
#'     \item{Events}{Number of events (for survival/logistic)}
#'     \item{Predictors}{Number of predictors}
#'     \item{Converged}{Whether model converged properly}
#'     \item{AIC}{Akaike Information Criterion}
#'     \item{BIC}{Bayesian Information Criterion}
#'     \item{Pseudo-R²}{McFadden's pseudo-R-squared (GLM)}
#'     \item{Concordance}{C-statistic (logistic/survival)}
#'     \item{Brier Score}{Brier accuracy score (logistic)}
#'     \item{Global p}{Overall model p-value}
#'   }
#'   
#'   Attributes include:
#'   \describe{
#'     \item{models}{List of fitted model objects}
#'     \item{coefficients}{Coefficient comparison table (if requested)}
#'     \item{best_model}{Name of recommended model}
#'   }
#'
#' @details
#' This function fits all specified models and computes comprehensive quality
#' metrics for comparison. It generates a composite score (Summata Score) that 
#' combines multiple metrics: lower AIC/BIC (information criteria), higher 
#' concordance (discrimination), and model convergence status.
#' 
#' For GLMs, McFadden's pseudo-R-squared is calculated as 1 - (logLik/logLik_null).
#' For survival models, the global p-value comes from the log-rank test.
#' 
#' Models that fail to converge are flagged and penalized in the composite score.
#' 
#' \strong{Interaction Terms:}
#' 
#' When \code{interactions_list} is provided, each element specifies the
#' interaction terms for the corresponding model in \code{model_list}. This is
#' particularly useful for testing whether adding interactions improves model fit:
#' \itemize{
#'   \item Use NULL for models without interactions
#'   \item Specify interactions using colon notation: c("age:treatment", "sex:stage")
#'   \item Main effects for all variables in interactions must be in the predictor list
#'   \item Common pattern: Compare main effects model vs model with interactions
#' }
#' 
#' Scoring weights can be customized based on model type:
#' \itemize{
#'   \item GLM: "convergence", "aic", "concordance", "pseudo_r2", "brier" 
#'   \item Cox: "convergence", "aic", "concordance", "global_p"
#'   \item Linear: "convergence", "aic", "pseudo_r2", "rmse"
#' }
#' Default weights emphasize discrimination (concordance) and model fit (AIC).
#'
#' The composite score is designed as a tool to quickly rank models by their 
#' quality metrics. It should be used alongside traditional model selection 
#' criteria rather than as a definitive selection method.
#'
#' @seealso
#' \code{\link{fit}} for individual model fitting,
#' \code{\link{fullfit}} for automated variable selection,
#' \code{\link{table2pdf}} for exporting results
#'
#' @examples
#' # Load example data
#' data(clintrial)
#' data(clintrial_labels)
#' 
#' # Example 1: Compare nested logistic regression models
#' models <- list(
#'     base = c("age", "sex"),
#'     clinical = c("age", "sex", "smoking", "diabetes"),
#'     full = c("age", "sex", "smoking", "diabetes", "stage", "ecog")
#' )
#' 
#' comparison <- compfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     model_list = models,
#'     model_names = c("Base", "Clinical", "Full")
#' )
#' comparison
#' 
#' # Example 2: Compare Cox survival models
#' library(survival)
#' surv_models <- list(
#'     simple = c("age", "sex"),
#'     clinical = c("age", "sex", "stage", "grade")
#' )
#' 
#' surv_comparison <- compfit(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     model_list = surv_models,
#'     model_type = "coxph"
#' )
#' surv_comparison
#' 
#' # Example 3: Test effect of adding interaction terms
#' interaction_models <- list(
#'     main = c("age", "treatment", "sex"),
#'     interact = c("age", "treatment", "sex")
#' )
#' 
#' interaction_comp <- compfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     model_list = interaction_models,
#'     model_names = c("Main Effects", "With Interaction"),
#'     interactions_list = list(
#'         NULL,
#'         c("treatment:sex")
#'     )
#' )
#' interaction_comp
#' 
#' # Example 4: Include coefficient comparison table
#' detailed <- compfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     model_list = models,
#'     include_coefficients = TRUE,
#'     labels = clintrial_labels
#' )
#' 
#' # Access coefficient table
#' coef_table <- attr(detailed, "coefficients")
#' coef_table
#' 
#' # Example 5: Access fitted model objects
#' fitted_models <- attr(comparison, "models")
#' names(fitted_models)
#' 
#' # Example 6: Get best model recommendation
#' best <- attr(comparison, "best_model")
#' cat("Recommended model:", best, "\n")
#'
#' @export
compfit <- function(data,
                    outcome,
                    model_list,
                    model_names = NULL,
                    interactions_list = NULL,
                    model_type = "auto",
                    family = "binomial",
                    conf_level = 0.95,
                    p_digits = 3,
                    include_coefficients = FALSE,
                    scoring_weights = NULL,
                    labels = NULL,
                    ...) {
    
    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }

    ## Check if any model has random effects (for auto-detection)
    ## Random effects are specified with | syntax, e.g., "(1|site)"
    has_any_random_effects <- any(sapply(model_list, function(preds) {
        any(grepl("\\|", preds))
    }))

    ## Auto-detect model type if requested
    if (model_type == "auto") {
        if (grepl("^Surv\\(", outcome)) {
            if (has_any_random_effects) {
                model_type <- "coxme"
                message("Auto-detected survival outcome with random effects, using coxme")
            } else {
                model_type <- "coxph"
                message("Auto-detected survival outcome, using Cox regression")
            }
        } else if (is.factor(data[[outcome]]) || length(unique(data[[outcome]])) == 2) {
            if (has_any_random_effects) {
                model_type <- "glmer"
                family <- "binomial"
                message("Auto-detected binary outcome with random effects, using glmer")
            } else {
                model_type <- "glm"
                family <- "binomial"
                message("Auto-detected binary outcome, using logistic regression")
            }
        } else if (is.numeric(data[[outcome]])) {
            if (has_any_random_effects) {
                model_type <- "lmer"
                message("Auto-detected continuous outcome with random effects, using lmer")
            } else {
                model_type <- "lm"
                message("Auto-detected continuous outcome, using linear regression")
            }
        } else {
            stop("Cannot auto-detect model type for outcome: ", outcome)
        }
    }

    ## Set model names if not provided
    if (is.null(model_names)) {
        ## Use list names if available, otherwise generate defaults
        if (!is.null(names(model_list)) && all(nzchar(names(model_list)))) {
            model_names <- names(model_list)
        } else {
            model_names <- paste("Model", seq_along(model_list))
        }
    }
    
    ## Validate interactions_list if provided
    if (!is.null(interactions_list)) {
        if (!is.list(interactions_list)) {
            stop("interactions_list must be a list")
        }
        if (length(interactions_list) != length(model_list)) {
            stop("interactions_list must have the same length as model_list (", 
                 length(model_list), " models)")
        }
        ## Check each element is either NULL or character vector
        for (i in seq_along(interactions_list)) {
            if (!is.null(interactions_list[[i]])) {
                if (!is.character(interactions_list[[i]])) {
                    stop("Each element of interactions_list must be NULL or a character vector")
                }
            }
        }
    }

    ## Initialize results storage
    n_models <- length(model_list)
    comparison_list <- vector("list", n_models)
    models <- vector("list", n_models)
    names(models) <- model_names
    coef_results <- list()

    ## Fit each model
    for (i in seq_along(model_list)) {
        ## Get interaction count for message
        n_interact <- if (!is.null(interactions_list) && !is.null(interactions_list[[i]])) length(interactions_list[[i]]) else 0
        if (n_interact > 0) {
            message(sprintf("Fitting %s with %d predictors + %d interaction%s...", model_names[i], length(model_list[[i]]), n_interact, if(n_interact > 1) "s" else ""))
        } else {
            message(sprintf("Fitting %s with %d predictors...", model_names[i], length(model_list[[i]])))
        }
        
        ## Extract interaction terms for this model BEFORE tryCatch
        current_interactions <- if (!is.null(interactions_list)) {
                                    interactions_list[[i]]
                                } else {
                                    NULL
                                }
        
        ## Fit model using fit() for consistency
        model_result <- tryCatch({
            fit(data = data,
                outcome = outcome,
                interactions = current_interactions,
                predictors = model_list[[i]],
                model_type = model_type,
                family = family,
                conf_level = conf_level,
                labels = labels,
                keep_qc_stats = TRUE,
                ...)
        }, error = function(e) {
            message("  Warning - model failed to fit: ", e$message)
            NULL
        })

        ## Build comparison row - conditional on model type
        if (is.null(model_result)) {
            ## Model failed to fit
            row <- data.table::data.table(
                                   Model = model_names[i],
                                   N = nrow(data),
                                   Events = NA_integer_,
                                   Predictors = length(model_list[[i]]),
                                   Converged = "Failed",
                                   AIC = NA_real_,
                                   BIC = NA_real_,
                                   `Pseudo-R^2` = NA_real_,
                                   Concordance = NA_real_,
                                   `Global p` = NA_real_
                               )
            ## Add model-type specific columns for failed models
            if (model_type == "glm") {
                row$`Brier Score` <- NA_real_
            } else if (model_type %in% c("lmer", "glmer", "coxme")) {
                row$Groups <- NA_integer_
                row$`Marginal R^2` <- NA_real_
                row$`Conditional R^2` <- NA_real_
                row$ICC <- NA_real_
                if (model_type == "glmer") {
                    row$`Brier Score` <- NA_real_
                }
            }
        } else {

            ## Extract model and raw data
            model <- attr(model_result, "model")
            raw_data <- attr(model_result, "raw_data")
            models[[model_names[i]]] <- model
            
            ## Store coefficients if requested
            if (include_coefficients) {
                coef_results[[model_names[i]]] <- model_result
            }
            
            ## Check convergence
            converged <- check_convergence(model)
            
            ## Extract metrics
            metrics <- extract_model_metrics(model, raw_data, model_type)

            ## Build comparison row
            row <- data.table::data.table(
                                   Model = model_names[i],
                                   N = metrics$n,
                                   Events = metrics$events,
                                   Predictors = length(model_list[[i]]),
                                   Converged = converged,
                                   AIC = if (!is.null(metrics$aic) && !is.na(metrics$aic)) round(metrics$aic, 1) else NA_real_,
                                   BIC = if (!is.null(metrics$bic) && !is.na(metrics$bic)) round(metrics$bic, 1) else NA_real_,
                                   `Pseudo-R^2` = if (!is.null(metrics$pseudo_r2) && !is.na(metrics$pseudo_r2)) round(metrics$pseudo_r2, 3) else NA_real_,
                                   Concordance = if (!is.null(metrics$concordance) && !is.na(metrics$concordance)) round(metrics$concordance, 3) else NA_real_,
                                   `Global p` = if (!is.null(metrics$global_p) && !is.na(metrics$global_p)) format_pvalues_fit(metrics$global_p, p_digits) else NA_character_
                               )
            
            ## Add model-type specific columns
            if (model_type == "glm") {
                row$`Brier Score` <- if (!is.null(metrics$brier_score)) round(metrics$brier_score, 3) else NA_real_
            } else if (model_type %in% c("lmer", "glmer", "coxme")) {
                ## Add mixed-effects specific columns
                row$Groups <- if (!is.null(metrics$n_groups)) as.integer(metrics$n_groups) else NA_integer_
                row$`Marginal R^2` <- if (!is.null(metrics$marginal_r2) && !is.na(metrics$marginal_r2)) round(metrics$marginal_r2, 3) else NA_real_
                row$`Conditional R^2` <- if (!is.null(metrics$conditional_r2) && !is.na(metrics$conditional_r2)) round(metrics$conditional_r2, 3) else NA_real_
                row$ICC <- if (!is.null(metrics$icc) && !is.na(metrics$icc)) round(metrics$icc, 3) else NA_real_
                
                if (model_type == "glmer") {
                    row$`Brier Score` <- if (!is.null(metrics$brier_score)) round(metrics$brier_score, 3) else NA_real_
                }
            }
        }
        
        ## Store row in pre-allocated list
        comparison_list[[i]] <- row
    }

    ## Combine all rows using rbindlist
    comparison <- data.table::rbindlist(comparison_list, fill = TRUE)

    ## Calculate scores and sort
    comparison <- calculate_model_scores(comparison, model_type, scoring_weights)

    ## Rename columns from ASCII to Unicode superscripts for display
    col_renames <- c(
        "Pseudo-R^2" = "Pseudo-R\u00b2",
        "Marginal R^2" = "Marginal R\u00b2",
        "Conditional R^2" = "Conditional R\u00b2",
        "Adjusted R^2" = "Adjusted R\u00b2"
    )
    for (old_name in names(col_renames)) {
        if (old_name %in% names(comparison)) {
            data.table::setnames(comparison, old_name, col_renames[[old_name]])
        }
    }

    ## Format for display (ensure proper column order with Score always last)
    ## Get all current column names
    all_cols <- names(comparison)
    
    ## Define the preferred order for known columns (Score will be moved to end)
    if (model_type == "glm") {
        preferred_order <- c("Model", "N", "Events", "Predictors", "Converged",
                             "AIC", "BIC", "Pseudo-R\u00b2", "Concordance", "Brier Score",
                             "Global p")
    } else if (model_type %in% c("lmer", "glmer")) {
        preferred_order <- c("Model", "N", "Events", "Predictors", "Groups", "Converged",
                             "AIC", "BIC", "Pseudo-R\u00b2", "Marginal R\u00b2", "Conditional R\u00b2",
                             "ICC", "Concordance", "Global p")
    } else if (model_type == "coxme") {
        preferred_order <- c("Model", "N", "Events", "Predictors", "Groups", "Converged",
                             "AIC", "BIC", "Concordance", "Global p")
    } else {
        preferred_order <- c("Model", "N", "Events", "Predictors", "Converged",
                             "AIC", "BIC", "Pseudo-R\u00b2", "Concordance", "Global p")
    }
    
    ## Build final column order: Model, Score, then preferred columns, then extras
    existing_preferred <- intersect(preferred_order, all_cols)
    existing_preferred <- setdiff(existing_preferred, "Model")
    other_cols <- setdiff(all_cols, c(preferred_order, "Summata Score"))
    final_order <- c("Model", "Summata Score", existing_preferred, other_cols)
    
    ## Only reorder columns that actually exist
    final_order <- intersect(final_order, all_cols)
    data.table::setcolorder(comparison, final_order)
    
    ## Attach attributes
    data.table::setattr(comparison, "models", models)
    data.table::setattr(comparison, "model_type", model_type)
    data.table::setattr(comparison, "outcome", outcome)

    if (include_coefficients && length(coef_results) > 0) {
        coef_table <- combine_coefficient_tables(coef_results, model_names)
        data.table::setattr(comparison, "coefficients", coef_table)
    }

    ## Identify best model (first one after sorting)
    if (nrow(comparison) > 0) {
        data.table::setattr(comparison, "best_model", comparison$Model[1])
    }

    class(comparison) <- c("compfit_result", class(comparison))

    return(comparison)
}

#' Print method showing scoring methodology
#' @keywords internal
#' @export
print.compfit_result <- function(x, ...) {
    cat("\nModel Comparison Results\n")
    cat("Outcome: ", attr(x, "outcome"), "\n", sep = "")
    cat("Model Type: ", attr(x, "model_type"), "\n", sep = "")
    
    ## Show scoring weights
    weights <- attr(x, "weights")
    if (!is.null(weights)) {
        cat("\nSummata Score Weights:\n")
        for (metric in names(weights)) {
            ## Formatting for metric names
            display_name <- switch(metric,
                                   "convergence" = "Convergence",
                                   "aic" = "AIC",
                                   "bic" = "BIC",
                                   "concordance" = "Concordance",
                                   "c_stat" = "C-statistic",
                                   "c_index" = "C-index",
                                   "pseudo_r2" = "Pseudo-R\u00b2",
                                   "global_p" = "Global p-value",
                                   "rmse" = "RMSE",
                                   "brier" = "Brier score",
                                   "adj_r2" = "Adjusted R\u00b2",
                                   "marginal_r2" = "Marginal R\u00b2",
                                   "conditional_r2" = "Conditional R\u00b2",
                                   "icc" = "ICC",
                                   metric  # Generic fallback
                                   )
            cat(sprintf("  %s: %.0f%%\n", display_name, weights[[metric]] * 100))
        }
    }
    
    ## Identify best model
    if (nrow(x) > 0) {
        best_model <- x$Model[1]  # Already sorted by score
        cat("\nRecommended Model: ", best_model, " (Summata Score: ", x$`Summata Score`[1], ")\n", sep = "")
    }
    
    cat("\nModels ranked by selection score:\n")
    NextMethod("print", x)

    cat("\nSummata Score interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor\n")
    
    if (!is.null(attr(x, "coefficients"))) {
        cat("\nNote: Coefficient comparison available via attr(result, 'coefficients')\n")
    }
    
    invisible(x)
}
