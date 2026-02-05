# Column name constants for Unicode superscripts
# Defined here to ensure consistent column naming throughout
.pseudo_r2_col <- "Pseudo-R\u00b2"
.marginal_r2_col <- "Marginal R\u00b2"
.conditional_r2_col <- "Conditional R\u00b2"
.adjusted_r2_col <- "Adjusted R\u00b2"

#' Check model convergence
#' 
#' Checks convergence status for various model types including standard
#' regression models and mixed-effects models.
#' 
#' @param model Fitted model object.
#' @return Character string: "Yes", "No", "Suspect", or "Failed"
#' @keywords internal
check_convergence <- function(model) {
    
    if (inherits(model, "glm")) {
        if (!model$converged) return("No")
        if (any(abs(coef(model)) > 10)) return("Suspect")
        return("Yes")
        
    } else if (inherits(model, c("coxph", "coxme"))) {
        
        ## Cox models do not always have an iter field
        ## Check for actual convergence issues
        if (!is.null(model$info) && !is.null(model$info$convergence)) {
            if (!model$info$convergence) return("No")
        }
        
        ## Check for extreme coefficients as a proxy for issues
        if (any(abs(coef(model)) > 10, na.rm = TRUE)) return("Suspect")
        
        ## Check if model failed (no coefficients)
        if (length(coef(model)) == 0) return("Failed")
        
        return("Yes")
        
    } else if (inherits(model, c("lmerMod", "glmerMod", "merMod"))) {
        ## Mixed-effects models from lme4
        
        ## Check for convergence warnings in the model object
        ## lme4 stores convergence info in the optinfo slot
        if (!is.null(model@optinfo$conv$lme4)) {
            conv_info <- model@optinfo$conv$lme4
            if (!is.null(conv_info$code) && conv_info$code != 0) {
                return("No")
            }
            if (!is.null(conv_info$messages) && length(conv_info$messages) > 0) {
                ## Check for specific convergence warnings
                msgs <- tolower(paste(conv_info$messages, collapse = " "))
                if (grepl("failed to converge|singular|gradient|hessian", msgs)) {
                    return("Suspect")
                }
            }
        }
        
        ## Check for singular fit (boundary estimates)
        if (lme4::isSingular(model)) {
            return("Suspect")
        }
        
        ## Check for extreme fixed effects coefficients
        fixed_coefs <- lme4::fixef(model)
        if (any(abs(fixed_coefs) > 10, na.rm = TRUE)) {
            return("Suspect")
        }
        
        ## Check that random effects variance is estimable
        var_comps <- lme4::VarCorr(model)
        for (i in seq_along(var_comps)) {
            var_vals <- diag(as.matrix(var_comps[[i]]))
            if (any(var_vals < 1e-10) || any(is.na(var_vals))) {
                return("Suspect")
            }
        }
        
        return("Yes")
        
    } else {
        
        return("Yes")
        
    }
}

#' Extract comprehensive model metrics based on academic consensus
#' 
#' Extracts quality control metrics from fitted models for comparison.
#' Supports GLM, Cox, linear, and mixed-effects models.
#' 
#' @param model Fitted model object.
#' @param raw_data Data.table with raw model information.
#' @param model_type Character string indicating model type.
#' @return Named list of metrics.
#' @keywords internal
extract_model_metrics <- function(model, raw_data, model_type) {
    
    metrics <- list(
        n = raw_data$n[1],
        events = if ("events" %chin% names(raw_data)) raw_data$events[1] else NA_integer_,
        predictors = length(coef(model)),
        aic = if ("AIC" %chin% names(raw_data)) raw_data$AIC[1] else AIC(model),
        bic = if ("BIC" %chin% names(raw_data)) raw_data$BIC[1] else BIC(model),
        ## Initialize all possible metrics as NA to avoid missing list elements
        pseudo_r2 = NA_real_,
        concordance = NA_real_,
        global_p = NA_real_,
        null_deviance = NA_real_,
        residual_deviance = NA_real_,
        deviance_ratio = NA_real_,
        mcfadden_r2 = NA_real_,
        nagelkerke_r2 = NA_real_,
        tjur_r2 = NA_real_,
        c_statistic = NA_real_,
        brier_score = NA_real_,
        hoslem_p = NA_real_,
        ## Mixed-effects-specific metrics
        marginal_r2 = NA_real_,
        conditional_r2 = NA_real_,
        icc = NA_real_,
        n_groups = NA_integer_,
        sigma = NA_real_
    )
    
    ## GLM-specific metrics
    if (model_type == "glm") {
        ## Deviance metrics
        metrics$null_deviance <- model$null.deviance
        metrics$residual_deviance <- model$deviance
        metrics$deviance_ratio <- 1 - (model$deviance / model$null.deviance)
        
        if (model$family$family == "binomial") {
            ## Multiple pseudo R-squared measures
            n <- length(model$y)
            
            ## McFadden R^2
            null_model <- stats::glm(model$y ~ 1, family = model$family)
            metrics$mcfadden_r2 <- as.numeric(1 - (logLik(model)/logLik(null_model)))
            metrics$pseudo_r2 <- metrics$mcfadden_r2  # Set generic pseudo_r2
            
            ## Nagelkerke R^2 (Cox-Snell adjusted)
            cox_snell <- 1 - exp((model$deviance - model$null.deviance)/n)
            metrics$nagelkerke_r2 <- cox_snell / (1 - exp(-model$null.deviance/n))
            
            ## Tjur R^2 (coefficient of discrimination)
            pred_probs <- fitted(model)
            metrics$tjur_r2 <- mean(pred_probs[model$y == 1]) - mean(pred_probs[model$y == 0])
            
            ## C-statistic (concordance)
            if (requireNamespace("pROC", quietly = TRUE)) {
                roc_obj <- pROC::roc(model$y, fitted(model), quiet = TRUE)
                metrics$c_statistic <- as.numeric(pROC::auc(roc_obj))
                metrics$concordance <- metrics$c_statistic  # Set generic concordance
                
                ## Optimal threshold metrics
                coords <- pROC::coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
                metrics$optimal_threshold <- coords$threshold
                metrics$sensitivity <- coords$sensitivity
                metrics$specificity <- coords$specificity
            }
            
            ## Brier score
            metrics$brier_score <- mean((fitted(model) - model$y)^2)
            
            ## Hosmer-Lemeshow test (if ResourceSelection available)
            if (requireNamespace("ResourceSelection", quietly = TRUE) && n >= 40) {
                tryCatch({
                    hl <- ResourceSelection::hoslem.test(model$y, fitted(model), g = 10)
                    metrics$hoslem_p <- hl$p.value
                }, error = function(e) {
                    metrics$hoslem_p <- NA_real_
                })
            }
            
            ## Global LR test
            lr_stat <- model$null.deviance - model$deviance
            df <- model$df.null - model$df.residual
            metrics$global_p <- pchisq(lr_stat, df, lower.tail = FALSE)
        }
        
    } else if (model_type == "coxph") {
        
        ## Initialize Cox-specific metrics
        metrics$c_index <- NA_real_
        metrics$rsq <- NA_real_
        metrics$lr_test_p <- NA_real_
        
        ## Cache summary to avoid repeated calls
        summ <- summary(model)
        
        ## Concordance
        if (!is.null(model$concordance)) {
            metrics$c_index <- model$concordance["concordance"]
            metrics$concordance <- metrics$c_index
            
        } else if (!is.null(summ$concordance)) {
            metrics$c_index <- summ$concordance["C"]
            metrics$concordance <- metrics$c_index
        }
        
        ## R-squared (different location in Cox models)
        if (!is.null(summ$rsq)) {
            metrics$rsq <- summ$rsq[1]
            metrics$pseudo_r2 <- metrics$rsq
        }
        
        ## Global tests
        if (!is.null(summ$logtest)) {
            metrics$lr_test_p <- summ$logtest["pvalue"]
            metrics$global_p <- metrics$lr_test_p
            
        } else if (!is.null(summ$waldtest)) {
            metrics$global_p <- summ$waldtest["pvalue"]
        }
        
    } else if (model_type == "lm") {
        
        ## Linear model metrics
        summ <- summary(model)
        metrics$pseudo_r2 <- summ$r.squared
        metrics$adj_r2 <- summ$adj.r.squared
        metrics$sigma <- summ$sigma
        metrics$rmse <- sqrt(mean(residuals(model)^2))
        
        ## F-test for global significance
        if (!is.null(summ$fstatistic)) {
            f_stat <- summ$fstatistic
            metrics$global_p <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
        }
        
    } else if (model_type == "lmer" || inherits(model, "lmerMod")) {
        
        ## Linear mixed-effects model metrics
        metrics <- extract_lmer_metrics(model, raw_data, metrics)
        
    } else if (model_type == "glmer" || inherits(model, "glmerMod")) {
        
        ## Generalized linear mixed-effects model metrics
        metrics <- extract_glmer_metrics(model, raw_data, metrics)
        
    } else if (model_type == "coxme") {
        
        ## Mixed-effects Cox model metrics
        metrics <- extract_coxme_metrics(model, raw_data, metrics)
    }
    
    return(metrics)
}


#' Extract metrics for linear mixed-effects models (lmer)
#' 
#' Extracts quality metrics specific to linear mixed-effects models including
#' R-squared measures, ICC, and global significance tests.
#' 
#' @param model Fitted lmerMod object from lme4.
#' @param raw_data Data.table with raw model information.
#' @param metrics Named list of initialized metrics to populate.
#' @return Updated metrics list with lmer-specific values.
#' @keywords internal
extract_lmer_metrics <- function(model, raw_data, metrics) {
    
    ## Number of observations and groups
    metrics$n <- nobs(model)
    n_grps <- lme4::ngrps(model)
    metrics$n_groups <- sum(n_grps)
    
    ## Fixed effects count (excluding intercept for predictor count)
    fixed_coefs <- lme4::fixef(model)
    metrics$predictors <- length(fixed_coefs) - 1  # Exclude intercept
    
    ## Residual standard deviation
    metrics$sigma <- sigma(model)
    
    ## RMSE (same as sigma for linear mixed models)
    metrics$rmse <- metrics$sigma
    
    ## R-squared measures using MuMIn if available
    if (requireNamespace("MuMIn", quietly = TRUE)) {
        r2_vals <- tryCatch({
            MuMIn::r.squaredGLMM(model)
        }, error = function(e) NULL)
        
        if (!is.null(r2_vals) && nrow(r2_vals) > 0) {
            metrics$marginal_r2 <- r2_vals[1, "R2m"]
            metrics$conditional_r2 <- r2_vals[1, "R2c"]
            metrics$pseudo_r2 <- metrics$marginal_r2  # Use marginal for comparison
        }
    }
    
    ## Alternative R-squared calculation if MuMIn not available
    if (is.na(metrics$pseudo_r2)) {
        ## Simple variance-based R^2
        var_fixed <- var(predict(model, re.form = NA))  # Fixed effects only
        var_total <- var(model@resp$y)
        metrics$pseudo_r2 <- var_fixed / var_total
        metrics$marginal_r2 <- metrics$pseudo_r2
    }
    
    ## ICC (Intraclass Correlation Coefficient)
    ## ICC = variance of random effects / total variance
    var_comps <- as.data.frame(lme4::VarCorr(model))
    random_var <- sum(var_comps$vcov[var_comps$grp != "Residual"])
    residual_var <- var_comps$vcov[var_comps$grp == "Residual"]
    if (!is.na(residual_var) && (random_var + residual_var) > 0) {
        metrics$icc <- random_var / (random_var + residual_var)
    }
    
    ## Global test using likelihood ratio test
    ## Compare to null model (intercept + random effects only)
    null_formula <- stats::update(stats::formula(model), . ~ 1 + (1 | .))
    
    ## Try to fit null model and compute LR test
    tryCatch({
        ## Get random effects structure
        re_terms <- lme4::findbars(formula(model))
        if (length(re_terms) > 0) {
            ## Build null formula with random effects
            re_string <- paste(sapply(re_terms, deparse), collapse = " + ")
            null_form <- stats::as.formula(paste(". ~ 1 +", re_string))
            null_model <- stats::update(model, formula = null_form)
            
            ## LR test
            lr_stat <- as.numeric(-2 * (logLik(null_model) - logLik(model)))
            df_diff <- attr(logLik(model), "df") - attr(logLik(null_model), "df")
            if (df_diff > 0) {
                metrics$global_p <- pchisq(lr_stat, df_diff, lower.tail = FALSE)
            }
        }
    }, error = function(e) {
        ## If null model fails, use Wald test approximation
        summ <- summary(model)
        coef_tab <- coef(summ)
        if (nrow(coef_tab) > 1) {
            ## Use minimum p-value as conservative estimate
            ## (excluding intercept)
            if (ncol(coef_tab) >= 5) {
                p_vals <- coef_tab[-1, ncol(coef_tab)]
                metrics$global_p <- min(p_vals, na.rm = TRUE)
            }
        }
    })
    
    return(metrics)
}


#' Extract metrics for generalized linear mixed-effects models (glmer)
#' 
#' Extracts quality metrics specific to generalized linear mixed-effects models
#' including concordance, R-squared measures, ICC, and Brier score for binomial.
#' 
#' @param model Fitted glmerMod object from lme4.
#' @param raw_data Data.table with raw model information.
#' @param metrics Named list of initialized metrics to populate.
#' @return Updated metrics list with glmer-specific values.
#' @keywords internal
extract_glmer_metrics <- function(model, raw_data, metrics) {
    
    summ <- summary(model)
    
    ## Number of observations and groups
    metrics$n <- nobs(model)
    n_grps <- lme4::ngrps(model)
    metrics$n_groups <- sum(n_grps)
    
    ## Fixed effects count
    fixed_coefs <- lme4::fixef(model)
    metrics$predictors <- length(fixed_coefs) - 1  # Exclude intercept
    
    ## For binomial models, get events
    if (summ$family == "binomial") {
        response_var <- model@resp$y
        if (!is.null(response_var)) {
            metrics$events <- sum(response_var, na.rm = TRUE)
        }
    }
    
    ## R-squared measures using MuMIn if available
    if (requireNamespace("MuMIn", quietly = TRUE)) {
        r2_vals <- tryCatch({
            MuMIn::r.squaredGLMM(model)
        }, error = function(e) NULL)
        
        if (!is.null(r2_vals) && nrow(r2_vals) > 0) {
            metrics$marginal_r2 <- r2_vals[1, "R2m"]
            metrics$conditional_r2 <- r2_vals[1, "R2c"]
            metrics$pseudo_r2 <- metrics$marginal_r2
        }
    }
    
    ## C-statistic for binomial models
    if (summ$family == "binomial") {
        if (requireNamespace("pROC", quietly = TRUE)) {
            c_stat <- tryCatch({
                fitted_vals <- fitted(model)
                response_vals <- model@resp$y
                roc_obj <- pROC::roc(response_vals, fitted_vals, quiet = TRUE)
                as.numeric(pROC::auc(roc_obj))
            }, error = function(e) NA_real_)
            
            if (!is.na(c_stat)) {
                metrics$c_statistic <- c_stat
                metrics$concordance <- c_stat
            }
        }
        
        ## Brier score
        fitted_vals <- fitted(model)
        response_vals <- model@resp$y
        metrics$brier_score <- mean((fitted_vals - response_vals)^2)
    }
    
    ## ICC for random effects
    var_comps <- as.data.frame(lme4::VarCorr(model))
    random_var <- sum(var_comps$vcov[var_comps$grp != "Residual"])
    
    ## For binomial, use pi^2/3 as residual variance on latent scale
    if (summ$family == "binomial") {
        latent_var <- pi^2 / 3
        metrics$icc <- random_var / (random_var + latent_var)
    }
    
    ## Global test using LR test
    tryCatch({
        re_terms <- lme4::findbars(formula(model))
        if (length(re_terms) > 0) {
            re_string <- paste(sapply(re_terms, deparse), collapse = " + ")
            null_form <- stats::as.formula(paste(". ~ 1 +", re_string))
            null_model <- stats::update(model, formula = null_form)
            
            lr_stat <- as.numeric(-2 * (logLik(null_model) - logLik(model)))
            df_diff <- attr(logLik(model), "df") - attr(logLik(null_model), "df")
            if (df_diff > 0) {
                metrics$global_p <- pchisq(lr_stat, df_diff, lower.tail = FALSE)
            }
        }
    }, error = function(e) {
        ## Fallback to Wald test
        coef_tab <- coef(summ)
        if (nrow(coef_tab) > 1 && ncol(coef_tab) >= 4) {
            p_vals <- coef_tab[-1, ncol(coef_tab)]
            metrics$global_p <- min(p_vals, na.rm = TRUE)
        }
    })
    
    return(metrics)
}


#' Extract metrics for mixed-effects Cox models (coxme)
#' 
#' Extracts quality metrics specific to mixed-effects Cox proportional hazards
#' models including concordance, pseudo-R-squared, and ICC.
#' 
#' @param model Fitted coxme object from coxme package.
#' @param raw_data Data.table with raw model information.
#' @param metrics Named list of initialized metrics to populate.
#' @return Updated metrics list with coxme-specific values.
#' @keywords internal
extract_coxme_metrics <- function(model, raw_data, metrics) {
    
    ## Sample size and events
    if (!is.null(model$n)) {
        metrics$events <- model$n[1]
        metrics$n <- model$n[2]
    }
    
    ## Number of groups
    if (!is.null(model$n)) {
        ## coxme stores n as c(events, n, ngroups)
        if (length(model$n) >= 3) {
            metrics$n_groups <- model$n[3]
        }
    }
    
    ## Fixed effects count
    fixed_coefs <- coxme::fixef(model)
    metrics$predictors <- length(fixed_coefs)
    
    ## Concordance/C-statistic
    if (requireNamespace("survival", quietly = TRUE)) {
        c_stat <- tryCatch({
            if (!is.null(model$y) && !is.null(model$linear.predictor)) {
                conc <- survival::concordance(model$y ~ model$linear.predictor)
                conc$concordance
            } else {
                NA_real_
            }
        }, error = function(e) NA_real_)
        
        if (!is.na(c_stat)) {
            metrics$concordance <- c_stat
            metrics$c_statistic <- c_stat
        }
    }
    
    ## Log-likelihood based pseudo-R^2
    loglik <- model$loglik
    if (!is.null(loglik) && length(loglik) >= 2) {
        ## McFadden-style pseudo-R^2
        metrics$pseudo_r2 <- 1 - (loglik[length(loglik)] / loglik[1])
    }
    
    ## Global LR test
    if (!is.null(loglik) && length(loglik) >= 2) {
        lr_stat <- 2 * (loglik[length(loglik)] - loglik[1])
        df <- length(fixed_coefs)
        if (df > 0) {
            metrics$global_p <- pchisq(lr_stat, df, lower.tail = FALSE)
        }
    }
    
    ## ICC from random effects variance
    var_comps <- coxme::VarCorr(model)
    if (!is.null(var_comps)) {
        ## Extract variance components
        total_var <- 0
        for (vc in var_comps) {
            if (is.matrix(vc)) {
                total_var <- total_var + sum(diag(vc))
            } else {
                total_var <- total_var + sum(vc)
            }
        }
        ## For Cox models, there's no residual variance, so ICC is relative
        ## to assumed baseline variance (often pi^2/6 for proportional hazards)
        baseline_var <- pi^2 / 6
        metrics$icc <- total_var / (total_var + baseline_var)
    }
    
    return(metrics)
}


#' Build comprehensive comparison table
#' 
#' Selects and renames columns based on model type following academic consensus
#' for reporting model fit statistics.
#' 
#' @param comparison Data.table with raw comparison metrics.
#' @param model_type Character string indicating model type.
#' @return Formatted data.table with appropriate columns for the model type.
#' @keywords internal
build_comparison_table <- function(comparison, model_type) {
    
    ## Select columns based on model type and academic consensus
    if (model_type == "glm") {
        ## Focus on discrimination, calibration, and information criteria
        key_cols <- c("Model", "N", "Events", "Predictors", "Converged",
                      "AIC", "BIC", "C-statistic", "Brier Score",
                      "McFadden R\u00b2", "Nagelkerke R\u00b2", 
                      "Hoslem p", "Global p")
        
        ## Rename for display
        old_names <- c("c_statistic", "brier_score", "mcfadden_r2", 
                       "nagelkerke_r2", "hoslem_p", "global_p")
        new_names <- c("C-statistic", "Brier Score", "McFadden R\u00b2", 
                       "Nagelkerke R\u00b2", "Hoslem p", "Global p")
        data.table::setnames(comparison, old_names, new_names, skip_absent = TRUE)
        
    } else if (model_type == "coxph") {
        key_cols <- c("Model", "N", "Events", "Predictors", "Converged",
                      "AIC", "BIC", "C-index", "R\u00b2", "R\u00b2 max",
                      "PH test p", "Global p")
        
        old_names <- c("c_index", "rsq", "rsq_max", "ph_global_p", "lr_test_p")
        new_names <- c("C-index", "R\u00b2", "R\u00b2 max", "PH test p", "Global p")
        data.table::setnames(comparison, old_names, new_names, skip_absent = TRUE)
        
    } else if (model_type == "lm") {
        key_cols <- c("Model", "N", "Predictors", "Converged",
                      "AIC", "BIC", "R\u00b2", "Adj R\u00b2", "RMSE",
                      "F-stat", "Global p")
        
        old_names <- c("r_squared", "adj_r_squared", "rmse", "f_statistic", "global_p")
        new_names <- c("R\u00b2", "Adj R\u00b2", "RMSE", "F-stat", "Global p")
        data.table::setnames(comparison, old_names, new_names, skip_absent = TRUE)
        
    } else if (model_type %in% c("lmer", "lmerMod")) {
        key_cols <- c("Model", "N", "Groups", "Predictors", "Converged",
                      "AIC", "BIC", "Marginal R\u00b2", "Conditional R\u00b2", 
                      "ICC", "RMSE", "Global p")
        
        old_names <- c("n_groups", "marginal_r2", "conditional_r2", "icc", "rmse", "global_p")
        new_names <- c("Groups", "Marginal R\u00b2", "Conditional R\u00b2", "ICC", "RMSE", "Global p")
        data.table::setnames(comparison, old_names, new_names, skip_absent = TRUE)
        
    } else if (model_type %in% c("glmer", "glmerMod")) {
        key_cols <- c("Model", "N", "Events", "Groups", "Predictors", "Converged",
                      "AIC", "BIC", "Concordance", "Marginal R\u00b2", "Conditional R\u00b2",
                      "ICC", "Brier Score", "Global p")
        
        old_names <- c("n_groups", "marginal_r2", "conditional_r2", "icc", "brier_score", "global_p")
        new_names <- c("Groups", "Marginal R\u00b2", "Conditional R\u00b2", "ICC", "Brier Score", "Global p")
        data.table::setnames(comparison, old_names, new_names, skip_absent = TRUE)
        
    } else if (model_type == "coxme") {
        key_cols <- c("Model", "N", "Events", "Groups", "Predictors", "Converged",
                      "AIC", "BIC", "Concordance", "Pseudo-R\u00b2", "ICC", "Global p")
        
        old_names <- c("n_groups", "icc", "global_p")
        new_names <- c("Groups", "ICC", "Global p")
        data.table::setnames(comparison, old_names, new_names, skip_absent = TRUE)
    }
    
    ## Keep only relevant columns that exist
    key_cols <- intersect(key_cols, names(comparison))
    comparison <- comparison[, ..key_cols]
    
    ## Format numeric columns appropriately
    format_model_comparison(comparison)
    
    return(comparison)
}


#' Format model comparison table
#' 
#' Rounds numeric columns to appropriate precision based on metric type.
#' 
#' @param comparison Data.table with comparison metrics.
#' @return Data.table with properly rounded numeric columns.
#' @keywords internal
format_model_comparison <- function(comparison) {
    
    ## Round numeric columns to appropriate precision
    numeric_cols <- names(comparison)[sapply(comparison, is.numeric)]
    
    for (col in numeric_cols) {
        if (col %chin% c("AIC", "BIC")) {
            comparison[, (col) := round(get(col), 1)]
        } else if (col %chin% c("C-statistic", "C-index", "Concordance", 
                                "McFadden R\u00b2", "Nagelkerke R\u00b2", "R\u00b2", "Adj R\u00b2",
                                "Marginal R\u00b2", "Conditional R\u00b2", "ICC",
                                "Brier Score", "RMSE", "Pseudo-R\u00b2")) {
            comparison[, (col) := round(get(col), 3)]
        }
    }
    
    return(comparison)
}


#' Calculate Composite Mean Scores (CMS) for model comparison
#' 
#' Computes composite Score based on weighted combination of model
#' quality metrics. Weights vary by model type to reflect academic consensus
#' on important metrics for each model class.
#' 
#' @param comparison Data.table with model comparison metrics.
#' @param model_type Character string indicating model type.
#' @param scoring_weights Optional named list of custom weights. If NULL,
#'   uses default weights for the model type.
#' @return Data.table with CMS column added, sorted by score.
#' @keywords internal
calculate_model_scores <- function(comparison, model_type, scoring_weights = NULL) {
    
    ## Define default weights based on model type
    default_weights <- list(
        glm = list(convergence = 0.15, aic = 0.25, concordance = 0.40, 
                   pseudo_r2 = 0.15, brier = 0.05),
        coxph = list(convergence = 0.15, aic = 0.30, concordance = 0.40, 
                     global_p = 0.15),
        lm = list(convergence = 0.15, aic = 0.25, pseudo_r2 = 0.45, 
                  rmse = 0.15),
        ## Mixed-effects models
        lmer = list(convergence = 0.20, aic = 0.25, marginal_r2 = 0.25,
                    conditional_r2 = 0.15, icc = 0.15),
        glmer = list(convergence = 0.15, aic = 0.25, concordance = 0.30,
                     marginal_r2 = 0.15, icc = 0.15),
        coxme = list(convergence = 0.15, aic = 0.30, concordance = 0.35,
                     pseudo_r2 = 0.10, icc = 0.10)
    )
    
    ## Use provided weights or defaults
    if (is.null(scoring_weights)) {
        weights <- default_weights[[model_type]]
        if (is.null(weights)) {
            weights <- list(convergence = 0.20, aic = 0.40, bic = 0.40)
        }
    } else {
        weights <- scoring_weights
    }
    
    ## Validate weights sum to 1
    weight_sum <- sum(unlist(weights))
    if (abs(weight_sum - 1) > 0.01) {
        warning("Scoring weights do not sum to 1, normalizing...")
        weights <- lapply(weights, function(x) x / weight_sum)
    }
    
    ## Initialize scores with pre-allocation
    n_models <- nrow(comparison)
    scores <- list(
        conv_score = numeric(n_models),
        aic_score = numeric(n_models),
        concordance_score = numeric(n_models),
        pseudo_r2_score = numeric(n_models),
        brier_score = numeric(n_models),
        global_score = numeric(n_models),
        rmse_score = numeric(n_models),
        bic_score = numeric(n_models),
        marginal_r2_score = numeric(n_models),
        conditional_r2_score = numeric(n_models),
        icc_score = numeric(n_models),
        total = numeric(n_models)
    )
    
    ## Convergence score (universal)
    scores$conv_score <- data.table::fcase(
                                         comparison$Converged == "Yes", 100,
                                         comparison$Converged == "Suspect", 70,
                                         comparison$Converged == "No", 30,
                                         default = 0
                                     )
    
    ## AIC score (universal, lower is better)
    if (!all(is.na(comparison$AIC))) {
        aic_values <- comparison$AIC[!is.na(comparison$AIC)]
        if (length(aic_values) > 1) {
            aic_best <- min(aic_values)
            aic_worst <- max(aic_values)
            if (aic_worst - aic_best > 0) {
                scores$aic_score <- data.table::fifelse(
                                                    is.na(comparison$AIC), 0,
                                                    100 * (1 - (comparison$AIC - aic_best)/(aic_worst - aic_best))
                                                )
            } else {
                scores$aic_score <- rep(100, n_models)
            }
        } else {
            scores$aic_score <- data.table::fifelse(is.na(comparison$AIC), 0, 100)
        }
        scores$aic_score[is.na(scores$aic_score)] <- 0
    } else {
        scores$aic_score <- rep(50, n_models)
    }
    
    ## Model-specific scores
    if (model_type == "glm") {
        ## Concordance/C-statistic
        if ("Concordance" %chin% names(comparison) && !all(is.na(comparison$Concordance))) {
            scores$concordance_score <- data.table::fcase(
                                                        is.na(comparison$Concordance), 0,
                                                        comparison$Concordance <= 0.5, 0,
                                                        comparison$Concordance <= 0.6, 40 * (comparison$Concordance - 0.5)/0.1,
                                                        comparison$Concordance <= 0.7, 40 + 20 * (comparison$Concordance - 0.6)/0.1,
                                                        comparison$Concordance <= 0.8, 60 + 20 * (comparison$Concordance - 0.7)/0.1,
                                                        comparison$Concordance <= 0.9, 80 + 10 * (comparison$Concordance - 0.8)/0.1,
                                                        default = 90 + 10 * (comparison$Concordance - 0.9)/0.1
                                                    )
            scores$concordance_score <- pmin(100, scores$concordance_score)
        } else {
            scores$concordance_score <- rep(50, n_models)
        }
        
        ## Pseudo-R^2
        if ("Pseudo-R^2" %chin% names(comparison) && !all(is.na(comparison[["Pseudo-R^2"]]))) {
            ## McFadden's R^2 rarely exceeds 0.4 for good models
            scores$pseudo_r2_score <- data.table::fifelse(
                                                      is.na(comparison[["Pseudo-R^2"]]), 0,
                                                      pmin(100, comparison[["Pseudo-R^2"]] * 250)  # 0.4 -> 100 points
                                                  )
        } else {
            scores$pseudo_r2_score <- rep(50, n_models)
        }

        ## Brier score - only if column exists and weight is non-zero
        if ("Brier Score" %chin% names(comparison) && "brier" %chin% names(weights) && weights$brier > 0) {
            ## Lower is better (0 = perfect, 0.25 = no skill)
            scores$brier_score <- data.table::fifelse(
                                                  is.na(comparison$`Brier Score`), 50,
                                                  100 * (1 - comparison$`Brier Score`/0.25)
                                              )
        } else {
            scores$brier_score <- rep(50, n_models)
            weights$brier <- 0  # Set weight to 0 if not using
        }
        
        ## Normalize weights if Brier is excluded
        if (weights$brier == 0) {
            remaining_weights <- weights[names(weights) != "brier"]
            weight_sum <- sum(unlist(remaining_weights))
            remaining_weights <- lapply(remaining_weights, function(x) x / weight_sum)
            weights[names(remaining_weights)] <- remaining_weights
        }
        
        ## Calculate weighted total
        scores$total <- scores$conv_score * weights$convergence +
            scores$aic_score * weights$aic +
            scores$concordance_score * weights$concordance +
            scores$pseudo_r2_score * weights$pseudo_r2
        
    } else if (model_type == "coxph") {
        ## Similar adjustments for Cox models
        if ("Concordance" %chin% names(comparison) && !all(is.na(comparison$Concordance))) {
            scores$concordance_score <- data.table::fcase(
                                                        is.na(comparison$Concordance), 0,
                                                        comparison$Concordance <= 0.5, 0,
                                                        default = pmin(100, 200 * (comparison$Concordance - 0.5))
                                                    )
        } else {
            scores$concordance_score <- rep(50, n_models)
        }
        
        ## Global p-value score
        if ("Global p" %chin% names(comparison)) {
            global_numeric <- suppressWarnings(as.numeric(gsub("< ", "", comparison$`Global p`)))
            scores$global_score <- data.table::fcase(
                                                   is.na(global_numeric), 50,
                                                   global_numeric > 0.05, 50,
                                                   default = 50 + 50 * (0.05 - global_numeric)/0.05
                                               )
        } else {
            scores$global_score <- rep(50, n_models)
        }
        
        scores$total <- scores$conv_score * weights$convergence +
            scores$aic_score * weights$aic +
            scores$concordance_score * weights$concordance +
            scores$global_score * weights$global_p
        
    } else if (model_type == "lm") {
        ## R^2 scoring
        if ("Pseudo-R^2" %chin% names(comparison) && !all(is.na(comparison[["Pseudo-R^2"]]))) {
            scores$pseudo_r2_score <- comparison[["Pseudo-R^2"]] * 100
            scores$pseudo_r2_score[is.na(scores$pseudo_r2_score)] <- 0
        } else {
            scores$pseudo_r2_score <- rep(50, n_models)
        }
        
        ## RMSE scoring would go here if column exists
        scores$rmse_score <- rep(50, n_models)  # Placeholder
        
        scores$total <- scores$conv_score * weights$convergence +
            scores$aic_score * weights$aic +
            scores$pseudo_r2_score * weights$pseudo_r2 +
            scores$rmse_score * weights$rmse
        
    } else if (model_type == "lmer") {
        ## Linear mixed-effects models
        scores <- calculate_lmer_scores(comparison, weights, scores, n_models)
        
    } else if (model_type == "glmer") {
        ## Generalized linear mixed-effects models
        scores <- calculate_glmer_scores(comparison, weights, scores, n_models)
        
    } else if (model_type == "coxme") {
        ## Mixed-effects Cox models
        scores <- calculate_coxme_scores(comparison, weights, scores, n_models)
        
    } else {
        ## Generic fallback
        scores$bic_score <- rep(50, n_models)
        if (!all(is.na(comparison$BIC))) {
            bic_values <- comparison$BIC[!is.na(comparison$BIC)]
            if (length(bic_values) > 1) {
                bic_best <- min(bic_values)
                bic_worst <- max(bic_values)
                if (bic_worst - bic_best > 0) {
                    scores$bic_score <- 100 * (1 - (comparison$BIC - bic_best)/(bic_worst - bic_best))
                } else {
                    scores$bic_score <- rep(100, n_models)
                }
            } else {
                scores$bic_score <- rep(100, n_models)
            }
            scores$bic_score[is.na(scores$bic_score)] <- 0
        }
        
        scores$total <- scores$conv_score * weights$convergence +
            scores$aic_score * weights$aic +
            scores$bic_score * weights$bic
    }
    
    ## Add final score to comparison table using := for efficiency
    comparison[, CMS := round(scores$total, 1)]
    
    ## Sort by score (highest first)
    data.table::setorder(comparison, -CMS)
    
    ## Store detailed scores as attribute
    data.table::setattr(comparison, "detailed_scores", scores)
    data.table::setattr(comparison, "weights", weights)
    
    return(comparison)
}


#' Calculate scores for linear mixed-effects models
#' 
#' Computes component scores for lmer models based on marginal R-squared,
#' conditional R-squared, and ICC metrics.
#' 
#' @param comparison Data.table with model comparison metrics.
#' @param weights Named list of scoring weights.
#' @param scores List of initialized score vectors.
#' @param n_models Integer number of models being compared.
#' @return Updated scores list with calculated values and total.
#' @keywords internal
calculate_lmer_scores <- function(comparison, weights, scores, n_models) {
    
    ## Marginal R^2 score (variance explained by fixed effects)
    ## Higher is better, typically ranges 0-0.8 for good models
    marg_r2_col <- if ("Marginal R^2" %chin% names(comparison)) "Marginal R^2" else "Pseudo-R^2"
    if (marg_r2_col %chin% names(comparison) && !all(is.na(comparison[[marg_r2_col]]))) {
        scores$marginal_r2_score <- data.table::fifelse(
                                                    is.na(comparison[[marg_r2_col]]), 0,
                                                    pmin(100, comparison[[marg_r2_col]] * 125)  # 0.8 -> 100 points
                                                )
    } else {
        scores$marginal_r2_score <- rep(50, n_models)
    }
    
    ## Conditional R^2 score (variance explained by full model)
    if ("Conditional R^2" %chin% names(comparison) && !all(is.na(comparison[["Conditional R^2"]]))) {
        scores$conditional_r2_score <- data.table::fifelse(
                                                       is.na(comparison[["Conditional R^2"]]), 0,
                                                       pmin(100, comparison[["Conditional R^2"]] * 110)  # 0.9 -> ~100 points
                                                   )
    } else {
        scores$conditional_r2_score <- rep(50, n_models)
    }
    
    ## ICC score - moderate ICC (0.05-0.5) is typically good
    ## Too low suggests random effects unnecessary, too high suggests problems
    if ("ICC" %chin% names(comparison) && !all(is.na(comparison$ICC))) {
        scores$icc_score <- data.table::fcase(
                                            is.na(comparison$ICC), 50,
                                            comparison$ICC < 0.01, 30,    # Very low - random effects may not be needed
                                            comparison$ICC < 0.05, 60,    # Low but present
                                            comparison$ICC < 0.25, 100,   # Ideal range
                                            comparison$ICC < 0.50, 90,    # Still good
                                            comparison$ICC < 0.75, 70,    # High - check for model issues
                                            default = 50                  # Very high - potential problems
                                        )
    } else {
        scores$icc_score <- rep(50, n_models)
    }
    
    ## Calculate weighted total
    scores$total <- scores$conv_score * weights$convergence +
        scores$aic_score * weights$aic +
        scores$marginal_r2_score * weights$marginal_r2 +
        scores$conditional_r2_score * weights$conditional_r2 +
        scores$icc_score * weights$icc
    
    return(scores)
}


#' Calculate scores for generalized linear mixed-effects models
#' 
#' Computes component scores for glmer models based on concordance,
#' marginal R-squared, and ICC metrics.
#' 
#' @param comparison Data.table with model comparison metrics.
#' @param weights Named list of scoring weights.
#' @param scores List of initialized score vectors.
#' @param n_models Integer number of models being compared.
#' @return Updated scores list with calculated values and total.
#' @keywords internal
calculate_glmer_scores <- function(comparison, weights, scores, n_models) {
    
    ## Concordance score (same as GLM)
    if ("Concordance" %chin% names(comparison) && !all(is.na(comparison$Concordance))) {
        scores$concordance_score <- data.table::fcase(
                                                    is.na(comparison$Concordance), 0,
                                                    comparison$Concordance <= 0.5, 0,
                                                    comparison$Concordance <= 0.6, 40 * (comparison$Concordance - 0.5)/0.1,
                                                    comparison$Concordance <= 0.7, 40 + 20 * (comparison$Concordance - 0.6)/0.1,
                                                    comparison$Concordance <= 0.8, 60 + 20 * (comparison$Concordance - 0.7)/0.1,
                                                    comparison$Concordance <= 0.9, 80 + 10 * (comparison$Concordance - 0.8)/0.1,
                                                    default = 90 + 10 * (comparison$Concordance - 0.9)/0.1
                                                )
        scores$concordance_score <- pmin(100, scores$concordance_score)
    } else {
        scores$concordance_score <- rep(50, n_models)
    }
    
    ## Marginal R^2 score
    marg_r2_col <- if ("Marginal R^2" %chin% names(comparison)) "Marginal R^2" else "Pseudo-R^2"
    if (marg_r2_col %chin% names(comparison) && !all(is.na(comparison[[marg_r2_col]]))) {
        ## For GLMMs, McFadden-style R^2 is lower, scale accordingly
        scores$marginal_r2_score <- data.table::fifelse(
                                                    is.na(comparison[[marg_r2_col]]), 0,
                                                    pmin(100, comparison[[marg_r2_col]] * 200)  # 0.5 -> 100 points
                                                )
    } else {
        scores$marginal_r2_score <- rep(50, n_models)
    }
    
    ## ICC score (same as lmer)
    if ("ICC" %chin% names(comparison) && !all(is.na(comparison$ICC))) {
        scores$icc_score <- data.table::fcase(
                                            is.na(comparison$ICC), 50,
                                            comparison$ICC < 0.01, 30,
                                            comparison$ICC < 0.05, 60,
                                            comparison$ICC < 0.25, 100,
                                            comparison$ICC < 0.50, 90,
                                            comparison$ICC < 0.75, 70,
                                            default = 50
                                        )
    } else {
        scores$icc_score <- rep(50, n_models)
    }
    
    ## Calculate weighted total
    scores$total <- scores$conv_score * weights$convergence +
        scores$aic_score * weights$aic +
        scores$concordance_score * weights$concordance +
        scores$marginal_r2_score * weights$marginal_r2 +
        scores$icc_score * weights$icc
    
    return(scores)
}


#' Calculate scores for mixed-effects Cox models
#' 
#' Computes component scores for coxme models based on concordance,
#' pseudo-R-squared, and ICC metrics.
#' 
#' @param comparison Data.table with model comparison metrics.
#' @param weights Named list of scoring weights.
#' @param scores List of initialized score vectors.
#' @param n_models Integer number of models being compared.
#' @return Updated scores list with calculated values and total.
#' @keywords internal
calculate_coxme_scores <- function(comparison, weights, scores, n_models) {
    
    ## Concordance score
    if ("Concordance" %chin% names(comparison) && !all(is.na(comparison$Concordance))) {
        scores$concordance_score <- data.table::fcase(
                                                    is.na(comparison$Concordance), 0,
                                                    comparison$Concordance <= 0.5, 0,
                                                    default = pmin(100, 200 * (comparison$Concordance - 0.5))
                                                )
    } else {
        scores$concordance_score <- rep(50, n_models)
    }
    
    ## Pseudo-R^2 score
    if ("Pseudo-R^2" %chin% names(comparison) && !all(is.na(comparison[["Pseudo-R^2"]]))) {
        scores$pseudo_r2_score <- data.table::fifelse(
                                                  is.na(comparison[["Pseudo-R^2"]]), 0,
                                                  pmin(100, comparison[["Pseudo-R^2"]] * 250)
                                              )
    } else {
        scores$pseudo_r2_score <- rep(50, n_models)
    }
    
    ## ICC score
    if ("ICC" %chin% names(comparison) && !all(is.na(comparison$ICC))) {
        scores$icc_score <- data.table::fcase(
                                            is.na(comparison$ICC), 50,
                                            comparison$ICC < 0.01, 30,
                                            comparison$ICC < 0.05, 60,
                                            comparison$ICC < 0.25, 100,
                                            comparison$ICC < 0.50, 90,
                                            comparison$ICC < 0.75, 70,
                                            default = 50
                                        )
    } else {
        scores$icc_score <- rep(50, n_models)
    }
    
    ## Calculate weighted total
    scores$total <- scores$conv_score * weights$convergence +
        scores$aic_score * weights$aic +
        scores$concordance_score * weights$concordance +
        scores$pseudo_r2_score * weights$pseudo_r2 +
        scores$icc_score * weights$icc
    
    return(scores)
}


#' Combine coefficient tables from multiple models
#' 
#' Merges coefficient tables from multiple fitted models into a single
#' data.table with a Model identifier column.
#' 
#' @param coef_list List of data.tables containing coefficient information.
#' @param model_names Character vector of model names corresponding to coef_list.
#' @return Combined data.table with Model column, or NULL if empty.
#' @keywords internal
combine_coefficient_tables <- function(coef_list, model_names) {
    if (length(coef_list) == 0) return(NULL)
    
    ## Add model identifier to each table
    for (i in seq_along(coef_list)) {
        if (!is.null(coef_list[[i]])) {
            coef_list[[i]]$Model <- model_names[i]
        }
    }
    
    ## Combine all tables
    combined <- data.table::rbindlist(coef_list, fill = TRUE)
    
    ## Reorder columns to put Model first
    if ("Model" %in% names(combined)) {
        cols <- names(combined)
        new_order <- c("Model", setdiff(cols, "Model"))
        combined <- combined[, ..new_order]
    }
    
    return(combined)
}


#' Auto-detect model type based on outcome and random effects
#' 
#' Determines the appropriate model type based on outcome variable
#' characteristics and presence of random effects.
#' 
#' @param data Data.frame or data.table containing the outcome variable.
#' @param outcome Character string specifying the outcome variable or Surv() expression.
#' @param has_random_effects Logical indicating if random effects are specified.
#' @param family Character string for GLM family (default "binomial").
#' @return Character string indicating detected model type.
#' @keywords internal
detect_model_type_auto <- function(data, outcome, has_random_effects, family = "binomial") {
    
    is_survival <- grepl("^Surv\\(", outcome)
    
    if (is_survival) {
        if (has_random_effects) {
            message("Auto-detected survival outcome with random effects, using coxme")
            return("coxme")
        } else {
            message("Auto-detected survival outcome, using Cox regression")
            return("coxph")
        }
    }
    
    ## Check outcome type
    is_binary <- is.factor(data[[outcome]]) || 
        (is.numeric(data[[outcome]]) && length(unique(data[[outcome]])) == 2)
    is_continuous <- is.numeric(data[[outcome]]) && length(unique(data[[outcome]])) > 2
    
    if (is_binary) {
        if (has_random_effects) {
            message("Auto-detected binary outcome with random effects, using glmer")
            return("glmer")
        } else {
            message("Auto-detected binary outcome, using logistic regression")
            return("glm")
        }
    } else if (is_continuous) {
        if (has_random_effects) {
            message("Auto-detected continuous outcome with random effects, using lmer")
            return("lmer")
        } else {
            message("Auto-detected continuous outcome, using linear regression")
            return("lm")
        }
    } else {
        stop("Cannot auto-detect model type for outcome: ", outcome)
    }
}


#' Normalize model type names
#' 
#' Converts model class names to standardized model type strings.
#' 
#' @param model_type Character string of model type or class name.
#' @return Normalized character string (e.g., "lmerMod" becomes "lmer").
#' @keywords internal
normalize_model_type <- function(model_type) {
    switch(model_type,
           "lmerMod" = "lmer",
           "glmerMod" = "glmer",
           model_type)
}


#' Check required packages for model type
#' 
#' Verifies that necessary packages are installed for the specified model type.
#' Stops with informative error if required packages are missing.
#' 
#' @param model_type Character string indicating model type.
#' @return NULL (invisibly). Stops execution if packages missing.
#' @keywords internal
check_required_packages <- function(model_type) {
    
    if (model_type %in% c("lmer", "glmer")) {
        if (!requireNamespace("lme4", quietly = TRUE)) {
            stop("Package 'lme4' is required for mixed-effects models. ",
                 "Install with: install.packages('lme4')")
        }
    }
    
    if (model_type == "coxme") {
        if (!requireNamespace("coxme", quietly = TRUE)) {
            stop("Package 'coxme' is required for mixed-effects Cox models. ",
                 "Install with: install.packages('coxme')")
        }
    }
    
    if (model_type %in% c("coxph", "clogit")) {
        if (!requireNamespace("survival", quietly = TRUE)) {
            stop("Package 'survival' is required for Cox models. ",
                 "Install with: install.packages('survival')")
        }
    }
}


#' Build row for failed model
#' 
#' Creates a comparison table row with appropriate NA values for a model
#' that failed to fit.
#' 
#' @param model_name Character string name of the model.
#' @param n Integer sample size.
#' @param n_predictors Integer number of predictors attempted.
#' @param model_type Character string indicating model type.
#' @return Data.table with single row of NA metrics and "Failed" convergence.
#' @keywords internal
build_failed_model_row <- function(model_name, n, n_predictors, model_type) {
    
    row <- data.table::data.table(
                           Model = model_name,
                           N = n,
                           Events = NA_integer_,
                           Predictors = n_predictors,
                           Converged = "Failed",
                           AIC = NA_real_,
                           BIC = NA_real_,
                           `Pseudo-R^2` = NA_real_,
                           Concordance = NA_real_,
                           `Global p` = NA_real_
                       )
    
    ## Add model-type specific columns
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
        if (model_type == "lmer") {
            row$RMSE <- NA_real_
        }
    }
    
    return(row)
}


#' Build comparison row for successfully fitted model
#' 
#' Creates a comparison table row with extracted metrics for a successfully
#' fitted model.
#' 
#' @param model_name Character string name of the model.
#' @param n_predictors Integer number of predictors in the model.
#' @param converged Character string convergence status.
#' @param metrics Named list of extracted model metrics.
#' @param model_type Character string indicating model type.
#' @return Data.table with single row of formatted metrics.
#' @keywords internal
build_model_row <- function(model_name, n_predictors, converged, metrics, model_type) {
    
    ## Base columns for all models
    row <- data.table::data.table(
                           Model = model_name,
                           N = metrics$n,
                           Events = metrics$events,
                           Predictors = n_predictors,
                           Converged = converged,
                           AIC = safe_round(metrics$aic, 1),
                           BIC = safe_round(metrics$bic, 1),
                           `Pseudo-R^2` = safe_round(metrics$pseudo_r2, 3),
                           Concordance = safe_round(metrics$concordance, 3),
                           `Global p` = if (!is.null(metrics$global_p) && !is.na(metrics$global_p)) {
                                            format_pvalues_fit(metrics$global_p, 3)
                                        } else {
                                            NA_character_
                                        }
                       )
    
    ## Add GLM-specific columns
    if (model_type == "glm") {
        row$`Brier Score` <- safe_round(metrics$brier_score, 3)
    }
    
    ## Add mixed-effects specific columns
    if (model_type %in% c("lmer", "glmer", "coxme")) {
        row$Groups <- metrics$n_groups
        row$`Marginal R^2` <- safe_round(metrics$marginal_r2, 3)
        row$`Conditional R^2` <- safe_round(metrics$conditional_r2, 3)
        row$ICC <- safe_round(metrics$icc, 3)
        
        if (model_type == "glmer") {
            row$`Brier Score` <- safe_round(metrics$brier_score, 3)
        }
        if (model_type == "lmer") {
            row$RMSE <- safe_round(metrics$rmse, 3)
        }
    }
    
    return(row)
}


#' Safe rounding that handles NULL and NA
#' 
#' Rounds numeric values while gracefully handling NULL, empty, and NA inputs.
#' 
#' @param x Numeric value to round.
#' @param digits Integer number of decimal places.
#' @return Rounded numeric value, or NA_real_ if input is NULL/NA/empty.
#' @keywords internal
safe_round <- function(x, digits) {
    if (is.null(x) || length(x) == 0) return(NA_real_)
    if (is.na(x)) return(NA_real_)
    round(x, digits)
}


#' Order comparison columns based on model type
#' 
#' Reorders columns in the comparison table to follow a logical sequence
#' appropriate for the model type.
#' 
#' @param comparison Data.table with comparison metrics.
#' @param model_type Character string indicating model type.
#' @return Data.table with reordered columns.
#' @keywords internal
order_comparison_columns <- function(comparison, model_type) {
    
    col_order <- switch(model_type,
                        "glm" = c("Model", "N", "Events", "Predictors", "Converged",
                                  "AIC", "BIC", "Pseudo-R\u00b2", "Concordance", "Brier Score",
                                  "Global p", "CMS"),
                        "coxph" = c("Model", "N", "Events", "Predictors", "Converged",
                                    "AIC", "BIC", "Pseudo-R\u00b2", "Concordance",
                                    "Global p", "CMS"),
                        "lm" = c("Model", "N", "Predictors", "Converged",
                                 "AIC", "BIC", "Pseudo-R\u00b2", "Global p", "CMS"),
                        "lmer" = c("Model", "N", "Groups", "Predictors", "Converged",
                                   "AIC", "BIC", "Marginal R\u00b2", "Conditional R\u00b2", "ICC",
                                   "RMSE", "Global p", "CMS"),
                        "glmer" = c("Model", "N", "Events", "Groups", "Predictors", "Converged",
                                    "AIC", "BIC", "Concordance", "Marginal R\u00b2", "Conditional R\u00b2",
                                    "ICC", "Brier Score", "Global p", "CMS"),
                        "coxme" = c("Model", "N", "Events", "Groups", "Predictors", "Converged",
                                    "AIC", "BIC", "Concordance", "Pseudo-R\u00b2", "ICC",
                                    "Global p", "CMS"),
                        ## Default fallback
                        c("Model", "N", "Events", "Predictors", "Converged",
                          "AIC", "BIC", "Pseudo-R\u00b2", "Concordance",
                          "Global p", "CMS")
                        )
    
    existing_cols <- intersect(col_order, names(comparison))
    extra_cols <- setdiff(names(comparison), col_order)
    
    data.table::setcolorder(comparison, c(existing_cols, extra_cols))
    
    return(comparison)
}
