#' Convert Model to Data Table
#'
#' Extracts coefficients, confidence intervals, and comprehensive model statistics 
#' from fitted regression models and converts them to a standardized data.table 
#' format suitable for further analysis or publication. This is a core utility 
#' function used internally by other \pkg{summata} regression functions.
#'
#' @param data Data frame or data.table containing the dataset used to fit the
#'   model. Required for computing group-level sample sizes and event counts.
#'
#' @param model Fitted model object. Supported classes include:
#'   \itemize{
#'     \item \code{glm} - Generalized linear models (logistic, Poisson, etc.)
#'     \item \code{lm} - Linear models
#'     \item \code{coxph} - Cox proportional hazards models
#'     \item \code{clogit} - Conditional logistic regression
#'     \item \code{coxme} - Mixed effects Cox models
#'     \item \code{lmerMod} - Linear mixed effects models
#'     \item \code{glmerMod} - Generalized linear mixed effects models
#'   }
#'   
#' @param conf_level Numeric confidence level for confidence intervals. Must be
#'   between 0 and 1. Default is 0.95 (95\% CI).
#'   
#' @param keep_qc_stats Logical. If \code{TRUE}, includes model quality statistics 
#'   such as AIC, BIC, R², concordance, and model fit tests. These appear 
#'   as additional columns in the output. Default is \code{TRUE}.
#'   
#' @param include_intercept Logical. If \code{TRUE}, includes the model intercept 
#'   in output. If \code{FALSE}, removes the intercept row from results. Useful 
#'   for creating cleaner presentation tables. Default is \code{TRUE}.
#'   
#' @param terms_to_exclude Character vector of term names to exclude from output.
#'   Useful for removing specific unwanted parameters (e.g., nuisance variables,
#'   spline terms). Default is \code{NULL}. Note: If \code{include_intercept = FALSE}, 
#'   "(Intercept)" is automatically added to this list.
#'   
#' @param reference_rows Logical. If \code{TRUE}, adds rows for reference 
#'   categories of factor variables with appropriate labels and baseline values 
#'   (OR/HR = 1, Coefficient = 0). This makes tables more complete and easier to 
#'   interpret. Default is \code{TRUE}.
#'   
#' @param reference_label Character string used to label reference category rows
#'   in the output. Appears in the \code{reference} column. Default is \code{"reference"}.
#'
#' @param skip_counts Logical. If \code{TRUE}, skips computation of group-level
#'   sample sizes and event counts (faster but less informative). Default is 
#'   \code{FALSE}.
#'
#' @return A \code{data.table} containing extracted model information with the 
#'   following standard columns:
#'   \describe{
#'     \item{model_scope}{Character. Either "Univariable" (unadjusted model with 
#'       single predictor) or "Multivariable" (adjusted model with multiple predictors)}
#'     \item{model_type}{Character. Type of regression (e.g., "Logistic", "Linear", 
#'       "Cox PH", "Poisson", etc.)}
#'     \item{variable}{Character. Variable name (for factor variables, the base 
#'       variable name without the level)}
#'     \item{group}{Character. Group/level name for factor variables; empty string 
#'       for continuous variables}
#'     \item{n}{Integer. Total sample size used in the model}
#'     \item{n_group}{Integer. Sample size for this specific variable level 
#'       (factor variables only)}
#'     \item{events}{Integer. Total number of events in the model (for survival 
#'       and logistic models)}
#'     \item{events_group}{Integer. Number of events for this specific variable 
#'       level (for survival and logistic models with factor variables)}
#'     \item{coefficient}{Numeric. Raw regression coefficient (log odds, log hazard, 
#'       etc.)}
#'     \item{se}{Numeric. Standard error of the coefficient}
#'     \item{OR/HR/RR/Coefficient}{Numeric. Effect estimate - column name depends on 
#'       model type:
#'       \itemize{
#'         \item \code{OR} for logistic regression (odds ratio)
#'         \item \code{HR} for Cox models (hazard ratio)
#'         \item \code{RR} for Poisson regression (rate/risk ratio)
#'         \item \code{Coefficient} for linear models or other GLMs
#'       }}
#'     \item{ci_lower}{Numeric. Lower bound of confidence interval for effect estimate}
#'     \item{ci_upper}{Numeric. Upper bound of confidence interval for effect estimate}
#'     \item{statistic}{Numeric. Test statistic (z-value for GLM/Cox, t-value for LM)}
#'     \item{p_value}{Numeric. P-value for coefficient test}
#'     \item{sig}{Character. Significance markers: "***" (p<0.001), "**" (p < 0.01), 
#'       "*" (p < 0.05), "." (p < 0.10), "" (p ≥ 0.10)}
#'     \item{sig_binary}{Logical. Binary indicator: \code{TRUE} if p < 0.05, 
#'       \code{FALSE} otherwise}
#'     \item{reference}{Character. Contains \code{reference_label} for reference 
#'       category rows when \code{reference_rows = TRUE}, empty string otherwise}
#'   }
#'
#' @details
#' This function is the core extraction utility used by \code{fit()} and other
#' regression functions. It handles the complexities of different model classes
#' and provides a consistent output format suitable for tables and forest plots.
#'
#' \strong{Model Type Detection:}
#' The function automatically detects model type and applies appropriate:
#' \itemize{
#'   \item Effect measure naming (OR, HR, RR, Coefficient)
#'   \item Confidence interval calculation method
#'   \item Event counting for binary/survival outcomes
#' }
#'
#' \strong{Mixed Effects Models:}
#' For lme4 models (glmer, lmer), the function extracts fixed effects only.
#' Random effects structure is not included in the output table.
#'
#' @seealso
#' \code{\link{fit}} for the main regression interface,
#' \code{\link{glmforest}}, \code{\link{coxforest}}, \code{\link{lmforest}} for
#' forest plot visualization
#'
#' @examples
#' # Load example data
#' data(clintrial)
#' 
#' # Example 1: Extract from logistic regression
#' glm_model <- glm(os_status ~ age + sex + treatment, 
#'                  data = clintrial, family = binomial)
#' 
#' glm_result <- m2dt(clintrial, glm_model)
#' glm_result
#' 
#' \donttest{
#' # Example 2: Extract from linear model
#' lm_model <- lm(los_days ~ age + sex + surgery, data = clintrial)
#' 
#' lm_result <- m2dt(clintrial, lm_model)
#' lm_result
#' 
#' # Example 3: Cox proportional hazards model
#' library(survival)
#' cox_model <- coxph(Surv(os_months, os_status) ~ age + sex + stage,
#'                    data = clintrial)
#' 
#' cox_result <- m2dt(clintrial, cox_model)
#' cox_result
#' 
#' # Example 4: Exclude intercept for cleaner tables
#' clean_result <- m2dt(clintrial, glm_model, include_intercept = FALSE)
#' clean_result
#' 
#' # Example 5: Change confidence level
#' result_90ci <- m2dt(clintrial, glm_model, conf_level = 0.90)
#' result_90ci
#'
#' }
#'
#' @export
m2dt <- function(data,
                 model,
                 conf_level = 0.95,
                 keep_qc_stats = TRUE,
                 include_intercept = TRUE,
                 terms_to_exclude = NULL,
                 reference_rows = TRUE,
                 reference_label = "reference",
                 skip_counts = FALSE) {

    ## Validate inputs
    if (missing(data)) {
        stop("data parameter is required. Usage: m2dt(data, model)")
    }
    
    if (!inherits(data, c("data.frame", "data.table"))) {
        stop("data must be a data.frame or data.table")
    }
    
    ## Store data for use throughout function
    model_data <- if (data.table::is.data.table(data)) data else data.table::as.data.table(data)
    
    ## Store as attribute for helper functions
    attr(model, "data") <- model_data

    ## Set model class - comprehensive detection for all model types
    model_classes <- class(model)
    
    ## Priority order for mixed models (check specific classes first)
    if ("lmerMod" %in% model_classes || inherits(model, "lmerMod")) {
        model_class <- "lmerMod"
    } else if ("glmerMod" %in% model_classes || inherits(model, "glmerMod")) {
        model_class <- "glmerMod"
    } else if ("lmerTest" %in% model_classes) {
        ## lmerTest package wraps lmer models
        model_class <- "lmerMod"
    } else if (inherits(model, "merMod")) {
        ## Generic merMod - determine if linear or generalized
        if ("lmer" %in% model_classes) {
            model_class <- "lmerMod"
        } else if ("glmer" %in% model_classes) {
            model_class <- "glmerMod"
        } else {
            ## Fallback for generic merMod
            model_class <- "lmerMod"  # Assume linear if not specified
        }
    } else if ("negbin" %in% model_classes) {
        ## MASS::glm.nb produces negbin class - treat as glm for extraction
        ## but mark it for special handling (rate ratios)
        model_class <- "negbin"
    } else {
        ## Standard detection for other models
        model_class <- class(model)[1]
    }
    
    ## Null operator for handling missing values
    `%||%` <- function(a, b) if (is.null(a)) b else a
    
    ## Add intercept to exclusion list if requested
    if (!include_intercept) {
        terms_to_exclude <- unique(c(terms_to_exclude, "(Intercept)"))
    }
    
    ## Get model type (Univariable vs Multivariable)
    model_scope <- detect_model_type(model)
    
    ## Get readable model type name
    model_type_name <- get_model_type_name(model)
    
    ## Extract results based on model class
    if (model_class %in% c("glm", "lm", "negbin")) {
        
        coef_summary <- summary(model)$coefficients
        conf_int <- stats::confint.default(model, level = conf_level)
        
        dt <- data.table::data.table(
                              model_scope = model_scope,
                              model_type = model_type_name,
                              term = rownames(coef_summary),
                              n = stats::nobs(model),
                              events = NA_real_,
                              coefficient = coef_summary[, "Estimate"],
                              se = coef_summary[, "Std. Error"]
                          )
        
        ## Calculate confidence intervals
        z_score <- stats::qnorm((1 + conf_level) / 2)
        dt[, `:=`(
            coef = coefficient,
            coef_lower = coefficient - z_score * se,
            coef_upper = coefficient + z_score * se
        )]
        
        ## Special handling for logistic regression (including quasibinomial)
        if (model_class == "glm" && !isS4(model)) {
            if (model$family$family %in% c("binomial", "quasibinomial", "poisson", "quasipoisson")) {
                if (!is.null(model$y)) {
                    dt[, events := sum(model$y, na.rm = TRUE)]
                } else if (!is.null(model$model)) {
                    outcome_col <- model$model[[1]]
                    if (is.factor(outcome_col)) {
                        dt[, events := sum(as.numeric(outcome_col) == 2, na.rm = TRUE)]
                    } else {
                        dt[, events := sum(outcome_col, na.rm = TRUE)]
                    }
                }
            }
        }
        
        ## Determine if should exponentiate
        should_exp <- FALSE
        is_logistic <- FALSE
        is_poisson <- FALSE
        is_negbin <- model_class == "negbin"
        
        if (model_class == "glm" && !isS4(model)) {
            family_name <- model$family$family
            link_name <- model$family$link
            
            ## binomial and quasibinomial both produce odds ratios
            is_logistic <- family_name %in% c("binomial", "quasibinomial")
            is_poisson <- family_name == "poisson" || 
                (family_name == "quasipoisson")
            
            should_exp <- is_logistic || is_poisson || link_name == "log"
        }
        
        ## Negative binomial uses log link - exponentiate for rate ratios
        if (is_negbin) {
            should_exp <- TRUE
        }
        
        ## Add exponentiated coefficients
        dt[, `:=`(
            exp_coef = if (should_exp) exp(coefficient) else coefficient,
            exp_lower = if (should_exp) exp(coef_lower) else coef_lower,
            exp_upper = if (should_exp) exp(coef_upper) else coef_upper
        )]
        
        if (is_logistic) {
            dt[, `:=`(
                OR = exp_coef,
                ci_lower = exp_lower,
                ci_upper = exp_upper
            )]
        } else if (is_poisson || is_negbin) {
            ## Poisson and negative binomial both produce rate ratios
            dt[, `:=`(
                RR = exp_coef,
                ci_lower = exp_lower,
                ci_upper = exp_upper
            )]
        } else if (should_exp) {
            ## Other log-link models
            dt[, `:=`(
                Coefficient = exp_coef,
                ci_lower = exp_lower,
                ci_upper = exp_upper
            )]
        } else {
            ## Linear models - use raw coefficients
            dt[, `:=`(
                Coefficient = coef,
                ci_lower = coef_lower,
                ci_upper = coef_upper
            )]
        }
        
        ## Add test statistics
        stat_col <- if ("z value" %in% colnames(coef_summary)) "z value" else "t value"
        dt[, `:=`(
            statistic = coef_summary[, stat_col],
            p_value = coef_summary[, ncol(coef_summary)]
        )]
        
        ## Add QC stats if requested
        if (keep_qc_stats) {
            if (model_class %in% c("glm", "negbin")) {
                dt[, `:=`(
                    AIC = stats::AIC(model),
                    BIC = stats::BIC(model),
                    deviance = stats::deviance(model),
                    null_deviance = if (!isS4(model)) model$null.deviance else NA,
                    df_residual = stats::df.residual(model)
                )]
                
                ## Add R-squared for non-binomial GLMs
                if (model_class == "glm" && !isS4(model) && model$family$family != "binomial") {
                    dt[, R2 := 1 - (deviance / null_deviance)]
                }
                
                ## Add theta for negative binomial
                if (model_class == "negbin") {
                    dt[, theta := model$theta]
                }
                
                ## For binomial, add discrimination/calibration metrics
                if (is_logistic && keep_qc_stats) {
                    ## C-statistic (if pROC available)
                    if (requireNamespace("pROC", quietly = TRUE)) {
                        if (!isS4(model)) {
                            roc_obj <- pROC::roc(model$y, stats::fitted(model), quiet = TRUE)
                            dt[, c_statistic := as.numeric(pROC::auc(roc_obj))]
                        }
                    }
                    
                    ## Hosmer-Lemeshow test (if ResourceSelection available)
                    if (requireNamespace("ResourceSelection", quietly = TRUE)) {
                        if (!isS4(model)) {
                            hl <- ResourceSelection::hoslem.test(model$y, stats::fitted(model), g = 10)
                            dt[, `:=`(
                                hoslem_chi2 = hl$statistic,
                                hoslem_p = hl$p.value
                            )]
                        }
                    }
                }
            } else if (model_class == "lm") {
                summ <- summary(model)
                dt[, `:=`(
                    R2 = summ$r.squared,
                    adj_R2 = summ$adj.r.squared,
                    AIC = stats::AIC(model),
                    BIC = stats::BIC(model),
                    sigma = summ$sigma,
                    df_residual = stats::df.residual(model)
                )]
            }
        }
        
    } else if (model_class %in% c("coxph", "clogit")) {
        
        if (!requireNamespace("survival", quietly = TRUE))
            stop("Package 'survival' required")
        
        summ <- summary(model)
        coef_summary <- summ$coefficients
        conf_int <- stats::confint(model, level = conf_level)
        
        dt <- data.table::data.table(
                              model_scope = model_scope %||% "Multivariable",
                              model_type = model_type_name,
                              term = rownames(coef_summary),
                              n = if (!is.null(model$n)) model$n[1] else summ$n,
                              events = if (!is.null(model$nevent)) model$nevent 
                                       else if (!is.null(model$n)) model$n[2] 
                                       else summ$nevent,
                              coefficient = coef_summary[, "coef"],
                              se = coef_summary[, "se(coef)"],
                              ## Store both versions
                              coef = coef_summary[, "coef"],
                              coef_lower = conf_int[, 1],
                              coef_upper = conf_int[, 2],
                              exp_coef = coef_summary[, "exp(coef)"],
                              exp_lower = exp(conf_int[, 1]),
                              exp_upper = exp(conf_int[, 2]),
                              ## Primary display columns
                              HR = coef_summary[, "exp(coef)"],
                              ci_lower = exp(conf_int[, 1]),
                              ci_upper = exp(conf_int[, 2]),
                              statistic = coef_summary[, "z"],
                              p_value = coef_summary[, "Pr(>|z|)"]
                          )

        ## Add QC stats
        if (keep_qc_stats) {
            dt[, `:=`(
                concordance = summ$concordance[1],
                concordance_se = summ$concordance[2],
                rsq = summ$rsq[1],
                rsq_max = summ$rsq[2],
                likelihood_ratio_test = summ$logtest[1],
                likelihood_ratio_df = summ$logtest[2],
                likelihood_ratio_p = summ$logtest[3],
                wald_test = summ$waldtest[1],
                wald_df = summ$waldtest[2],
                wald_p = summ$waldtest[3],
                score_test = summ$sctest[1],
                score_df = summ$sctest[2],
                score_p = summ$sctest[3]
            )]
        }

    } else if (model_class %in% c("coxme", "lme", "lmer", "lmerMod", "glmer", "glmerMod", "merMod") || 
               inherits(model, c("lmerMod", "glmerMod", "merMod", "lmer", "glmer"))) {
        
        ## Mixed effects models
        summ <- summary(model)
        
        if (model_class == "coxme") {
            coef_vec <- coxme::fixef(model)
            vcov_mat <- as.matrix(stats::vcov(model))
            se_vec <- sqrt(diag(vcov_mat))
            z_vec <- coef_vec / se_vec
            p_vec <- 2 * (1 - stats::pnorm(abs(z_vec)))
            
            dt <- data.table::data.table(
                                  model_scope = model_scope %||% "Multivariable",
                                  model_type = model_type_name,
                                  term = names(coef_vec),
                                  n = model$n[2],
                                  events = model$n[1],
                                  n_group = NA_real_,
                                  events_group = NA_real_, 
                                  coefficient = coef_vec,
                                  se = se_vec,
                                  coef = coef_vec,
                                  coef_lower = coef_vec - stats::qnorm((1 + conf_level) / 2) * se_vec,
                                  coef_upper = coef_vec + stats::qnorm((1 + conf_level) / 2) * se_vec,
                                  exp_coef = exp(coef_vec),
                                  exp_lower = exp(coef_vec - stats::qnorm((1 + conf_level) / 2) * se_vec),
                                  exp_upper = exp(coef_vec + stats::qnorm((1 + conf_level) / 2) * se_vec),
                                  HR = exp(coef_vec),
                                  ci_lower = exp(coef_vec - stats::qnorm((1 + conf_level) / 2) * se_vec),
                                  ci_upper = exp(coef_vec + stats::qnorm((1 + conf_level) / 2) * se_vec),
                                  statistic = z_vec,
                                  p_value = p_vec
                              )
            
            ## Add QC stats for coxme
            if (keep_qc_stats) {
                ## AIC and BIC
                dt[, `:=`(
                    AIC = stats::AIC(model),
                    BIC = stats::BIC(model)
                )]
                
                ## Log-likelihood
                loglik <- model$loglik
                if (!is.null(loglik) && length(loglik) >= 2) {
                    dt[, `:=`(
                        loglik_null = loglik[1],
                        loglik_model = loglik[length(loglik)]
                    )]
                }
                
                ## Integrated log-likelihood (penalized)
                if (!is.null(model$loglik)) {
                    dt[, integrated_loglik := model$loglik[length(model$loglik)]]
                }
                
                ## Initialize c_statistic with NA
                dt[, c_statistic := NA_real_]
                
                ## Concordance/C-statistic for coxme
                ## coxme stores both the survival object (y) and linear predictor directly
                if (requireNamespace("survival", quietly = TRUE)) {
                    c_stat <- tryCatch({
                        ## Both y (Surv object) and linear.predictor are stored in coxme models
                        if (!is.null(model$y) && !is.null(model$linear.predictor)) {
                            conc <- survival::concordance(model$y ~ model$linear.predictor)
                            conc$concordance
                        } else {
                            NA_real_
                        }
                    }, error = function(e) NA_real_)
                    
                    if (!is.na(c_stat)) {
                        dt[, c_statistic := c_stat]
                    }
                }
            }
            
        } else {
            
            ## lmer/glmer from lme4
            if (!requireNamespace("lme4", quietly = TRUE))
                stop("Package 'lme4' required")
            
            coef_summary <- stats::coef(summ)
            
            ## Determine if should exponentiate
            if (model_class %in% c("lmer", "lmerMod")) {
                should_exp <- FALSE  ## Linear mixed models: no exponentiation
            } else if (model_class %in% c("glmer", "glmerMod")) {
                should_exp <- (summ$family == "binomial" || summ$link == "log")
            } else {
                should_exp <- FALSE  # Default for other mixed models
            }
            
            dt <- data.table::data.table(
                                  model_scope = model_scope %||% "Multivariable",
                                  model_type = model_type_name,
                                  term = rownames(coef_summary),
                                  n = stats::nobs(model),
                                  events = NA_real_,
                                  coefficient = coef_summary[, "Estimate"],
                                  se = coef_summary[, "Std. Error"]
                              )
            
            ## For glmer binomial, calculate total events from the response
            if (model_class %in% c("glmer", "glmerMod") && inherits(model, "merMod")) {
                if (summ$family == "binomial") {
                    response_var <- model@resp$y  
                    if (!is.null(response_var)) {
                        dt[, events := sum(response_var, na.rm = TRUE)]
                    }
                }
            }
            
            z_score <- stats::qnorm((1 + conf_level) / 2)
            dt[, `:=`(
                coef = coefficient,
                coef_lower = coefficient - z_score * se,
                coef_upper = coefficient + z_score * se,
                exp_coef = if (should_exp) exp(coefficient) else coefficient,
                exp_lower = if (should_exp) exp(coefficient - z_score * se) else coefficient - z_score * se,
                exp_upper = if (should_exp) exp(coefficient + z_score * se) else coefficient + z_score * se
            )]
            
            ## Add appropriate effect column
            if (model_class %in% c("glmer", "glmerMod") && summ$family == "binomial") {
                dt[, `:=`(
                    OR = exp_coef,
                    ci_lower = exp_lower,
                    ci_upper = exp_upper
                )]
            } else if (should_exp) {
                dt[, `:=`(
                    RR = exp_coef,
                    ci_lower = exp_lower,
                    ci_upper = exp_upper
                )]
            } else {
                dt[, `:=`(
                    Coefficient = coef,
                    ci_lower = coef_lower,
                    ci_upper = coef_upper
                )]
            }
            
            ## Add test statistics
            if (ncol(coef_summary) >= 3) {
                ## Use z-values if available
                stat_col <- if ("z value" %in% colnames(coef_summary)) {
                                "z value"
                            } else if ("t value" %in% colnames(coef_summary)) {
                                "t value"
                            } else {
                                NULL
                            }
                
                if (!is.null(stat_col)) {
                    dt[, statistic := coef_summary[, stat_col]]
                } else {
                    ## Calculate z-statistics manually
                    dt[, statistic := coefficient / se]
                }
                
                ## Check if p-values are provided
                p_col <- grep("^Pr\\(", colnames(coef_summary))
                if (length(p_col) > 0) {
                    dt[, p_value := coef_summary[, p_col[1]]]
                } else {
                    ## Calculate p-values from z/t statistics
                    dt[, p_value := 2 * (1 - stats::pnorm(abs(statistic)))]
                }
            } else {
                ## No test statistics - calculate manually
                dt[, `:=`(
                    statistic = coefficient / se,
                    p_value = 2 * (1 - stats::pnorm(abs(coefficient / se)))
                )]
            }
            
            ## Add QC stats for lmer/glmer
            if (keep_qc_stats) {
                ## AIC and BIC
                dt[, `:=`(
                    AIC = stats::AIC(model),
                    BIC = stats::BIC(model)
                )]
                
                ## Log-likelihood
                loglik_val <- as.numeric(stats::logLik(model))
                dt[, loglik := loglik_val]
                
                ## Deviance - handle REML vs ML fits
                dev_val <- tryCatch({
                    if (model_class %in% c("lmer", "lmerMod") && lme4::isREML(model)) {
                        ## For REML fits, use REMLcrit or deviance with REML=FALSE
                        stats::deviance(model, REML = FALSE)
                    } else {
                        stats::deviance(model)
                    }
                }, error = function(e) NA_real_)
                dt[, deviance := dev_val]
                
                ## For lmer models, add R-squared measures
                if (model_class %in% c("lmer", "lmerMod")) {
                    ## Initialize R-squared columns with NA
                    dt[, `:=`(
                        r_squared = NA_real_,
                        r_squared_marginal = NA_real_,
                        r_squared_conditional = NA_real_
                    )]
                    
                    ## Try to get R-squared using MuMIn if available
                    if (requireNamespace("MuMIn", quietly = TRUE)) {
                        r2_vals <- tryCatch({
                            MuMIn::r.squaredGLMM(model)
                        }, error = function(e) NULL)
                        
                        if (!is.null(r2_vals)) {
                            dt[, `:=`(
                                r_squared = r2_vals[1, "R2m"],
                                r_squared_marginal = r2_vals[1, "R2m"],
                                r_squared_conditional = r2_vals[1, "R2c"]
                            )]
                        }
                    }
                    
                    ## Residual standard deviation (sigma)
                    dt[, sigma := stats::sigma(model)]
                }
                
                ## For glmer binomial, try to compute C-statistic
                if (model_class %in% c("glmer", "glmerMod") && inherits(model, "merMod")) {
                    ## Initialize c_statistic with NA
                    dt[, c_statistic := NA_real_]
                    
                    if (summ$family == "binomial") {
                        if (requireNamespace("pROC", quietly = TRUE)) {
                            c_stat <- tryCatch({
                                fitted_vals <- stats::fitted(model)
                                response_vals <- model@resp$y
                                roc_obj <- pROC::roc(response_vals, fitted_vals, quiet = TRUE)
                                as.numeric(pROC::auc(roc_obj))
                            }, error = function(e) NA_real_)
                            
                            if (!is.na(c_stat)) {
                                dt[, c_statistic := c_stat]
                            }
                        }
                    }
                }
            }
        }
        
    } else {
        stop("Unsupported model class: ", model_class)
    }
    
    ## Parse terms into variable and group AFTER creating dt
    xlevels <- get_model_xlevels(model)
    
    ## Parse terms - pass model for coxme special handling
    parsed <- parse_term(dt$term, xlevels, model)
    dt[, `:=`(variable = parsed$variable, group = parsed$group)]
    
    ## Initialize n_group and events_group columns if absent
    if (!"n_group" %in% names(dt)) {
        dt[, n_group := NA_real_]
    }
    if (!"events_group" %in% names(dt)) {
        dt[, events_group := NA_real_]
    }

    ## Group counts calculation
    all_counts <- NULL
    data_dt <- NULL
    event_var <- NULL
    outcome_var <- NULL

    if (!skip_counts) {
        
        ## Prepare data for group counting
        if (!is.null(xlevels) || model_class == "coxme") {
            data_source <- model_data
            data_dt <- if (data.table::is.data.table(data_source)) data_source else data.table::as.data.table(data_source)
            
            ## For coxme, reconstruct xlevels from the data if needed
            if (is.null(xlevels) && model_class == "coxme") {
                ## Use formulaList$fixed for coxme
                formula_to_use <- model$formulaList$fixed
                
                ## Get variables
                all_formula_vars <- all.vars(formula_to_use)
                term_vars <- all_formula_vars[-1]  # Exclude response
                
                xlevels <- list()
                
                for (var_name in term_vars) {
                    if (var_name %in% names(data_dt) && is.factor(data_dt[[var_name]])) {
                        xlevels[[var_name]] <- levels(data_dt[[var_name]])
                    }
                }
            }
            
            ## Determine outcome and event variables
            outcome_var <- NULL
            event_var <- NULL
            
            if (model_class %in% c("glm", "negbin") && !isS4(model)) {
                outcome_var <- all.vars(model$formula)[1]
            } else if (model_class %in% c("lmer", "lmerMod", "glmer", "glmerMod") && inherits(model, "merMod")) {
                ## For mixed models, get the response variable name from the formula
                formula_obj <- model@call$formula
                if (!is.null(formula_obj)) {
                    outcome_var <- all.vars(formula_obj)[1]
                }
            } else if (model_class %in% c("coxph", "clogit", "coxme")) {
                ## For all survival models, use the unified helper function
                event_var <- get_event_variable(model, model_class)
            }
            
            ## Calculate n_group and events_group for all factor variables at once
            if (length(xlevels) > 0 && !is.null(data_dt)) {
                
                ## Get factor variables that exist in the data
                factor_vars <- names(xlevels)[names(xlevels) %in% names(data_dt)]
                
                if (length(factor_vars) > 0) {
                    
                    ## Use melt for efficient long-format conversion (faster than lapply+rbindlist)
                    ## Only select the columns we need to minimize memory usage
                    cols_needed <- factor_vars
                    if (!is.null(event_var) && event_var %in% names(data_dt)) {
                        cols_needed <- c(cols_needed, event_var)
                    } else if (!is.null(outcome_var) && outcome_var %in% names(data_dt)) {
                        cols_needed <- c(cols_needed, outcome_var)
                    }
                    
                    ## Melt to long format - much faster than lapply for multiple variables
                    long_dt <- data.table::melt(
                                               data_dt[, ..cols_needed],
                                               measure.vars = factor_vars,
                                               variable.name = "variable",
                                               value.name = "group",
                                               na.rm = TRUE,
                                               variable.factor = FALSE
                                           )
                    long_dt[, group := as.character(group)]
                    
                    if (!is.null(event_var) && event_var %in% names(data_dt)) {
                        ## For survival models
                        all_counts <- long_dt[, .(
                            n_group = .N,
                            events_group = sum(get(event_var), na.rm = TRUE)
                        ), by = .(variable, group)]
                        
                        ## Single join to update all counts at once
                        dt[all_counts, `:=`(
                                           n_group = i.n_group,
                                           events_group = i.events_group
                                       ), on = .(variable, group)]
                        
                    } else if (!is.null(outcome_var) && outcome_var %in% names(data_dt)) {
                        ## For GLM/GLMER models
                        outcome_col <- long_dt[[outcome_var]]
                        if (is.factor(outcome_col)) {
                            long_dt[, .events_calc := as.numeric(outcome_col) == 2]
                        } else {
                            long_dt[, .events_calc := outcome_col]
                        }
                        
                        all_counts <- long_dt[, .(
                            n_group = .N,
                            events_group = sum(.events_calc, na.rm = TRUE)
                        ), by = .(variable, group)]
                        
                        ## Single join to update all counts
                        dt[all_counts, `:=`(
                                           n_group = i.n_group,
                                           events_group = i.events_group
                                       ), on = .(variable, group)]
                        
                    } else {
                        ## No outcome/event variable: just count n
                        all_counts <- long_dt[, .(n_group = .N), by = .(variable, group)]
                        
                        ## Single join to update counts
                        dt[all_counts, n_group := i.n_group, on = .(variable, group)]
                    }
                }
            }
        }
        
        ## Calculate n_group and events_group for interaction terms
        ## Interactions are identified by ":" in the variable name
        interaction_rows <- grep(":", dt$variable, fixed = TRUE)
        
        if (length(interaction_rows) > 0 && !is.null(model_data)) {
            data_dt <- if (data.table::is.data.table(data_source)) data_dt else data.table::as.data.table(data_dt)
            
            ## Get xlevels if not already available
            if (is.null(xlevels)) {
                xlevels <- get_model_xlevels(model)
            }
            
            ## Determine event variable for survival models
            event_var <- NULL
            outcome_var <- NULL
            
            if (model_class %in% c("coxph", "clogit", "coxme")) {
                event_var <- get_event_variable(model, model_class)
            } else if (model_class %in% c("glm", "negbin") && !isS4(model)) {
                outcome_var <- all.vars(model$formula)[1]
            } else if (model_class %in% c("lmer", "lmerMod", "glmer", "glmerMod") && inherits(model, "merMod")) {
                formula_obj <- model@call$formula
                if (!is.null(formula_obj)) {
                    outcome_var <- all.vars(formula_obj)[1]
                }
            }
            
            ## Calculate counts for each interaction term
            for (row_idx in interaction_rows) {
                term_name <- dt$term[row_idx]
                
                ## Split interaction term into components (e.g., "treatmentDrug A:stageII" -> c("treatmentDrug A", "stageII"))
                components <- strsplit(term_name, ":", fixed = TRUE)[[1]]
                
                ## Parse each component into variable and level
                conditions <- list()
                valid_interaction <- TRUE
                
                for (comp in components) {
                    ## Try to match against xlevels to find variable name and level
                    found <- FALSE
                    
                    if (!is.null(xlevels)) {
                        for (var_name in names(xlevels)) {
                            if (startsWith(comp, var_name)) {
                                level_val <- substring(comp, nchar(var_name) + 1)
                                if (nchar(level_val) > 0 && var_name %in% names(data_dt)) {
                                    conditions[[length(conditions) + 1]] <- list(var = var_name, level = level_val)
                                    found <- TRUE
                                    break
                                }
                            }
                        }
                    }
                    
                    ## If not found in xlevels, treat as continuous variable
                    if (!found) {
                        ## Check if the component is a known variable name (continuous)
                        if (comp %in% names(data_dt)) {
                            ## For continuous variables, cannot filter by level
                            ## Just mark it as continuous (no filtering needed for count)
                            conditions[[length(conditions) + 1]] <- list(var = comp, level = NULL, continuous = TRUE)
                            found <- TRUE
                        }
                    }
                    
                    if (!found) {
                        valid_interaction <- FALSE
                        break
                    }
                }
                
                ## Calculate counts if all components were parsed successfully
                if (valid_interaction && length(conditions) > 0) {
                    ## Build filter expression for all factor conditions
                    filter_expr_parts <- character()
                    
                    for (cond in conditions) {
                        if (!isTRUE(cond$continuous) && !is.null(cond$level)) {
                            ## Factor variable - filter by level
                            filter_expr_parts <- c(filter_expr_parts, 
                                                   sprintf("get('%s') == '%s'", cond$var, cond$level))
                        }
                    }
                    
                    if (length(filter_expr_parts) > 0) {
                        ## Create combined filter expression
                        filter_expr <- paste(filter_expr_parts, collapse = " & ")
                        
                        ## Calculate n_group
                        n_val <- tryCatch({
                            filtered_dt <- data_dt[eval(parse(text = filter_expr))]
                            nrow(filtered_dt)
                        }, error = function(e) NA_real_)
                        
                        ## Calculate events_group if applicable
                        events_val <- NA_real_
                        
                        if (!is.null(event_var) && event_var %in% names(data_dt)) {
                            events_val <- tryCatch({
                                filtered_dt <- data_dt[eval(parse(text = filter_expr))]
                                sum(filtered_dt[[event_var]], na.rm = TRUE)
                            }, error = function(e) NA_real_)
                        } else if (!is.null(outcome_var) && outcome_var %in% names(data_dt)) {
                            events_val <- tryCatch({
                                filtered_dt <- data_dt[eval(parse(text = filter_expr))]
                                outcome_col <- filtered_dt[[outcome_var]]
                                if (is.factor(outcome_col)) {
                                    sum(as.numeric(outcome_col) == 2, na.rm = TRUE)
                                } else {
                                    sum(outcome_col, na.rm = TRUE)
                                }
                            }, error = function(e) NA_real_)
                        }
                        
                        ## Update the dt row
                        if (!is.na(n_val)) {
                            dt[row_idx, n_group := n_val]
                        }
                        if (!is.na(events_val)) {
                            dt[row_idx, events_group := events_val]
                        }
                    }
                }
            }
        }
    }
    
    ## Add reference rows for factor variables while maintaining original order
    ## Create all reference rows at once (vectorized approach)
    xlevels_ref <- xlevels
    
    if (!is.null(xlevels_ref) && reference_rows && length(xlevels_ref) > 0) {

        ## Add reference column to existing data
        dt[, reference := ""]
        
        ## Build reference counts data.table from all_counts if available
        ref_vars <- names(xlevels_ref)
        ref_levels <- vapply(xlevels_ref, `[`, character(1), 1)
        
        ## Create a lookup table for reference level counts
        ## Check all_counts more carefully - it may not exist or may be NULL
        have_all_counts <- !is.null(all_counts) && 
            inherits(all_counts, "data.table") && 
            nrow(all_counts) > 0 &&
            all(c("variable", "group") %in% names(all_counts))
        
        if (have_all_counts) {
            ## Build reference lookup from all_counts in one operation
            ref_lookup <- data.table::data.table(
                                          variable = ref_vars,
                                          group = ref_levels
                                      )
            ref_lookup <- all_counts[ref_lookup, on = .(variable, group)]
        } else if (!skip_counts && !is.null(data_dt) && inherits(data_dt, "data.table")) {
            ## Fallback: calculate reference counts directly
            ref_counts_list <- lapply(seq_along(ref_vars), function(i) {
                var <- ref_vars[i]
                ref_level <- ref_levels[i]
                if (var %in% names(data_dt)) {
                    if (!is.null(event_var) && event_var %in% names(data_dt)) {
                        data_dt[get(var) == ref_level & !is.na(get(var)), .(
                                                                              variable = var,
                                                                              group = ref_level,
                                                                              n_group = .N,
                                                                              events_group = sum(get(event_var), na.rm = TRUE)
                                                                          )]
                    } else {
                        data_dt[get(var) == ref_level & !is.na(get(var)), .(
                                                                              variable = var,
                                                                              group = ref_level,
                                                                              n_group = .N
                                                                          )]
                    }
                } else {
                    NULL
                }
            })
            ref_lookup <- data.table::rbindlist(ref_counts_list[!vapply(ref_counts_list, is.null, logical(1))], fill = TRUE)
        } else {
            ref_lookup <- data.table::data.table(
                                          variable = ref_vars,
                                          group = ref_levels,
                                          n_group = NA_real_
                                      )
        }
        
        ## Ensure ref_lookup has required columns
        if (!inherits(ref_lookup, "data.table") || nrow(ref_lookup) == 0) {
            ref_lookup <- data.table::data.table(
                                          variable = ref_vars,
                                          group = ref_levels,
                                          n_group = NA_real_
                                      )
        }
        
        ## Get the column structure from dt
        dt_cols <- names(dt)
        
        ## Find which variables have matching rows in dt
        vars_in_dt <- ref_vars[ref_vars %in% unique(dt$variable)]
        
        if (length(vars_in_dt) > 0 && nrow(dt) > 0) {
            ## Get template values from first matching row of each variable
            template_dt <- dt[, .SD[1], by = variable][variable %in% vars_in_dt]
            
            ## Build reference rows data.table directly
            ref_rows_dt <- data.table::data.table(
                                           variable = vars_in_dt,
                                           group = ref_levels[match(vars_in_dt, ref_vars)],
                                           term = paste0(vars_in_dt, ref_levels[match(vars_in_dt, ref_vars)]),
                                           coefficient = 0,
                                           se = NA_real_,
                                           coef = 0,
                                           coef_lower = NA_real_,
                                           coef_upper = NA_real_,
                                           exp_coef = 1,
                                           exp_lower = NA_real_,
                                           exp_upper = NA_real_,
                                           statistic = NA_real_,
                                           p_value = NA_real_,
                                           reference = reference_label
                                       )
            
            ## Add counts from ref_lookup
            if (nrow(ref_lookup) > 0 && "n_group" %in% names(ref_lookup)) {
                ref_rows_dt <- ref_lookup[ref_rows_dt, on = .(variable, group)]
            } else {
                ref_rows_dt[, n_group := NA_real_]
                if ("events_group" %in% dt_cols) {
                    ref_rows_dt[, events_group := NA_real_]
                }
            }
            
            ## Copy metadata columns from template using a proper join
            ## This avoids issues with match() returning wrong-length vectors
            meta_cols <- c("model_scope", "model_type", "n", "events")
            meta_cols_to_add <- meta_cols[meta_cols %in% dt_cols & !meta_cols %in% names(ref_rows_dt)]
            
            if (length(meta_cols_to_add) > 0) {
                ## Select only the columns we need from template
                template_subset <- template_dt[, c("variable", meta_cols_to_add), with = FALSE]
                ## Ensure no duplicates in template
                template_subset <- unique(template_subset, by = "variable")
                ## Join to add metadata columns
                ref_rows_dt <- template_subset[ref_rows_dt, on = "variable"]
            }
            
            ## Update n and events with group-specific values where available
            if ("n_group" %in% names(ref_rows_dt) && "n" %in% names(ref_rows_dt)) {
                ref_rows_dt[!is.na(n_group), n := n_group]
            }
            if ("events_group" %in% names(ref_rows_dt) && "events" %in% names(ref_rows_dt)) {
                ref_rows_dt[!is.na(events_group), events := events_group]
            }
            
            ## Add model-specific effect columns
            if ("HR" %in% dt_cols) {
                ref_rows_dt[, `:=`(HR = 1, ci_lower = NA_real_, ci_upper = NA_real_)]
            }
            if ("OR" %in% dt_cols) {
                ref_rows_dt[, `:=`(OR = 1, ci_lower = NA_real_, ci_upper = NA_real_)]
            }
            if ("RR" %in% dt_cols) {
                ref_rows_dt[, `:=`(RR = 1, ci_lower = NA_real_, ci_upper = NA_real_)]
            }
            if ("Coefficient" %in% dt_cols) {
                ref_rows_dt[, `:=`(Coefficient = 0, ci_lower = NA_real_, ci_upper = NA_real_)]
            }
            
            ## Add any QC stat columns that exist in dt
            qc_cols <- intersect(dt_cols, c("AIC", "BIC", "deviance", "null_deviance", 
                                            "loglik", "c_statistic", "rsq", "rsq_adj"))
            for (col in qc_cols) {
                if (!(col %in% names(ref_rows_dt))) {
                    ref_rows_dt[, (col) := template_dt[[col]][match(variable, template_dt$variable)]]
                }
            }
            
            ## Combine main dt and reference rows
            dt <- data.table::rbindlist(list(dt, ref_rows_dt), use.names = TRUE, fill = TRUE)
            
            ## Sort to put reference rows in correct positions
            ## Within each variable, reference row comes first
            dt[, .is_ref := reference == reference_label]
            data.table::setorder(dt, variable, -.is_ref, group)
            dt[, .is_ref := NULL]
        }
    }

    ## Update n and events with group-specific counts where available  
    dt[!is.na(n_group), n := n_group]
    if ("events_group" %in% names(dt)) {
        dt[!is.na(events_group), events := events_group]
    }
    
    ## Filter excluded terms
    if (!is.null(terms_to_exclude)) {
        dt <- dt[!term %in% terms_to_exclude]
    }
    
    ## Add significance indicators
    dt[, `:=`(
        sig = data.table::fcase(
                              is.na(p_value), "",
                              p_value < 0.001, "***",
                              p_value < 0.01, "**",
                              p_value < 0.05, "*",
                              p_value < 0.1, ".",
                              default = ""
                          ),
        sig_binary = !is.na(p_value) & p_value < 0.05
    )]

    ## Reorder variables to match original predictor order if available
    predictors_order <- attr(model, "predictors")
    if (!is.null(predictors_order) && "variable" %in% names(dt)) {
        
        ## Remove random effects notation - vectorized
        clean_predictors <- predictors_order[!grepl("\\|", predictors_order)]
        
        ## Get unique variables in dt
        dt_vars <- unique(dt$variable)
        
        ## Separate interaction terms from main effects
        interaction_vars <- dt_vars[grepl(":", dt_vars, fixed = TRUE)]
        main_effect_vars <- setdiff(dt_vars, interaction_vars)
        
        ## Build ordered_vars using vectorized matching
        ## First, direct matches
        direct_matches <- clean_predictors[clean_predictors %in% main_effect_vars]
        
        ## Then, factor variable matches (where dt_vars start with predictor name)
        remaining_predictors <- setdiff(clean_predictors, direct_matches)
        factor_matches <- character()
        if (length(remaining_predictors) > 0) {
            ## For each remaining predictor, find variables that start with it
            for (pred in remaining_predictors) {
                matches <- main_effect_vars[startsWith(main_effect_vars, pred)]
                factor_matches <- c(factor_matches, setdiff(matches, c(direct_matches, factor_matches)))
            }
        }
        
        ## Combine: direct matches first, then factor matches, preserving order
        ordered_vars <- c(direct_matches, factor_matches)
        
        ## Get main effects that aren't already ordered
        unordered_main <- setdiff(main_effect_vars, ordered_vars)
        
        ## Final order: ordered variables, unordered main effects, interactions
        final_order <- unique(c(ordered_vars, unordered_main, interaction_vars))
        
        ## Use factor levels for efficient sorting
        dt[, variable := factor(variable, levels = final_order)]
        data.table::setkey(dt, variable)
        dt[, variable := as.character(variable)]
        
        ## If reference rows exist, ensure they come first within each variable
        if ("reference" %in% names(dt)) {
            dt[, .temp_var_order := match(variable, final_order)]
            dt[, .is_ref := reference != ""]
            data.table::setorder(dt, .temp_var_order, -.is_ref, group)
            dt[, c(".temp_var_order", ".is_ref") := NULL]
        }
    }
    
    ## Add attributes
    data.table::setattr(dt, "model_class", model_class)
    data.table::setattr(dt, "formula_str", deparse(stats::formula(model)))
    
    if (model_class == "glm") {
        data.table::setattr(dt, "model_family", model$family$family)
        data.table::setattr(dt, "model_link", model$family$link)
    }
    
    if (model_class == "negbin") {
        data.table::setattr(dt, "model_family", "Negative Binomial")
        data.table::setattr(dt, "model_link", "log")
        data.table::setattr(dt, "theta", model$theta)
    }
    
    dt[]
    return(dt)
}
