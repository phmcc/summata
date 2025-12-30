#' Multivariate Regression Analysis
#'
#' Performs regression analyses of a single predictor (exposure) across multiple 
#' outcomes. This function is designed for studies where a single exposure variable 
#' is tested against multiple endpoints, such as complication screening, biomarker 
#' associations, or phenome-wide association studies. Returns publication-ready 
#' formatted results with optional covariate adjustment. Supports interactions,
#' mixed effects models, stratification, and clustered standard errors.
#'
#' @param data A data.frame or data.table containing the analysis dataset. The 
#'   function automatically converts data.frames to data.tables for efficient 
#'   processing.
#'   
#' @param outcomes Character vector of outcome variable names to analyze. Each 
#'   outcome is tested in its own model with the predictor. For survival analysis, 
#'   use \code{Surv()} syntax (e.g., \code{c("Surv(time1, status1)", "Surv(time2, status2)")}).
#'   
#' @param predictor Character string specifying the predictor (exposure) variable 
#'   name. This variable is tested against each outcome. Can be continuous or 
#'   categorical (factor).
#'   
#' @param covariates Optional character vector of covariate variable names to 
#'   include in adjusted models. When specified, models are fit as 
#'   \code{outcome ~ predictor + covariate1 + covariate2 + ...}, and only the 
#'   predictor effect is reported. Default is \code{NULL} (unadjusted models).
#'
#' @param interactions Optional character vector of interaction terms to include
#'   in adjusted models, using colon notation (e.g., \code{c("predictor:sex")}).
#'   Interactions involving the predictor will have their effects extracted and
#'   reported. Default is \code{NULL}.
#'
#' @param random Optional character string specifying random effects formula for
#'   mixed effects models (e.g., \code{"(1|hospital)"} or \code{"(1|site/patient)"}).
#'   Required when \code{model_type} is \code{"glmer"}, \code{"lmer"}, or 
#'   \code{"coxme"}. Default is \code{NULL}.
#'
#' @param strata Optional character string naming the stratification variable for
#'   Cox or conditional logistic models. Creates separate baseline hazards for 
#'   each stratum. Default is \code{NULL}.
#'
#' @param cluster Optional character string naming the clustering variable for
#'   Cox models. Computes robust clustered standard errors. Default is \code{NULL}.
#'   
#' @param model_type Character string specifying the type of regression model to 
#'   fit. Options include:
#'   \itemize{
#'     \item \code{"glm"} - Generalized linear model [default]
#'     \item \code{"lm"} - Linear regression
#'     \item \code{"coxph"} - Cox proportional hazards (survival analysis)
#'     \item \code{"clogit"} - Conditional logistic regression
#'     \item \code{"glmer"} - Generalized linear mixed effects model (requires \code{random})
#'     \item \code{"lmer"} - Linear mixed effects model (requires \code{random})
#'     \item \code{"coxme"} - Mixed effects Cox model (requires \code{random})
#'   }
#'   
#' @param family For GLM models, the error distribution and link function. Common 
#'   options include:
#'   \itemize{
#'     \item \code{"binomial"} - Logistic regression for binary outcomes [default]
#'     \item \code{"poisson"} - Poisson regression for count data
#'     \item \code{"gaussian"} - Normal linear regression via GLM
#'   }
#'   See \code{\link[stats]{family}} for all options. Ignored for non-GLM models.
#'   
#' @param columns Character string specifying which result columns to display when 
#'   both unadjusted and adjusted models are fit (i.e., when \code{covariates} is 
#'   specified):
#'   \itemize{
#'     \item \code{"adjusted"} - Show only adjusted (covariate-controlled) results [default]
#'     \item \code{"unadjusted"} - Show only unadjusted (crude) results
#'     \item \code{"both"} - Show both unadjusted and adjusted results side-by-side
#'   }
#'   Ignored when \code{covariates = NULL}.
#'   
#' @param p_threshold Numeric value between 0 and 1 specifying a p-value threshold 
#'   for filtering results. Only outcomes with p-value <= threshold are included 
#'   in the output. Default is 1 (no filtering, all outcomes returned).
#'   
#' @param conf_level Numeric confidence level for confidence intervals. Must be 
#'   between 0 and 1. Default is 0.95 (95\% confidence intervals).
#'   
#' @param show_n Logical. If \code{TRUE}, includes the sample size column in 
#'   the output table. Default is \code{TRUE}.
#'   
#' @param show_events Logical. If \code{TRUE}, includes the events column in the 
#'   output table (relevant for survival and logistic regression). Default is 
#'   \code{TRUE}.
#'   
#' @param digits Integer specifying the number of decimal places for effect 
#'   estimates (OR, HR, RR, coefficients). Default is 2.
#'   
#' @param p_digits Integer specifying the number of decimal places for p-values. 
#'   P-values smaller than \code{10^(-p_digits)} are displayed as "< 0.001", 
#'   etc. Default is 3.
#'   
#' @param labels Named character vector or list providing custom display 
#'   labels for variables. Can include labels for outcomes, predictors, and 
#'   covariates. Names should match variable names, values are the display labels.
#'   Labels are applied to: (1) outcome names in the Outcome column, (2) predictor
#'   variable name when displayed, and (3) variable names in formatted interaction
#'   terms. Variables not in \code{labels} use their original names. 
#'   Default is \code{NULL}.
#'   
#' @param predictor_label Optional character string providing a custom display 
#'   label for the predictor variable. Takes precedence over \code{labels} for
#'   the predictor. Default is \code{NULL} (uses label from \code{labels} or 
#'   original name).
#'
#' @param include_predictor Logical. If \code{TRUE} (default), includes the 
#'   Predictor column in the output table showing which level of a factor 
#'   predictor is being compared. If \code{FALSE}, omits the Predictor column,
#'   which may be useful when the predictor information will be explained in 
#'   a table caption or figure legend.
#'   
#' @param keep_models Logical. If \code{TRUE}, stores all fitted model objects 
#'   in the output as an attribute. Models are accessible via 
#'   \code{attr(result, "models")}. Default is \code{FALSE}.
#'   
#' @param exponentiate Logical. Whether to exponentiate coefficients (display 
#'   OR/HR/RR instead of log odds/log hazards). Default is \code{NULL}, which 
#'   automatically exponentiates for logistic, Poisson, and Cox models, and 
#'   displays raw coefficients for linear models.
#'   
#' @param parallel Logical. If \code{TRUE}, fits models in parallel for improved 
#'   performance with many outcomes. Default is \code{TRUE}.
#'   
#' @param n_cores Integer. Number of cores for parallel processing. Default is 
#'   \code{NULL} (auto-detect: uses number of available cores - 1).
#'   
#' @param ... Additional arguments passed to the underlying model fitting functions.
#'
#' @return A data.table with S3 class \code{"multifit_result"} containing formatted 
#'   multivariate results. The table structure includes:
#'   \describe{
#'     \item{Outcome}{Character. Outcome variable name or custom label}
#'     \item{Predictor}{Character. For factor predictors: formatted as 
#'       "Variable (Level)" showing the level being compared to reference.
#'       For binary variables where the non-reference level is an affirmative
#'       value (Yes, 1, True, Present, Positive, +), shows just "Variable".
#'       For continuous predictors: the variable name. For interactions: the
#'       formatted interaction term (e.g., "Treatment (Drug A) × Sex (Male)")}
#'     \item{n}{Integer. Sample size used in the model (if \code{show_n = TRUE})}
#'     \item{Events}{Integer. Number of events (if \code{show_events = TRUE})}
#'     \item{OR/HR/RR/Coefficient (95\% CI)}{Character. Unadjusted effect 
#'       estimate with CI (if \code{columns = "unadjusted"} or \code{"both"})}
#'     \item{aOR/aHR/aRR/Adj. Coefficient (95\% CI)}{Character. Adjusted 
#'       effect estimate with CI (if \code{columns = "adjusted"} or \code{"both"})}
#'     \item{Uni p / Multi p / p-value}{Character. Formatted p-value(s). Column 
#'       names depend on \code{columns} setting}
#'   }
#'   
#'   The returned object includes the following attributes accessible via \code{attr()}:
#'   \describe{
#'     \item{raw_data}{data.table. Unformatted numeric results with separate 
#'       columns for effect estimates, standard errors, confidence intervals, 
#'       and p-values. Suitable for custom analysis or visualization}
#'     \item{models}{list (if \code{keep_models = TRUE}). Named list of fitted 
#'       model objects, with outcome names as list names. Each element contains
#'       \code{$unadjusted} and/or \code{$adjusted} models depending on settings}
#'     \item{predictor}{Character. The predictor variable name}
#'     \item{outcomes}{Character vector. The outcome variable names}
#'     \item{covariates}{Character vector or NULL. The covariate variable names}
#'     \item{interactions}{Character vector or NULL. The interaction terms}
#'     \item{random}{Character or NULL. The random effects formula}
#'     \item{strata}{Character or NULL. The stratification variable}
#'     \item{cluster}{Character or NULL. The clustering variable}
#'     \item{model_type}{Character. The regression model type used}
#'     \item{columns}{Character. Which columns were displayed}
#'     \item{analysis_type}{Character. "multi_outcome" to identify analysis type}
#'   }
#'
#' @details
#' \strong{Analysis Approach:}
#' 
#' The function implements a multivariate (multi-outcome) screening workflow that
#' inverts the typical regression paradigm:
#' \enumerate{
#'   \item For each outcome in \code{outcomes}, fits a separate model with the 
#'     predictor as the main exposure
#'   \item If \code{covariates} specified, fits adjusted model: 
#'     \code{outcome ~ predictor + covariates + interactions}
#'   \item Extracts only the predictor effect(s) from each model, ignoring 
#'     covariate coefficients
#'   \item Combines results into a single table for comparison across outcomes
#'   \item Optionally filters by p-value threshold
#' }
#' 
#' This is conceptually opposite to \code{\link{uniscreen}}, which tests multiple 

#' predictors against a single outcome. Use \code{multifit} when you have one 
#' exposure of interest and want to screen across multiple endpoints.
#' 
#' \strong{When to Use Multivariate Analysis:}
#' \itemize{
#'   \item \strong{Complication screening}: Test one exposure (e.g., operative time, 
#'     BMI, biomarker level) against multiple postoperative complications
#'   \item \strong{Treatment effects}: Test one treatment against multiple efficacy 
#'     and safety endpoints simultaneously
#'   \item \strong{Biomarker studies}: Test one biomarker against multiple clinical 
#'     outcomes to understand its prognostic value
#'   \item \strong{Phenome-wide association studies (PheWAS)}: Test genetic variants 
#'     or exposures against many phenotypes
#'   \item \strong{Risk factor profiling}: Understand how one risk factor relates 
#'     to a spectrum of outcomes
#' }
#' 
#' \strong{Handling Categorical Predictors:}
#' 
#' When the predictor is a factor variable with multiple levels:
#' \itemize{
#'   \item Each non-reference level gets its own row for each outcome
#'   \item Reference category is determined by factor level ordering
#'   \item The Predictor column shows "Variable (Level)" format 
#'     (e.g., "Treatment (Drug A)", "Treatment (Drug B)")
#'   \item For binary variables with affirmative non-reference levels 
#'     (Yes, 1, True, Present, Positive, +), shows just "Variable" 
#'     (e.g., "Diabetes" instead of "Diabetes (Yes)")
#'   \item Effect estimates compare each level to the reference
#' }
#' 
#' \strong{Adjusted vs. Unadjusted Results:}
#' 
#' When \code{covariates} is specified, the function fits both models but only 
#' extracts predictor effects:
#' \itemize{
#'   \item \code{columns = "adjusted"}: Reports only covariate-adjusted effects.
#'     Column labeled "Multivariable aOR/aHR" etc.
#'   \item \code{columns = "unadjusted"}: Reports only crude effects. Column 
#'     labeled "Univariable OR/HR" etc.
#'   \item \code{columns = "both"}: Reports both side-by-side. Useful for 
#'     identifying confounding (large change in effect) or independent effects 
#'     (similar estimates)
#' }
#' 
#' \strong{Interaction Terms:}
#' 
#' When \code{interactions} includes terms involving the predictor:
#' \itemize{
#'   \item Main effect of predictor is always reported
#'   \item Interaction effects are extracted and displayed with formatted names
#'   \item Format: "Variable (Level) * Variable (Level)" using asterisk notation
#'   \item Useful for testing effect modification (e.g., does treatment effect 
#'     differ by sex?)
#' }
#' 
#' \strong{Mixed Effects Models:}
#' 
#' For clustered or hierarchical data (e.g., patients within hospitals):
#' \itemize{
#'   \item Use \code{model_type = "glmer"} with \code{random = "(1|cluster)"} for 
#'     random intercept models
#'   \item Nested random effects: \code{random = "(1|site/patient)"}
#'   \item Crossed random effects: \code{random = "(1|site) + (1|doctor)"}
#'   \item For survival outcomes, use \code{model_type = "coxme"}
#' }
#' 
#' \strong{Stratification and Clustering (Cox models):}
#' 
#' For Cox proportional hazards models:
#' \itemize{
#'   \item \code{strata}: Creates separate baseline hazards for each stratum level.
#'     Use when hazards are non-proportional across strata but you don't want to 
#'     estimate stratum effects
#'   \item \code{cluster}: Computes robust (sandwich) standard errors accounting 
#'     for within-cluster correlation. Alternative to mixed effects when only 
#'     robust SEs are needed
#' }
#' 
#' \strong{P-value Filtering:}
#' 
#' The \code{p_threshold} parameter filters results after fitting all models:
#' \itemize{
#'   \item Only outcomes with p <= threshold are retained in output
#'   \item For factor predictors, outcome is kept if any level is significant
#'   \item Useful for focusing on significant associations in exploratory analyses
#'   \item Default is 1 (no filtering) - recommended for confirmatory analyses
#' }
#' 
#' \strong{Effect Measures by Model Type:}
#' \itemize{
#'   \item \strong{Logistic} (\code{model_type = "glm"}, \code{family = "binomial"}): 
#'     Odds ratios (OR/aOR)
#'   \item \strong{Cox} (\code{model_type = "coxph"}): Hazard ratios (HR/aHR)
#'   \item \strong{Poisson} (\code{model_type = "glm"}, \code{family = "poisson"}): 
#'     Rate ratios (RR/aRR)
#'   \item \strong{Linear} (\code{model_type = "lm"}): Coefficient estimates
#'   \item \strong{Mixed effects}: Same as fixed effects counterparts
#' }
#' 
#' \strong{Memory and Performance:}
#' 
#' \itemize{
#'   \item \code{parallel = TRUE} (default) uses multiple cores for faster fitting
#'   \item \code{keep_models = FALSE} (default) discards model objects to save memory
#'   \item For many outcomes, parallel processing provides substantial speedup
#'   \item Set \code{keep_models = TRUE} only when you need model diagnostics
#' }
#'
#' @seealso 
#' \code{\link{uniscreen}} for screening multiple predictors against one outcome,
#' \code{\link{multiforest}} for creating forest plots from multifit results,
#' \code{\link{fit}} for single-outcome regression with full coefficient output,
#' \code{\link{fullfit}} for complete univariable-to-multivariable workflow
#'
#' @examples
#' # Load example data
#' data(clintrial)
#' data(clintrial_labels)
#' 
#' # Example 1: Basic multivariate analysis (unadjusted)
#' # Test treatment effect on multiple binary outcomes
#' result1 <- multifit(
#'     data = clintrial,
#'     outcomes = c("surgery", "pfs_status", "os_status"),
#'     predictor = "treatment",
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(result1)
#' # Shows odds ratios comparing Drug A and Drug B to Control
#' 
#' # Example 2: Adjusted analysis with covariates
#' # Adjust for age, sex, and disease stage
#' result2 <- multifit(
#'     data = clintrial,
#'     outcomes = c("surgery", "pfs_status", "os_status"),
#'     predictor = "treatment",
#'     covariates = c("age", "sex", "stage"),
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(result2)
#' # Shows adjusted odds ratios (aOR)
#' 
#' # Example 3: Compare unadjusted and adjusted results
#' result3 <- multifit(
#'     data = clintrial,
#'     outcomes = c("surgery", "pfs_status", "os_status"),
#'     predictor = "treatment",
#'     covariates = c("age", "sex", "stage"),
#'     columns = "both",
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(result3)
#' # Useful for identifying confounding effects
#' 
#' # Example 4: Continuous predictor across outcomes
#' # Test age effect on multiple outcomes
#' result4 <- multifit(
#'     data = clintrial,
#'     outcomes = c("surgery", "pfs_status", "os_status"),
#'     predictor = "age",
#'     covariates = c("sex", "treatment", "stage"),
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(result4)
#' # One row per outcome for continuous predictor
#' 
#' # Example 5: Cox regression for survival outcomes
#' library(survival)
#' cox_result <- multifit(
#'     data = clintrial,
#'     outcomes = c("Surv(pfs_months, pfs_status)", 
#'                  "Surv(os_months, os_status)"),
#'     predictor = "treatment",
#'     covariates = c("age", "sex", "stage"),
#'     model_type = "coxph",
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(cox_result)
#' # Returns hazard ratios (HR/aHR)
#' 
#' # Example 6: Cox with stratification by site
#' cox_strat <- multifit(
#'     data = clintrial,
#'     outcomes = c("Surv(os_months, os_status)"),
#'     predictor = "treatment",
#'     covariates = c("age", "sex"),
#'     strata = "site",
#'     model_type = "coxph",
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(cox_strat)
#' 
#' # Example 7: Cox with clustered standard errors
#' cox_cluster <- multifit(
#'     data = clintrial,
#'     outcomes = c("Surv(os_months, os_status)"),
#'     predictor = "treatment",
#'     covariates = c("age", "sex", "stage"),
#'     cluster = "site",
#'     model_type = "coxph",
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(cox_cluster)
#' 
#' # Example 8: Interaction between predictor and covariate
#' # Test if treatment effect differs by sex
#' result_int <- multifit(
#'     data = clintrial,
#'     outcomes = c("surgery", "os_status"),
#'     predictor = "treatment",
#'     covariates = c("age", "sex", "stage"),
#'     interactions = c("treatment:sex"),
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(result_int)
#' # Shows main effects and interaction terms with * notation
#' 
#' # Example 9: Linear model for continuous outcomes
#' linear_result <- multifit(
#'     data = clintrial,
#'     outcomes = c("los_days", "biomarker_x"),
#'     predictor = "treatment",
#'     covariates = c("age", "sex"),
#'     model_type = "lm",
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(linear_result)
#' # Returns coefficient estimates, not ratios
#' 
#' # Example 10: Poisson regression for count outcomes
#' poisson_result <- multifit(
#'     data = clintrial,
#'     outcomes = c("los_days"),
#'     predictor = "treatment",
#'     covariates = c("age", "sex", "stage"),
#'     model_type = "glm",
#'     family = "poisson",
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(poisson_result)
#' # Returns rate ratios (RR)
#' 
#' # Example 11: Filter to significant results only
#' sig_results <- multifit(
#'     data = clintrial,
#'     outcomes = c("surgery", "pfs_status", "os_status"),
#'     predictor = "stage",
#'     p_threshold = 0.05,
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(sig_results)
#' # Only outcomes with significant associations shown
#' 
#' # Example 12: Custom outcome labels
#' result_labeled <- multifit(
#'     data = clintrial,
#'     outcomes = c("surgery", "pfs_status", "os_status"),
#'     predictor = "treatment",
#'     labels = c(
#'         surgery = "Surgical Resection",
#'         pfs_status = "Disease Progression",
#'         os_status = "Death",
#'         treatment = "Treatment Group"
#'     ),
#'     parallel = FALSE
#' )
#' print(result_labeled)
#' 
#' # Example 13: Keep models for diagnostics
#' result_models <- multifit(
#'     data = clintrial,
#'     outcomes = c("surgery", "os_status"),
#'     predictor = "treatment",
#'     covariates = c("age", "sex"),
#'     keep_models = TRUE,
#'     parallel = FALSE
#' )
#' 
#' # Access stored models
#' models <- attr(result_models, "models")
#' names(models)
#' 
#' # Get adjusted model for surgery outcome
#' surgery_model <- models$surgery$adjusted
#' summary(surgery_model)
#' 
#' # Example 14: Access raw numeric data
#' result <- multifit(
#'     data = clintrial,
#'     outcomes = c("surgery", "os_status"),
#'     predictor = "age",
#'     parallel = FALSE
#' )
#' 
#' # Get unformatted results for custom analysis
#' raw_data <- attr(result, "raw_data")
#' print(raw_data)
#' # Contains exp_coef, ci_lower, ci_upper, p_value, etc.
#' 
#' # Example 15: Hide sample size and event columns
#' result_minimal <- multifit(
#'     data = clintrial,
#'     outcomes = c("surgery", "os_status"),
#'     predictor = "treatment",
#'     show_n = FALSE,
#'     show_events = FALSE,
#'     parallel = FALSE
#' )
#' print(result_minimal)
#' 
#' # Example 16: Customize decimal places
#' result_digits <- multifit(
#'     data = clintrial,
#'     outcomes = c("surgery", "os_status"),
#'     predictor = "age",
#'     digits = 3,
#'     p_digits = 4,
#'     parallel = FALSE
#' )
#' print(result_digits)
#' 
#' # Example 17: Force coefficient display (no exponentiation)
#' result_coef <- multifit(
#'     data = clintrial,
#'     outcomes = c("surgery"),
#'     predictor = "age",
#'     exponentiate = FALSE,
#'     parallel = FALSE
#' )
#' print(result_coef)
#' 
#' # Example 18: Complete publication workflow
#' final_table <- multifit(
#'     data = clintrial,
#'     outcomes = c("surgery", "pfs_status", "os_status"),
#'     predictor = "treatment",
#'     covariates = c("age", "sex", "stage", "grade"),
#'     columns = "both",
#'     labels = clintrial_labels,
#'     digits = 2,
#'     p_digits = 3,
#'     parallel = FALSE
#' )
#' print(final_table)
#' 
#' # Export to various formats for publication
#' # table2pdf(final_table, "treatment_effects.pdf", 
#' #         caption = "Treatment effects across outcomes")
#' # table2docx(final_table, "treatment_effects.docx")
#' 
#' # Example 19: Create forest plot from results
#' result_forest <- multifit(
#'     data = clintrial,
#'     outcomes = c("surgery", "pfs_status", "os_status"),
#'     predictor = "treatment",
#'     covariates = c("age", "sex", "stage"),
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' 
#' # Create forest plot (requires ggplot2)
#' # p <- multiforest(result_forest, 
#' #                  title = "Treatment Effects Across Outcomes")
#' # print(p)
#'
#' @export
multifit <- function(data,
                     outcomes,
                     predictor,
                     covariates = NULL,
                     interactions = NULL,
                     random = NULL,
                     strata = NULL,
                     cluster = NULL,
                     model_type = "glm",
                     family = "binomial",
                     columns = "adjusted",
                     p_threshold = 1,
                     conf_level = 0.95,
                     show_n = TRUE,
                     show_events = TRUE,
                     digits = 2,
                     p_digits = 3,
                     labels = NULL,
                     predictor_label = NULL,
                     include_predictor = TRUE,
                     keep_models = FALSE,
                     exponentiate = NULL,
                     parallel = TRUE,
                     n_cores = NULL,
                     ...) {
    
    ## Validate inputs
    if (missing(data)) stop("'data' is required")
    if (missing(outcomes)) stop("'outcomes' is required")
    if (missing(predictor)) stop("'predictor' is required")
    if (length(predictor) != 1) stop("'predictor' must be a single variable name")
    if (length(outcomes) < 1) stop("'outcomes' must contain at least one outcome")
    
    ## Validate columns parameter
    columns <- match.arg(columns, c("adjusted", "unadjusted", "both"))
    
    ## Validate model_type for mixed effects
    mixed_types <- c("glmer", "lmer", "coxme")
    if (!is.null(random) && !(model_type %in% mixed_types)) {
        stop("'random' effects require model_type to be one of: ", 
             paste(mixed_types, collapse = ", "))
    }
    
    ## Validate strata/cluster usage
    if (!is.null(strata) && !(model_type %in% c("coxph", "clogit"))) {
        stop("'strata' is only supported for coxph and clogit models")
    }
    if (!is.null(cluster) && model_type != "coxph") {
        stop("'cluster' is only supported for coxph models")
    }
    
    ## If no covariates/interactions/random, columns doesn't matter - treat as unadjusted
    has_adjustment <- !is.null(covariates) || !is.null(interactions) || !is.null(random)
    if (!has_adjustment) {
        columns <- "unadjusted"
    }
    
    ## Convert to data.table once at the start
    ## Use .data to avoid any potential conflicts with base::data
    .data <- if (data.table::is.data.table(data)) {
                 data.table::copy(data)
             } else {
                 data.table::as.data.table(data)
             }
    
    ## Determine if predictor is a factor (for handling multiple levels)
    is_factor_predictor <- is.factor(.data[[predictor]]) || 
        is.character(.data[[predictor]])
    
    ## Build the terms to extract from models
    ## For interactions involving the predictor, we need to extract those too
    terms_to_extract <- predictor
    if (!is.null(interactions)) {
        ## Find interactions that involve the predictor
        predictor_interactions <- interactions[grepl(predictor, interactions, fixed = TRUE)]
        terms_to_extract <- c(terms_to_extract, predictor_interactions)
    }
    
    ## Define model fitting function for a single outcome
    fit_one_outcome <- function(outcome) {
        
        results_list <- list()
        models_list <- list()
        
        ## Build unadjusted formula (predictor only, no interactions in univariable)
        unadj_formula_str <- paste(outcome, "~", predictor)
        
        ## Add strata for Cox models in unadjusted
        if (!is.null(strata) && model_type %in% c("coxph", "clogit")) {
            unadj_formula_str <- paste(unadj_formula_str, "+ strata(", strata, ")")
        }
        
        unadj_formula <- stats::as.formula(unadj_formula_str)
        
        ## Fit unadjusted model if needed
        if (columns %in% c("unadjusted", "both")) {
            unadj_model <- tryCatch({
                switch(model_type,
                       "glm" = stats::glm(unadj_formula, data = .data, family = family, ...),
                       "lm" = stats::lm(unadj_formula, data = .data, ...),
                       "coxph" = {
                           if (!requireNamespace("survival", quietly = TRUE)) 
                               stop("Package 'survival' required for Cox models")
                           if (!is.null(cluster)) {
                               survival::coxph(unadj_formula, data = .data, 
                                               cluster = .data[[cluster]], ...)
                           } else {
                               survival::coxph(unadj_formula, data = .data, ...)
                           }
                       },
                       "clogit" = {
                           if (!requireNamespace("survival", quietly = TRUE))
                               stop("Package 'survival' required for conditional logistic regression")
                           survival::clogit(unadj_formula, data = .data, ...)
                       },
                       ## Mixed effects models - univariable still includes random effects
                       "glmer" = {
                           if (!requireNamespace("lme4", quietly = TRUE))
                               stop("Package 'lme4' required for mixed effects models")
                           if (is.null(random)) stop("'random' is required for glmer models")
                           unadj_me_formula <- stats::as.formula(
                                                          paste(outcome, "~", predictor, "+", random)
                                                      )
                           lme4::glmer(unadj_me_formula, data = .data, family = family, ...)
                       },
                       "lmer" = {
                           if (!requireNamespace("lme4", quietly = TRUE))
                               stop("Package 'lme4' required for mixed effects models")
                           if (is.null(random)) stop("'random' is required for lmer models")
                           unadj_me_formula <- stats::as.formula(
                                                          paste(outcome, "~", predictor, "+", random)
                                                      )
                           lme4::lmer(unadj_me_formula, data = .data, ...)
                       },
                       "coxme" = {
                           if (!requireNamespace("coxme", quietly = TRUE))
                               stop("Package 'coxme' required for mixed effects Cox models")
                           if (is.null(random)) stop("'random' is required for coxme models")
                           unadj_me_formula <- stats::as.formula(
                                                          paste(outcome, "~", predictor, "+", random)
                                                      )
                           coxme::coxme(unadj_me_formula, data = .data, ...)
                       },
                       stop("Unsupported model type: ", model_type)
                       )
            }, error = function(e) {
                warning("Model failed for outcome '", outcome, "': ", e$message)
                return(NULL)
            })
            
            if (!is.null(unadj_model)) {
                ## Extract predictor rows from model
                unadj_raw <- extract_predictor_effects(
                    model = unadj_model,
                    predictor = predictor,
                    outcome = outcome,
                    conf_level = conf_level,
                    adjusted = FALSE,
                    terms_to_extract = predictor  # Only main effect for univariable
                )
                results_list$unadjusted <- unadj_raw
                
                if (keep_models) {
                    models_list$unadjusted <- unadj_model
                }
            }
        }
        
        ## Fit adjusted model if needed
        if (has_adjustment && columns %in% c("adjusted", "both")) {
            ## Build adjusted formula
            adj_terms <- predictor
            if (!is.null(covariates)) {
                adj_terms <- c(adj_terms, covariates)
            }
            if (!is.null(interactions)) {
                adj_terms <- c(adj_terms, interactions)
            }
            
            adj_formula_str <- paste(outcome, "~", paste(adj_terms, collapse = " + "))
            
            ## Add random effects for mixed models
            if (!is.null(random) && model_type %in% mixed_types) {
                adj_formula_str <- paste(adj_formula_str, "+", random)
            }
            
            ## Add strata for Cox models
            if (!is.null(strata) && model_type %in% c("coxph", "clogit")) {
                adj_formula_str <- paste(adj_formula_str, "+ strata(", strata, ")")
            }
            
            adj_formula <- stats::as.formula(adj_formula_str)
            
            adj_model <- tryCatch({
                switch(model_type,
                       "glm" = stats::glm(adj_formula, data = .data, family = family, ...),
                       "lm" = stats::lm(adj_formula, data = .data, ...),
                       "coxph" = {
                           if (!requireNamespace("survival", quietly = TRUE)) 
                               stop("Package 'survival' required for Cox models")
                           if (!is.null(cluster)) {
                               survival::coxph(adj_formula, data = .data, 
                                               cluster = .data[[cluster]], ...)
                           } else {
                               survival::coxph(adj_formula, data = .data, ...)
                           }
                       },
                       "clogit" = {
                           if (!requireNamespace("survival", quietly = TRUE))
                               stop("Package 'survival' required for conditional logistic regression")
                           survival::clogit(adj_formula, data = .data, ...)
                       },
                       "glmer" = {
                           if (!requireNamespace("lme4", quietly = TRUE))
                               stop("Package 'lme4' required for mixed effects models")
                           lme4::glmer(adj_formula, data = .data, family = family, ...)
                       },
                       "lmer" = {
                           if (!requireNamespace("lme4", quietly = TRUE))
                               stop("Package 'lme4' required for mixed effects models")
                           lme4::lmer(adj_formula, data = .data, ...)
                       },
                       "coxme" = {
                           if (!requireNamespace("coxme", quietly = TRUE))
                               stop("Package 'coxme' required for mixed effects Cox models")
                           coxme::coxme(adj_formula, data = .data, ...)
                       },
                       stop("Unsupported model type: ", model_type)
                       )
            }, error = function(e) {
                warning("Adjusted model failed for outcome '", outcome, "': ", e$message)
                return(NULL)
            })
            
            if (!is.null(adj_model)) {
                adj_raw <- extract_predictor_effects(
                    model = adj_model,
                    predictor = predictor,
                    outcome = outcome,
                    conf_level = conf_level,
                    adjusted = TRUE,
                    terms_to_extract = terms_to_extract  # Include interaction terms
                )
                results_list$adjusted <- adj_raw
                
                if (keep_models) {
                    models_list$adjusted <- adj_model
                }
            }
        }
        
        list(
            results = results_list,
            models = if (keep_models) models_list else NULL
        )
    }
    
    ## Fit models (parallel or sequential)
    if (parallel && length(outcomes) > 1) {
        ## Determine number of cores
        if (is.null(n_cores)) {
            n_cores <- max(1L, parallel::detectCores() - 1L)
        }
        
        if (.Platform$OS.type == "windows") {
            ## Windows: use parLapply with PSOCK cluster
            cl <- parallel::makeCluster(n_cores)
            on.exit(parallel::stopCluster(cl), add = TRUE)
            
            ## Export required objects to workers
            parallel::clusterExport(cl, 
                                    varlist = c(".data", "predictor", "covariates", "interactions",
                                                "random", "strata", "cluster", "model_type", 
                                                "family", "conf_level", "columns", "keep_models",
                                                "has_adjustment", "mixed_types", "terms_to_extract",
                                                "extract_predictor_effects"),
                                    envir = environment()
                                    )
            
            ## Load required packages on workers
            parallel::clusterEvalQ(cl, {
                library(data.table)
                library(stats)
            })
            
            if (model_type %in% c("coxph", "clogit")) {
                parallel::clusterEvalQ(cl, library(survival))
            }
            if (model_type %in% c("glmer", "lmer")) {
                parallel::clusterEvalQ(cl, library(lme4))
            }
            if (model_type == "coxme") {
                parallel::clusterEvalQ(cl, library(coxme))
            }
            
            all_results <- parallel::parLapply(cl, outcomes, fit_one_outcome)
            
        } else {
            ## Unix/Mac: use mclapply (fork-based)
            all_results <- parallel::mclapply(
                                         outcomes, 
                                         fit_one_outcome,
                                         mc.cores = n_cores
                                     )
        }
    } else {
        ## Sequential processing
        all_results <- lapply(outcomes, fit_one_outcome)
    }
    names(all_results) <- outcomes
    
    ## Combine results based on columns parameter
    combined_raw <- combine_multifit_results(all_results, columns)
    
    ## Store models if requested
    if (keep_models) {
        models <- lapply(all_results, `[[`, "models")
        models <- models[!vapply(models, is.null, logical(1))]
    }
    
    ## Filter by p-value if requested
    if (p_threshold < 1 && nrow(combined_raw) > 0) {
        ## Determine which p-value column to use for filtering
        p_col <- if ("p_adj" %in% names(combined_raw)) "p_adj" else "p_value"
        combined_raw <- combined_raw[get(p_col) <= p_threshold | is.na(get(p_col))]
    }
    
    ## Format results for publication
    formatted <- format_multifit_table(
        combined_raw,
        columns = columns,
        show_n = show_n,
        show_events = show_events,
        digits = digits,
        p_digits = p_digits,
        labels = labels,
        predictor_label = predictor_label,
        include_predictor = include_predictor,
        exponentiate = exponentiate
    )
    
    ## Attach attributes
    data.table::setattr(formatted, "raw_data", combined_raw)
    
    if (keep_models) {
        data.table::setattr(formatted, "models", models)
    }
    
    data.table::setattr(formatted, "predictor", predictor)
    data.table::setattr(formatted, "outcomes", outcomes)
    data.table::setattr(formatted, "covariates", covariates)
    data.table::setattr(formatted, "interactions", interactions)
    data.table::setattr(formatted, "random", random)
    data.table::setattr(formatted, "strata", strata)
    data.table::setattr(formatted, "cluster", cluster)
    data.table::setattr(formatted, "model_type", model_type)
    data.table::setattr(formatted, "columns", columns)
    data.table::setattr(formatted, "include_predictor", include_predictor)
    data.table::setattr(formatted, "analysis_type", "multi_outcome")
    
    class(formatted) <- c("multifit_result", class(formatted))
    
    return(formatted)
}


#' Format interaction term for display
#' 
#' Converts R's internal interaction term format (e.g., "treatmentDrug A:stageII")
#' to a more readable format (e.g., "Treatment (Drug A) * Stage (II)").
#' 
#' @param term Character string of the interaction term from model coefficients.
#' @param labels Optional named vector of labels for variable names.
#' @return Formatted interaction term string.
#' @keywords internal
format_interaction_term <- function(term, labels = NULL) {
    ## Split by colon
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    
    formatted_parts <- vapply(parts, function(part) {
        ## Try to separate variable name from level
        ## Common patterns: varLevel, var.Level, varLevel1Level2
        
        ## First, check if this looks like a factor level (contains uppercase after lowercase)
        ## or ends with a number pattern
        
        ## Try to find where variable name ends and level begins
        ## Look for transition from lowercase to uppercase, or number after letter
        
        ## Strategy: find the longest matching variable name prefix
        ## For "treatmentDrug A", we want "treatment" and "Drug A"
        ## For "stageII", we want "stage" and "II"
        ## For "age" (continuous in interaction), we want just "age"
        
        ## Check for common variable name patterns
        ## Match: lowercase letters, possibly followed by underscore/lowercase
        var_match <- regexpr("^[a-z][a-z0-9_]*", part, perl = TRUE)
        
        if (var_match > 0) {
            var_len <- attr(var_match, "match.length")
            var_name <- substr(part, 1, var_len)
            level <- substr(part, var_len + 1, nchar(part))
            
            ## Apply label if available
            if (!is.null(labels) && var_name %in% names(labels)) {
                var_display <- labels[[var_name]]
            } else {
                ## Convert to title case
                var_display <- paste0(toupper(substr(var_name, 1, 1)), 
                                      substr(var_name, 2, nchar(var_name)))
            }
            
            if (nchar(level) > 0) {
                ## Has a level - format as "Variable (Level)"
                paste0(var_display, " (", level, ")")
            } else {
                ## Continuous variable - just the name
                var_display
            }
        } else {
            ## Fallback - return as-is
            part
        }
    }, character(1))
    
    ## Join with asterisk (consistent with other fullfit functions)
    paste(formatted_parts, collapse = " * ")
}


#' Extract predictor effects from a fitted model
#' 
#' Internal helper function that extracts only the predictor variable's 
#' coefficient(s) from a fitted model, ignoring intercept and covariates.
#' Supports standard models (glm, lm, coxph) and mixed effects models
#' (glmer, lmer, coxme).
#' 
#' @param model A fitted model object.
#' @param predictor Character string of the predictor variable name.
#' @param outcome Character string of the outcome variable name.
#' @param conf_level Numeric confidence level.
#' @param adjusted Logical indicating if this is an adjusted model.
#' @param terms_to_extract Character vector of terms to extract (predictor and 
#'   optionally interaction terms involving the predictor).
#' @return data.table with predictor effect information.
#' @keywords internal
extract_predictor_effects <- function(model, predictor, outcome, 
                                      conf_level = 0.95, adjusted = FALSE,
                                      terms_to_extract = NULL) {
    
    if (is.null(terms_to_extract)) {
        terms_to_extract <- predictor
    }
    
    model_classes <- class(model)
    model_class <- model_classes[1]
    
    ## Detect mixed effects models
    is_merMod <- inherits(model, "merMod") || any(c("glmerMod", "lmerMod") %in% model_classes)
    is_coxme <- inherits(model, "coxme")
    is_mixed <- is_merMod || is_coxme
    
    ## Get coefficient summary based on model type
    if (is_coxme) {
        ## coxme models
        coefs <- coxme::fixef(model)
        vcov_mat <- as.matrix(vcov(model))
        se_values <- sqrt(diag(vcov_mat))
        coef_names <- names(coefs)
        
        ## Build coefficient summary matrix
        z_vals <- coefs / se_values
        p_vals <- 2 * stats::pnorm(-abs(z_vals))
        coef_summary <- cbind(coef = coefs, `se(coef)` = se_values, 
                              z = z_vals, `Pr(>|z|)` = p_vals)
        rownames(coef_summary) <- coef_names
        
    } else if (is_merMod) {
        ## lme4 models (glmer, lmer)
        coef_summary <- summary(model)$coefficients
        coef_names <- rownames(coef_summary)
        
    } else {
        ## Standard models
        coef_summary <- summary(model)$coefficients
        coef_names <- rownames(coef_summary)
    }
    
    ## Find which coefficients belong to the terms we want to extract
    ## Build pattern to match predictor and any interaction terms
    predictor_rows <- integer(0)
    for (term in terms_to_extract) {
        if (grepl(":", term, fixed = TRUE)) {
            ## Interaction term - match exactly or with factor level suffixes
            ## Split interaction components
            components <- strsplit(term, ":", fixed = TRUE)[[1]]
            ## Build regex that matches the interaction in either order
            ## e.g., "age:sexMale" or "sexMale:age"
            pattern <- paste0("(^|:)", components[1], ".*:", components[2], "|",
                              "(^|:)", components[2], ".*:", components[1])
            matches <- grep(pattern, coef_names)
        } else {
            ## Main effect - match term at start of coefficient name
            pattern <- paste0("^", term)
            matches <- grep(pattern, coef_names)
            ## Exclude matches that are actually interaction terms (contain ":")
            ## unless we're specifically looking for interactions
            matches <- matches[!grepl(":", coef_names[matches], fixed = TRUE)]
        }
        predictor_rows <- c(predictor_rows, matches)
    }
    predictor_rows <- unique(predictor_rows)
    
    if (length(predictor_rows) == 0) {
        warning("Predictor '", predictor, "' not found in model for outcome '", outcome, "'")
        return(NULL)
    }
    
    ## Extract predictor coefficients
    pred_coefs <- coef_summary[predictor_rows, , drop = FALSE]
    
    ## Get sample size and events
    if (is_merMod) {
        n_obs <- nrow(model@frame)
    } else if (is_coxme) {
        n_obs <- model$n[1]
    } else {
        n_obs <- stats::nobs(model)
    }
    
    events <- NA_integer_
    if (model_class == "glm") {
        if (model$family$family == "binomial") {
            ## Binary outcome - count events (1s)
            if (!is.null(model$y)) {
                events <- sum(model$y, na.rm = TRUE)
            }
        } else if (model$family$family %in% c("poisson", "quasipoisson", "negbin")) {
            ## Count outcome - sum total events
            if (!is.null(model$y)) {
                events <- sum(model$y, na.rm = TRUE)
            }
        }
    } else if (model_class == "coxph") {
        events <- model$nevent
    } else if (is_coxme) {
        events <- sum(model$y[, "status"], na.rm = TRUE)
    } else if (inherits(model, "glmerMod")) {
        ## For glmer models - check family
        resp <- model@resp$y
        if (!is.null(resp)) {
            fam <- family(model)$family
            if (fam %in% c("binomial", "poisson", "quasipoisson")) {
                events <- sum(resp, na.rm = TRUE)
            }
        }
    }
    
    ## Build coefficient names (handle factor levels and interactions)
    term_names <- rownames(pred_coefs)
    group_names <- character(length(term_names))
    
    for (i in seq_along(term_names)) {
        term <- term_names[i]
        if (grepl(":", term, fixed = TRUE)) {
            ## Interaction term - store raw term, will be formatted later with labels
            group_names[i] <- term
        } else {
            ## Main effect - extract level name
            group_names[i] <- sub(paste0("^", predictor), "", term)
            if (group_names[i] == "") group_names[i] <- "-"  # Continuous predictor
        }
    }
    
    ## Get appropriate column names based on model class
    if (model_class == "coxph" || is_coxme) {
        coef_col <- "coef"
        se_col <- "se(coef)"
        stat_col <- "z"
        p_col <- if ("Pr(>|z|)" %in% colnames(pred_coefs)) "Pr(>|z|)" else "p"
    } else if (model_class == "glm" || inherits(model, "glmerMod")) {
        coef_col <- "Estimate"
        se_col <- "Std. Error"
        stat_col <- "z value"
        p_col <- "Pr(>|z|)"
    } else if (model_class == "lm" || inherits(model, "lmerMod")) {
        coef_col <- "Estimate"
        se_col <- "Std. Error"
        stat_col <- if ("t value" %in% colnames(pred_coefs)) "t value" else "z value"
        p_col <- if ("Pr(>|t|)" %in% colnames(pred_coefs)) "Pr(>|t|)" else "Pr(>|z|)"
    } else {
        coef_col <- colnames(pred_coefs)[1]
        se_col <- colnames(pred_coefs)[2]
        stat_col <- colnames(pred_coefs)[3]
        p_col <- colnames(pred_coefs)[4]
    }
    
    ## Calculate confidence intervals
    z_score <- stats::qnorm((1 + conf_level) / 2)
    coefficients <- pred_coefs[, coef_col]
    se_values <- pred_coefs[, se_col]
    ci_lower <- coefficients - z_score * se_values
    ci_upper <- coefficients + z_score * se_values
    
    ## Determine effect type and calculate exponentiated values
    if (model_class == "coxph" || is_coxme) {
        effect_type <- "HR"
        exp_coef <- exp(coefficients)
        exp_lower <- exp(ci_lower)
        exp_upper <- exp(ci_upper)
    } else if (model_class == "glm" || inherits(model, "glmerMod")) {
        family_name <- if (is_merMod) model@resp$family$family else model$family$family
        if (family_name == "binomial") {
            effect_type <- "OR"
            exp_coef <- exp(coefficients)
            exp_lower <- exp(ci_lower)
            exp_upper <- exp(ci_upper)
        } else if (family_name %in% c("poisson", "quasipoisson")) {
            effect_type <- "RR"
            exp_coef <- exp(coefficients)
            exp_lower <- exp(ci_lower)
            exp_upper <- exp(ci_upper)
        } else {
            effect_type <- "Coefficient"
            exp_coef <- coefficients
            exp_lower <- ci_lower
            exp_upper <- ci_upper
        }
    } else {
        effect_type <- "Coefficient"
        exp_coef <- coefficients
        exp_lower <- ci_lower
        exp_upper <- ci_upper
    }
    
    ## Handle p-values for lmer models (which may not have p-values)
    p_values <- if (p_col %in% colnames(pred_coefs)) {
                    pred_coefs[, p_col]
                } else {
                    rep(NA_real_, nrow(pred_coefs))
                }
    
    ## Handle statistics column
    stat_values <- if (stat_col %in% colnames(pred_coefs)) {
                       pred_coefs[, stat_col]
                   } else {
                       rep(NA_real_, nrow(pred_coefs))
                   }
    
    ## Build result data.table
    dt <- data.table::data.table(
                          outcome = outcome,
                          predictor = predictor,
                          group = group_names,
                          n = n_obs,
                          events = events,
                          coefficient = coefficients,
                          se = se_values,
                          exp_coef = exp_coef,
                          ci_lower = exp_lower,
                          ci_upper = exp_upper,
                          statistic = stat_values,
                          p_value = p_values,
                          effect_type = effect_type,
                          adjusted = adjusted
                      )
    
    return(dt)
}


#' Combine multivariate results
#' 
#' Internal helper to combine unadjusted and/or adjusted results from 
#' multiple outcomes into a single data.table.
#' 
#' @param all_results List of results from fit_one_outcome calls.
#' @param columns Character specifying which columns to include.
#' @return Combined data.table.
#' @keywords internal
combine_multifit_results <- function(all_results, columns) {
    
    ## Extract results based on columns parameter
    if (columns == "unadjusted") {
        results_list <- lapply(all_results, function(x) x$results$unadjusted)
        results_list <- results_list[!vapply(results_list, is.null, logical(1))]
        
        if (length(results_list) == 0) return(data.table::data.table())
        
        combined <- data.table::rbindlist(results_list, fill = TRUE)
        
    } else if (columns == "adjusted") {
        results_list <- lapply(all_results, function(x) x$results$adjusted)
        results_list <- results_list[!vapply(results_list, is.null, logical(1))]
        
        if (length(results_list) == 0) return(data.table::data.table())
        
        combined <- data.table::rbindlist(results_list, fill = TRUE)
        
    } else {
        ## columns == "both" - merge unadjusted and adjusted
        unadj_list <- lapply(all_results, function(x) x$results$unadjusted)
        adj_list <- lapply(all_results, function(x) x$results$adjusted)
        
        unadj_list <- unadj_list[!vapply(unadj_list, is.null, logical(1))]
        adj_list <- adj_list[!vapply(adj_list, is.null, logical(1))]
        
        if (length(unadj_list) == 0 && length(adj_list) == 0) {
            return(data.table::data.table())
        }
        
        unadj_combined <- if (length(unadj_list) > 0) {
                              data.table::rbindlist(unadj_list, fill = TRUE)
                          } else {
                              data.table::data.table()
                          }
        
        adj_combined <- if (length(adj_list) > 0) {
                            data.table::rbindlist(adj_list, fill = TRUE)
                        } else {
                            data.table::data.table()
                        }
        
        ## Rename columns for merging
        if (nrow(unadj_combined) > 0) {
            data.table::setnames(unadj_combined, 
                                 c("exp_coef", "ci_lower", "ci_upper", "p_value"),
                                 c("exp_coef_unadj", "ci_lower_unadj", "ci_upper_unadj", "p_unadj"),
                                 skip_absent = TRUE
                                 )
        }
        
        if (nrow(adj_combined) > 0) {
            data.table::setnames(adj_combined,
                                 c("exp_coef", "ci_lower", "ci_upper", "p_value"),
                                 c("exp_coef_adj", "ci_lower_adj", "ci_upper_adj", "p_adj"),
                                 skip_absent = TRUE
                                 )
        }
        
        ## Merge on outcome, predictor, group
        if (nrow(unadj_combined) > 0 && nrow(adj_combined) > 0) {
            merge_cols <- c("outcome", "predictor", "group")
            keep_from_unadj <- c(merge_cols, "n", "events", "effect_type",
                                 "exp_coef_unadj", "ci_lower_unadj", "ci_upper_unadj", "p_unadj")
            keep_from_adj <- c(merge_cols, "exp_coef_adj", "ci_lower_adj", "ci_upper_adj", "p_adj")
            
            unadj_subset <- unadj_combined[, intersect(keep_from_unadj, names(unadj_combined)), with = FALSE]
            adj_subset <- adj_combined[, intersect(keep_from_adj, names(adj_combined)), with = FALSE]
            
            combined <- merge(unadj_subset, adj_subset, by = merge_cols, all = TRUE)
        } else if (nrow(unadj_combined) > 0) {
            combined <- unadj_combined
        } else {
            combined <- adj_combined
        }
    }
    
    return(combined)
}


#' Format multifit results for publication
#' 
#' Internal helper that formats raw multivariate results into 
#' publication-ready table format.
#' 
#' @param data Raw combined data.table from combine_multifit_results.
#' @param columns Character specifying column layout.
#' @param show_n Logical for sample size column.
#' @param show_events Logical for events column.
#' @param digits Integer decimal places for effects.
#' @param p_digits Integer decimal places for p-values.
#' @param outcome_labels Named vector of outcome labels.
#' @param predictor_label Label for predictor variable.
#' @param exponentiate Logical or NULL for coefficient handling.
#' @return Formatted data.table.
#' @keywords internal
format_multifit_table <- function(data, 
                                  columns,
                                  show_n = TRUE,
                                  show_events = TRUE,
                                  digits = 2,
                                  p_digits = 3,
                                  labels = NULL,
                                  predictor_label = NULL,
                                  include_predictor = TRUE,
                                  exponentiate = NULL) {
    
    if (nrow(data) == 0) {
        return(data.table::data.table())
    }
    
    result <- data.table::copy(data)
    
    ## Apply outcome labels
    if (!is.null(labels) && length(labels) > 0) {
        for (orig_name in names(labels)) {
            result[outcome == orig_name, outcome := labels[[orig_name]]]
        }
    }
    
    ## Rename outcome column
    data.table::setnames(result, "outcome", "Outcome")
    
    ## Handle predictor group display
    ## If only one group value ("-"), omit Predictor column
    has_multiple_groups <- length(unique(result$group)) > 1 || 
        !all(result$group == "-")
    
    if (has_multiple_groups) {
        ## Get predictor name for labeling
        pred_name <- result$predictor[1]
        pred_display <- if (!is.null(predictor_label)) {
                            predictor_label
                        } else if (!is.null(labels) && pred_name %in% names(labels)) {
                            labels[[pred_name]]
                        } else {
                            pred_name
                        }
        
        data.table::setnames(result, "group", "Predictor")
        
        ## Check if this is a binary variable with affirmative non-reference level
        ## Common affirmative values that do not need "(Level)" suffix
        affirmative_values <- c("Yes", "YES", "yes", "1", "True", "TRUE", "true", 
                                "Present", "Positive", "+")
        unique_levels <- unique(result$Predictor[result$Predictor != "-"])
        is_binary_affirmative <- length(unique_levels) == 1 && 
            unique_levels[1] %in% affirmative_values
        
        ## Format predictor values as "Variable (Level)" or just "Variable" for binary
        for (i in seq_len(nrow(result))) {
            pred_val <- result$Predictor[i]
            if (grepl(":", pred_val, fixed = TRUE)) {
                ## This is an interaction term - format nicely
                result[i, Predictor := format_interaction_term(pred_val, labels)]
            } else if (pred_val != "-") {
                if (is_binary_affirmative) {
                    ## Binary variable with affirmative level - just show variable name
                    result[i, Predictor := pred_display]
                } else {
                    ## Regular factor level - format as "Variable (Level)"
                    result[i, Predictor := paste0(pred_display, " (", pred_val, ")")]
                }
            }
        }
        
        ## For continuous predictors (marked with "-"), just show variable name
        result[Predictor == "-", Predictor := pred_display]
    }
    
    ## Determine effect type for column naming
    effect_type <- result$effect_type[1]
    if (is.na(effect_type)) effect_type <- "Coefficient"
    
    ## Disable show_events for linear models only
    ## Events make sense for binomial (OR), survival (HR), and count (RR) models
    if (effect_type == "Coefficient") {
        show_events <- FALSE
    }
    
    ## Format effect columns based on columns parameter
    if (columns == "both") {
        ## Effect type label is already correct (Coefficient for linear models)
        effect_display <- effect_type
        
        ## Format unadjusted effects - simple label without "Univariable" prefix
        if ("exp_coef_unadj" %in% names(result)) {
            unadj_label <- paste0(effect_display, " (95% CI)")
            
            if (effect_type %in% c("OR", "HR", "RR")) {
                result[, (unadj_label) := sprintf("%.*f (%.*f-%.*f)", 
                                                  digits, exp_coef_unadj, digits, ci_lower_unadj, digits, ci_upper_unadj)]
            } else {
                result[, (unadj_label) := sprintf("%.*f (%.*f, %.*f)", 
                                                  digits, exp_coef_unadj, digits, ci_lower_unadj, digits, ci_upper_unadj)]
            }
            result[is.na(exp_coef_unadj), (unadj_label) := "-"]
        }
        
        ## Format unadjusted p-values
        if ("p_unadj" %in% names(result)) {
            result[, `Uni p` := format_pvalues_multifit(p_unadj, p_digits)]
        }
        
        ## Format adjusted effects - use aOR/aHR/aRR/Adj. Coefficient without "Multivariable" prefix
        if ("exp_coef_adj" %in% names(result)) {
            adj_effect <- switch(effect_type,
                                 "OR" = "aOR",
                                 "HR" = "aHR", 
                                 "RR" = "aRR",
                                 "Adj. Coefficient"
                                 )
            adj_label <- paste0(adj_effect, " (95% CI)")
            
            if (effect_type %in% c("OR", "HR", "RR")) {
                result[, (adj_label) := sprintf("%.*f (%.*f-%.*f)", 
                                                digits, exp_coef_adj, digits, ci_lower_adj, digits, ci_upper_adj)]
            } else {
                result[, (adj_label) := sprintf("%.*f (%.*f, %.*f)", 
                                                digits, exp_coef_adj, digits, ci_lower_adj, digits, ci_upper_adj)]
            }
            result[is.na(exp_coef_adj), (adj_label) := "-"]
        }
        
        ## Format adjusted p-values
        if ("p_adj" %in% names(result)) {
            result[, `Multi p` := format_pvalues_multifit(p_adj, p_digits)]
        }
        
    } else {
        ## Single column mode (unadjusted or adjusted)
        ## Note: In single-column mode, columns are NOT renamed with _adj/_unadj suffix
        effect_col <- "exp_coef"
        ci_lower_col <- "ci_lower"
        ci_upper_col <- "ci_upper"
        p_col <- "p_value"
        
        ## Effect display name is already correct (Coefficient for linear models)
        effect_display <- effect_type
        
        ## Determine label based on adjustment status - no Univariable/Multivariable prefix
        if (columns == "adjusted") {
            adj_effect <- switch(effect_type,
                                 "OR" = "aOR",
                                 "HR" = "aHR",
                                 "RR" = "aRR",
                                 "Coefficient" = "Adj. Coefficient",
                                 "Adj. Coefficient"
                                 )
            effect_label <- paste0(adj_effect, " (95% CI)")
        } else {
            effect_label <- paste0(effect_display, " (95% CI)")
        }
        
        ## Format effect with CI
        if (effect_col %in% names(result)) {
            if (effect_type %in% c("OR", "HR", "RR")) {
                result[, (effect_label) := sprintf("%.*f (%.*f-%.*f)", 
                                                   digits, get(effect_col), digits, get(ci_lower_col), digits, get(ci_upper_col))]
            } else {
                result[, (effect_label) := sprintf("%.*f (%.*f, %.*f)", 
                                                   digits, get(effect_col), digits, get(ci_lower_col), digits, get(ci_upper_col))]
            }
            result[is.na(get(effect_col)), (effect_label) := "-"]
        }
        
        ## Format p-value
        if (p_col %in% names(result)) {
            result[, `p-value` := format_pvalues_multifit(get(p_col), p_digits)]
        }
    }
    
    ## Format n and events
    if ("n" %in% names(result)) {
        result[, n := data.table::fcase(
                                      is.na(n), NA_character_,
                                      n >= 1000, format(n, big.mark = ","),
                                      default = as.character(n)
                                  )]
    }
    
    if ("events" %in% names(result)) {
        result[, events := data.table::fcase(
                                           is.na(events), NA_character_,
                                           events >= 1000, format(events, big.mark = ","),
                                           default = as.character(events)
                                       )]
        data.table::setnames(result, "events", "Events")
    }
    
    ## Select display columns
    display_cols <- "Outcome"
    if (has_multiple_groups && include_predictor) display_cols <- c(display_cols, "Predictor")
    if (show_n && "n" %in% names(result)) display_cols <- c(display_cols, "n")
    if (show_events && "Events" %in% names(result)) display_cols <- c(display_cols, "Events")
    
    ## Add effect and p-value columns based on layout
    if (columns == "both") {
        effect_display <- effect_type
        unadj_label <- paste0(effect_display, " (95% CI)")
        adj_effect <- switch(effect_type, "OR" = "aOR", "HR" = "aHR", "RR" = "aRR", "Coefficient" = "Adj. Coefficient", "Adj. Coefficient")
        adj_label <- paste0(adj_effect, " (95% CI)")
        
        if (unadj_label %in% names(result)) display_cols <- c(display_cols, unadj_label)
        if ("Uni p" %in% names(result)) display_cols <- c(display_cols, "Uni p")
        if (adj_label %in% names(result)) display_cols <- c(display_cols, adj_label)
        if ("Multi p" %in% names(result)) display_cols <- c(display_cols, "Multi p")
    } else {
        effect_display <- effect_type
        if (columns == "adjusted") {
            adj_effect <- switch(effect_type, "OR" = "aOR", "HR" = "aHR", "RR" = "aRR", "Coefficient" = "Adj. Coefficient", "Adj. Coefficient")
            effect_label <- paste0(adj_effect, " (95% CI)")
        } else {
            effect_label <- paste0(effect_display, " (95% CI)")
        }
        
        if (effect_label %in% names(result)) display_cols <- c(display_cols, effect_label)
        if ("p-value" %in% names(result)) display_cols <- c(display_cols, "p-value")
    }
    
    ## Select only display columns
    formatted <- result[, intersect(display_cols, names(result)), with = FALSE]
    
    return(formatted)
}


#' Format p-values for multifit display
#' 
#' Converts numeric p-values to formatted character strings. Values below the 
#' threshold (determined by digits parameter) display as "< 0.001" (for digits=3),
#' "< 0.0001" (for digits=4), etc. NA values display as "-".
#' 
#' @param p Numeric vector of p-values.
#' @param digits Integer number of decimal places. Also determines the threshold
#'   for "less than" display: threshold = 10^(-digits). Default is 3.
#' @return Character vector of formatted p-values.
#' @keywords internal
format_pvalues_multifit <- function(p, digits = 3) {
    ## Calculate threshold based on digits (e.g., digits=3 -> 0.001, digits=4 -> 0.0001)
    threshold <- 10^(-digits)
    less_than_string <- paste0("< ", format(threshold, scientific = FALSE))
    
    data.table::fcase(
                    is.na(p), "-",
                    p < threshold, less_than_string,
                    default = sprintf(paste0("%.", digits, "f"), p)
                )
}


#' Print method for multifit results
#' @export
print.multifit_result <- function(x, ...) {
    cat("\nMultivariate Analysis Results\n")
    cat("Predictor: ", attr(x, "predictor"), "\n", sep = "")
    cat("Outcomes: ", length(attr(x, "outcomes")), "\n", sep = "")
    cat("Model Type: ", attr(x, "model_type"), "\n", sep = "")
    
    covariates <- attr(x, "covariates")
    if (!is.null(covariates)) {
        cat("Covariates: ", paste(covariates, collapse = ", "), "\n", sep = "")
    }
    
    interactions <- attr(x, "interactions")
    if (!is.null(interactions)) {
        cat("Interactions: ", paste(interactions, collapse = ", "), "\n", sep = "")
    }
    
    random <- attr(x, "random")
    if (!is.null(random)) {
        cat("Random Effects: ", random, "\n", sep = "")
    }
    
    strata <- attr(x, "strata")
    if (!is.null(strata)) {
        cat("Strata: ", strata, "\n", sep = "")
    }
    
    cluster <- attr(x, "cluster")
    if (!is.null(cluster)) {
        cat("Cluster: ", cluster, "\n", sep = "")
    }
    
    cat("Display: ", attr(x, "columns"), "\n", sep = "")
    
    if (!is.null(attr(x, "models"))) {
        cat("Models stored: Yes\n")
    }
    
    cat("\n")
    NextMethod("print", x)
    invisible(x)
}
