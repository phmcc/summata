#' Fit Regression Model with Publication-Ready Output
#'
#' Provides a unified interface for fitting various types of regression models 
#' with automatic formatting of results for publication. Supports generalized 
#' linear models, linear models, survival models, and mixed effects models with 
#' consistent syntax and output formatting. Handles both univariable and 
#' multivariable models automatically.
#' 
#' Can be used in two ways:
#' \enumerate{
#'   \item \strong{Formula-based}: Provide \code{data}, \code{outcome}, and 
#'     \code{predictors} to fit a new model
#'   \item \strong{Model-based}: Provide a pre-fitted \code{model} object to 
#'     format existing results (e.g., from \code{MASS::glm.nb()})
#' }
#'
#' @param data Data frame or data.table containing the analysis dataset. 
#'   Required for formula-based workflow; optional for model-based workflow
#'   (extracted from model if not provided).
#'   
#' @param outcome Character string specifying the outcome variable name. For 
#'   survival analysis, use \code{Surv()} syntax from the \pkg{survival} package 
#'   (e.g., \code{"Surv(time, status)"} or \code{"Surv(os_months, os_status)"}).
#'   Required for formula-based workflow; ignored if \code{model} is provided.
#'   
#' @param predictors Character vector of predictor variable names to include in 
#'   the model. All predictors are included simultaneously (multivariable model). 
#'   For univariable models, provide a single predictor. Can include continuous, 
#'   categorical (factor), or binary variables. Required for formula-based
#'   workflow; ignored if \code{model} is provided.
#'
#' @param model Optional pre-fitted model object to format. When provided,
#'   \code{outcome} and \code{predictors} are ignored and the model is formatted
#'   directly. Supported model classes include:
#'   \itemize{
#'     \item \code{glm} - Generalized linear models
#'     \item \code{negbin} - Negative binomial models from \code{MASS::glm.nb()}
#'     \item \code{lm} - Linear models
#'     \item \code{coxph} - Cox proportional hazards models
#'     \item \code{lmerMod}, \code{glmerMod} - Mixed-effects models from lme4
#'     \item \code{coxme} - Mixed-effects Cox models
#'   }
#'   
#' @param model_type Character string specifying the type of regression model. 
#'   Ignored if \code{model} is provided. Options include:
#'   \itemize{
#'     \item \code{"glm"} - Generalized linear model (default). Supports multiple 
#'       distributions via the \code{family} parameter including logistic, Poisson, 
#'       Gamma, Gaussian, and quasi-likelihood models.
#'     \item \code{"lm"} - Linear regression for continuous outcomes with normally 
#'       distributed errors. Equivalent to \code{glm} with \code{family = "gaussian"}.
#'     \item \code{"coxph"} - Cox proportional hazards model for time-to-event 
#'       survival analysis. Requires \code{Surv()} outcome syntax.
#'     \item \code{"clogit"} - Conditional logistic regression for matched 
#'       case-control studies or stratified analyses.
#'     \item \code{"negbin"} - Negative binomial regression for overdispersed count 
#'       data (requires MASS package). Estimates an additional dispersion parameter 
#'       compared to Poisson regression.
#'     \item \code{"lmer"} - Linear mixed-effects model for hierarchical or clustered 
#'       data with continuous outcomes (requires \pkg{lme4} package).
#'     \item \code{"glmer"} - Generalized linear mixed-effects model for hierarchical 
#'       or clustered data with non-normal outcomes (requires \pkg{lme4} package).
#'     \item \code{"coxme"} - Cox mixed-effects model for clustered survival data 
#'       (requires coxme package).
#'   }
#'   
#' @param family For GLM and GLMER models, specifies the error distribution and link 
#'   function. Can be a character string, a family function, or a family object.
#'   Ignored for non-GLM/GLMER models.
#'   
#'   \strong{Binary/Binomial outcomes:}
#'   \itemize{
#'     \item \code{"binomial"} or \code{binomial()} - Logistic regression for binary 
#'       outcomes (0/1, TRUE/FALSE). Returns odds ratios (OR). Default.
#'     \item \code{"quasibinomial"} or \code{quasibinomial()} - Logistic regression 
#'       with overdispersion. Use when residual deviance >> residual df.
#'     \item \code{binomial(link = "probit")} - Probit regression (normal CDF link).
#'     \item \code{binomial(link = "cloglog")} - Complementary log-log link for 
#'       asymmetric binary outcomes.
#'   }
#'   
#'   \strong{Count outcomes:}
#'   \itemize{
#'     \item \code{"poisson"} or \code{poisson()} - Poisson regression for count 
#'       data. Returns rate ratios (RR). Assumes mean = variance.
#'     \item \code{"quasipoisson"} or \code{quasipoisson()} - Poisson regression 
#'       with overdispersion. Use when variance > mean.
#'   }
#'   
#'   \strong{Continuous outcomes:}
#'   \itemize{
#'     \item \code{"gaussian"} or \code{gaussian()} - Normal/Gaussian distribution 
#'       for continuous outcomes. Equivalent to linear regression.
#'     \item \code{gaussian(link = "log")} - Log-linear model for positive continuous 
#'       outcomes. Returns multiplicative effects.
#'     \item \code{gaussian(link = "inverse")} - Inverse link for specific applications.
#'   }
#'   
#'   \strong{Positive continuous outcomes:}
#'   \itemize{
#'     \item \code{"Gamma"} or \code{Gamma()} - Gamma distribution for positive, 
#'       right-skewed continuous data (e.g., costs, lengths of stay). Default log link.
#'     \item \code{Gamma(link = "inverse")} - Gamma with inverse (canonical) link.
#'     \item \code{Gamma(link = "identity")} - Gamma with identity link for additive 
#'       effects on positive outcomes.
#'     \item \code{"inverse.gaussian"} or \code{inverse.gaussian()} - Inverse Gaussian 
#'       for positive, highly right-skewed data.
#'   }
#'   
#'   For negative binomial regression (overdispersed counts), use 
#'   \code{model_type = "negbin"} instead of the \code{family} parameter.
#'   
#'   See \code{\link[stats]{family}} for additional details and options.
#'   
#' @param interactions Character vector of interaction terms using colon 
#'   notation (e.g., \code{c("age:sex", "treatment:stage")}). Interaction terms 
#'   are added to the model in addition to main effects. Default is \code{NULL} 
#'   (no interactions).
#'   
#' @param strata For Cox or conditional logistic models, character string naming 
#'   the stratification variable. Creates separate baseline hazards for each 
#'   stratum level without estimating stratum effects. Default is \code{NULL}.
#'   
#' @param cluster For Cox models, character string naming the variable for 
#'   robust clustered standard errors. Accounts for within-cluster correlation 
#'   (e.g., patients within hospitals). Default is \code{NULL}.
#'   
#' @param weights Character string naming the weights variable in \code{data}. 
#'   The specified column should contain numeric values for observation weights. 
#'   Used for weighted regression, survey data, or inverse probability weighting. 
#'   Default is \code{NULL}.
#'   
#' @param conf_level Numeric confidence level for confidence intervals. Must be 
#'   between 0 and 1. Default is 0.95 (95\% confidence intervals).
#'   
#' @param reference_rows Logical. If \code{TRUE}, adds rows for reference 
#'   categories of factor variables with baseline values (OR/HR/RR = 1, 
#'   coefficient = 0). Default is \code{TRUE}.
#'   
#' @param show_n Logical. If \code{TRUE}, includes the sample size column in 
#'   the output. Default is \code{TRUE}.
#'   
#' @param show_events Logical. If \code{TRUE}, includes the events column in 
#'   the output (for survival and logistic regression). Default is \code{TRUE}.
#'   
#' @param digits Integer specifying the number of decimal places for effect 
#'   estimates (OR, HR, RR, coefficients). Default is 2.
#'   
#' @param p_digits Integer specifying the number of decimal places for p-values. 
#'   P-values smaller than \code{10^(-p_digits)} are displayed as "< 0.001", 
#'   etc. Default is 3.
#'   
#' @param labels Named character vector or list providing custom display 
#'   labels for variables. Names should match variable names, values are display 
#'   labels. Default is \code{NULL}.
#'   
#' @param keep_qc_stats Logical. If \code{TRUE}, includes model quality statistics 
#'   (AIC, BIC, R-squared, concordance, etc.) in the raw data attribute for 
#'   model diagnostics and comparison. Default is \code{TRUE}.
#'   
#' @param exponentiate Logical. Whether to exponentiate coefficients. Default 
#'   is \code{NULL}, which automatically exponentiates for logistic, Poisson, 
#'   and Cox models, and displays raw coefficients for linear models. Set to 
#'   \code{TRUE} to force exponentiation or \code{FALSE} to force coefficients.
#'   
#' @param ... Additional arguments passed to the underlying model fitting 
#'   function (\code{\link[stats]{glm}}, \code{\link[stats]{lm}}, 
#'   \code{\link[survival]{coxph}}, \code{\link[lme4]{glmer}}, etc.). Common 
#'   options include \code{subset}, \code{na.action}, and model-specific control 
#'   parameters.
#'
#' @return A data.table with S3 class \code{"fit_result"} containing formatted 
#'   regression results. The table structure includes:
#'   \describe{
#'     \item{Variable}{Character. Predictor name or custom label}
#'     \item{Group}{Character. For factor variables: category level. For 
#'       interactions: interaction term. For continuous: typically empty}
#'     \item{n}{Integer. Total sample size (if \code{show_n = TRUE})}
#'     \item{n_group}{Integer. Sample size for this factor level}
#'     \item{events}{Integer. Total number of events (if \code{show_events = TRUE})}
#'     \item{events_group}{Integer. Events for this factor level}
#'     \item{OR/HR/RR/Coefficient or aOR/aHR/aRR/Adj. Coefficient (95\% CI)}{Character. 
#'       Formatted effect estimate with confidence interval. Column name depends on 
#'       model type and scope. Univariable models use: OR, HR, RR, Coefficient.
#'       Multivariable models use adjusted notation: aOR, aHR, aRR, Adj. Coefficient}
#'     \item{p-value}{Character. Formatted p-value from Wald test}
#'   }
#'   
#'   The returned object includes the following attributes accessible via \code{attr()}:
#'   \describe{
#'     \item{model}{The fitted model object (glm, lm, coxph, etc.). Access for 
#'       diagnostics, predictions, or further analysis}
#'     \item{raw_data}{data.table. Unformatted numeric results with columns for 
#'       coefficients, standard errors, confidence bounds, quality statistics, etc.}
#'     \item{outcome}{Character. The outcome variable name}
#'     \item{predictors}{Character vector. The predictor variable names}
#'     \item{formula_str}{Character. The complete model formula as a string}
#'     \item{model_scope}{Character. "Univariable" (one predictor) or 
#'       "Multivariable" (multiple predictors)}
#'     \item{model_type}{Character. The regression model type used}
#'     \item{interactions}{Character vector (if interactions specified). 
#'       The interaction terms included}
#'     \item{strata}{Character (if stratification used). The stratification variable}
#'     \item{cluster}{Character (if clustering used). The cluster variable}
#'     \item{weights}{Character (if weighting used). The weights variable}
#'   }
#'
#' @details
#' \strong{Model Scope Detection:}
#' 
#' The function automatically detects whether the model is:
#' \itemize{
#'   \item \strong{Univariable}: Single predictor (e.g., \code{predictors = "age"})
#'     - Effect estimates labeled as "Univariable OR", "Univariable HR", etc.
#'     - Represents crude (unadjusted) association
#'   \item \strong{Multivariable}: Multiple predictors (e.g., 
#'     \code{predictors = c("age", "sex", "treatment")})
#'     - Effect estimates labeled as "Multivariable aOR", "Multivariable aHR", etc.
#'     - "a" prefix indicates "adjusted" for other variables in the model
#'     - Represents associations adjusted for confounding
#' }
#' 
#' \strong{Interaction Terms:}
#' 
#' Interactions are specified using colon notation and added to the model:
#' \itemize{
#'   \item \code{interactions = c("age:treatment")} creates interaction 
#'     between age and treatment
#'   \item Main effects for both variables are automatically included
#'   \item Multiple interactions can be specified: 
#'     \code{c("age:sex", "treatment:stage")}
#'   \item For interactions between categorical variables, separate terms are 
#'     created for each combination of levels
#' }
#' 
#' \strong{Stratification (Cox/Conditional Logistic):}
#' 
#' The \code{strata} parameter creates separate baseline hazards:
#' \itemize{
#'   \item Allows baseline hazard to vary across strata without estimating 
#'     stratum effects
#'   \item Useful when proportional hazards assumption violated across strata
#'   \item Example: \code{strata = "center"} for multicenter studies
#'   \item Stratification variable is not included as a predictor
#' }
#' 
#' \strong{Clustering (Cox Models):}
#' 
#' The \code{cluster} parameter computes robust standard errors:
#' \itemize{
#'   \item Accounts for within-cluster correlation (e.g., multiple observations 
#'     per patient)
#'   \item Uses sandwich variance estimator
#'   \item Does not change point estimates, only standard errors and p-values
#' }
#' 
#' \strong{Weighting:}
#' 
#' The \code{weights} parameter enables weighted regression:
#' \itemize{
#'   \item For survey data with sampling weights
#'   \item Inverse probability weighting for causal inference
#'   \item Frequency weights for aggregated data
#'   \item Weights should be in a column of \code{data}
#' }
#' 
#' \strong{Mixed Effects Models (lmer/glmer/coxme):}
#' 
#' Mixed effects models handle hierarchical or clustered data:
#' \itemize{
#'   \item Use \code{model_type = "lmer"} for continuous/normal outcomes
#'   \item Use \code{model_type = "glmer"} with appropriate \code{family} for GLM outcomes
#'   \item Use \code{model_type = "coxme"} for survival outcomes with clustering
#'   \item Random effects are specified in predictors using lme4 syntax:
#'     \itemize{
#'       \item \code{"(1|site)"} - Random intercepts by site
#'       \item \code{"(treatment|site)"} - Random slopes for treatment by site
#'       \item \code{"(1 + treatment|site)"} - Both random intercepts and slopes
#'     }
#'   \item Include random effects as part of the predictors vector
#'   \item Example: \code{predictors = c("age", "treatment", "(1|site)")}
#' }
#' 
#' \strong{Effect Measures by Model Type:}
#' \itemize{
#'   \item \strong{Logistic} (\code{family = "binomial"/"quasibinomial"}): Odds ratios (OR/aOR)
#'   \item \strong{Cox} (\code{model_type = "coxph"}): Hazard ratios (HR/aHR)
#'   \item \strong{Poisson/Count} (\code{family = "poisson"/"quasipoisson"}): Rate ratios (RR/aRR)
#'   \item \strong{Negative binomial} (\code{model_type = "negbin"}): Rate ratios (RR/aRR)
#'   \item \strong{Gamma/Log-link}: Ratios (multiplicative effects)
#'   \item \strong{Linear/Gaussian}: Raw coefficient estimates (additive effects)
#' }
#' 
#' \strong{Confidence Intervals:}
#' 
#' Confidence intervals are computed using the Wald method:
#' \deqn{\hat{\beta} \pm z_{1-\alpha/2} \times \mathrm{SE}(\hat{\beta})}
#' 
#' This is the standard approach for GLM, Cox, and mixed-effects regression, 
#' and is appropriate for standard sample sizes. The Wald interval is computed
#' directly from the coefficient and standard error, avoiding redundant matrix
#' operations.
#' 
#' For small samples, sparse data, or parameters near boundary values, profile 
#' likelihood confidence intervals may be preferred. These can be obtained from 
#' the underlying model object:
#' 
#' \preformatted{
#' result <- fit(data, outcome, predictors)
#' confint(attr(result, "model"))
#' }
#'
#' @seealso 
#' \code{\link{uniscreen}} for univariable screening of multiple predictors,
#' \code{\link{fullfit}} for complete univariable-to-multivariable workflow,
#' \code{\link{compfit}} for comparing multiple models,
#' \code{\link{m2dt}} for model-to-table conversion
#'
#' @examples
#' # Load example data
#' data(clintrial)
#' data(clintrial_labels)
#' library(survival)
#' 
#' # Example 1: Univariable logistic regression
#' uni_model <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = "age"
#' )
#' print(uni_model)
#' # Labeled as "Univariable OR"
#' 
#' \donttest{
#' # Example 2: Multivariable logistic regression
#' multi_model <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "treatment"),
#'     labels = clintrial_labels
#' )
#' print(multi_model)
#' # Labeled as "Multivariable aOR" (adjusted OR)
#' 
#' # Example 3: Cox proportional hazards model
#' cox_model <- fit(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     predictors = c("age", "sex", "treatment", "stage"),
#'     model_type = "coxph",
#'     labels = clintrial_labels
#' )
#' print(cox_model)
#' 
#' # Example 4: Model with interaction terms
#' interact_model <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "treatment", "sex"),
#'     interactions = c("age:treatment"),
#'     labels = clintrial_labels
#' )
#' print(interact_model)
#' 
#' # Example 5: Cox model with stratification
#' strat_model <- fit(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     predictors = c("age", "sex", "treatment"),
#'     model_type = "coxph",
#'     strata = "site",  # Separate baseline hazards by site
#'     labels = clintrial_labels
#' )
#' print(strat_model)
#' 
#' # Example 6: Cox model with clustering
#' cluster_model <- fit(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     predictors = c("age", "treatment"),
#'     model_type = "coxph",
#'     cluster = "site",  # Robust SEs accounting for site clustering
#'     labels = clintrial_labels
#' )
#' print(cluster_model)
#' 
#' # Example 7: Linear regression
#' linear_model <- fit(
#'     data = clintrial,
#'     outcome = "bmi",
#'     predictors = c("age", "sex", "smoking"),
#'     model_type = "lm",
#'     labels = clintrial_labels
#' )
#' print(linear_model)
#' 
#' # Example 8: Poisson regression for count data
#' poisson_model <- fit(
#'     data = clintrial,
#'     outcome = "los_days",
#'     predictors = c("age", "treatment", "surgery", "stage"),
#'     model_type = "glm",
#'     family = "poisson",
#'     labels = clintrial_labels
#' )
#' print(poisson_model)
#' # Returns rate ratios (RR/aRR)
#' 
#' # Example 9: Negative binomial regression for overdispersed counts
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'     nb_result <- fit(
#'         data = clintrial,
#'         outcome = "los_days",
#'         predictors = c("age", "sex", "treatment", "stage"),
#'         model_type = "negbin",
#'         labels = clintrial_labels
#'     )
#'     print(nb_result)
#' }
#' 
#' # Example 10: Gamma regression for positive continuous outcomes
#' gamma_model <- fit(
#'     data = clintrial,
#'     outcome = "los_days",
#'     predictors = c("age", "treatment", "surgery"),
#'     model_type = "glm",
#'     family = Gamma(link = "log"),
#'     labels = clintrial_labels
#' )
#' print(gamma_model)
#' 
#' # Example 11: Access the underlying fitted model
#' result <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi")
#' )
#' 
#' # Get the model object
#' model_obj <- attr(result, "model")
#' summary(model_obj)
#' 
#' # Model diagnostics
#' plot(model_obj)
#' 
#' # Predictions
#' preds <- predict(model_obj, type = "response")
#' 
#' # Example 12: Access raw numeric data
#' raw_data <- attr(result, "raw_data")
#' print(raw_data)
#' # Contains unformatted coefficients, SEs, CIs, AIC, BIC, etc.
#' 
#' # Example 13: Multiple interactions
#' complex_model <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment", "bmi"),
#'     interactions = c("age:treatment", "sex:bmi"),
#'     labels = clintrial_labels
#' )
#' print(complex_model)
#' 
#' # Example 14: Customize output columns
#' minimal <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment"),
#'     show_n = FALSE,
#'     show_events = FALSE,
#'     reference_rows = FALSE
#' )
#' print(minimal)
#' 
#' # Example 15: Different confidence levels
#' ci90 <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "treatment"),
#'     conf_level = 0.90  # 90% confidence intervals
#' )
#' print(ci90)
#' 
#' # Example 16: Force coefficient display instead of OR
#' coef_model <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "bmi"),
#'     exponentiate = FALSE  # Show log odds instead of OR
#' )
#' print(coef_model)
#' 
#' # Example 17: Check model quality statistics
#' result <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment", "stage"),
#'     keep_qc_stats = TRUE
#' )
#' 
#' raw <- attr(result, "raw_data")
#' cat("AIC:", raw$AIC[1], "\n")
#' cat("BIC:", raw$BIC[1], "\n")
#' cat("C-statistic:", raw$c_statistic[1], "\n")
#' }
#' 
#' # Example 18: Linear mixed effects model with site random effects
#' # (Requires \pkg{lme4} package)
#' \dontrun{
#' library(lme4)
#' lmer_model <- fit(
#'     data = clintrial,
#'     outcome = "los_days",
#'     predictors = c("age", "treatment", "stage", "(1|site)"),
#'     model_type = "lmer"
#' )
#' print(lmer_model)
#' 
#' # With random slopes for treatment effect by site
#' lmer_slopes <- fit(
#'     data = clintrial,
#'     outcome = "los_days", 
#'     predictors = c("age", "treatment", "stage", "(treatment|site)"),
#'     model_type = "lmer"
#' )
#' print(lmer_slopes)
#'
#' # Example 19: Format a pre-fitted model (model-based workflow)
#' # Useful for models fitted outside of fit(), like MASS::glm.nb()
#' pre_fitted <- glm(os_status ~ age + sex + treatment,
#'                   family = binomial, data = clintrial)
#' result <- fit(model = pre_fitted, labels = clintrial_labels)
#' print(result)
#' }
#'
#' @export
fit <- function(data = NULL, 
                outcome = NULL,
                predictors = NULL,
                model = NULL,
                model_type = "glm",
                family = "binomial",
                interactions = NULL,
                strata = NULL,
                cluster = NULL,
                weights = NULL,
                conf_level = 0.95,
                reference_rows = TRUE,
                show_n = TRUE,
                show_events = TRUE,
                digits = 2,
                p_digits = 3,
                labels = NULL,
                keep_qc_stats = TRUE,
                exponentiate = NULL,
                ...) {
    
    ## Check if we're using model-based workflow
    if (!is.null(model)) {
        ## Model-based workflow: format existing model
        
        ## Get data from model if not provided
        if (is.null(data)) {
            data <- get_model_data(model)
            if (is.null(data)) {
                stop("Cannot extract data from model. Please provide 'data' parameter.")
            }
        }
        
        ## Convert to data.table
        if (!data.table::is.data.table(data)) {
            data <- data.table::as.data.table(data)
        }
        
        ## Store data in model for m2dt
        if (isS4(model)) {
            attr(model, "data") <- data
        } else {
            model$data <- data
            attr(model, "data") <- data
        }
        
        ## Convert to readable format using m2dt
        raw_data <- m2dt(data = data,
                         model = model, 
                         conf_level = conf_level,
                         keep_qc_stats = keep_qc_stats,
                         include_intercept = FALSE,
                         reference_rows = reference_rows,
                         skip_counts = (!show_n && !show_events))
        
        ## Format results for publication
        formatted <- format_model_table(raw_data,
                                        show_n = show_n,
                                        show_events = show_events,
                                        digits = digits,
                                        p_digits = p_digits,
                                        labels = labels,
                                        exponentiate = exponentiate)
        
        ## Extract formula and predictors from model
        formula_str <- deparse(stats::formula(model))
        
        ## Attach metadata
        data.table::setattr(formatted, "model", model)
        data.table::setattr(formatted, "raw_data", raw_data)
        data.table::setattr(formatted, "formula_str", formula_str)
        data.table::setattr(formatted, "model_scope", raw_data$model_scope[1])
        data.table::setattr(formatted, "model_type", raw_data$model_type[1])
        
        class(formatted) <- c("fit_result", class(formatted))
        
        return(formatted)
    }
    
    ## Formula-based workflow: fit new model
    
    ## Validate required parameters for formula-based workflow
    if (is.null(data)) {
        stop("'data' parameter is required when not providing a pre-fitted model.")
    }
    if (is.null(outcome)) {
        stop("'outcome' parameter is required when not providing a pre-fitted model.")
    }
    if (is.null(predictors)) {
        stop("'predictors' parameter is required when not providing a pre-fitted model.")
    }
    
    ## Store original data name for reference
    data_name <- deparse(substitute(data))
    
    ## Convert to data.table once at the start
    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }
    
    ## Validate inputs and auto-correct model type if needed
    validation <- validate_fit_inputs(
        data = data,
        outcome = outcome,
        predictors = predictors,
        model_type = model_type,
        family = if (model_type %in% c("glm", "glmer")) family else NULL,
        conf_level = conf_level,
        digits = digits,
        p_digits = p_digits,
        auto_correct_model = TRUE
    )
    
    ## Apply any auto-corrections
    if (validation$auto_corrected) {
        model_type <- validation$model_type
    }
    
    ## Build formula string
    if (!is.null(interactions)) {
        formula_str <- paste(outcome, "~", 
                             paste(c(predictors, interactions), collapse = " + "))
    } else {
        formula_str <- paste(outcome, "~", paste(predictors, collapse = " + "))
    }
    
    ## Add strata if provided (for Cox models)
    if (!is.null(strata) && model_type %in% c("coxph", "clogit")) {
        formula_str <- paste(formula_str, "+ strata(", strata, ")")
    }
    
    formula <- stats::as.formula(formula_str)
    
    ## Fit model based on type - use 'data' directly
    if (model_type == "glm") {
        if (!is.null(weights)) {
            model <- stats::glm(formula, data = data, family = family, 
                                weights = data[[weights]], ...)
        } else {
            model <- stats::glm(formula, data = data, family = family, ...)
        }
        
    } else if (model_type == "negbin") {
        if (!requireNamespace("MASS", quietly = TRUE))
            stop("Package 'MASS' required for negative binomial models")
        if (!is.null(weights)) {
            model <- MASS::glm.nb(formula, data = data, 
                                   weights = data[[weights]], ...)
        } else {
            model <- MASS::glm.nb(formula, data = data, ...)
        }
        
    } else if (model_type == "lm") {
        if (!is.null(weights)) {
            model <- stats::lm(formula, data = data, weights = data[[weights]], ...)
        } else {
            model <- stats::lm(formula, data = data, ...)
        }
        
    } else if (model_type == "coxph") {
        if (!requireNamespace("survival", quietly = TRUE)) 
            stop("Package 'survival' required for Cox models")
        
        if (!is.null(cluster)) {
            model <- survival::coxph(formula, data = data, 
                                     cluster = data[[cluster]], ...)
        } else {
            model <- survival::coxph(formula, data = data, ...)
        }
        
    } else if (model_type == "clogit") {
        if (!requireNamespace("survival", quietly = TRUE)) 
            stop("Package 'survival' required for conditional logistic regression")
        model <- survival::clogit(formula, data = data, ...)
        
    } else if (model_type == "coxme") {
        if (!requireNamespace("coxme", quietly = TRUE))
            stop("Package 'coxme' required for mixed effects Cox models")
        model <- coxme::coxme(formula, data = data, ...)
        
    } else if (model_type == "lmer") {
        if (!requireNamespace("lme4", quietly = TRUE))
            stop("Package 'lme4' required for linear mixed effects models")
        model <- lme4::lmer(formula, data = data, ...)
        
    } else if (model_type == "glmer") {
        if (!requireNamespace("lme4", quietly = TRUE))
            stop("Package 'lme4' required for mixed effects models")
        model <- lme4::glmer(formula, data = data, family = family, ...)
        
    } else {
        stop("Unsupported model type: ", model_type)
    }
    
    ## Store the data directly in the model
    if (isS4(model)) {
        attr(model, "data") <- data
    } else {
        model$data <- data
        attr(model, "data") <- data
    }
    
    ## Attach metadata as attributes
    data.table::setattr(model, "formula_str", formula_str)
    data.table::setattr(model, "predictors", predictors)
    data.table::setattr(model, "model_type", model_type)
    data.table::setattr(model, "data_name", data_name)
    
    if (!is.null(interactions)) {
        data.table::setattr(model, "interactions", interactions)
    }
    if (!is.null(strata)) {
        data.table::setattr(model, "strata", strata)
    }
    if (!is.null(cluster)) {
        data.table::setattr(model, "cluster", cluster)
    }
    if (!is.null(weights)) {
        data.table::setattr(model, "weights", weights)
    }
    
    ## Convert to readable format using m2dt
    raw_data <- m2dt(data = data,
                     model = model, 
                     conf_level = conf_level,
                     keep_qc_stats = keep_qc_stats,
                     include_intercept = FALSE,
                     reference_rows = reference_rows,
                     skip_counts = (!show_n && !show_events))

    ## Format results for publication
    formatted <- format_model_table(raw_data,
                                    show_n = show_n,
                                    show_events = show_events,
                                    digits = digits,
                                    p_digits = p_digits,
                                    labels = labels,
                                    exponentiate = exponentiate)

    ## Attach metadata
    data.table::setattr(formatted, "model", model)
    data.table::setattr(formatted, "raw_data", raw_data)
    data.table::setattr(formatted, "outcome", outcome)
    data.table::setattr(formatted, "predictors", predictors)
    data.table::setattr(formatted, "formula_str", formula_str)
    data.table::setattr(formatted, "model_scope", raw_data$model_scope[1])
    data.table::setattr(formatted, "model_type", raw_data$model_type[1])

    ## Add metadata from model fitting
    if (!is.null(interactions)) {
        data.table::setattr(formatted, "interactions", interactions)
    }
    if (!is.null(strata)) {
        data.table::setattr(formatted, "strata", strata)
    }
    if (!is.null(cluster)) {
        data.table::setattr(formatted, "cluster", cluster)
    }
    if (!is.null(weights)) {
        data.table::setattr(formatted, "weights", weights)
    }

    class(formatted) <- c("fit_result", class(formatted))

    return(formatted)
}

#' Print method for fit results
#' 
#' Displays a summary header with model scope (Univariable/Multivariable),
#' model type, formula, sample size, and event count before printing the
#' formatted results table.
#' 
#' @param x fit_result object.
#' @param ... Additional arguments passed to print methods.
#' @return Invisibly returns the input object.
#' @keywords internal
#' @export
print.fit_result <- function(x, ...) {
    cat("\n", attr(x, "model_scope"), " ", attr(x, "model_type"), " Model\n", sep = "")
    cat("Formula: ", attr(x, "formula_str"), "\n", sep = "")
    
    ## Get sample size from raw results
    raw <- attr(x, "raw_data")
    if (!is.null(raw)) {
        cat("N = ", raw$n[1], sep = "")
        if (!is.na(raw$events[1])) cat(", Events = ", raw$events[1], sep = "")
        cat("\n")
    }
    cat("\n")
    
    NextMethod("print", x)
    invisible(x)
}
