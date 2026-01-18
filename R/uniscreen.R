#' Univariable Screening for Multiple Predictors
#'
#' Performs comprehensive univariable (unadjusted) regression analyses by fitting 
#' separate models for each predictor against a single outcome. This function is 
#' designed for initial variable screening, hypothesis generation, and understanding 
#' crude associations before multivariable modeling. Returns publication-ready 
#' formatted results with optional p-value filtering.
#'
#' @param data Data frame or data.table containing the analysis dataset. The 
#'   function automatically converts data frames to data.tables for efficient 
#'   processing.
#'   
#' @param outcome Character string specifying the outcome variable name. For 
#'   survival analysis, use \code{Surv()} syntax from the \pkg{survival} package 
#'   (e.g., \code{"Surv(time, status)"} or \code{"Surv(os_months, os_status)"}).
#'   
#' @param predictors Character vector of predictor variable names to screen. Each 
#'   predictor is tested independently in its own univariable model. Can include 
#'   continuous, categorical (factor), or binary variables.
#'   
#' @param model_type Character string specifying the type of regression model to 
#'   fit. Options include:
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
#'     \item \code{"glmer"} - Generalized linear mixed-effects model for hierarchical 
#'       or clustered data with non-normal outcomes (requires \pkg{lme4} package and 
#'       \code{random} parameter).
#'     \item \code{"lmer"} - Linear mixed-effects model for hierarchical or clustered 
#'       data with continuous outcomes (requires \pkg{lme4} package and \code{random} 
#'       parameter).
#'     \item \code{"coxme"} - Cox mixed-effects model for clustered survival data 
#'       (requires coxme package and \code{random} parameter).
#'   }
#'
#' @param random Character string specifying the random effects formula for 
#'   mixed-effects models (\code{glmer}, \code{lmer}, \code{coxme}). Use standard
#'   lme4/coxme syntax, e.g., \code{"(1|site)"} for random intercepts by site,
#'   \code{"(1|site/patient)"} for nested random effects. Required when 
#'   \code{model_type} is a mixed-effects model type. Default is \code{NULL}.
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
#' @param p_threshold Numeric value between 0 and 1 specifying the p-value threshold 
#'   used to count significant predictors in the printed summary. All predictors
#'   are always included in the output table. Default is 0.05.
#'   
#' @param conf_level Numeric confidence level for confidence intervals. Must be 
#'   between 0 and 1. Default is 0.95 (95\% confidence intervals).
#'   
#' @param reference_rows Logical. If \code{TRUE}, adds rows for reference 
#'   categories of factor variables with baseline values (OR/HR/RR = 1, 
#'   coefficient = 0). Makes tables complete and easier to interpret. 
#'   Default is \code{TRUE}.
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
#'   "< 0.01", etc. Default is 3.
#'   
#' @param labels Named character vector or list providing custom display 
#'   labels for variables. Names should match predictor names, values are the 
#'   display labels. Predictors not in \code{labels} use their original names. 
#'   Default is \code{NULL}.
#'   
#' @param keep_models Logical. If \code{TRUE}, stores all fitted model objects 
#'   in the output as an attribute. This allows access to models for diagnostics, 
#'   predictions, or further analysis, but can consume significant memory for 
#'   large datasets or many predictors. Models are accessible via 
#'   \code{attr(result, "models")}. Default is \code{FALSE}.
#'   
#' @param exponentiate Logical. Whether to exponentiate coefficients (display 
#'   OR/HR/RR instead of log odds/log hazards). Default is \code{NULL}, which 
#'   automatically exponentiates for logistic, Poisson, and Cox models, and 
#'   displays raw coefficients for linear models and other GLM families. Set 
#'   to \code{TRUE} to force exponentiation or \code{FALSE} to force coefficients.
#'
#' @param parallel Logical. If \code{TRUE} (default), fits models in parallel 
#'   using multiple CPU cores for improved performance with many predictors. 
#'   On Unix/Mac systems, uses fork-based parallelism via \code{mclapply}; 
#'   on Windows, uses socket clusters via \code{parLapply}. Set to \code{FALSE} 
#'   for sequential processing.
#'
#' @param n_cores Integer specifying the number of CPU cores to use for parallel 
#'   processing. Default is \code{NULL}, which automatically detects available 
#'   cores and uses \code{detectCores() - 1}. During R CMD check, the number 
#'   of cores is automatically limited to 2 per CRAN policy. Ignored when 
#'   \code{parallel = FALSE}.
#'   
#' @param ... Additional arguments passed to the underlying model fitting functions 
#'   (\code{\link[stats]{glm}}, \code{\link[stats]{lm}}, 
#'   \code{\link[survival]{coxph}}, etc.). Common options include \code{weights}, 
#'   \code{subset}, \code{na.action}, and model-specific control parameters.
#'
#' @return A data.table with S3 class \code{"uniscreen_result"} containing formatted 
#'   univariable screening results. The table structure includes:
#'   \describe{
#'     \item{Variable}{Character. Predictor name or custom label (from \code{labels})}
#'     \item{Group}{Character. For factor variables: category level. For continuous 
#'       variables: typically empty or descriptive statistic label}
#'     \item{n}{Integer. Sample size used in the model (if \code{show_n = TRUE})}
#'     \item{n_group}{Integer. Sample size for this specific factor level 
#'       (factor variables only)}
#'     \item{events}{Integer. Total number of events in the model for survival 
#'       or logistic regression (if \code{show_events = TRUE})}
#'     \item{events_group}{Integer. Number of events for this specific factor 
#'       level (factor variables only)}
#'     \item{OR/HR/RR/Coefficient (95\% CI)}{Character. Formatted effect 
#'       estimate with confidence interval. Column name depends on model type:
#'       "OR (95\% CI)" for logistic, "HR (95\% CI)" for Cox, 
#'       "RR (95\% CI)" for Poisson, "Coefficient (95\% CI)" for linear models}
#'     \item{p-value}{Character. Formatted p-value from the Wald test}
#'   }
#'   
#'   The returned object includes the following attributes accessible via \code{attr()}:
#'   \describe{
#'     \item{raw_data}{data.table. Unformatted numeric results with separate 
#'       columns for coefficients, standard errors, confidence interval bounds, 
#'       etc. Suitable for further statistical analysis or custom formatting}
#'     \item{models}{list (if \code{keep_models = TRUE}). Named list of fitted 
#'       model objects, with predictor names as list names. Access specific models 
#'       via \code{attr(result, "models")[["predictor_name"]]}}
#'     \item{outcome}{Character. The outcome variable name used}
#'     \item{model_type}{Character. The regression model type used}
#'     \item{model_scope}{Character. Always "Univariable" for screening results}
#'     \item{screening_type}{Character. Always "univariable" to identify the 
#'       analysis type}
#'   }
#'
#' @details
#' \strong{Analysis Approach:}
#' 
#' The function implements a comprehensive univariable screening workflow:
#' \enumerate{
#'   \item For each predictor in \code{predictors}, fits a separate model: 
#'     \code{outcome ~ predictor}
#'   \item Extracts coefficients, confidence intervals, and p-values from each model
#'   \item Combines results into a single table for easy comparison
#'   \item Formats output for publication with appropriate effect measures
#' }
#' 
#' Each predictor is tested \emph{independently} - these are crude (unadjusted) 
#' associations that do not account for confounding or interaction effects.
#' 
#' \strong{When to Use Univariable Screening:}
#' \itemize{
#'   \item \strong{Initial variable selection}: Identify predictors associated 
#'     with the outcome before building multivariable models
#'   \item \strong{Hypothesis generation}: Explore potential associations in 
#'     exploratory analyses
#'   \item \strong{Understanding crude associations}: Report unadjusted effects 
#'     alongside adjusted estimates
#'   \item \strong{Variable reduction}: Use p-value thresholds (e.g., p < 0.20) 
#'     to identify candidates for multivariable modeling
#'   \item \strong{Checking multicollinearity}: Compare univariable and 
#'     multivariable effects to identify potential collinearity
#' }
#' 
#' \strong{Factor Variables and Reference Categories:}
#' 
#' When \code{reference_rows = TRUE} (default):
#' \itemize{
#'   \item Reference categories are explicitly shown with OR/HR/RR = 1.00
#'   \item The reference row displays "(Reference)" instead of an effect estimate
#'   \item P-values are shown only for non-reference categories
#'   \item Group-specific sample sizes and event counts are calculated
#' }
#' 
#' \strong{P-value Threshold:}
#' 
#' The \code{p_threshold} parameter controls the significance threshold used 
#' in the printed summary to count how many predictors are significant. All
#' predictors are always included in the output table regardless of this setting.
#' 
#' \strong{Effect Measures by Model Type:}
#' \itemize{
#'   \item \strong{Logistic regression} (\code{model_type = "glm"}, 
#'     \code{family = "binomial"}): Odds ratios (OR)
#'   \item \strong{Cox regression} (\code{model_type = "coxph"}): Hazard ratios (HR)
#'   \item \strong{Poisson regression} (\code{model_type = "glm"}, 
#'     \code{family = "poisson"}): Rate/risk ratios (RR)
#'   \item \strong{Negative binomial} (\code{model_type = "negbin"}): Rate ratios (RR)
#'   \item \strong{Linear regression} (\code{model_type = "lm"} or GLM with 
#'     identity link): Raw coefficient estimates
#'   \item \strong{Gamma regression} (\code{model_type = "glm"}, 
#'     \code{family = "Gamma"}): Multiplicative effects (with default log link)
#' }
#' 
#' \strong{Memory Considerations:}
#' 
#' When \code{keep_models = FALSE} (default), fitted models are discarded after 
#' extracting results to conserve memory. Set \code{keep_models = TRUE} only when 
#' you need:
#' \itemize{
#'   \item Model diagnostic plots
#'   \item Predictions from individual models
#'   \item Additional model statistics not extracted by default
#'   \item Further analysis of specific models
#' }
#'
#' @seealso 
#' \code{\link{fit}} for fitting a single multivariable model,
#' \code{\link{fullfit}} for complete univariable-to-multivariable workflow,
#' \code{\link{compfit}} for comparing multiple models,
#' \code{\link{m2dt}} for converting individual models to tables
#'
#' @examples
#' # Load example data
#' data(clintrial)
#' data(clintrial_labels)
#' 
#' # Example 1: Basic logistic regression screening
#' screen1 <- uniscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "smoking", "hypertension"),
#'     model_type = "glm",
#'     family = "binomial",
#'     parallel = FALSE
#' )
#' print(screen1)
#' 
#' \donttest{
#' # Example 2: With custom variable labels
#' screen2 <- uniscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "treatment"),
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(screen2)
#' 
#' # Example 3: Filter by p-value threshold
#' # Only keep predictors with p < 0.20 (common for screening)
#' screen3 <- uniscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "smoking", "hypertension", 
#'                   "diabetes", "ecog", "stage"),
#'     p_threshold = 0.20,
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(screen3)
#' # Only significant predictors are shown
#' 
#' # Example 4: Cox proportional hazards screening
#' library(survival)
#' cox_screen <- uniscreen(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     predictors = c("age", "sex", "treatment", "stage", "grade"),
#'     model_type = "coxph",
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(cox_screen)
#' # Returns hazard ratios (HR) instead of odds ratios
#' 
#' # Example 5: Keep models for diagnostics
#' screen5 <- uniscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "bmi", "creatinine"),
#'     keep_models = TRUE,
#'     parallel = FALSE
#' )
#' 
#' # Access stored models
#' models <- attr(screen5, "models")
#' summary(models[["age"]])
#' plot(models[["age"]])  # Diagnostic plots
#' 
#' # Example 6: Linear regression screening
#' linear_screen <- uniscreen(
#'     data = clintrial,
#'     outcome = "bmi",
#'     predictors = c("age", "sex", "smoking", "creatinine", "hemoglobin"),
#'     model_type = "lm",
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(linear_screen)
#' 
#' # Example 7: Poisson regression for equidispersed count outcomes
#' # fu_count has variance ≈ mean, appropriate for standard Poisson
#' poisson_screen <- uniscreen(
#'     data = clintrial,
#'     outcome = "fu_count",
#'     predictors = c("age", "stage", "treatment", "surgery"),
#'     model_type = "glm",
#'     family = "poisson",
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(poisson_screen)
#' # Returns rate ratios (RR)
#' 
#' # Example 8: Negative binomial for overdispersed counts
#' # ae_count has variance > mean (overdispersed), use negbin
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'     nb_screen <- uniscreen(
#'         data = clintrial,
#'         outcome = "ae_count",
#'         predictors = c("age", "treatment", "diabetes", "surgery"),
#'         model_type = "negbin",
#'         labels = clintrial_labels,
#'         parallel = FALSE
#'     )
#'     print(nb_screen)
#' }
#' 
#' # Example 9: Gamma regression for positive continuous outcomes (e.g., costs)
#' gamma_screen <- uniscreen(
#'     data = clintrial,
#'     outcome = "los_days",
#'     predictors = c("age", "sex", "treatment", "surgery"),
#'     model_type = "glm",
#'     family = Gamma(link = "log"),
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' print(gamma_screen)
#' 
#' # Example 10: Hide reference rows for factor variables
#' screen10 <- uniscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("treatment", "stage", "grade"),
#'     reference_rows = FALSE,
#'     parallel = FALSE
#' )
#' print(screen10)
#' # Reference categories not shown
#' 
#' # Example 11: Customize decimal places
#' screen11 <- uniscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "bmi", "creatinine"),
#'     digits = 3,      # 3 decimal places for OR
#'     p_digits = 4     # 4 decimal places for p-values
#' )
#' print(screen11)
#' 
#' # Example 12: Hide sample size and event columns
#' screen12 <- uniscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi"),
#'     show_n = FALSE,
#'     show_events = FALSE,
#'     parallel = FALSE
#' )
#' print(screen12)
#' 
#' # Example 13: Access raw numeric data
#' screen13 <- uniscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment"),
#'     parallel = FALSE
#' )
#' raw_data <- attr(screen13, "raw_data")
#' print(raw_data)
#' # Contains unformatted coefficients, SEs, CIs, etc.
#' 
#' # Example 14: Force coefficient display instead of OR
#' screen14 <- uniscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "bmi"),
#'     model_type = "glm",
#'     family = "binomial",
#'     parallel = FALSE,
#'     exponentiate = FALSE  # Show log odds instead of OR
#' )
#' print(screen14)
#' 
#' # Example 15: Screening with weights
#' screen15 <- uniscreen(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     predictors = c("age", "sex", "bmi"),
#'     model_type = "coxph",
#'     weights = runif(nrow(clintrial), min = 0.5, max = 2),  # Random numbers for example
#'     parallel = FALSE
#' )
#' 
#' # Example 16: Strict significance filter (p < 0.05)
#' sig_only <- uniscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "smoking", "hypertension", 
#'                   "diabetes", "ecog", "treatment", "stage", "grade"),
#'     p_threshold = 0.05,
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' 
#' # Check how many predictors passed the filter
#' n_significant <- length(unique(sig_only$Variable[sig_only$Variable != ""]))
#' cat("Significant predictors:", n_significant, "\n")
#' 
#' # Example 17: Complete workflow - screen then use in multivariable
#' # Step 1: Screen with liberal threshold
#' candidates <- uniscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "smoking", "hypertension",
#'                   "diabetes", "treatment", "stage", "grade"),
#'     p_threshold = 0.20,
#'     parallel = FALSE
#' )
#' 
#' # Step 2: Extract significant predictor names from raw data
#' sig_predictors <- unique(attr(candidates, "raw_data")$variable)
#' 
#' # Step 3: Fit multivariable model with selected predictors
#' multi_model <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = sig_predictors,
#'     labels = clintrial_labels
#' )
#' print(multi_model)
#'
#' # Example 18: Mixed-effects logistic regression (glmer)
#' # Accounts for clustering by site
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'     glmer_screen <- uniscreen(
#'         data = clintrial,
#'         outcome = "os_status",
#'         predictors = c("age", "sex", "treatment", "stage"),
#'         model_type = "glmer",
#'         random = "(1|site)",
#'         family = "binomial",
#'         labels = clintrial_labels,
#'         parallel = FALSE
#'     )
#'     print(glmer_screen)
#' }
#'
#' # Example 19: Mixed-effects linear regression (lmer)
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'     lmer_screen <- uniscreen(
#'         data = clintrial,
#'         outcome = "biomarker_x",
#'         predictors = c("age", "sex", "treatment", "smoking"),
#'         model_type = "lmer",
#'         random = "(1|site)",
#'         labels = clintrial_labels,
#'         parallel = FALSE
#'     )
#'     print(lmer_screen)
#' }
#'
#' # Example 20: Mixed-effects Cox model (coxme)
#' if (requireNamespace("coxme", quietly = TRUE)) {
#'     coxme_screen <- uniscreen(
#'         data = clintrial,
#'         outcome = "Surv(os_months, os_status)",
#'         predictors = c("age", "sex", "treatment", "stage"),
#'         model_type = "coxme",
#'         random = "(1|site)",
#'         labels = clintrial_labels,
#'         parallel = FALSE
#'     )
#'     print(coxme_screen)
#' }
#'
#' # Example 21: Mixed-effects with nested random effects
#' # Patients nested within sites
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'     nested_screen <- uniscreen(
#'         data = clintrial,
#'         outcome = "os_status",
#'         predictors = c("age", "treatment"),
#'         model_type = "glmer",
#'         random = "(1|site/patient_id)",
#'         family = "binomial",
#'         parallel = FALSE
#'     )
#' }
#' }
#'
#' @export
uniscreen <- function(data,
                      outcome,
                      predictors,
                      model_type = "glm",
                      family = "binomial",
                      random = NULL,
                      p_threshold = 0.05,
                      conf_level = 0.95,
                      reference_rows = TRUE,
                      show_n = TRUE,
                      show_events = TRUE,
                      digits = 2,
                      p_digits = 3,
                      labels = NULL,
                      keep_models = FALSE,
                      exponentiate = NULL,
                      parallel = TRUE,
                      n_cores = NULL,
                      ...) {
    
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
        p_threshold = p_threshold,
        auto_correct_model = TRUE
    )
    
    ## Apply any auto-corrections
    if (validation$auto_corrected) {
        model_type <- validation$model_type
    }
    
    ## Validate random effects for mixed-effects models
    mixed_types <- c("glmer", "lmer", "coxme")
    if (model_type %in% mixed_types && is.null(random)) {
        stop("'random' parameter is required for mixed-effects models (", model_type, ").\n",
             "Example: random = \"(1|site)\" for random intercepts by site.")
    }
    if (!is.null(random) && !model_type %in% mixed_types) {
        warning("'random' parameter is ignored for non-mixed-effects models (", model_type, ").")
    }
    
    ## Define model fitting function (used by lapply/mclapply)
    fit_one_predictor <- function(pred) {
        ## Build formula - for mixed-effects, append random effects
        formula_str <- paste(outcome, "~", pred)
        if (model_type %in% c("glmer", "lmer")) {
            formula_str <- paste(formula_str, "+", random)
        }
        formula <- stats::as.formula(formula_str)
        
        ## Fit model based on type
        ## Use model=FALSE for glm/lm when not keeping models to save memory
        model <- switch(model_type,
                        "glm" = stats::glm(formula, data = data, family = family, 
                                           model = keep_models, x = FALSE, y = TRUE, ...),
                        "negbin" = {
                            if (!requireNamespace("MASS", quietly = TRUE))
                                stop("Package 'MASS' required for negative binomial models")
                            MASS::glm.nb(formula, data = data, 
                                          model = keep_models, x = FALSE, y = TRUE, ...)
                        },
                        "lm" = stats::lm(formula, data = data, 
                                         model = keep_models, x = FALSE, y = TRUE, ...),
                        "coxph" = {
                            if (!requireNamespace("survival", quietly = TRUE)) 
                                stop("Package 'survival' required for Cox models")
                            survival::coxph(formula, data = data, 
                                            model = keep_models, x = FALSE, y = TRUE, ...)
                        },
                        "clogit" = {
                            if (!requireNamespace("survival", quietly = TRUE))
                                stop("Package 'survival' required for conditional logistic regression")
                            survival::clogit(formula, data = data, ...)
                        },
                        "glmer" = {
                            if (!requireNamespace("lme4", quietly = TRUE))
                                stop("Package 'lme4' required for glmer models")
                            lme4::glmer(formula, data = data, family = family, ...)
                        },
                        "lmer" = {
                            if (!requireNamespace("lme4", quietly = TRUE))
                                stop("Package 'lme4' required for lmer models")
                            lme4::lmer(formula, data = data, ...)
                        },
                        "coxme" = {
                            if (!requireNamespace("coxme", quietly = TRUE))
                                stop("Package 'coxme' required for coxme models")
                            if (!requireNamespace("survival", quietly = TRUE))
                                stop("Package 'survival' required for coxme models")
                            ## coxme uses different formula syntax - random in formula
                            coxme_formula <- stats::as.formula(paste(outcome, "~", pred, "+", random))
                            coxme::coxme(coxme_formula, data = data, ...)
                        },
                        stop("Unsupported model type: ", model_type)
                        )
        
        raw_result <- m2dt(
            data = data,
            model = model,
            conf_level = conf_level,
            keep_qc_stats = FALSE,
            include_intercept = FALSE,
            reference_rows = reference_rows,
            skip_counts = (!show_n && !show_events)
        )
        
        ## Add predictor name for tracking
        raw_result[, predictor := pred]
        
        ## Return both raw result and optionally the model
        list(
            raw = raw_result,
            model = if (keep_models) model else NULL
        )
    }
    
    ## Fit models (parallel or sequential)
    can_use_windows_parallel <- .Platform$OS.type == "windows" && 
        interactive() && 
        !nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_", ""))
    
    use_parallel <- parallel && 
        length(predictors) > 1 && 
        (.Platform$OS.type != "windows" || can_use_windows_parallel)
    
    if (use_parallel) {
        ## Determine number of cores
        if (is.null(n_cores)) {
            n_cores <- max(1L, parallel::detectCores() - 1L)
        }
        ## CRAN policy: respect _R_CHECK_LIMIT_CORES_ environment variable
        chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
        if (nzchar(chk) && chk == "true") {
            n_cores <- min(n_cores, 2L)
        }
        
        if (.Platform$OS.type == "windows") {
            ## Windows: use parLapply() with PSOCK cluster
            ## Only reaches here in interactive sessions with package installed
            cl <- parallel::makeCluster(n_cores)
            on.exit(parallel::stopCluster(cl), add = TRUE)
            
            ## Export the fit_one_predictor closure and all required local objects
            parallel::clusterExport(cl, 
                                    varlist = c("fit_one_predictor", "data", "outcome", 
                                                "model_type", "family", "random", 
                                                "conf_level", "reference_rows", 
                                                "keep_models", "show_n", "show_events"),
                                    envir = environment()
                                    )
            
            ## Export m2dt from summata namespace
            parallel::clusterExport(cl, 
                                    varlist = "m2dt",
                                    envir = asNamespace("summata")
                                    )
            
            ## Load required packages on workers
            parallel::clusterEvalQ(cl, library(data.table))
            
            ## Load survival if needed
            if (model_type %in% c("coxph", "clogit", "coxme")) {
                parallel::clusterEvalQ(cl, library(survival))
            }
            
            ## Load lme4 if needed
            if (model_type %in% c("glmer", "lmer")) {
                parallel::clusterEvalQ(cl, library(lme4))
            }
            
            ## Load coxme if needed
            if (model_type == "coxme") {
                parallel::clusterEvalQ(cl, library(coxme))
            }
            
            results <- parallel::parLapply(cl, predictors, fit_one_predictor)
            
        } else {
            ## Unix/Mac: use mclapply (fork-based, shares memory so no export needed)
            results <- parallel::mclapply(
                                     predictors, 
                                     fit_one_predictor,
                                     mc.cores = n_cores
                                 )
        }
    } else {
        ## Sequential processing with lapply()
        ## Used when: parallel = FALSE, single predictor, or Windows in non-interactive context
        results <- lapply(predictors, fit_one_predictor)
    }
    
    ## Extract raw results and combine
    raw_results <- lapply(results, `[[`, "raw")
    combined_raw <- data.table::rbindlist(raw_results, fill = TRUE)
    
    ## Extract models if kept
    if (keep_models) {
        models <- lapply(results, `[[`, "model")
        names(models) <- predictors
        models <- models[!vapply(models, is.null, logical(1))]
    }
    
    ## Note: p_threshold is used for reporting "Significant" count, not filtering
    
    ## Format results
    formatted <- format_model_table(
        combined_raw,
        show_n = show_n,
        show_events = show_events,
        digits = digits,
        p_digits = p_digits,
        labels = labels,
        exponentiate = exponentiate,
        conf_level = conf_level
    )
    
    ## Attach attributes
    data.table::setattr(formatted, "raw_data", combined_raw)
    
    if (keep_models) {
        data.table::setattr(formatted, "models", models)
    }
    
    data.table::setattr(formatted, "outcome", outcome)
    data.table::setattr(formatted, "model_type", unique(combined_raw$model_type)[1])
    data.table::setattr(formatted, "model_scope", "Univariable")
    data.table::setattr(formatted, "screening_type", "univariable")
    data.table::setattr(formatted, "p_threshold", p_threshold)
    
    class(formatted) <- c("uniscreen_result", class(formatted))
    
    return(formatted)
}


#' Print method for uniscreen results
#' @keywords internal
#' @export
print.uniscreen_result <- function(x, ...) {
    cat("\nUnivariable Screening Results\n")
    cat("Outcome: ", attr(x, "outcome"), "\n", sep = "")
    cat("Model Type: ", attr(x, "model_type"), "\n", sep = "")
    
    n_predictors <- length(unique(x$Variable[x$Variable != ""]))
    cat("Predictors Screened: ", n_predictors, "\n", sep = "")
    
    raw <- attr(x, "raw_data")
    p_thresh <- attr(x, "p_threshold")
    if (is.null(p_thresh)) p_thresh <- 0.05
    
    if (!is.null(raw) && "p_value" %in% names(raw)) {
        sig_predictors <- unique(raw[p_value < p_thresh]$predictor)
        n_sig <- length(sig_predictors)
        cat("Significant (p < ", p_thresh, "): ", n_sig, "\n", sep = "")
    }
    
    if (!is.null(attr(x, "models"))) {
        cat("Models stored: Yes (", length(attr(x, "models")), ")\n", sep = "")
    }
    
    cat("\n")
    NextMethod("print", x)
    invisible(x)
}
