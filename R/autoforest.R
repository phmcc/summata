#' Create Forest Plot with Automatic Model Detection
#'
#' A convenience wrapper that automatically detects the input type and routes to
#' the appropriate specialized forest plot function. This eliminates the need to
#' remember which forest function to call for different model types or analysis
#' objects, making it ideal for exploratory analysis and rapid prototyping.
#'
#' @param x One of the following:
#'   \itemize{
#'     \item A fitted model object: \code{glm}, \code{lm}, \code{coxph}, or \code{clogit}
#'     \item A \code{fit_result} object from \code{\link{fit}}
#'     \item A \code{fullfit_result} object from \code{\link{fullfit}}
#'     \item A \code{uniscreen_result} object from \code{\link{uniscreen}}
#'     \item A \code{multifit_result} object from \code{\link{multifit}}
#'   }
#'
#' @param data Data frame or data.table containing the original data. Required 
#'   when \code{x} is a raw model object. Ignored when \code{x} is a result object
#'   from \code{fit()}, \code{fullfit()}, \code{uniscreen()}, or \code{multifit()}
#'   since these contain embedded data.
#'
#' @param title Character string for plot title. If \code{NULL} (default), an 
#'   appropriate title is generated based on the detected model type:
#'   \itemize{
#'     \item Cox models: "Cox Proportional Hazards Model"
#'     \item Logistic regression: "Logistic Regression Model"
#'     \item Poisson regression: "Poisson Regression Model"
#'     \item Linear regression: "Linear Regression Model"
#'     \item Uniscreen results: "Univariable [Type] Screening"
#'     \item Multifit results: "Multi-Outcome [Type] Analysis"
#'   }
#'
#' @param ... Additional arguments passed to the specific forest plot function.
#'   Common arguments include:
#'   \describe{
#'     \item{labels}{Named character vector for variable labels}
#'     \item{digits}{Number of decimal places for estimates (default 2)}
#'     \item{p_digits}{Number of decimal places for p-values (default 3)}
#'     \item{conf_level}{Confidence level for intervals (default 0.95)}
#'     \item{show_n}{Logical, show sample sizes (default TRUE)}
#'     \item{show_events}{Logical, show event counts (default TRUE for survival/binomial)}
#'     \item{qc_footnote}{Logical, show model QC stats in footer (default TRUE)}
#'     \item{zebra_stripes}{Logical, alternating row shading (default TRUE)}
#'     \item{indent_groups}{Logical, indent factor levels (default FALSE)}
#'     \item{color}{Color for point estimates}
#'     \item{table_width}{Proportion of width for table (default 0.6)}
#'     \item{plot_width, plot_height}{Explicit dimensions}
#'     \item{units}{Dimension units: "in", "cm", or "mm"}
#'   }
#'   See the documentation for the specific forest function for all available options.
#'
#' @return A \code{ggplot} object containing the forest plot. The returned object
#'   includes an attribute \code{"recommended_dims"} accessible via 
#'   \code{attr(plot, "recommended_dims")} containing recommended width and height.
#'
#' @details
#' \strong{Automatic Detection Logic:}
#' 
#' The function uses the following priority order for detection:
#' \enumerate{
#'   \item \strong{uniscreen results}: Detected by class \code{"uniscreen_result"} or
#'     presence of attributes \code{outcome}, \code{predictors}, \code{model_type},
#'     and \code{model_scope = "Univariable"}. Routes to \code{\link{uniforest}}.
#'   \item \strong{multifit results}: Detected by presence of attributes 
#'     \code{predictor}, \code{outcomes}, \code{model_type}, and \code{raw_data}.
#'     Routes to \code{\link{multiforest}}.
#'   \item \strong{Cox models}: Classes \code{coxph} or \code{clogit}. Routes to 
#'     \code{\link{coxforest}}.
#'   \item \strong{GLM models}: Class \code{glm}. Routes to \code{\link{glmforest}}.
#'   \item \strong{Linear models}: Class \code{lm} (but not \code{glm}). Routes to 
#'     \code{\link{lmforest}}.
#' }
#' 
#' \strong{When to Use autoforest() vs Specific Functions:}
#' 
#' Use \code{autoforest()} when:
#' \itemize{
#'   \item You want quick visualization without remembering function names
#'   \item You're working interactively and exploring different model types
#'   \item You're writing generic code that handles multiple model types
#' }
#' 
#' Use the specific functions (\code{coxforest}, \code{glmforest}, etc.) when:
#' \itemize{
#'   \item You need model-specific options not available in autoforest
#'   \item You want explicit control over the visualization
#'   \item You're writing production code and want clarity
#' }
#' 
#' \strong{Title Generation:}
#' 
#' When \code{title = NULL}, titles are generated based on model type:
#' \tabular{ll}{
#'   \strong{Input Type} \tab \strong{Generated Title} \cr
#'   coxph \tab "Cox Proportional Hazards Model" \cr
#'   clogit \tab "Conditional Logistic Regression" \cr
#'   glm (binomial) \tab "Logistic Regression Model" \cr
#'   glm (poisson) \tab "Poisson Regression Model" \cr
#'   glm (other) \tab "Generalized Linear Model" \cr
#'   lm \tab "Linear Regression Model" \cr
#'   uniscreen (glm) \tab "Univariable Logistic Regression Screening" \cr
#'   uniscreen (coxph) \tab "Univariable Survival Analysis Screening" \cr
#'   uniscreen (lm) \tab "Univariable Linear Regression Screening" \cr
#'   uniscreen (glmer) \tab "Univariable Mixed Effects Logistic Screening" \cr
#'   uniscreen (lmer) \tab "Univariable Mixed Effects Linear Screening" \cr
#'   uniscreen (coxme) \tab "Univariable Mixed Effects Survival Screening" \cr
#'   multifit (glm) \tab "Multi-Outcome Logistic Regression" \cr
#'   multifit (coxph) \tab "Multi-Outcome Survival Analysis" \cr
#'   multifit (lm) \tab "Multi-Outcome Linear Regression" \cr
#'   multifit (glmer) \tab "Multi-Outcome Mixed Effects Logistic Regression" \cr
#'   multifit (lmer) \tab "Multi-Outcome Mixed Effects Linear Regression" \cr
#'   multifit (coxme) \tab "Multi-Outcome Mixed Effects Survival Analysis" \cr
#' }
#'
#' @seealso 
#' \code{\link{glmforest}} for GLM forest plots,
#' \code{\link{coxforest}} for Cox model forest plots,
#' \code{\link{lmforest}} for linear model forest plots,
#' \code{\link{uniforest}} for univariable screening forest plots,
#' \code{\link{multiforest}} for multi-outcome forest plots,
#' \code{\link{fit}} for single-model regression,
#' \code{\link{fullfit}} for combined univariable/multivariable regression,
#' \code{\link{uniscreen}} for univariable screening,
#' \code{\link{multifit}} for multi-outcome analysis
#'
#' @examples
#' # Load example data
#' data(clintrial)
#' data(clintrial_labels)
#' library(survival)
#' 
#' # Example 1: Logistic regression model
#' glm_model <- glm(surgery ~ age + sex + bmi + smoking,
#'                  family = binomial, data = clintrial)
#' 
#' plot1 <- autoforest(glm_model, data = clintrial)
#' print(plot1)
#' # Automatically detects GLM and routes to glmforest()
#' 
#' # Example 2: Cox proportional hazards model
#' cox_model <- coxph(Surv(os_months, os_status) ~ age + sex + treatment + stage,
#'                    data = clintrial)
#' 
#' plot2 <- autoforest(cox_model, data = clintrial)
#' print(plot2)
#' # Automatically detects coxph and routes to coxforest()
#' 
#' # Example 3: Linear regression model
#' lm_model <- lm(biomarker ~ age + sex + bmi + treatment, data = clintrial)
#' 
#' plot3 <- autoforest(lm_model, data = clintrial)
#' print(plot3)
#' # Automatically detects lm and routes to lmforest()
#' 
#' # Example 4: With custom labels
#' plot4 <- autoforest(
#'     glm_model,
#'     data = clintrial,
#'     labels = clintrial_labels,
#'     title = "Risk Factors for Surgical Intervention"
#' )
#' print(plot4)
#' 
#' # Example 5: Pass additional formatting options
#' plot5 <- autoforest(
#'     cox_model,
#'     data = clintrial,
#'     labels = clintrial_labels,
#'     zebra_stripes = TRUE,
#'     indent_groups = TRUE,
#'     color = "#E74C3C"
#' )
#' print(plot5)
#' 
#' # Example 6: From fit() result - data and labels extracted automatically
#' fit_result <- fit(
#'     data = clintrial,
#'     outcome = "surgery",
#'     predictors = c("age", "sex", "bmi", "treatment"),
#'     labels = clintrial_labels
#' )
#' 
#' plot6 <- autoforest(fit_result)
#' print(plot6)
#' # No need to pass data or labels - extracted from fit_result
#' 
#' # Example 7: From fullfit() result
#' ff_result <- fullfit(
#'     data = clintrial,
#'     outcome = "surgery",
#'     predictors = c("age", "sex", "bmi", "smoking", "treatment"),
#'     labels = clintrial_labels
#' )
#' 
#' plot7 <- autoforest(ff_result)
#' print(plot7)
#' 
#' # Example 8: From uniscreen() result
#' uni_result <- uniscreen(
#'     data = clintrial,
#'     outcome = "surgery",
#'     predictors = c("age", "sex", "bmi", "smoking", "treatment", "stage"),
#'     labels = clintrial_labels
#' )
#' 
#' plot8 <- autoforest(uni_result)
#' print(plot8)
#' # Automatically detects uniscreen_result and routes to uniforest()
#' 
#' # Example 9: From multifit() result
#' mf_result <- multifit(
#'     data = clintrial,
#'     outcomes = c("surgery", "complication", "readmission"),
#'     predictor = "treatment",
#'     labels = clintrial_labels
#' )
#' 
#' plot9 <- autoforest(mf_result)
#' print(plot9)
#' # Automatically detects multifit result and routes to multiforest()
#' 
#' # Example 10: Survival uniscreen
#' surv_uni <- uniscreen(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     predictors = c("age", "sex", "treatment", "stage", "grade"),
#'     model_type = "coxph",
#'     labels = clintrial_labels
#' )
#' 
#' plot10 <- autoforest(surv_uni)
#' print(plot10)
#' # Title automatically set to "Univariable Survival Analysis Screening"
#' 
#' # Example 11: Override auto-generated title
#' plot11 <- autoforest(
#'     surv_uni,
#'     title = "Univariable Predictors of Overall Survival"
#' )
#' print(plot11)
#' 
#' # Example 12: Save with recommended dimensions
#' plot12 <- autoforest(cox_model, data = clintrial, labels = clintrial_labels)
#' dims <- attr(plot12, "recommended_dims")
#' 
#' # Save to file
#' # ggsave("forest_plot.pdf", plot12, width = dims$width, height = dims$height)
#' # ggsave("forest_plot.png", plot12, width = dims$width, height = dims$height, dpi = 300)
#' 
#' # Example 13: Poisson regression
#' pois_model <- glm(n_events ~ age + sex + treatment,
#'                   family = poisson, data = clintrial)
#' 
#' plot13 <- autoforest(pois_model, data = clintrial)
#' print(plot13)
#' # Automatically detects Poisson GLM, title set to "Poisson Regression Model"
#' 
#' # Example 14: Mixed effects uniscreen
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'     me_uni <- uniscreen(
#'         data = clintrial,
#'         outcome = "surgery",
#'         predictors = c("age", "sex", "treatment"),
#'         model_type = "glmer",
#'         random = "(1 | site)",
#'         labels = clintrial_labels
#'     )
#'     
#'     plot14 <- autoforest(me_uni)
#'     print(plot14)
#'     # Title: "Univariable Mixed Effects Logistic Screening"
#' }
#' 
#' # Example 15: Quick comparison workflow
#' # Fit multiple model types and visualize each
#' models <- list(
#'     logistic = glm(surgery ~ age + sex + treatment, binomial, clintrial),
#'     survival = coxph(Surv(os_months, os_status) ~ age + sex + treatment, clintrial),
#'     linear = lm(biomarker ~ age + sex + treatment, clintrial)
#' )
#' 
#' # autoforest handles each appropriately
#' plots <- lapply(models, function(m) autoforest(m, data = clintrial))
#' 
#' # Each plot uses the correct forest function automatically
#' print(plots$logistic)
#' print(plots$survival)
#' print(plots$linear)
#'
#' @export
autoforest <- function(x, data = NULL, title = NULL, ...) {
    
    ## First check if this is a uniscreen result
    if (is_uniscreen_result(x)) {
        ## Generate appropriate title if not provided
        if (is.null(title)) {
            model_type <- attr(x, "model_type")
            title <- switch(model_type,
                            "glm" = "Univariable Logistic Regression Screening",
                            "coxph" = "Univariable Survival Analysis Screening",
                            "lm" = "Univariable Linear Regression Screening",
                            "glmer" = "Univariable Mixed Effects Logistic Screening",
                            "lmer" = "Univariable Mixed Effects Linear Screening",
                            "coxme" = "Univariable Mixed Effects Survival Screening",
                            "Univariable Screening"
                            )
        }
        
        return(uniforest(x = x, title = title, ...))
    }
    
    ## Check if this is a multifit result
    if (is_multifit_result(x)) {
        ## Generate appropriate title if not provided
        if (is.null(title)) {
            model_type <- attr(x, "model_type")
            title <- switch(model_type,
                            "glm" = "Multi-Outcome Logistic Regression",
                            "coxph" = "Multi-Outcome Survival Analysis",
                            "lm" = "Multi-Outcome Linear Regression",
                            "glmer" = "Multi-Outcome Mixed Effects Logistic Regression",
                            "lmer" = "Multi-Outcome Mixed Effects Linear Regression",
                            "coxme" = "Multi-Outcome Mixed Effects Survival Analysis",
                            "Multi-Outcome Analysis"
                            )
        }
        
        return(multiforest(x = x, title = title, ...))
    }
    
    ## Otherwise, treat as a model object
    model <- x
    model_class <- class(model)[1]
    
    ## Check for mixed-effects model classes
    is_glmer <- inherits(model, c("glmerMod", "glmerMod"))
    is_lmer <- inherits(model, "lmerMod") && !is_glmer
    is_coxme <- inherits(model, "coxme")
    
    ## Generate appropriate title if not provided
    if (is.null(title)) {
        if (is_glmer) {
            ## Check family for glmer
            family_name <- model@resp$family$family
            if (family_name == "binomial") {
                title <- "Mixed Effects Logistic Regression"
            } else if (family_name == "poisson") {
                title <- "Mixed Effects Poisson Regression"
            } else {
                title <- "Mixed Effects Generalized Linear Model"
            }
        } else if (is_lmer) {
            title <- "Mixed Effects Linear Regression"
        } else if (is_coxme) {
            title <- "Mixed Effects Cox Proportional Hazards"
        } else {
            title <- switch(model_class,
                            "coxph" = "Cox Proportional Hazards Model",
                            "clogit" = "Conditional Logistic Regression",
                            "glm" = {
                                if (model$family$family == "binomial") {
                                    "Logistic Regression Model"
                                } else if (model$family$family == "poisson") {
                                    "Poisson Regression Model" 
                                } else {
                                    "Generalized Linear Model"
                                }
                            },
                            "lm" = "Linear Regression Model",
                            "Model Results"  # Generic fallback
                            )
        }
    }
    
    ## Route to appropriate function
    ## Mixed-effects models first
    if (is_glmer) {
        glmforest(x = model, data = data, title = title, ...)
        
    } else if (is_lmer) {
        lmforest(x = model, data = data, title = title, ...)
        
    } else if (is_coxme) {
        coxforest(x = model, data = data, title = title, ...)
        
    } else if (model_class %in% c("coxph", "clogit")) {
        coxforest(x = model, data = data, title = title, ...)
        
    } else if (model_class == "glm") {
        glmforest(x = model, data = data, title = title, ...)
        
    } else if (model_class == "lm") {
        ## Ensure not actually a GLM
        if (inherits(model, "glm")) {
            glmforest(x = model, data = data, title = title, ...)
        } else {
            lmforest(x = model, data = data, title = title, ...)
        }
        
    } else {
        ## Check if it might be a supported class with different name
        if (inherits(model, "coxph")) {
            coxforest(x = model, data = data, title = title, ...)
        } else if (inherits(model, "glm")) {
            glmforest(x = model, data = data, title = title, ...)
        } else if (inherits(model, "lm")) {
            lmforest(x = model, data = data, title = title, ...)
        } else {
            stop(paste("Input class", model_class, 
                       "is not supported. Supported classes are:",
                       "lm, glm, coxph, clogit, glmerMod, lmerMod, coxme,",
                       "uniscreen_result, multifit_result"))
        }
    }
}


#' Check if object is a multifit result
#' 
#' Internal helper to detect multifit output objects.
#' 
#' @param x Object to check.
#' @return Logical indicating if x is a multifit result.
#' @keywords internal
is_multifit_result <- function(x) {
    ## Check for data.table with multifit-specific attributes
    if (!data.table::is.data.table(x)) {
        return(FALSE)
    }
    
    ## Check for required multifit attributes
    required_attrs <- c("predictor", "outcomes", "model_type", "raw_data")
    has_attrs <- vapply(required_attrs, function(a) !is.null(attr(x, a)), logical(1))
    
    all(has_attrs)
}


#' Check if object is a uniscreen result
#' 
#' Internal helper to detect uniscreen output objects.
#' 
#' @param x Object to check.
#' @return Logical indicating if x is a uniscreen result.
#' @keywords internal
is_uniscreen_result <- function(x) {
    ## Check for uniscreen_result class
    if ("uniscreen_result" %in% class(x)) {
        return(TRUE)
    }
    
    ## Also check for data.table with uniscreen-specific attributes
    if (!data.table::is.data.table(x)) {
        return(FALSE)
    }
    
    ## Check for required uniscreen attributes
    required_attrs <- c("outcome", "predictors", "model_type", "model_scope")
    has_attrs <- vapply(required_attrs, function(a) !is.null(attr(x, a)), logical(1))
    
    ## Also verify model_scope is "Univariable"
    if (all(has_attrs)) {
        return(attr(x, "model_scope") == "Univariable")
    }
    
    FALSE
}
