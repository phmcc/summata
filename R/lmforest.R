#' Create Forest Plot for Linear Models
#'
#' Generates a publication-ready forest plot that combines a formatted data table 
#' with a graphical representation of regression coefficients from a linear model. 
#' The plot integrates variable names, group levels, sample sizes, coefficients 
#' with confidence intervals, p-values, and model diagnostics (R^2, F-statistic, 
#' AIC) in a single comprehensive visualization designed for manuscripts and 
#' presentations.
#'
#' @param x Either a fitted linear model object (class \code{lm} or \code{lmerMod}), 
#'   a \code{fit_result} object from \code{\link{fit}}, or a \code{fullfit_result}
#'   object from \code{\link{fullfit}}. When a \code{fit_result} or \code{fullfit_result}
#'   is provided, the model, data, and labels are automatically extracted.
#'   Note: GLM objects are not accepted - use \code{\link{glmforest}} instead.
#'   
#' @param data A data.frame or data.table containing the original data used to 
#'   fit the model. If \code{NULL} (default) and \code{x} is a model, the function 
#'   attempts to extract data from the model object. If \code{x} is a \code{fit_result},
#'   data is extracted automatically. Providing data explicitly is recommended when
#'   passing a model directly.
#'   
#' @param title Character string specifying the plot title displayed at the top. 
#'   Default is \code{"Linear Model"}. Use descriptive titles for publication.
#'   
#' @param effect_label Character string for the effect measure label on the 
#'   forest plot axis. Default is \code{"Coefficient"}. Could be customized to 
#'   reflect units (e.g., "Change in BMI (kg/m^2)").
#'   
#' @param digits Integer specifying the number of decimal places for coefficients 
#'   and confidence intervals. Default is 2.
#'
#' @param p_digits Integer specifying the number of decimal places for p-values.
#'   P-values smaller than \code{10^(-p_digits)} are displayed as "< 0.001" 
#'   (for \code{p_digits = 3}), "< 0.0001" (for \code{p_digits = 4}), etc. 
#'   Default is 3.
#'
#' @param conf_level Numeric confidence level for confidence intervals. Must be
#'   between 0 and 1. Default is 0.95 (95\% confidence intervals). The CI
#'   percentage is automatically displayed in column headers (e.g., "90\% CI"
#'   when \code{conf_level = 0.90}).
#'   
#' @param font_size Numeric multiplier controlling the base font size for all 
#'   text elements. Default is 1.0.
#'   
#' @param annot_size Numeric value controlling the relative font size for data 
#'   annotations. Default is 3.88.
#'   
#' @param header_size Numeric value controlling the relative font size for column 
#'   headers. Default is 5.82.
#'   
#' @param title_size Numeric value controlling the relative font size for the 
#'   main plot title. Default is 23.28.
#'   
#' @param table_width Numeric value between 0 and 1 specifying the proportion of 
#'   total plot width allocated to the data table. Default is 0.6.
#'   
#' @param plot_width Numeric value specifying the intended output width in 
#'   specified \code{units}. Default is \code{NULL} (automatic).
#'   
#' @param plot_height Numeric value specifying the intended output height in 
#'   specified \code{units}. Default is \code{NULL} (automatic).
#'   
#' @param show_n Logical. If \code{TRUE}, includes a column showing group-specific 
#'   sample sizes. Note: Linear models don't have "events" so there's no 
#'   \code{show_events} parameter. Default is \code{TRUE}.
#'   
#' @param indent_groups Logical. If \code{TRUE}, indents factor levels under 
#'   their parent variable name, creating hierarchical structure. The "Group" 
#'   column is hidden when \code{TRUE}. Default is \code{FALSE}.
#'   
#' @param condense_table Logical. If \code{TRUE}, condenses binary categorical 
#'   variables into single rows. Automatically sets \code{indent_groups = TRUE}. 
#'   Default is \code{FALSE}.
#'
#' @param bold_variables Logical. If \code{TRUE}, variable names are displayed
#'   in bold. If \code{FALSE} (default), variable names are displayed in plain
#'   text.
#'   
#' @param center_padding Numeric value specifying horizontal spacing between 
#'   table and forest plot. Default is 4.
#'   
#' @param zebra_stripes Logical. If \code{TRUE}, applies alternating gray 
#'   background shading to different variables. Default is \code{TRUE}.
#'   
#' @param ref_label Character string to display for reference categories of 
#'   factor variables. Default is \code{"reference"}.
#'   
#' @param labels Named character vector providing custom display labels for 
#'   variables. Example: \code{c(age = "Age (years)", height = "Height (cm)")}. 
#'   Default is \code{NULL}.
#'   
#' @param units Character string specifying units for plot dimensions: 
#'   \code{"in"} (inches), \code{"cm"}, or \code{"mm"}. Default is \code{"in"}.
#'   
#' @param color Character string specifying the color for coefficient point 
#'   estimates in the forest plot. Default is \code{"#5A8F5A"} (green). Use 
#'   hex codes or R color names.
#'
#' @param qc_footer Logical. If \code{TRUE}, displays model quality control
#'   statistics in the footer (observations analyzed, R-squared, adjusted R-squared,
#'   F-statistic, AIC). Default is \code{TRUE}.
#'
#' @return A \code{ggplot} object containing the complete forest plot. The plot 
#'   can be displayed, saved, or further customized.
#'   
#'   The returned object includes an attribute \code{"recommended_dims"} 
#'   accessible via \code{attr(plot, "recommended_dims")}, containing:
#'   \describe{
#'     \item{width}{Numeric. Recommended plot width in specified units}
#'     \item{height}{Numeric. Recommended plot height in specified units}
#'     \item{units}{Character. The units used}
#'   }
#'   
#'   Recommendations are printed to console if dimensions are not specified.
#'
#' @details
#' \strong{Linear Model-Specific Features:}
#' 
#' The linear model forest plot differs from logistic and Cox plots in several ways:
#' \itemize{
#'   \item \strong{Coefficients}: Raw regression coefficients shown (not exponentiated)
#'   \item \strong{Reference line}: At coefficient = 0 (not at 1)
#'   \item \strong{Linear scale}: Forest plot uses linear scale (not log scale)
#'   \item \strong{No events column}: Only sample sizes shown (no event counts)
#'   \item \strong{R^2 statistics}: Model fit assessed by R² and adjusted R²
#'   \item \strong{F-test}: Overall model significance from F-statistic
#' }
#' 
#' \strong{Plot Components:}
#' 
#' \enumerate{
#'   \item \strong{Title}: Centered at top
#'   \item \strong{Data Table} (left): Contains:
#'     \itemize{
#'       \item Variable: Predictor names
#'       \item Group: Factor levels (if applicable)
#'       \item n: Sample sizes by group
#'       \item Coefficient (95\% CI); \emph{p}-value: Raw coefficients with CIs and p-values
#'     }
#'   \item \strong{Forest Plot} (right):
#'     \itemize{
#'       \item Point estimates (squares sized by sample size)
#'       \item 95\% confidence intervals (error bars)
#'       \item Reference line at coefficient = 0
#'       \item Linear scale
#'     }
#'   \item \strong{Model Statistics} (footer):
#'     \itemize{
#'       \item Observations analyzed (with percentage of total data)
#'       \item R^2 and adjusted R^2
#'       \item F-statistic with degrees of freedom and p-value
#'       \item AIC
#'     }
#' }
#' 
#' \strong{Interpreting Coefficients:}
#' 
#' Linear regression coefficients represent the change in the outcome variable 
#' for a one-unit change in the predictor:
#' \itemize{
#'   \item \strong{Continuous predictors}: Coefficient = change in Y per unit of X
#'   \item \strong{Binary predictors}: Coefficient = difference in Y between groups
#'   \item \strong{Factor predictors}: Coefficients = differences from reference 
#'     category
#'   \item \strong{Sign matters}: Positive = increase in Y, Negative = decrease in Y
#'   \item \strong{Zero crossing}: CI crossing zero suggests no significant effect
#' }
#' 
#' Example: If the coefficient for "age" is 0.50 when predicting BMI, 
#' BMI increases by 0.50 kg/m^2 for each additional year of age.
#' 
#' \strong{Model Fit Statistics:}
#' 
#' The footer displays key diagnostics:
#' \itemize{
#'   \item \strong{R^2}: Proportion of variance explained (0 to 1)
#'     \itemize{
#'       \item 0.0-0.3: Weak explanatory power
#'       \item 0.3-0.5: Moderate
#'       \item 0.5-0.7: Good
#'       \item >0.7: Strong (rare in social/biological sciences)
#'     }
#'   \item \strong{Adjusted R^2}: R^2 penalized for number of predictors
#'     \itemize{
#'       \item Always ≤ R^2
#'       \item Preferred for model comparison
#'       \item Accounts for model complexity
#'     }
#'   \item \strong{F-statistic}: Tests null hypothesis that all coefficients = 0
#'     \itemize{
#'       \item Degrees of freedom: df1 = # predictors, df2 = # observations - # predictors - 1
#'       \item Significant p-value indicates model explains variance better than intercept-only
#'     }
#'   \item \strong{AIC}: For model comparison (lower is better)
#' }
#' 
#' \strong{Assumptions:}
#' 
#' Linear regression assumes:
#' \enumerate{
#'   \item Linearity of relationships
#'   \item Independence of observations
#'   \item Homoscedasticity (constant variance)
#'   \item Normality of residuals
#'   \item No multicollinearity
#' }
#' 
#' Check assumptions using:
#' \itemize{
#'   \item \code{plot(model)} for diagnostic plots
#'   \item \code{car::vif(model)} for multicollinearity
#'   \item \code{lmtest::bptest(model)} for heteroscedasticity
#'   \item \code{shapiro.test(residuals(model))} for normality
#' }
#' 
#' \strong{Reference Categories:}
#' 
#' For factor variables:
#' \itemize{
#'   \item First level is the reference (coefficient = 0)
#'   \item Other levels show difference from reference
#'   \item Reference displayed with \code{ref_label}
#'   \item Relevel factors before modeling if needed: 
#'     \code{factor(x, levels = c("desired_ref", ...))}
#' }
#' 
#' \strong{Sample Size Reporting:}
#' 
#' The "n" column shows:
#' \itemize{
#'   \item For continuous variables: Total observations with non-missing data
#'   \item For factor variables: Number of observations in each category
#'   \item Footer shows total observations analyzed and percentage of original 
#'     data (accounting for missing values)
#' }
#'
#' @seealso 
#' \code{\link{glmforest}} for logistic/GLM forest plots,
#' \code{\link{coxforest}} for Cox model forest plots,
#' \code{\link[stats]{lm}} for fitting linear models,
#' \code{\link{fit}} for fullfit regression modeling,
#' \code{\link[stats]{plot.lm}} for diagnostic plots
#'
#' @examples
#' # Load example data
#' data(clintrial)
#' data(clintrial_labels)
#' 
#' # Example 1: Basic linear model forest plot
#' model1 <- lm(bmi ~ age + sex + smoking,
#'              data = clintrial)
#' 
#' plot1 <- lmforest(model1, data = clintrial)
#' print(plot1)
#' 
#' # Example 2: With custom labels and title
#' plot2 <- lmforest(
#'     model = model1,
#'     data = clintrial,
#'     title = "Predictors of Body Mass Index",
#'     effect_label = "Change in BMI (kg/m^2)",
#'     labels = clintrial_labels
#' )
#' print(plot2)
#' 
#' # Example 3: More comprehensive model
#' model3 <- lm(
#'     bmi ~ age + sex + smoking + hypertension + diabetes + creatinine,
#'     data = clintrial
#' )
#' 
#' plot3 <- lmforest(
#'     model = model3,
#'     data = clintrial,
#'     labels = clintrial_labels,
#'     indent_groups = TRUE
#' )
#' print(plot3)
#' 
#' # Example 4: Check model diagnostics
#' summary(model3)
#' plot(model3)  # Diagnostic plots
#' 
#' # Example 5: Condensed layout
#' model5 <- lm(
#'     bmi ~ age + sex + smoking + hypertension + diabetes,
#'     data = clintrial
#' )
#' 
#' plot5 <- lmforest(
#'     model = model5,
#'     data = clintrial,
#'     condense_table = TRUE,
#'     labels = clintrial_labels
#' )
#' print(plot5)
#' 
#' # Example 6: Custom color
#' plot6 <- lmforest(
#'     model = model1,
#'     data = clintrial,
#'     color = "#2ECC71",  # Bright green
#'     labels = clintrial_labels
#' )
#' print(plot6)
#' 
#' # Example 7: Hide sample sizes
#' plot7 <- lmforest(
#'     model = model1,
#'     data = clintrial,
#'     show_n = FALSE,
#'     labels = clintrial_labels
#' )
#' print(plot7)
#' 
#' # Example 8: Adjust table width
#' plot8 <- lmforest(
#'     model = model3,
#'     data = clintrial,
#'     table_width = 0.55,
#'     labels = clintrial_labels
#' )
#' print(plot8)
#' 
#' # Example 9: Get recommended dimensions
#' plot9 <- lmforest(model1, data = clintrial)
#' dims <- attr(plot9, "recommended_dims")
#' cat("Recommended:", dims$width, "x", dims$height, dims$units, "\n")
#' 
#' # Example 10: Specify exact dimensions
#' plot10 <- lmforest(
#'     model = model1,
#'     data = clintrial,
#'     plot_width = 12,
#'     plot_height = 7,
#'     labels = clintrial_labels
#' )
#' 
#' # Example 11: Different units
#' plot11 <- lmforest(
#'     model = model1,
#'     data = clintrial,
#'     plot_width = 30,
#'     plot_height = 20,
#'     units = "cm",
#'     labels = clintrial_labels
#' )
#' 
#' # Example 12: No zebra stripes
#' plot12 <- lmforest(
#'     model = model1,
#'     data = clintrial,
#'     zebra_stripes = FALSE,
#'     labels = clintrial_labels
#' )
#' print(plot12)
#' 
#' # Example 13: Custom reference label
#' plot13 <- lmforest(
#'     model = model1,
#'     data = clintrial,
#'     ref_label = "0.00 (ref)",
#'     labels = clintrial_labels
#' )
#' print(plot13)
#' 
#' # Example 14: Larger fonts for presentation
#' plot14 <- lmforest(
#'     model = model1,
#'     data = clintrial,
#'     font_size = 1.3,
#'     title_size = 26,
#'     labels = clintrial_labels
#' )
#' print(plot14)
#' 
#' # Example 15: More decimal places
#' plot15 <- lmforest(
#'     model = model1,
#'     data = clintrial,
#'     digits = 3,
#'     labels = clintrial_labels
#' )
#' print(plot15)
#' 
#' # Example 16: Hemoglobin as outcome
#' model16 <- lm(
#'     hemoglobin ~ age + sex + bmi + smoking + creatinine,
#'     data = clintrial
#' )
#' 
#' plot16 <- lmforest(
#'     model = model16,
#'     data = clintrial,
#'     title = "Predictors of Baseline Hemoglobin",
#'     effect_label = "Change in Hemoglobin (g/dL)",
#'     labels = clintrial_labels
#' )
#' print(plot16)
#' 
#' # Check assumptions
#' par(mfrow = c(2, 2))
#' plot(model16)
#' par(mfrow = c(1, 1))
#' 
#' # Example 17: Publication-ready final plot
#' final_model <- lm(
#'     bmi ~ age + sex + race + smoking + hypertension + 
#'         diabetes + creatinine + hemoglobin,
#'     data = clintrial
#' )
#' 
#' # Check model fit
#' summary(final_model)
#' 
#' final_plot <- lmforest(
#'     model = final_model,
#'     data = clintrial,
#'     title = "Multivariable Linear Regression: Predictors of Body Mass Index",
#'     effect_label = "Change in BMI (kg/m^2)",
#'     labels = clintrial_labels,
#'     indent_groups = TRUE,
#'     zebra_stripes = TRUE,
#'     show_n = TRUE,
#'     color = "#5A8F5A",
#'     digits = 2
#' )
#' 
#' # Save for publication
#' dims <- attr(final_plot, "recommended_dims")
#' # ggsave("figure3_linear.pdf", final_plot,
#' #        width = dims$width, height = dims$height)
#' # ggsave("figure3_linear.png", final_plot,
#' #        width = dims$width, height = dims$height, dpi = 300)
#'
#' @export
lmforest <- function(x, data = NULL,
                     title = "Linear Model",
                     effect_label = "Coefficient",
                     digits = 2,
                     p_digits = 3,
                     conf_level = 0.95,
                     font_size = 1.0,
                     annot_size = 3.88,
                     header_size = 5.82,
                     title_size = 23.28,
                     plot_width = NULL,
                     plot_height = NULL,
                     table_width = 0.6,
                     show_n = TRUE,
                     indent_groups = FALSE,
                     condense_table = FALSE,
                     bold_variables = FALSE,
                     center_padding = 4,
                     zebra_stripes = TRUE,
                     ref_label = "reference",
                     labels = NULL,
                     units = "in",
                     color = "#5A8F5A",
                     qc_footer = TRUE) {
    
    ## Check for required packages
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required but not installed.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required but not installed.")
    }
    if (!requireNamespace("grid", quietly = TRUE)) {
        stop("Package 'grid' is required but not installed.")
    }
    
    ## Handle input: accept either a model object or a fit_result/fullfit_result from fit()/fullfit()
    if (inherits(x, "fit_result") || inherits(x, "fullfit_result")) {
        ## Extract model, data, and labels from fit_result or fullfit_result
        model <- attr(x, "model")
        
        if (is.null(model)) {
            stop("The fit_result/fullfit_result does not contain a model.\n",
                 "This may occur if fullfit() was run with columns='uni' only.")
        }
        
        ## Use data from fit_result if not provided
        if (is.null(data)) {
            ## Try to get data from model
            if (inherits(model, c("lmerMod", "lmerModLmerTest", "merMod"))) {
                data <- model@frame
            } else if (!is.null(model$data)) {
                data <- model$data
            } else if (!is.null(model$model)) {
                data <- model$model
            } else {
                tryCatch({
                    data <- eval(model$call$data)
                }, error = function(e) {
                    data <- NULL
                })
            }
        }
        
        ## Use labels from fit_result if not provided
        if (is.null(labels)) {
            labels <- attr(x, "labels")
        }
        
        ## Validate that the model is a linear model (not GLM)
        is_glm <- inherits(model, "glm")
        is_valid_lm <- inherits(model, c("lm", "lmerMod", "lmerModLmerTest", "merMod")) && !is_glm
        
        if (!is_valid_lm) {
            stop("fit_result does not contain a linear model (lm or lmerMod).\n",
                 "Model type: ", class(model)[1], "\n",
                 "Use glmforest() for GLM models or coxforest() for Cox models.")
        }
    } else if (inherits(x, c("lm", "lmerMod", "lmerModLmerTest", "merMod")) && !inherits(x, "glm")) {
        ## Direct model object (lm but not glm)
        model <- x
    } else {
        stop("x must be either:\n",
             "  - An lm or lmerMod model object (not glm), or\n",
             "  - A fit_result from fit(), or\n",
             "  - A fullfit_result from fullfit()\n",
             "
Received class: ", paste(class(x), collapse = ", "))
    }
    
    ## Check model class - support both LM and LMER
    is_lme4 <- inherits(model, c("lmerMod", "lmerModLmerTest", "merMod"))
    is_lm <- inherits(model, "lm") && !inherits(model, "glm")

    ## Internally work in inches
    if (!is.null(plot_width) && units != "in") {
        plot_width <- convert_units(plot_width, from = units, to = "in")
    }
    
    ## Get model data
    if (is_lme4) {
        ## For lme4 models, the model frame contains the complete cases
        model_data <- model@frame
    } else {
        ## For regular lm, use model$model which contains complete cases
        model_data <- model$model
    }
    
    ## Get the original data for reference (if provided)
    if(is.null(data)){
        warning("The `data` argument is not provided. Data will be extracted from model fit.")
        if (is_lme4) {
            data <- model@frame
        } else {
            data <- model$data
            if (is.null(data))
                data <- model$model
        }
        if (is.null(data))
            stop("The `data` argument should be provided either to lmforest or lm/lmer.")
    }
    
    ## Convert to data.table if not already
    if(!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }
    if(!data.table::is.data.table(model_data)) {
        model_data <- data.table::as.data.table(model_data)
    }
    
    ## Extract terms based on model type
    if (is_lme4) {
        ## For lme4 models, we need to work with the fixed effects only
        ## Get the fixed effects formula (removing random effects)
        fixed_formula <- lme4::nobars(stats::formula(model))
        
        ## Get variable names from the fixed formula (excluding outcome)
        var_names <- all.vars(fixed_formula)[-1]  # Remove outcome variable
        
        ## Get the model frame for determining variable types
        model_frame <- model@frame
        
        ## Create terms object similar to LM structure
        terms <- sapply(var_names, function(v) {
            if (v %in% names(model_frame)) {
                if (is.factor(model_frame[[v]])) "factor"
                else if (is.character(model_frame[[v]])) "character"
                else "numeric"
            } else if (v %in% names(data)) {
                if (is.factor(data[[v]])) "factor"
                else if (is.character(data[[v]])) "character"
                else "numeric"
            } else "numeric"
        }, USE.NAMES = TRUE)
        
    } else {
        ## For regular lm
        terms <- attr(model$terms, "dataClasses")[-1]
    }
    
    ## Filter out interaction terms (contain ":") 
    terms <- terms[!grepl(":", names(terms), fixed = TRUE)]
    
    ## Extract coefficients and confidence intervals
    if (is_lme4) {
        ## For lme4 models
        coef_summary <- summary(model)$coefficients
        
        ## Calculate confidence intervals using Wald method with specified conf_level
        coef_vals <- lme4::fixef(model)
        se_vals <- sqrt(diag(as.matrix(vcov(model))))
        t_crit <- stats::qt((1 + conf_level) / 2, df = summary(model)$devcomp$dims["n"] - 
                                                      summary(model)$devcomp$dims["p"])
        
        ## Use profile CI if available, else fall back to Wald
        tryCatch({
            conf_int <- confint(model, method = "Wald", level = conf_level)
            ## Remove random effect rows if present
            fixed_names <- names(coef_vals)
            if (any(rownames(conf_int) %in% fixed_names)) {
                conf_int <- conf_int[fixed_names, , drop = FALSE]
            }
        }, error = function(e) {
            ci_pct_low <- paste0(round((1 - conf_level) / 2 * 100, 1), " %")
            ci_pct_high <- paste0(round((1 + conf_level) / 2 * 100, 1), " %")
            conf_int <<- cbind(
                coef_vals - t_crit * se_vals,
                coef_vals + t_crit * se_vals
            )
            colnames(conf_int) <<- c(ci_pct_low, ci_pct_high)
            rownames(conf_int) <<- names(coef_vals)
        })
        
        ## Extract t values and p values
        ## Note: lme4 may or may not provide p-values depending on the package
        if ("Pr(>|t|)" %in% colnames(coef_summary)) {
            t_values <- coef_summary[, "t value"]
            p_values <- coef_summary[, "Pr(>|t|)"]
        } else {
            ## Calculate p-values using normal approximation if not provided
            t_values <- coef_summary[, "t value"]
            ## Approximate df for p-value calculation
            df_resid <- summary(model)$devcomp$dims["n"] - 
                                     summary(model)$devcomp$dims["p"]
            p_values <- 2 * stats::pt(abs(t_values), df = df_resid, lower.tail = FALSE)
        }
        
    } else {
        ## For regular lm
        coef_summary <- summary(model)$coefficients
        conf_int <- stats::confint(model, level = conf_level)
        t_values <- coef_summary[, "t value"]
        p_values <- coef_summary[, "Pr(>|t|)"]
    }
    
    ## Remove intercept if present
    if("(Intercept)" %in% rownames(coef_summary)) {
        intercept_idx <- which(rownames(coef_summary) == "(Intercept)")
        coef_summary <- coef_summary[-intercept_idx, , drop = FALSE]
        conf_int <- conf_int[-intercept_idx, , drop = FALSE]
        t_values <- t_values[-intercept_idx]
        p_values <- p_values[-intercept_idx]
    }
    
    ## Create coefficient data table
    coef <- data.table::data.table(
                            term = rownames(coef_summary),
                            estimate = coef_summary[, "Estimate"],
                            std_error = coef_summary[, "Std. Error"],
                            statistic = t_values,
                            p_value = p_values,
                            conf_low = conf_int[, 1],
                            conf_high = conf_int[, 2]
                        )
    
    ## Get model statistics
    if (is_lme4) {
        ## For lme4 models
        model_summary <- summary(model)
        
        gmodel <- list(
            nobs = nrow(model@frame),
            ## lme4 doesn't have traditional R-squared; use conditional/marginal if available
            r_squared = NA,
            adj_r_squared = NA,
            f_statistic = NA,
            f_df1 = NA,
            f_df2 = NA,
            AIC = stats::AIC(model),
            ## Store info about random effects structure
            n_groups = sapply(model_summary$ngrps, identity)
        )
        
        ## Try to get pseudo R-squared using MuMIn if available
        tryCatch({
            if (requireNamespace("MuMIn", quietly = TRUE)) {
                r2_vals <- MuMIn::r.squaredGLMM(model)
                gmodel$r_squared <- r2_vals[1, "R2m"]  # Marginal R^2
                gmodel$conditional_r2 <- r2_vals[1, "R2c"]  # Conditional R^2
            }
        }, error = function(e) {
            ## Keep as NA
        })
        
        ## F-test not directly available for lmer
        gmodel$f_pvalue <- NA
        
    } else {
        ## For regular lm
        model_summary <- summary(model)
        gmodel <- list(
            nobs = nobs(model),
            r_squared = model_summary$r.squared,
            adj_r_squared = model_summary$adj.r.squared,
            f_statistic = model_summary$fstatistic[1],
            f_df1 = model_summary$fstatistic[2],
            f_df2 = model_summary$fstatistic[3],
            AIC = stats::AIC(model)
        )
        
        ## Calculate F-test p-value
        gmodel$f_pvalue <- pf(gmodel$f_statistic, gmodel$f_df1, gmodel$f_df2, 
                              lower.tail = FALSE)
    }
    
    ## Calculate total observations and percentage analyzed
    total_obs <- nrow(data)
    gmodel$pct_analyzed <- (gmodel$nobs / total_obs) * 100
    
    ## Format for display
    gmodel$nobs_formatted <- format(gmodel$nobs, big.mark = ",", scientific = FALSE)
    gmodel$nobs_with_pct <- paste0(gmodel$nobs_formatted, " (", 
                                   sprintf("%.1f%%", gmodel$pct_analyzed), ")")
    gmodel$AIC_formatted <- format(round(gmodel$AIC, 2), big.mark = ",", 
                                   scientific = FALSE, nsmall = 2)
    
    ## Extract xlevels (factor levels) based on model type
    if (is_lme4) {
        ## For lme4 models, extract from the model frame
        frame_data <- model@frame
        xlevels <- list()
        
        for (var_name in names(terms)) {
            if (terms[var_name] %in% c("factor", "character")) {
                if (var_name %in% names(frame_data) && is.factor(frame_data[[var_name]])) {
                    xlevels[[var_name]] <- levels(frame_data[[var_name]])
                } else if (var_name %in% names(model_data) && is.factor(model_data[[var_name]])) {
                    xlevels[[var_name]] <- levels(model_data[[var_name]])
                }
            }
        }
        
        if (length(xlevels) == 0) xlevels <- NULL
        
    } else {
        ## For regular lm
        xlevels <- model$xlevels
    }

    ## Extract statistics for every variable, preserving order
    all_terms <- lapply(seq_along(terms), function(i){
        var <- names(terms)[i]
        
        if (terms[i] %in% c("factor", "character")) {
            ## Get the factor levels from xlevels (proper order)
            if(!is.null(xlevels) && var %in% names(xlevels)) {
                factor_levels <- xlevels[[var]]
                
                ## Create data table with proper levels from model data
                level_counts <- model_data[!is.na(get(var)), .N, by = var]
                data.table::setnames(level_counts, c("level", "Freq"))
                
                ## Ensure all levels are present
                all_levels_dt <- data.table::data.table(
                                                 level = factor_levels,
                                                 Freq = 0
                                             )
                
                ## Update with actual counts
                all_levels_dt[level_counts, Freq := i.Freq, on = "level"]
                
                all_levels_dt[, var := var]
                all_levels_dt[, pos := .I]
                all_levels_dt[, var_order := i]
                return(all_levels_dt[, .(var, level, Freq, pos, var_order)])
            } else {
                ## Fallback for variables not in xlevels
                adf <- model_data[!is.na(get(var)), .N, by = var]
                data.table::setnames(adf, old = c(var, "N"), new = c("level", "Freq"))
                adf[, var := var]
                adf[, pos := .I]
                adf[, var_order := i]
                return(adf[, .(var, level, Freq, pos, var_order)])
            }
        }
        else if (terms[i] == "numeric") {
            ## For numeric variables, return a single row
            return(data.table::data.table(
                                   var = var, 
                                   level = "-", 
                                   Freq = sum(!is.na(model_data[[var]])), 
                                   pos = 1, 
                                   var_order = i
                               ))
        }
        else {
            ## Other cases
            vars = grep(paste0("^", var, "*."), coef$term, value=TRUE)
            return(data.table::data.table(
                                   var = vars, 
                                   level = "", 
                                   Freq = nrow(model_data),  # Use model_data row count
                                   pos = seq_along(vars), 
                                   var_order = i
                               ))
        }
    })
    
    all_terms_df <- data.table::rbindlist(all_terms)
    data.table::setnames(all_terms_df, c("var", "level", "N", "pos", "var_order"))

    ## Determine interaction terms via ":" strings
    interaction_terms <- coef$term[grepl(":", coef$term, fixed = TRUE)]
    
    if (length(interaction_terms) > 0) {
        ## Get the next var_order value
        next_var_order <- max(all_terms_df$var_order) + 1L
        
        ## Process each interaction term
        interaction_rows <- lapply(seq_along(interaction_terms), function(idx) {
            int_term <- interaction_terms[idx]
            
            ## Parse the interaction term to create a readable label
            int_parts <- strsplit(int_term, ":", fixed = TRUE)[[1]]
            
            ## Build readable variable name from parts
            display_parts <- sapply(int_parts, function(part) {
                ## Try to match against known factor variables to extract var and level
                matched_result <- NULL
                for (var_name in names(terms)) {
                    if (startsWith(part, var_name)) {
                        level_str <- sub(paste0("^", var_name), "", part)
                        if (nchar(level_str) > 0) {
                            ## Apply label if available
                            var_label <- if (!is.null(labels) && var_name %in% names(labels)) {
                                             labels[[var_name]]
                                         } else {
                                             var_name
                                         }
                            matched_result <- paste0(var_label, " (", level_str, ")")
                            break
                        }
                    }
                }
                ## If no match found (continuous variable in interaction or exact match)
                if (is.null(matched_result)) {
                    var_label <- if (!is.null(labels) && part %in% names(labels)) {
                                     labels[[part]]
                                 } else {
                                     part
                                 }
                    matched_result <- var_label
                }
                matched_result
            })
            
            display_name <- paste(display_parts, collapse = " * ")
            
            ## Calculate N for this interaction (observations where both conditions are met)
            int_n <- tryCatch({
                n_count <- NA_integer_
                
                ## Try to extract the subset count from the model data
                ## For each part of the interaction, build a condition
                conditions <- list()
                for (j in seq_along(int_parts)) {
                    part <- int_parts[j]
                    for (var_name in names(terms)) {
                        if (startsWith(part, var_name) && terms[var_name] %in% c("factor", "character")) {
                            level_str <- sub(paste0("^", var_name), "", part)
                            if (nchar(level_str) > 0 && var_name %in% names(model_data)) {
                                conditions[[var_name]] <- level_str
                            }
                            break
                        }
                    }
                }
                
                ## If conditions found, count matching rows
                if (length(conditions) > 0) {
                    ## Build subset expression safely
                    subset_dt <- data.table::copy(model_data)
                    for (cond_var in names(conditions)) {
                        subset_dt <- subset_dt[get(cond_var) == conditions[[cond_var]]]
                    }
                    n_count <- nrow(subset_dt)
                }
                
                n_count
            }, error = function(e) NA_integer_)
            
            data.table::data.table(
                            var = display_name,
                            level = "-",
                            N = int_n,
                            pos = 1L,
                            var_order = next_var_order + idx - 1L,
                            .int_term = int_term  # Store raw term for later inds assignment
                        )
        })
        
        interaction_df <- data.table::rbindlist(interaction_rows)
        
        ## Bind interaction rows to main terms
        all_terms_df <- data.table::rbindlist(list(all_terms_df, interaction_df), 
                                              use.names = TRUE, fill = TRUE)
    }

    ## Apply condensing and indenting if requested
    if (condense_table || indent_groups) {
        if (condense_table) {
            indent_groups <- TRUE
        }
        
        ## For interaction rows, use the stored raw term; for others use var+level pattern
        if (".int_term" %in% names(all_terms_df)) {
            all_terms_df[, inds := data.table::fifelse(
                                                   !is.na(.int_term), 
                                                   .int_term,
                                                   data.table::fifelse(level == "-", var, paste0(var, level))
                                               )]
        } else {
            all_terms_df[, inds := data.table::fifelse(level == "-", var, paste0(var, level))]
        }
        orig_inds_map <- data.table::copy(all_terms_df[, .(var, level, inds, N)])
        
        processed_rows <- list()
        row_counter <- 1
        unique_vars <- unique(all_terms_df[, var])
        
        for (v in unique_vars) {
            var_rows <- all_terms_df[var == v]
            
            if (nrow(var_rows) == 1) {
                ## Single-level variable (continuous)
                processed_rows[[row_counter]] <- var_rows
                row_counter <- row_counter + 1
            } else {
                ## Multi-level variable (categorical)
                is_binary <- nrow(var_rows) == 2
                
                if (condense_table && is_binary) {
                    ## Only condense if both flags are set
                    non_ref_row <- var_rows[2]
                    condensed_row <- data.table::copy(non_ref_row)
                    condensed_row[, var := paste0(v, " (", level, ")")]
                    condensed_row[, level := "-"]
                    processed_rows[[row_counter]] <- condensed_row
                    row_counter <- row_counter + 1
                } else if (indent_groups) {
                    ## Add header for any multi-level variable when indenting
                    header_row <- data.table::data.table(
                                                  var = v,
                                                  level = "-",
                                                  N = NA_integer_,
                                                  pos = var_rows$pos[1],
                                                  var_order = var_rows$var_order[1],
                                                  inds = NA_character_
                                              )
                    processed_rows[[row_counter]] <- header_row
                    row_counter <- row_counter + 1
                    
                    ## Add indented levels
                    for (i in 1:nrow(var_rows)) {
                        group_row <- data.table::copy(var_rows[i])
                        group_row[, var := paste0("    ", level)]
                        group_row[, level := "-"]
                        processed_rows[[row_counter]] <- group_row
                        row_counter <- row_counter + 1
                    }
                } else {
                    ## Standard layout
                    for (i in 1:nrow(var_rows)) {
                        processed_rows[[row_counter]] <- var_rows[i]
                        row_counter <- row_counter + 1
                    }
                }
            }
        }
        
        all_terms_df <- data.table::rbindlist(processed_rows, fill = TRUE)
        
        for (i in 1:nrow(all_terms_df)) {
            current_var <- all_terms_df$var[i]
            
            if (is.na(all_terms_df$inds[i]) && !grepl("^    ", current_var) && !grepl("\\(", current_var)) {
                next
            }
            
            if (grepl("\\(", current_var)) {
                orig_var <- gsub(" \\(.*\\)", "", current_var)
                orig_level <- gsub(".*\\((.*)\\)", "\\1", current_var)
                
                matching <- orig_inds_map[var == orig_var & level == orig_level]
                if (nrow(matching) > 0) {
                    all_terms_df[i, `:=`(inds = matching$inds[1],
                                         N = matching$N[1])]
                }
            }
            else if (grepl("^    ", current_var)) {
                clean_level <- gsub("^    ", "", current_var)
                
                parent_var <- NA_character_
                for (j in (i-1):1) {
                    if (!grepl("^    ", all_terms_df$var[j]) && all_terms_df$var[j] != "") {
                        parent_var <- gsub(" \\(.*\\)", "", all_terms_df$var[j])
                        break
                    }
                }
                
                if (!is.na(parent_var)) {
                    matching <- orig_inds_map[var == parent_var & level == clean_level]
                    if (nrow(matching) > 0) {
                        all_terms_df[i, `:=`(inds = matching$inds[1],
                                             N = matching$N[1])]
                    }
                }
            }
        }
    } else {
        ## For non-condensed/non-indented display, set inds for coefficient matching
        ## For interaction rows, use the stored raw term; for others use var+level pattern
        if (".int_term" %in% names(all_terms_df)) {
            all_terms_df[, inds := data.table::fifelse(
                                                   !is.na(.int_term), 
                                                   .int_term,
                                                   data.table::fifelse(level == "-", var, paste0(var, level))
                                               )]
        } else {
            all_terms_df[, inds := data.table::fifelse(level == "-", var, paste0(var, level))]
        }
    }
    
    ## Clean up temporary column if it exists
    if (".int_term" %in% names(all_terms_df)) {
        all_terms_df[, .int_term := NULL]
    }
    
    ## Process coefficients
    coef[, term := gsub(term, pattern = "`", replacement = "")]
    coef[, inds := term]
    
    ## Merge data
    to_show <- merge(all_terms_df, coef, by.x = "inds", by.y = "inds", all.x = TRUE, sort = FALSE)
    
    ## Sort by variable order first, then position within variable
    data.table::setorder(to_show, var_order, pos)
    
    ## Add variable-based shading indicator (zebra stripes)
    if (zebra_stripes) {
        to_show[, shade_group := var_order %% 2]
        shade_colors <- c("#FFFFFF", "#EEEEEE")
    } else {
        to_show[, shade_group := 0]
        shade_colors <- c("#FFFFFF", "#FFFFFF")
    }
    
    ## Select columns
    to_show <- to_show[, .(var, level, N, p_value, estimate, conf_low, conf_high, pos, var_order, shade_group)]

    ## Format the values
    to_show_exp_clean <- data.table::copy(to_show)
    
    ## Create formatted columns for display
    to_show_exp_clean[, effect := estimate]
    to_show_exp_clean[, effect_formatted := data.table::fifelse(is.na(N) & is.na(estimate),
                                                                "",
                                                                data.table::fifelse(is.na(estimate), 
                                                                                    ref_label,
                                                                                    format(round(estimate, digits), nsmall = digits)))]
    to_show_exp_clean[, conf_low_formatted := data.table::fifelse(is.na(conf_low), 
                                                                  NA_character_,
                                                                  format(round(conf_low, digits), nsmall = digits))]
    to_show_exp_clean[, conf_high_formatted := data.table::fifelse(is.na(conf_high), 
                                                                   NA_character_,
                                                                   format(round(conf_high, digits), nsmall = digits))]
    
    ## Format CI percentage for display in headers
    ci_pct <- round(conf_level * 100)

    ## Format p-values using p_digits parameter
    p_threshold <- 10^(-p_digits)
    p_threshold_str <- paste0("< ", format(p_threshold, scientific = FALSE))
    
    to_show_exp_clean[, p_formatted := data.table::fifelse(is.na(p_value), 
                                                           NA_character_,
                                                           data.table::fifelse(p_value < p_threshold, 
                                                                               p_threshold_str,
                                                                               format(round(p_value, p_digits), nsmall = p_digits)))]
    
    ## Create the combined effect string with expression for italic p
    to_show_exp_clean[, effect_string_expr := data.table::fifelse(
                                                              is.na(N) & is.na(estimate),
                                                              "''",
                                                              data.table::fcase(
                                                                  is.na(estimate), paste0("'", ref_label, "'"),
                                                                  
                                                                  p_value < p_threshold & (conf_low < 0 | conf_high < 0),
                                                                  paste0("'", effect_formatted, " (", conf_low_formatted, ", ", 
                                                                         conf_high_formatted, "); '*~italic(p)~'", p_threshold_str, "'"),
                                                                  
                                                                  p_value < p_threshold,
                                                                  paste0("'", effect_formatted, " (", conf_low_formatted, "-", 
                                                                         conf_high_formatted, "); '*~italic(p)~'", p_threshold_str, "'"),
                                                                  
                                                                  conf_low < 0 | conf_high < 0,
                                                                  paste0("'", effect_formatted, " (", conf_low_formatted, ", ", 
                                                                         conf_high_formatted, "); '*~italic(p)~'= ", p_formatted, "'"),
                                                                  
                                                                  default = paste0("'", effect_formatted, " (", conf_low_formatted, "-", 
                                                                                   conf_high_formatted, "); '*~italic(p)~'= ", p_formatted, "'")
                                                              )
                                                          )]
    
    ## Format N with thousands separator
    to_show_exp_clean[, n_formatted := data.table::fifelse(is.na(N), "", format(N, big.mark = ",", scientific = FALSE))]
    
    ## Clean up variable names for display
    to_show_exp_clean[, var_display := as.character(var)]
    
    if (indent_groups || condense_table) {
        to_show_exp_clean[, var_display := var]
        
        for (v in unique(to_show_exp_clean$var)) {
            if (v != "" && !grepl("^    ", v)) {
                ## Skip interaction terms (already formatted) - identified by " * " in the name
                if (grepl(" \\* ", v)) next
                
                clean_v <- gsub(" \\(.*\\)", "", v)

                label <- if (!is.null(labels) && clean_v %in% names(labels)) {
                             labels[clean_v]
                         } else if (clean_v %in% names(data) && !is.null(attr(data[[clean_v]], "label"))) {
                             attr(data[[clean_v]], "label")
                         } else {
                             NULL
                         }
                
                if (!is.null(label)) {
                    if (grepl("\\(", v)) {
                        category <- gsub(".*\\((.*)\\)", "\\1", v)
                        if (category %in% c("Yes", "YES", "yes", "1", "True", "TRUE", "true", "Present", "Positive", "+")) {
                            to_show_exp_clean[var == v, var_display := label]
                        } else {
                            to_show_exp_clean[var == v, var_display := paste0(label, " (", category, ")")]
                        }
                    } else {
                        to_show_exp_clean[var == v, var_display := label]
                    }
                }
            }
        }
        
        to_show_exp_clean[, level := ""]
    } else {
        for(v in unique(to_show_exp_clean$var)) {
            if(v %in% to_show_exp_clean$var) {
                ## Skip interaction terms (already formatted) - identified by " * " in the name
                if (grepl(" \\* ", v)) next
                
                if(!is.null(labels) && v %in% names(labels)) {
                    to_show_exp_clean[var == v, var_display := labels[v]]
                }
                else if(v %in% names(data) && !is.null(attr(data[[v]], "label"))) {
                    to_show_exp_clean[var == v, var_display := attr(data[[v]], "label")]
                }
            }
        }
    }

    if (!indent_groups) {
        to_show_exp_clean[duplicated(var), var_display := ""]
    }
    
    ## Handle missing estimates for plotting
    to_show_exp_clean[is.na(estimate), estimate := 0]
    
    ## Reorder (flip), but maintain the variable grouping
    to_show_exp_clean <- to_show_exp_clean[order(nrow(to_show_exp_clean):1)]
    to_show_exp_clean[, x_pos := .I]
    
    ## Add shade_group based on var_order for alternating by variable
    to_show_exp_clean[, shade_group := var_order %% 2]
    
    ## Number of rows for later adjustments
    n_rows <- nrow(to_show_exp_clean)
    
    ## Calculate plot ranges
    rangeb <- range(to_show$conf_low, to_show$conf_high, na.rm = TRUE)
    breaks <- pretty(rangeb, n = 7)
    breaks <- breaks[breaks >= rangeb[1] & breaks <= rangeb[2]]
    
    ## Check for one-sided case
    is_one_sided <- (min(rangeb) > 0) || (max(rangeb) < 0)
    
    ## Ensure 0 is included for coefficients
    if(!0 %in% breaks) {
        breaks <- sort(unique(c(breaks, 0)))
    }
    
    ## Extend range if needed
    if(min(rangeb) > 0) {
        rangeb[1] <- -0.1 * abs(max(rangeb))
    } else if(max(rangeb) < 0) {
        rangeb[2] <- 0.1 * abs(min(rangeb))
    }
    
    reference_value <- 0

    ## Calculate plot ranges
    rangeb <- range(to_show$conf_low, to_show$conf_high, na.rm = TRUE)
    breaks <- pretty(rangeb, n = 7)
    breaks <- breaks[breaks >= rangeb[1] & breaks <= rangeb[2]]

    ## Check for one-sided case
    is_one_sided <- (min(rangeb) > 0) || (max(rangeb) < 0)

    ## Ensure 0 is included for coefficients
    if(!0 %in% breaks) {
        breaks <- sort(unique(c(breaks, 0)))
    }

    ## Extend range if needed
    if(min(rangeb) > 0) {
        rangeb[1] <- -0.1 * abs(max(rangeb))
    } else if(max(rangeb) < 0) {
        rangeb[2] <- 0.1 * abs(min(rangeb))
    }

    reference_value <- 0

    ## Calculate layout using helper function
    layout <- calculate_forest_layout(
        to_show_exp_clean = to_show_exp_clean,
        show_n = show_n,
        show_events = FALSE,  # linear models do not have events
        indent_groups = indent_groups,
        condense_table = condense_table,
        effect_label = effect_label,
        ref_label = ref_label,
        font_size = font_size,
        table_width = data.table::fifelse(is.null(table_width), 0.6, table_width),
        rangeb = rangeb,
        center_padding = center_padding
    )

    ## Set up the extended range for plotting
    rangeplot <- c(layout$rangeplot_start, rangeb[2] + diff(rangeb) * 0.05)

    ## Extract positions
    y_variable <- layout$positions$var
    if (!(indent_groups || condense_table)) {
        y_level <- layout$positions$level
    }
    if (show_n) {
        y_n <- layout$positions$n
    }
    y_coef <- layout$positions$effect

    ## Use the effect abbreviation from layout
    effect_abbrev <- layout$effect_abbrev

    ## Calculate recommended dimensions
    rec_height <- max(5, min(20, 3 + nrow(to_show_exp_clean) * 0.25))

    if(!is.null(plot_width)) {
        rec_width <- plot_width
        if(!is.null(plot_height)) {
            rec_height <- plot_height
        }
    } else {
        ## Use the calculated total width
        rec_width <- layout$total_width + 1.0  # Add margins
        
        ## Apply reasonable bounds based on layout
        if (condense_table || indent_groups) {
            if (!show_n) {
                rec_width <- max(7.5, min(11, rec_width))
            } else {
                rec_width <- max(9, min(14, rec_width))
            }
        } else {
            if (!show_n) {
                rec_width <- max(12, min(15, rec_width))
            } else {
                rec_width <- max(13, min(16, rec_width))
            }
        }
    }

    ## Font sizes
    annot_font <- font_size * annot_size
    header_font <- font_size * header_size

    ## Custom ticks data
    ticks_df <- data.frame(
        x = -0.5,
        xend = -0.7,
        y = breaks,
        yend = breaks
    )
    
    ## Create the plot
    p <- ggplot2::ggplot(to_show_exp_clean, ggplot2::aes(x_pos, estimate)) +
        
        ## Shading rectangles
        ggplot2::geom_rect(ggplot2::aes(xmin = x_pos - .5, xmax = x_pos + .5,
                                        ymin = rangeplot[1], ymax = rangeplot[2],
                                        fill = ordered(shade_group))) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.05))) +
        ggplot2::scale_size_area(max_size = 6, guide = "none") +
        ggplot2::scale_fill_manual(values = shade_colors, guide = "none") +
        
        ## Forest plot elements
        ggplot2::geom_point(ggplot2::aes(size = N), pch = 22, color = "#000000", fill = color, na.rm = TRUE) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = conf_low, ymax = conf_high), width = 0.15) +
        
        ## Y-axis for forest plot
        ggplot2::annotate(geom = "segment",
                          x = -0.5, xend = -0.5,
                          y = rangeb[1],
                          yend = rangeb[2],
                          color = "#000000", linewidth = 1) +
        
        ## Reference line at 0
        ggplot2::annotate(geom = "segment", 
                          x = -0.5, xend = max(to_show_exp_clean$x_pos) + 0.5, 
                          y = reference_value, yend = reference_value, 
                          linetype = "longdash") +
        
        ## Ticks
        ggplot2::geom_segment(data = ticks_df,
                              ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                              inherit.aes = FALSE,
                              color = "#000000",
                              linewidth = 1) +
        ggplot2::geom_text(data = ticks_df,
                           ggplot2::aes(x = xend - 0.05, y = y, label = sprintf("%g", y)),
                           inherit.aes = FALSE,
                           hjust = 0.5,
                           vjust = 1.3,
                           size = annot_font * 1.5) +
        ## Set coordinate system
        ggplot2::coord_flip(ylim = rangeplot) +
        ggplot2::ggtitle(title) +
        ggplot2::scale_y_continuous(name = effect_label,
                                    labels = sprintf("%g", breaks),
                                    expand = c(0.02, 0.02),
                                    breaks = breaks) +
        ggplot2::theme_light() +
        ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0),
                       panel.grid.minor.y = ggplot2::element_blank(),
                       panel.grid.minor.x = ggplot2::element_blank(),
                       panel.grid.major.y = ggplot2::element_blank(),
                       panel.grid.major.x = ggplot2::element_blank(),
                       legend.position = "none",
                       panel.border = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(size = font_size * title_size, face = "bold", hjust = 0.5)) +
        ggplot2::xlab("") +
        
        ## Variable column
        ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = y_variable,
                          label = "Variable", fontface = "bold", hjust = 0,
                          size = header_font) +
        
        {if (indent_groups || condense_table) {
             fontfaces <- if (bold_variables) {
                 data.table::fifelse(grepl("^    ", to_show_exp_clean$var_display), "plain", "bold")
             } else {
                 rep("plain", nrow(to_show_exp_clean))
             }
             ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_variable,
                               label = to_show_exp_clean$var_display, 
                               fontface = fontfaces, 
                               hjust = 0,
                               size = annot_font)
         } else {
             ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_variable,
                               label = to_show_exp_clean$var_display, 
                               fontface = if (bold_variables) "bold" else "plain", 
                               hjust = 0,
                               size = annot_font)
         }} +
        
        ## Group/level column
        {if (!(indent_groups || condense_table)) {
             list(
                 ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = y_level,
                                   label = "Group", fontface = "bold", hjust = 0,
                                   size = header_font),
                 ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_level,
                                   label = to_show_exp_clean$level, hjust = 0,
                                   size = annot_font)
             )
         }} +
        
        ## N column (conditional)
        {if (show_n) {
             list(
                 ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = y_n,
                                   label = "n", fontface = "bold.italic", hjust = 0.5,
                                   size = header_font),
                 ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_n,
                                   label = to_show_exp_clean$n_formatted, hjust = 0.5,
                                   size = annot_font)
             )
         }} +
        
        ## Coefficient column - use dynamic CI percentage
        ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.4, y = y_coef,
                          label = paste0("bold('", effect_abbrev, " (", ci_pct, "% CI); '*bolditalic(p)*'-value')"),
                          hjust = 0, size = header_font, parse = TRUE) +
        
        ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_coef,
                          label = to_show_exp_clean$effect_string_expr, hjust = 0,
                          size = annot_font, parse = TRUE) +
        
        ## X-axis label
        ggplot2::annotate(geom = "text", x = -1.5, y = reference_value,
                          label = effect_label, fontface = "bold",
                          hjust = 0.5, vjust = 2, size = annot_font * 1.5) +
        
        ## Model statistics at bottom (conditional)
        {if (qc_footer) {
             ggplot2::annotate(geom = "text", x = 0.5, y = y_variable,
                               label = paste0("Observations analyzed: ", gmodel$nobs_with_pct,
                                              "\nR^2: ", round(gmodel$r_squared, 3),
                                              " (Adjusted: ", round(gmodel$adj_r_squared, 3), ")",
                                              "\nF-statistic: ", round(gmodel$f_statistic, 2),
                                              " (df1 = ", gmodel$f_df1, ", df2 = ", gmodel$f_df2,
                                              "); p ", data.table::fifelse(gmodel$f_pvalue < p_threshold, p_threshold_str, 
                                                                           paste0("= ", format(round(gmodel$f_pvalue, p_digits), nsmall = p_digits))),
                                              "\nAIC: ", gmodel$AIC_formatted),
                               size = annot_font * 0.8, hjust = 0, vjust = 1.2, fontface = "italic")
        }}

    ## Convert units back for output if needed
    if (units != "in") {
        rec_width <- convert_units(rec_width, from = "in", to = units)
        rec_height <- convert_units(rec_height, from = "in", to = units)
    }

    ## Provide dimension recommendations
    if(is.null(plot_width) || is.null(plot_height)) {
        message(sprintf("Recommended plot dimensions: width = %.1f %s, height = %.1f %s",
                        rec_width, units, rec_height, units))
    }
    
    ## Add recommended dimensions as an attribute
    attr(p, "recommended_dims") <- list(width = rec_width, height = rec_height, units = units)

    ## Return the plot
    return(p)
}
