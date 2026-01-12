#' Create Forest Plot for Univariable Screening Results
#'
#' Generates a publication-ready forest plot from a \code{\link{uniscreen}} output
#' object. The plot displays effect estimates (OR, HR, RR, or coefficients) with
#' confidence intervals for each predictor tested in univariable analysis against
#' a single outcome.
#'
#' @param x A uniscreen result object (data.table with class \code{uniscreen_result}
#'   from \code{\link{uniscreen}}).
#'   
#' @param title Character string specifying the plot title. Default is 
#'   \code{"Univariable Screening"}. Use descriptive titles for publication.
#'   
#' @param effect_label Character string for the effect measure label on the 
#'   forest plot axis. Default is \code{NULL}, which auto-detects based on 
#'   model type (e.g., "Odds Ratio", "Hazard Ratio", "Rate Ratio", "Coefficient").
#'   
#' @param digits Integer specifying the number of decimal places for effect 
#'   estimates and confidence intervals. Default is 2.
#'
#' @param p_digits Integer specifying the number of decimal places for p-values.
#'   P-values smaller than \code{10^(-p_digits)} are displayed as "< 0.001" 
#'   (for \code{p_digits = 3}). Default is 3.
#'
#' @param conf_level Numeric confidence level for confidence intervals. Must be
#'   between 0 and 1. Default is 0.95 (95\% confidence intervals). The CI
#'   percentage is automatically displayed in column headers (e.g., "90\% CI"
#'   when \code{conf_level = 0.90}). Note: This parameter affects display only;
#'   the underlying CIs come from the uniscreen result.
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
#'   total plot width allocated to the data table. Default is 0.6 (60\% table, 
#'   40\% forest plot).
#'   
#' @param plot_width Numeric value specifying the intended output width in 
#'   specified \code{units}. Used for optimizing layout. Default is \code{NULL} 
#'   (automatic). Recommended: 10-16 inches.
#'   
#' @param plot_height Numeric value specifying the intended output height in 
#'   specified \code{units}. Default is \code{NULL} (automatic based on number 
#'   of rows).
#'   
#' @param show_n Logical. If \code{TRUE}, includes a column showing sample sizes. 
#'   Default is \code{TRUE}.
#'   
#' @param show_events Logical. If \code{TRUE}, includes a column showing the 
#'   number of events for each row. Default is \code{NULL}, which auto-detects
#'   based on model type (TRUE for binomial/survival, FALSE for linear).
#'   
#' @param indent_groups Logical. If \code{TRUE}, indents factor levels
#'   under their parent variable name, creating a hierarchical display. If 
#'   \code{FALSE} (default), shows variable and level in separate columns.
#'
#' @param condense_table Logical. If \code{TRUE}, condenses binary categorical 
#'   variables into single rows by showing only the non-reference category. 
#'   Automatically sets \code{indent_groups = TRUE}. Useful for tables with 
#'   many binary variables. Default is \code{FALSE}.
#'
#' @param bold_variables Logical. If \code{TRUE}, variable names are displayed
#'   in bold. If \code{FALSE} (default), variable names are displayed in plain
#'   text.
#'   
#' @param center_padding Numeric value specifying horizontal spacing between 
#'   table and forest plot. Default is 4.
#'   
#' @param zebra_stripes Logical. If \code{TRUE}, applies alternating gray 
#'   background shading to different variables for improved readability. 
#'   Default is \code{TRUE}.
#'   
#' @param color Character string specifying the color for point estimates in 
#'   the forest plot. Default is \code{NULL}, which auto-selects based on
#'   model type (purple for Cox, teal for GLM, blue for Poisson, green for LM).
#'   Use hex codes or R color names for custom colors.
#'   
#' @param null_line Numeric value for the reference line position. Default is 
#'   \code{NULL}, which uses 1 for ratio measures (OR, HR, RR) and 0 for 
#'   coefficients.
#'   
#' @param log_scale Logical. If \code{TRUE}, uses log scale for the x-axis.
#'   Default is \code{NULL}, which auto-detects (TRUE for OR/HR/RR, FALSE for
#'   coefficients).
#'   
#' @param labels Named character vector providing custom display labels for 
#'   variables. Applied to predictor names in the plot.
#'   Default is \code{NULL} (uses original variable names).
#'
#' @param show_footer Logical. If \code{TRUE}, displays a footer with the
#'   outcome variable name. Default is \code{TRUE}.
#'   
#' @param units Character string specifying units for plot dimensions: 
#'   \code{"in"} (inches), \code{"cm"}, or \code{"mm"}. Default is \code{"in"}.
#'
#' @return A \code{ggplot} object containing the complete forest plot.
#'   
#'   The returned object includes an attribute \code{"recommended_dims"} 
#'   accessible via \code{attr(plot, "recommended_dims")}, containing:
#'   \describe{
#'     \item{width}{Numeric. Recommended plot width in specified units}
#'     \item{height}{Numeric. Recommended plot height in specified units}
#'   }
#'
#' @details
#' The forest plot displays univariable (unadjusted) associations between each
#' predictor and the outcome. This is useful for:
#' \itemize{
#'   \item Visualizing results of initial variable screening
#'   \item Identifying potential predictors for multivariable modeling
#'   \item Presenting crude associations alongside adjusted results
#'   \item Quick visual assessment of effect sizes and significance
#' }
#'
#' The plot automatically handles:
#' \itemize{
#'   \item Different effect types (OR, HR, RR, coefficients) with appropriate
#'     axis scaling (log vs linear)
#'   \item Factor variables with multiple levels (grouped under variable name)
#'   \item Continuous variables (single row per predictor)
#'   \item Reference categories for categorical variables
#' }
#'
#' @seealso 
#' \code{\link{uniscreen}} for generating univariable screening results,
#' \code{\link{multiforest}} for multi-outcome forest plots,
#' \code{\link{coxforest}}, \code{\link{glmforest}}, \code{\link{lmforest}} for
#' single-model forest plots
#'
#' @examples
#' # Load example data
#' data(clintrial)
#' data(clintrial_labels)
#' 
#' # Example 1: Basic logistic regression screening
#' uni_results <- uniscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "smoking", "treatment", "stage"),
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' 
#' uniforest(uni_results, title = "Univariable Associations with Mortality")
#' 
#' \donttest{
#' # Example 2: Survival analysis
#' library(survival)
#' surv_results <- uniscreen(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     predictors = c("age", "sex", "treatment", "stage"),
#'     model_type = "coxph",
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' 
#' uniforest(surv_results, title = "Univariable Survival Analysis")
#' 
#' # Example 3: Linear regression
#' lm_results <- uniscreen(
#'     data = clintrial,
#'     outcome = "los_days",
#'     predictors = c("age", "sex", "surgery", "diabetes"),
#'     model_type = "lm",
#'     labels = clintrial_labels,
#'     parallel = FALSE
#' )
#' 
#' uniforest(lm_results, title = "Predictors of Length of Stay")
#' 
#' # Example 4: Customize appearance
#' uniforest(
#'     uni_results,
#'     title = "Crude Associations with Mortality",
#'     color = "#E74C3C",
#'     indent_groups = TRUE,
#'     zebra_stripes = TRUE,
#'     bold_variables = TRUE
#' )
#' 
#' # Example 5: Hide footer
#' uniforest(uni_results, 
#'           title = "Predictors of Overall Survival Status",
#'           show_footer = FALSE)
#' 
#' # Example 6: Save with recommended dimensions
#' p <- uniforest(uni_results)
#' dims <- attr(p, "recommended_dims")
#' # ggsave("univariable_forest.pdf", p, width = dims$width, height = dims$height)
#' }
#'
#' @export
uniforest <- function(x,
                      title = "Univariable Screening",
                      effect_label = NULL,
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
                      show_events = NULL,
                      indent_groups = FALSE,
                      condense_table = FALSE,
                      bold_variables = FALSE,
                      center_padding = 4,
                      zebra_stripes = TRUE,
                      color = NULL,
                      null_line = NULL,
                      log_scale = NULL,
                      labels = NULL,
                      show_footer = TRUE,
                      units = "in") {
    
    ## Input validation
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required but not installed.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required but not installed.")
    }
    
    ## Check if user passed correct input type
    if (!inherits(x, "uniscreen_result")) {
        if (inherits(x, c("coxph", "glm", "lm", "glmerMod", "lmerMod", "coxme"))) {
            stop("uniforest() expects a uniscreen() result, not a model object.\n",
                 "\nFor single-model forest plots, use:\n",
                 "  - coxforest() for Cox models\n",
                 "  - glmforest() for GLM models\n",
                 "  - lmforest() for linear models\n",
                 "\nFor univariable screening visualization, first run uniscreen(), then pass the result to uniforest().")
        }
        if (inherits(x, "fit_result") || inherits(x, "fullfit_result")) {
            stop("uniforest() expects a uniscreen() result, not a fit()/fullfit() result.\n",
                 "\nFor fit_result/fullfit_result objects, use the appropriate forest function:\n",
                 "  - coxforest() for Cox models\n",
                 "  - glmforest() for GLM models\n",
                 "  - lmforest() for linear models")
        }
        stop("uniforest() expects a uniscreen() result.\n",
             "Input is of class: ", paste(class(x), collapse = ", "))
    }
    
    ## Extract raw data from uscreen result
    raw_data <- attr(x, "raw_data")
    
    if (is.null(raw_data) || !data.table::is.data.table(raw_data)) {
        stop("Input does not appear to be a valid uniscreen() result (missing raw_data attribute).")
    }
    
    ## Get attributes from uscreen result
    us_outcome <- attr(x, "outcome")
    us_model_type <- attr(x, "model_type")
    
    ## Determine effect type from raw data column names
    ## m2dt creates columns named OR, HR, RR, or Estimate
    if ("OR" %in% names(raw_data)) {
        effect_type <- "OR"
    } else if ("HR" %in% names(raw_data)) {
        effect_type <- "HR"
    } else if ("RR" %in% names(raw_data)) {
        effect_type <- "RR"
    } else {
        effect_type <- "Estimate"
    }
    
    ## Set defaults based on effect type
    if (is.null(null_line)) {
        null_line <- if (effect_type %in% c("OR", "HR", "RR")) 1 else 0
    }
    
    if (is.null(log_scale)) {
        log_scale <- effect_type %in% c("OR", "HR", "RR")
    }
    
    ## Set effect label based on effect type
    if (is.null(effect_label)) {
        effect_label <- switch(effect_type,
                               "OR" = "Odds Ratio",
                               "HR" = "Hazard Ratio",
                               "RR" = "Rate Ratio",
                               "Coefficient")
    }
    
    ## Set effect abbreviation for column header
    ## For univariable screening, use non-adjusted abbreviations
    effect_abbrev <- switch(effect_type,
                            "OR" = "OR",
                            "HR" = "HR",
                            "RR" = "RR",
                            "Coefficient")
    
    ## Format CI percentage for display in headers
    ci_pct <- round(conf_level * 100)
    
    ## Create header expression with italic p
    effect_header_expr <- paste0("bold('", effect_abbrev, " (", ci_pct, "% CI); '*bolditalic(p)*'-value')")
    
    ## Set default show_events based on model type
    if (is.null(show_events)) {
        show_events <- effect_type %in% c("OR", "HR", "RR")
    }
    
    ## Set default color based on effect type
    if (is.null(color)) {
        color <- switch(effect_type,
                        "HR" = "#8A61D8",
                        "OR" = "#3C8D9C",
                        "RR" = "#3064A6",
                        "#5A8F5A")
    }
    
    ## Internally work in inches (uses convert_units from forest_utils.R)
    if (!is.null(plot_width) && units != "in") {
        plot_width <- convert_units(plot_width, from = units, to = "in")
    }
    
    ## Prepare data for plotting
    to_show <- data.table::copy(raw_data)
    
    ## Handle estimate column
    if ("coefficient" %in% names(to_show) && !"estimate" %in% names(to_show)) {
        to_show[, estimate := coefficient]
    }
    
    ## Handle CI columns - m2dt() uses ci_lower/ci_upper (exponentiated) or coef_lower/coef_upper (raw)
    ## For forest plots, we need the raw (log scale) values for log_scale models
    ## and raw values for linear models
    if (!"conf_low" %in% names(to_show)) {
        if ("coef_lower" %in% names(to_show)) {
            ## Use raw coefficient bounds if available (these are on log scale for GLM/Cox)
            to_show[, conf_low := coef_lower]
            to_show[, conf_high := coef_upper]
        } else if ("ci_lower" %in% names(to_show)) {
            ## ci_lower/ci_upper are exponentiated for GLM/Cox, raw for LM
            ## Use the effect_type detected earlier
            if (effect_type %in% c("OR", "HR", "RR")) {
                ## Log scale model - need to log the CIs back
                to_show[, conf_low := log(ci_lower)]
                to_show[, conf_high := log(ci_upper)]
            } else {
                ## Linear model - CIs are already on raw scale
                to_show[, conf_low := ci_lower]
                to_show[, conf_high := ci_upper]
            }
        }
    }
    
    ## Handle level column - m2dt might use 'group' instead
    if ("group" %in% names(to_show) && !"level" %in% names(to_show)) {
        to_show[, level := group]
    }
    
    ## Handle n columns - m2dt uses n_group/events_group for categorical, n/events for continuous
    ## First try n_group, then fall back to n for continuous variables
    if (!"n_obs" %in% names(to_show)) {
        if ("n_group" %in% names(to_show) && "n" %in% names(to_show)) {
            ## Use n_group where available, fall back to n where n_group is NA
            to_show[, n_obs := data.table::fifelse(is.na(n_group), n, n_group)]
        } else if ("n_group" %in% names(to_show)) {
            to_show[, n_obs := n_group]
        } else if ("n" %in% names(to_show)) {
            to_show[, n_obs := n]
        }
    }
    if (!"n_events" %in% names(to_show)) {
        if ("events_group" %in% names(to_show) && "events" %in% names(to_show)) {
            ## Use events_group where available, fall back to events where events_group is NA
            to_show[, n_events := data.table::fifelse(is.na(events_group), events, events_group)]
        } else if ("events_group" %in% names(to_show)) {
            to_show[, n_events := events_group]
        } else if ("events" %in% names(to_show)) {
            to_show[, n_events := events]
        }
    }
    
    ## Handle is_reference - m2dt uses 'reference' column with text, not boolean
    if ("reference" %in% names(to_show) && !"is_reference" %in% names(to_show)) {
        to_show[, is_reference := (reference != "" & !is.na(reference))]
    }
    
    ## Keep only relevant columns that exist
    keep_cols <- c("predictor", "term", "level", "estimate", "conf_low", "conf_high", 
                   "p_value", "n_obs", "n_events", "is_reference", "effect_type")
    keep_cols <- keep_cols[keep_cols %in% names(to_show)]
    to_show <- to_show[, ..keep_cols]
    
    ## Ensure required columns exist even if empty
    if (!"level" %in% names(to_show)) {
        to_show[, level := NA_character_]
    }
    if (!"n_obs" %in% names(to_show)) {
        to_show[, n_obs := NA_integer_]
    }
    if (!"n_events" %in% names(to_show)) {
        to_show[, n_events := NA_integer_]
    }
    if (!"estimate" %in% names(to_show)) {
        to_show[, estimate := NA_real_]
    }
    if (!"conf_low" %in% names(to_show)) {
        to_show[, conf_low := NA_real_]
    }
    if (!"conf_high" %in% names(to_show)) {
        to_show[, conf_high := NA_real_]
    }
    if (!"is_reference" %in% names(to_show)) {
        to_show[, is_reference := FALSE]
    }
    
    ## Apply labels to predictor names
    to_show[, var_display := predictor]
    if (!is.null(labels)) {
        for (pred in unique(to_show$predictor)) {
            if (pred %in% names(labels)) {
                to_show[predictor == pred, var_display := labels[[pred]]]
            }
        }
    }
    
    ## Create display columns - handle case where level might be missing or NA
    to_show[, level_display := data.table::fifelse(
                                               is.na(level) | level == "" | level == "-", 
                                               "-", 
                                               as.character(level)
                                           )]
    
    ## Assign variable order for grouping and shading
    predictor_order <- unique(to_show$predictor)
    to_show[, var_order := match(predictor, predictor_order)]
    
    ## Calculate N and events formatting
    to_show[, N := n_obs]
    to_show[, n_formatted := data.table::fifelse(is.na(n_obs), "", format(n_obs, big.mark = ","))]
    to_show[, events_formatted := data.table::fifelse(is.na(n_events), "", format(n_events, big.mark = ","))]
    
    ## Format the values for display based on log_scale
    to_show_exp_clean <- data.table::copy(to_show)
    
    if (log_scale) {
        ## Log scale (OR, HR, RR): exponentiate values
        to_show_exp_clean[, hr := data.table::fifelse(is.na(estimate), NA_real_, exp(estimate))]
        to_show_exp_clean[, hr_formatted := data.table::fifelse(
                                                            is.na(estimate), "",
                                                            format_number(exp(estimate), digits)
                                                        )]
        to_show_exp_clean[, conf_low_formatted := data.table::fifelse(
                                                                  is.na(conf_low), "",
                                                                  format_number(exp(conf_low), digits)
                                                              )]
        to_show_exp_clean[, conf_high_formatted := data.table::fifelse(
                                                                   is.na(conf_high), "",
                                                                   format_number(exp(conf_high), digits)
                                                               )]
    } else {
        ## Linear scale: use raw coefficients
        to_show_exp_clean[, hr := data.table::fifelse(is.na(estimate), NA_real_, estimate)]
        to_show_exp_clean[, hr_formatted := data.table::fifelse(
                                                            is.na(estimate), "",
                                                            format_number(estimate, digits)
                                                        )]
        to_show_exp_clean[, conf_low_formatted := data.table::fifelse(
                                                                  is.na(conf_low), "",
                                                                  format_number(conf_low, digits)
                                                              )]
        to_show_exp_clean[, conf_high_formatted := data.table::fifelse(
                                                                   is.na(conf_high), "",
                                                                   format_number(conf_high, digits)
                                                               )]
    }
    
    ## Format p-values
    p_threshold <- 10^(-p_digits)
    p_threshold_str <- paste0("< ", format(p_threshold, scientific = FALSE))
    
    to_show_exp_clean[, p_formatted := data.table::fifelse(
                                                       is.na(p_value), "",
                                                       data.table::fifelse(p_value < p_threshold, p_threshold_str,
                                                                           format_number(p_value, p_digits))
                                                   )]
    
    ## Determine if ANY coefficient or CI bound is negative (only relevant for linear scale)
    ## If so, use comma notation throughout for consistency
    use_comma_notation <- !log_scale && any(to_show_exp_clean$conf_low < 0 | to_show_exp_clean$conf_high < 0, na.rm = TRUE)
    ci_separator <- if (use_comma_notation) ", " else "-"
    
    ## Create effect string with italic p
    to_show_exp_clean[, hr_string_expr := data.table::fifelse(
                                                          is.na(estimate) | is_reference == TRUE,
                                                          "''",
                                                          data.table::fcase(
                                                              p_value < p_threshold,
                                                              paste0("'", hr_formatted, " (", conf_low_formatted, ci_separator, conf_high_formatted, 
                                                                     "); '*~italic(p)~'", p_threshold_str, "'"),
                                                              
                                                              default = paste0("'", hr_formatted, " (", conf_low_formatted, ci_separator, conf_high_formatted, 
                                                                               "); '*~italic(p)~'= ", p_formatted, "'")
                                                          )
                                                      )]
    
    ## Handle reference rows
    to_show_exp_clean[is_reference == TRUE, hr_string_expr := "'reference'"]
    
    ## Handle condense_table: collapse binary variables to single rows
    if (condense_table) {
        indent_groups <- TRUE  # condense_table implies indent_groups
        
        ## Identify binary predictors (exactly 2 rows)
        predictor_counts <- to_show_exp_clean[, .N, by = predictor]
        binary_predictors <- predictor_counts[N == 2]$predictor
        
        if (length(binary_predictors) > 0) {
            ## Process each binary predictor
            rows_to_remove <- integer()
            
            for (pred in binary_predictors) {
                pred_rows <- to_show_exp_clean[predictor == pred]
                
                ## Find reference and non-reference rows
                ref_idx <- which(pred_rows$is_reference == TRUE)
                non_ref_idx <- which(pred_rows$is_reference == FALSE)
                
                if (length(ref_idx) == 1 && length(non_ref_idx) == 1) {
                    ref_category <- pred_rows$level[ref_idx]
                    non_ref_category <- pred_rows$level[non_ref_idx]
                    
                    ## Get label for this predictor
                    var_label <- pred_rows$var_display[1]
                    
                    ## Use greedy approach to determine if we should condense
                    if (should_condense_binary(ref_category, non_ref_category, var_label)) {
                        ## Keep the non-reference row but update display
                        to_show_exp_clean[predictor == pred & is_reference == FALSE, 
                                          var_display := var_label]
                    } else {
                        ## Append category name
                        to_show_exp_clean[predictor == pred & is_reference == FALSE, 
                                          var_display := paste0(var_label, " (", non_ref_category, ")")]
                    }
                    
                    ## Clear the level for condensed row
                    to_show_exp_clean[predictor == pred & is_reference == FALSE, level := "-"]
                    to_show_exp_clean[predictor == pred & is_reference == FALSE, level_display := "-"]
                    
                    ## Mark reference row for removal
                    ref_row_idx <- which(to_show_exp_clean$predictor == pred & 
                                         to_show_exp_clean$is_reference == TRUE)
                    rows_to_remove <- c(rows_to_remove, ref_row_idx)
                }
            }
            
            ## Remove reference rows of binary predictors
            if (length(rows_to_remove) > 0) {
                to_show_exp_clean <- to_show_exp_clean[-rows_to_remove]
            }
        }
    }
    
    ## Zebra stripe shading by variable
    if (zebra_stripes) {
        shade_colors <- c("0" = "#EEEEEE", "1" = "#FFFFFF")
        to_show_exp_clean[, shade_group := (var_order + 1) %% 2]
    } else {
        shade_colors <- c("0" = "#FFFFFF", "1" = "#FFFFFF")
        to_show_exp_clean[, shade_group := 0]
    }
    
    ## Handle indented display vs separate columns
    if (indent_groups || condense_table) {
        ## Create combined display with variable name as header, levels indented
        ## First, identify header rows (first occurrence of each predictor)
        to_show_exp_clean[, is_header := !duplicated(predictor)]
        
        ## Determine if each row is for a continuous variable (level is "-" or NA)
        to_show_exp_clean[, is_continuous := level_display == "-" | is.na(level_display)]
        
        ## For header rows of multi-level variables, show variable name (bold)
        ## For subsequent rows of multi-level variables, indent the level
        ## For continuous variables, show just the variable name
        to_show_exp_clean[, display_text := data.table::fifelse(
                                                            is_continuous,
                                                            var_display,  # Continuous variable: just name
                                                            data.table::fifelse(
                                                                            is_header,
                                                                            var_display,  # First row of categorical: variable name
                                                                            paste0("    ", level_display)  # Subsequent rows: indented level
                                                                        )
                                                        )]
        
        ## For multi-level variables, insert a header row if needed
        ## Check which predictors have multiple rows (categorical variables)
        multi_level_predictors <- to_show_exp_clean[, .N, by = predictor][N > 1]$predictor
        
        ## Create new data for header rows
        if (length(multi_level_predictors) > 0) {
            header_rows <- to_show_exp_clean[predictor %in% multi_level_predictors & is_header == TRUE]
            header_rows <- data.table::copy(header_rows)
            header_rows[, `:=`(
                display_text = var_display,
                estimate = NA_real_,
                conf_low = NA_real_,
                conf_high = NA_real_,
                p_value = NA_real_,
                hr = NA_real_,
                hr_string_expr = "''",
                n_formatted = "",
                events_formatted = "",
                is_header = TRUE,
                is_continuous = FALSE,
                row_type = "header"
            )]
            
            ## Mark data rows
            to_show_exp_clean[, row_type := "data"]
            to_show_exp_clean[predictor %in% multi_level_predictors, display_text := paste0("    ", level_display)]
            to_show_exp_clean[predictor %in% multi_level_predictors, is_continuous := FALSE]
            
            ## Combine and sort
            to_show_exp_clean <- rbind(header_rows, to_show_exp_clean, fill = TRUE)
            
            ## Sort by var_order and then by row_type (header first)
            to_show_exp_clean[, sort_key := data.table::fifelse(row_type == "header", 0, 1)]
            data.table::setorder(to_show_exp_clean, var_order, sort_key, level_display)
        } else {
            to_show_exp_clean[, row_type := "data"]
        }
        
        ## Determine which rows should be bold: headers and continuous variables
        to_show_exp_clean[, should_bold := (row_type == "header") | (is_continuous == TRUE)]
    } else {
        ## Non-indented display: variable and level in separate columns
        to_show_exp_clean[, display_text := var_display]
        to_show_exp_clean[duplicated(predictor), display_text := ""]
        to_show_exp_clean[, is_header := FALSE]
        to_show_exp_clean[, is_continuous := FALSE]
        to_show_exp_clean[, should_bold := TRUE]  # All variable names bold in non-indented mode
        to_show_exp_clean[, row_type := "data"]
    }
    
    ## Handle missing estimates for plotting
    to_show_exp_clean[is.na(estimate), estimate := 0]
    
    ## Reorder (flip) - maintain variable grouping
    to_show_exp_clean <- to_show_exp_clean[order(nrow(to_show_exp_clean):1)]
    to_show_exp_clean[, x_pos := .I]
    
    ## Calculate plot ranges
    if (log_scale) {
        rangeb <- range(to_show$conf_low, to_show$conf_high, na.rm = TRUE)
        min_ci <- min(to_show$conf_low, na.rm = TRUE)
        max_ci <- max(to_show$conf_high, na.rm = TRUE)
        
        ## Intelligent tick selection for log scale
        range_magnitude <- diff(rangeb)
        
        if (exp(min_ci) < 0.01 && exp(max_ci) > 2) {
            breaks <- c(0.01, 0.1, 0.5, 1, 2, 5)
        } else if (range_magnitude > 3) {
            all_breaks <- grDevices::axisTicks(rangeb/2, log = TRUE, nint = 7)
            if (length(all_breaks) > 7) {
                important <- c(1)
                if (min(all_breaks) < 0.5) important <- c(min(all_breaks), important)
                if (max(all_breaks) > 2) important <- c(important, max(all_breaks))
                other_breaks <- setdiff(all_breaks, important)
                if (length(other_breaks) > 3) {
                    keep_idx <- round(seq(1, length(other_breaks), length.out = 3))
                    other_breaks <- other_breaks[keep_idx]
                }
                breaks <- sort(unique(c(important, other_breaks)))
            } else {
                breaks <- all_breaks
            }
        } else {
            breaks <- grDevices::axisTicks(rangeb/2, log = TRUE, nint = 7)
        }
        
        if (!1 %in% breaks) breaks <- sort(unique(c(breaks, 1)))
        
        if (min_ci > 0) {
            rangeb[1] <- log(0.9)
        } else if (max_ci < 0) {
            rangeb[2] <- log(1.1)
        }
        
        breaks <- breaks[breaks >= exp(rangeb[1]) & breaks <= exp(rangeb[2])]
        reference_value <- 1
    } else {
        ## Linear scale
        rangeb <- range(to_show$conf_low, to_show$conf_high, na.rm = TRUE)
        range_width <- diff(rangeb)
        rangeb[1] <- rangeb[1] - range_width * 0.05
        rangeb[2] <- rangeb[2] + range_width * 0.05
        
        if (rangeb[1] > 0) rangeb[1] <- -0.1 * abs(rangeb[2])
        if (rangeb[2] < 0) rangeb[2] <- 0.1 * abs(rangeb[1])
        
        breaks <- pretty(rangeb, n = 7)
        breaks <- breaks[breaks >= rangeb[1] & breaks <= rangeb[2]]
        
        if (!0 %in% breaks) breaks <- sort(unique(c(breaks, 0)))
        
        reference_value <- 0
    }
    
    ## Calculate layout (note: condense_table implies indent_groups)
    layout <- calculate_uniforest_layout(
        to_show_exp_clean = to_show_exp_clean,
        show_n = show_n,
        show_events = show_events,
        indent_groups = (indent_groups || condense_table),
        table_width = table_width,
        center_padding = center_padding,
        effect_abbrev = effect_abbrev,
        font_size = font_size,
        log_scale = log_scale,
        rangeb = rangeb,
        ci_pct = ci_pct
    )

    ## Set up column positions
    y_variable <- layout$positions$variable
    if (!(indent_groups || condense_table)) {
        y_level <- layout$positions$level
    }
    if (show_n) {
        y_n <- layout$positions$n
    }
    if (show_events) {
        y_events <- layout$positions$events
    }
    y_effect <- layout$positions$effect

    rangeplot <- c(layout$rangeplot_start, rangeb[2] + diff(rangeb) * 0.05)

    ## Calculate recommended dimensions
    rec_height <- max(5, min(20, 3 + nrow(to_show_exp_clean) * 0.25))

    if (!is.null(plot_width)) {
        rec_width <- plot_width
        if (!is.null(plot_height)) {
            rec_height <- plot_height
        }
    } else {
        rec_width <- layout$total_width + 1.0
        rec_width <- max(10, min(20, rec_width))
    }

    ## Font sizes
    annot_font <- font_size * annot_size
    header_font <- font_size * header_size

    ## Custom ticks - for log scale, breaks are on exponentiated scale (0.5, 1, 2)
    ## For linear scale, breaks are on raw scale (-2, -1, 0, 1, 2)
    ticks_df <- data.frame(
        x = -0.5,
        xend = -0.7,
        y = breaks,  # Position on the plot (already transformed for log scale)
        yend = breaks,
        label = sprintf("%g", breaks)
    )

    ## Create filtered data for forest plot elements (exclude header rows)
    plot_data <- if ("row_type" %in% names(to_show_exp_clean)) {
                     to_show_exp_clean[row_type == "data" | is.na(row_type)]
                 } else {
                     to_show_exp_clean
                 }

    ## Define coordinate transformation
    tfm <- if (log_scale) exp else identity

    ## Create the plot
    p <- ggplot2::ggplot(to_show_exp_clean, ggplot2::aes(x_pos, tfm(estimate))) +
        
        ## Shading rectangles
        ggplot2::geom_rect(ggplot2::aes(xmin = x_pos - .5, xmax = x_pos + .5,
                                        ymin = tfm(rangeplot[1]), ymax = tfm(rangeplot[2]),
                                        fill = ordered(shade_group))) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.05))) +
        ggplot2::scale_size_area(max_size = 6, guide = "none") +
        ggplot2::scale_fill_manual(values = shade_colors, guide = "none") +
        
        ## Forest plot elements
        ggplot2::geom_point(data = plot_data[!is.na(hr) & is_reference != TRUE], 
                            ggplot2::aes(size = N), pch = 22, color = "#000000", fill = color, na.rm = TRUE) +
        ggplot2::geom_errorbar(data = plot_data[!is.na(hr) & is_reference != TRUE], 
                               ggplot2::aes(ymin = tfm(conf_low), ymax = tfm(conf_high)), 
                               width = 0.15, na.rm = TRUE) +
        
        ## Y-axis for forest plot
        ggplot2::annotate(geom = "segment",
                          x = -0.5, xend = -0.5,
                          y = tfm(rangeb[1]), yend = tfm(rangeb[2]),
                          color = "#000000", linewidth = 1) +
        
        ## Reference line
        ggplot2::annotate(geom = "segment",
                          x = -0.5, xend = max(to_show_exp_clean$x_pos) + 0.5,
                          y = reference_value, yend = reference_value,
                          linetype = "longdash") +
        
        ## Ticks
        ggplot2::geom_segment(data = ticks_df,
                              ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                              inherit.aes = FALSE, color = "#000000", linewidth = 1) +
        ggplot2::geom_text(data = ticks_df,
                           ggplot2::aes(x = xend - 0.05, y = y, label = label),
                           inherit.aes = FALSE, hjust = 0.5, vjust = 1.3,
                           size = annot_font * 1.5) +
        
        ## Coordinate system
        ggplot2::coord_flip(ylim = tfm(rangeplot)) +
        ggplot2::ggtitle(title) +
        
        ## Y-axis scale
        {if (log_scale) {
             ggplot2::scale_y_log10(name = effect_label,
                                    labels = sprintf("%g", breaks),
                                    expand = c(0.02, 0.02),
                                    breaks = breaks)
         } else {
             ggplot2::scale_y_continuous(name = effect_label,
                                         labels = sprintf("%g", breaks),
                                         expand = c(0.02, 0.02),
                                         breaks = breaks)
         }} +
        
        ggplot2::theme_light(base_family = detect_plot_font()) +
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
                       plot.title = ggplot2::element_text(size = font_size * title_size, 
                                                          face = "bold", hjust = 0.5)) +
        ggplot2::xlab("") +
        
        ## Variable column header
        ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = tfm(y_variable),
                          label = "Variable", fontface = "bold", hjust = 0, size = header_font) +
        
        ## Variable column content
        {if (indent_groups) {
             ## Use should_bold column, but respect bold_variables parameter
             fontfaces <- if (bold_variables) {
                              data.table::fifelse(to_show_exp_clean$should_bold, "bold", "plain")
                          } else {
                              rep("plain", nrow(to_show_exp_clean))
                          }
             ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = tfm(y_variable),
                               label = to_show_exp_clean$display_text,
                               fontface = fontfaces,
                               hjust = 0, size = annot_font)
         } else {
             ## Non-indented: bold all variable names if bold_variables is TRUE
             ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = tfm(y_variable),
                               label = to_show_exp_clean$display_text, 
                               fontface = if (bold_variables) "bold" else "plain",
                               hjust = 0, size = annot_font)
         }} +
        
        ## Level column (only when not indented or condensed)
        {if (!(indent_groups || condense_table)) {
             list(
                 ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = tfm(y_level),
                                   label = "Group", fontface = "bold", hjust = 0, size = header_font),
                 ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = tfm(y_level),
                                   label = to_show_exp_clean$level_display, hjust = 0, size = annot_font)
             )
         }} +
        
        ## N column
        {if (show_n) {
             list(
                 ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = tfm(y_n),
                                   label = "n", fontface = "bold.italic", hjust = 0.5, size = header_font),
                 ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = tfm(y_n),
                                   label = to_show_exp_clean$n_formatted, hjust = 0.5, size = annot_font)
             )
         }} +
        
        ## Events column
        {if (show_events) {
             list(
                 ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = tfm(y_events),
                                   label = "Events", fontface = "bold", hjust = 0.5, size = header_font),
                 ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = tfm(y_events),
                                   label = to_show_exp_clean$events_formatted, hjust = 0.5, size = annot_font)
             )
         }} +
        
        ## Effect column header (with italic p)
        ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = tfm(y_effect),
                          label = effect_header_expr, hjust = 0, size = header_font, parse = TRUE) +
        
        ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = tfm(y_effect),
                          label = to_show_exp_clean$hr_string_expr, hjust = 0,
                          size = annot_font, parse = TRUE) +
        
        ## X-axis label
        ggplot2::annotate(geom = "text", x = -1.5, y = reference_value,
                          label = effect_label, fontface = "bold",
                          hjust = 0.5, vjust = 2, size = annot_font * 1.5) +
        
        ## Outcome footer (conditional)
        {if (show_footer && !is.null(us_outcome)) {
             ## Apply label to outcome if available
             outcome_display <- us_outcome
             if (!is.null(labels)) {
                 ## For simple outcomes, check directly
                 if (us_outcome %in% names(labels)) {
                     outcome_display <- labels[[us_outcome]]
                 } else {
                     ## For Surv() outcomes, try to extract and label component variables
                     ## e.g., "Surv(os_months, os_status)" -> check for os_months, os_status labels
                     if (grepl("^Surv\\(", us_outcome)) {
                         ## Extract variable names from Surv()
                         surv_vars <- gsub("^Surv\\(|\\)$", "", us_outcome)
                         surv_vars <- trimws(strsplit(surv_vars, ",")[[1]])
                         
                         ## Apply labels to each component
                         labeled_vars <- sapply(surv_vars, function(v) {
                             if (v %in% names(labels)) labels[[v]] else v
                         })
                         outcome_display <- paste0("Surv(", paste(labeled_vars, collapse = ", "), ")")
                     }
                 }
             }
             outcome_text <- paste0("Outcome: ", outcome_display)
             ggplot2::annotate(geom = "text", x = 0.25, y = tfm(y_variable),
                               label = outcome_text,
                               size = annot_font, hjust = 0, vjust = 1.2, fontface = "italic")
         }}

    ## Convert units back for output if needed
    if (units != "in") {
        rec_width <- convert_units(rec_width, from = "in", to = units)
        rec_height <- convert_units(rec_height, from = "in", to = units)
    }

    ## Provide dimension recommendations
    if (is.null(plot_width) || is.null(plot_height)) {
        message(sprintf("Recommended plot dimensions: width = %.1f %s, height = %.1f %s",
                        rec_width, units, rec_height, units))
    }

    ## Add recommended dimensions as attribute
    attr(p, "recommended_dims") <- list(width = rec_width, height = rec_height)

    return(p)
}


#' Calculate table layout for uniforest plots
#' 
#' Internal function to determine column positions and widths for forest plot
#' table section. Positions are calculated in the same units as the data
#' (log scale for OR/HR/RR, linear for coefficients).
#' 
#' @param to_show_exp_clean Data.table with formatted data for plotting.
#' @param show_n Logical whether to include n column.
#' @param show_events Logical whether to include events column.
#' @param indent_groups Logical whether levels are indented.
#' @param table_width Proportion of width for table.
#' @param center_padding Padding between table and forest.
#' @param effect_abbrev Effect type abbreviation.
#' @param font_size Font size multiplier.
#' @param log_scale Logical whether using log scale.
#' @param rangeb Numeric vector with plot range bounds (in data units).
#' 
#' @return List with column positions, widths, and layout parameters.
#' 
#' @keywords internal
calculate_uniforest_layout <- function(to_show_exp_clean,
                                       show_n,
                                       show_events,
                                       indent_groups,
                                       table_width,
                                       center_padding,
                                       effect_abbrev,
                                       font_size,
                                       log_scale,
                                       rangeb,
                                       ci_pct = 95) {
    
    ## Character to inch conversion (approximate)
    char_to_inch <- 0.08 * font_size
    margin <- 0.15  # Increased margin for better separation between columns
    
    ## Calculate column widths in inches
    columns <- list()
    
    ## Variable column width
    var_width_chars <- max(nchar(as.character(to_show_exp_clean$display_text)), nchar("Variable"), na.rm = TRUE)
    columns$variable <- var_width_chars * char_to_inch
    
    ## Level column (only if not indented)
    if (!indent_groups) {
        level_width_chars <- max(nchar(as.character(to_show_exp_clean$level_display)), nchar("Group"), na.rm = TRUE)
        columns$level <- level_width_chars * char_to_inch
    }
    
    ## N column
    if (show_n) {
        n_width_chars <- max(nchar(to_show_exp_clean$n_formatted), nchar("n"), na.rm = TRUE)
        columns$n <- n_width_chars * char_to_inch
    }
    
    ## Events column
    if (show_events) {
        events_width_chars <- max(nchar(to_show_exp_clean$events_formatted), nchar("Events"), na.rm = TRUE)
        columns$events <- events_width_chars * char_to_inch
    }
    
    ## Effect column - use dynamic CI percentage
    effect_header <- paste0(effect_abbrev, " (", ci_pct, "% CI); p-value")
    effect_lengths <- nchar(paste0(
        to_show_exp_clean$hr_formatted, " (",
        to_show_exp_clean$conf_low_formatted, "-",
        to_show_exp_clean$conf_high_formatted, "); p = ",
        to_show_exp_clean$p_formatted
    ))
    effect_width_chars <- max(effect_lengths, nchar(effect_header), na.rm = TRUE) + center_padding
    columns$effect <- effect_width_chars * char_to_inch
    
    ## Calculate total table width in inches
    calc_table_width <- sum(unlist(columns)) + length(columns) * margin * 2
    
    ## Calculate forest width based on table_width proportion
    calc_forest_width <- calc_table_width * ((1 - table_width) / table_width)
    
    ## Convert table width to data scale units
    ## rangeb is in data units (log scale for OR/HR/RR, linear for coefficients)
    calc_range_width <- diff(rangeb)
    calc_table_width_data_units <- (calc_table_width / calc_forest_width) * calc_range_width
    
    ## Calculate positions in data scale units
    ## Start from the left edge (rangeb[1] minus table width)
    rangeplot_start <- rangeb[1] - calc_table_width_data_units
    
    ## Convert inch measurements to data scale units
    inch_to_data <- calc_range_width / calc_forest_width
    
    positions <- list()
    current_pos <- rangeplot_start
    
    for (name in names(columns)) {
        if (name == "events") {
            current_pos <- current_pos + margin * inch_to_data * 3.5
        } else {
            current_pos <- current_pos + margin * inch_to_data
        }
        positions[[name]] <- current_pos
        current_pos <- current_pos + columns[[name]] * inch_to_data + margin * inch_to_data
    }
    
    return(list(
        calc_table_width = calc_table_width,
        calc_forest_width = calc_forest_width,
        positions = positions,
        rangeplot_start = rangeplot_start,
        total_width = calc_table_width + calc_forest_width
    ))
}
