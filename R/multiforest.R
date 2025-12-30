#' Create Forest Plot for Multivariate Analysis
#'
#' Generates a publication-ready forest plot from a \code{\link{multifit}} output
#' object. The plot displays effect estimates (OR, HR, RR, or coefficients) with
#' confidence intervals across multiple outcomes, organized by outcome with the
#' predictor levels shown for each.
#'
#' @param x A multifit result object (data.table with class attributes from 
#'   \code{\link{multifit}}).
#'   
#' @param title Character string specifying the plot title. Default is 
#'   \code{"Multivariate Analysis"}. Use descriptive titles for publication.
#'   
#' @param effect_label Character string for the effect measure label on the 
#'   forest plot axis. Default is \code{NULL}, which auto-detects based on 
#'   model type (e.g., "Odds Ratio", "Hazard Ratio", "Rate Ratio", "Estimate").
#'   
#' @param column Character string specifying which column to plot when both 
#'   univariable and multivariable results are present. Options are 
#'   \code{"adjusted"} (default), \code{"unadjusted"}, or \code{"both"}.
#'   When \code{"both"}, shows side-by-side forest plots.
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
#'   number of events for each row. Default is \code{TRUE} for binomial and 
#'   survival models, \code{FALSE} for linear models.
#'
#' @param show_predictor Logical. If \code{TRUE}, includes the Predictor column 
#'   showing which level of a factor predictor is being compared. If \code{FALSE},
#'   omits the column (useful when predictor info is in the caption). Default is 
#'   \code{NULL}, which uses the \code{include_predictor} setting from multifit() 
#'   if available, otherwise \code{TRUE}.
#'
#' @param covariates_footer Logical. If \code{TRUE} (default), displays a footer
#'   listing the covariates used in adjusted models. Covariate names are formatted
#'   using the \code{labels} parameter if provided. Only shown when displaying
#'   adjusted results.
#'   
#' @param group_by_outcome Logical. If \code{TRUE} (default), groups rows by 
#'   outcome with outcome names as headers. If \code{FALSE}, shows flat list.
#'
#' @param bold_variables Logical. If \code{TRUE}, variable names are displayed
#'   in bold. If \code{FALSE} (default), variable names are displayed in plain
#'   text.
#'   
#' @param center_padding Numeric value specifying horizontal spacing between 
#'   table and forest plot. Default is 4.
#'   
#' @param zebra_stripes Logical. If \code{TRUE}, applies alternating gray 
#'   background shading to different outcomes for improved readability. 
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
#'   outcomes and variables. Applied to outcome names in the plot.
#'   Default is \code{NULL} (uses labels already applied in multifit, or 
#'   original names).
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
#' \strong{Plot Layout:}
#' 
#' The forest plot is organized with outcomes as grouping headers and predictor
#' levels (or interaction terms) as rows within each outcome. This provides a
#' clear visual comparison of how a single predictor affects multiple outcomes.
#' 
#' \enumerate{
#'   \item \strong{Title}: Centered at top
#'   \item \strong{Data Table} (left): Contains:
#'     \itemize{
#'       \item Outcome column (or grouped headers)
#'       \item Predictor/Group column
#'       \item n: Sample sizes (optional)
#'       \item Events: Event counts (optional, for applicable models)
#'       \item Effect (95\% CI); \emph{p}-value
#'     }
#'   \item \strong{Forest Plot} (right):
#'     \itemize{
#'       \item Point estimates (squares)
#'       \item 95\% confidence intervals
#'       \item Reference line at null value (1 or 0)
#'       \item Log scale for ratio measures
#'     }
#' }
#' 
#' \strong{Data Source:}
#' 
#' The function extracts effect estimates directly from the multifit output
#' object's \code{raw_data} attribute, which contains the numeric values
#' needed for plotting. This approach is efficient and ensures consistency
#' with the formatted table output.
#'
#' @seealso 
#' \code{\link{multifit}} for multivariate regression analysis,
#' \code{\link{glmforest}} for single GLM forest plots,
#' \code{\link{coxforest}} for single Cox model forest plots,
#' \code{\link{autoforest}} for automatic model detection
#'
#' @examples
#' # Load example data
#' data(clintrial)
#' data(clintrial_labels)
#' 
#' # Example 1: Basic multivariate forest plot
#' outcomes <- c("surgery", "pfs_status", "os_status")
#' result <- multifit(
#'     data = clintrial,
#'     outcomes = outcomes,
#'     predictor = "treatment",
#'     covariates = c("age", "sex", "stage")
#' )
#' 
#' plot1 <- multiforest(result)
#' print(plot1)
#' 
#' # Example 2: With custom title and labels
#' plot2 <- multiforest(
#'     result,
#'     title = "Treatment Effects Across Clinical Outcomes",
#'     labels = clintrial_labels
#' )
#' print(plot2)
#' 
#' # Example 3: Cox survival outcomes
#' surv_result <- multifit(
#'     data = clintrial,
#'     outcomes = c("Surv(pfs_months, pfs_status)", "Surv(os_months, os_status)"),
#'     predictor = "treatment",
#'     covariates = c("age", "sex", "stage"),
#'     model_type = "coxph"
#' )
#' 
#' plot3 <- multiforest(
#'     surv_result,
#'     title = "Treatment Effects on Survival Outcomes"
#' )
#' print(plot3)
#' 
#' # Example 4: Customize appearance
#' plot4 <- multiforest(
#'     result,
#'     color = "#E74C3C",
#'     zebra_stripes = TRUE,
#'     show_n = FALSE,
#'     font_size = 1.2
#' )
#' print(plot4)
#' 
#' # Example 5: Save with recommended dimensions
#' plot5 <- multiforest(result)
#' dims <- attr(plot5, "recommended_dims")
#' # ggsave("multioutcome_forest.pdf", plot5,
#' #        width = dims$width, height = dims$height)
#'
#' @export
multiforest <- function(x,
                        title = "Multivariate Analysis",
                        effect_label = NULL,
                        column = "adjusted",
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
                        show_predictor = NULL,
                        covariates_footer = TRUE,
                        indent_predictor = FALSE,
                        bold_variables = TRUE,
                        center_padding = 4,
                        zebra_stripes = TRUE,
                        color = NULL,
                        null_line = NULL,
                        log_scale = NULL,
                        labels = NULL,
                        units = "in") {
    
    ## Input validation
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required but not installed.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required but not installed.")
    }
    
    ## Check if user passed a model object instead of a multifit result
    ## Common mistake since other forest functions (coxforest, glmforest) take models
    if (inherits(x, c("coxph", "glm", "lm", "glmerMod", "lmerMod", "coxme"))) {
        stop("multiforest() expects a multifit() result, not a model object.\n",
             "\nFor single-model forest plots, use:\n",
             "  - coxforest() for Cox models\n",
             "  - glmforest() for GLM models\n",
             "  - lmforest() for linear models\n",
             "\nFor multi-outcome analysis, first run multifit(), then pass the result to multiforest().")
    }
    
    if (!data.table::is.data.table(x)) {
        stop("multiforest() expects a multifit() result (data.table with raw_data attribute).\n",
             "\nInput is of class: ", paste(class(x), collapse = ", "))
    }
    
    raw_data <- attr(x, "raw_data")
    if (is.null(raw_data)) {
        stop("Input does not appear to be a multifit() result (missing raw_data attribute).\n",
             "\nmultiforest() requires output from multifit(). For single-model forest plots, use:\n",
             "  - coxforest() for Cox models\n",
             "  - glmforest() for GLM models\n",
             "  - lmforest() for linear models")
    }
    
    ## Extract multifit attributes
    mf_predictor <- attr(x, "predictor")
    mf_outcomes <- attr(x, "outcomes")
    mf_model_type <- attr(x, "model_type")
    mf_columns <- attr(x, "columns")
    
    ## Resolve show_predictor: use multifit attribute if not specified
    if (is.null(show_predictor)) {
        mf_include_predictor <- attr(x, "include_predictor")
        show_predictor <- if (!is.null(mf_include_predictor)) mf_include_predictor else TRUE
    }
    
    ## Validate column parameter
    column <- match.arg(column, c("adjusted", "unadjusted", "both"))
    
    ## If multifit only has one column type, use that
    if (mf_columns %in% c("adjusted", "unadjusted")) {
        column <- mf_columns
    }
    
    ## Determine effect type and defaults
    effect_type <- raw_data$effect_type[1]
    if (is.na(effect_type)) effect_type <- "Estimate"
    
    if (is.null(null_line)) {
        null_line <- if (effect_type %in% c("OR", "HR", "RR")) 1 else 0
    }
    
    if (is.null(log_scale)) {
        log_scale <- effect_type %in% c("OR", "HR", "RR")
    }
    
    ## Determine adjustment prefix based on column parameter
    ## "adjusted" -> "a" prefix for OR/HR/RR, "Adj." for Coefficient
    ## "unadjusted" -> no prefix
    is_adjusted <- (column == "adjusted") || (mf_columns == "adjusted")
    
    if (effect_type %in% c("OR", "HR", "RR")) {
        adj_prefix <- if (is_adjusted) "a" else ""
        effect_abbrev <- paste0(adj_prefix, effect_type)  # e.g., "aOR", "aHR", "HR", "OR"
    } else {
        adj_prefix <- if (is_adjusted) "Adj. " else ""
        effect_abbrev <- paste0(adj_prefix, "Coefficient")  # e.g., "Adj. Coefficient" or "Coefficient"
    }
    
    ## Set effect_label for axis (full name)
    if (is.null(effect_label)) {
        effect_label <- switch(effect_type,
                               "OR" = "Odds Ratio",
                               "HR" = "Hazard Ratio",
                               "RR" = "Rate Ratio",
                               "Coefficient")
    }
    
    ## Format CI percentage for display in headers
    ci_pct <- round(conf_level * 100)
    
    ## Create header label for the effect column with italic p (e.g., "aOR (95% CI); p-value")
    ## Format: bold('aOR (95% CI); '*bolditalic(p)*'-value')
    effect_header_expr <- paste0("bold('", effect_abbrev, " (", ci_pct, "% CI); '*bolditalic(p)*'-value')")
    
    if (is.null(show_events)) {
        show_events <- effect_type %in% c("OR", "HR", "RR")
    }
    
    ## Set default color based on effect type to match other forest plot conventions
    if (is.null(color)) {
        color <- switch(effect_type,
                        "HR" = "#8A61D8",
                        "OR" = "#3C8D9C",
                        "RR" = "#3064A6",
                        "#5A8F5A")  # Default for Estimate/Coefficient
    }
    
    ## Prepare data - handle "both" columns case
    dt <- data.table::copy(raw_data)
    
    ## For "both" column multifit results, select the appropriate columns
    if ("exp_coef_adj" %in% names(dt)) {
        if (column == "adjusted") {
            dt[, `:=`(
                exp_coef = exp_coef_adj,
                ci_lower = ci_lower_adj,
                ci_upper = ci_upper_adj,
                p_value = p_adj
            )]
        } else {
            dt[, `:=`(
                exp_coef = exp_coef_unadj,
                ci_lower = ci_lower_unadj,
                ci_upper = ci_upper_unadj,
                p_value = p_unadj
            )]
        }
    }
    
    ## Apply labels to outcomes if provided
    if (!is.null(labels)) {
        for (orig_name in names(labels)) {
            dt[outcome == orig_name, outcome := labels[[orig_name]]]
        }
    }
    
    ## Convert to log scale using conf_level for CI calculation
    if ("coefficient" %in% names(dt)) {
        dt[, estimate := coefficient]
        if ("se" %in% names(dt)) {
            z_crit <- qnorm((1 + conf_level) / 2)
            dt[, conf_low := coefficient - z_crit * se]
            dt[, conf_high := coefficient + z_crit * se]
        } else {
            dt[, conf_low := log(ci_lower)]
            dt[, conf_high := log(ci_upper)]
        }
    } else {
        dt[, estimate := log(exp_coef)]
        dt[, conf_low := log(ci_lower)]
        dt[, conf_high := log(ci_upper)]
    }
    
    ## Rename columns to match glmforest format
    ## var = outcome (the "variable" in multifit context is the outcome)
    dt[, var := outcome]
    
    ## level = group (predictor level like "Drug A", "Drug B")
    dt[, level := data.table::fifelse(is.na(group) | group == "-", "-", group)]
    
    ## N and Events (uppercase to match glmforest)
    if ("n" %in% names(dt)) {
        data.table::setnames(dt, "n", "N")
    }
    if ("events" %in% names(dt)) {
        data.table::setnames(dt, "events", "Events")
    }
    
    ## Assign row positions (pos, var_order)
    outcomes <- unique(dt$var)
    
    rows_list <- list()
    row_counter <- 0
    
    for (i in seq_along(outcomes)) {
        outcome_name <- outcomes[i]
        outcome_rows <- dt[var == outcome_name]
        
        for (j in seq_len(nrow(outcome_rows))) {
            row_counter <- row_counter + 1
            data_row <- data.table::copy(outcome_rows[j])
            data_row[, `:=`(
                pos = j,
                var_order = i
            )]
            rows_list[[row_counter]] <- data_row
        }
    }
    
    dt <- data.table::rbindlist(rows_list, fill = TRUE)
    
    ## Add zebra stripe shading
    if (zebra_stripes) {
        ## Use (var_order + 1) %% 2 so first group (var_order=1) gets shade_group=0 (gray)
        dt[, shade_group := (var_order + 1) %% 2]
        shade_colors <- c("#EEEEEE", "#FFFFFF")  # 0 = gray, 1 = white
    } else {
        dt[, shade_group := 0]
        shade_colors <- c("#FFFFFF", "#FFFFFF")
    }
    
    ## Select final columns and sort
    cols_to_keep <- c("var", "level", "N", "Events", "p_value", "estimate", "conf_low", "conf_high", "pos", "var_order", "shade_group")
    cols_available <- intersect(cols_to_keep, names(dt))
    
    to_show <- dt[, ..cols_available]
    
    data.table::setorder(to_show, var_order, pos)
    
    ## Format the values for display
    to_show_exp_clean <- data.table::copy(to_show)
    
    ## Create formatted columns for display
    ## For log_scale = TRUE (OR, HR, RR): exponentiate values
    ## For log_scale = FALSE (linear models): use raw coefficients
    if (log_scale) {
        to_show_exp_clean[, hr := data.table::fifelse(is.na(estimate), 
                                                      NA_real_,
                                                      exp(estimate))]
        
        ## Format effect estimate (exponentiated)
        to_show_exp_clean[, hr_formatted := data.table::fifelse(is.na(estimate),
                                                                "",
                                                                format(round(exp(estimate), digits), nsmall = digits))]
        
        to_show_exp_clean[, conf_low_formatted := data.table::fifelse(is.na(conf_low), 
                                                                      "",
                                                                      format(round(exp(conf_low), digits), nsmall = digits))]
        to_show_exp_clean[, conf_high_formatted := data.table::fifelse(is.na(conf_high), 
                                                                       "",
                                                                       format(round(exp(conf_high), digits), nsmall = digits))]
    } else {
        ## Linear scale - use raw coefficients
        to_show_exp_clean[, hr := data.table::fifelse(is.na(estimate), 
                                                      NA_real_,
                                                      estimate)]
        
        ## Format effect estimate (raw coefficient)
        to_show_exp_clean[, hr_formatted := data.table::fifelse(is.na(estimate),
                                                                "",
                                                                format(round(estimate, digits), nsmall = digits))]
        
        to_show_exp_clean[, conf_low_formatted := data.table::fifelse(is.na(conf_low), 
                                                                      "",
                                                                      format(round(conf_low, digits), nsmall = digits))]
        to_show_exp_clean[, conf_high_formatted := data.table::fifelse(is.na(conf_high), 
                                                                       "",
                                                                       format(round(conf_high, digits), nsmall = digits))]
    }

    ## Format p-values using p_digits parameter
    p_threshold <- 10^(-p_digits)
    p_threshold_str <- paste0("< ", format(p_threshold, scientific = FALSE))
    
    to_show_exp_clean[, p_formatted := data.table::fifelse(is.na(p_value), 
                                                           "",
                                                           data.table::fifelse(p_value < p_threshold, 
                                                                               p_threshold_str,
                                                                               format(round(p_value, p_digits), nsmall = p_digits)))]

    ## Create the combined effect string with expression for italic p
    ## Use comma separator when CI bounds are negative (for linear models)
    to_show_exp_clean[, hr_string_expr := data.table::fifelse(
                                                          is.na(estimate),
                                                          "''",  # Empty string for header rows
                                                          data.table::fcase(
                                                              p_value < p_threshold & (conf_low < 0 | conf_high < 0),
                                                              paste0("'", hr_formatted, " (", conf_low_formatted, ", ", 
                                                                     conf_high_formatted, "); '*~italic(p)~'", p_threshold_str, "'"),
                                                              
                                                              p_value < p_threshold,
                                                              paste0("'", hr_formatted, " (", conf_low_formatted, "-", 
                                                                     conf_high_formatted, "); '*~italic(p)~'", p_threshold_str, "'"),
                                                              
                                                              conf_low < 0 | conf_high < 0,
                                                              paste0("'", hr_formatted, " (", conf_low_formatted, ", ", 
                                                                     conf_high_formatted, "); '*~italic(p)~'= ", p_formatted, "'"),
                                                              
                                                              default = paste0("'", hr_formatted, " (", conf_low_formatted, "-", 
                                                                               conf_high_formatted, "); '*~italic(p)~'= ", p_formatted, "'")
                                                          )
                                                      )]

    ## Format N and events with thousands separator
    to_show_exp_clean[, n_formatted := data.table::fifelse(is.na(N), "", format(N, big.mark = ",", scientific = FALSE))]
    to_show_exp_clean[, events_formatted := data.table::fifelse(is.na(Events), "", format(as.integer(Events), big.mark = ",", scientific = FALSE))]
    
    ## Clean up variable names for display
    to_show_exp_clean[, var_display := as.character(var)]

    if (indent_predictor) {
        ## Apply labels to outcome names
        for (v in unique(to_show_exp_clean$var)) {
            if (v != "" && !grepl("^    ", v)) {
                if (!is.null(labels) && v %in% names(labels)) {
                    to_show_exp_clean[var == v, var_display := labels[v]]
                }
            }
        }
        
        ## Create indented display: outcome header rows + indented predictor rows
        outcomes <- unique(to_show_exp_clean$var)
        new_rows <- list()
        row_idx <- 0
        
        for (i in seq_along(outcomes)) {
            outcome_name <- outcomes[i]
            outcome_rows <- to_show_exp_clean[var == outcome_name]
            outcome_label <- outcome_rows$var_display[1]
            
            ## shade_group: first outcome (i = 1) should be 0 (gray)
            current_shade <- (i + 1) %% 2
            
            ## Add header row for outcome (bold, no estimate)
            row_idx <- row_idx + 1
            header_row <- data.table::data.table(
                                          var = outcome_name,
                                          level = "",
                                          N = NA_integer_,
                                          Events = NA_real_,
                                          p_value = NA_real_,
                                          estimate = NA_real_,
                                          conf_low = NA_real_,
                                          conf_high = NA_real_,
                                          pos = 0L,
                                          var_order = as.integer(i),
                                          shade_group = current_shade,
                                          var_display = outcome_label,
                                          is_header = TRUE,
                                          hr_formatted = "",
                                          conf_low_formatted = "",
                                          conf_high_formatted = "",
                                          p_formatted = "",
                                          hr_string_expr = "''",
                                          n_formatted = "",
                                          events_formatted = ""
                                      )
            new_rows[[row_idx]] <- header_row
            
            ## Add indented predictor level rows
            for (j in seq_len(nrow(outcome_rows))) {
                row_idx <- row_idx + 1
                data_row <- data.table::copy(outcome_rows[j])
                data_row[, `:=`(
                    var_display = paste0("    ", level),  # Indent predictor level
                    is_header = FALSE,
                    shade_group = current_shade  # Match header row shading
                )]
                new_rows[[row_idx]] <- data_row
            }
        }
        
        to_show_exp_clean <- data.table::rbindlist(new_rows, fill = TRUE)
        to_show_exp_clean[, level := ""]  # Clear level column since it's now in var_display
        
    } else {
        ## Without indent_predictor: show outcome and predictor in separate columns
        ## Apply labels to outcome names
        for (v in unique(to_show_exp_clean$var)) {
            if (v %in% to_show_exp_clean$var) {
                if (!is.null(labels) && v %in% names(labels)) {
                    to_show_exp_clean[var == v, var_display := labels[v]]
                }
            }
        }
        
        ## For non-indented display, avoid repeating outcome name for same outcome
        to_show_exp_clean[duplicated(var), var_display := ""]
        to_show_exp_clean[, is_header := FALSE]
    }
    
    ## Handle missing estimates for plotting
    to_show_exp_clean[is.na(estimate), estimate := 0]
    
    ## Reorder (flip) - but maintain the variable grouping
    to_show_exp_clean <- to_show_exp_clean[order(nrow(to_show_exp_clean):1)]
    to_show_exp_clean[, x_pos := .I]
    
    ## Calculate plot ranges - different for log vs linear scale
    if (log_scale) {
        ## Log scale (OR, HR, RR)
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
        
        if (!1 %in% breaks) {
            breaks <- sort(unique(c(breaks, 1)))
        }
        
        if (min_ci > 0) {
            rangeb[1] <- log(0.9)
        } else if(max_ci < 0) {
            rangeb[2] <- log(1.1)
        }
        
        breaks <- breaks[breaks >= exp(rangeb[1]) & breaks <= exp(rangeb[2])]
        reference_value <- 1
    } else {
        ## Linear scale (coefficients)
        rangeb <- range(to_show$conf_low, to_show$conf_high, na.rm = TRUE)
        
        ## Add some padding
        range_width <- diff(rangeb)
        rangeb[1] <- rangeb[1] - range_width * 0.05
        rangeb[2] <- rangeb[2] + range_width * 0.05
        
        ## Ensure 0 is included in the range for reference line
        if (rangeb[1] > 0) rangeb[1] <- -0.1 * abs(rangeb[2])
        if (rangeb[2] < 0) rangeb[2] <- 0.1 * abs(rangeb[1])
        
        ## Generate breaks for linear scale
        breaks <- pretty(rangeb, n = 7)
        breaks <- breaks[breaks >= rangeb[1] & breaks <= rangeb[2]]
        
        if (!0 %in% breaks) {
            breaks <- sort(unique(c(breaks, 0)))
        }
        
        reference_value <- 0
    }

    ## Calculate layout
    layout <- calculate_multiforest_layout(
        to_show_exp_clean = to_show_exp_clean,
        show_n = show_n,
        show_events = show_events,
        indent_predictor = indent_predictor,
        show_predictor = show_predictor,
        effect_abbrev = effect_abbrev,
        font_size = font_size,
        table_width = data.table::fifelse(is.null(table_width), 0.6, table_width),
        rangeb = rangeb,
        center_padding = center_padding,
        ci_pct = ci_pct
    )

    ## Set up the extended range for plotting
    rangeplot <- c(layout$rangeplot_start, rangeb[2] + diff(rangeb) * 0.05)

    ## Extract positions
    y_outcome <- layout$positions$outcome
    if (!indent_predictor && show_predictor) {
        y_predictor <- layout$positions$predictor
    }
    if (show_n) {
        y_n <- layout$positions$n
    }
    if (show_events) {
        y_events <- layout$positions$events
    }
    y_effect <- layout$positions$effect

    ## Calculate recommended dimensions
    rec_height <- max(5, min(20, 3 + nrow(to_show_exp_clean) * 0.25))

    if(!is.null(plot_width)) {
        rec_width <- plot_width
        if(!is.null(plot_height)) {
            rec_height <- plot_height
        }
    } else {
        rec_width <- layout$total_width + 1.0
        rec_width <- max(10, min(20, rec_width))
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

    ## Create filtered data for forest plot elements (exclude header rows)
    plot_data <- if ("is_header" %in% names(to_show_exp_clean)) {
                     to_show_exp_clean[is_header == FALSE | is.na(is_header)]
                 } else {
                     to_show_exp_clean
                 }

    ## Define coordinate transformation based on scale type
    ## For log scale: use exp() to transform from log to original scale
    ## For linear scale: use identity (no transformation)
    tfm <- if (log_scale) exp else identity

    ## Create the plot
    p <- ggplot2::ggplot(to_show_exp_clean, ggplot2::aes(x_pos, tfm(estimate))) +
        
        ## Shading rectangles
        ggplot2::geom_rect(ggplot2::aes(xmin = x_pos - .5, xmax = x_pos + .5,
                                        ymin = tfm(rangeplot[1]), ymax = tfm(rangeplot[2]),
                                        fill = ordered(shade_group))) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.05))) +
        ## Size scaling: use sqrt transformation with reasonable limits
        ## This ensures points are visible even when N values are similar
        ggplot2::scale_size_area(max_size = 6, guide = "none") +
        ggplot2::scale_fill_manual(values = shade_colors, guide = "none") +
        
        ## Forest plot elements (use filtered data without header rows)
        ggplot2::geom_point(data = plot_data, ggplot2::aes(size = N), pch = 22, color = "#000000", fill = color, na.rm = TRUE) +
        ggplot2::geom_errorbar(data = plot_data, ggplot2::aes(ymin = tfm(conf_low), ymax = tfm(conf_high)), width = 0.15, na.rm = TRUE) +
        
        ## Y-axis for forest plot
        ggplot2::annotate(geom = "segment",
                          x = -0.5, xend = -0.5,
                          y = tfm(rangeb[1]),
                          yend = tfm(rangeb[2]),
                          color = "#000000", linewidth = 1) +
        
        ## Reference line (1 for log scale, 0 for linear)
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
        
        ## Set coordinate system with extended limits
        ggplot2::coord_flip(ylim = tfm(rangeplot)) +
        ggplot2::ggtitle(title) +
        
        ## Y-axis scale (log10 or continuous)
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
        
        ## Outcome column header
        ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = tfm(y_outcome),
                          label = "Outcome", fontface = "bold", hjust = 0,
                          size = header_font) +

    {if (indent_predictor) {
         ## When indented, use conditional formatting for headers vs data rows
         fontfaces <- if (bold_variables) {
                          data.table::fifelse(to_show_exp_clean$is_header, "bold", "plain")
                      } else {
                          rep("plain", nrow(to_show_exp_clean))
                      }
         ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = tfm(y_outcome),
                           label = to_show_exp_clean$var_display, 
                           fontface = fontfaces, 
                           hjust = 0,
                           size = annot_font)
     } else {
         ## Non-indented: outcome names are bold if bold_variables is TRUE
         ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = tfm(y_outcome),
                           label = to_show_exp_clean$var_display, 
                           fontface = if (bold_variables) "bold" else "plain", 
                           hjust = 0,
                           size = annot_font)
     }} +
    
    ## Predictor level column (only shown when not indented AND show_predictor is TRUE)
    {if (!indent_predictor && show_predictor) {
         list(
             ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = tfm(y_predictor),
                               label = "Predictor", fontface = "bold", hjust = 0,
                               size = header_font),
             ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = tfm(y_predictor),
                               label = to_show_exp_clean$level, hjust = 0,
                               size = annot_font)
         )
     }} +
    
    ## N column (conditional)
    {if (show_n) {
         list(
             ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = tfm(y_n),
                               label = "n", fontface = "bold.italic", hjust = 0.5,
                               size = header_font),
             ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = tfm(y_n),
                               label = to_show_exp_clean$n_formatted, hjust = 0.5,
                               size = annot_font)
         )
     }} +
    
    ## Events column (conditional)
    {if (show_events) {
         list(
             ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = tfm(y_events),
                               label = "Events", fontface = "bold", hjust = 0.5,
                               size = header_font),
             ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = tfm(y_events),
                               label = to_show_exp_clean$events_formatted, hjust = 0.5,
                               size = annot_font)
         )
     }} +
    
    ## Effect column header (with italic p)
    ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = tfm(y_effect),
                      label = effect_header_expr,
                      hjust = 0, size = header_font, parse = TRUE) +
    
    ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = tfm(y_effect),
                      label = to_show_exp_clean$hr_string_expr, hjust = 0,
                      size = annot_font, parse = TRUE) +
    
    ## X-axis label
    ggplot2::annotate(geom = "text", x = -1.5, y = reference_value,
                      label = effect_label, fontface = "bold",
                      hjust = 0.5, vjust = 2, size = annot_font * 1.5) +
    
    ## Covariates footer (if adjusted model with covariates and covariates_footer = TRUE)
    {
        mf_covariates <- attr(x, "covariates")
        if (covariates_footer && !is.null(mf_covariates) && length(mf_covariates) > 0 && is_adjusted) {
            ## Apply labels to covariate names if available
            covar_display <- sapply(mf_covariates, function(cv) {
                if (!is.null(labels) && cv %in% names(labels)) {
                    labels[[cv]]
                } else {
                    cv
                }
            })
            
            ## Format covariates list, truncating if very long
            if (length(covar_display) > 8) {
                covar_text <- paste0("Covariates: ", 
                                     paste(covar_display[1:8], collapse = ", "),
                                     ", ... (", length(covar_display), " total)")
            } else {
                covar_text <- paste0("Covariates: ", paste(covar_display, collapse = ", "))
            }
            ggplot2::annotate(geom = "text", x = 0.25, y = tfm(y_outcome),
                              label = covar_text,
                              size = annot_font, hjust = 0, vjust = 1.2, fontface = "italic")
        }
    }
    
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
    attr(p, "recommended_dims") <- list(width = rec_width, height = rec_height)

    return(p)
}

#' Calculate table layout for multiforest plots
#' 
#' Computes column widths and positions for the table portion of a multiforest plot.
#' 
#' @keywords internal
calculate_multiforest_layout <- function(to_show_exp_clean, show_n, show_events, 
                                         indent_predictor, show_predictor = TRUE,
                                         effect_abbrev,
                                         font_size, table_width = 0.6, rangeb, 
                                         center_padding, ci_pct = 95) {
    
    char_to_inch <- 0.08 * font_size
    margin <- 0.3  # inches between columns
    
    ## Build list of active columns with their widths
    columns <- list()
    
    ## Always have outcome column
    outcome_width_chars <- max(nchar(to_show_exp_clean$var_display), nchar("Outcome"), na.rm = TRUE)
    columns$outcome <- outcome_width_chars * char_to_inch
    
    ## Conditionally add predictor column (only when not indented AND show_predictor is TRUE)
    if (!indent_predictor && show_predictor) {
        predictor_width_chars <- max(nchar(to_show_exp_clean$level), nchar("Predictor"), na.rm = TRUE)
        columns$predictor <- predictor_width_chars * char_to_inch + 0.2
    }
    
    ## Conditionally add n column
    if (show_n) {
        n_width_chars <- max(nchar(to_show_exp_clean$n_formatted), nchar("____"), na.rm = TRUE)
        columns$n <- n_width_chars * char_to_inch
    }
    
    ## Conditionally add events column
    if (show_events) {
        events_width_chars <- max(nchar(to_show_exp_clean$events_formatted), nchar("Events"), na.rm = TRUE)
        columns$events <- events_width_chars * char_to_inch
    }
    
    ## Always have effect column - use dynamic CI percentage
    effect_header <- paste0(effect_abbrev, " (", ci_pct, "% CI); p-value")
    
    ## Calculate max effect string length
    effect_lengths <- nchar(paste0(
        to_show_exp_clean$hr_formatted, " (",
        to_show_exp_clean$conf_low_formatted, "-",
        to_show_exp_clean$conf_high_formatted, "); p = ",
        to_show_exp_clean$p_formatted
    ))
    
    effect_width_chars <- max(effect_lengths, nchar(effect_header), na.rm = TRUE) + center_padding
    columns$effect <- effect_width_chars * char_to_inch
    
    ## Calculate total table width
    calc_table_width <- sum(unlist(columns)) + length(columns) * margin * 2
    
    ## Calculate forest width based on table_width proportion
    calc_forest_width <- calc_table_width * ((1 - table_width) / table_width)
    
    ## Convert table width to log scale units
    calc_range_width <- diff(rangeb)
    calc_table_width_log_units <- (calc_table_width / calc_forest_width) * calc_range_width
    
    ## Calculate positions in log scale
    rangeplot_start <- rangeb[1] - calc_table_width_log_units
    
    ## Convert inch measurements to log scale units
    inch_to_log <- calc_range_width / calc_forest_width
    
    positions <- list()
    current_pos <- rangeplot_start
    
    for (name in names(columns)) {
        current_pos <- current_pos + margin * inch_to_log
        positions[[name]] <- current_pos
        current_pos <- current_pos + columns[[name]] * inch_to_log + margin * inch_to_log
    }
    
    return(list(
        calc_table_width = calc_table_width,
        calc_forest_width = calc_forest_width,
        positions = positions,
        rangeplot_start = rangeplot_start,
        total_width = calc_table_width + calc_forest_width
    ))
}
