#' Convert between units
#' 
#' Converts measurements between different unit systems commonly used in
#' graphics (inches, centimeters, millimeters, pixels, points).
#' 
#' @param value Numeric value to convert.
#' @param from Character string specifying source unit ("in", "cm", "mm", "px", "pt").
#' @param to Character string specifying target unit ("in", "cm", "mm", "px", "pt").
#' @param dpi Integer dots per inch for pixel conversions (default 96).
#' @return Numeric value in target units.
#' @keywords internal
convert_units <- function(value, from = "in", to = "in", dpi = 96) {
    to_inches <- c(
        "in" = 1,
        "cm" = 1/2.54,    # 1 cm = 1/2.54 inches
        "mm" = 1/25.4,    # 1 mm = 1/25.4 inches
        "px" = 1/dpi,     # 1 px = 1/dpi inches
        "pt" = 1/72       # 1 pt = 1/72 inches
    )
    
    value_in_inches <- value * to_inches[[from]]
    value_in_target <- value_in_inches / to_inches[[to]]
    
    return(value_in_target)
}

#' Calculate table layout for forest plots
#' 
#' Computes column widths and positions for the table portion of a forest plot.
#' Determines spacing based on content width, font size, and desired table/forest
#' proportion. Returns positions in log-scale units for plot coordinate system.
#' 
#' @param to_show_exp_clean Data.table with formatted display data for the plot.
#' @param show_n Logical whether to include sample size column.
#' @param show_events Logical whether to include events column.
#' @param indent_groups Logical whether groups are indented (affects level column).
#' @param condense_table Logical whether table is condensed (affects level column).
#' @param effect_label Character string describing effect measure type.
#' @param ref_label Character string label for reference categories.
#' @param font_size Numeric font size for width calculations.
#' @param table_width Numeric proportion of total width for table (0-1).
#' @param rangeb Numeric vector of length 2 with plot x-axis range.
#' @param center_padding Numeric additional padding for effect column.
#' @return List with table_width, forest_width, positions, rangeplot_start,
#'   total_width, and effect_abbrev components.
#' @keywords internal
calculate_forest_layout <- function(to_show_exp_clean, show_n, show_events, 
                                    indent_groups, condense_table, 
                                    effect_label, ref_label, font_size, 
                                    table_width = 0.6, rangeb, center_padding) {
    
    char_to_inch <- 0.08 * font_size
    margin <- 0.3  # inches between columns
    
    ## Build list of active columns with their widths
    columns <- list()
    
    ## Always have variable column
    var_width_chars <- max(nchar(to_show_exp_clean$var_display), nchar("Variable"), na.rm = TRUE)
    columns$var <- var_width_chars * char_to_inch
    
    ## Conditionally add level column
    if (!(indent_groups || condense_table)) {
        level_width_chars <- max(nchar(to_show_exp_clean$level), nchar("Group"), na.rm = TRUE)
        columns$level <- level_width_chars * char_to_inch
    }
    
    ## Conditionally add n column
    if (show_n) {
        n_width_chars <- max(nchar(to_show_exp_clean$n_formatted), nchar("____"), na.rm = TRUE)
        columns$n <- n_width_chars * char_to_inch
    }
    
    ## Conditionally add events column
    if (show_events) {
        events_width_chars <- max(nchar(to_show_exp_clean$events_formatted), nchar("___"), na.rm = TRUE)
        columns$events <- events_width_chars * char_to_inch
    }
    
    ## Always have effect column
    ## Create abbreviation from effect_label - handle common cases and custom labels
    ## For multivariable models, use adjusted abbreviations
    effect_abbrev <- if(effect_label == "Odds Ratio") "aOR" 
                     else if(effect_label == "Risk Ratio") "aRR" 
                     else if(effect_label == "Rate Ratio") "aRR"
                     else if(effect_label == "Hazard Ratio") "aHR"
                     else if(effect_label == "Coefficient") "Adj. Coefficient"
                     else if(effect_label == "Exp(Coefficient)") "Exp(Adj. Coef)"
                     else effect_label  ## Use custom label as-is for the header
    
    effect_header <- paste0(effect_abbrev, " (95% CI); p-value")
    
    ## Calculate max effect string length (without expression parsing)
    effect_lengths <- nchar(data.table::fifelse(
                                            is.na(to_show_exp_clean$estimate),
                                            ref_label,
                                            paste0(to_show_exp_clean$effect_formatted, " (",
                                                   to_show_exp_clean$conf_low_formatted, "-",
                                                   to_show_exp_clean$conf_high_formatted, "); p = ",
                                                   to_show_exp_clean$p_formatted)
                                        )) + center_padding
    
    effect_width_chars <- max(effect_lengths, nchar(effect_header), na.rm = TRUE)
    columns$effect <- effect_width_chars * char_to_inch
    
    ## Calculate total table width
    calc_table_width <- sum(unlist(columns)) + length(columns) * margin * 2
    
    ## Calculate forest width based on table_width proportion
    calc_forest_width <- calc_table_width * ((1 - table_width) / table_width)
    
    ## Convert table width to log scale units
    calc_range_width <- diff(rangeb)
    calc_table_width_log_units <- (calc_table_width / calc_forest_width) * calc_range_width
    
    ## Calculate positions in log scale
    ## Start from the extended left edge
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
        total_width = calc_table_width + calc_forest_width,
        effect_abbrev = effect_abbrev
    ))
}
