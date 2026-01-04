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

#' Detect available sans-serif font for plots
#' 
#' Checks for commonly available sans-serif fonts in order of preference
#' (Helvetica, Arial, Helvetica Neue) and returns the first available one.
#' Falls back to "sans" if none are found or if systemfonts is unavailable.
#' 
#' @return Character string with the font family name to use.
#' @keywords internal
detect_plot_font <- function() {
    if (requireNamespace("systemfonts", quietly = TRUE)) {
        available_fonts <- systemfonts::system_fonts()$family
        ## Check fonts in order of preference
        if ("Helvetica" %in% available_fonts) return("Helvetica")
        if ("Nimbus Sans L" %in% available_fonts) return("Nimbus Sans L")
        if ("Nimbus Sans" %in% available_fonts) return("Nimbus Sans")
        if ("Liberation Sans" %in% available_fonts) return("Liberation Sans")
        if ("Arial" %in% available_fonts) return("Arial")
        if ("Helvetica Neue" %in% available_fonts) return("Helvetica Neue")
    }
    ## Fallback to generic sans-serif
    return("sans")
}

#' Format numeric value with fixed decimal places
#' 
#' Formats a numeric value to a specified number of decimal places, removing
#' leading/trailing whitespace and fixing negative zero display (e.g., "-0.00"
#' becomes "0.00").
#' 
#' @param x Numeric value to format.
#' @param digits Integer number of decimal places.
#' @return Character string with formatted value.
#' @keywords internal
format_number <- function(x, digits) {
    result <- trimws(format(round(x, digits), nsmall = digits))
    ## Fix negative zero: "-0.00" -> "0.00" (handles embedded occurrences too)
    gsub("(?<![0-9])-0(\\.0+)(?![0-9])", "0\\1", result, perl = TRUE)
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


#' Find non-reference row for binary variable condensing
#' 
#' Identifies the non-reference row in a binary categorical variable by checking
#' for NA estimates (reference rows have NA). This is more robust than assuming
#' row position or matching specific strings like "Yes" or "Positive".
#' 
#' @param var_rows Data.table containing rows for a single variable.
#' @param estimate_col Character string naming the estimate column 
#'   (e.g., "estimate", "coef"). Default is "estimate".
#' @return Integer index of the non-reference row within var_rows, or NULL if
#'   cannot be determined (e.g., no NA estimates found, or multiple non-NA rows).
#' @keywords internal
find_non_reference_row <- function(var_rows, estimate_col = "estimate") {
    if (!estimate_col %in% names(var_rows)) {
        ## Fallback: try common column names
        possible_cols <- c("estimate", "coef", "coefficient", "hr", "or", "rr")
        estimate_col <- intersect(tolower(names(var_rows)), possible_cols)[1]
        if (is.na(estimate_col)) return(NULL)
        ## Match actual case
        estimate_col <- names(var_rows)[tolower(names(var_rows)) == estimate_col][1]
    }
    
    estimates <- var_rows[[estimate_col]]
    
    ## Reference row has NA estimate; non-reference has actual value
    non_ref_idx <- which(!is.na(estimates))
    
    if (length(non_ref_idx) == 1) {
        return(non_ref_idx)
    } else if (length(non_ref_idx) == 0) {
        ## All NA - shouldn't happen for binary, but fallback to row 2
        return(2L)
    } else {
        ## Multiple non-NA values - not a simple binary reference pattern
        ## Return NULL to signal condensing should not occur
        return(NULL)
    }
}


#' Check if category name is a standard reference/negative value
#' 
#' Determines whether a category name represents a standard reference or negative
#' value that indicates absence. Used to suppress redundant category names when
#' condensing binary variables.
#' 
#' @param category Character string with the category name.
#' @param label Optional character string with the variable label. If provided,
#'   checks if category is "No [label]" or similar patterns.
#' @return Logical indicating whether category is a reference/negative value.
#' @keywords internal
is_reference_category <- function(category, label = NULL) {
    if (is.null(category) || is.na(category) || category == "") {
        return(FALSE)
    }
    
    norm_category <- tolower(trimws(category))
    
    ## Standard reference/negative values
    reference_values <- c(
        "no", "none",
        "0",
        "false",
        "absent",
        "negative",
        "normal",
        "-"
    )
    
    if (norm_category %in% reference_values) {
        return(TRUE)
    }
    
    ## Check for "No [something]" pattern
    if (grepl("^no\\s+", norm_category)) {
        return(TRUE)
    }
    
    ## Check for "Non-[something]" or "Non [something]" pattern
    if (grepl("^non[-\\s]", norm_category)) {
        return(TRUE)
    }
    
    ## Check for "Without [something]" pattern
    if (grepl("^without\\s+", norm_category)) {
        return(TRUE)
    }
    
    ## Check for "[something] absent" or "[something] negative" suffix
    if (grepl("\\s+(absent|negative)$", norm_category)) {
        return(TRUE)
    }
    
    ## If label provided, check if category is "No [label]" pattern
    if (!is.null(label) && !is.na(label) && label != "") {
        norm_label <- tolower(trimws(label))
        
        ## "No [label]" pattern (e.g., label="DVT", category="No DVT")
        if (norm_category == paste("no", norm_label)) {
            return(TRUE)
        }
        
        ## "No [partial label]" pattern
        ## e.g., label="30-Day Readmission", category="No 30-day readmission"
        if (grepl("^no\\s+", norm_category)) {
            category_suffix <- sub("^no\\s+", "", norm_category)
            if (grepl(category_suffix, norm_label, fixed = TRUE) ||
                grepl(norm_label, category_suffix, fixed = TRUE)) {
                return(TRUE)
            }
        }
    }
    
    return(FALSE)
}


#' Check if category name should be suppressed in condensed label
#' 
#' Determines whether a category name should be suppressed when condensing
#' binary variables. Returns TRUE for standard affirmative values (e.g., "Yes", 
#' "1", "Positive"), standard reference values (e.g., "No", "Absent", "None"),
#' or when the category name essentially matches the variable label 
#' (case-insensitive comparison).
#' 
#' @param category Character string with the category name.
#' @param label Optional character string with the variable label. If provided,
#'   returns TRUE when category is a case-insensitive match or substring.
#' @return Logical indicating whether category should be suppressed.
#' @keywords internal
is_affirmative_category <- function(category, label = NULL) {
    if (is.null(category) || is.na(category) || category == "") {
        return(FALSE)
    }
    
    norm_category <- tolower(trimws(category))
    
    ## Standard affirmative values
    affirmative_values <- c(
        "yes",
        "1",
        "true",
        "present",
        "positive",
        "+"
    )
    
    if (norm_category %in% affirmative_values) {
        return(TRUE)
    }
    
    ## Check if it's a reference category (also suppress these)
    if (is_reference_category(category, label)) {
        return(TRUE)
    }
    
    ## Check if category matches or is contained in the label (case-insensitive)
    if (!is.null(label) && !is.na(label) && label != "") {
        norm_label <- tolower(trimws(label))
        
        ## Check for exact match
        if (norm_category == norm_label) {
            return(TRUE)
        }
        
        ## Check if category is a substantial substring of label
        ## (e.g., "30-day readmission" in "30-Day Readmission")
        if (nchar(norm_category) >= 3 && grepl(norm_category, norm_label, fixed = TRUE)) {
            return(TRUE)
        }
        
        ## Check reverse: label is substring of category
        if (nchar(norm_label) >= 3 && grepl(norm_label, norm_category, fixed = TRUE)) {
            return(TRUE)
        }
    }
    
    return(FALSE)
}


#' Check if a binary variable should be condensed without category suffix
#' 
#' Uses a greedy/liberal approach to determine if a binary variable's condensed
#' display should omit the category name. Returns TRUE if EITHER level of the
#' binary variable is a standard reference/affirmative value, OR if either level
#' matches/contains the variable label.
#' 
#' This function is designed for binary (2-level) categorical variables where
#' one level is a reference and one is the "event" or "condition" level.
#' 
#' @param ref_category Character string with the reference category name
#'   (the level with NA estimate).
#' @param non_ref_category Character string with the non-reference category name
#'   (the level with the actual estimate).
#' @param label Optional character string with the variable label. Used for
#'   intelligent matching (e.g., "30-Day Readmission" label with 
#'   "30-day readmission" / "No 30-day readmission" levels).
#' @return Logical indicating whether the binary variable should be condensed
#'   without appending the category name.
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' # Standard Yes/No - should condense
#' should_condense_binary("No", "Yes", "Diabetes")
#' # Returns TRUE (both are standard values)
#' 
#' # Descriptive levels matching label - should condense
#' should_condense_binary("No 30-day readmission", "30-day readmission", 
#'                        "30-Day Readmission")
#' # Returns TRUE (levels match/contain label)
#' 
#' # Custom levels not matching - should NOT condense
#' should_condense_binary("Group A", "Group B", "Treatment Arm")
#' # Returns FALSE (neither is standard nor matches label)
#' }
should_condense_binary <- function(ref_category, non_ref_category, label = NULL) {
    ## Greedy approach: if EITHER category is recognizable, condense
    
    ## Check if reference category is a standard reference value
    ref_is_standard <- is_reference_category(ref_category, label)
    
    ## Check if non-reference category is a standard affirmative value
    ## or matches the label
    non_ref_is_standard <- is_affirmative_category(non_ref_category, label)
    
    ## Return TRUE if either is recognized
    return(ref_is_standard || non_ref_is_standard)
}
