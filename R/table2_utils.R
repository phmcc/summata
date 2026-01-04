## Note: is_reference_category(), is_affirmative_category(), and should_condense_binary()
## are defined in forest_utils.R as the canonical location. They are also used by 
## condense_table_rows() below. Since R packages load all files in the R/ directory, 
## these functions will be available when the package is loaded. For standalone use 
## of this file, ensure forest_utils.R is sourced first.


#' Add padding to exported table headers
#' 
#' Adds LaTeX vertical spacing rules to column headers for proper vertical
#' alignment in PDF/LaTeX exports.
#' 
#' @param col_names Character vector of column names.
#' @return Character vector with LaTeX padding rules added.
#' @keywords internal
add_header_padding <- function (col_names) {
    padded_names <- paste0("\\rule{0pt}{3ex}", col_names, "\\rule[-1.5ex]{0pt}{0pt}")
    return(as.character(padded_names))
}

#' Add padding to exported table variables
#' 
#' Inserts blank padding rows between variable groups in exported tables
#' for improved visual separation.
#' 
#' @param df Data.table with Variable column.
#' @return Data.table with padding rows inserted between variable groups.
#' @keywords internal
add_variable_padding <- function (df) {
    var_col <- if ("Variable" %chin% names(df)) 
                   "Variable"
               else if ("variable" %chin% names(df)) 
                   "variable"
               else NULL
    if (is.null(var_col)) 
        return(df)
    var_starts <- which(df[[var_col]] != "" & !grepl("^\\\\hspace", 
                                                     df[[var_col]]) & !grepl("^\\s+", df[[var_col]]))
    if (length(var_starts) == 0) 
        return(df)
    new_df <- data.table()
    padding_row <- df[1, ]
    padding_row[1, ] <- ""
    new_df <- rbind(new_df, padding_row)
    for (i in seq_along(var_starts)) {
        if (i < length(var_starts)) {
            var_rows <- var_starts[i]:(var_starts[i + 1] - 1)
        }
        else {
            var_rows <- var_starts[i]:nrow(df)
        }
        new_df <- rbind(new_df, df[var_rows, ])
        padding_row <- df[1, ]
        padding_row[1, ] <- ""
        new_df <- rbind(new_df, padding_row)
    }
    return(new_df)
}

#' Check LaTeX installation
#' 
#' Verifies that a LaTeX distribution (pdflatex or xelatex) is available
#' on the system for PDF compilation.
#' 
#' @return Logical TRUE if LaTeX is available, FALSE otherwise.
#' @keywords internal
check_latex <- function () {
    pdflatex_check <- Sys.which("pdflatex")
    if (pdflatex_check != "") 
        return(TRUE)
    xelatex_check <- Sys.which("xelatex")
    if (xelatex_check != "") 
        return(TRUE)
    return(FALSE)
}

#' Determine alignment for exported tables
#' 
#' Creates column alignment string for LaTeX tables. Variable and Group
#' columns are left-aligned; all others are centered.
#' 
#' @param df Data.frame or data.table to determine alignment for.
#' @return Character string with alignment codes (e.g., "rlcc").
#' @keywords internal
determine_alignment <- function(df) {
    align <- "r"  # Start with row numbers column (if any)
    for (col in names(df)) {
        if (col %chin% c("Variable", "Group")) {
            align <- paste0(align, "l")
        } else {
            align <- paste0(align, "c")
        }
    }
    return(align)
}

#' Apply formatting to column headers in exported tables (PDF/LaTeX)
#' 
#' Formats column headers for LaTeX output by escaping special characters,
#' italicizing 'n' and 'p', and optionally adding vertical spacing.
#' 
#' @param col_names Character vector of column names.
#' @param add_header_space Logical whether to add vertical padding.
#' @return Character vector with LaTeX-formatted column names.
#' @keywords internal
format_column_headers <- function (col_names, add_header_space = TRUE) {
    if (!is.character(col_names)) {
        col_names <- as.character(col_names)
    }
    col_names <- gsub("%", "\\\\%", col_names)
    col_names <- gsub("^n$", "\\\\textit{n}", col_names, ignore.case = TRUE)
    col_names <- gsub("p-value", "\\\\textit{p}-value", col_names, 
                      ignore.case = TRUE)
    col_names <- gsub("Uni p", "Uni \\\\textit{p}", col_names)
    col_names <- gsub("Multi p", "Multi \\\\textit{p}", col_names)
    if (add_header_space && !is.null(col_names)) {
        col_names <- add_header_padding(col_names)
    }
    return(as.character(col_names))
}

#' Apply formatting to column headers in exported tables (HTML)
#' 
#' Formats column headers for HTML output by italicizing 'n' and 'p'
#' using HTML tags.
#' 
#' @param col_names Character vector of column names.
#' @return Character vector with HTML-formatted column names.
#' @keywords internal
format_column_headers_html <- function (col_names) {
    col_names <- gsub("^n$", "<i>n</i>", col_names, ignore.case = TRUE)
    col_names <- gsub("p-value", "<i>p</i>-value", col_names, 
                      ignore.case = TRUE)
    return(col_names)
}

#' Apply formatting to indented groups
#' 
#' Transforms tables with Variable/Group columns into indented format where
#' group levels appear as indented rows under variable names. Handles both
#' regression and descriptive tables with appropriate p-value placement.
#' 
#' @param df Data.table with Variable and Group columns.
#' @param indent_string Character string to use for indentation.
#' @return Data.table with Group column removed and levels indented under Variables.
#' @keywords internal
format_indented_groups <- function (df, indent_string = "    ") {
    if (!("Variable" %chin% names(df) && "Group" %chin% names(df))) {
        return(df)
    }
    
    var_rows <- which(df$Variable != "")
    p_cols <- grep("p-value|Uni p|Multi p", names(df), value = TRUE)
    effect_cols <- c("OR (95% CI)", "HR (95% CI)", "RR (95% CI)", "Estimate (95% CI)",
                     "Coefficient (95% CI)",  # ADD THIS
                     "Univariable OR (95% CI)", "Multivariable aOR (95% CI)", 
                     "Univariable HR (95% CI)", "Multivariable aHR (95% CI)",
                     "Univariable RR (95% CI)", "Multivariable aRR (95% CI)",
                     "Univariable Estimate (95% CI)", "Multivariable Estimate (95% CI)",
                     "Univariable Coefficient (95% CI)", "Multivariable Coefficient (95% CI)")
    is_regression_table <- any(effect_cols %chin% names(df))
    is_fullfit_table <- any(c("Uni p", "Multi p") %chin% names(df))
    new_df <- data.table()
    
    for (i in seq_along(var_rows)) {
        current <- var_rows[i]
        next_var <- if (i < length(var_rows)) 
                        var_rows[i + 1]
                    else nrow(df) + 1
        var_name <- df$Variable[current]
        has_group <- df$Group[current] != "" && !is.na(df$Group[current]) && 
            df$Group[current] != "-"
        
        if (has_group) {
            if (is_regression_table) {
                var_row <- df[current, ]
                var_row$Variable <- var_name
                var_row$Group <- ""
                data_cols <- setdiff(names(df), c("Variable", "Group", p_cols))
                for (col in data_cols) {
                    var_row[[col]] <- ""
                }
                new_df <- rbind(new_df, var_row)
                
                for (j in current:(next_var - 1)) {
                    group_row <- df[j, ]
                    group_row$Variable <- paste0(indent_string, 
                                                 df$Group[j])
                    if (is_fullfit_table) {
                        for (ec in effect_cols) {
                            if (ec %chin% names(group_row)) {
                                val <- group_row[[ec]]
                                if (!is.na(val) && (val == "-" || grepl("Reference", 
                                                                        val))) {
                                    if (grepl("Univariable", ec) && "Uni p" %chin% 
                                        names(group_row)) {
                                        group_row[["Uni p"]] <- "-"
                                    }
                                    else if (grepl("Multivariable", ec) && 
                                             "Multi p" %chin% names(group_row)) {
                                        group_row[["Multi p"]] <- "-"
                                    }
                                }
                            }
                        }
                    }
                    else {
                        is_reference <- FALSE
                        for (ec in effect_cols) {
                            if (ec %chin% names(group_row)) {
                                val <- group_row[[ec]]
                                if (!is.na(val) && (val == "-" || grepl("Reference", 
                                                                        val))) {
                                    is_reference <- TRUE
                                    break
                                }
                            }
                        }
                        if (is_reference) {
                            for (p_col in p_cols) {
                                if (p_col %chin% names(group_row)) {
                                    group_row[[p_col]] <- ""
                                }
                            }
                        }
                    }
                    new_df <- rbind(new_df, group_row)
                }
            }
            else {
                var_row <- df[current, ]
                var_row$Variable <- var_name
                data_cols <- setdiff(names(df), c("Variable", 
                                                  "Group", p_cols))
                for (col in data_cols) {
                    var_row[[col]] <- ""
                }
                new_df <- rbind(new_df, var_row)
                
                for (j in current:(next_var - 1)) {
                    group_row <- df[j, ]
                    group_row$Variable <- paste0(indent_string, 
                                                 df$Group[j])
                    for (p_col in p_cols) {
                        if (p_col %chin% names(group_row)) {
                            group_row[[p_col]] <- ""
                        }
                    }
                    new_df <- rbind(new_df, group_row)
                }
            }
        }
        else {
            new_df <- rbind(new_df, df[current, ])
        }
    }
    new_df$Group <- NULL
    rownames(new_df) <- NULL
    return(new_df)
}

#' Format p-values for exported tables
#' 
#' Applies bold formatting to significant p-values in LaTeX tables using
#' the textbf command.
#' 
#' @param df Data.table containing p-value columns.
#' @param p_threshold Numeric threshold for significance (default 0.05).
#' @return Data.table with significant p-values wrapped in LaTeX bold commands.
#' @keywords internal
format_pvalues_export_tex <- function (df, p_threshold = 0.05) {
    for (col in names(df)) {
        if (col == "p-value" || col == "Uni p" || col == "Multi p" || 
            grepl("p.value|pvalue", col, ignore.case = TRUE)) {
            for (i in seq_len(nrow(df))) {
                cell_value <- as.character(df[[col]][i])
                if (is.na(cell_value) || cell_value == "" || 
                    cell_value == "NA" || grepl("\\\\textbf", cell_value)) {
                    next
                }
                is_significant <- FALSE
                if (grepl("^<\\s*0\\.001", cell_value)) {
                    is_significant <- TRUE
                }
                else if (grepl("^[0-9]\\.[0-9]", cell_value) || 
                         grepl("^0\\.[0-9]", cell_value)) {
                    p_numeric <- suppressWarnings(as.numeric(gsub("[^0-9.]", 
                                                                  "", cell_value)))
                    if (!is.na(p_numeric) && p_numeric < p_threshold) {
                        is_significant <- TRUE
                    }
                }
                if (is_significant) {
                    df[[col]][i] <- paste0("\\textbf{", cell_value, 
                                           "}")
                }
            }
        }
    }
    return(df)
}

#' Format p-values for exported tables (HTML)
#' 
#' Applies bold formatting to significant p-values in HTML tables using
#' the b tag.
#' 
#' @param df Data.table containing p-value columns.
#' @param p_threshold Numeric threshold for significance (default 0.05).
#' @return Data.table with significant p-values wrapped in HTML bold tags.
#' @keywords internal
format_pvalues_export_html <- function (df, p_threshold = 0.05) {
    for (col in names(df)) {
        if (col == "p-value" || col == "Uni p" || col == "Multi p" || 
            grepl("p.value|pvalue", col, ignore.case = TRUE)) {
            for (i in seq_len(nrow(df))) {
                cell_value <- as.character(df[[col]][i])
                if (is.na(cell_value) || cell_value == "" || 
                    cell_value == "NA") {
                    next
                }
                is_significant <- FALSE
                if (grepl("^<\\s*0\\.001", cell_value)) {
                    is_significant <- TRUE
                }
                else if (grepl("^[0-9]\\.[0-9]", cell_value) || 
                         grepl("^0\\.[0-9]", cell_value)) {
                    p_numeric <- suppressWarnings(as.numeric(gsub("[^0-9.]", 
                                                                  "", cell_value)))
                    if (!is.na(p_numeric) && p_numeric < p_threshold) {
                        is_significant <- TRUE
                    }
                }
                if (is_significant) {
                    df[[col]][i] <- paste0("<b>", cell_value, "</b>")
                }
            }
        }
    }
    return(df)
}

#' Get paper size for PDF/LaTeX export
#' 
#' Returns paper dimensions and margin settings for the specified paper size.
#' 
#' @param paper Character string: "letter", "a4", or "auto".
#' @param margins Optional numeric vector of margins (length 1 or 4).
#' @return List with latex_paper, width, height, and margins components.
#' @keywords internal
get_paper_settings <- function (paper, margins = NULL) {
    paper <- match.arg(paper, c("letter", "a4", "auto"))
    settings <- switch(paper,
                       letter = list(latex_paper = "letterpaper", 
                                     width = 8.5,
                                     height = 11,
                                     default_margins = c(1, 1, 1, 1)),
                       a4 = list(latex_paper = "a4paper",
                                 width = 8.27, 
                                 height = 11.69,
                                 default_margins = c(1, 1, 1, 1)),
                       auto = list(latex_paper = "letterpaper", 
                                   width = NULL, height = NULL, default_margins = c(0.5, 0.5, 0.5, 0.5)))
    if (!is.null(margins)) {
        if (length(margins) == 1) {
            margins <- rep(margins, 4)
        }
        else if (length(margins) != 4) {
            stop("margins must be length 1 or 4")
        }
        settings$margins <- margins
    }
    else {
        settings$margins <- settings$default_margins
    }
    return(settings)
}

#' Sanitize certain symbols for LaTeX
#' 
#' Escapes special LaTeX characters (%, &, #, _, $, ^, ~, {, }) while
#' preserving existing LaTeX commands. Uses negative lookbehind to avoid
#' double-escaping already escaped characters.
#' 
#' @param x Character vector to sanitize.
#' @return Character vector with special characters escaped for LaTeX.
#' @keywords internal
sanitize_for_latex <- function(x) {
    if (is.null(x) || length(x) == 0) 
        return(x)
    
    ## Process each element
    result <- x
    
    for (i in seq_along(x)) {
        if (is.na(x[i])) next
        
        val <- x[i]
        
        ## Check for already-formatted LaTeX commands
        has_latex <- grepl("\\\\(text(bf|it|tt|sc|sl|rm)|hspace|vspace|rule|begin|end|[a-zA-Z]+\\{)", val)
        
        if (has_latex) {
            ## String has LaTeX commands - only escape unescaped special chars
            ## Escape % that is not already escaped (not preceded by \)
            val <- gsub("(?<!\\\\)%", "\\\\%", val, perl = TRUE)
            ## Escape & that is not already escaped
            val <- gsub("(?<!\\\\)&", "\\\\&", val, perl = TRUE)
            ## Escape # that is not already escaped
            val <- gsub("(?<!\\\\)#", "\\\\#", val, perl = TRUE)
            ## Escape _ that is not already escaped
            val <- gsub("(?<!\\\\)_", "\\\\_", val, perl = TRUE)
            ## Escape $ that is not already escaped
            val <- gsub("(?<!\\\\)\\$", "\\\\$", val, perl = TRUE)
        } else {
            ## No LaTeX commands - check if already escaped
            already_escaped <- grepl("\\\\[%&#_$]", val)
            
            if (!already_escaped) {
                ## Full sanitization needed
                val <- gsub("\\\\", "\\\\textbackslash{}", val)
                val <- gsub("%", "\\\\%", val)
                val <- gsub("&", "\\\\&", val)
                val <- gsub("#", "\\\\#", val)
                val <- gsub("_", "\\\\_", val)
                val <- gsub("\\$", "\\\\$", val)
                val <- gsub("\\^", "\\\\textasciicircum{}", val)
                val <- gsub("~", "\\\\textasciitilde{}", val)
                val <- gsub("\\{", "\\\\{", val)
                val <- gsub("\\}", "\\\\}", val)
            }
        }
        
        result[i] <- val
    }
    
    return(result)
}

#' Format column headers with n counts (TeX)
#' 
#' Creates LaTeX-formatted column headers with sample size counts displayed
#' below the column name in a stacked format.
#' 
#' @param col_names Character vector of column names.
#' @param n_row_data Named list or data.table row with n values for each column.
#' @return Character vector with LaTeX-formatted headers including n counts.
#' @keywords internal
format_column_headers_with_n_tex <- function(col_names, n_row_data) {
    new_names <- character(length(col_names))
    
    for (i in seq_along(col_names)) {
        col <- col_names[i]
        
        ## Sanitize the column name first
        col_sanitized <- sanitize_for_latex(col)
        
        has_n <- col %chin% names(n_row_data) && 
            !is.na(n_row_data[[col]]) && 
            n_row_data[[col]] != "" && 
            n_row_data[[col]] != "0"
        
        if (col %chin% c("Variable", "Group")) {
            new_names[i] <- format_column_headers(col_sanitized)
        } else if (col %chin% c("p-value", "p value")) {
            new_names[i] <- format_column_headers(col_sanitized)
        } else if (has_n) {
            n_value <- n_row_data[[col]]
            ## Build the complex structure with already-sanitized name
            new_names[i] <- paste0("\\begin{tabular}{@{}c@{}}\\rule{0pt}{2.5ex}", 
                                   col_sanitized,  # Use sanitized version
                                   "\\\\[-0ex] (\\textit{n} = ", 
                                   format(as.numeric(gsub(",", "", n_value)), big.mark = ","), 
                                   ")\\rule[-1ex]{0pt}{0pt}\\end{tabular}")
        } else {
            new_names[i] <- format_column_headers(col_sanitized)
        }
    }
    
    return(new_names)
}

#' Format column headers with n counts (HTML)
#' 
#' Creates HTML-formatted column headers with sample size counts displayed
#' below the column name using line breaks.
#' 
#' @param col_names Character vector of column names.
#' @param n_row_data Named list or data.table row with n values for each column.
#' @return Character vector with HTML-formatted headers including n counts.
#' @keywords internal
format_column_headers_with_n_html <- function(col_names, n_row_data) {
    new_names <- col_names
    for (i in seq_along(col_names)) {
        col <- col_names[i]
        
        ## Skip Variable, Group, and p-value columns
        if (col %chin% c("Variable", "Group", "p-value", "p value")) {
            new_names[i] <- format_column_headers_html(col)
            next
        }
        
        ## Get n value from the N row
        if (col %chin% names(n_row_data)) {
            n_value <- n_row_data[[col]]
            if (!is.na(n_value) && n_value != "" && n_value != "0") {
                ## Format header with only n italicized
                clean_col <- format_column_headers_html(col)
                new_names[i] <- paste0(clean_col, "<br>(<i>n</i> = ", n_value, ")")
            } else {
                new_names[i] <- format_column_headers_html(col)
            }
        } else {
            new_names[i] <- format_column_headers_html(col)
        }
    }
    return(new_names)
}

#' Condense table rows for more compact display
#' 
#' Collapses multi-row variables into single rows for compact tables.
#' Continuous variables show only the first statistic row, binary categorical
#' variables show only the non-reference category, and survival variables
#' show only the median row.
#' 
#' @param df Data.table with Variable and Group columns.
#' @param indent_groups Logical whether indentation will be applied (affects processing).
#' @return Data.table with condensed rows.
#' @keywords internal
condense_table_rows <- function(df, indent_groups = TRUE) {
    
    ## Work with a copy
    result <- data.table::copy(as.data.table(df))
    
    ## Detect table type
    is_descriptive <- !any(grepl("(OR|HR|RR|Coefficient|Estimate).*\\(95% CI\\)", names(result)))
    
    ## Track rows to delete
    rows_to_delete <- integer()
    
    ## Process each variable
    vars_to_process <- which(result$Variable != "" & !is.na(result$Variable))
    
    for (i in seq_along(vars_to_process)) {
        var_start <- vars_to_process[i]
        var_end <- if (i < length(vars_to_process)) {
                       vars_to_process[i + 1] - 1
                   } else {
                       nrow(result)
                   }
        
        if (var_start %in% rows_to_delete) next
        
        var_name <- result$Variable[var_start]
        var_rows <- result[var_start:var_end]
        n_rows <- nrow(var_rows)
        
        if ("Group" %chin% names(var_rows)) {
            groups <- var_rows$Group
            
            ## Check for continuous variable
            is_continuous <- any(grepl("^(Mean|Median|Range|SD|IQR)", groups[1], ignore.case = TRUE))
            
            ## Check for binary categorical
            non_empty_groups <- groups[groups != "" & !is.na(groups)]
            is_binary <- length(non_empty_groups) == 2 && !is_continuous
            
            ## Check for survival
            is_survival <- any(grepl("Median.*\\(.*CI.*\\)", groups, ignore.case = TRUE))
            
            if (is_continuous || is_survival) {
                ## For continuous/survival, keep only first row
                if (is_descriptive) {
                    stat_type <- gsub("\\s*\\(.*\\)", "", groups[1])
                    stat_type <- gsub("\\s*\\+/-.*", "", stat_type)
                    if (stat_type != "" && !is.na(stat_type) && stat_type != "-") {
                        if (grepl("Mean \\+/- SD", groups[1])) {
                            result[var_start, Variable := paste0(var_name, ", mean +/- SD")]
                        } else if (grepl("Median", groups[1]) && grepl("IQR", groups[1])) {
                            result[var_start, Variable := paste0(var_name, ", median [IQR]")]
                        } else if (grepl("Median", groups[1]) && is_survival) {
                            result[var_start, Variable := paste0(var_name, ", median (95% CI)")]
                        } else {
                            result[var_start, Variable := paste0(var_name, ", ", tolower(stat_type))]
                        }
                    }
                    ## Clear Group column
                    result[var_start, Group := ""]
                }
                
                ## Mark extra rows for deletion
                if (n_rows > 1) {
                    rows_to_delete <- c(rows_to_delete, (var_start + 1):var_end)
                }
                
            } else if (is_binary) {
                
                ## For binary categorical
                data_cols <- setdiff(names(result), c("Variable", "Group", "p-value", "p.value"))
                
                if (length(data_cols) > 0) {
                    first_data_col <- data_cols[1]
                    ## Find non-reference row by checking for actual data (not "-", "reference", etc.)
                    non_ref_idx <- which(!var_rows[[first_data_col]] %chin% c("-", "reference", "Reference", ""))
                    
                    if (length(non_ref_idx) > 1) {
                        ## Multiple non-ref rows - take the first one
                        non_ref_idx <- non_ref_idx[1]
                    } else if (length(non_ref_idx) == 0) {
                        ## All rows are reference-like - skip condensing
                        next
                    }
                    
                    if (non_ref_idx <= n_rows) {
                        non_ref_row <- var_start + non_ref_idx - 1
                        ref_idx <- setdiff(1:n_rows, non_ref_idx)[1]
                        ref_row <- var_start + ref_idx - 1
                        
                        non_ref_category <- result$Group[non_ref_row]
                        ref_category <- result$Group[ref_row]
                        
                        if (!is.na(non_ref_category) && non_ref_category != "") {
                            ## Use greedy helper function for condensing detection
                            ## Pass var_name as label for case-insensitive matching
                            if (should_condense_binary(ref_category, non_ref_category, var_name)) {
                                result[var_start, Variable := paste0(var_name)]
                            } else {
                                result[var_start, Variable := paste0(var_name, " (", non_ref_category, ")")]                                
                            }
                        }
                        
                        ## Copy statistics
                        for (col in data_cols) {
                            if (col %chin% names(result)) {
                                result[var_start, (col) := result[non_ref_row, get(col)]]
                            }
                        }
                        
                        ## Copy p-value
                        if ("p-value" %chin% names(result)) {
                            pval <- result[var_start:var_end, `p-value`]
                            pval <- pval[!is.na(pval) & pval != ""]
                            if (length(pval) > 0) {
                                result[var_start, `p-value` := pval[1]]
                            }
                        }
                        
                        ## Clear Group column
                        result[var_start, Group := ""]
                        
                        ## Mark other rows for deletion
                        rows_to_delete <- c(rows_to_delete, (var_start + 1):var_end)
                    }
                }
            } 
        }
    }
    
    ## Remove marked rows
    if (length(rows_to_delete) > 0) {
        rows_to_delete <- sort(unique(rows_to_delete))
        rows_to_delete <- rows_to_delete[rows_to_delete <= nrow(result)]
        if (length(rows_to_delete) > 0) {
            result <- result[-rows_to_delete]
        }
    }
    
    return(result)
}


#' Condense quantitative variable rows only
#' 
#' Collapses multi-row continuous and survival variables into single rows
#' while preserving all categorical variable rows (including binary).
#' Only applies to descriptive tables from desctable().
#' 
#' @param df A data.table or data.frame
#' @param indent_groups Logical. Whether to apply indentation formatting.
#' @return A data.table with condensed continuous/survival rows
#' @keywords internal
condense_quantitative_rows <- function(df, indent_groups = TRUE) {
    
    ## Work with a copy
    result <- data.table::copy(as.data.table(df))
    
    ## Detect table type - only apply to descriptive tables
    is_descriptive <- !any(grepl("(OR|HR|RR|Coefficient|Estimate).*\\(95% CI\\)", names(result)))
    
    if (!is_descriptive) {
        warning("condense_quantitative only applies to descriptive tables. Returning unchanged.")
        return(result)
    }
    
    ## Track rows to delete
    rows_to_delete <- integer()
    
    ## Process each variable
    vars_to_process <- which(result$Variable != "" & !is.na(result$Variable))
    
    for (i in seq_along(vars_to_process)) {
        var_start <- vars_to_process[i]
        var_end <- if (i < length(vars_to_process)) {
                       vars_to_process[i + 1] - 1
                   } else {
                       nrow(result)
                   }
        
        if (var_start %in% rows_to_delete) next
        
        var_name <- result$Variable[var_start]
        var_rows <- result[var_start:var_end]
        n_rows <- nrow(var_rows)
        
        if ("Group" %chin% names(var_rows)) {
            groups <- var_rows$Group
            
            ## Check for continuous variable
            is_continuous <- any(grepl("^(Mean|Median|Range|SD|IQR)", groups[1], ignore.case = TRUE))
            
            ## Check for survival
            is_survival <- any(grepl("Median.*\\(.*CI.*\\)", groups, ignore.case = TRUE))
            
            ## Only condense continuous and survival - NOT categorical (including binary)
            if (is_continuous || is_survival) {
                ## For continuous/survival, keep only first row
                stat_type <- gsub("\\s*\\(.*\\)", "", groups[1])
                stat_type <- gsub("\\s*\\+/-.*", "", stat_type)
                if (stat_type != "" && !is.na(stat_type) && stat_type != "-") {
                    if (grepl("Mean \\+/- SD", groups[1])) {
                        result[var_start, Variable := paste0(var_name, ", mean +/- SD")]
                    } else if (grepl("Median", groups[1]) && grepl("IQR", groups[1])) {
                        result[var_start, Variable := paste0(var_name, ", median [IQR]")]
                    } else if (grepl("Median", groups[1]) && is_survival) {
                        result[var_start, Variable := paste0(var_name, ", median (95% CI)")]
                    } else {
                        result[var_start, Variable := paste0(var_name, ", ", tolower(stat_type))]
                    }
                }
                ## Clear Group column
                result[var_start, Group := ""]
                
                ## Mark extra rows for deletion
                if (n_rows > 1) {
                    rows_to_delete <- c(rows_to_delete, (var_start + 1):var_end)
                }
            }
            ## Note: Categorical variables (including binary) are not processed here
        }
    }
    
    ## Remove marked rows
    if (length(rows_to_delete) > 0) {
        rows_to_delete <- sort(unique(rows_to_delete))
        rows_to_delete <- rows_to_delete[rows_to_delete <= nrow(result)]
        if (length(rows_to_delete) > 0) {
            result <- result[-rows_to_delete]
        }
    }
    
    return(result)
}


#' Core flextable processing function
#' 
#' Central processing function for creating flextable objects from data tables.
#' Handles N row extraction, condensing, indentation, zebra stripes, formatting,
#' and styling. Used by table2docx, table2pptx, and table2rtf.
#' 
#' @param table Data.frame or data.table to process.
#' @param caption Optional character string for table caption.
#' @param font_size Numeric font size in points.
#' @param font_family Character string font family name.
#' @param format_headers Logical whether to format headers.
#' @param bold_significant Logical whether to bold significant p-values.
#' @param p_threshold Numeric p-value threshold for significance.
#' @param indent_groups Logical whether to indent group levels.
#' @param condense_table Logical whether to condense all variable types.
#' @param condense_quantitative Logical whether to condense only continuous/survival.
#' @param zebra_stripes Logical whether to apply alternating row shading.
#' @param dark_header Logical whether to use dark header style.
#' @param bold_variables Logical whether to bold variable names (non-indented rows).
#' @param paper Character string paper size.
#' @param orientation Character string page orientation.
#' @param width Optional numeric table width in inches.
#' @param align Optional alignment specification.
#' @return List with ft (flextable object) and caption components.
#' @keywords internal
process_table_for_flextable <- function(table,
                                        caption = NULL,
                                        font_size = 10,
                                        font_family = "Arial",
                                        format_headers = TRUE,
                                        bold_significant = TRUE,
                                        p_threshold = 0.05,
                                        indent_groups = FALSE,
                                        condense_table = FALSE,
                                        condense_quantitative = FALSE,
                                        zebra_stripes = FALSE,
                                        dark_header = FALSE,
                                        bold_variables = TRUE,
                                        paper = "letter",
                                        orientation = "portrait",
                                        width = NULL,
                                        align = NULL) {
    
    ## Convert to data.table
    df <- data.table::as.data.table(table)
    
    ## Handle N row if present
    has_n_row <- FALSE
    n_row_data <- NULL
    if (nrow(df) > 0 && "Variable" %chin% names(df) && 
        !is.na(df$Variable[1]) && df$Variable[1] == "N") {
        has_n_row <- TRUE
        n_row_data <- df[1, ]
        df <- df[-1, ]
    }
    
    ## Track variable groups BEFORE any transformation
    var_groups <- NULL
    if (zebra_stripes && "Variable" %chin% names(df)) {
        var_groups <- identify_variable_groups(df)
    }
    
    ## Apply condensing if requested
    if (condense_table) {
        indent_groups <- TRUE
        df <- condense_table_rows(df, indent_groups = indent_groups)
        
        ## Update variable groups after condensing
        if (zebra_stripes) {
            var_groups <- identify_variable_groups(df)
        }
        
        df <- format_indented_groups(df, indent_string = "    ")
    } else if (condense_quantitative) {
        ## Only condense continuous/survival variables (not categorical)
        ## Also set indent_groups = TRUE to avoid awkward empty Group column
        indent_groups <- TRUE
        df <- condense_quantitative_rows(df, indent_groups = indent_groups)
        
        ## Update variable groups after condensing
        if (zebra_stripes) {
            var_groups <- identify_variable_groups(df)
        }
        
        df <- format_indented_groups(df, indent_string = "    ")
    } else if (indent_groups) {
        df <- format_indented_groups(df, indent_string = "    ")
    }
    
    ## Replace empty cells with "-" for consistency
    df <- replace_empty_cells(df)
    
    ## Create flextable
    ft <- flextable::flextable(df)
    
    ## Set font
    ft <- flextable::font(ft, fontname = font_family, part = "all")
    ft <- flextable::fontsize(ft, size = font_size, part = "all")
    
    ## Format headers
    if (format_headers) {
        ft <- format_headers_ft(ft, has_n_row, n_row_data)
    }

    ## Apply dark header if requested
    if (dark_header) {
        ft <- flextable::bg(ft, bg = "#000000", part = "header")
        ft <- flextable::color(ft, color = "#FFFFFF", part = "header")
        ft <- flextable::bold(ft, bold = TRUE, part = "header")
    }
    
    ## Bold significant p-values
    if (bold_significant) {
        ft <- bold_pvalues_ft(ft, df, p_threshold)
    }
    
    ## Bold variable names (non-indented rows in Variable column)
    if (bold_variables && "Variable" %chin% names(df)) {
        ## Find rows where Variable is not empty and not indented (doesn't start with spaces)
        var_rows <- which(!is.na(df$Variable) & 
                          df$Variable != "" & 
                          df$Variable != "-" &
                          !grepl("^\\s+", df$Variable))
        if (length(var_rows) > 0) {
            ft <- flextable::bold(ft, i = var_rows, j = "Variable", part = "body")
        }
    }
    
    ## Set alignment
    if (is.null(align)) {
        for (col in names(df)) {
            if (col %chin% c("Variable", "Group")) {
                ft <- flextable::align(ft, j = col, align = "left", part = "all")
            } else {
                ft <- flextable::align(ft, j = col, align = "center", part = "all")
            }
        }
    } else {
        if (length(align) == 1) {
            ft <- flextable::align(ft, align = align, part = "all")
        } else if (length(align) == ncol(df)) {
            for (i in seq_along(align)) {
                ft <- flextable::align(ft, j = i, align = align[i], part = "all")
            }
        }
    }
    
    ## Add borders
    ft <- flextable::border_remove(ft)
    ft <- flextable::hline_top(ft, border = officer::fp_border(width = 2), part = "header")
    ft <- flextable::hline_bottom(ft, border = officer::fp_border(width = 1), part = "header")
    ft <- flextable::hline_bottom(ft, border = officer::fp_border(width = 2), part = "body")
    
    ## Reduce line spacing and padding
    ft <- flextable::line_spacing(ft, space = 1)
    ft <- flextable::padding(ft, i = NULL, j = NULL,
                             padding.top = 1, padding.bottom = 1,
                             padding.left = 1, padding.right = 1)
    
    ## Add zebra stripes by variable group if requested
    if (zebra_stripes && !is.null(var_groups)) {
        ft <- apply_zebra_stripes_ft(ft, df, var_groups)
    }
    
    ## Calculate width based on paper and orientation if not specified
    if (is.null(width)) {
        width <- calculate_table_width(paper, orientation)
    }
    
    ## Set width
    ft <- flextable::width(ft, width = width / ncol(df))
    
    return(list(ft = ft, caption = caption))
}

#' Identify variable groups before indentation
#' 
#' Detects variable group boundaries by finding rows where Variable column
#' is non-empty. Returns row indices for each group for zebra stripe application.
#' 
#' @param df Data.table with Variable column.
#' @return List of integer vectors, each containing row indices for one variable group.
#' @keywords internal
identify_variable_groups <- function(df) {
    if (!"Variable" %in% names(df)) return(NULL)
    
    var_starts <- which(df$Variable != "" & !is.na(df$Variable))
    if (length(var_starts) == 0) return(NULL)
    
    ## Vectorized group creation
    var_ends <- c(var_starts[-1] - 1, nrow(df))
    groups <- mapply(seq, var_starts, var_ends, SIMPLIFY = FALSE)
    
    return(groups)
}

#' Replace empty cells with "-"
#' 
#' Converts empty strings and NA values to "-" for consistent display
#' in exported tables. Preserves Variable column values.
#' 
#' @param df Data.frame or data.table to process.
#' @return Data.table with empty cells replaced by "-".
#' @keywords internal
replace_empty_cells <- function(df) {
    ## Convert to data.table for efficient in-place modification
    dt <- data.table::as.data.table(df)
    
    ## Get columns to process (excluding Variable)
    cols_to_process <- setdiff(names(dt), "Variable")
    
    ## Vectorized replacement using data.table
    for (col in cols_to_process) {
        dt[is.na(get(col)) | get(col) == "", (col) := "-"]
    }
    
    return(data.table::as.data.table(dt))
}

#' Apply zebra stripes with proper variable group detection for indented tables
#' 
#' Applies alternating background colors to variable groups in flextable objects.
#' Handles both indented tables (detects groups by leading whitespace) and
#' non-indented tables (uses pre-identified groups).
#' 
#' @param ft A flextable object.
#' @param df The source data.table used to create the flextable.
#' @param var_groups List of row index vectors for variable groups.
#' @return Flextable object with zebra stripe formatting applied.
#' @keywords internal
apply_zebra_stripes_ft <- function(ft, df, var_groups) {
    ## Check if table has been indented (look for leading spaces in Variable column)
    is_indented <- any(grepl("^\\s{2,}", df$Variable))
    
    if (is_indented) {
        ## For indented tables, identify variable groups by non-indented rows
        var_starts <- which(!grepl("^\\s", df$Variable) & df$Variable != "")
        
        for (i in seq_along(var_starts)) {
            start_row <- var_starts[i]
            end_row <- if (i < length(var_starts)) {
                           var_starts[i + 1] - 1
                       } else {
                           nrow(df)
                       }
            
            if (i %% 2 == 1) {  ## Odd variable groups get gray shading
                ft <- flextable::bg(ft, i = start_row:end_row, 
                                    bg = "#EEEEEE", part = "body")
            } else {  ## Even variable groups get white background
                ft <- flextable::bg(ft, i = start_row:end_row, 
                                    bg = "#FFFFFF", part = "body")
            }
        }
    } else if (!is.null(var_groups)) {
        ## Use pre-identified groups for non-indented tables
        for (i in seq_along(var_groups)) {
            rows <- var_groups[[i]]
            rows <- rows[rows <= nrow(df)]
            if (length(rows) > 0) {
                if (i %% 2 == 1) {  ## Odd variable groups get gray shading
                    ft <- flextable::bg(ft, i = rows, bg = "#EEEEEE", part = "body")
                } else {  ## Even variable groups get white background
                    ft <- flextable::bg(ft, i = rows, bg = "#FFFFFF", part = "body")
                }
            }
        }
    } else {
        ## Fallback to row-based striping
        odd_rows <- seq(1, nrow(df), 2)
        even_rows <- seq(2, nrow(df), 2)
        ft <- flextable::bg(ft, i = odd_rows, bg = "#EEEEEE", part = "body")
        ft <- flextable::bg(ft, i = even_rows, bg = "#FFFFFF", part = "body")
    }
    
    return(ft)
}

#' Calculate table width based on paper size and orientation
#' 
#' Computes usable table width in inches based on paper dimensions and
#' orientation, accounting for standard 1-inch margins.
#' 
#' @param paper Character string paper size ("letter", "a4", "legal").
#' @param orientation Character string page orientation ("portrait", "landscape").
#' @return Numeric usable width in inches.
#' @keywords internal
calculate_table_width <- function(paper, orientation) {
    ## Define paper sizes (in inches) with margins
    paper_sizes <- list(
        letter = c(width = 8.5, height = 11),
        a4 = c(width = 8.27, height = 11.69),
        legal = c(width = 8.5, height = 14)
    )
    
    if (!paper %chin% names(paper_sizes)) {
        paper <- "letter"
    }
    
    dims <- paper_sizes[[paper]]
    
    ## Swap for landscape
    if (orientation == "landscape") {
        dims <- c(width = dims["height"], height = dims["width"])
    }
    
    ## Subtract margins (1 inch on each side)
    usable_width <- dims["width"] - 2
    
    return(as.numeric(usable_width))
}

#' Format headers for flextable
#' 
#' Applies formatting to flextable headers including italicizing 'n',
#' adding sample size counts from N row data, and bolding all headers.
#' 
#' @param ft A flextable object.
#' @param has_n_row Logical whether source data had an N row.
#' @param n_row_data Data from the N row for adding counts to headers.
#' @return Formatted flextable object.
#' @keywords internal
format_headers_ft <- function(ft, has_n_row, n_row_data) {
    col_names <- names(ft$body$dataset)
    
    for (i in seq_along(col_names)) {
        col <- col_names[i]
        
        ## Skip Variable column for n count addition
        if (col == "Variable") {
            ## Just keep the original label without adding (n = N)
            next
        }
        
        ## Italicize 'n' in headers
        if (col == "n") {
            ft <- flextable::italic(ft, j = i, part = "header")
        }
        
        ## Add n counts for data columns only (not Variable)
        if (has_n_row && col %chin% names(n_row_data) && col != "Variable") {
            n_val <- n_row_data[[col]]
            if (!is.na(n_val) && n_val != "" && n_val != "0") {
                new_label <- paste0(col, "\n(n = ", n_val, ")")
                ft <- flextable::set_header_labels(ft, 
                                                   values = setNames(list(new_label), col))
            }
        }
    }
    
    ## Bold all headers
    ft <- flextable::bold(ft, part = "header")
    
    return(ft)
}

#' Bold significant p-values in DOCX
#' 
#' Applies bold formatting to significant p-values in flextable objects
#' by detecting values below threshold or "< 0.001" patterns.
#' 
#' @param ft A flextable object.
#' @param df The source data.table.
#' @param p_threshold Numeric p-value threshold for significance.
#' @return Flextable object with significant p-values bolded.
#' @keywords internal
bold_pvalues_ft <- function(ft, df, p_threshold = 0.05) {
    p_cols <- grep("p-value|p value|Uni p|Multi p|pvalue", names(df), 
                   ignore.case = TRUE, value = TRUE)
    
    if (length(p_cols) == 0) return(ft)
    
    for (p_col in p_cols) {
        if (p_col %in% names(df)) {
            vals <- df[[p_col]]
            
            ## Significance check
            is_very_small <- grepl("^<\\s*0\\.001", vals)
            p_numeric <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", vals)))
            is_small_numeric <- !is.na(p_numeric) & p_numeric < p_threshold
            
            sig_rows <- which((is_very_small | is_small_numeric) & 
                              vals != "" & !is.na(vals))
            
            ## Bulk bold operation
            if (length(sig_rows) > 0) {
                ft <- flextable::bold(ft, i = sig_rows, j = p_col, part = "body")
            }
        }
    }
    
    return(ft)
}
