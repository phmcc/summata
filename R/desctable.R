#' Create Publication-Ready Descriptive Statistics Tables
#'
#' Generates comprehensive descriptive statistics tables with automatic variable
#' type detection, group comparisons, and appropriate statistical testing. This
#' function is designed to create "Table 1"-style summaries commonly used in
#' clinical and epidemiological research, with full support for continuous,
#' categorical, and time-to-event variables.
#'
#' @param data Data frame or data.table containing the dataset to summarize.
#'   Automatically converted to a data.table for efficient processing.
#'
#' @param by Character string specifying the column name of the grouping
#'   variable for stratified analysis (\emph{e.g.}, treatment arm, exposure
#'   status). When \code{NULL} (default), produces overall summaries only
#'   without group comparisons or statistical tests.
#'
#' @param variables Character vector of variable names to summarize. Can
#'   include standard column names for continuous or categorical variables,
#'   and survival expressions using \code{Surv()} syntax (\emph{e.g.},
#'   \code{"Surv(os_months, os_status)"}). Variables are processed in the
#'   order provided.
#'
#' @param stats_continuous Character vector specifying which statistics to
#'   compute for continuous variables. Multiple values create separate rows
#'   for each variable. Options:
#'   \itemize{
#'     \item \code{"mean_sd"} - Mean \eqn{\pm} standard deviation
#'     \item \code{"median_iqr"} - Median [interquartile range]
#'     \item \code{"median_range"} - Median (minimum-maximum)
#'     \item \code{"range"} - Minimum-maximum only
#'   }
#'   Default is \code{"median_iqr"}.
#'
#' @param stats_categorical Character string specifying the format for
#'   categorical variable summaries:
#'   \itemize{
#'     \item \code{"n"} - Count only
#'     \item \code{"percent"} - Percentage only
#'     \item \code{"n_percent"} - Count (percentage) [default]
#'   }
#'
#' @param digits Integer specifying the number of decimal places for
#'   continuous statistics. Default is 1.
#'
#' @param p_digits Integer specifying the number of decimal places for
#'   \emph{p}-values. Values smaller than \code{10^(-p_digits)} are displayed
#'   as \code{"< 0.001"} (for \code{p_digits = 3}), \code{"< 0.0001"} (for
#'   \code{p_digits = 4}), etc. Default is 3.
#'
#' @param conf_level Numeric confidence level for confidence intervals in
#'   survival variable summaries (median survival time with CI). Must be
#'   between 0 and 1. Default is 0.95 (95\% confidence intervals).
#'
#' @param p_per_stat Logical. If \code{TRUE}, displays \emph{p}-values on each
#'   row (per statistic) rather than only on the first row of each variable.
#'   Useful when different statistics within a variable warrant separate
#'   significance testing. Default is \code{FALSE}.
#'
#' @param na_include Logical. If \code{TRUE}, missing values (NAs) are
#'   displayed as a separate category/row for each variable. If \code{FALSE},
#'   missing values are silently excluded from calculations. Default is
#'   \code{FALSE}.
#'
#' @param na_label Character string used to label the missing values row when
#'   \code{na_include = TRUE}. Default is \code{"Unknown"}.
#'
#' @param na_percent Logical. Controls how percentages are calculated for
#'   categorical variables when \code{na_include = TRUE}:
#'   \itemize{
#'     \item If \code{TRUE}, percentages include NAs in the denominator (all
#'       categories sum to 100\%)
#'     \item If \code{FALSE}, percentages exclude NAs from the denominator
#'       (non-missing categories sum to 100\%, missing shown separately)
#'   }
#'   Only affects categorical variables. Default is \code{FALSE}.
#'
#' @param test Logical. If \code{TRUE}, performs appropriate statistical tests
#'   for group comparisons and adds a \emph{p}-value column. Requires
#'   \code{by} to be specified. Tests are automatically selected based on
#'   variable type and test parameters. Default is \code{TRUE}.
#'
#' @param test_continuous Character string specifying the statistical test for
#'   continuous variables:
#'   \itemize{
#'     \item \code{"auto"} - Automatic selection: \emph{t}-test/ANOVA for means,
#'       Wilcoxon/Kruskal-Wallis for medians [default]
#'     \item \code{"t"} - Independent samples \emph{t}-test (2 groups only)
#'     \item \code{"aov"} - One-way ANOVA (2+ groups)
#'     \item \code{"wrs"} - Wilcoxon rank-sum test (2 groups only)
#'     \item \code{"kwt"} - Kruskal-Wallis test (2+ groups)
#'   }
#'
#' @param test_categorical Character string specifying the statistical test
#'   for categorical variables:
#'   \itemize{
#'     \item \code{"auto"} - Automatic selection: Fisher exact test if any
#'       expected cell frequency < 5, otherwise \eqn{\chi^2} test [default]
#'     \item \code{"fisher"} - Fisher exact test
#'     \item \code{"chisq"} - \eqn{\chi^2} test
#'   }
#'
#' @param total Logical or character string controlling the total column:
#'   \itemize{
#'     \item \code{TRUE} or \code{"first"} - Include total column as first
#'       column after Variable/Group [default]
#'     \item \code{"last"} - Include total column as last column before
#'       \emph{p}-value
#'     \item \code{FALSE} - Exclude total column
#'   }
#'
#' @param total_label Character string for the total column header.
#'   Default is \code{"Total"}.
#'
#' @param labels Named character vector or list providing custom display
#'   labels for variables. Names should match variable names (or \code{Surv()}
#'   expressions), values are the display labels. Variables not in
#'   \code{labels} use their original names. Can also label the grouping
#'   variable specified in \code{by}. Default is \code{NULL}.
#'
#' @param number_format Character string or two-element character vector
#'   controlling thousand and decimal separators in formatted output. Named
#'   presets:
#'   \itemize{
#'     \item \code{"us"} - Comma thousands, period decimal: \code{1,234.56} [default]
#'     \item \code{"eu"} - Period thousands, comma decimal: \code{1.234,56}
#'     \item \code{"space"} - Thin-space thousands, period decimal: \code{1 234.56}
#'       (SI/ISO 31-0)
#'     \item \code{"none"} - No thousands separator: \code{1234.56}
#'   }
#'   Or provide a custom two-element vector \code{c(big.mark, decimal.mark)},
#'   \emph{e.g.}, \code{c("'", ".")} for Swiss-style: \verb{1'234.56}.
#'
#'   When \code{NULL} (default), uses
#'   \code{getOption("summata.number_format", "us")}. Set the global option
#'   once per session to avoid passing this argument repeatedly:
#'   \preformatted{
#'     options(summata.number_format = "eu")
#'   }
#'
#' @param ... Additional arguments passed to the underlying statistical test
#'   functions (\emph{e.g.}, \code{var.equal = TRUE} for \emph{t}-tests,
#'   \code{simulate.p.value = TRUE} for Fisher test).
#'
#' @return A data.table with S3 class \code{"desctable"} containing formatted
#'   descriptive statistics. The table structure includes:
#'   \describe{
#'     \item{Variable}{Variable name or label (from \code{labels})}
#'     \item{Group}{For continuous variables: statistic type (\emph{e.g.},
#'       "Mean \eqn{\pm} SD", "Median [IQR]"). For categorical variables:
#'       category level. Empty for variable name rows.}
#'     \item{Total}{Statistics for the total sample (if
#'       \code{total = TRUE})}
#'     \item{Group columns}{Statistics for each group level (when \code{by}
#'       is specified). Column names match group levels.}
#'     \item{\emph{p}-value}{Formatted \emph{p}-values from statistical tests
#'       (when \code{test = TRUE} and \code{by} is specified)}
#'   }
#'
#'   The first row always shows sample sizes for each column. All numeric
#'   output (counts, statistics, \emph{p}-values) respects the
#'   \code{number_format} setting for locale-appropriate formatting.
#'
#'   The returned object includes the following attributes accessible via
#'   \code{attr()}:
#'   \describe{
#'     \item{raw_data}{A data.table containing unformatted numeric values
#'       suitable for further statistical analysis or custom formatting.
#'       Includes additional columns for standard deviations, quartiles,
#'       etc.}
#'     \item{by_variable}{The grouping variable name used (value of
#'       \code{by})}
#'     \item{variables}{The variables analyzed (value of
#'       \code{variables})}
#'   }
#'
#' @details
#' \strong{Variable Type Detection:}
#'
#' The function automatically detects variable types and applies appropriate
#' summaries:
#' \itemize{
#'   \item \strong{Continuous}: Numeric variables (integer or double) receive
#'     statistics specified in \code{stats_continuous}
#'   \item \strong{Categorical}: Character, factor, or logical variables
#'     receive frequency counts and percentages
#'   \item \strong{Time-to-Event}: Variables specified as
#'     \code{Surv(time, event)} display median survival with confidence
#'     intervals (level controlled by \code{conf_level})
#' }
#'
#' \strong{Statistical Testing:}
#'
#' When \code{test = TRUE} and \code{by} is specified:
#' \itemize{
#'   \item \strong{Continuous with "auto"}: Parametric tests (\emph{t}-test, ANOVA)
#'     for mean-based statistics; non-parametric tests (Wilcoxon,
#'     Kruskal-Wallis) for median-based statistics
#'   \item \strong{Categorical with "auto"}: Fisher exact test when any
#'     expected cell frequency < 5; \eqn{\chi^2} test otherwise
#'   \item \strong{Survival}: Log-rank test for comparing survival curves
#'   \item \strong{Range statistics}: No \emph{p}-value computed (ranges
#'     are descriptive)
#' }
#'
#' \strong{Missing Data Handling:}
#'
#' Missing values are handled differently by variable type:
#' \itemize{
#'   \item \strong{Continuous}: NAs excluded from calculations; optionally
#'     shown as count when \code{na_include = TRUE}
#'   \item \strong{Categorical}: NAs can be included as a category when
#'     \code{na_include = TRUE}. The \code{na_percent} parameter controls
#'     whether percentages are calculated with or without NAs in the
#'     denominator
#'   \item \strong{Survival}: NAs in time or event excluded from analysis
#' }
#'
#' \strong{Formatting Conventions:}
#'
#' All numeric output respects the \code{number_format} parameter. Separators
#' within ranges and confidence intervals adapt automatically to avoid
#' ambiguity:
#' \itemize{
#'   \item Mean \eqn{\pm} SD: \code{"45.2 \eqn{\pm} 12.3"} (US) or
#'     \code{"45,2 \eqn{\pm} 12,3"} (EU)
#'   \item Median [IQR]: \code{"38.0 [28.0-52.0]"} (US) or
#'     \code{"38,0 [28,0-52,0]"} (EU, en-dash separator)
#'   \item Range: \code{"18.0-75.0"} (positive, US),
#'     \code{"-5.0 to 10.0"} (when bounds are negative)
#'   \item Survival: \code{"24.5 (21.2-28.9)"} (US) or
#'     \code{"24,5 (21,2-28,9)"} (EU)
#'   \item Counts \eqn{\ge} 1000: \code{"1,234"} (US) or \code{"1.234"} (EU)
#'   \item \emph{p}-values: \code{"< 0.001"} (US) or \code{"< 0,001"} (EU)
#' }
#'
#' @seealso
#' \code{\link{survtable}} for detailed survival summary tables,
#' \code{\link{fit}} for regression modeling,
#' \code{\link{table2pdf}} for PDF export,
#' \code{\link{table2docx}} for Word export,
#' \code{\link{table2html}} for HTML export
#'
#' @examples
#' # Load example clinical trial data
#' data(clintrial)
#'
#' # Example 1: Basic descriptive table without grouping
#' desctable(clintrial,
#'         variables = c("age", "sex", "bmi"))
#'
#'
#' \donttest{
#' # Example 2: Grouped comparison with default tests
#' desctable(clintrial,
#'         by = "treatment",
#'         variables = c("age", "sex", "race", "bmi"))
#'
#' # Example 3: Customize continuous statistics
#' desctable(clintrial,
#'         by = "treatment",
#'         variables = c("age", "bmi", "creatinine"),
#'         stats_continuous = c("median_iqr", "range"))
#'
#' # Example 4: Change categorical display format
#' desctable(clintrial,
#'         by = "treatment",
#'         variables = c("sex", "race", "smoking"),
#'         stats_categorical = "n")  # Show counts only
#'
#' # Example 5: Include missing values
#' desctable(clintrial,
#'         by = "treatment",
#'         variables = c("age", "smoking", "hypertension"),
#'         na_include = TRUE,
#'         na_label = "Missing")
#'
#' # Example 6: Disable statistical testing
#' desctable(clintrial,
#'         by = "treatment",
#'         variables = c("age", "sex", "bmi"),
#'         test = FALSE)
#'
#' # Example 7: Force specific tests
#' desctable(clintrial,
#'         by = "surgery",
#'         variables = c("age", "sex"),
#'         test_continuous = "t",      # t-test instead of auto
#'         test_categorical = "fisher") # Fisher test instead of auto
#'
#' # Example 8: Adjust decimal places
#' desctable(clintrial,
#'         by = "treatment",
#'         variables = c("age", "bmi"),
#'         digits = 2,    # 2 decimals for continuous
#'         p_digits = 4)  # 4 decimals for p-values
#'
#' # Example 9: Custom variable labels
#' labels <- c(
#'     age = "Age (years)",
#'     sex = "Sex",
#'     bmi = "Body Mass Index (kg/m\u00b2)",
#'     treatment = "Treatment Arm"
#' )
#'
#' desctable(clintrial,
#'         by = "treatment",
#'         variables = c("age", "sex", "bmi"),
#'         labels = labels)
#'
#' # Example 10: Position total column last
#' desctable(clintrial,
#'         by = "treatment",
#'         variables = c("age", "sex"),
#'         total = "last")
#'
#' # Example 11: Exclude total column
#' desctable(clintrial,
#'         by = "treatment",
#'         variables = c("age", "sex"),
#'         total = FALSE)
#'
#' # Example 12: Survival analysis
#' desctable(clintrial,
#'         by = "treatment",
#'         variables = "Surv(os_months, os_status)")
#'
#' # Example 13: Multiple survival endpoints
#' desctable(clintrial,
#'         by = "treatment",
#'         variables = c(
#'             "Surv(pfs_months, pfs_status)",
#'             "Surv(os_months, os_status)"
#'         ),
#'         labels = c(
#'             "Surv(pfs_months, pfs_status)" = "Progression-Free Survival",
#'             "Surv(os_months, os_status)" = "Overall Survival"
#'         ))
#'
#' # Example 14: Mixed variable types
#' desctable(clintrial,
#'         by = "treatment",
#'         variables = c(
#'             "age", "sex", "race",           # Demographics
#'             "bmi", "creatinine",            # Labs
#'             "smoking", "hypertension",      # Risk factors
#'             "Surv(os_months, os_status)"    # Survival
#'         ))
#'
#' # Example 15: Export table
#' table1 <- desctable(clintrial,
#'                   by = "treatment",
#'                   variables = c("age", "sex", "bmi"))
#'
#' # Can export directly to PDF/LaTeX/HTML for publication
#' # table2pdf(table1, "table1.pdf")
#' # table2docx(table1, "table1.docx")
#'
#' # Example 16: Three or more groups
#' desctable(clintrial,
#'         by = "stage",  # Assuming stage has 3+ levels
#'         variables = c("age", "sex", "bmi"))
#' # Automatically uses ANOVA/Kruskal-Wallis and chi-squared
#'
#' # Example 17: Access raw unformatted data
#' result <- desctable(clintrial,
#'                   by = "treatment",
#'                   variables = c("age", "bmi"))
#' raw_data <- attr(result, "raw_data")
#' print(raw_data)
#' # Raw data includes unformatted numbers, SDs, quartiles, etc.
#'
#' # Example 18: Check which grouping variable was used
#' result <- desctable(clintrial,
#'                   by = "treatment",
#'                   variables = c("age", "sex"))
#' attr(result, "by_variable")  # "treatment"
#'
#' # Example 19: NA percentage calculation options
#' # Include NAs in percentage denominator (all sum to 100%)
#' desctable(clintrial,
#'         by = "treatment",
#'         variables = "smoking",
#'         na_include = TRUE,
#'         na_percent = TRUE)
#'
#' # Exclude NAs from denominator (non-missing sum to 100%)
#' desctable(clintrial,
#'         by = "treatment",
#'         variables = "smoking",
#'         na_include = TRUE,
#'         na_percent = FALSE)
#'
#' # Example 20: Passing additional test arguments
#' # Equal variance t-test
#' desctable(clintrial,
#'         by = "sex",
#'         variables = "age",
#'         test_continuous = "t",
#'         var.equal = TRUE)
#'
#' # Example 21: European number formatting
#' desctable(clintrial,
#'         by = "treatment",
#'         variables = c("age", "sex", "bmi"),
#'         number_format = "eu")
#'
#' # Example 22: Complete Table 1 for publication
#' table1 <- desctable(
#'     data = clintrial,
#'     by = "treatment",
#'     variables = c(
#'         "age", "sex", "race", "ethnicity", "bmi",
#'         "smoking", "hypertension", "diabetes",
#'         "ecog", "creatinine", "hemoglobin",
#'         "site", "stage", "grade",
#'         "Surv(os_months, os_status)"
#'     ),
#'     labels = clintrial_labels,
#'     stats_continuous = c("median_iqr", "range"),
#'     total = TRUE,
#'     na_include = FALSE
#' )
#' print(table1)
#' }
#' @family descriptive functions
#' @export
desctable <- function(data,
                      by = NULL,
                    variables,
                    stats_continuous = c("median_iqr"),
                    stats_categorical = "n_percent",
                    digits = 1,
                    p_digits = 3,
                    conf_level = 0.95,
                    p_per_stat = FALSE,
                    na_include = FALSE,
                    na_label = "Unknown",
                    na_percent = FALSE,
                    test = TRUE,
                    test_continuous = "auto",
                    test_categorical = "auto",
                    total = TRUE,
                    total_label = "Total",
                    labels = NULL,
                    number_format = NULL,
                    ...) {
    
    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }
    
    ## Resolve number formatting marks
    validate_number_format(number_format)
    marks <- resolve_number_marks(number_format)
    
    ## Set group_var from 'by' parameter
    group_var <- by
    group_var_label <- NULL
    
    ## Apply label to group variable if provided
    if (!is.null(group_var) && !is.null(labels) && group_var %chin% names(labels)) {
        group_var_label <- labels[group_var]
    } else if (!is.null(group_var)) {
        group_var_label <- group_var
    }
    
    ## Variables are already provided as a vector
    vars <- variables
    
    ## Pre-allocate result lists
    n_vars <- length(vars)
    result_list <- vector("list", n_vars)
    raw_result_list <- vector("list", n_vars)
    
    ## Process each variable
    for (i in seq_along(vars)) {
        var_data <- process_variable(
            data = data,
            var = vars[i],
            group_var = group_var,
            stats_continuous = stats_continuous,
            stats_categorical = stats_categorical,
            digits = digits,
            p_digits = p_digits,
            conf_level = conf_level,
            p_per_stat = p_per_stat,
            na_include = na_include,
            na_label = na_label,
            na_percent = na_percent,
            test = test,
            test_continuous = test_continuous,
            test_categorical = test_categorical,
            total = total,
            total_label = total_label,
            labels = labels,
            marks = marks,
            ...
        )
        
        result_list[[i]] <- var_data$formatted
        raw_result_list[[i]] <- var_data$raw
    }
    
    ## Combine all results using rbindlist
    result <- data.table::rbindlist(result_list, fill = TRUE)
    raw_result <- data.table::rbindlist(raw_result_list, fill = TRUE)

    ## Standardize column names
    old_names <- c("variable", "level")
    new_names <- c("Variable", "Group")
    
    if ("variable" %chin% names(result)) {
        data.table::setnames(result, old_names, new_names, skip_absent = TRUE)
        data.table::setnames(raw_result, old_names, new_names, skip_absent = TRUE)
    }

    ## Add p-value column if tests requested (after standardization)
    if (test && !is.null(group_var)) {
        result <- format_pvalues_desctable(result, p_digits, marks)
    }

    ## Add N row as first row (for both grouped and ungrouped tables)
    if (!is.null(group_var)) {
        ## Get the group values in the correct order
        if (is.factor(data[[group_var]])) {
            ## Use factor levels for proper ordering
            groups <- levels(data[[group_var]])
        } else {
            ## Fall back to unique values for non-factors
            groups <- unique(data[[group_var]])
            groups <- groups[!is.na(groups)]
        }
        
        ## Create N row
        n_row <- data.table::data.table(
                                 Variable = "N",
                                 Group = ""
                             )
        
        ## Calculate and add total if present
        if (total_label %chin% names(result)) {
            n_total <- nrow(data)
            n_row[[total_label]] <- format_count(n_total, marks)
        }
        
        ## Calculate for each group in the correct order
        for (grp in groups) {
            grp_col <- as.character(grp)
            if (grp_col %chin% names(result)) {
                n_group <- sum(data[[group_var]] == grp, na.rm = TRUE)
                n_row[[grp_col]] <- format_count(n_group, marks)
            }
        }
        
        ## Add empty p-value column if it exists
        if ("p-value" %chin% names(result)) {
            n_row[["p-value"]] <- ""
        }
        
        ## Prepend N row
        result <- rbind(n_row, result, fill = TRUE)
        
    } else if (total && total_label %chin% names(result)) {
        ## Add N row for ungrouped tables with Total column
        n_total <- nrow(data)
        n_row <- data.table::data.table(
                                 Variable = "N",
                                 Group = ""
                             )
        n_row[[total_label]] <- format_count(n_total, marks)
        
        ## Prepend N row
        result <- rbind(n_row, result, fill = TRUE)
    }

    ## Reorder columns if total position specified
    if (!isFALSE(total) && !is.null(group_var)) {
        result <- reorder_total_column(result, total, total_label)
    }

    ## Attach raw data and metadata as attributes
    data.table::setattr(result, "raw_data", raw_result)
    data.table::setattr(result, "by_variable", group_var)
    data.table::setattr(result, "variables", variables)

    result[]
    return(result)
}
