### * Main function

#' Save a Table to File
#'
#' Writes a table to file in the format indicated by the file extension,
#' dispatching to the appropriate specialized export function. Provides a
#' unified interface for table export across all supported formats, and is the
#' counterpart to \code{forestsave()} for forest plots.
#'
#' @param table A table produced by \code{desctable()}, \code{survtable()},
#'   \code{fit()}, \code{uniscreen()}, \code{fullfit()}, \code{compfit()}, or
#'   \code{multifit()}. Any data frame, data.table, or matrix is accepted.
#'   
#' @param file Character string specifying the output filename. The file extension
#'   determines the export format:
#'   \itemize{
#'     \item \code{.csv} - Comma-separated values (uses \code{data.table::fwrite()})
#'     \item \code{.tsv} - Tab-separated values (uses \code{data.table::fwrite()})
#'     \item \code{.pdf} - PDF via LaTeX (uses \code{table2pdf()})
#'     \item \code{.docx} - Microsoft Word (uses \code{table2docx()})
#'     \item \code{.html} or \code{.htm} - HTML (uses \code{table2html()})
#'     \item \code{.pptx} - Microsoft PowerPoint (uses \code{table2pptx()})
#'     \item \code{.tex} - LaTeX source (uses \code{table2tex()})
#'     \item \code{.rtf} - Rich Text Format (uses \code{table2rtf()})
#'   }
#'   
#' @param quiet Logical. Suppress progress and confirmation messages. The
#'   setting is forwarded to the format-specific export function, so the
#'   LaTeX compilation notices from \code{table2pdf()} are suppressed with
#'   it. Default is \code{FALSE}.
#'
#' @param ... Additional arguments passed to the format-specific function. See
#'   the documentation for individual functions for available parameters:
#'   \describe{
#'     \item{PDF}{\code{table2pdf()} - \code{orientation}, \code{paper}, \code{margins},
#'                \code{fit_to_page}, \emph{etc.}}
#'     \item{DOCX}{\code{table2docx()} - \code{font_size}, \code{font_family},
#'                \code{caption}, \emph{etc.}}
#'     \item{HTML}{\code{table2html()} - \code{format_headers}, \code{zebra_stripes},
#'                \emph{etc.}}
#'     \item{PPTX}{\code{table2pptx()} - \code{font_size}, \code{font_family},
#'                \code{caption}, \emph{etc.}}
#'     \item{TEX}{\code{table2tex()} - \code{caption}, \code{format_headers},
#'                \code{align}, \emph{etc.}}
#'     \item{RTF}{\code{table2rtf()} - \code{font_size}, \code{font_family},
#'                \code{caption}, \emph{etc.}}
#'     \item{CSV, TSV}{\code{data.table::fwrite()} - \code{sep}, \code{quote},
#'                \code{na}, \emph{etc.}}
#'   }
#'   
#'   Common parameters across formats include:
#'   \describe{
#'     \item{\code{caption}}{Table caption (supported by most formats)}
#'     \item{\code{font_size}}{Base font size in points (PDF, DOCX, PPTX, RTF)}
#'     \item{\code{format_headers}}{Format column headers (all formats)}
#'     \item{\code{bold_significant}}{Bold significant \emph{p}-values (all formats)}
#'     \item{\code{p_threshold}}{Threshold for \emph{p}-value bolding (all formats)}
#'     \item{\code{indent_groups}}{Indent factor levels (all formats)}
#'     \item{\code{condense_table}}{Condense to essential rows (all formats)}
#'     \item{\code{zebra_stripes}}{Alternating background colors (most formats)}
#'   }
#'
#' @return Invisibly returns \code{file}.
#'
#' @details
#' This function provides a convenient wrapper around format-specific export 
#' functions, automatically routing to the appropriate function based on the 
#' file extension. All parameters are passed through to the underlying function,
#' so the full range of format-specific options remains available.
#' 
#' For format-specific advanced features, individual export functions may be
#' called directly:
#' \itemize{
#'   \item PDF exports support orientation, paper size, margins, and auto-sizing
#'   \item DOCX/PPTX/RTF support font customization and \pkg{flextable} formatting
#'   \item HTML supports CSS styling, responsive design, and custom themes
#'   \item TeX generates standalone LaTeX source with booktabs styling
#' }
#'
#' @examples
#' # Create example data
#' data(clintrial)
#' data(clintrial_labels)
#' tbl <- desctable(clintrial, by = "treatment",
#'     variables = c("age", "sex"), labels = clintrial_labels)
#'
#' # Example 1: The format follows the file extension
#' if (requireNamespace("xtable", quietly = TRUE)) {
#'   tablesave(tbl, file.path(tempdir(), "example.html"))
#' }
#'
#' # Example 2: Delimited output needs no additional packages
#' tablesave(tbl, file.path(tempdir(), "example.csv"))
#'
#' # Example 3: Suppress the message reporting the file written
#' tablesave(tbl, file.path(tempdir(), "example.tsv"), quiet = TRUE)
#'
#' \donttest{
#' # Load example data
#' data(clintrial)
#' data(clintrial_labels)
#' 
#' # Create a regression table
#' results <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment"),
#'     labels = clintrial_labels
#' )
#' 
#' # Test that LaTeX can actually compile (needed for PDF export)
#' has_latex <- local({
#'   if (!nzchar(Sys.which("pdflatex"))) return(FALSE)
#'   test_tex <- file.path(tempdir(), "summata_latex_test.tex")
#'   writeLines(c("\\documentclass{article}",
#'                "\\usepackage{booktabs}",
#'                "\\begin{document}", "test",
#'                "\\end{document}"), test_tex)
#'   result <- tryCatch(
#'     system2("pdflatex", c("-interaction=nonstopmode",
#'             paste0("-output-directory=", tempdir()), test_tex),
#'             stdout = FALSE, stderr = FALSE),
#'     error = function(e) 1L)
#'   result == 0L
#' })
#' 
#' # Example 4: The format is taken from the file extension
#' tablesave(results, file.path(tempdir(), "results.html"))  # Creates HTML file
#' tablesave(results, file.path(tempdir(), "results.docx"))  # Creates Word document
#' tablesave(results, file.path(tempdir(), "results.pptx"))  # Creates PowerPoint slide
#' tablesave(results, file.path(tempdir(), "results.tex"))   # Creates LaTeX source
#' tablesave(results, file.path(tempdir(), "results.rtf"))   # Creates RTF document
#' if (has_latex) {
#'   tablesave(results, file.path(tempdir(), "results.pdf")) # Creates PDF
#' }
#' 
#' # Example 5: Format-specific parameters are passed through
#' if (has_latex) {
#'   tablesave(results, file.path(tempdir(), "results.pdf"), 
#'              orientation = "landscape",
#'              paper = "a4",
#'              font_size = 10)
#' }
#' 
#' tablesave(results, file.path(tempdir(), "results.docx"),
#'            caption = "Table 1: Logistic Regression Results",
#'            font_family = "Times New Roman",
#'            condense_table = TRUE)
#' 
#' tablesave(results, file.path(tempdir(), "results.html"),
#'            zebra_stripes = TRUE,
#'            dark_header = TRUE,
#'            bold_significant = TRUE)
#' 
#' # Example 6: Any summata table may be saved
#' desc <- desctable(clintrial,
#'                   by = "treatment",
#'                   variables = c("age", "sex", "bmi"))
#' if (has_latex) {
#'   tablesave(desc, file.path(tempdir(), "demographics.pdf"))
#' }
#' 
#' # Example 7: Model comparison table
#' # Information criteria assume a common sample, so the candidate
#' # predictors are restricted to complete cases before comparison
#' comparison_data <- na.omit(
#'     clintrial[, c("os_status", "age", "sex", "treatment", "stage")]
#' )
#' comparison <- compfit(
#'     data = comparison_data,
#'     outcome = "os_status",
#'     model_list = list(
#'         base = c("age", "sex"),
#'         full = c("age", "sex", "treatment", "stage")
#'     )
#' )
#' tablesave(comparison, file.path(tempdir(), "model_comparison.docx"))
#'
#' }
#'
#' @seealso
#' \code{\link{table2pdf}}, \code{\link{table2docx}}, \code{\link{table2pptx}},
#' \code{\link{table2html}}, \code{\link{table2rtf}}, \code{\link{table2tex}}
#'
#' @family export functions
#' @export
tablesave <- function(table, file, quiet = FALSE, ...) {

    if (missing(table) || missing(file)) {
        stop("Both 'table' and 'file' arguments are required.", call. = FALSE)
    }

    if (!is.data.frame(table) && !is.matrix(table)) {
        stop("'table' must be a data frame, data.table, or matrix.",
             call. = FALSE)
    }

    if (!is.character(file) || length(file) != 1 || !nzchar(file)) {
        stop("'file' must be a single non-empty character string.",
             call. = FALSE)
    }

    file_ext <- tolower(tools::file_ext(file))

    supported <- paste0(
        "Supported formats are:\n",
        "  .pdf   - PDF via LaTeX\n",
        "  .docx  - Microsoft Word\n",
        "  .pptx  - Microsoft PowerPoint\n",
        "  .html  - HTML\n",
        "  .htm   - HTML\n",
        "  .rtf   - Rich Text Format\n",
        "  .tex   - LaTeX source\n",
        "  .csv   - Comma-separated values\n",
        "  .tsv   - Tab-separated values")

    if (file_ext == "") {
        stop("File '", file, "' has no extension, so no format can be ",
             "determined.\n", supported, call. = FALSE)
    }

    ## The delimited formats are written here rather than delegated, as the
    ## package has no table2csv() to dispatch to. data.table::fwrite() is used
    ## for consistency with the export routes documented elsewhere.
    write_delimited <- function(sep) {
        data.table::fwrite(table, file = file, sep = sep, ...)
        if (!quiet) message(sprintf("Table saved to %s", file))
    }

    ## The export functions accept quiet themselves, so it is forwarded rather
    ## than applied by suppressing messages here
    switch(
        file_ext,
        "pdf"  = table2pdf(table = table, file = file, quiet = quiet, ...),
        "docx" = table2docx(table = table, file = file, quiet = quiet, ...),
        "pptx" = table2pptx(table = table, file = file, quiet = quiet, ...),
        "html" = ,
        "htm"  = table2html(table = table, file = file, quiet = quiet, ...),
        "rtf"  = table2rtf(table = table, file = file, quiet = quiet, ...),
        "tex"  = table2tex(table = table, file = file, quiet = quiet, ...),
        "csv"  = write_delimited(","),
        "tsv"  = write_delimited("\t"),
        stop("Unsupported file extension: '.", file_ext, "'\n", supported,
             call. = FALSE)
    )

    invisible(file)
}



#' @rdname tablesave
#' @export
autotable <- function(table, file, quiet = FALSE, ...) {

    .Deprecated("tablesave", package = "summata",
                msg = paste0("'autotable()' is deprecated and will be removed ",
                             "in a future release.\n",
                             "Use 'tablesave()' instead, which is the same ",
                             "function under a name that pairs with ",
                             "'forestsave()'."))

    tablesave(table = table, file = file, quiet = quiet, ...)
}
