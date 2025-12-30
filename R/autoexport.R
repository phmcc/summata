#' Export Table with Automatic Format Detection
#'
#' Automatically detects the output format based on file extension and exports
#' the table using the appropriate specialized function (table2pdf, table2docx, 
#' table2pptx, table2html, table2rtf, or table2tex). Provides a unified interface for
#' table export across all supported formats.
#'
#' @param table A data.frame, data.table, or matrix to export. Can be output from 
#'   \code{\link{desctable}}, \code{\link{fit}}, \code{\link{uniscreen}}, 
#'   \code{\link{fullfit}}, \code{\link{compfit}}, or any tabular data structure.
#'   
#' @param file Character string specifying the output filename. The file extension
#'   determines the export format:
#'   \itemize{
#'     \item \code{.pdf} - PDF via LaTeX (uses \code{table2pdf})
#'     \item \code{.docx} - Microsoft Word (uses \code{table2docx})
#'     \item \code{.pptx} - Microsoft PowerPoint (uses \code{table2pptx})
#'     \item \code{.html} or \code{.htm} - HTML (uses \code{table2html})
#'     \item \code{.rtf} - Rich Text Format (uses \code{table2rtf})
#'     \item \code{.tex} - LaTeX source (uses \code{table2tex})
#'   }
#'   
#' @param ... Additional arguments passed to the format-specific function. See
#'   the documentation for individual functions for available parameters:
#'   \itemize{
#'     \item PDF: \code{\link{table2pdf}} - orientation, paper, margins, fit_to_page, etc.
#'     \item DOCX: \code{\link{table2docx}} - font_size, font_family, caption, etc.
#'     \item PPTX: \code{\link{table2pptx}} - font_size, font_family, caption, etc.
#'     \item HTML: \code{\link{table2html}} - format_headers, zebra_stripes, etc.
#'     \item RTF: \code{\link{table2rtf}} - font_size, font_family, caption, etc.
#'     \item TEX: \code{\link{table2tex}} - caption, format_headers, align, etc.
#'   }
#'   
#'   Common parameters across formats include:
#'   \itemize{
#'     \item \code{caption} - Table caption (supported by most formats)
#'     \item \code{font_size} - Base font size in points (PDF, DOCX, PPTX, RTF)
#'     \item \code{format_headers} - Format column headers (all formats)
#'     \item \code{bold_significant} - Bold significant p-values (all formats)
#'     \item \code{p_threshold} - P-value threshold for bolding (all formats)
#'     \item \code{indent_groups} - Indent factor levels (all formats)
#'     \item \code{condense_table} - Condense to essential rows (all formats)
#'     \item \code{zebra_stripes} - Alternating background colors (most formats)
#'   }
#'
#' @return Invisibly returns the file path. Called primarily for its side effect
#'   of creating the output file.
#'
#' @details
#' This function provides a convenient wrapper around format-specific export 
#' functions, automatically routing to the appropriate function based on the 
#' file extension. All parameters are passed through to the underlying function,
#' so the full range of format-specific options remains available.
#' 
#' For format-specific advanced features, you may prefer to use the individual
#' export functions directly:
#' \itemize{
#'   \item PDF exports support orientation, paper size, margins, and auto-sizing
#'   \item DOCX/PPTX/RTF support font customization and flextable formatting
#'   \item HTML supports CSS styling, responsive design, and custom themes
#'   \item TEX generates standalone LaTeX source with booktabs styling
#' }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(clintrial)
#' 
#' # Create a regression table
#' mod <- glm(outcome ~ age + gender + treatment, 
#'            family = binomial, 
#'            data = clintrial)
#' table <- fit(mod, data = clintrial)
#' 
#' # Export automatically detects format from extension
#' autoexport(table, "results.pdf")                 # Creates PDF
#' autoexport(table, "results.docx")                # Creates Word document
#' autoexport(table, "results.pptx")                # Creates PowerPoint slide
#' autoexport(table, "results.html")                # Creates HTML file
#' autoexport(table, "results.rtf")                 # Creates RTF document
#' autoexport(table, "results.tex")                 # Creates LaTeX source
#' 
#' # Pass format-specific parameters
#' autoexport(table, "results.pdf", 
#'            orientation = "landscape",
#'            paper = "a4",
#'            font_size = 10)
#' 
#' autoexport(table, "results.docx",
#'            caption = "Table 1: Logistic Regression Results",
#'            font_family = "Times New Roman",
#'            condense_table = TRUE)
#' 
#' autoexport(table, "results.html",
#'            zebra_stripes = TRUE,
#'            style = "modern",
#'            bold_significant = TRUE)
#' 
#' # Works with any tabular output
#' desc <- desctable(clintrial, strata = "treatment")
#' autoexport(desc, "demographics.pdf")
#' 
#' comp <- compfit(
#'   fit(glm(outcome ~ age, family = binomial, data = clintrial)),
#'   fit(glm(outcome ~ age + gender, family = binomial, data = clintrial))
#' )
#' autoexport(comp, "model_comparison.docx")
#' }
#'
#' @seealso
#' \code{\link{table2pdf}}, \code{\link{table2docx}}, \code{\link{table2pptx}},
#' \code{\link{table2html}}, \code{\link{table2rtf}}, \code{\link{table2tex}}
#'
#' @export
autoexport <- function(table, file, ...) {
    
    ## Validate inputs
    if (missing(table) || missing(file)) {
        stop("Both 'table' and 'file' arguments are required")
    }
    
    if (!is.character(file) || length(file) != 1) {
        stop("'file' must be a single character string")
    }
    
    ## Extract file extension
    file_ext <- tolower(tools::file_ext(file))
    
    ## Check if extension is empty
    if (file_ext == "") {
        stop(paste0("File '", file, "' has no extension. ",
                    "Please specify one of: .pdf, .docx, .pptx, .html, .htm, .rtf, .tex"))
    }
    
    ## Route to appropriate export function based on extension
    result <- switch(file_ext,
                     "pdf" = {
                         table2pdf(table = table, file = file, ...)
                     },
                     "docx" = {
                         table2docx(table = table, file = file, ...)
                     },
                     "pptx" = {
                         table2pptx(table = table, file = file, ...)
                     },
                     "html" = ,
                     "htm" = {
                         table2html(table = table, file = file, ...)
                     },
                     "rtf" = {
                         table2rtf(table = table, file = file, ...)
                     },
                     "tex" = {
                         table2tex(table = table, file = file, ...)
                     },
                     {
                         ## Unknown extension
                         stop(paste0("Unsupported file extension: '.", file_ext, "'\n",
                                     "Supported formats are:\n",
                                     "  .pdf   - PDF via LaTeX (table2pdf)\n",
                                     "  .docx  - Microsoft Word (table2docx)\n",
                                     "  .pptx  - Microsoft PowerPoint (table2pptx)\n",
                                     "  .html  - HTML (table2html)\n",
                                     "  .htm   - HTML (table2html)\n",
                                     "  .rtf   - Rich Text Format (table2rtf)\n",
                                     "  .tex   - LaTeX source (table2tex)"))
                     }
                     )
    
    ## Return file path invisibly (matches pattern of export functions)
    invisible(file)
}
