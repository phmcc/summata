# Save a Table to File

Writes a table to file in the format indicated by the file extension,
dispatching to the appropriate specialized export function. Provides a
unified interface for table export across all supported formats, and is
the counterpart to
[`forestsave()`](https://phmcc.codeberg.page/summata/reference/forestsave.md)
for forest plots.

## Usage

``` r
tablesave(table, file, quiet = FALSE, ...)

autotable(table, file, quiet = FALSE, ...)
```

## Arguments

- table:

  A table produced by
  [`desctable()`](https://phmcc.codeberg.page/summata/reference/desctable.md),
  [`survtable()`](https://phmcc.codeberg.page/summata/reference/survtable.md),
  [`fit()`](https://phmcc.codeberg.page/summata/reference/fit.md),
  [`uniscreen()`](https://phmcc.codeberg.page/summata/reference/uniscreen.md),
  [`fullfit()`](https://phmcc.codeberg.page/summata/reference/fullfit.md),
  [`compfit()`](https://phmcc.codeberg.page/summata/reference/compfit.md),
  or
  [`multifit()`](https://phmcc.codeberg.page/summata/reference/multifit.md).
  Any data frame, data.table, or matrix is accepted.

- file:

  Character string specifying the output filename. The file extension
  determines the export format:

  - `.csv` - Comma-separated values (uses
    [`data.table::fwrite()`](https://rdrr.io/pkg/data.table/man/fwrite.html))

  - `.tsv` - Tab-separated values (uses
    [`data.table::fwrite()`](https://rdrr.io/pkg/data.table/man/fwrite.html))

  - `.pdf` - PDF via LaTeX (uses
    [`table2pdf()`](https://phmcc.codeberg.page/summata/reference/table2pdf.md))

  - `.docx` - Microsoft Word (uses
    [`table2docx()`](https://phmcc.codeberg.page/summata/reference/table2docx.md))

  - `.html` or `.htm` - HTML (uses
    [`table2html()`](https://phmcc.codeberg.page/summata/reference/table2html.md))

  - `.pptx` - Microsoft PowerPoint (uses
    [`table2pptx()`](https://phmcc.codeberg.page/summata/reference/table2pptx.md))

  - `.tex` - LaTeX source (uses
    [`table2tex()`](https://phmcc.codeberg.page/summata/reference/table2tex.md))

  - `.rtf` - Rich Text Format (uses
    [`table2rtf()`](https://phmcc.codeberg.page/summata/reference/table2rtf.md))

- quiet:

  Logical. Suppress progress and confirmation messages. The setting is
  forwarded to the format-specific export function, so the LaTeX
  compilation notices from
  [`table2pdf()`](https://phmcc.codeberg.page/summata/reference/table2pdf.md)
  are suppressed with it. Default is `FALSE`.

- ...:

  Additional arguments passed to the format-specific function. See the
  documentation for individual functions for available parameters:

  PDF

  :   [`table2pdf()`](https://phmcc.codeberg.page/summata/reference/table2pdf.md) -
      `orientation`, `paper`, `margins`, `fit_to_page`, *etc.*

  DOCX

  :   [`table2docx()`](https://phmcc.codeberg.page/summata/reference/table2docx.md) -
      `font_size`, `font_family`, `caption`, *etc.*

  HTML

  :   [`table2html()`](https://phmcc.codeberg.page/summata/reference/table2html.md) -
      `format_headers`, `zebra_stripes`, *etc.*

  PPTX

  :   [`table2pptx()`](https://phmcc.codeberg.page/summata/reference/table2pptx.md) -
      `font_size`, `font_family`, `caption`, *etc.*

  TEX

  :   [`table2tex()`](https://phmcc.codeberg.page/summata/reference/table2tex.md) -
      `caption`, `format_headers`, `align`, *etc.*

  RTF

  :   [`table2rtf()`](https://phmcc.codeberg.page/summata/reference/table2rtf.md) -
      `font_size`, `font_family`, `caption`, *etc.*

  CSV, TSV

  :   [`data.table::fwrite()`](https://rdrr.io/pkg/data.table/man/fwrite.html) -
      `sep`, `quote`, `na`, *etc.*

  Common parameters across formats include:

  `caption`

  :   Table caption (supported by most formats)

  `font_size`

  :   Base font size in points (PDF, DOCX, PPTX, RTF)

  `format_headers`

  :   Format column headers (all formats)

  `bold_significant`

  :   Bold significant *p*-values (all formats)

  `p_threshold`

  :   Threshold for *p*-value bolding (all formats)

  `indent_groups`

  :   Indent factor levels (all formats)

  `condense_table`

  :   Condense to essential rows (all formats)

  `zebra_stripes`

  :   Alternating background colors (most formats)

## Value

Invisibly returns `file`.

## Details

This function provides a convenient wrapper around format-specific
export functions, automatically routing to the appropriate function
based on the file extension. All parameters are passed through to the
underlying function, so the full range of format-specific options
remains available.

For format-specific advanced features, individual export functions may
be called directly:

- PDF exports support orientation, paper size, margins, and auto-sizing

- DOCX/PPTX/RTF support font customization and flextable formatting

- HTML supports CSS styling, responsive design, and custom themes

- TeX generates standalone LaTeX source with booktabs styling

## See also

[`table2pdf`](https://phmcc.codeberg.page/summata/reference/table2pdf.md),
[`table2docx`](https://phmcc.codeberg.page/summata/reference/table2docx.md),
[`table2pptx`](https://phmcc.codeberg.page/summata/reference/table2pptx.md),
[`table2html`](https://phmcc.codeberg.page/summata/reference/table2html.md),
[`table2rtf`](https://phmcc.codeberg.page/summata/reference/table2rtf.md),
[`table2tex`](https://phmcc.codeberg.page/summata/reference/table2tex.md)

Other export functions:
[`table2docx()`](https://phmcc.codeberg.page/summata/reference/table2docx.md),
[`table2html()`](https://phmcc.codeberg.page/summata/reference/table2html.md),
[`table2pdf()`](https://phmcc.codeberg.page/summata/reference/table2pdf.md),
[`table2pptx()`](https://phmcc.codeberg.page/summata/reference/table2pptx.md),
[`table2rtf()`](https://phmcc.codeberg.page/summata/reference/table2rtf.md),
[`table2tex()`](https://phmcc.codeberg.page/summata/reference/table2tex.md)

## Examples

``` r
# Create example data
data(clintrial)
data(clintrial_labels)
tbl <- desctable(clintrial, by = "treatment",
    variables = c("age", "sex"), labels = clintrial_labels)

# Example 1: The format follows the file extension
if (requireNamespace("xtable", quietly = TRUE)) {
  tablesave(tbl, file.path(tempdir(), "example.html"))
}
#> Table saved to /tmp/RtmpPQeTCd/example.html

# Example 2: Delimited output needs no additional packages
tablesave(tbl, file.path(tempdir(), "example.csv"))
#> Table saved to /tmp/RtmpPQeTCd/example.csv

# Example 3: Suppress the message reporting the file written
tablesave(tbl, file.path(tempdir(), "example.tsv"), quiet = TRUE)

# \donttest{

# Load example data
data(clintrial)
data(clintrial_labels)

# Create a regression table
results <- fit(
    data = clintrial,
    outcome = "os_status",
    predictors = c("age", "sex", "treatment"),
    labels = clintrial_labels
)

# Test that LaTeX can actually compile (needed for PDF export)
has_latex <- local({
  if (!nzchar(Sys.which("pdflatex"))) return(FALSE)
  test_tex <- file.path(tempdir(), "summata_latex_test.tex")
  writeLines(c("\\documentclass{article}",
               "\\usepackage{booktabs}",
               "\\begin{document}", "test",
               "\\end{document}"), test_tex)
  result <- tryCatch(
    system2("pdflatex", c("-interaction=nonstopmode",
            paste0("-output-directory=", tempdir()), test_tex),
            stdout = FALSE, stderr = FALSE),
    error = function(e) 1L)
  result == 0L
})

# Example 4: The format is taken from the file extension
tablesave(results, file.path(tempdir(), "results.html"))  # Creates HTML file
#> Table saved to /tmp/RtmpPQeTCd/results.html
tablesave(results, file.path(tempdir(), "results.docx"))  # Creates Word document
#> Table saved to /tmp/RtmpPQeTCd/results.docx
tablesave(results, file.path(tempdir(), "results.pptx"))  # Creates PowerPoint slide
#> Table saved to /tmp/RtmpPQeTCd/results.pptx
tablesave(results, file.path(tempdir(), "results.tex"))   # Creates LaTeX source
#> Table saved to /tmp/RtmpPQeTCd/results.tex
tablesave(results, file.path(tempdir(), "results.rtf"))   # Creates RTF document
#> Table saved to /tmp/RtmpPQeTCd/results.rtf
if (has_latex) {
  tablesave(results, file.path(tempdir(), "results.pdf")) # Creates PDF
}
#> Compiling PDF...
#> Table saved to /tmp/RtmpPQeTCd/results.pdf

# Example 5: Format-specific parameters are passed through
if (has_latex) {
  tablesave(results, file.path(tempdir(), "results.pdf"), 
             orientation = "landscape",
             paper = "a4",
             font_size = 10)
}
#> Compiling PDF...
#> Table saved to /tmp/RtmpPQeTCd/results.pdf

tablesave(results, file.path(tempdir(), "results.docx"),
           caption = "Table 1: Logistic Regression Results",
           font_family = "Times New Roman",
           condense_table = TRUE)
#> Table saved to /tmp/RtmpPQeTCd/results.docx

tablesave(results, file.path(tempdir(), "results.html"),
           zebra_stripes = TRUE,
           dark_header = TRUE,
           bold_significant = TRUE)
#> Table saved to /tmp/RtmpPQeTCd/results.html

# Example 6: Any summata table may be saved
desc <- desctable(clintrial,
                  by = "treatment",
                  variables = c("age", "sex", "bmi"))
if (has_latex) {
  tablesave(desc, file.path(tempdir(), "demographics.pdf"))
}
#> Compiling PDF...
#> Table saved to /tmp/RtmpPQeTCd/demographics.pdf

# Example 7: Model comparison table
# Information criteria assume a common sample, so the candidate
# predictors are restricted to complete cases before comparison
comparison_data <- na.omit(
    clintrial[, c("os_status", "age", "sex", "treatment", "stage")]
)
comparison <- compfit(
    data = comparison_data,
    outcome = "os_status",
    model_list = list(
        base = c("age", "sex"),
        full = c("age", "sex", "treatment", "stage")
    )
)
#> Auto-detected binary outcome, using logistic regression
#> Fitting base with 2 predictors...
#> Fitting full with 4 predictors...
tablesave(comparison, file.path(tempdir(), "model_comparison.docx"))
#> Table saved to /tmp/RtmpPQeTCd/model_comparison.docx

# }
```
