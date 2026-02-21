# Export Table with Automatic Format Detection

Automatically detects the output format based on file extension and
exports the table using the appropriate specialized function. Provides a
unified interface for table export across all supported formats.

## Usage

``` r
autotable(table, file, ...)
```

## Arguments

- table:

  Data frame, data.table, or matrix to export. Can be output from
  [`desctable()`](https://phmcc.github.io/summata/reference/desctable.md),
  [`survtable()`](https://phmcc.github.io/summata/reference/survtable.md),
  [`fit()`](https://phmcc.github.io/summata/reference/fit.md),
  [`uniscreen()`](https://phmcc.github.io/summata/reference/uniscreen.md),
  [`fullfit()`](https://phmcc.github.io/summata/reference/fullfit.md),
  [`compfit()`](https://phmcc.github.io/summata/reference/compfit.md),
  [`multifit()`](https://phmcc.github.io/summata/reference/multifit.md),
  or any tabular data structure.

- file:

  Character string specifying the output filename. The file extension
  determines the export format:

  - `.pdf` - PDF via LaTeX (uses
    [`table2pdf()`](https://phmcc.github.io/summata/reference/table2pdf.md))

  - `.docx` - Microsoft Word (uses
    [`table2docx()`](https://phmcc.github.io/summata/reference/table2docx.md))

  - `.html` or `.htm` - HTML (uses
    [`table2html()`](https://phmcc.github.io/summata/reference/table2html.md))

  - `.pptx` - Microsoft PowerPoint (uses
    [`table2pptx()`](https://phmcc.github.io/summata/reference/table2pptx.md))

  - `.tex` - LaTeX source (uses
    [`table2tex()`](https://phmcc.github.io/summata/reference/table2tex.md))

  - `.rtf` - Rich Text Format (uses
    [`table2rtf()`](https://phmcc.github.io/summata/reference/table2rtf.md))

- ...:

  Additional arguments passed to the format-specific function. See the
  documentation for individual functions for available parameters:

  PDF

  :   [`table2pdf()`](https://phmcc.github.io/summata/reference/table2pdf.md) -
      `orientation`, `paper`, `margins`, `fit_to_page`, *etc.*

  DOCX

  :   [`table2docx()`](https://phmcc.github.io/summata/reference/table2docx.md) -
      `font_size`, `font_family`, `caption`, *etc.*

  HTML

  :   [`table2html()`](https://phmcc.github.io/summata/reference/table2html.md) -
      `format_headers`, `zebra_stripes`, *etc.*

  PPTX

  :   [`table2pptx()`](https://phmcc.github.io/summata/reference/table2pptx.md) -
      `font_size`, `font_family`, `caption`, *etc.*

  TEX

  :   [`table2tex()`](https://phmcc.github.io/summata/reference/table2tex.md) -
      `caption`, `format_headers`, `align`, *etc.*

  RTF

  :   [`table2rtf()`](https://phmcc.github.io/summata/reference/table2rtf.md) -
      `font_size`, `font_family`, `caption`, *etc.*

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

Invisibly returns the file path. Called primarily for its side effect of
creating the output file.

## Details

This function provides a convenient wrapper around format-specific
export functions, automatically routing to the appropriate function
based on the file extension. All parameters are passed through to the
underlying function, so the full range of format-specific options
remains available.

For format-specific advanced features, you may prefer to use the
individual export functions directly:

- PDF exports support orientation, paper size, margins, and auto-sizing

- DOCX/PPTX/RTF support font customization and flextable formatting

- HTML supports CSS styling, responsive design, and custom themes

- TeX generates standalone LaTeX source with booktabs styling

## See also

[`table2pdf`](https://phmcc.github.io/summata/reference/table2pdf.md),
[`table2docx`](https://phmcc.github.io/summata/reference/table2docx.md),
[`table2pptx`](https://phmcc.github.io/summata/reference/table2pptx.md),
[`table2html`](https://phmcc.github.io/summata/reference/table2html.md),
[`table2rtf`](https://phmcc.github.io/summata/reference/table2rtf.md),
[`table2tex`](https://phmcc.github.io/summata/reference/table2tex.md)

Other export functions:
[`table2docx()`](https://phmcc.github.io/summata/reference/table2docx.md),
[`table2html()`](https://phmcc.github.io/summata/reference/table2html.md),
[`table2pdf()`](https://phmcc.github.io/summata/reference/table2pdf.md),
[`table2pptx()`](https://phmcc.github.io/summata/reference/table2pptx.md),
[`table2rtf()`](https://phmcc.github.io/summata/reference/table2rtf.md),
[`table2tex()`](https://phmcc.github.io/summata/reference/table2tex.md)

## Examples

``` r
if (FALSE) {
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

# Export automatically detects format from extension
autotable(results, "results.pdf")   # Creates PDF
autotable(results, "results.docx")  # Creates Word document
autotable(results, "results.html")  # Creates HTML file
autotable(results, "results.pptx")  # Creates PowerPoint slide
autotable(results, "results.tex")   # Creates LaTeX source
autotable(results, "results.rtf")   # Creates RTF document

# Pass format-specific parameters
autotable(results, "results.pdf", 
           orientation = "landscape",
           paper = "a4",
           font_size = 10)

autotable(results, "results.docx",
           caption = "Table 1: Logistic Regression Results",
           font_family = "Times New Roman",
           condense_table = TRUE)

autotable(results, "results.html",
           zebra_stripes = TRUE,
           dark_header = TRUE,
           bold_significant = TRUE)

# Works with any summata table output
desc <- desctable(clintrial,
                  by = "treatment",
                  variables = c("age", "sex", "bmi"))
autotable(desc, "demographics.pdf")

comparison <- compfit(
    data = clintrial,
    outcome = "os_status",
    model_list = list(
        base = c("age", "sex"),
        full = c("age", "sex", "treatment", "stage")
    )
)
autotable(comparison, "model_comparison.docx")
}
```
