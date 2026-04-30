# Export Table to PDF Format

Converts a data frame, data.table, or matrix to a professionally
formatted PDF document using LaTeX as an intermediate format. Provides
extensive control over page layout, typography, and formatting for
publication-ready output. Particularly well-suited for tables from
regression analyses, descriptive statistics, and model comparisons.
Requires xtable for export.

## Usage

``` r
table2pdf(
  table,
  file,
  orientation = "portrait",
  paper = "letter",
  margins = NULL,
  fit_to_page = TRUE,
  font_size = 8,
  caption = NULL,
  caption_size = NULL,
  format_headers = TRUE,
  variable_padding = FALSE,
  cell_padding = "normal",
  bold_significant = TRUE,
  bold_variables = FALSE,
  p_threshold = 0.05,
  align = NULL,
  indent_groups = FALSE,
  condense_table = FALSE,
  condense_quantitative = FALSE,
  zebra_stripes = FALSE,
  stripe_color = "gray!20",
  dark_header = FALSE,
  show_logs = FALSE,
  ...
)
```

## Arguments

- table:

  Data frame, data.table, or matrix to export. Can be output from
  [`desctable()`](https://phmcc.codeberg.page/summata/reference/desctable.md),
  [`survtable()`](https://phmcc.codeberg.page/summata/reference/survtable.md),
  [`fit()`](https://phmcc.codeberg.page/summata/reference/fit.md),
  [`uniscreen()`](https://phmcc.codeberg.page/summata/reference/uniscreen.md),
  [`fullfit()`](https://phmcc.codeberg.page/summata/reference/fullfit.md),
  [`compfit()`](https://phmcc.codeberg.page/summata/reference/compfit.md),
  or any tabular data structure.

- file:

  Character string specifying the output PDF filename. Must have `.pdf`
  extension. Example: `"results.pdf"`, `"Table1.pdf"`.

- orientation:

  Character string specifying page orientation:

  - `"portrait"` - Vertical orientation \[default\]

  - `"landscape"` - Horizontal orientation (recommended for wide tables)

- paper:

  Character string specifying paper size:

  - `"letter"` - US Letter (8.5" x 11") \[default\]

  - `"a4"` - A4 (210 mm x 297 mm)

  - `"auto"` - Auto-size to content (no margins, crops to fit)

- margins:

  Numeric vector of length 4 specifying margins in inches as
  `c(top, right, bottom, left)`. Default is `c(1, 1, 1, 1)`. Ignored
  when `paper = "auto"`.

- fit_to_page:

  Logical. If `TRUE`, scales table to fit within the text width
  (respects margins). Useful for wide tables that would otherwise
  overflow. Default is `TRUE`.

- font_size:

  Numeric. Base font size in points. Default is 8. Smaller values
  accommodate more content; larger values improve readability. Typical
  range: 6-12 points.

- caption:

  Character string. Optional caption displayed below the table. Supports
  LaTeX formatting for multi-line captions, superscripts, italics,
  *etc.* See Details for formatting guidance. Default is `NULL`.

- caption_size:

  Numeric. Caption font size in points. If `NULL` (default), uses the
  base `font_size`. Set to a specific value (*e.g.,* 6, 7, 8) to control
  caption size independently of table font size. Useful for fitting
  captions on constrained page sizes. Typical range: 6-10 points.

- format_headers:

  Logical. If `TRUE`, applies automatic formatting to column headers:
  converts underscores to spaces, italicizes statistical notation
  ("*n*", "*p*"), and improves readability. Default is `TRUE`.

- variable_padding:

  Logical. If `TRUE`, adds vertical spacing between different variables
  in the table, creating visual grouping. Particularly useful for
  regression tables with multiple predictors. Default is `FALSE`.

- cell_padding:

  Character string or numeric specifying vertical padding within table
  cells:

  - `"none"` - No extra padding (most compact)

  - `"normal"` - Standard padding \[default\]

  - `"relaxed"` - Increased padding

  - `"loose"` - Maximum padding

  - Numeric value - Custom multiplier (*e.g.,* `1.5`)

  Adjusts `\arraystretch` in LaTeX.

- bold_significant:

  Logical. If `TRUE`, applies bold formatting to p-values below the
  significance threshold, making significant results stand out visually.
  Default is `TRUE`.

- bold_variables:

  Logical. If `TRUE`, variable names are displayed in bold. Default is
  `FALSE`.

- p_threshold:

  Numeric. Threshold for bold *p*-value formatting. Only used when
  `bold_significant = TRUE`. Default is 0.05.

- align:

  Character string or vector specifying column alignment. Options:

  - `"l"` - Left aligned

  - `"c"` - Center aligned

  - `"r"` - Right aligned

  If `NULL` (default), automatically determines alignment based on
  content (text left, numbers right). Can specify per-column:
  `c("l", "l", "r", "r")`.

- indent_groups:

  Logical. If `TRUE`, indents factor levels/groups under their parent
  variable using horizontal space, creating a hierarchical display.
  Useful for factor variables in regression tables. Default is `FALSE`.

- condense_table:

  Logical. If `TRUE`, condenses the table by:

  - Showing only one row per continuous variable (estimate + CI)

  - Showing only non-reference categories for binary factors

  - Automatically setting `indent_groups = TRUE`

  Significantly reduces table height. Default is `FALSE`.

- condense_quantitative:

  Logical. If `TRUE`, condenses continuous and survival variables into
  single rows while preserving all categorical variable rows (including
  binary). Only applies to descriptive tables from
  [`desctable()`](https://phmcc.codeberg.page/summata/reference/desctable.md).
  Automatically sets `indent_groups = TRUE`. Unlike `condense_table`,
  this does not collapse binary categorical variables. Default is
  `FALSE`.

- zebra_stripes:

  Logical. If `TRUE`, applies alternating gray background shading to
  different variables (not individual rows) for improved visual grouping
  and readability. Default is `FALSE`.

- stripe_color:

  Character string. LaTeX color specification for zebra stripes. Default
  is `"gray!20"` (20% gray). Can use other colors like `"blue!10"`,
  `"red!15"`. Requires `zebra_stripes = TRUE`.

- dark_header:

  Logical. If `TRUE`, creates a dark (black) background with white text
  for the header row. Provides strong visual contrast. Default is
  `FALSE`.

- show_logs:

  Logical. If `TRUE`, retains LaTeX log and auxiliary files after PDF
  compilation for troubleshooting. If `FALSE`, deletes these files.
  Default is `FALSE`.

- ...:

  Additional arguments passed to
  [`xtable`](https://rdrr.io/pkg/xtable/man/xtable.html) for advanced
  LaTeX table customization.

## Value

Invisibly returns `NULL`. Creates a PDF file at the specified location.
If compilation fails, check the `.log` file (if `show_logs = TRUE`) for
error details.

## Details

**LaTeX Requirements:**

This function requires a working LaTeX installation. The function checks
for LaTeX availability and provides installation guidance if missing.

*Recommended LaTeX distributions:*

- **TinyTeX** (lightweight, R-integrated): Install via
  [`tinytex::install_tinytex()`](https://rdrr.io/pkg/tinytex/man/install_tinytex.html)

- **TeX Live** (comprehensive, cross-platform)

- **MiKTeX** (Windows)

- **MacTeX** (macOS)

*Required LaTeX packages* (auto-installed with most distributions):

- `fontenc`, `inputenc` - Character encoding

- `array`, `booktabs`, `longtable` - Table formatting

- `graphicx` - Scaling tables

- `geometry` - Page layout

- `pdflscape`, `lscape` - Landscape orientation

- `helvet` - Sans-serif fonts

- `standalone`, `varwidth` - Auto-sizing (for `paper = "auto"`)

- `float`, `caption` - Floats and captions

- `xcolor`, `colortable` - Colors (for `zebra_stripes` or `dark_header`)

**Caption Formatting:**

Captions support LaTeX commands for rich formatting:

    # Multi-line caption with line breaks
    caption = "Table 1: Multivariable Analysis\\
              OR = odds ratio; CI = confidence interval"

    # With superscripts (using LaTeX syntax)
    caption = "Table 1: Results\\
              Adjusted for age and sex\\
              p-values from Wald tests"

    # With special characters (must escape percent signs)
    caption = "Results for income (in thousands)"

**Auto-Sizing (`paper = "auto"`):**

When `paper = "auto"`, the function attempts to create a minimal PDF
sized exactly to the table content:

1.  Using the `standalone` LaTeX class (cleanest output)

2.  Fallback to `pdfcrop` utility if standalone unavailable

3.  Fallback to minimal margins if neither available

**Table Width Management:**

For wide tables that don't fit on the page:

1.  Use `orientation = "landscape"`

2.  Use `fit_to_page = TRUE` (default) to auto-scale

3.  Reduce `font_size` (*e.g.,* 7 or 6)

4.  Consider `paper = "auto"` for maximum flexibility

**Troubleshooting:**

If PDF compilation fails:

1.  Check that LaTeX is installed: Run `Sys.which("pdflatex")`

2.  Set `show_logs = TRUE` and examine the .log file

3.  Common issues:

    - Missing LaTeX packages: Install via package manager

    - Special characters in text: Escape properly

    - Very wide tables: Use landscape or reduce font size

    - Caption formatting: Check LaTeX syntax

## See also

[`autotable`](https://phmcc.codeberg.page/summata/reference/autotable.md)
for automatic format detection,
[`table2tex`](https://phmcc.codeberg.page/summata/reference/table2tex.md)
for LaTeX source files,
[`table2html`](https://phmcc.codeberg.page/summata/reference/table2html.md)
for HTML output,
[`table2docx`](https://phmcc.codeberg.page/summata/reference/table2docx.md)
for Microsoft Word,
[`table2pptx`](https://phmcc.codeberg.page/summata/reference/table2pptx.md)
for PowerPoint,
[`table2rtf`](https://phmcc.codeberg.page/summata/reference/table2rtf.md)
for Rich Text Format,
[`desctable`](https://phmcc.codeberg.page/summata/reference/desctable.md)
for descriptive tables,
[`fit`](https://phmcc.codeberg.page/summata/reference/fit.md) for
regression tables

Other export functions:
[`autotable()`](https://phmcc.codeberg.page/summata/reference/autotable.md),
[`table2docx()`](https://phmcc.codeberg.page/summata/reference/table2docx.md),
[`table2html()`](https://phmcc.codeberg.page/summata/reference/table2html.md),
[`table2pptx()`](https://phmcc.codeberg.page/summata/reference/table2pptx.md),
[`table2rtf()`](https://phmcc.codeberg.page/summata/reference/table2rtf.md),
[`table2tex()`](https://phmcc.codeberg.page/summata/reference/table2tex.md)

## Examples

``` r
data(clintrial)
data(clintrial_labels)

# Create example table
results <- fit(
    data = clintrial,
    outcome = "os_status",
    predictors = c("age", "sex", "treatment", "stage"),
    labels = clintrial_labels
)

# Test that LaTeX can compile (needed for all PDF examples)
has_latex <- local({
  if (!nzchar(Sys.which("pdflatex"))) return(FALSE)
  test_tex <- file.path(tempdir(), "summata_latex_test.tex")
  writeLines(c("\\documentclass{article}", "\\usepackage{booktabs}",
               "\\begin{document}", "test", "\\end{document}"), test_tex)
  tryCatch(
    system2("pdflatex", c("-interaction=nonstopmode",
            paste0("-output-directory=", tempdir()), test_tex),
            stdout = FALSE, stderr = FALSE),
    error = function(e) 1L) == 0L
})

# Example 1: Basic PDF export
if(has_latex){
  table2pdf(results, file.path(tempdir(), "basic_results.pdf"))
}
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/basic_results.pdf

# \donttest{

if(has_latex){

# Example 2: Landscape orientation for wide tables
table2pdf(results, file.path(tempdir(), "wide_results.pdf"),
         orientation = "landscape")

# Example 3: With caption
table2pdf(results, file.path(tempdir(), "captioned.pdf"),
         caption = "Table 1: Multivariable logistic regression results")

# Example 4: Multi-line caption with formatting
table2pdf(results, file.path(tempdir(), "formatted_caption.pdf"),
         caption = "Table 1: Risk Factors for Mortality\\\\
                   aOR = adjusted odds ratio; CI = confidence interval")

# Example 5: Auto-sized PDF (no fixed page dimensions)
table2pdf(results, file.path(tempdir(), "autosize.pdf"),
         paper = "auto")

# Example 6: A4 paper with custom margins
table2pdf(results, file.path(tempdir(), "a4_custom.pdf"),
         paper = "a4",
         margins = c(0.75, 0.75, 0.75, 0.75))

# Example 7: Larger font for readability
table2pdf(results, file.path(tempdir(), "large_font.pdf"),
         font_size = 11)

# Example 8: Indented hierarchical display
table2pdf(results, file.path(tempdir(), "indented.pdf"),
         indent_groups = TRUE)

# Example 9: Condensed table (reduced height)
table2pdf(results, file.path(tempdir(), "condensed.pdf"),
         condense_table = TRUE)

# Example 10: With zebra stripes
table2pdf(results, file.path(tempdir(), "striped.pdf"),
         zebra_stripes = TRUE,
         stripe_color = "gray!15")

# Example 11: Dark header style
table2pdf(results, file.path(tempdir(), "dark_header.pdf"),
         dark_header = TRUE)

# Example 12: Combination of formatting options
table2pdf(results, file.path(tempdir(), "publication_ready.pdf"),
         orientation = "portrait",
         paper = "letter",
         font_size = 9,
         caption = "Table 2: Multivariable Analysis\\\\
                   Model adjusted for age, sex, and clinical factors",
         indent_groups = TRUE,
         zebra_stripes = TRUE,
         bold_significant = TRUE,
         p_threshold = 0.05)

# Example 13: Adjust cell padding
table2pdf(results, file.path(tempdir(), "relaxed_padding.pdf"),
         cell_padding = "relaxed")  # More spacious

# Example 14: No scaling (natural table width)
table2pdf(results, file.path(tempdir(), "no_scale.pdf"),
         fit_to_page = FALSE,
         font_size = 10)

# Example 15: Hide significance bolding
table2pdf(results, file.path(tempdir(), "no_bold.pdf"),
         bold_significant = FALSE)

# Example 16: Custom column alignment
table2pdf(results, file.path(tempdir(), "custom_align.pdf"),
         align = c("c", "c", "c", "c", "c", "c", "c"))

# Example 17: Descriptive statistics table
desc_table <- desctable(clintrial, by = "treatment",
     variables = c("age", "sex", "bmi", "stage"), labels = clintrial_labels)

table2pdf(desc_table, file.path(tempdir(), "descriptive.pdf"),
         caption = "Table 1: Baseline Characteristics by Treatment Group",
         orientation = "landscape")

# Example 18: Model comparison table
models <- list(
     base = c("age", "sex"),
     full = c("age", "sex", "bmi", "treatment")
)

comparison <- compfit(
     data = clintrial,
     outcome = "os_status",
     model_list = models
)

table2pdf(comparison, file.path(tempdir(), "model_comparison.pdf"),
         caption = "Table 3: Model Comparison Statistics")

# Example 19: Very wide table with aggressive fitting
wide_model <- fit(
     data = clintrial,
     outcome = "os_status",
     predictors = c("age", "sex", "race", "bmi", "smoking", 
                   "hypertension", "diabetes", "treatment", "stage")
)

table2pdf(wide_model, file.path(tempdir(), "very_wide.pdf"),
         orientation = "landscape",
         font_size = 7,
         fit_to_page = TRUE,
         condense_table = TRUE)

# Example 20: With caption size control
table2pdf(results, file.path(tempdir(), "caption_size.pdf"),
         font_size = 8,
         caption_size = 6,
         caption = "Table 4: Results with Compact Caption\\\\
                   Smaller caption fits better on constrained pages")

# Example 21: Troubleshooting - keep logs
table2pdf(results, file.path(tempdir(), "debug.pdf"),
         show_logs = TRUE)
# If it fails, check debug.log for error messages

}
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/wide_results.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/captioned.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/formatted_caption.pdf
#> Using standalone class for auto-sized output
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/autosize.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/a4_custom.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/large_font.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/indented.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/condensed.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/striped.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/dark_header.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/publication_ready.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/relaxed_padding.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/no_scale.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/no_bold.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/custom_align.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/descriptive.pdf
#> Auto-detected binary outcome, using logistic regression
#> Fitting base with 2 predictors...
#> Fitting full with 4 predictors...
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/model_comparison.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/very_wide.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/caption_size.pdf
#> Compiling PDF...
#> Table exported to /tmp/RtmpsIjHNN/debug.pdf
# }
```
