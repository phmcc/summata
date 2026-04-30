# Export Table to Microsoft Word Format (DOCX)

Converts a data frame, data.table, or matrix to a fully editable
Microsoft Word document (`.docx`) using the flextable and officer
packages. Creates publication-ready tables with extensive formatting
options including typography, alignment, colors, and page layout. Tables
can be further edited in Microsoft Word after creation.

## Usage

``` r
table2docx(
  table,
  file,
  caption = NULL,
  font_size = 8,
  font_family = "Arial",
  format_headers = TRUE,
  bold_significant = TRUE,
  bold_variables = FALSE,
  p_threshold = 0.05,
  indent_groups = FALSE,
  condense_table = FALSE,
  condense_quantitative = FALSE,
  zebra_stripes = FALSE,
  dark_header = FALSE,
  paper = "letter",
  orientation = "portrait",
  width = NULL,
  align = NULL,
  return_ft = FALSE,
  ...
)
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
  or any tabular data.

- file:

  Character string specifying the output DOCX filename. Must have
  `.docx` extension. Example: `"results.docx"`, `"Table1.docx"`.

- caption:

  Character string. Optional caption displayed above the table in the
  Word document. Default is `NULL`.

- font_size:

  Numeric. Base font size in points for table content. Default is 8.
  Typical range: 8-12 points. Headers use slightly larger size.

- font_family:

  Character string. Font family name for the table. Must be a font
  installed on the system. Default is `"Arial"`. Common options:
  `"Times New Roman"`, `"Calibri"`, `"Helvetica"`.

- format_headers:

  Logical. If `TRUE`, formats column headers by italicizing statistical
  notation ("*n*", "*p*"), converting underscores to spaces, and
  improving readability. Default is `TRUE`.

- bold_significant:

  Logical. If `TRUE`, applies bold formatting to *p*-values below the
  significance threshold. Makes significant results stand out. Default
  is `TRUE`.

- bold_variables:

  Logical. If `TRUE`, variable names are displayed in bold. Default is
  `FALSE`.

- p_threshold:

  Numeric. Threshold for bold *p*-value formatting. Only used when
  `bold_significant = TRUE`. Default is 0.05.

- indent_groups:

  Logical. If `TRUE`, indents factor levels under their parent variable
  using horizontal spacing, creating hierarchical display. Useful for
  categorical variables in regression tables. Default is `FALSE`.

- condense_table:

  Logical. If `TRUE`, condenses table by showing only essential rows
  (single row for continuous, non-reference for binary). Automatically
  sets `indent_groups = TRUE`. Significantly reduces table height.
  Default is `FALSE`.

- condense_quantitative:

  Logical. If `TRUE`, condenses continuous and survival variables into
  single rows while preserving all categorical variable rows (including
  binary). Only applies to descriptive tables from
  [`desctable()`](https://phmcc.github.io/summata/reference/desctable.md).
  Automatically sets `indent_groups = TRUE`. Unlike `condense_table`,
  this does not collapse binary categorical variables. Default is
  `FALSE`.

- zebra_stripes:

  Logical. If `TRUE`, applies alternating row shading to different
  variables (not individual rows) for visual grouping. Default is
  `FALSE`.

- dark_header:

  Logical. If `TRUE`, creates a dark background with light text for the
  header row, providing strong visual contrast. Default is `FALSE`.

- paper:

  Character string specifying paper size:

  - `"letter"` - US Letter (8.5" × 11") \[default\]

  - `"a4"` - A4 (210 mm × 297 mm)

  - `"legal"` - US Legal (8.5" × 14")

- orientation:

  Character string specifying page orientation:

  - `"portrait"` - Vertical \[default\]

  - `"landscape"` - Horizontal (for wide tables)

- width:

  Numeric. Table width in inches. If `NULL` (default), automatically
  fits to content and page width. Specify to control exactly.

- align:

  Character vector specifying column alignment for each column. Options:
  `"left"`, `"center"`, or `"right"`. If `NULL` (default), automatically
  determines based on content (text left, numbers right). Example:
  `c("left", "left", "center", "right", "right")`.

- return_ft:

  Logical. If `TRUE`, returns the flextable object directly for further
  customization. If `FALSE` (default), returns invisibly with flextable
  object as attribute. See Details for usage. Default is `FALSE`.

- ...:

  Additional arguments passed to
  [`read_docx`](https://davidgohel.github.io/officer/reference/read_docx.html)
  for document initialization.

## Value

Behavior depends on `return_ft`:

- `return_ft = FALSE`:

  Invisibly returns a list with components:

  - `file` - Path to created file

  - `caption` - Caption text (if provided)

  The flextable object is accessible via `attr(result, "flextable")`

- `return_ft = TRUE`:

  Directly returns the flextable object for immediate further
  customization

In both cases, creates a `.docx` file at the specified location.

## Details

**Package Requirements:**

This function requires:

- **flextable** - For creating formatted tables

- **officer** - For Word document manipulation

Install if needed:

    install.packages(c("flextable", "officer"))

**Output Features:**

The generated Word document contains:

- Fully editable table (native Word table, not image)

- Professional typography and spacing

- Proper page setup (size, orientation, margins)

- Caption (if provided) as separate paragraph above table

- All formatting preserved but editable

- Compatible with Word 2007 and later

**Further Customization:**

For programmatic customization beyond the built-in options, access the
`flextable` object:

*Method 1: Via attribute (default)*

    result <- table2docx(table, "output.docx")
    ft <- attr(result, "flextable")

    # Customize flextable
    ft <- flextable::bold(ft, i = 1, j = 1, part = "body")
    ft <- flextable::color(ft, i = 2, j = 3, color = "red")

    # Re-save if needed
    doc <- officer::read_docx()
    doc <- flextable::body_add_flextable(doc, ft)
    print(doc, target = "customized.docx")

*Method 2: Direct return*

    ft <- table2docx(table, "output.docx", return_ft = TRUE)

    # Customize immediately
    ft <- flextable::bg(ft, bg = "yellow", part = "header")
    ft <- flextable::autofit(ft)

    # Save to new document
    doc <- officer::read_docx()
    doc <- flextable::body_add_flextable(doc, ft)
    print(doc, target = "custom.docx")

**Page Layout:**

The function automatically sets up the Word document with:

- Specified paper size and orientation

- Standard margins (1 inch by default)

- Continuous section (no page breaks before table)

- Left-aligned table placement

For landscape orientation:

- Automatically swaps page width and height

- Applies landscape property to section

- Useful for wide tables with many columns

**Table Width Management:**

Width behavior:

- `width = NULL` - Auto-fits to content and page width

- `width = 6` - Exactly 6 inches wide

- Width distributed evenly across columns by default

- Can adjust individual column widths in Word after creation

For very wide tables:

1.  Use `orientation = "landscape"`

2.  Use `paper = "legal"` for extra width

3.  Reduce `font_size`

4.  Use `condense_table = TRUE`

5.  Consider breaking across multiple tables

**Typography:**

The function applies professional typography:

- Column headers: Bold, slightly larger font

- Body text: Regular weight, specified font size

- Numbers: Right-aligned for easy comparison

- Text: Left-aligned for readability

- Consistent spacing: Adequate padding in cells

Font family must be installed on the system where Word opens the
document. Common cross-platform choices:

- Arial - Sans-serif, highly readable

- Times New Roman - Serif, traditional

- Calibri - Microsoft default, modern

- Helvetica - Sans-serif, professional

**Zebra Striping:**

When `zebra_stripes = TRUE`:

- Alternating variables receive light gray background

- All rows of same variable share same shading

- Improves visual grouping

- Particularly useful for tables with many factor variables

- Color can be changed in Word after creation

**Dark Header:**

When `dark_header = TRUE`:

- Header row: Dark gray/black background

- Header text: White for high contrast

- Modern, professional appearance

- Draws attention to column names

**Integration with R Markdown/Quarto:**

For R Markdown/Quarto Word output:

    # Create flextable for inline display
    ft <- table2docx(results, "temp.docx", return_ft = TRUE)

    # Display in R Markdown chunk
    ft  # Renders in Word output

Or use `flextable` directly in chunks:

    flextable::flextable(results)

## See also

[`autotable`](https://phmcc.github.io/summata/reference/autotable.md)
for automatic format detection,
[`table2pptx`](https://phmcc.github.io/summata/reference/table2pptx.md)
for PowerPoint slides,
[`table2pdf`](https://phmcc.github.io/summata/reference/table2pdf.md)
for PDF output,
[`table2html`](https://phmcc.github.io/summata/reference/table2html.md)
for HTML tables,
[`table2rtf`](https://phmcc.github.io/summata/reference/table2rtf.md)
for Rich Text Format,
[`table2tex`](https://phmcc.github.io/summata/reference/table2tex.md)
for LaTeX output,
[`flextable`](https://davidgohel.github.io/flextable/reference/flextable.html)
for the underlying table object,
[`read_docx`](https://davidgohel.github.io/officer/reference/read_docx.html)
for Word document manipulation

Other export functions:
[`autotable()`](https://phmcc.github.io/summata/reference/autotable.md),
[`table2html()`](https://phmcc.github.io/summata/reference/table2html.md),
[`table2pdf()`](https://phmcc.github.io/summata/reference/table2pdf.md),
[`table2pptx()`](https://phmcc.github.io/summata/reference/table2pptx.md),
[`table2rtf()`](https://phmcc.github.io/summata/reference/table2rtf.md),
[`table2tex()`](https://phmcc.github.io/summata/reference/table2tex.md)

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

# Example 1: Basic Word export
if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  table2docx(results, file.path(tempdir(), "results.docx"))
}
#> Table exported to /tmp/RtmpLSeNVB/results.docx

# \donttest{
old_width <- options(width = 180)
# Example 2: With caption
table2docx(results, file.path(tempdir(), "captioned.docx"),
        caption = "Table 1: Multivariable Logistic Regression Results")
#> Table exported to /tmp/RtmpLSeNVB/captioned.docx

# Example 3: Landscape orientation for wide tables
table2docx(results, file.path(tempdir(), "wide.docx"),
        orientation = "landscape")
#> Table exported to /tmp/RtmpLSeNVB/wide.docx

# Example 4: Custom font and size
table2docx(results, file.path(tempdir(), "custom_font.docx"),
        font_family = "Times New Roman",
        font_size = 11)
#> Table exported to /tmp/RtmpLSeNVB/custom_font.docx

# Example 5: Hierarchical display
table2docx(results, file.path(tempdir(), "indented.docx"),
        indent_groups = TRUE)
#> Table exported to /tmp/RtmpLSeNVB/indented.docx

# Example 6: Condensed table
table2docx(results, file.path(tempdir(), "condensed.docx"),
        condense_table = TRUE)
#> Table exported to /tmp/RtmpLSeNVB/condensed.docx

# Example 7: With zebra stripes
table2docx(results, file.path(tempdir(), "striped.docx"),
        zebra_stripes = TRUE)
#> Table exported to /tmp/RtmpLSeNVB/striped.docx

# Example 8: Dark header style
table2docx(results, file.path(tempdir(), "dark.docx"),
        dark_header = TRUE)
#> Table exported to /tmp/RtmpLSeNVB/dark.docx

# Example 9: A4 paper for international journals
table2docx(results, file.path(tempdir(), "a4.docx"),
        paper = "a4")
#> Table exported to /tmp/RtmpLSeNVB/a4.docx

# Example 10: Get flextable for customization
result <- table2docx(results, file.path(tempdir(), "base.docx"))
#> Table exported to /tmp/RtmpLSeNVB/base.docx
ft <- attr(result, "flextable")

# Customize the flextable
ft <- flextable::bold(ft, i = 1, part = "body")
ft <- flextable::color(ft, j = "p-value", color = "blue")

# Example 11: Direct flextable return
ft <- table2docx(results, file.path(tempdir(), "direct.docx"), return_ft = TRUE)
#> Table exported to /tmp/RtmpLSeNVB/direct.docx
ft <- flextable::bg(ft, bg = "yellow", part = "header")

# Example 12: Publication-ready table
table2docx(results, file.path(tempdir(), "publication.docx"),
        caption = "Table 2: Adjusted Odds Ratios for Mortality",
        font_family = "Times New Roman",
        font_size = 10,
        indent_groups = TRUE,
        zebra_stripes = FALSE,
        bold_significant = TRUE)
#> Table exported to /tmp/RtmpLSeNVB/publication.docx

# Example 13: Custom column alignment
table2docx(results, file.path(tempdir(), "aligned.docx"),
        align = c("left", "left", "center", "right", "right"))
#> Table exported to /tmp/RtmpLSeNVB/aligned.docx

# Example 14: Disable significance bolding
table2docx(results, file.path(tempdir(), "no_bold.docx"),
        bold_significant = FALSE)
#> Table exported to /tmp/RtmpLSeNVB/no_bold.docx

# Example 15: Stricter significance threshold
table2docx(results, file.path(tempdir(), "strict.docx"),
        bold_significant = TRUE,
        p_threshold = 0.01)
#> Table exported to /tmp/RtmpLSeNVB/strict.docx

options(old_width)
# }
```
