# Export Table to Rich Text Format (RTF)

Converts a data frame, data.table, or matrix to a Rich Text Format
(`.rtf`) document using the flextable and officer packages. Creates
widely compatible tables with extensive formatting options. RTF files
can be opened and edited in Microsoft Word, LibreOffice, WordPad, and
many other word processors. Particularly useful for regulatory
submissions, cross-platform compatibility, and when maximum editability
is required.

## Usage

``` r
table2rtf(
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

  Character string specifying the output RTF filename. Must have `.rtf`
  extension. Example: `"results.rtf"`, `"Table1.rtf"`.

- caption:

  Character string. Optional caption displayed above the table in the
  RTF document. Default is `NULL`.

- font_size:

  Numeric. Base font size in points for table content. Default is 8.
  Typical range: 8-12 points. Headers use slightly larger size.

- font_family:

  Character string. Font family name for the table. Must be a font
  installed on the system. Default is `"Arial"`. Common options:
  `"Times New Roman"`, `"Calibri"`, `"Courier New"`.

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
  as attribute. See Details for usage. Default is `FALSE`.

- ...:

  Additional arguments (currently unused, reserved for future
  extensions).

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

In both cases, creates a .rtf file at the specified location.

## Details

**Package Requirements:**

This function requires:

- **flextable** - For creating formatted tables

- **officer** - For RTF document generation

Install if needed:

    install.packages(c("flextable", "officer"))

**RTF Format Advantages:**

RTF (Rich Text Format) is a universal document format with several
advantages:

- **Maximum compatibility** - Opens in virtually all word processors

- **Cross-platform** - Works on Windows, Mac, Linux without conversion

- **Fully editable** - Native text format, not embedded objects

- **Lightweight** - Smaller file sizes than DOCX

- **Regulatory compliance** - Widely accepted for submissions (FDA, EMA)

- **Long-term accessibility** - Simple text-based format

- **Version control friendly** - Text-based, works with diff tools

Applications that can open RTF files:

- Microsoft Word (Windows, Mac)

- LibreOffice Writer

- Apache OpenOffice Writer

- WordPad (Windows built-in)

- TextEdit (Mac built-in)

- Google Docs (with import)

- Pages (Mac)

- Many other word processors

**Output Features:**

The generated RTF document contains:

- Fully editable table (native RTF table, not image)

- Professional typography and spacing

- Proper page setup (size, orientation, margins)

- Caption (if provided) as separate paragraph above table

- All formatting preserved but editable

- Compatible with RTF 1.5 specification

**Further Customization:**

For programmatic customization beyond the built-in options, access the
`flextable` object:

*Method 1: Via attribute (default)*

    result <- table2rtf(table, "output.rtf")
    ft <- attr(result, "flextable")

    # Customize flextable
    ft <- flextable::bold(ft, i = 1, j = 1, part = "body")
    ft <- flextable::color(ft, i = 2, j = 3, color = "red")

    # Re-save if needed
    flextable::save_as_rtf(ft, path = "customized.rtf")

*Method 2: Direct return*

    ft <- table2rtf(table, "output.rtf", return_ft = TRUE)

    # Customize immediately
    ft <- flextable::bg(ft, bg = "yellow", part = "header")
    ft <- flextable::autofit(ft)

    # Save to new file
    flextable::save_as_rtf(ft, path = "custom.rtf")

**Page Layout:**

The function automatically sets up the RTF document with:

- Specified paper size and orientation

- Standard margins (1 inch by default)

- Table positioned at document start

- Left-aligned table placement

For landscape orientation:

- Automatically swaps page dimensions

- Applies landscape property

- Useful for wide tables with many columns

**Table Width Management:**

Width behavior:

- `width = NULL` - Auto-fits to content and page width

- `width = 6` - Exactly 6 inches wide

- Width distributed evenly across columns by default

- Can adjust individual column widths in word processor after creation

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

- Statistical notation: Italicized appropriately

## See also

[`autotable`](https://phmcc.github.io/summata/reference/autotable.md)
for automatic format detection,
[`table2docx`](https://phmcc.github.io/summata/reference/table2docx.md)
for Word documents,
[`table2pptx`](https://phmcc.github.io/summata/reference/table2pptx.md)
for PowerPoint slides,
[`table2pdf`](https://phmcc.github.io/summata/reference/table2pdf.md)
for PDF output,
[`table2html`](https://phmcc.github.io/summata/reference/table2html.md)
for HTML tables,
[`table2tex`](https://phmcc.github.io/summata/reference/table2tex.md)
for LaTeX output,
[`flextable`](https://davidgohel.github.io/flextable/reference/flextable.html)
for the underlying table object,
[`save_as_rtf`](https://davidgohel.github.io/flextable/reference/save_as_rtf.html)
for direct RTF export

Other export functions:
[`autotable()`](https://phmcc.github.io/summata/reference/autotable.md),
[`table2docx()`](https://phmcc.github.io/summata/reference/table2docx.md),
[`table2html()`](https://phmcc.github.io/summata/reference/table2html.md),
[`table2pdf()`](https://phmcc.github.io/summata/reference/table2pdf.md),
[`table2pptx()`](https://phmcc.github.io/summata/reference/table2pptx.md),
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

# Example 1: Basic RTF export
if (requireNamespace("flextable", quietly = TRUE)) {
  table2rtf(results, file.path(tempdir(), "results.rtf"))
}
#> Table exported to /tmp/RtmpLSeNVB/results.rtf

# \donttest{
old_width <- options(width = 180)
# Example 2: With caption
table2rtf(results, file.path(tempdir(), "captioned.rtf"),
       caption = "Table 1: Multivariable Logistic Regression Results")
#> Table exported to /tmp/RtmpLSeNVB/captioned.rtf

# Example 3: Landscape orientation for wide tables
table2rtf(results, file.path(tempdir(), "wide.rtf"),
       orientation = "landscape")
#> Table exported to /tmp/RtmpLSeNVB/wide.rtf

# Example 4: Custom font and size
table2rtf(results, file.path(tempdir(), "custom_font.rtf"),
       font_family = "Times New Roman",
       font_size = 11)
#> Table exported to /tmp/RtmpLSeNVB/custom_font.rtf

# Example 5: Hierarchical display
table2rtf(results, file.path(tempdir(), "indented.rtf"),
       indent_groups = TRUE)
#> Table exported to /tmp/RtmpLSeNVB/indented.rtf

# Example 6: Condensed table
table2rtf(results, file.path(tempdir(), "condensed.rtf"),
       condense_table = TRUE)
#> Table exported to /tmp/RtmpLSeNVB/condensed.rtf

# Example 7: With zebra stripes
table2rtf(results, file.path(tempdir(), "striped.rtf"),
       zebra_stripes = TRUE)
#> Table exported to /tmp/RtmpLSeNVB/striped.rtf

# Example 8: Dark header style
table2rtf(results, file.path(tempdir(), "dark.rtf"),
       dark_header = TRUE)
#> Table exported to /tmp/RtmpLSeNVB/dark.rtf

# Example 9: A4 paper for international submissions
table2rtf(results, file.path(tempdir(), "a4.rtf"),
       paper = "a4")
#> Table exported to /tmp/RtmpLSeNVB/a4.rtf

# Example 10: Get flextable for customization
result <- table2rtf(results, file.path(tempdir(), "base.rtf"))
#> Table exported to /tmp/RtmpLSeNVB/base.rtf
ft <- attr(result, "flextable")

# Customize the flextable
ft <- flextable::bold(ft, i = 1, part = "body")
ft <- flextable::color(ft, j = "p-value", color = "blue")

# Re-save
flextable::save_as_rtf(ft, path = file.path(tempdir(), "customized.rtf"))

# Example 11: Direct flextable return
ft <- table2rtf(results, file.path(tempdir(), "direct.rtf"), return_ft = TRUE)
#> Table exported to /tmp/RtmpLSeNVB/direct.rtf
ft <- flextable::bg(ft, bg = "yellow", part = "header")

# Example 12: Regulatory submission table
table2rtf(results, file.path(tempdir(), "submission.rtf"),
       caption = "Table 2: Adjusted Odds Ratios for Mortality",
       font_family = "Times New Roman",
       font_size = 10,
       indent_groups = TRUE,
       zebra_stripes = FALSE,
       bold_significant = TRUE)
#> Table exported to /tmp/RtmpLSeNVB/submission.rtf

# Example 13: Custom column alignment
table2rtf(results, file.path(tempdir(), "aligned.rtf"),
       align = c("left", "left", "center", "right", "right"))
#> Table exported to /tmp/RtmpLSeNVB/aligned.rtf

# Example 14: Disable significance bolding
table2rtf(results, file.path(tempdir(), "no_bold.rtf"),
       bold_significant = FALSE)
#> Table exported to /tmp/RtmpLSeNVB/no_bold.rtf

# Example 15: Stricter significance threshold
table2rtf(results, file.path(tempdir(), "strict.rtf"),
       bold_significant = TRUE,
       p_threshold = 0.01)
#> Table exported to /tmp/RtmpLSeNVB/strict.rtf

# Example 16: Descriptive statistics for baseline characteristics
desc <- desctable(clintrial, by = "treatment",
   variables = c("age", "sex", "bmi", "stage"), labels = clintrial_labels)

table2rtf(desc, file.path(tempdir(), "baseline.rtf"),
       caption = "Table 1: Baseline Patient Characteristics",
       zebra_stripes = TRUE)
#> Table exported to /tmp/RtmpLSeNVB/baseline.rtf

# Example 17: Clinical trial efficacy table
table2rtf(results, file.path(tempdir(), "efficacy.rtf"),
       caption = "Table 3: Primary Efficacy Analysis - Intent to Treat Population",
       font_family = "Courier New",  # Monospace for alignment
       paper = "letter",
       orientation = "landscape",
       condense_table = TRUE)
#> Table exported to /tmp/RtmpLSeNVB/efficacy.rtf

options(old_width)
# }
```
