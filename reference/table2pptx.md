# Export Table to Microsoft PowerPoint Format (PPTX)

Converts a data frame, data.table, or matrix to a Microsoft PowerPoint
slide (`.pptx`) with a formatted table using the flextable and officer
packages. Creates presentation-ready slides with extensive control over
table formatting, positioning, and layout. Tables can be further edited
in PowerPoint after creation. Ideal for creating data-driven
presentations and conference talks.

## Usage

``` r
table2pptx(
  table,
  file,
  caption = NULL,
  font_size = 10,
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
  width = NULL,
  align = NULL,
  template = NULL,
  layout = "Title and Content",
  master = "Office Theme",
  left = 0.5,
  top = 1.5,
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

  Character string specifying the output PPTX filename. Must have
  `.pptx` extension. Example: `"results.pptx"`, `"Slide1.pptx"`.

- caption:

  Character string. Optional title displayed in the slide's title
  placeholder or as text box above the table. Default is `NULL`.

- font_size:

  Numeric. Base font size in points for table content. Default is 10.
  Typical range for presentations: 10-14 points. Larger than print
  documents for visibility at distance.

- font_family:

  Character string. Font family name for the table. Must be installed on
  the system. Default is `"Arial"`. Common presentation fonts:
  `"Calibri"`, `"Helvetica"`, `"Arial"`.

- format_headers:

  Logical. If `TRUE`, formats column headers by italicizing statistical
  notation ("*n*", "*p*") and improving readability. Default is `TRUE`.

- bold_significant:

  Logical. If `TRUE`, applies bold formatting to *p*-values below the
  significance threshold. Makes important results stand out in
  presentations. Default is `TRUE`.

- bold_variables:

  Logical. If `TRUE`, variable names are displayed in bold. Default is
  `FALSE`.

- p_threshold:

  Numeric. Threshold for bold *p*-value formatting. Only used when
  `bold_significant = TRUE`. Default is 0.05.

- indent_groups:

  Logical. If `TRUE`, indents factor levels under their parent variable,
  creating hierarchical display. Useful for categorical variables.
  Default is `FALSE`.

- condense_table:

  Logical. If `TRUE`, condenses table to essential rows only.
  Automatically sets `indent_groups = TRUE`. Crucial for fitting content
  on slides. Default is `FALSE`.

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
  variables for visual grouping. Improves readability during
  presentations. Default is `FALSE`.

- dark_header:

  Logical. If `TRUE`, creates dark background with light text for header
  row. Provides strong contrast visible from distance. Default is
  `FALSE`.

- width:

  Numeric. Table width in inches. If `NULL` (default), auto-fits to
  slide width (approximately 9 inches for standard 10-inch slide with
  margins). Specify for exact control.

- align:

  Character vector specifying column alignment. Options: `"left"`,
  `"center"`, or `"right"`. If `NULL` (default), automatically
  determines based on content.

- template:

  Character string. Path to custom PPTX template file. If `NULL`
  (default), uses officer's default blank template. Use to match
  corporate branding or conference themes.

- layout:

  Character string. Name of slide layout to use from template. Default
  is `"Title and Content"`. Common layouts:

  - `"Title and Content"` - Standard layout \[default\]

  - `"Blank"` - Empty slide for maximum control

  - `"Title Only"` - Title area only

  - `"Two Content"` - Title with two content areas

- master:

  Character string. Name of slide master to use. Default is
  `"Office Theme"`. Varies by template. Check template for available
  masters.

- left:

  Numeric. Horizontal position from left edge of slide in inches.
  Default is 0.5. Standard slide is 10 inches wide.

- top:

  Numeric. Vertical position from top edge of slide in inches. Default
  is 1.5 (leaves room for title). Standard slide is 7.5 inches tall.
  Adjust based on table size and layout.

- return_ft:

  Logical. If `TRUE`, returns the flextable object directly for further
  customization. If `FALSE` (default), returns invisibly with flextable
  object as attribute. See Details for usage. Default is `FALSE`.

- ...:

  Additional arguments passed to
  [`read_pptx`](https://davidgohel.github.io/officer/reference/read_pptx.html).

## Value

Behavior depends on `return_ft`:

- `return_ft = FALSE`:

  Invisibly returns a list with:

  - `file` - Path to created file

  - `caption` - Caption/title text

  - `layout` - Layout name used

  - `master` - Master name used

  - `template` - Template path (if provided)

  - `position` - List with `left` and `top` coordinates

  Flextable accessible via `attr(result, "flextable")`

- `return_ft = TRUE`:

  Directly returns the flextable object

Always creates a `.pptx` file at the specified location.

## Details

**Package Requirements:**

Requires:

- **flextable** - Table creation and formatting

- **officer** - PowerPoint manipulation

Install: `install.packages(c("flextable", "officer"))`

**Slide Dimensions:**

Standard PowerPoint slide:

- Width: 10 inches (25.4 cm)

- Height: 7.5 inches (19.05 cm)

- Aspect ratio: 4:3 (standard) or 16:9 (widescreen)

Safe content area (with margins):

- Width: ~9 inches

- Height: ~6 inches (accounting for title)

**Positioning:**

The `left` and `top` parameters control table placement:

- (0, 0) = Top-left corner of slide

- Default (0.5, 1.5) = Standard position with title room

- Center: `left = (10 - table_width) / 2`

When caption is provided:

- Attempts to use title placeholder (if layout supports)

- Falls back to text box above table

- Automatically adjusts table position downward

**Slide Layouts:**

Different layouts serve different purposes:

*Title and Content (default):*

- Has title and content placeholders

- Caption goes in title area

- Table in content area

- Most common for data slides

*Blank:*

- No predefined areas

- Maximum flexibility

- Use absolute positioning (`left`, `top`)

- Good for custom layouts

*Title-Only:*

- Title area only

- Large space for table

- Good for data-heavy slides

**Custom Templates:**

Use organizational or conference templates:

    table2pptx(table, "branded.pptx",
             template = "company_template.pptx",
             layout = "Content Layout",  # Name from template
             master = "Company Theme")   # Name from template

To find layout and master names in template:

    pres <- officer::read_pptx("template.pptx")
    officer::layout_summary(pres)

**Multiple Slides:**

Creating presentations with multiple tables:

    # Each call creates new presentation - combine after
    table2pptx(table1, "slide1.pptx", caption = "Results Part 1")
    table2pptx(table2, "slide2.pptx", caption = "Results Part 2")

    # Then manually combine in PowerPoint, or:
    # Use officer to create multi-slide presentation
    pres <- officer::read_pptx()

    # Add first table
    ft1 <- table2pptx(table1, "temp1.pptx", return_ft = TRUE)
    pres <- officer::add_slide(pres)
    pres <- officer::ph_with(pres, ft1,
                             location = officer::ph_location(left = 0.5, top = 1.5))

    # Add second table
    ft2 <- table2pptx(table2, "temp2.pptx", return_ft = TRUE)
    pres <- officer::add_slide(pres)
    pres <- officer::ph_with(pres, ft2,
                             location = officer::ph_location(left = 0.5, top = 1.5))

    print(pres, target = "combined.pptx")

**Further Customization:**

Access the `flextable` object for advanced formatting:

    ft <- table2pptx(table, "base.pptx", return_ft = TRUE)

    # Customize
    ft <- flextable::color(ft, j = "p-value", color = "red")
    ft <- flextable::bg(ft, i = 1, bg = "yellow")
    ft <- flextable::bold(ft, i = ~ estimate > 0, j = "estimate")

    # Save to new slide
    pres <- officer::read_pptx()
    pres <- officer::add_slide(pres)
    pres <- officer::ph_with(pres, ft,
                             location = officer::ph_location(left = 0.5, top = 1.5))
    print(pres, target = "custom.pptx")

## See also

[`autotable`](https://phmcc.github.io/summata/reference/autotable.md)
for automatic format detection,
[`table2docx`](https://phmcc.github.io/summata/reference/table2docx.md)
for Word documents,
[`table2pdf`](https://phmcc.github.io/summata/reference/table2pdf.md)
for PDF output,
[`table2html`](https://phmcc.github.io/summata/reference/table2html.md)
for HTML tables,
[`table2rtf`](https://phmcc.github.io/summata/reference/table2rtf.md)
for Rich Text Format,
[`table2tex`](https://phmcc.github.io/summata/reference/table2tex.md)
for LaTeX output,
[`flextable`](https://davidgohel.github.io/flextable/reference/flextable.html)
for table customization,
[`read_pptx`](https://davidgohel.github.io/officer/reference/read_pptx.html)
for PowerPoint manipulation

Other export functions:
[`autotable()`](https://phmcc.github.io/summata/reference/autotable.md),
[`table2docx()`](https://phmcc.github.io/summata/reference/table2docx.md),
[`table2html()`](https://phmcc.github.io/summata/reference/table2html.md),
[`table2pdf()`](https://phmcc.github.io/summata/reference/table2pdf.md),
[`table2rtf()`](https://phmcc.github.io/summata/reference/table2rtf.md),
[`table2tex()`](https://phmcc.github.io/summata/reference/table2tex.md)

## Examples

``` r
# Create example data
data(clintrial)
data(clintrial_labels)
tbl <- desctable(clintrial, by = "treatment",
    variables = c("age", "sex"), labels = clintrial_labels)

# Basic PowerPoint export
if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  table2pptx(tbl, file.path(tempdir(), "example.pptx"))
}
#> Table exported to /tmp/Rtmpj56AYI/example.pptx

# \donttest{
old_width <- options(width = 180)
# Load data
data(clintrial)
data(clintrial_labels)

# Create regression table
results <- fit(
   data = clintrial,
   outcome = "os_status",
   predictors = c("age", "sex", "treatment"),
   labels = clintrial_labels
)

# Example 1: Basic PowerPoint slide
table2pptx(results, file.path(tempdir(), "results.pptx"))
#> Table exported to /tmp/Rtmpj56AYI/results.pptx

# Example 2: With title
table2pptx(results, file.path(tempdir(), "titled.pptx"),
        caption = "Multivariable Regression Results")
#> Table exported to /tmp/Rtmpj56AYI/titled.pptx

# Example 3: Larger font for visibility
table2pptx(results, file.path(tempdir(), "large_font.pptx"),
        font_size = 12,
        caption = "Main Findings")
#> Table exported to /tmp/Rtmpj56AYI/large_font.pptx

# Example 4: Condensed for slide space
table2pptx(results, file.path(tempdir(), "condensed.pptx"),
        condense_table = TRUE,
        caption = "Key Results")
#> Table exported to /tmp/Rtmpj56AYI/condensed.pptx

# Example 5: Dark header for emphasis
table2pptx(results, file.path(tempdir(), "dark.pptx"),
        dark_header = TRUE,
        caption = "Risk Factors")
#> Table exported to /tmp/Rtmpj56AYI/dark.pptx

# Example 6: With zebra stripes
table2pptx(results, file.path(tempdir(), "striped.pptx"),
        zebra_stripes = TRUE)
#> Table exported to /tmp/Rtmpj56AYI/striped.pptx

# Example 7: Blank layout with custom positioning
table2pptx(results, file.path(tempdir(), "blank.pptx"),
        layout = "Blank",
        left = 1,
        top = 1.5,
        width = 8)
#> Table exported to /tmp/Rtmpj56AYI/blank.pptx

# Example 8: Get flextable for customization
ft <- table2pptx(results, file.path(tempdir(), "base.pptx"), return_ft = TRUE)
#> Table exported to /tmp/Rtmpj56AYI/base.pptx

# Customize the returned flextable object
ft <- flextable::color(ft, j = "p-value", color = "darkred")

# Example 9: Presentation-optimized table
table2pptx(results, file.path(tempdir(), "presentation.pptx"),
        caption = "Main Analysis Results",
        font_size = 11,
        condense_table = TRUE,
        zebra_stripes = TRUE,
        dark_header = TRUE,
        bold_significant = TRUE)
#> Table exported to /tmp/Rtmpj56AYI/presentation.pptx

# Example 10: Descriptive statistics slide
desc <- desctable(
   data = clintrial,
   by = "treatment",
   variables = c("age", "sex", "bmi"),
   labels = clintrial_labels
)

table2pptx(desc, file.path(tempdir(), "baseline.pptx"),
        caption = "Baseline Characteristics",
        font_size = 10)
#> Table exported to /tmp/Rtmpj56AYI/baseline.pptx

# Example 11: Conference presentation style
table2pptx(results, file.path(tempdir(), "conference.pptx"),
        caption = "Study Outcomes",
        font_family = "Calibri",
        font_size = 14,  # Large for big rooms
        dark_header = TRUE,
        condense_table = TRUE)
#> Table exported to /tmp/Rtmpj56AYI/conference.pptx

options(old_width)
# }
```
