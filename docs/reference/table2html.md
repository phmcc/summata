# Export Table to HTML Format

Converts a data frame, data.table, or matrix to HTML format with
optional CSS styling for web display, HTML documents, or embedding in
web applications. Generates clean, standards-compliant HTML with
professional styling options including responsive design support, color
schemes, and interactive features. Requires xtable for export.

## Usage

``` r
table2html(
  table,
  file,
  caption = NULL,
  format_headers = TRUE,
  variable_padding = FALSE,
  bold_significant = TRUE,
  bold_variables = FALSE,
  p_threshold = 0.05,
  indent_groups = FALSE,
  condense_table = FALSE,
  condense_quantitative = FALSE,
  zebra_stripes = FALSE,
  stripe_color = "#EEEEEE",
  dark_header = FALSE,
  include_css = TRUE,
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

  Character string specifying the output HTML filename. Must have
  `.html` or `.htm` extension. Example: `"results.html"`.

- caption:

  Character string. Optional caption displayed below the table. Supports
  basic HTML formatting. Default is `NULL`.

- format_headers:

  Logical. If `TRUE`, formats column headers by converting underscores
  to spaces and applying proper casing. Wraps statistical notation
  ("*n*", "*p*") in `<i>` tags. Default is `TRUE`.

- variable_padding:

  Logical. If `TRUE`, adds vertical spacing around variable groups for
  improved readability. Particularly useful for multi-row factor
  variables. Default is `FALSE`.

- bold_significant:

  Logical. If `TRUE`, wraps significant *p*-values in `<b>` tags for
  bold display. Makes important results stand out visually. Default is
  `TRUE`.

- bold_variables:

  Logical. If `TRUE`, variable names are displayed in bold. Default is
  `FALSE`.

- p_threshold:

  Numeric. Threshold for bold *p*-value formatting. Only used when
  `bold_significant = TRUE`. Default is 0.05.

- indent_groups:

  Logical. If `TRUE`, indents grouped rows using non-breaking spaces
  (`&nbsp;`) for hierarchical display. Useful for factor variables in
  regression output. Default is `FALSE`.

- condense_table:

  Logical. If `TRUE`, condenses table by showing only essential rows.
  Automatically sets `indent_groups = TRUE`. Default is `FALSE`.

- condense_quantitative:

  Logical. If `TRUE`, condenses continuous and survival variables into
  single rows while preserving all categorical variable rows (including
  binary). Only applies to descriptive tables from
  [`desctable()`](https://phmcc.github.io/summata/reference/desctable.md).
  Automatically sets `indent_groups = TRUE`. Unlike `condense_table`,
  this does not collapse binary categorical variables. Default is
  `FALSE`.

- zebra_stripes:

  Logical. If `TRUE`, applies alternating background shading to
  different variables (not individual rows) for visual grouping. Default
  is `FALSE`.

- stripe_color:

  Character string. HTML color specification for zebra stripes. Can use
  hex codes (`"#EEEEEE"`), RGB (`"rgb(238,238,238)"`), or color names
  (`"lightgray"`). Default is `"#EEEEEE"`.

- dark_header:

  Logical. If `TRUE`, creates black background with white text for the
  header row. Provides strong visual contrast. Default is `FALSE`.

- include_css:

  Logical. If `TRUE`, includes embedded CSS styling in the output file
  for standalone HTML. Set to `FALSE` when embedding in existing HTML
  with its own stylesheet. Default is `TRUE`.

- ...:

  Additional arguments passed to
  [`xtable`](https://rdrr.io/pkg/xtable/man/xtable.html).

## Value

Invisibly returns `NULL`. Creates an HTML file at the specified location
that can be opened in web browsers or embedded in HTML documents.

## Details

**Output Format:**

The function generates standards-compliant HTML5 markup with:

- Semantic `<table>` structure

- Proper `<thead>` and `<tbody>` sections

- Accessible header cells (`<th>`)

- Clean, readable markup

- Optional embedded CSS styling

**Standalone vs. Embedded:**

*Standalone HTML (`include_css = TRUE`):*

- Can be opened directly in web browsers

- Includes all necessary styling

- Self-contained, portable

- Suitable for sharing via email or web hosting

*Embedded HTML (`include_css = FALSE`):*

- For inclusion in existing HTML documents

- No CSS included (use parent document's styles)

- Smaller file size

- Integrates with web frameworks (Shiny, R Markdown, Quarto)

**CSS Styling:**

When `include_css = TRUE`, the function applies professional styling:

- **Table:** Border-collapse, sans-serif font (Arial), 20px margin

- **Cells:** 8px vertical × 12px horizontal padding, left-aligned text

- **Borders:** 1px solid `#DDD` (light gray)

- **Headers:** Bold text, light gray background (`#F2F2F2`)

- **Numeric columns:** Center-aligned (auto-detected)

- **Caption:** Bold, 1.1em font, positioned below table

*With `dark_header = TRUE`:*

- Header background: Black (`#000000`)

- Header text: White (`#FFFFFF`)

- Creates high contrast, modern appearance

*With `zebra_stripes = TRUE`:*

- Alternating variable groups receive background color

- Default color: `#EEEEEE` (light gray)

- Applied via CSS class `.zebra-stripe`

- Groups entire variable (all factor levels together)

**Hierarchical Display:**

The `indent_groups` option creates visual hierarchy using HTML
non-breaking spaces:

    <td><b>Treatment</b></td>  <!-- Variable name -->
    <td>&nbsp;&nbsp;&nbsp;&nbsp;Control</td>  <!-- Indented level -->
    <td>&nbsp;&nbsp;&nbsp;&nbsp;Active</td>   <!-- Indented level -->

**Integration with R Markdown/Quarto:**

For R Markdown or Quarto documents:

    # Generate HTML fragment (no CSS)
    table2html(results, "table.html", include_css = FALSE)

Then include in your document chunk with `results='asis'`:

    cat(readLines("table.html"), sep = "\n")

Or directly render without file:

    # For inline display
    htmltools::HTML(
      capture.output(
        print(xtable::xtable(results), type = "html")
      )
    )

**Integration with Shiny:**

For Shiny applications:

    # In server function
    output$results_table <- renderUI({
      table2html(results_data(), "temp.html", include_css = FALSE)
      HTML(readLines("temp.html"))
    })

    # Or use directly with DT package for interactive tables
    output$interactive_table <- DT::renderDT({
      results_data()
    })

**Accessibility:**

The generated HTML follows accessibility best practices:

- Semantic table structure

- Proper header cells (`<th>`) with scope attributes

- Clear visual hierarchy

- Adequate color contrast (when using default styles)

- Screen reader friendly markup

## See also

[`autotable`](https://phmcc.github.io/summata/reference/autotable.md)
for automatic format detection,
[`table2pdf`](https://phmcc.github.io/summata/reference/table2pdf.md)
for PDF output,
[`table2tex`](https://phmcc.github.io/summata/reference/table2tex.md)
for LaTeX output,
[`table2docx`](https://phmcc.github.io/summata/reference/table2docx.md)
for Word documents,
[`table2pptx`](https://phmcc.github.io/summata/reference/table2pptx.md)
for PowerPoint,
[`table2rtf`](https://phmcc.github.io/summata/reference/table2rtf.md)
for Rich Text Format,
[`fit`](https://phmcc.github.io/summata/reference/fit.md) for regression
tables,
[`desctable`](https://phmcc.github.io/summata/reference/desctable.md)
for descriptive tables

Other export functions:
[`autotable()`](https://phmcc.github.io/summata/reference/autotable.md),
[`table2docx()`](https://phmcc.github.io/summata/reference/table2docx.md),
[`table2pdf()`](https://phmcc.github.io/summata/reference/table2pdf.md),
[`table2pptx()`](https://phmcc.github.io/summata/reference/table2pptx.md),
[`table2rtf()`](https://phmcc.github.io/summata/reference/table2rtf.md),
[`table2tex()`](https://phmcc.github.io/summata/reference/table2tex.md)

## Examples

``` r
if (FALSE) {
# Load data and create table
data(clintrial)
data(clintrial_labels)

results <- fit(
     data = clintrial,
     outcome = "os_status",
     predictors = c("age", "sex", "treatment", "stage"),
     labels = clintrial_labels
)

# Example 1: Basic HTML export (standalone)
table2html(results, "results.html")
# Open results.html in web browser

# Example 2: With caption
table2html(results, "captioned.html",
          caption = "Table 1: Multivariable Logistic Regression Results")

# Example 3: For embedding (no CSS)
table2html(results, "embed.html",
          include_css = FALSE)
# Include in your HTML document

# Example 4: Hierarchical display
table2html(results, "indented.html",
          indent_groups = TRUE)

# Example 5: Condensed table
table2html(results, "condensed.html",
          condense_table = TRUE)

# Example 6: With zebra stripes
table2html(results, "striped.html",
          zebra_stripes = TRUE,
          stripe_color = "#F0F0F0")

# Example 7: Dark header style
table2html(results, "dark.html",
          dark_header = TRUE)

# Example 8: Combination styling
table2html(results, "styled.html",
          zebra_stripes = TRUE,
          dark_header = TRUE,
          bold_significant = TRUE)

# Example 9: Custom stripe color
table2html(results, "blue_stripes.html",
          zebra_stripes = TRUE,
          stripe_color = "#E3F2FD")  # Light blue

# Example 10: Disable significance bolding
table2html(results, "no_bold.html",
          bold_significant = FALSE)

# Example 11: Stricter significance threshold
table2html(results, "strict.html",
          bold_significant = TRUE,
          p_threshold = 0.01)

# Example 12: No header formatting
table2html(results, "raw_headers.html",
          format_headers = FALSE)

# Example 13: Descriptive statistics table
desc_table <- desctable(
     data = clintrial,
     by = "treatment",
     variables = c("age", "sex", "bmi"),
     labels = clintrial_labels
)

table2html(desc_table, "baseline.html",
          caption = "Table 1: Baseline Characteristics by Treatment Group")

# Example 14: For R Markdown (no CSS, for inline display)
table2html(results, "rmd_table.html",
          include_css = FALSE,
          indent_groups = TRUE)

# Then in R Markdown:
# ```{r results='asis', echo=FALSE}
# cat(readLines("rmd_table.html"), sep = "\n")
# ```

# Example 15: Email-friendly version
table2html(results, "email.html",
          include_css = TRUE,  # Self-contained
          zebra_stripes = TRUE,
          caption = "Regression Results - See Attached")
# Can be directly included in HTML emails

# Example 16: Publication-ready web version
table2html(results, "publication.html",
          caption = "Table 2: Multivariable Analysis of Risk Factors",
          indent_groups = TRUE,
          zebra_stripes = FALSE,  # Clean look
          bold_significant = TRUE,
          dark_header = FALSE)

# Example 17: Modern dark theme
table2html(results, "dark_theme.html",
          dark_header = TRUE,
          stripe_color = "#2A2A2A",  # Dark gray stripes
          zebra_stripes = TRUE)

# Example 18: Minimal styling for custom CSS
table2html(results, "minimal.html",
          include_css = FALSE,
          format_headers = FALSE,
          bold_significant = FALSE)
# Apply your own CSS classes and styling

# Example 19: Model comparison table
models <- list(
     base = c("age", "sex"),
     full = c("age", "sex", "treatment", "stage")
)

comparison <- compfit(
     data = clintrial,
     outcome = "os_status",
     model_list = models
)

table2html(comparison, "comparison.html",
          caption = "Model Comparison Statistics")
}
```
