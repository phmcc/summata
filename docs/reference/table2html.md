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
data(clintrial)
data(clintrial_labels)

# Create example table
results <- fit(
     data = clintrial,
     outcome = "os_status",
     predictors = c("age", "sex", "treatment", "stage"),
     labels = clintrial_labels
)

# Example 1: Basic HTML export (standalone)
if (requireNamespace("xtable", quietly = TRUE)) {
  table2html(results, file.path(tempdir(), "results.html"))
}
#> Table exported to /tmp/RtmphwDdUy/results.html

# \donttest{
# Example 2: With caption
table2html(results, file.path(tempdir(), "captioned.html"),
          caption = "Table 1: Multivariable Logistic Regression Results")
#> Table exported to /tmp/RtmphwDdUy/captioned.html

# Example 3: For embedding (no CSS)
table2html(results, file.path(tempdir(), "embed.html"),
          include_css = FALSE)
#> Table exported to /tmp/RtmphwDdUy/embed.html
# Include in your HTML document

# Example 4: Hierarchical display
table2html(results, file.path(tempdir(), "indented.html"),
          indent_groups = TRUE)
#> Table exported to /tmp/RtmphwDdUy/indented.html

# Example 5: Condensed table
table2html(results, file.path(tempdir(), "condensed.html"),
          condense_table = TRUE)
#> Table exported to /tmp/RtmphwDdUy/condensed.html

# Example 6: With zebra stripes
table2html(results, file.path(tempdir(), "striped.html"),
          zebra_stripes = TRUE,
          stripe_color = "#F0F0F0")
#> Table exported to /tmp/RtmphwDdUy/striped.html

# Example 7: Dark header style
table2html(results, file.path(tempdir(), "dark.html"),
          dark_header = TRUE)
#> Table exported to /tmp/RtmphwDdUy/dark.html

# Example 8: Combination styling
table2html(results, file.path(tempdir(), "styled.html"),
          zebra_stripes = TRUE,
          dark_header = TRUE,
          bold_significant = TRUE)
#> Table exported to /tmp/RtmphwDdUy/styled.html

# Example 9: Custom stripe color
table2html(results, file.path(tempdir(), "blue_stripes.html"),
          zebra_stripes = TRUE,
          stripe_color = "#E3F2FD")  # Light blue
#> Table exported to /tmp/RtmphwDdUy/blue_stripes.html

# Example 10: Disable significance bolding
table2html(results, file.path(tempdir(), "no_bold.html"),
          bold_significant = FALSE)
#> Table exported to /tmp/RtmphwDdUy/no_bold.html

# Example 11: Stricter significance threshold
table2html(results, file.path(tempdir(), "strict.html"),
          bold_significant = TRUE,
          p_threshold = 0.01)
#> Table exported to /tmp/RtmphwDdUy/strict.html

# Example 12: No header formatting
table2html(results, file.path(tempdir(), "raw_headers.html"),
          format_headers = FALSE)
#> Table exported to /tmp/RtmphwDdUy/raw_headers.html

# Example 13: Descriptive statistics table
desc_table <- desctable(clintrial, by = "treatment",
     variables = c("age", "sex", "bmi"), labels = clintrial_labels)

table2html(desc_table, file.path(tempdir(), "baseline.html"),
          caption = "Table 1: Baseline Characteristics by Treatment Group")
#> Table exported to /tmp/RtmphwDdUy/baseline.html

# Example 14: For R Markdown (no CSS, for inline display)
table2html(results, file.path(tempdir(), "rmd_table.html"),
          include_css = FALSE,
          indent_groups = TRUE)
#> Table exported to /tmp/RtmphwDdUy/rmd_table.html

# Then in R Markdown, use a chunk with results='asis' to display inline:
cat(readLines(file.path(tempdir(), "rmd_table.html")), sep = "\n")
#> <!-- html table generated in R 4.5.2 by xtable 1.8-4 package -->
#> <!-- Mon Mar 23 20:52:53 2026 -->
#> <table border=1>
#> <tr> <th> Variable </th> <th> <i>n</i> </th> <th> Events </th> <th> aOR (95% CI) </th> <th> <i>p</i>-value </th>  </tr>
#>   <tr> <td> Age (years) </td> <td> 847 </td> <td> 606 </td> <td> 1.05 (1.04-1.07) </td> <td> <b>< 0.001</b> </td> </tr>
#>   <tr> <td> Sex </td> <td>  </td> <td>  </td> <td>  </td> <td> - </td> </tr>
#>   <tr> <td> &nbsp;&nbsp;&nbsp;&nbsp;Female </td> <td> 449 </td> <td> 297 </td> <td> reference </td> <td>  </td> </tr>
#>   <tr> <td> &nbsp;&nbsp;&nbsp;&nbsp;Male </td> <td> 398 </td> <td> 309 </td> <td> 2.02 (1.45-2.83) </td> <td>  </td> </tr>
#>   <tr> <td> Treatment Group </td> <td>  </td> <td>  </td> <td>  </td> <td> - </td> </tr>
#>   <tr> <td> &nbsp;&nbsp;&nbsp;&nbsp;Control </td> <td> 194 </td> <td> 149 </td> <td> reference </td> <td>  </td> </tr>
#>   <tr> <td> &nbsp;&nbsp;&nbsp;&nbsp;Drug A </td> <td> 292 </td> <td> 184 </td> <td> 0.44 (0.28-0.68) </td> <td>  </td> </tr>
#>   <tr> <td> &nbsp;&nbsp;&nbsp;&nbsp;Drug B </td> <td> 361 </td> <td> 273 </td> <td> 0.74 (0.47-1.16) </td> <td>  </td> </tr>
#>   <tr> <td> Disease Stage </td> <td>  </td> <td>  </td> <td>  </td> <td> - </td> </tr>
#>   <tr> <td> &nbsp;&nbsp;&nbsp;&nbsp;I </td> <td> 211 </td> <td> 127 </td> <td> reference </td> <td>  </td> </tr>
#>   <tr> <td> &nbsp;&nbsp;&nbsp;&nbsp;II </td> <td> 263 </td> <td> 172 </td> <td> 1.32 (0.88-1.97) </td> <td>  </td> </tr>
#>   <tr> <td> &nbsp;&nbsp;&nbsp;&nbsp;III </td> <td> 241 </td> <td> 186 </td> <td> 2.70 (1.74-4.21) </td> <td>  </td> </tr>
#>   <tr> <td> &nbsp;&nbsp;&nbsp;&nbsp;IV </td> <td> 132 </td> <td> 121 </td> <td> 9.08 (4.66-19.28) </td> <td>  </td> </tr>
#>    </table>

# Example 15: Email-friendly version
table2html(results, file.path(tempdir(), "email.html"),
          include_css = TRUE,  # Self-contained
          zebra_stripes = TRUE,
          caption = "Regression Results - See Attached")
#> Table exported to /tmp/RtmphwDdUy/email.html
# Can be directly included in HTML emails

# Example 16: Publication-ready web version
table2html(results, file.path(tempdir(), "publication.html"),
          caption = "Table 2: Multivariable Analysis of Risk Factors",
          indent_groups = TRUE,
          zebra_stripes = FALSE,  # Clean look
          bold_significant = TRUE,
          dark_header = FALSE)
#> Table exported to /tmp/RtmphwDdUy/publication.html

# Example 17: Modern dark theme
table2html(results, file.path(tempdir(), "dark_theme.html"),
          dark_header = TRUE,
          stripe_color = "#2A2A2A",  # Dark gray stripes
          zebra_stripes = TRUE)
#> Table exported to /tmp/RtmphwDdUy/dark_theme.html

# Example 18: Minimal styling for custom CSS
table2html(results, file.path(tempdir(), "minimal.html"),
          include_css = FALSE,
          format_headers = FALSE,
          bold_significant = FALSE)
#> Table exported to /tmp/RtmphwDdUy/minimal.html
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
#> Auto-detected binary outcome, using logistic regression
#> Fitting base with 2 predictors...
#> Fitting full with 4 predictors...

table2html(comparison, file.path(tempdir(), "comparison.html"),
          caption = "Model Comparison Statistics")
#> Table exported to /tmp/RtmphwDdUy/comparison.html

# }
```
