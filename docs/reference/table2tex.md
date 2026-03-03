# Export Table to LaTeX Format

Converts a data frame, data.table, or matrix to LaTeX source code
suitable for inclusion in LaTeX documents. Generates publication-quality
table markup with extensive formatting options including booktabs
styling, color schemes, and hierarchical displays. Output can be
directly `\input{}` or `\include{}` into LaTeX manuscripts. Requires
xtable for export.

## Usage

``` r
table2tex(
  table,
  file,
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
  booktabs = FALSE,
  zebra_stripes = FALSE,
  stripe_color = "gray!20",
  dark_header = FALSE,
  caption = NULL,
  caption_size = NULL,
  label = NULL,
  show_logs = FALSE,
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

  Character string specifying the output `.tex` filename. Must have
  `.tex` extension. Example: `"results.tex"`, `"table1.tex"`.

- format_headers:

  Logical. If `TRUE`, formats column headers by converting underscores
  to spaces, italicizing statistical notation ("*n*", "*p*"), and
  applying title case. Default is `TRUE`.

- variable_padding:

  Logical. If `TRUE`, adds vertical spacing around variable groups using
  `\addlinespace` for improved readability. Default is `FALSE`.

- cell_padding:

  Character string or numeric. Vertical padding within cells:

  - `"none"` - No extra padding

  - `"normal"` - Standard padding \[default\]

  - `"relaxed"` - Increased padding

  - `"loose"` - Maximum padding

  - Numeric - Custom `\arraystretch` value

- bold_significant:

  Logical. If `TRUE`, wraps significant *p*-values in `textbf` commands
  for bold display. Default is `TRUE`.

- bold_variables:

  Logical. If `TRUE`, variable names are displayed in bold. Default is
  `FALSE`.

- p_threshold:

  Numeric. Threshold for bold *p*-value formatting. Only used when
  `bold_significant = TRUE`. Default is 0.05.

- align:

  Character string or vector specifying column alignment:

  - `"l"` - Left

  - `"c"` - Center

  - `"r"` - Right

  - Paragraph column with specified width (p-type)

  If `NULL`, automatically determines based on content. Can specify
  per-column as vector. Default is `NULL`.

- indent_groups:

  Logical. If `TRUE`, uses hspace to indent grouped rows, creating
  hierarchical display. Useful for factor variables in regression
  tables. Default is `FALSE`.

- condense_table:

  Logical. If `TRUE`, condenses table by showing only essential rows
  (single row for continuous, non-reference for binary). Automatically
  sets `indent_groups = TRUE`. Default is `FALSE`.

- condense_quantitative:

  Logical. If `TRUE`, condenses continuous and survival variables into
  single rows while preserving all categorical variable rows (including
  binary). Only applies to descriptive tables from
  [`desctable()`](https://phmcc.github.io/summata/reference/desctable.md).
  Automatically sets `indent_groups = TRUE`. Unlike `condense_table`,
  this does not collapse binary categorical variables. Default is
  `FALSE`.

- booktabs:

  Logical. If `TRUE`, uses booktabs package commands (toprule, midrule,
  bottomrule) for professional table rules. Requires booktabs package in
  LaTeX preamble. Default is `FALSE`.

- zebra_stripes:

  Logical. If `TRUE`, adds alternating row colors for variable groups
  using rowcolor command. Requires `xcolor` package with table option in
  preamble. Default is `FALSE`.

- stripe_color:

  Character string. LaTeX color specification for zebra stripes (*e.g.,*
  `"gray!20"`, `"blue!10"`). Only used when `zebra_stripes = TRUE`.
  Default is `"gray!20"`.

- dark_header:

  Logical. If `TRUE`, creates white text on black background for header
  row. Requires `xcolor` package with `table` option. Default is
  `FALSE`.

- caption:

  Character string. Table caption for LaTeX caption command. Supports
  multi-line captions using double backslash. Default is `NULL`.

- caption_size:

  Numeric. Caption font size in points. If `NULL` (default), caption
  will use the document's default caption size (typically slightly
  smaller than body text). Set to a specific value (*e.g.,* 6, 7, 8, 9)
  to control caption size explicitly. This generates a LaTeX comment
  that you can use when wrapping the table. Typical range: 6-10 points.

- label:

  Character string. LaTeX label for cross-references. Example:
  `"tab:regression"`. Default is `NULL`.

- show_logs:

  Logical. If `TRUE`, displays informational messages about required
  LaTeX packages and formatting options applied. If `FALSE`, suppresses
  these messages. Default is `FALSE`.

- ...:

  Additional arguments passed to
  [`xtable`](https://rdrr.io/pkg/xtable/man/xtable.html).

## Value

Invisibly returns `NULL`. Creates a `.tex` file at the specified
location containing a LaTeX tabular environment.

## Details

**Output Format:**

The function generates a standalone LaTeX tabular environment that can
be:

1.  Included in documents with `\input` command

2.  Embedded in table/figure environments

3.  Used in manuscript classes (`article`, `report`, *etc.*)

The output includes:

- Complete tabular environment with proper alignment

- Horizontal rules (`\hline` or `booktabs` rules)

- Column headers with optional formatting

- Data rows with automatic escaping of special characters

- Optional caption and label commands

**Required LaTeX Packages:**

Add these to your LaTeX document preamble:

*Always required:*

    \usepackage[T1]{fontenc}
    \usepackage[utf8]{inputenc}
    \usepackage{array}
    \usepackage{graphicx}  

*Optional (based on parameters):*

    \usepackage{booktabs}  
    \usepackage[table]{xcolor}  

**Booktabs Style:**

When `booktabs = TRUE`, the table uses publication-quality rules:

- `\toprule` - Heavy rule at top

- `\midrule` - Medium rule below headers

- `\bottomrule` - Heavy rule at bottom

- No vertical rules (`booktabs` style)

- Better spacing around rules

This is the preferred style for most academic journals.

**Color Features:**

*Zebra Stripes:* Creates alternating background colors for visual
grouping:

    zebra_stripes = TRUE
    stripe_color = "gray!20"  # 20% gray
    stripe_color = "blue!10"  # 10% blue

*Dark Header:* Creates high-contrast header row:

    dark_header = TRUE  # Black background, white text

Both require the xcolor package with table option in your document.

**Integration with LaTeX Documents:**

*Basic inclusion:*

    \begin{table}[htbp]
      \centering
      \caption{Regression Results}
      \label{tab:regression}
      \input{results.tex}
    \end{table}

*With resizing:*

    \begin{table}[htbp]
      \centering
      \caption{Results}
      \resizebox{\textwidth}{!}{\input{results.tex}}
    \end{table}

*Landscape orientation:*

    \usepackage{pdflscape}
    \begin{landscape}
      \begin{table}[htbp]
        \centering
        \input{wide_results.tex}
      \end{table}
    \end{landscape}

**Caption Formatting:**

Captions in the `caption` parameter are written as LaTeX comments in the
output file for reference. For actual LaTeX captions, wrap the table in
a table environment (see examples above).

**Special Characters:**

The function automatically escapes LaTeX special characters in your
data:

- Ampersand, percent, dollar sign, hash, underscore

- Left and right braces

- Tilde and caret (using `textasciitilde` and `textasciicircum`)

Variable names and labels should not include these characters unless
intentionally using LaTeX commands.

## See also

[`autotable`](https://phmcc.github.io/summata/reference/autotable.md)
for automatic format detection,
[`table2pdf`](https://phmcc.github.io/summata/reference/table2pdf.md)
for direct PDF output,
[`table2html`](https://phmcc.github.io/summata/reference/table2html.md)
for HTML tables,
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
[`table2html()`](https://phmcc.github.io/summata/reference/table2html.md),
[`table2pdf()`](https://phmcc.github.io/summata/reference/table2pdf.md),
[`table2pptx()`](https://phmcc.github.io/summata/reference/table2pptx.md),
[`table2rtf()`](https://phmcc.github.io/summata/reference/table2rtf.md)

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

# Example 1: Basic LaTeX export
if (requireNamespace("xtable", quietly = TRUE)) {
  table2tex(results, file.path(tempdir(), "basic.tex"))
}
#> Table exported to /tmp/Rtmpox9B2N/basic.tex

# \donttest{
# Example 2: With booktabs for publication
table2tex(results, file.path(tempdir(), "publication.tex"),
       booktabs = TRUE,
       caption = "Multivariable logistic regression results",
       label = "tab:regression")
#> Table exported to /tmp/Rtmpox9B2N/publication.tex

# Example 3: Multi-line caption with abbreviations
table2tex(results, file.path(tempdir(), "detailed.tex"),
       booktabs = TRUE,
       caption = "Table 1: Risk Factors for Mortality\\\\
                 aOR = adjusted odds ratio; CI = confidence interval\\\\
                 Model adjusted for age, sex, treatment, and disease stage",
       label = "tab:mortality")
#> Table exported to /tmp/Rtmpox9B2N/detailed.tex

# Example 4: Hierarchical display with indentation
table2tex(results, file.path(tempdir(), "indented.tex"),
       indent_groups = TRUE,
       booktabs = TRUE)
#> Table exported to /tmp/Rtmpox9B2N/indented.tex

# Example 5: Condensed table (reduced height)
table2tex(results, file.path(tempdir(), "condensed.tex"),
       condense_table = TRUE,
       booktabs = TRUE)
#> Table exported to /tmp/Rtmpox9B2N/condensed.tex

# Example 6: With zebra stripes
table2tex(results, file.path(tempdir(), "striped.tex"),
       zebra_stripes = TRUE,
       stripe_color = "gray!15",
       booktabs = TRUE)
#> Table exported to /tmp/Rtmpox9B2N/striped.tex
# Remember to add \usepackage[table]{xcolor} to the LaTeX document

# Example 7: Dark header style
table2tex(results, file.path(tempdir(), "dark_header.tex"),
       dark_header = TRUE,
       booktabs = TRUE)
#> Table exported to /tmp/Rtmpox9B2N/dark_header.tex
# Requires \usepackage[table]{xcolor}

# Example 8: Custom cell padding
table2tex(results, file.path(tempdir(), "relaxed.tex"),
       cell_padding = "relaxed",
       booktabs = TRUE)
#> Table exported to /tmp/Rtmpox9B2N/relaxed.tex

# Example 9: Custom column alignment (auto-detected by default)
table2tex(results, file.path(tempdir(), "custom_align.tex"),
       align = c("c", "c", "c", "c", "c", "c", "c"))
#> Table exported to /tmp/Rtmpox9B2N/custom_align.tex

# Example 10: No header formatting (keep original names)
table2tex(results, file.path(tempdir(), "raw_headers.tex"),
       format_headers = FALSE)
#> Table exported to /tmp/Rtmpox9B2N/raw_headers.tex

# Example 11: Disable significance bolding
table2tex(results, file.path(tempdir(), "no_bold.tex"),
       bold_significant = FALSE,
       booktabs = TRUE)
#> Table exported to /tmp/Rtmpox9B2N/no_bold.tex

# Example 12: Stricter significance threshold
table2tex(results, file.path(tempdir(), "strict_sig.tex"),
       bold_significant = TRUE,
       p_threshold = 0.01,  # Bold only if p < 0.01
       booktabs = TRUE)
#> Table exported to /tmp/Rtmpox9B2N/strict_sig.tex

# Example 13: With caption size control
table2tex(results, file.path(tempdir(), "caption_size.tex"),
       caption_size = 6,
       caption = "Table 1 - Results with Compact Caption\\\\
                 Smaller caption fits better on constrained pages")
#> Table exported to /tmp/Rtmpox9B2N/caption_size.tex

# Example 14: Complete publication-ready table
table2tex(results, file.path(tempdir(), "final_table1.tex"),
       booktabs = TRUE,
       caption = "Table 1: Multivariable Analysis of Mortality Risk Factors",
       label = "tab:main_results",
       indent_groups = TRUE,
       zebra_stripes = FALSE,  # Many journals prefer no stripes
       bold_significant = TRUE,
       cell_padding = "normal")
#> Table exported to /tmp/Rtmpox9B2N/final_table1.tex

# Example 15: Descriptive statistics table
desc_table <- desctable(clintrial, by = "treatment",
   variables = c("age", "sex", "bmi"), labels = clintrial_labels)

table2tex(desc_table, file.path(tempdir(), "table1_descriptive.tex"),
       booktabs = TRUE,
       caption = "Table 1: Baseline Characteristics",
       label = "tab:baseline")
#> Table exported to /tmp/Rtmpox9B2N/table1_descriptive.tex

# Example 16: Model comparison table
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

table2tex(comparison, file.path(tempdir(), "model_comparison.tex"),
       booktabs = TRUE,
       caption = "Model Comparison Statistics",
       label = "tab:models")
#> Table exported to /tmp/Rtmpox9B2N/model_comparison.tex

# }
```
