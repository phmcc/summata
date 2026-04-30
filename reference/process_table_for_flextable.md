# Core flextable processing function

Central processing function for creating flextable objects from data
tables. Handles N row extraction, condensing, indentation, zebra
stripes, formatting, and styling. Used by table2docx, table2pptx, and
table2rtf.

## Usage

``` r
process_table_for_flextable(
  table,
  caption = NULL,
  font_size = 10,
  font_family = "Arial",
  format_headers = TRUE,
  bold_significant = TRUE,
  p_threshold = 0.05,
  indent_groups = FALSE,
  condense_table = FALSE,
  condense_quantitative = FALSE,
  zebra_stripes = FALSE,
  dark_header = FALSE,
  bold_variables = TRUE,
  paper = "letter",
  orientation = "portrait",
  width = NULL,
  align = NULL
)
```

## Arguments

- table:

  Data.frame or data.table to process.

- caption:

  Optional character string for table caption.

- font_size:

  Numeric font size in points.

- font_family:

  Character string font family name.

- format_headers:

  Logical whether to format headers.

- bold_significant:

  Logical whether to bold significant *p*-values.

- p_threshold:

  Numeric *p*-value threshold for significance.

- indent_groups:

  Logical whether to indent group levels.

- condense_table:

  Logical whether to condense all variable types.

- condense_quantitative:

  Logical whether to condense only continuous/survival.

- zebra_stripes:

  Logical whether to apply alternating row shading.

- dark_header:

  Logical whether to use dark header style.

- bold_variables:

  Logical whether to bold variable names (non-indented rows).

- paper:

  Character string paper size.

- orientation:

  Character string page orientation.

- width:

  Optional numeric table width in inches.

- align:

  Optional alignment specification.

## Value

List with ft (flextable object) and caption components.
