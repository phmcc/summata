# Check if category name is a standard reference/negative value

Determines whether a category name represents a standard reference or
negative value that indicates absence. Used to suppress redundant
category names when condensing binary variables.

## Usage

``` r
is_reference_category(
  category,
  label = NULL,
  norm_category = NULL,
  norm_label = NULL
)
```

## Arguments

- category:

  Character string with the category name.

- label:

  Optional character string with the variable label. If provided, checks
  if category is "No \[label\]" or similar patterns.

- norm_category:

  Optional pre-normalized category (lowercase, trimmed). If provided,
  skips normalization for performance.

- norm_label:

  Optional pre-normalized label (lowercase, trimmed). If provided, skips
  normalization for performance.

## Value

Logical indicating whether category is a reference/negative value.
