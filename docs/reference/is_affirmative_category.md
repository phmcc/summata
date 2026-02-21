# Check if category name should be suppressed in condensed label

Determines whether a category name should be suppressed when condensing
binary variables. Returns `TRUE` for standard affirmative values
(*e.g.,* "Yes", "1", "Positive"), standard reference values (*e.g.,*
"No", "Absent", "None"), or when the category name essentially matches
the variable label (case-insensitive comparison).

## Usage

``` r
is_affirmative_category(
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

  Optional character string with the variable label. If provided,
  returns `TRUE` when category is a case-insensitive match or substring.

- norm_category:

  Optional pre-normalized category (lowercase, trimmed). If provided,
  skips normalization for performance.

- norm_label:

  Optional pre-normalized label (lowercase, trimmed). If provided, skips
  normalization for performance.

## Value

Logical indicating whether category should be suppressed.
