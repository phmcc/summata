# Check if a binary variable should be condensed without category suffix

Uses a greedy/liberal approach to determine if a binary variable's
condensed display should omit the category name. Returns `TRUE` if
EITHER level of the binary variable is a standard reference/affirmative
value, OR if either level matches/contains the variable label.

## Usage

``` r
should_condense_binary(ref_category, non_ref_category, label = NULL)
```

## Arguments

- ref_category:

  Character string with the reference category name (the level with NA
  estimate).

- non_ref_category:

  Character string with the non-reference category name (the level with
  the actual estimate).

- label:

  Optional character string with the variable label. Used for
  intelligent matching (*e.g.,* "30-Day Readmission" label with "30-day
  readmission" / "No 30-day readmission" levels).

## Value

Logical indicating whether the binary variable should be condensed
without appending the category name.

## Details

This function is designed for binary (2-level) categorical variables
where one level is a reference and one is the "event" or "condition"
level.
