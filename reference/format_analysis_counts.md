# Describe the analyzed sample for console output

Renders the counts produced by
[`get_analysis_counts()`](https://phmcc.codeberg.page/summata/reference/get_analysis_counts.md)
as a single line for the [`print()`](https://rdrr.io/r/base/print.html)
methods. The line is produced whenever the counts are available,
including when every supplied observation entered the analysis: a
complete sample is itself worth stating, and reporting it only on loss
would leave the reader to infer the difference between no exclusions and
no disclosure.

## Usage

``` r
format_analysis_counts(counts, label = "Observations analyzed", marks = NULL)
```

## Arguments

- counts:

  List as returned by
  [`get_analysis_counts()`](https://phmcc.codeberg.page/summata/reference/get_analysis_counts.md).

- label:

  Character label preceding the counts.

- marks:

  List of number marks as returned by
  [`resolve_number_marks()`](https://phmcc.codeberg.page/summata/reference/resolve_number_marks.md).
  Counts and percentages are separated as the accompanying table
  separates them, so that a locale setting applies to the whole of the
  output rather than to the table alone. Resolved from the global option
  when not supplied.

## Value

Character string, or `NULL` when the counts are unavailable.
