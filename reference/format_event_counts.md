# Describe the analyzed events for console output

Renders the event counts held in
[`get_analysis_counts()`](https://phmcc.codeberg.page/summata/reference/get_analysis_counts.md)
using the same formatter as the observation counts, so that the two
lines cannot drift apart in rounding, separators or range handling.

## Usage

``` r
format_event_counts(counts, label = "Events analyzed", marks = NULL)
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

## Value

Character string, or `NULL` when the model carries no events.
