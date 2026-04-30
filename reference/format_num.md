# Format a numeric value with locale-aware separators

General-purpose number formatter used by all display functions. For
values \\\ge\\ 1000 (in absolute value), inserts the appropriate
thousands separator. Fixes negative-zero display artefacts.

## Usage

``` r
format_num(x, fmt_str, marks)
```

## Arguments

- x:

  Numeric value to format.

- fmt_str:

  Character string `sprintf` format specification (*e.g.,* `"%.1f"`).
  Used only for \|x\| \< 1000.

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.github.io/summata/reference/resolve_number_marks.md).

## Value

Character string with the formatted number.
