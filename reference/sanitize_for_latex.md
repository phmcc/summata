# Sanitize certain symbols for LaTeX

Escapes special LaTeX characters ( preserving existing LaTeX commands.
Uses negative lookbehind to avoid double-escaping already escaped
characters.

## Usage

``` r
sanitize_for_latex(x)
```

## Arguments

- x:

  Character vector to sanitize.

## Value

Character vector with special characters escaped for LaTeX.
