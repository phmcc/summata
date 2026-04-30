# Detect available sans-serif font for plots

Checks for commonly available sans-serif fonts in order of preference
(Helvetica, Arial, Helvetica Neue) and returns the first available one.
Falls back to "sans" if none are found or if systemfonts is unavailable.

## Usage

``` r
detect_plot_font()
```

## Value

Character string with the font family name to use.

## Details

When ragg is being used as the graphics device (detected via options or
knitr settings), font detection works in non-interactive sessions since
ragg handles font rendering independently of the R graphics system.
