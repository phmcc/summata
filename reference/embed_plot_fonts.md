# Embed fonts in a PDF after it is written

Post-processes a PDF so that its fonts are embedded, through
[`grDevices::embedFonts()`](https://rdrr.io/r/grDevices/embedFonts.html).
This is offered as an alternative to drawing with a Cairo device, which
embeds fonts but renders italics unreliably. If Cairo is not selected,
embedding is performed by Ghostscript (external).

## Usage

``` r
embed_plot_fonts(file, quiet = FALSE)
```

## Arguments

- file:

  Character string giving the path to the written file.

- quiet:

  Logical. Suppress the confirmation message.

## Value

Invisibly returns `TRUE` where the fonts were embedded and `FALSE`
otherwise.
