# Embed fonts in a PDF after it is written

Post-processes a PDF so that its fonts are embedded, through
[`grDevices::embedFonts()`](https://rdrr.io/r/grDevices/embedFonts.html).
This is offered as an alternative to drawing with a Cairo device, which
embeds fonts but selects its own italic face.

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

## Details

Ghostscript performs the embedding and is not part of R, so its absence
is reported rather than allowed to fail obscurely. A failure leaves the
original file in place: the plot is already written by the time this
runs, and an unembedded figure is more useful than none.
