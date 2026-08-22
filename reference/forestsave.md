# Save a Forest Plot to File

Writes a forest plot to file, sized by default from the recommendations
the plotting function computed from the plot's own content, and rendered
through a graphics device chosen to reproduce the fonts seen on screen.

## Usage

``` r
forestsave(
  plot,
  file,
  width = NULL,
  height = NULL,
  units = NULL,
  dpi = 300,
  device = NULL,
  embed_fonts = FALSE,
  quiet = FALSE,
  ...
)
```

## Arguments

- plot:

  A forest plot produced by
  [`lmforest()`](https://phmcc.codeberg.page/summata/reference/lmforest.md),
  [`glmforest()`](https://phmcc.codeberg.page/summata/reference/glmforest.md),
  [`coxforest()`](https://phmcc.codeberg.page/summata/reference/coxforest.md),
  [`uniforest()`](https://phmcc.codeberg.page/summata/reference/uniforest.md),
  [`multiforest()`](https://phmcc.codeberg.page/summata/reference/multiforest.md),
  or
  [`autoforest()`](https://phmcc.codeberg.page/summata/reference/autoforest.md).
  Any ggplot2 object is accepted, but one that carries no recommended
  dimensions requires `width` and `height` to be supplied.

- file:

  Character string giving the output path. The file extension determines
  the format: `.pdf`, `.png`, `.tiff`, `.jpeg`, `.svg`, `.eps`, and the
  other formats
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  supports.

- width, height:

  Numeric plot dimensions. Taken from the plot's recommended dimensions
  when not supplied, and converted where these differ from `units`.

- units:

  Character string giving the units for `width` and `height`: `"in"`
  (inches), `"cm"`, or `"mm"`. Defaults to the units the recommended
  dimensions were expressed in, and to `"in"` where none were recorded.

- dpi:

  Numeric resolution in dots per inch for raster formats. Ignored for
  vector formats. Default is `300`.

- device:

  Graphics device. Selected automatically when `NULL`: see Details.
  Supply `"cairo"` to use the Cairo device for the output format, which
  embeds fonts, or any device function or name recognized by
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

- embed_fonts:

  Logical. Embed the fonts in a PDF after it is written, using
  [`grDevices::embedFonts()`](https://rdrr.io/r/grDevices/embedFonts.html).
  Requires Ghostscript, and applies only to PDF output. Default is
  `FALSE`.

- quiet:

  Logical. Suppress the message reporting the dimensions used. Default
  is `FALSE`.

- ...:

  Additional arguments passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html),
  such as `bg` or `limitsize`.

## Value

Invisibly returns `file`.

## Details

Forest plots are sized by their content rather than by the page: a plot
with more variables, longer labels, or more factor levels requires a
larger canvas, and the plotting functions compute one and attach it to
the returned object as a `"rec_dims"` attribute. `forestsave()` reads
that attribute, so the following are equivalent:

    dims <- attr(p, "rec_dims")
    ggplot2::ggsave("forest.pdf", p, width = dims$width,
                    height = dims$height, units = dims$units)

    forestsave(p, "forest.pdf")

**Device selection.** Fonts are measured differently by different
graphics devices, so a figure written through one device may not
reproduce the spacing seen on screen or in the package documentation.
Unless `device` is supplied, `forestsave()` selects:

- PDF:

  R's internal PDF device,
  [`grDevices::pdf()`](https://rdrr.io/r/grDevices/pdf.html), reached by
  deferring to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).
  Cairo is not used by default: although it embeds fonts, its face
  selection can substitute an incorrect italic, which matters here
  because forest plot footers and statistical notation are set in
  italic. Pass `device = "cairo"` to use it where font embedding is
  required.

- PNG, TIFF, JPEG:

  The corresponding ragg device where that package is installed, whose
  text metrics match those used to render the package vignettes. Falls
  back to the standard raster devices otherwise.

- Other formats:

  Left to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

**Font embedding.** R's PDF device draws on the standard fonts that
every PDF reader is required to provide, so its output is portable
without embedding. Embedding matters where the plot uses a font outside
that set, and where a publisher requires it. Setting
`embed_fonts = TRUE` embeds the fonts after the file is written, through
[`grDevices::embedFonts()`](https://rdrr.io/r/grDevices/embedFonts.html),
which requires a Ghostscript installation. This is offered as an
alternative to `device = "cairo"`, which embeds fonts as it draws but
selects its own italic face.

`device = "cairo"` is a shorthand of this package's own, resolving to
[`grDevices::cairo_pdf()`](https://rdrr.io/r/grDevices/cairo.html),
[`grDevices::cairo_ps()`](https://rdrr.io/r/grDevices/cairo.html), or
[`grDevices::svg()`](https://rdrr.io/r/grDevices/cairo.html) according
to the output format:

    forestsave(p, "forest.pdf", device = "cairo")

It applies to vector formats only, since raster output already uses the
ragg devices, which supersede Cairo's raster backends. Cairo is absent
from some installations even where `capabilities("cairo")` reports
otherwise, so the device is opened once before use and its
unavailability reported. Any other value of `device` is passed to
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
unchanged.

The dimensions used are reported through a
[`message()`](https://rdrr.io/r/base/message.html) unless
`quiet = TRUE`, so that a figure written at an unexpected size is
apparent at the point it is written.

## See also

[`autoforest`](https://phmcc.codeberg.page/summata/reference/autoforest.md),
[`coxforest`](https://phmcc.codeberg.page/summata/reference/coxforest.md),
[`glmforest`](https://phmcc.codeberg.page/summata/reference/glmforest.md),
[`lmforest`](https://phmcc.codeberg.page/summata/reference/lmforest.md),
[`uniforest`](https://phmcc.codeberg.page/summata/reference/uniforest.md),
[`multiforest`](https://phmcc.codeberg.page/summata/reference/multiforest.md)
for producing forest plots;
[`tablesave`](https://phmcc.codeberg.page/summata/reference/tablesave.md)
for exporting formatted tables

Other visualization functions:
[`autoforest()`](https://phmcc.codeberg.page/summata/reference/autoforest.md),
[`coxforest()`](https://phmcc.codeberg.page/summata/reference/coxforest.md),
[`glmforest()`](https://phmcc.codeberg.page/summata/reference/glmforest.md),
[`lmforest()`](https://phmcc.codeberg.page/summata/reference/lmforest.md),
[`multiforest()`](https://phmcc.codeberg.page/summata/reference/multiforest.md),
[`recdims()`](https://phmcc.codeberg.page/summata/reference/recdims.md),
[`uniforest()`](https://phmcc.codeberg.page/summata/reference/uniforest.md)

## Examples

``` r
data(clintrial)
data(clintrial_labels)

model <- stats::glm(readmission_30d ~ age + sex + stage,
                    data = clintrial, family = stats::binomial)

p <- glmforest(model, data = clintrial, labels = clintrial_labels)
#> Recommended plot dimensions: width = 13.9 in, height = 5.0 in

# Example 1: Save at the recommended dimensions
forestsave(p, file.path(tempdir(), "forest.pdf"))
#> Forest plot saved to /tmp/RtmppNr8bf/forest.pdf (width = 13.9 in, height = 5.0 in)

# \donttest{

# Example 2: Raster output at publication resolution
forestsave(p, file.path(tempdir(), "forest.png"), dpi = 600)
#> Forest plot saved to /tmp/RtmppNr8bf/forest.png (width = 13.9 in, height = 5.0 in)

# Example 3: Dimensions given explicitly, in any supported units
forestsave(p, file.path(tempdir(), "forest_sized.pdf"),
           width = 24, height = 16, units = "cm")
#> Forest plot saved to /tmp/RtmppNr8bf/forest_sized.pdf (width = 24.0 cm, height = 16.0 cm)

# Example 4: Embed fonts after writing, where Ghostscript is available
forestsave(p, file.path(tempdir(), "forest_embedded.pdf"),
           embed_fonts = TRUE)
#> Forest plot saved to /tmp/RtmppNr8bf/forest_embedded.pdf (width = 13.9 in, height = 5.0 in)
#> Fonts embedded in /tmp/RtmppNr8bf/forest_embedded.pdf

# }
```
