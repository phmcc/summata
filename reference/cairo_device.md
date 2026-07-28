# Resolve the Cairo device for an output format

Maps the `"cairo"` shorthand accepted by
[`forestsave()`](https://phmcc.codeberg.page/summata/reference/forestsave.md)
onto the Cairo-backed device for the requested format. Cairo's raster
backends are superseded by ragg and are not offered, so the shorthand
applies to vector formats only.

## Usage

``` r
cairo_device(file)
```

## Arguments

- file:

  Character string giving the output path.

## Value

A device function.
