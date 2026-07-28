# Select a graphics device for forest plot output

Chooses a device whose text metrics reproduce those used to render the
package documentation, falling back to
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)'s
own selection where the preferred device is unavailable.

## Usage

``` r
select_forest_device(file)
```

## Arguments

- file:

  Character string giving the output path.

## Value

A device function, or `NULL` to defer to
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).
