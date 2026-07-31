# Recommended Figure Dimensions

Returns the figure dimensions recommended for a forest plot, computed
from the plot's own content by the function that produced it.

## Usage

``` r
recdims(x, units = NULL)
```

## Arguments

- x:

  A forest plot produced by
  [`lmforest()`](https://phmcc.codeberg.page/summata/reference/lmforest.md),
  [`glmforest()`](https://phmcc.codeberg.page/summata/reference/glmforest.md),
  [`coxforest()`](https://phmcc.codeberg.page/summata/reference/coxforest.md),
  [`uniforest()`](https://phmcc.codeberg.page/summata/reference/uniforest.md),
  [`multiforest()`](https://phmcc.codeberg.page/summata/reference/multiforest.md),
  or
  [`autoforest()`](https://phmcc.codeberg.page/summata/reference/autoforest.md).

- units:

  Character string giving the units the dimensions are returned in:
  `"in"` (inches), `"cm"`, or `"mm"`. Defaults to the units the plot was
  created in, so that the dimensions pass through unconverted unless a
  change is requested.

## Value

A named numeric vector with elements `width` and `height`, expressed in
`units`. The units are recorded on the result as a `"units"` attribute,
so that a value carried between functions remains self-describing.

## Details

Forest plots are sized by their content rather than by the page: a plot
with more variables, longer labels, or more factor levels requires a
larger canvas. Each plotting function computes a suitable size and
attaches it to the object it returns, and `recdims()` reports it.

The values returned are those
[`forestsave()`](https://phmcc.codeberg.page/summata/reference/forestsave.md)
applies, and are reported exactly as computed rather than rounded, so
that a figure written with
[`forestsave()`](https://phmcc.codeberg.page/summata/reference/forestsave.md)
and one written with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
at these dimensions are identical. Explicit use is therefore needed only
where the dimensions themselves are wanted, as when sizing a knitr chunk
or arranging several plots on one page.

Where `units` differs from the units the plot was created in, the
dimensions are converted. A plot produced before the units were recorded
carries none, and is treated as having been created in inches.

## See also

[`forestsave`](https://phmcc.codeberg.page/summata/reference/forestsave.md)
for saving at these dimensions,
[`autoforest`](https://phmcc.codeberg.page/summata/reference/autoforest.md)
for producing the plot

Other visualization functions:
[`autoforest()`](https://phmcc.codeberg.page/summata/reference/autoforest.md),
[`coxforest()`](https://phmcc.codeberg.page/summata/reference/coxforest.md),
[`forestsave()`](https://phmcc.codeberg.page/summata/reference/forestsave.md),
[`glmforest()`](https://phmcc.codeberg.page/summata/reference/glmforest.md),
[`lmforest()`](https://phmcc.codeberg.page/summata/reference/lmforest.md),
[`multiforest()`](https://phmcc.codeberg.page/summata/reference/multiforest.md),
[`uniforest()`](https://phmcc.codeberg.page/summata/reference/uniforest.md)

## Examples

``` r
data(clintrial)
data(clintrial_labels)

model <- stats::glm(readmission_30d ~ age + sex + stage,
                    data = clintrial, family = stats::binomial)

p <- glmforest(model, data = clintrial, labels = clintrial_labels)
#> Recommended plot dimensions: width = 13.9 in, height = 5.0 in

# Example 1: Dimensions in the units the plot was created in
recdims(p)
#>    width   height 
#> 13.86667  5.00000 
#> attr(,"units")
#> [1] "in"

# Example 2: Journals commonly specify figure widths in millimeters
recdims(p, units = "mm")
#>    width   height 
#> 352.2133 127.0000 
#> attr(,"units")
#> [1] "mm"

# Example 3: Sizing an output canvas from the individual elements
dims <- recdims(p)
dims[["width"]]
#> [1] 13.86667
dims[["height"]]
#> [1] 5
attr(dims, "units")
#> [1] "in"
```
