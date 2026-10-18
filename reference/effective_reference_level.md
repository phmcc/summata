# Level that serves as the reference of a factor whose first level is empty

Cox models keep factor levels without observations, where
[`lm()`](https://rdrr.io/r/stats/lm.html) and
[`glm()`](https://rdrr.io/r/stats/glm.html) drop them. When the first
level, the reference of the treatment contrasts, has no observations,
the indicators of the other levels sum to one, which the baseline hazard
absorbs. One of them is then aliased, its coefficient fixed at zero, and
the others are estimated against it: that level, not the empty first
level, is the reference of the estimates. A level may also be aliased
because other terms of the model determine it, as when a covariate is
the indicator of that level. Its effect is then carried by those terms
rather than fixed at zero, so that it is not the reference. Where
several levels are aliased, the levels that other terms determine (see
[`determined_levels()`](https://phmcc.codefloe.page/summata/reference/determined_levels.md))
are set aside.

## Usage

``` r
effective_reference_level(levels, counts, aliased, determined = NULL)
```

## Arguments

- levels:

  Character vector of the levels of a factor, in order.

- counts:

  Numeric vector of the number of observations at each level.

- aliased:

  Logical vector indicating, for each level, whether its coefficient is
  aliased.

- determined:

  Optional logical vector indicating, for each level, whether other
  terms of the model determine its indicator.

## Value

The level that serves as the reference, or `NA` when the first level has
observations or no single aliased level with observations can serve.
