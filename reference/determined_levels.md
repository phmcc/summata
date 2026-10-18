# Levels of a factor that the other terms of a model determine

Tests, for each level of a factor, whether its indicator lies in the
span of the columns of the model's other terms and of the constant, or
of the strata of a stratified model, whose baseline hazards absorb
constants within strata. The columns are centered within strata (or
about their means), which takes the place of the constant or the strata,
and the span is read from the model matrix, whose rows are the
observations analyzed. The columns of the interactions that the factor
enters are set aside with its own, as they are built from its levels.

## Usage

``` r
determined_levels(model, var, levels)
```

## Arguments

- model:

  Fitted `coxph` or `clogit` model.

- var:

  Character string naming the factor.

- levels:

  Character vector of its levels.

## Value

Logical vector, one per level, `TRUE` where the other terms determine
the level, or `NULL` when the model matrix cannot be rebuilt for the
observations analyzed.
