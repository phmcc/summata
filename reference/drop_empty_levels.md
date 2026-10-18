# Drop factor levels without observations before fitting a Cox model

[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html),
[`survival::clogit()`](https://rdrr.io/pkg/survival/man/clogit.html) and
[`coxme::coxme()`](https://rdrr.io/pkg/coxme/man/coxme.html) keep the
levels of a factor that have no observations among the complete cases of
the model, where [`lm()`](https://rdrr.io/r/stats/lm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html) and lme4 drop them from the
model frame. An empty level is then aliased; when it is the first level,
the other levels are estimated against whichever of them is aliased in
turn (see
[`effective_reference_level()`](https://phmcc.codefloe.page/summata/reference/effective_reference_level.md)),
and `coxme()` fails to fit. The empty levels of the factors that a model
uses are therefore dropped before a Cox model is fitted, so that it
takes the first level with observations as its reference, as the other
models do.

## Usage

``` r
drop_empty_levels(data, terms)
```

## Arguments

- data:

  Data table.

- terms:

  Character vector of the terms of the model, such as the outcome,
  predictors, random effects, strata and cluster variable, from which
  the variables are read.

## Value

`data` when none of its factors among the variables has an empty level,
and otherwise a copy in which they have none.

## Details

A level counts as empty when no complete case of the model variables
takes it, as for the model frame of
[`glm()`](https://rdrr.io/r/stats/glm.html), so that a level that occurs
only in observations excluded for missing values is dropped as well.
Those observations, which the model does not use, are then missing for
the factor; the attributes of the factor, such as a variable label, are
kept.
