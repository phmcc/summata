# C-statistic of predictions for a binary outcome

Computes the area under the empirical ROC curve as the Mann-Whitney
statistic, with tied predictions counted as one half, from a single sort
(see
[`average_ranks()`](https://phmcc.codefloe.page/summata/reference/average_ranks.md)).
The orientation matches that of
[`pROC::roc()`](https://rdrr.io/pkg/pROC/man/roc.html) with
`direction = "auto"`, which this function replaces: higher predictions
indicate events unless the median prediction of the events is below that
of the non-events. The result therefore no longer depends on whether the
suggested package is installed.

## Usage

``` r
c_statistic(y, p)
```

## Arguments

- y:

  Binary outcome coded 0 and 1, or logical.

- p:

  Numeric vector of predictions of the same length.

## Value

Numeric scalar; `NA` when either outcome is absent.
