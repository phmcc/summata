# Resolve the number of worker processes

Uses all but one of the detected cores when no number is requested.
[`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
returns `NA` when the number of cores cannot be determined, which would
otherwise propagate into `mc.cores`; one core is used then. Under
`R CMD check`, the limit of two cores signaled by
`_R_CHECK_LIMIT_CORES_` applies to any number, requested or detected.

## Usage

``` r
resolve_core_count(n_cores = NULL)
```

## Arguments

- n_cores:

  Requested number of cores, or `NULL`.

## Value

Positive integer.
