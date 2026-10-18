# Find the reference level of a factor whose first level is empty

Applies
[`effective_reference_level()`](https://phmcc.codefloe.page/summata/reference/effective_reference_level.md),
setting aside the levels that other terms determine where several
aliased levels could serve. Some aliasing leaves no single reference:
the levels may each be estimated against a different aliased level. No
reference is then named, and a warning says so, as when the levels that
other terms determine cannot be found because the model matrix cannot be
rebuilt.

## Usage

``` r
factor_reference(model, var, levels, counts, aliased)
```

## Arguments

- model:

  Fitted model, or `NULL`.

- var:

  Character string naming the factor.

- levels, counts, aliased:

  As for
  [`effective_reference_level()`](https://phmcc.codefloe.page/summata/reference/effective_reference_level.md).

## Value

List with `reference`, the reference level or `NA`, and `ambiguous`,
`TRUE` when the first level is empty but no single level serves as the
reference.
