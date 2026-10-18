# Name the arguments that each descriptive test accepts

Returns the names of the arguments that
[`desctable()`](https://phmcc.codefloe.page/summata/reference/desctable.md)
forwards from `...` to a given statistical test. A single call can
summarize continuous, categorical, and time-to-event variables, so an
argument intended for one test (*e.g.,* `var.equal` for the *t*-test)
must not reach a test that does not accept it. The lists are explicit
because the test functions are S3 generics whose formal arguments do not
include the arguments of their methods.

## Usage

``` r
test_arg_names(test_type = NULL)
```

## Arguments

- test_type:

  Character string naming the test: `"t"`, `"wrs"`, `"aov"`, `"kwt"`,
  `"chisq"`, `"fisher"`, or `"survival"`. When `NULL`, the names
  accepted by any of the tests are returned.

## Value

Character vector of argument names.
