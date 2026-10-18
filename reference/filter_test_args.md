# Select the test arguments that a given test accepts

Subsets the arguments collected from `...` in
[`desctable()`](https://phmcc.codefloe.page/summata/reference/desctable.md)
to those accepted by one statistical test, as listed by
[`test_arg_names()`](https://phmcc.codefloe.page/summata/reference/test_arg_names.md).

## Usage

``` r
filter_test_args(test_args, test_type)
```

## Arguments

- test_args:

  Named list of arguments collected from `...`.

- test_type:

  Character string naming the test, as in
  [`test_arg_names()`](https://phmcc.codefloe.page/summata/reference/test_arg_names.md).

## Value

Named list of the arguments accepted by the test (possibly empty).
