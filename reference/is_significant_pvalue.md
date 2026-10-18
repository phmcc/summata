# Identify significant formatted *p*-values

Reads *p*-values as they appear in formatted tables, with the decimal
mark of any locale and with or without a threshold prefix (`"0.023"`,
`"0,023"`, `"< 0.001"`, `"< 0,0001"`). A value given as below a
threshold is significant when that threshold does not exceed
`p_threshold`; any other value is significant when it is below
`p_threshold`. Empty cells, placeholders such as `"-"`, and other text
are not significant.

## Usage

``` r
is_significant_pvalue(x, p_threshold = 0.05)
```

## Arguments

- x:

  Character vector of formatted *p*-values.

- p_threshold:

  Numeric threshold for significance.

## Value

Logical vector of the same length as `x`.
