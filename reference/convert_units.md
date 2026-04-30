# Convert between units

Converts measurements between different unit systems commonly used in
graphics (inches, centimeters, millimeters, pixels, points).

## Usage

``` r
convert_units(value, from = "in", to = "in", dpi = 96)
```

## Arguments

- value:

  Numeric value to convert.

- from:

  Character string specifying source unit ("in", "cm", "mm", "px",
  "pt").

- to:

  Character string specifying target unit ("in", "cm", "mm", "px",
  "pt").

- dpi:

  Integer dots per inch for pixel conversions (default 96).

## Value

Numeric value in target units.
