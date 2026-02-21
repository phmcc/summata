# Get paper size for PDF/LaTeX export

Returns paper dimensions and margin settings for the specified paper
size.

## Usage

``` r
get_paper_settings(paper, margins = NULL)
```

## Arguments

- paper:

  Character string: "letter", "a4", or "auto".

- margins:

  Optional numeric vector of margins (length 1 or 4).

## Value

List with latex_paper, width, height, and margins components.
