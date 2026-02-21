# Installation and Setup

This document provides instructions for installing the `summata` package
and configuring the dependencies required for complete functionality.
The package requires R ≥ 4.2.0 and several additional components for
document generation.

------------------------------------------------------------------------

## Package Installation

### From Repository

This package is not yet available on CRAN. The stable release may be
installed from GitHub:

``` r
install.packages("remotes")
remotes::install_github("phmcc/summata")
```

The development version is available from Codeberg:

``` r
remotes::install_git("https://codeberg.org/phmcc/summata.git")
```

### Verification

Confirm successful installation by loading the package:

``` r
library(summata)
packageVersion("summata")
```

------------------------------------------------------------------------

## Dependencies

### Required Dependencies

The following packages are installed automatically as dependencies:

| Package      | Purpose                                                     |
|:-------------|:------------------------------------------------------------|
| `data.table` | High-performance data manipulation and core data operations |
| `survival`   | Survival and time-to-event analysis                         |
| `ggplot2`    | Visualization system for forest plots and graphics          |
| `stats`      | Statistical functions and model fitting                     |
| `grDevices`  | Graphics devices and color palettes                         |

### Optional Dependencies

Additional packages extend functionality for specific use-cases:

| Package | Purpose | Required For |
|:---|:---|:---|
| `lme4` | Linear mixed-effects models | `lmer()`, `glmer()` model types |
| `coxme` | Mixed-effects Cox models | `coxme()` model type |
| `xtable` | LaTeX table generation | [`table2pdf()`](https://phmcc.github.io/summata/reference/table2pdf.md), [`table2tex()`](https://phmcc.github.io/summata/reference/table2tex.md), [`table2html()`](https://phmcc.github.io/summata/reference/table2html.md) |
| `flextable` | Office document tables | [`table2docx()`](https://phmcc.github.io/summata/reference/table2docx.md), [`table2pptx()`](https://phmcc.github.io/summata/reference/table2pptx.md), [`table2rtf()`](https://phmcc.github.io/summata/reference/table2rtf.md) |
| `officer` | Office document creation | [`table2docx()`](https://phmcc.github.io/summata/reference/table2docx.md), [`table2pptx()`](https://phmcc.github.io/summata/reference/table2pptx.md) |
| `knitr` | Dynamic document generation | R Markdown integration, vignettes |
| `rmarkdown` | R Markdown documents | Vignette rendering |
| `ragg` | High-quality graphics | Enhanced PNG rendering with better fonts |
| `systemfonts` | Font management | Font detection for graphics |
| `tinytex` | R-integrated LaTeX distribution | PDF table export via [`table2pdf()`](https://phmcc.github.io/summata/reference/table2pdf.md) |
| `MASS` | Statistical methods | Model diagnostics and testing |
| `MuMIn` | Multi-model inference | Model selection and averaging |
| `pROC` | ROC curve analysis | Diagnostic performance evaluation |
| `ResourceSelection` | Goodness-of-fit tests | Hosmer-Lemeshow tests for logistic models |
| `stringr` | String manipulation | Text processing utilities |
| `utils` | Utility functions | Data import/export helpers |
| `withr` | Temporary state changes | Testing and safe state modifications |

Install optional dependencies as needed:

``` r
# For mixed-effects models
install.packages(c("lme4", "coxme"))

# For table export (all formats)
install.packages(c("xtable", "flextable", "officer", "knitr", "tinytex"))

# For enhanced graphics
install.packages(c("ragg", "systemfonts"))

# For model diagnostics and selection
install.packages(c("MASS", "MuMIn", "pROC", "ResourceSelection"))
```

Or install all suggested packages at once:

``` r
install.packages("summata", dependencies = TRUE)
```

------------------------------------------------------------------------

## LaTeX Configuration

PDF and LaTeX export through
[`table2pdf()`](https://phmcc.github.io/summata/reference/table2pdf.md)
and
[`table2tex()`](https://phmcc.github.io/summata/reference/table2tex.md)
require a LaTeX distribution.

### TinyTeX (Recommended)

TinyTeX provides a lightweight, R-integrated LaTeX distribution:

``` r
install.packages("tinytex")
tinytex::install_tinytex()
```

Verify the installation:

``` r
tinytex::is_tinytex()
tinytex::tlmgr_version()
```

TinyTeX installs required LaTeX packages automatically on first use.

### Alternative LaTeX Distributions

Full LaTeX distributions may be installed system-wide:

| Distribution | Platform       | URL                            |
|:-------------|:---------------|:-------------------------------|
| TeX Live     | Cross-platform | <https://www.tug.org/texlive/> |
| MiKTeX       | Windows        | <https://miktex.org/>          |
| MacTeX       | macOS          | <https://www.tug.org/mactex/>  |

For Debian/Ubuntu systems:

``` bash
sudo apt install texlive-latex-base texlive-latex-extra texlive-fonts-recommended
```

### Required LaTeX Packages

The export functions utilize the following LaTeX packages:

| Category    | Packages                             |
|:------------|:-------------------------------------|
| Typography  | fontenc, inputenc, helvet            |
| Tables      | array, booktabs, longtable, colortbl |
| Layout      | geometry, pdflscape, lscape          |
| Graphics    | graphicx, xcolor                     |
| Specialized | standalone, varwidth, float, caption |

With TinyTeX, these packages install automatically. For standard
distributions:

``` bash
tlmgr install fontenc inputenc array booktabs longtable graphicx geometry \
  pdflscape lscape helvet standalone varwidth float caption xcolor colortbl
```

------------------------------------------------------------------------

## Microsoft Office Export

Export to Word, PowerPoint, and RTF formats requires:

``` r
install.packages(c("flextable", "officer"))
```

These packages enable creation of `.docx`, `.pptx`, and `.rtf` files
that can be opened and edited in Microsoft Office or compatible
applications.

------------------------------------------------------------------------

## Export Verification Procedure

Execute the following to verify complete functionality:

``` r
library(summata)

# Create test data
test_data <- data.frame(
  Variable = c("Sample Size", "Mean", "Standard Deviation"),
  Value = c("150", "45.3", "12.1")
)

# Test Word export (requires flextable, officer)
table2docx(test_data, file = "verification_test.docx")

# Test PDF export (requires xtable, LaTeX)
table2pdf(test_data, file = "verification_test.pdf")

# Clean up
file.remove("verification_test.docx", "verification_test.pdf")
```

Successful execution confirms proper configuration of all export-related
dependencies.

------------------------------------------------------------------------

## Troubleshooting

### LaTeX Not Found

If PDF generation fails, verify LaTeX accessibility:

``` r
Sys.which("pdflatex")
```

An empty return value indicates LaTeX is not in the system PATH.
Solutions include:

1.  Restart the R session after LaTeX installation
2.  Add the LaTeX binary directory to PATH manually
3.  Reinstall TinyTeX:
    [`tinytex::reinstall_tinytex()`](https://rdrr.io/pkg/tinytex/man/install_tinytex.html)

### Missing LaTeX Packages

For TinyTeX, install missing packages directly:

``` r
tinytex::tlmgr_install("package_name")
```

For compilation errors, examine the log file:

``` r
table2pdf(test_data, file = "debug.pdf", show_logs = TRUE)
```

### Package Installation Failures

Ensure repository access is configured:

``` r
options(repos = c(CRAN = "https://cloud.r-project.org"))
```

For packages requiring system libraries (Linux):

``` bash
sudo apt install r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev
```

### Write Permission Errors

Verify write access to the target directory:

``` r
getwd()
file.access(getwd(), mode = 2)  # Returns 0 if writable
```

------------------------------------------------------------------------

## System Requirements

| Component | Requirement |
|:---|:---|
| R | ≥ 4.2.0 |
| Operating System | Windows, macOS, or Linux |
| RAM | 4 GB minimum; 8 GB recommended |
| Disk Space | Base package (\< 5 MB), TinyTeX (~150 MB), Full LaTeX distribution (2–6 GB) |
