# Development Roadmap

This article outlines planned features for future versions of `summata`,
informed by capabilities available in related packages and user
feedback. The roadmap is tentative and subject to change based on
development priorities and community input.

------------------------------------------------------------------------

## Planned Features

The following features are under consideration for future releases,
organized by priority level.

### Higher Priority

Features expected in near-term releases:

| Feature | Description | Reference |
|:---|:---|:---|
| Paired/longitudinal tables | Summary tables for paired measurements across time points | `arsenal` |
| Subgroup analysis | Stratified analyses with syntax similar to `tbl_strata()` | `gtsummary` |
| Missing data visualization | Visual patterns of missingness for data quality assessment | `finalfit` |

### Medium Priority

Features planned for subsequent releases:

| Feature | Description | Reference |
|:---|:---|:---|
| Inline text reporting | Extract statistics for R Markdown inline reporting | `gtsummary` |
| Survey-weighted analysis | Support for `svyglm()`, `svycoxph()` from the survey package | `gtsummary`, `survey` |
| Competing risks regression | Fine-Gray models via the `cmprsk` package | `finalfit` |
| Multinomial regression | Support for [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) with multi-level categorical outcomes | `gtsummary`, `nnet` |
| Ordinal regression | Support for [`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) and [`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) for ordered categorical outcomes | `MASS`, `ordinal` |
| Custom statistics | User-defined summary statistics functions for [`desctable()`](https://phmcc.codeberg.page/summata/reference/desctable.md) | `gtsummary` |

### Lower Priority

Features under consideration for future releases:

| Feature | Description | Reference |
|:---|:---|:---|
| Bootstrap CIs | Bootstrapped confidence intervals for model predictions | `finalfit` |
| GAM models | Generalized additive model support via `mgcv` | `gtsummary` |
| Bayesian models | Support for `brmsfit` and `stanreg` model objects | `gtsummary` |

------------------------------------------------------------------------

## Contributing

Contributions to `summata` are welcome. If interested in implementing a
roadmap feature or proposing a new one, please consult the contributing
guidelines in the [package
repository](https://codeberg.org/phmcc/summata).

### Bug Reports

Bug reports and feature requests may be submitted via the issue tracker,
either on [Codeberg](https://codeberg.org/phmcc/summata/issues) or
[GitHub](https://github.com/phmcc/summata/issues).

### Development Repository

- **Primary development**:
  [codeberg.org/phmcc/summata](https://codeberg.org/phmcc/summata)
- **GitHub mirror**:
  [github.com/phmcc/summata](https://github.com/phmcc/summata)

------------------------------------------------------------------------

## Version History

See the
[Changelog](https://phmcc.codeberg.page/summata/articles/news/index.md)
for a detailed history of changes in each release.

------------------------------------------------------------------------

## Additional Resources

- [Feature
  Comparison](https://phmcc.codeberg.page/summata/articles/feature_comparison.md):
  Comparison with related packages
- [Benchmarks](https://phmcc.codeberg.page/summata/articles/benchmarks.md):
  Performance comparisons
- [gtsummary documentation](https://www.danieldsjoberg.com/gtsummary/)
- [finalfit documentation](https://finalfit.org/)
- [arsenal documentation](https://mayoverse.github.io/arsenal/)
