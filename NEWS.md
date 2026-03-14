# *summata* 0.11.4 (2026-03-14)

* Use profile likelihood CIs for GLM/negbin and exact *t*-distribution CIs for linear models, replacing Wald approximation
* Add `conf_method` parameter (`"profile"` / `"wald"`) to `fit()`, `uniscreen()`, `fullfit()`, `multifit()`, and `m2dt()`, with global option `summata.conf_method`
* Report complete-case *n*/Events in `fullfit()` multivariable rows (STROBE item 12)
* Cache profile likelihood CIs from `fit()` for reuse in forest plot functions
* Resolve `family = "Gamma"` string to `Gamma(link = "log")` for interpretable multiplicative effects
* Fix `survival::strata()` namespace in conditional logistic regression test
* Update benchmarks feature comparison articles

# *summata* 0.11.3 (2026-03-03)

* CRAN initial release re-submission
* Fix `table2pdf()` to specify output directory for files
* Minor WORDLIST changes

# *summata* 0.11.2 (2026-03-02)

* CRAN initial release re-submission
* Minimal executable examples added to all exported functions, with lengthier demonstrations wrapped in `\donttest{}`
* Modify file-writing examples in documentation headers to write to `tempdir()`

# *summata* 0.11.1 (2026-02-22)

* CRAN initial release re-submission
* Documentation LaTeX fixes
* Broken link fixes

# *summata* 0.11.0 (2026-02-21)

* CRAN initial release submission
* Upload documentation

# *summata* 0.10.9 (2026-02-20)

* Documentation headers fine editing and revisions
* Minor bugfixes

# *summata* 0.10.8 (2026-02-09)

* Add compatibility with US/EU locales and different big/small number marks
* Add quiet vs. verbose settings to regression functions for warnings
* Documentation and vignette updates

# *summata* 0.10.7 (2026-02-05)

* Add R-hub compatability
* Documentation and vignette updates

# *summata* 0.10.6 (2026-01-20)

* Fix "n" and "Events" counts in regression summary tables for extended GLM families
* Improve clogit example in "Regression Modeling" vignette
* Fix errors in `condense_table` logic in table export and forest plot functions
* Fix `variable_padding` in table export functions, ensure consistent behavior with `zebra_stripes`
* Lighten forest plot colors
* Change *n*/*N* labeling for consistency
* Expanded documentation examples for interaction effects, mixed effects, and extended GLM families
* Diversify use of `random` parameter in regression models
* Documentation cleaning

# *summata* 0.10.5 (2026-01-18)

* Add count outcomes (`ae_count` and `fu_count`) to `clintrial` mock dataset
* Fix color coding in forest plots
* Fix `glmforest()` to correctly extract values from `MASS::glm.nb()`
* Fix regression functions to respect level order in categorical functions
* Fix bug in `uniscreen()` where specified `p_threshold` values were not generating screened outputs
* Fixes to print outputs in regression functions
* Multiple `multifit()` revisions, including improved "n" and "Events" column handling and interaction effect formatting

# *summata* 0.10.4 (2026-01-15)

* Revise and standardize documentation
* Add linting fixes

# *summata* 0.10.3 (2026-01-12)

* Add vignettes
* Add GitHub Actions CI/CD workflows

# *summata* 0.10.2 (2026-01-12)

* Add support for various GLM families
* Add `MASS` dependency for `MASS::glm.nb()`
* Testing suite bugfixes
* Update documentation

# *summata* 0.10.1 (2026-01-05)

* Streamline intelligent `condense_table` string finding
* Update documentation

# *summata* 0.10.0 (2026-01-04)

* Add automatic font detection to forest plot functions
* Fix formatting for negative coefficients in forest plots
* Fix "negative zero" formatting issues
* Add safeguards to avoid multiple model types in `multifit()` calls
* Add intelligent handling of binary categorical variables for `condense_table`

# *summata* 0.9.4 (2026-01-02)

* Add `show_logs` parameter to `table2tex()`
* Consistent Unicode formatting for files
* Increased R version requirement to >= 4.2.0 for better Unicode compatability
* Add `survtable()` function with utilities
* Minor optimizations and function cleanup, with notation standardization

# *summata* 0.9.3 (2025-12-29)

* Package refresh, with revised naming scheme
* README rewrite

# *summata* 0.9.2 (2025-12-27)

* Expand "affirmative level" string criteria for binary categorical variables
* Add option to remove "Predictors" column in `multifit()` and `multiforest()`
* Update documentation

# *summata* 0.9.1 (2025-12-27)

* Ensure user-modifiable `p_digits` and `conf_level` parameters across all regression and forest plot functions
* Add parameter for toggling the footer in forest plot functions
* Update `uniscreen()` to accept mixed-effect models
* Update `autoforest()` to accept lmer and glmer objects

# *summata* 0.9.0 (2025-12-27)

* Add multivariate analysis support (`multifit()` and `multiforest()`)
* New vignette: "Multivariate Analysis"
* Modify core forest plot functions (`lmforest()`, `glmforest()`, and `coxforest()`) to accept *Summata* objects or models
* Add univariable screening forest plot (`uniforest()`)
* Rename `uscreen()` to `uniscreen()` for consistency
* Update `autoforest()` to handle new forest plot functions

# *summata* 0.8.3 (2025-12-24)

* Add optimizations to descriptive statistics workflows
* Add `p_per_stat` parameter to `desctable()`
* Change `add_reference_rows` parameter to just `reference_rows` in regression functions

# *summata* 0.8.2 (2025-12-23)

* Add optimizations to regression workflows

# *summata* 0.8.1 (2025-12-21)

* Add cleanup to documentation headers

# *summata* 0.8.0 (2025-12-20)

* Add pkgdown website
* Add testing suite
* Add NEWS
* Add README
* Add vignettes
* Add citation function

# *summata* 0.7.14 (2025-12-20)

* Add `autotable()` function

# *summata* 0.7.13 (2025-12-19)

* Rename `digits_p` parameter to `p_digits` and add to forest plot functions
* Fix `effect_label` parameter in forest plot functions
* Update/clean comments and documentation headers

# *summata* 0.7.12 (2025-12-15)

* Add `condense_quantitative` parameter to table export functions
* Modify `fit()` and `glmforest()` to correctly display Poisson model statistics
* Fix LaTeX/PDF export bugs

# *summata* 0.7.11 (2025-12-13)

* Update R version dependency (>= 4.4.0) for built-in null coalescing operator

# *summata* 0.7.10 (2025-12-08)

* Update MuMIn package dependency and QC statistics in lmer models
* Ensure that lme4 and coxme models work with `compfit()`

# *summata* 0.7.9 (2025-12-06)

* Fix "n" and "Events" columns for interaction effects
* Fix lmerMod compatibility with `fit()` and `m2dt()`
* Fix lmer compatibility for `lmforest()`
* Get interaction effects to show on forest plots

# *summata* 0.7.8 (2025-11-24)

* Add GLM and Cox mixed-effect compatibility with `m2dt()` and `fit()`
* Add GLM and Cox mixed-effect compatibility with `glmforest()` and `coxforest()`

# *summata* 0.7.7 (2025-11-20)

* Add data requirement to `m2dt()` to allow for accurate per-group "n" and "events" columns for all models
* Multiple `R CMD check` fixes

# *summata* 0.7.6 (2025-11-12)

* Rename `var_labels` parameter to just `labels`

# *summata* 0.7.5 (2025-11-12)

* Add compatibility with interaction terms

# *summata* 0.7.4 (2025-11-03)

* Forest plot unit conversion fixes

# *summata* 0.7.3 (2025-11-03)

* Add coefficient table combination function for `compfit()`
* Global and imported variable fixes

# *summata* 0.7.2 (2025-11-03)

* Multiple `R CMD check` fixes

# *summata* 0.7.1 (2025-11-03)

* Fix global variables

# *summata* 0.7.0 (2025-11-02)

* Add global variables

# *summata* 0.6.3 (2025-11-03)

* Add suggested dependencies for vignettes and documentation

# *summata* 0.6.2 (2025-11-03)

* Add GitHub Actions workflow

# *summata* 0.6.1 (2025-11-03)

* Update internal documentation for `flextable` export

# *summata* 0.6.0 (2025-11-02)

* Add `flextable` table export functions with helpers
* Documentation updates for existing table export functions

# *summata* 0.5.12 (2025-10-31)

* Add `data.table`-specific performance enhancements for core functions

# *summata* 0.5.11 (2025-10-29)

* Update documentation for core functions

# *summata* 0.5.10 (2025-10-29)

* Indent groups and condense tables logic refinements

# *summata* 0.5.9 (2025-10-26)

* Add zebra stripes and dark header options for all table outputs

# *summata* 0.5.8 (2025-10-23)

* Update forest plot functions: additive columns, zebra stripes, and units conversion

# *summata* 0.5.7 (2025-10-18)

* Add condensing/indenting to forest plots

# *summata* 0.5.6 (2025-10-12)

* Add table condensing to table export functions

# *summata* 0.5.5 (2025-10-11)

* Add `lmforest()` and `autoforest()` functions

# *summata* 0.5.4 (2025-10-11)

* Add `clintrial` sample data

# *summata* 0.5.3 (2025-10-10)

* Fix `desctable()` ordering to follow variable levels
* Add "events" columns in `glmforest()` and `coxforest()`
* Fix outcome factor handling in `glmforest()` and `coxforest()`
* Implement exponentiation options for "fit" functions

# *summata* 0.5.2 (2025-10-09)

* Fix "n" and "events" columns in `uscreen()`, `fit()`, and `*summata*()`
* Add hyphen space fillers for "p-value", "Uni p", and "Multi p" columns

# *summata* 0.5.1 (2025-10-08)

* Add N rows to `desctable()`
* Fix errors with "rolling" p-values in `desctable()`
* Fix `m2dt()` edge cases

# *summata* 0.5.0 (2025-10-06)

* Add export functions `table2pdf()`, `table2tex()`, and `table2html()`

# *summata* 0.4.1 (2025-10-06)

* Added `compfit()` function
* Reorganized helper functions
* Function name cleanup

# *summata* 0.4.0 (2025-10-05)

* Replacement of `mmodel()` function with `fit()`
* Standardization of output formats

# *summata* 0.3.1 (2025-10-02)

* Fixes for unknown/missing rows, reporting of ranges, and non-grouped tables in `desctable()`

# *summata* 0.3.0 (2025-09-30)

* Edits to `desctable()` and `fit()`
* Addition of "raw data" outputs
* Expand internal documentation

# *summata* 0.2.1 (2025-09-26)

* Output formatting edits

# *summata* 0.2.0 (2025-09-26)

* Add `desctable()` function
* Add `fit()` function

# *summata* 0.1.0 (2025-09-23)

* Initial commit
