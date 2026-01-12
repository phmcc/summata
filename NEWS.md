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

* Add `autoexport()` function

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

* Fixed "n" and "events" columns in `uscreen()`, `fit()`, and `*summata*()`
* Added hyphen space fillers for "p-value", "Uni p", and "Multi p" columns

# *summata* 0.5.1 (2025-10-08)

* Added N rows to `desctable()`
* Fixed errors with "rolling" p-values in `desctable()`
* Fixed `m2dt()` edge cases

# *summata* 0.5.0 (2025-10-06)

* Added export functions `table2pdf()`, `table2tex()`, and `table2html()`

# *summata* 0.4.1 (2025-10-06)

* Added `compfit()` function
* Reorganized helper functions
* Cleaned up function names

# *summata* 0.4.0 (2025-10-05)

* Replacement of `mmodel()` function with `fit()`
* Standardization of output formats

# *summata* 0.3.1 (2025-10-02)

* Fixes for unknown/missing rows, reporting of ranges, and non-grouped tables in `desctable()`

# *summata* 0.3.0 (2025-09-30)

* Edits to `desctable()` and `*summata*()`
* Addition of "raw data" outputs
* Expanded internal documentation

# *summata* 0.2.1 (2025-09-26)

* Output formatting edits

# *summata* 0.2.0 (2025-09-26)

* Added `desctable()` function
* Added `*summata*()` function

# *summata* 0.1.0 (2025-09-23)

* Initial commit
