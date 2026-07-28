## Shared expectation helpers.
##
## Helpers used by more than one test file live here so that a single definition
## governs the contract. Helpers specific to one file remain in that file.
##
## Objects defined here are visible to every test file. A test file that assigns
## to the same name shadows this definition for that file only.


## Check that a formatted result carries an effect column, optionally of a
## specific type.
##
## Applies to the display tables produced by fit(), fullfit(), uniscreen() and
## multifit(), whose effect columns combine the estimate with its confidence
## interval (for example, "OR (95% CI)"). For the unformatted columns returned
## by m2dt(), see expect_raw_effect_column() in test-m2dt.R.
expect_effect_column <- function(result, effect_type = NULL) {
    col_names <- names(result)

    ## Find effect column (contains OR, HR, RR, Coefficient, or Estimate with CI)
    effect_cols <- grep("(OR|HR|RR|Coefficient|Estimate).*CI", col_names, value = TRUE)
    expect_true(length(effect_cols) >= 1,
                info = paste("Expected effect column, found:", paste(col_names, collapse = ", ")))

    if (!is.null(effect_type)) {
        ## For linear models, both "Coefficient" and "Estimate" are valid
        if (effect_type == "Coefficient") {
            expect_true(any(grepl("Coefficient|Estimate", effect_cols)),
                        info = paste("Expected Coefficient or Estimate in columns, found:",
                                     paste(effect_cols, collapse = ", ")))
        } else {
            expect_true(any(grepl(effect_type, effect_cols)),
                        info = paste("Expected", effect_type, "in columns, found:",
                                     paste(effect_cols, collapse = ", ")))
        }
    }
}


## Check that a formatted result carries a p-value column and that its entries
## are formatted as expected.
##
## Accepted forms are a decimal value, a threshold below the reporting precision,
## an em-dash placeholder for reference rows, and an empty string for rows that
## carry no test. The threshold pattern is deliberately written as 0\\.0*1 rather
## than 0\\.001, because the reporting precision is configurable through
## p_digits and the threshold moves with it.
expect_pvalue_column <- function(result) {
    expect_true("p-value" %in% names(result))

    pvals <- result$`p-value`
    valid_pvals <- grepl("^[0-9]\\.[0-9]+$|^< 0\\.0*1$|^-$|^$", pvals)
    expect_true(all(valid_pvals),
                info = paste("Invalid p-values found:",
                             paste(pvals[!valid_pvals], collapse = ", ")))
}


## Load a pristine copy of the packaged dataset.
##
## Reading the object directly picks up whatever copy is nearest in scope, which
## may be a test file's own augmented version or, in an interactive session, a
## modified copy in the global environment. Tests that depend on the dataset
## being unmodified should use this instead.
fresh_clintrial <- function() {
    env <- new.env(parent = emptyenv())
    utils::data("clintrial", package = "summata", envir = env)
    as.data.frame(env$clintrial)
}


## Return clintrial with the first n_missing values of the named outcome set to
## NA, as a data.table. Used by the tests covering complete-case counting.
with_missing_outcome <- function(outcome, n_missing = 100L) {
    df <- fresh_clintrial()
    df[[outcome]][seq_len(n_missing)] <- NA
    data.table::as.data.table(df)
}
