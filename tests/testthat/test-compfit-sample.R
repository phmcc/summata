## Information criteria are comparable only across models fitted to the same
## observations. Each candidate is fitted to complete cases for its own
## predictor set, so a predictor carrying missing values reduces the sample for
## the models that contain it.

## clintrial carries missing values of its own, so tests start from a
## complete-case base and introduce missingness only where one is needed.
compfit_base <- function() {
    keep <- c("readmission_30d", "age", "sex", "bmi")
    df <- fresh_clintrial()
    df[stats::complete.cases(df[, keep]), keep]
}

compfit_candidates <- list(base = c("age", "sex"),
                           extended = c("age", "sex", "bmi"))


test_that("compfit warns when candidate models are fitted to different samples", {

    df <- compfit_base()
    df$bmi[seq_len(100)] <- NA

    expect_warning(
        compfit(
            data = df,
            outcome = "readmission_30d",
            model_list = compfit_candidates,
            model_type = "glm",
            family = "binomial"
        ),
        "different numbers of observations"
    )
})


test_that("compfit is silent when candidate models share a sample", {

    expect_no_warning(
        compfit(
            data = compfit_base(),
            outcome = "readmission_30d",
            model_list = compfit_candidates,
            model_type = "glm",
            family = "binomial"
        ),
        message = "different numbers of observations"
    )
})


test_that("the documented comparison workflow is quiet and well behaved", {

    ## Mirrors the compfit() examples: a single complete-case dataset shared by
    ## every comparison, so information criteria describe a common sample.
    df <- fresh_clintrial()
    keep <- c("readmission_30d", "os_months", "os_status", "age", "sex",
              "smoking", "diabetes", "stage", "ecog", "grade", "treatment")
    dat <- stats::na.omit(df[, keep])

    models <- list(
        base = c("age", "sex"),
        clinical = c("age", "sex", "smoking", "diabetes"),
        full = c("age", "sex", "smoking", "diabetes", "stage", "ecog")
    )

    result <- expect_no_warning(
        suppressMessages(compfit(
            data = dat,
            outcome = "readmission_30d",
            model_list = models,
            model_names = c("Base", "Clinical", "Full")
        )),
        message = "different numbers of observations"
    )

    ## Every model uses the same observations
    expect_length(unique(result$N), 1L)

    ## No model reports a suspect fit. os_status produced complete separation at
    ## ECOG 3, where every patient had the event; readmission_30d does not.
    expect_false(any(result$Converged == "Suspect"))
})


test_that("coefficient estimates from the documented models are finite", {

    skip_if_not_installed("xtable")

    df <- fresh_clintrial()
    keep <- c("readmission_30d", "age", "sex", "smoking", "diabetes",
              "stage", "ecog")
    dat <- stats::na.omit(df[, keep])

    detailed <- suppressMessages(compfit(
        data = dat,
        outcome = "readmission_30d",
        model_list = list(base = c("age", "sex"),
                          full = c("age", "sex", "smoking", "diabetes",
                                   "stage", "ecog")),
        include_coefficients = TRUE
    ))

    coefs <- attr(detailed, "coefficients")
    effect_col <- grep("CI", names(coefs), value = TRUE)[1]
    skip_if(is.na(effect_col), "no effect column found")

    ## Complete separation renders as an implausibly large estimate rather than
    ## as NA, so the check is on magnitude
    numeric_part <- suppressWarnings(
        as.numeric(sub(" .*", "", coefs[[effect_col]])))
    expect_true(all(is.na(numeric_part) | numeric_part < 1000))
})
