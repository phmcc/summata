## get_model_type_name() labels generalized linear models by family. Families
## without a fixed label are capitalized in base R, so that no suggested
## package is required.

test_that("GLM families receive readable labels", {

    df <- fresh_clintrial()
    label <- function(family) {
        model <- stats::glm(recovery_days ~ age, data = df, family = family)
        summata:::get_model_type_name(model)
    }

    expect_identical(label(stats::Gamma(link = "log")), "Gamma")
    expect_identical(label(stats::inverse.gaussian(link = "log")), "Inverse Gaussian")
    expect_identical(label(stats::quasi(link = "log", variance = "mu")), "Quasi GLM")
})
