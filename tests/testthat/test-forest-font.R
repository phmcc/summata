## The forest plots choose their default font in detect_plot_font(), which
## checks the knitr graphics device so that a plot rendered by knitr with a
## ragg device can use a system font. knitr is only suggested, so the check
## must not require it; when knitr was not installed, every forest plot
## failed with "there is no package called 'knitr'".

test_that("forest plots do not load knitr", {

    skip_if(isNamespaceLoaded("knitr"), "knitr is already loaded")
    df <- fresh_clintrial()

    expect_type(summata:::detect_plot_font(), "character")

    model <- stats::glm(os_status ~ age + sex, data = df, family = stats::binomial)
    p <- suppressMessages(suppressWarnings(glmforest(model, data = df)))
    expect_s3_class(p, "ggplot")

    expect_false(isNamespaceLoaded("knitr"))
})
