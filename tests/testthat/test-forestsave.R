## forestsave() writes a forest plot at the dimensions the plotting function
## recommended, so that a figure is not silently sized to the device default.

forest_fixture <- function() {
    dt <- data.table::as.data.table(fresh_clintrial())
    model <- stats::glm(readmission_30d ~ age + sex + stage,
                        data = dt, family = stats::binomial)
    suppressMessages(suppressWarnings(glmforest(model, data = dt)))
}


test_that("forest plots record the units their dimensions are expressed in", {

    dt <- data.table::as.data.table(fresh_clintrial())
    model <- stats::glm(readmission_30d ~ age + sex,
                        data = dt, family = stats::binomial)

    p_in <- suppressMessages(suppressWarnings(
        glmforest(model, data = dt, units = "in")))
    p_cm <- suppressMessages(suppressWarnings(
        glmforest(model, data = dt, units = "cm")))

    expect_equal(attr(p_in, "rec_dims")$units, "in")
    expect_equal(attr(p_cm, "rec_dims")$units, "cm")

    ## The centimetre figures are the inch figures scaled, so a consumer that
    ## assumed inches would mis-size the output by the conversion factor
    expect_equal(attr(p_cm, "rec_dims")$width,
                 attr(p_in, "rec_dims")$width * 2.54,
                 tolerance = 1e-6)
})


test_that("forestsave writes a file at the recommended dimensions", {

    p <- forest_fixture()
    dims <- attr(p, "rec_dims")
    out <- tempfile(fileext = ".pdf")

    expect_message(forestsave(p, out),
                   sprintf("width = %.1f", dims$width))
    expect_true(file.exists(out))
    expect_gt(file.size(out), 0)
})


test_that("forestsave returns its path invisibly and can be silenced", {

    p <- forest_fixture()
    out <- tempfile(fileext = ".pdf")

    expect_no_message(result <- forestsave(p, out, quiet = TRUE))
    expect_identical(result, out)
})


test_that("forestsave converts when the requested units differ", {

    dt <- data.table::as.data.table(fresh_clintrial())
    model <- stats::glm(readmission_30d ~ age + sex,
                        data = dt, family = stats::binomial)

    ## Recommended in inches, requested in centimetres
    p <- suppressMessages(suppressWarnings(
        glmforest(model, data = dt, units = "in")))
    expected <- attr(p, "rec_dims")$width * 2.54

    expect_message(forestsave(p, tempfile(fileext = ".pdf"), units = "cm"),
                   sprintf("width = %.1f cm", expected))
})


test_that("forestsave requires dimensions when none were recommended", {

    plain <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::geom_point()

    ## Writing at the device default is the mis-sizing the function prevents
    expect_error(forestsave(plain, tempfile(fileext = ".pdf")),
                 "no recommended dimensions")

    ## Explicit dimensions are accepted
    out <- tempfile(fileext = ".pdf")
    expect_silent(forestsave(plain, out, width = 6, height = 4, quiet = TRUE))
    expect_true(file.exists(out))
})


test_that("forestsave validates its arguments", {

    p <- forest_fixture()

    expect_error(forestsave("not a plot", tempfile(fileext = ".pdf")),
                 "must be a forest plot")
    expect_error(forestsave(p, character(0)), "non-empty character")
    expect_error(forestsave(p, tempfile(fileext = ".pdf"), units = "furlongs"),
                 "'units' must be one of")
})


test_that("PDF output defers to R's internal device rather than Cairo", {

    ## Cairo embeds fonts but can substitute an incorrect italic face, which
    ## the italic footer would expose. It remains available as an override.
    expect_null(summata:::select_forest_device("x.pdf"))
})


test_that("raster output prefers ragg where available", {

    skip_if_not_installed("ragg")

    expect_identical(summata:::select_forest_device("x.png"), ragg::agg_png)
    expect_identical(summata:::select_forest_device("x.tiff"), ragg::agg_tiff)
    expect_identical(summata:::select_forest_device("x.jpeg"), ragg::agg_jpeg)
})


test_that("formats with no preferred device defer to ggsave", {

    expect_null(summata:::select_forest_device("x.svg"))
    expect_null(summata:::select_forest_device("x"))
})


## capabilities("cairo") reports that the library was compiled in, not that the
## device can be opened. On a headless macOS runner the two disagree: the
## capability is TRUE while cairo_pdf() fails for want of the supporting
## infrastructure. The device is therefore probed rather than asked.
cairo_pdf_works <- function() {
    if (!isTRUE(unname(capabilities("cairo")))) {
        return(FALSE)
    }
    probe <- tempfile(fileext = ".pdf")
    ok <- tryCatch({
        grDevices::cairo_pdf(filename = probe)
        grDevices::dev.off()
        file.exists(probe) && file.size(probe) > 0
    }, error = function(e) FALSE, warning = function(w) FALSE)
    unlink(probe)
    return(isTRUE(ok))
}


test_that("an explicit device overrides the automatic selection", {

    skip_if_not(cairo_pdf_works(), "Cairo PDF device is not usable here")

    p <- forest_fixture()

    ## Both the device function and the shorthand reach the same device
    out_fn <- tempfile(fileext = ".pdf")
    expect_silent(forestsave(p, out_fn, device = grDevices::cairo_pdf,
                             quiet = TRUE))
    expect_true(file.exists(out_fn))

    out_str <- tempfile(fileext = ".pdf")
    expect_silent(forestsave(p, out_str, device = "cairo", quiet = TRUE))
    expect_true(file.exists(out_str))
})


test_that("the cairo shorthand resolves per format and rejects raster", {

    ## This test only inspects the returned device function, so the capability
    ## flag is the correct guard; the device is never opened.
    skip_if(!isTRUE(unname(capabilities("cairo"))), "Cairo not compiled in")

    expect_identical(summata:::cairo_device("x.pdf"), grDevices::cairo_pdf)
    expect_identical(summata:::cairo_device("x.eps"), grDevices::cairo_ps)
    expect_identical(summata:::cairo_device("x.svg"), grDevices::svg)

    ## ragg supersedes Cairo's raster backends, so the shorthand does not
    ## silently fall through to something else
    expect_error(summata:::cairo_device("x.png"), "vector formats")
})


test_that("embed_fonts is rejected for non-PDF output", {

    p <- forest_fixture()
    out <- tempfile(fileext = ".png")

    expect_warning(forestsave(p, out, embed_fonts = TRUE, quiet = TRUE),
                   "PDF output only")
    ## The plot is still written; only the embedding step is skipped
    expect_true(file.exists(out))
})


test_that("embed_fonts degrades to a warning when Ghostscript is absent", {

    p <- forest_fixture()
    out <- tempfile(fileext = ".pdf")

    if (nzchar(tools::find_gs_cmd())) {
        expect_silent(forestsave(p, out, embed_fonts = TRUE, quiet = TRUE))
    } else {
        expect_warning(forestsave(p, out, embed_fonts = TRUE, quiet = TRUE),
                       "Ghostscript")
    }

    ## Either way the figure exists, since embedding runs after the write
    expect_true(file.exists(out))
    expect_gt(file.size(out), 0)
})
