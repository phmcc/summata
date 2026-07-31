## recdims() reports the dimensions the plotting function computed, and must
## agree with what forestsave() applies. A discrepancy between the two would
## mean a figure sized one way and reported another.

recdims_fixture <- function(units = "in") {
    dt <- data.table::as.data.table(fresh_clintrial())
    model <- stats::glm(readmission_30d ~ age + sex + stage,
                        data = dt, family = stats::binomial)
    suppressMessages(suppressWarnings(
        glmforest(model, data = dt, units = units)))
}


test_that("recdims reports the attached dimensions unchanged", {

    p <- recdims_fixture()
    dims <- attr(p, "rec_dims")

    result <- recdims(p)

    expect_named(result, c("width", "height"))
    expect_equal(unname(result[["width"]]), dims$width)
    expect_equal(unname(result[["height"]]), dims$height)
    expect_identical(attr(result, "units"), dims$units)
})


test_that("recdims converts where the requested units differ", {

    p <- recdims_fixture(units = "in")
    inches <- recdims(p)

    mm <- recdims(p, units = "mm")

    expect_identical(attr(mm, "units"), "mm")
    expect_equal(unname(mm[["width"]]),
                 unname(inches[["width"]]) * 25.4,
                 tolerance = 1e-6)

    ## A round trip returns the original values
    back <- recdims(p, units = "in")
    expect_equal(unname(back[["width"]]), unname(inches[["width"]]))
})


test_that("recdims defaults to the units the plot was created in", {

    p_cm <- recdims_fixture(units = "cm")

    expect_identical(attr(recdims(p_cm), "units"), "cm")

    ## Requesting inches from a centimetre plot converts rather than relabels
    in_dims <- recdims(p_cm, units = "in")
    expect_equal(unname(in_dims[["width"]]),
                 unname(recdims(p_cm)[["width"]]) / 2.54,
                 tolerance = 1e-6)
})


test_that("recdims and forestsave agree on the dimensions used", {

    p <- recdims_fixture()
    dims <- recdims(p)

    out <- tempfile(fileext = ".pdf")

    ## forestsave() reports the dimensions it applies, which must be the
    ## figures recdims() returns
    expect_message(forestsave(p, out),
                   sprintf("width = %.1f", dims[["width"]]))
})


test_that("recdims validates its arguments", {

    p <- recdims_fixture()

    expect_error(recdims("not a plot"), "must be a forest plot")
    expect_error(recdims(p, units = "furlongs"), "'units' must be one of")

    ## A plot from another source carries no recommendation
    plain <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::geom_point()
    expect_error(recdims(plain), "no recommended dimensions")
})
