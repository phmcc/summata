### * Main function

#' Save a Forest Plot at Its Recommended Dimensions
#'
#' Writes a forest plot to file, sized by default from the recommendations the
#' plotting function computed from the plot's own content, and rendered through
#' a graphics device chosen to reproduce the fonts seen on screen.
#'
#' @param plot A forest plot produced by \code{lmforest()}, \code{glmforest()},
#'   \code{coxforest()}, \code{uniforest()}, \code{multiforest()}, or
#'   \code{autoforest()}. Any \pkg{ggplot2} object is accepted, but one that
#'   carries no recommended dimensions requires \code{width} and \code{height}
#'   to be supplied.
#' @param file Character string giving the output path. The file extension
#'   determines the format: \code{.pdf}, \code{.png}, \code{.tiff},
#'   \code{.jpeg}, \code{.svg}, \code{.eps}, and the other formats
#'   \code{ggplot2::ggsave()} supports.
#' @param width,height Numeric plot dimensions. Taken from the plot's
#'   recommended dimensions when not supplied, and converted where these differ
#'   from \code{units}.
#' @param units Character string giving the units for \code{width} and
#'   \code{height}: \code{"in"} (inches), \code{"cm"}, or \code{"mm"}. Defaults
#'   to the units the recommended dimensions were expressed in, and to
#'   \code{"in"} where none were recorded.
#' @param dpi Numeric resolution in dots per inch for raster formats. Ignored
#'   for vector formats. Default is \code{300}.
#' @param device Graphics device. Selected automatically when \code{NULL}: see
#'   Details. Supply \code{"cairo"} to use the Cairo device for the output
#'   format, which embeds fonts, or any device function or name recognized by
#'   \code{ggplot2::ggsave()}.
#' @param embed_fonts Logical. Embed the fonts in a PDF after it is written,
#'   using \code{grDevices::embedFonts()}. Requires Ghostscript, and applies
#'   only to PDF output. Default is \code{FALSE}.
#' @param quiet Logical. Suppress the message reporting the dimensions used.
#'   Default is \code{FALSE}.
#' @param ... Additional arguments passed to \code{ggplot2::ggsave()}, such as
#'   \code{bg} or \code{limitsize}.
#'
#' @return Invisibly returns \code{file}.
#'
#' @details
#' Forest plots are sized by their content rather than by the page: a plot with
#' more variables, longer labels, or more factor levels requires a larger
#' canvas, and the plotting functions compute one and attach it to the returned
#' object as a \code{"rec_dims"} attribute. \code{forestsave()} reads that
#' attribute, so the following are equivalent:
#'
#' \preformatted{
#' dims <- attr(p, "rec_dims")
#' ggplot2::ggsave("forest.pdf", p, width = dims$width,
#'                 height = dims$height, units = dims$units)
#'
#' forestsave(p, "forest.pdf")
#' }
#'
#' \strong{Device selection.} Fonts are measured differently by different
#' graphics devices, so a figure written through one device may not reproduce
#' the spacing seen on screen or in the package documentation. Unless
#' \code{device} is supplied, \code{forestsave()} selects:
#'
#' \describe{
#'   \item{PDF}{R's internal PDF device, \code{grDevices::pdf()}, reached by
#'     deferring to \code{ggplot2::ggsave()}. Cairo is not used by default:
#'     although it embeds fonts, its face selection can substitute an
#'     incorrect italic, which matters here because forest plot footers and
#'     statistical notation are set in italic. Pass \code{device = "cairo"}
#'     to use it where font embedding is required.}
#'   \item{PNG, TIFF, JPEG}{The corresponding \pkg{ragg} device where that
#'     package is installed, whose text metrics match those used to render the
#'     package vignettes. Falls back to the standard raster devices
#'     otherwise.}
#'   \item{Other formats}{Deferred to \code{ggplot2::ggsave()}.}
#' }
#'
#' \strong{Font embedding.} R's PDF device draws on the standard fonts that every
#' PDF reader is required to provide, so its output is portable without
#' embedding. Embedding matters where the plot uses a font outside that set, and
#' where a publisher requires it. Setting \code{embed_fonts = TRUE} embeds the
#' fonts after the file is written, through \code{grDevices::embedFonts()}, which
#' requires a Ghostscript installation. This is offered as an alternative to
#' \code{device = "cairo"}, which embeds fonts as it draws but selects its own
#' italic face.
#'
#' \code{device = "cairo"} is a shorthand of this package's own, resolving to
#' \code{grDevices::cairo_pdf()}, \code{grDevices::cairo_ps()}, or
#' \code{grDevices::svg()} according to the output format. It applies to vector
#' formats only: raster output already uses the \pkg{ragg} devices, which
#' supersede Cairo's raster backends. Any other value of \code{device} is
#' passed to \code{ggplot2::ggsave()} unchanged.
#'
#' The dimensions used are reported through a \code{message()} unless
#' \code{quiet = TRUE}, so that a figure written at an unexpected size is
#' apparent at the point it is written.
#'
#' @seealso \code{\link{autoforest}}, \code{\link{lmforest}},
#'   \code{\link{glmforest}}, \code{\link{coxforest}}, \code{\link{uniforest}},
#'   \code{\link{multiforest}} for producing the plots;
#'   \code{\link{tablesave}} for the equivalent operation on tables
#'
#' @examples
#' data(clintrial)
#' data(clintrial_labels)
#'
#' model <- stats::glm(readmission_30d ~ age + sex + stage,
#'                     data = clintrial, family = stats::binomial)
#'
#' p <- glmforest(model, data = clintrial, labels = clintrial_labels)
#'
#' \donttest{
#' # Example 1: Save at the recommended dimensions
#' # Written under tempdir() so the example respects CRAN's no-write policy;
#' # in practice any desired path may be supplied.
#' forestsave(p, file.path(tempdir(), "forest.pdf"))
#'
#' # Example 2: Raster output at publication resolution
#' forestsave(p, file.path(tempdir(), "forest.png"), dpi = 600)
#'
#' # Example 3: Dimensions given explicitly, in any supported units
#' forestsave(p, file.path(tempdir(), "forest_sized.pdf"),
#'            width = 24, height = 16, units = "cm")
#'
#' # Example 4: Embed fonts after writing, where Ghostscript is available
#' forestsave(p, file.path(tempdir(), "forest_embedded.pdf"),
#'            embed_fonts = TRUE)
#'
#' # Example 5: Cairo embeds fonts as it draws, at the cost of its own
#' # italic face selection
#' forestsave(p, file.path(tempdir(), "forest_cairo.pdf"), device = "cairo")
#' }
#'
#' @family visualization functions
#' @export
forestsave <- function(plot, file, width = NULL, height = NULL,
                       units = NULL, dpi = 300, device = NULL,
                       embed_fonts = FALSE, quiet = FALSE, ...) {

    if (!inherits(plot, "ggplot")) {
        stop("'plot' must be a forest plot or other ggplot2 object.",
             call. = FALSE)
    }

    if (!is.character(file) || length(file) != 1 || !nzchar(file)) {
        stop("'file' must be a single non-empty character string.",
             call. = FALSE)
    }

    dims <- attr(plot, "rec_dims")

    ## The units the recommended dimensions are expressed in. Plots produced
    ## before this element was recorded, and plots from other sources, are
    ## treated as inches.
    dims_units <- if (!is.null(dims$units)) dims$units else "in"

    ## Output units default to those of the recommendation, so that the
    ## dimensions pass through unconverted in the common case
    if (is.null(units)) {
        units <- if (!is.null(dims)) dims_units else "in"
    }

    valid_units <- c("in", "cm", "mm", "px", "pt")

    if (!units %in% valid_units) {
        stop("'units' must be one of: ",
             paste(sQuote(valid_units), collapse = ", "), ".",
             call. = FALSE)
    }

    ## Fill either dimension from the recommendation, converting where the
    ## requested units differ from those the recommendation was computed in
    if (is.null(width) && !is.null(dims$width)) {
        width <- convert_units(dims$width, from = dims_units, to = units)
    }

    if (is.null(height) && !is.null(dims$height)) {
        height <- convert_units(dims$height, from = dims_units, to = units)
    }

    ## A plot carrying no recommendation and given no dimensions would be
    ## written at the device default, which is the mis-sizing this function
    ## exists to prevent
    if (is.null(width) || is.null(height)) {
        stop("'plot' carries no recommended dimensions, so 'width' and ",
             "'height' must be supplied.", call. = FALSE)
    }

    if (is.null(device)) {
        device <- select_forest_device(file)
    } else if (identical(device, "cairo")) {
        ## A shorthand of this package's own, resolved here because
        ## ggplot2::ggsave() does not recognize it
        device <- cairo_device(file)
    }

    if (!quiet) {
        ## Phrased to match the table export functions, with the dimensions
        ## appended because an unexpected size is the failure this reports
        message(sprintf(
            "Forest plot saved to %s (width = %.1f %s, height = %.1f %s)",
            file, width, units, height, units))
    }

    ## ggsave() accepts NULL for device, in which case it infers one from the
    ## file extension
    ggplot2::ggsave(filename = file, plot = plot,
                    width = width, height = height, units = units,
                    dpi = dpi, device = device, ...)

    ## A graphics device that cannot be opened may warn rather than raise an
    ## error, leaving no file behind. Returning the path in that case would
    ## report a success that did not occur.
    if (!file.exists(file)) {
        stop("The graphics device wrote no file to '", file, "'.",
             if (!is.null(device)) {
                 paste0(" The device supplied through 'device' may be ",
                        "unavailable on this system.")
             } else {
                 ""
             }, call. = FALSE)
    }

    if (isTRUE(embed_fonts)) {
        embed_plot_fonts(file, quiet = quiet)
    }

    invisible(file)
}


### * Graphics device selection

#' Select a graphics device for forest plot output
#'
#' Chooses a device whose text metrics reproduce those used to render the
#' package documentation, falling back to \code{ggplot2::ggsave()}'s own
#' selection where the preferred device is unavailable.
#'
#' @param file Character string giving the output path.
#' @return A device function, or \code{NULL} to defer to
#'   \code{ggplot2::ggsave()}.
#' @keywords internal
select_forest_device <- function(file) {

    ext <- tolower(tools::file_ext(file))

    has_ragg <- requireNamespace("ragg", quietly = TRUE)

    device <- switch(
        ext,
        ## PDF defers to ggsave, and so to R's internal device. Cairo
        ## embeds fonts but can substitute an incorrect italic face, which
        ## the italic footer and statistical notation would expose.
        pdf  = NULL,
        png  = if (has_ragg) ragg::agg_png else NULL,
        tif  = ,
        tiff = if (has_ragg) ragg::agg_tiff else NULL,
        jpg  = ,
        jpeg = if (has_ragg) ragg::agg_jpeg else NULL,
        NULL
    )

    return(device)
}


#' Resolve the Cairo device for an output format
#'
#' Maps the \code{"cairo"} shorthand accepted by \code{forestsave()} onto the
#' Cairo-backed device for the requested format. Cairo's raster backends are
#' superseded by \pkg{ragg} and are not offered, so the shorthand applies to
#' vector formats only.
#'
#' @param file Character string giving the output path.
#' @return A device function.
#' @keywords internal
cairo_device <- function(file) {

    ## capabilities() reports whether the Cairo graphics library was compiled
    ## in, which is not guaranteed on all platforms. Reporting that here is
    ## clearer than the device failing when it is opened.
    if (!isTRUE(unname(capabilities("cairo")))) {
        stop("device = \"cairo\" requires an R installation built with Cairo ",
             "support, which capabilities(\"cairo\") reports as unavailable.",
             call. = FALSE)
    }

    ext <- tolower(tools::file_ext(file))

    device <- switch(
        ext,
        pdf = grDevices::cairo_pdf,
        ps  = ,
        eps = grDevices::cairo_ps,
        svg = grDevices::svg,
        NULL
    )

    if (is.null(device)) {
        stop("device = \"cairo\" applies to vector formats ('.pdf', '.ps', ",
             "'.eps', '.svg'), not '", ext, "'. Raster output already uses the ",
             "ragg devices where available.", call. = FALSE)
    }

    return(device)
}


### * Font embedding

#' Embed fonts in a PDF after it is written
#'
#' Post-processes a PDF so that its fonts are embedded, through
#' \code{grDevices::embedFonts()}. This is offered as an alternative to drawing
#' with a Cairo device, which embeds fonts but selects its own italic face.
#' Embedding is performed by Ghostscript (external).
#'
#' @param file Character string giving the path to the written file.
#' @param quiet Logical. Suppress the confirmation message.
#' @return Invisibly returns \code{TRUE} where the fonts were embedded and
#'   \code{FALSE} otherwise.
#' @keywords internal
embed_plot_fonts <- function(file, quiet = FALSE) {

    if (!identical(tolower(tools::file_ext(file)), "pdf")) {
        warning("'embed_fonts' applies to PDF output only; '",
                tools::file_ext(file), "' left unchanged.", call. = FALSE)
        return(invisible(FALSE))
    }

    ## embedFonts() calls Ghostscript, which ships separately from R.
    ## find_gs_cmd() consults R_GSCMD and the PATH, and returns an empty string
    ## when nothing is found.
    if (!nzchar(tools::find_gs_cmd())) {
        warning("Font embedding requires Ghostscript, which was not found. ",
                "The file was written without embedded fonts.", call. = FALSE)
        return(invisible(FALSE))
    }

    ## The plot is already on disk, so a failure here leaves a usable file
    ok <- tryCatch({
        grDevices::embedFonts(file)
        TRUE
    }, error = function(e) {
        warning("Font embedding failed: ", conditionMessage(e),
                ". The file was written without embedded fonts.",
                call. = FALSE)
        FALSE
    })

    if (ok && !quiet) {
        message("Fonts embedded in ", file)
    }

    return(invisible(ok))
}
