#' Recommended Figure Dimensions
#'
#' Returns the figure dimensions recommended for a forest plot, computed from
#' the plot's own content by the function that produced it.
#'
#' @param x A forest plot produced by \code{lmforest()}, \code{glmforest()},
#'   \code{coxforest()}, \code{uniforest()}, \code{multiforest()}, or
#'   \code{autoforest()}.
#' @param units Character string giving the units the dimensions are returned
#'   in: \code{"in"} (inches), \code{"cm"}, or \code{"mm"}. Defaults to the
#'   units the plot was created in, so that the dimensions pass through
#'   unconverted unless a change is requested.
#'
#' @return A named numeric vector with elements \code{width} and \code{height},
#'   expressed in \code{units}. The units are recorded on the result as a
#'   \code{"units"} attribute, so that a value carried between functions
#'   remains self-describing.
#'
#' @details
#' Forest plots are sized by their content rather than by the page: a plot with
#' more variables, longer labels, or more factor levels requires a larger
#' canvas. Each plotting function computes a suitable size and attaches it to
#' the object it returns, and \code{recdims()} reports it.
#'
#' The values returned are those \code{forestsave()} applies, and are reported
#' exactly as computed rather than rounded, so that a figure written with
#' \code{forestsave()} and one written with \code{ggplot2::ggsave()} at these
#' dimensions are identical. Explicit use is therefore needed only where the
#' dimensions themselves are wanted, as when sizing a \pkg{knitr} chunk or
#' arranging several plots on one page.
#'
#' Where \code{units} differs from the units the plot was created in, the
#' dimensions are converted. A plot produced before the units were recorded
#' carries none, and is treated as having been created in inches.
#'
#' @seealso \code{\link{forestsave}} for saving at these dimensions,
#'   \code{\link{autoforest}} for producing the plot
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
#' # Example 1: Dimensions in the units the plot was created in
#' recdims(p)
#'
#' # Example 2: Journals commonly specify figure widths in millimeters
#' recdims(p, units = "mm")
#'
#' # Example 3: Sizing an output canvas from the individual elements
#' dims <- recdims(p)
#' dims[["width"]]
#' dims[["height"]]
#' attr(dims, "units")
#'
#' @family visualization functions
#' @export
recdims <- function(x, units = NULL) {

    if (!inherits(x, "ggplot")) {
        stop("'x' must be a forest plot produced by one of the forest plot ",
             "functions.", call. = FALSE)
    }

    dims <- attr(x, "rec_dims")

    if (is.null(dims) || is.null(dims$width) || is.null(dims$height)) {
        stop("'x' carries no recommended dimensions. These are attached by ",
             "the forest plot functions, so a plot from another source has ",
             "none.", call. = FALSE)
    }

    valid_units <- c("in", "cm", "mm", "px", "pt")

    ## Plots produced before the units were recorded carry none, and the
    ## dimensions were expressed in inches at that time
    from_units <- if (!is.null(dims$units) && dims$units %in% valid_units) {
                      dims$units
                  } else {
                      "in"
                  }

    ## Output units default to those of the recommendation, so the dimensions
    ## pass through unconverted in the common case
    if (is.null(units)) {
        units <- from_units
    }

    if (!is.character(units) || length(units) != 1 ||
        !units %in% valid_units) {
        stop("'units' must be one of: ",
             paste(sQuote(valid_units), collapse = ", "), ".",
             call. = FALSE)
    }

    result <- c(
        width  = unname(convert_units(dims$width,
                                      from = from_units, to = units)),
        height = unname(convert_units(dims$height,
                                      from = from_units, to = units))
    )

    ## Recorded so that a value carried into a saving function cannot be
    ## mis-sized by the conversion factor
    attr(result, "units") <- units

    return(result)
}
