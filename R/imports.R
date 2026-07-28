#' Package Imports
#'
#' Declares the imports that are called without a namespace prefix, so that
#' R CMD check does not report them as undefined globals. Functions called as
#' \code{package::function()} elsewhere in the package need no entry here, but
#' their package must still appear under Imports in DESCRIPTION.
#'
#' @name summata-imports
#' @import data.table
#' @importFrom stats sd var median quantile IQR
#' @importFrom stats logLik AIC BIC coef confint vcov
#' @importFrom stats chisq.test fisher.test t.test wilcox.test kruskal.test
#' @importFrom stats lm glm anova
#' @importFrom stats family binomial gaussian poisson
#' @importFrom stats predict fitted residuals
#' @importFrom stats formula terms model.frame model.matrix
#' @importFrom stats na.omit complete.cases
#' @importFrom stats cor pnorm qnorm
#' @importFrom stats nobs pchisq pf
#' @importFrom grDevices axisTicks
#' @importFrom utils head tail str capture.output
#' @importFrom survival Surv strata cluster coxph
#' @keywords internal
NULL

#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
