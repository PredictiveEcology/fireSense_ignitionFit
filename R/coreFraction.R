#' The exposure of each row of the fit data
#'
#' The ignitions are counted in the core of the ELF, so a pixel with a core share of `coreFraction` can hold only
#' that share of the pixel's expected ignitions. For a Poisson deviance, `weight = coreFraction` with
#' `label = ignitions / coreFraction` is the offset `log(coreFraction)`: the model learns the rate of a whole
#' coarse pixel, and its prediction stays that of the whole pixel. Only the ignition fit is weighted; the
#' escape fit (rows with an ignition, escapes per ignition) is not, and without the column every weight is 1.
#'
#' @param dat data.table of the fit data, with `coreFraction` if fireSense_dataPrepFit made it.
#' @param type `"ignition"` or `"escape"`.
#' @return Numeric, one per row: `coreFraction` for `"ignition"` with the column, otherwise 1.
exposureOf <- function(dat, type) {
  w <- if (identical(type, "ignition")) dat[[fireSenseUtils::coreFractionTxt]]
  if (is.null(w)) rep(1, nrow(dat)) else w
}
