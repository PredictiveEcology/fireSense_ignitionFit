#' How concentrated the predicted events are, against the observed ones
#'
#' A boosted model with few, shallow trees and shrinkage pulls predictions toward the mean, so it can spread the
#' expected events more evenly over the pixel-years than they happen. This measures it on out-of-fold (or
#' out-of-year) predictions, so a model is not credited for what it has seen. Pixel-years are ranked by the
#' prediction, highest first; tied predictions share their events equally.
#'
#' @param observed The observed count of each pixel-year (in the core, for ignitions).
#' @param expected The out-of-fold prediction of each pixel-year, as a count in the same units
#'   (`outOfFoldPrediction()`).
#' @param basis `"out-of-year"` or `"out-of-fold"`: what the predictions are (`runXGBOOST()`).
#' @param top The shares of the pixel-years of the lift table (the top 1, 5, 10 and 25%).
#' @return A list: `basis`, `n` (pixel-years), `observed` and `expected` (the totals), `lift` (a `data.frame`
#'   with `top`, `observed`, the share of the observed events in the top `top` of the pixel-years, and `predicted`,
#'   the share of the predicted events there: `predicted` below `observed` means the model spreads them too
#'   evenly), `gini` (the Gini index of the observed and of the predicted counts across pixel-years: 0 is
#'   even, 1 is all in one; the model's below the observed means too even), and `deviance`, the mean Poisson
#'   deviance per pixel-year, `2 * (y * log(y / mu) - (y - mu))`.
concentrationOf <- function(observed, expected, basis = "out-of-fold", top = c(0.01, 0.05, 0.10, 0.25)) {
  ok <- !is.na(expected) & !is.na(observed)
  observed <- as.numeric(observed[ok]); expected <- expected[ok]
  share <- function(x) { # cumulative share of `x` at each `top` of the rows ranked by `expected`; linear within ties
    o <- order(expected, decreasing = TRUE)
    tie <- cumsum(c(TRUE, diff(expected[o]) != 0))
    rows <- c(0, cumsum(tabulate(tie)) / length(x))
    events <- c(0, cumsum(rowsum(x[o], tie)[, 1]))
    stats::approx(rows, events / sum(x), xout = top)$y
  }
  mu <- pmax(expected, 1e-10)
  dev <- 2 * (ifelse(observed > 0, observed * log(pmax(observed, 1e-300) / mu), 0) - (observed - mu))
  list(basis = basis, n = length(observed), observed = sum(observed), expected = sum(expected),
       lift = data.frame(top = top, observed = share(observed), predicted = share(expected)),
       gini = c(observed = giniIndex(observed), predicted = giniIndex(expected)),
       deviance = mean(dev))
}

## The Gini index of non-negative `x`: 0 when all equal, 1 when one has everything; NA when the total is 0.
giniIndex <- function(x) {
  n <- length(x)
  if (!n || sum(x) <= 0) return(NA_real_)
  x <- sort(x)
  2 * sum(seq_len(n) * x) / (n * sum(x)) - (n + 1) / n
}

#' The rows of the ignition table that each row of the escape table is the pixel-year of
#'
#' @inheritParams chainedEscapeFactor
#' @return Integer, one per row of `escapes`; stops when a row has no ignition row of the same pixel and year.
matchedIgnitionRows <- function(ignitions, escapes) {
  rows <- match(paste(escapes$pixelID, escapes$year), paste(ignitions$pixelID, ignitions$year))
  if (anyNA(rows))
    stop("fittedYearEffects(): ", sum(is.na(rows)), " escape rows have no ignition row of the same pixel and year")
  rows
}

#' The factor by which the escapes that follow model-placed ignitions differ from the observed ones
#'
#' When `fireSense_burn` places ignitions with the ignition model, each ignition escapes with the escape model's
#' probability at the pixel it lands in (`fireSenseUtils::escapeProbPerIgnition()`), and the expected escapes of
#' the fit period are `sum(w * lambda * q)` over the pixel-years (`w` the share in the core, `exposureOf()`,
#' `lambda` the expected ignitions of the whole pixel, `q` the escape probability per ignition), the quantities
#' `fittedYearEffects()` expects the years' escapes with, in-sample (fold-mean predictions). `c` is the
#' observed escapes over these. It is split exactly into three factors,
#' `c = ignitionRatio * escapeCalibration * placement`: the observed over the expected ignitions; the observed
#' escapes over those expected from the observed ignitions (`sum(ignitions * q)`: the escape model's calibration
#' where ignitions happened); and `qObserved / qModel`, the mean `q` where ignitions were observed over the mean
#' `q` where the model places them. A placement factor above 1 means the ignition model places ignitions where
#' escaping is less likely than where they happened: it is spread too evenly.
#'
#' @param ignitions `data.table` of the ignition fit data, one row per pixel and year, as for
#'   `fittedYearEffects()`: `pixelID`, `year`, `observed` (ignitions), `mean` (predicted ignitions of the whole
#'   pixel) and `weight` (its share in the core).
#' @param escapes `data.table` of the escape fit data for the same pixel-years: `pixelID`, `year`, `observed`
#'   (escapes) and `mean` (the predicted escapes if ignited).
#' @return A list: `c`, `ignitionRatio`, `escapeCalibration`, `placement`, `qObserved`, `qModel`, and the totals
#'   `observedIgnitions`, `expectedIgnitions`, `observedEscapes`, `expectedEscapes` (`c = observedEscapes / expectedEscapes`).
chainedEscapeFactor <- function(ignitions, escapes) {
  rows <- matchedIgnitionRows(ignitions, escapes)
  lambda <- ignitions$mean[rows]
  w <- ignitions$weight[rows]
  q <- fireSenseUtils::escapeProbPerIgnition(escapes$mean, lambda)
  obsIgn <- ignitions$observed[rows]
  expIgn <- sum(w * lambda)
  expEsc <- sum(w * lambda * q)
  qObserved <- sum(obsIgn * q) / sum(obsIgn)
  qModel <- expEsc / expIgn
  list(c = sum(escapes$observed) / expEsc, ignitionRatio = sum(obsIgn) / expIgn,
       escapeCalibration = sum(escapes$observed) / (sum(obsIgn) * qObserved), placement = qObserved / qModel,
       qObserved = qObserved, qModel = qModel,
       observedIgnitions = sum(obsIgn), expectedIgnitions = expIgn,
       observedEscapes = sum(escapes$observed), expectedEscapes = expEsc)
}

#' One line of the log for a fit's concentration
#'
#' @param concentration `concentrationOf()` of the fit.
#' @param process `"ignition"` or `"escape"`.
#' @return Character, one string.
concentrationMessage <- function(concentration, process) {
  lift <- concentration$lift
  pct <- function(x) paste0(format(round(100 * x, 1)), collapse = "/")
  paste0("fireSense_ignitionFit ", process, " concentration (", concentration$basis, ", ", concentration$n,
         " pixel-years): share of events in the top ", paste0(format(100 * lift$top), collapse = "/"), "%, observed ",
         pct(lift$observed), "% vs predicted ", pct(lift$predicted), "%; Gini observed ",
         format(round(concentration$gini[["observed"]], 3)), " vs predicted ",
         format(round(concentration$gini[["predicted"]], 3)), "; deviance ", format(round(concentration$deviance, 4)))
}

#' One line of the log for the chained-escape factor
#'
#' @param chain `chainedEscapeFactor()`.
#' @return Character, one string.
chainMessage <- function(chain) {
  paste0("fireSense_ignitionFit chained-escape factor c = ", format(round(chain$c, 3)), " = ignitions ",
         format(round(chain$ignitionRatio, 3)), " x escape model ", format(round(chain$escapeCalibration, 3)),
         " x placement ", format(round(chain$placement, 3)), " (observed ", round(chain$observedEscapes, 1),
         " escapes, expected ", round(chain$expectedEscapes, 1), " when the ignition model places the ignitions)")
}
