#' The yearly mean of a covariate over the fit data, weighted by the core share
#'
#' The z basis of the year climate term: `fireSenseUtils::yearClimateMeans()` of the raw (unscaled) ignition
#' covariates, each pixel weighted by `exposureOf(raw, "ignition")`. (`fireSense_ignitionPredict` uses the plain mean
#' over the ELF's pixels.)
#'
#' @param raw The raw ignition covariates, `sim$fireSense_ignitionCovariates`.
#' @param variable Covariate name, with or without underscores (`fireSenseUtils::yearClimateColumn()`).
#' @return A `data.frame` with `year` (calendar year) and `mean`.
yearClimateMeansOf <- function(raw, variable) {
  fireSenseUtils::yearClimateMeans(raw[[fireSenseUtils::yearClimateColumn(variable, names(raw))]],
                                   raw[[fireSenseUtils::yearTxt]], exposureOf(raw, "ignition"))
}

#' The year climate term of an ignition fit
#'
#' Fits `E[O_y] = P_y * exp(beta * z_y)`, `beta >= 0`, where `P_y` and `O_y` are the out-of-year predicted and
#' observed ignitions of year `y` (`yearTotals`) and `z_y` is the standardised yearly mean of `variable`
#' (`yearClimateMeansOf()`: centre is the mean over the fit years of the yearly means, scale their sd).
#' `beta` is the coefficient of `glm(observed ~ 0 + z, offset = log(predicted), family = poisson)`, floored at 0;
#' it is 0 if the glm fails or nothing is expected.
#'
#' @param variable Covariate name; see `yearClimateMeansOf()`. The column must be in `raw`.
#' @param raw The raw ignition covariates.
#' @param yearTotals The `data.frame` of `outOfFoldYearTotals()`: `year` (calendar year; the fit's is standardised
#'   when `rescaleVars`, so map it through the covariates' rows first), `observed`, `predicted`; or `NULL`.
#' @return `list(variable, center, scale, beta)`, `variable` being the covariate's column name.
fitYearClimate <- function(variable, raw, yearTotals) {
  variable <- fireSenseUtils::yearClimateColumn(variable, names(raw))
  m <- yearClimateMeansOf(raw, variable)
  center <- mean(m$mean)
  scale <- stats::sd(m$mean)
  if (!is.finite(scale) || scale <= 0) scale <- 1   # one year, or none differs: z is 0 and beta unidentified
  yearClimate <- list(variable = variable, center = center, scale = scale, beta = 0)
  if (is.null(yearTotals)) return(yearClimate)
  d <- data.frame(observed = yearTotals$observed, predicted = yearTotals$predicted,
                  z = (m$mean[match(yearTotals$year, m$year)] - center) / scale)
  d <- d[stats::complete.cases(d) & d$predicted > 0, ]
  if (nrow(d) == 0L || sum(d$observed) <= 0) return(yearClimate)
  b <- tryCatch(suppressWarnings(unname(stats::coef(
    stats::glm(observed ~ 0 + z, offset = log(predicted), family = stats::poisson, data = d))["z"])),
    error = function(e) 0)
  yearClimate$beta <- if (is.finite(b)) max(0, b) else 0
  yearClimate
}

#' The year climate multiplier of each row
#'
#' @param yearClimate `modelList$yearClimate` of the ignition fit, or `NULL` (no term).
#' @param raw The raw ignition covariates.
#' @param year The calendar year of each element of the result; defaults to the year of each row of `raw`.
#' @return Numeric, `fireSenseUtils::yearClimateMultiplier()` of each year's `yearClimateMeansOf()`; 1 if `yearClimate` is `NULL`.
yearClimateRowMultiplier <- function(yearClimate, raw, year = raw[[fireSenseUtils::yearTxt]]) {
  if (is.null(yearClimate)) return(rep(1, length(year)))
  m <- yearClimateMeansOf(raw, yearClimate$variable)
  fireSenseUtils::yearClimateMultiplier(yearClimate, m$mean[match(year, m$year)])
}

#' The calendar year of a year as fitted
#'
#' With `rescaleVars` the fit's `year` column is standardised, and so is the `year` of `outOfFoldYearTotals()`
#' (whose years went through `rownames()`, so they are right to about 15 digits, not exactly). The covariates keep
#' the order of the raw table, so a year as fitted is found among the covariates' and read from the raw table;
#' a year with no rows is not assumed to exist.
#'
#' @param year Years as fitted, e.g. `yearTotals$year`.
#' @param fitYear The `year` column of the covariates as fitted.
#' @param rawYear The `year` column of the raw covariates, same rows as `fitYear`.
#' @return The calendar year of each `year` (`NA` if it is not among `fitYear`).
calendarYearOf <- function(year, fitYear, rawYear) {
  rawYear[match(signif(year, 10), signif(fitYear, 10))]
}
