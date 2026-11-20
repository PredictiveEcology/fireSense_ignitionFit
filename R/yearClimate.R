#' The year climate term of an ignition fit
#'
#' Fits `E[O_y] = P_y * exp(beta * z_y) / norm`, `beta >= 0`, where `P_y` and `O_y` are the out-of-year predicted and
#' observed ignitions of year `y` (`yearTotals`) and `z_y` is the standardised yearly mean of `variable`
#' (`fireSenseUtils::yearClimateMeans()` of the raw covariates: centre is the mean over the fit years of the yearly
#' means, scale their sd). `beta` is the coefficient of
#' `glm(observed ~ 0 + z, offset = log(predicted * sum(observed) / sum(predicted)), family = poisson)`, floored at 0:
#' the offset removes the overall level gap, as `fireSenseUtils::yearEffectSD()` does, so `beta` is only the
#' year-to-year response. `norm` is the mean of `exp(beta * z_y)` over the fit years, so the multiplier averages 1
#' over them. `beta` is 0 (with a message) if there are fewer than 2 years with a mean, the means do not vary, nothing
#' is observed or expected, or the glm fails; `center` and `scale` are always finite.
#'
#' @param variable Covariate name, with or without underscores. The column must be in `raw`.
#' @param raw The raw ignition covariates, `sim$fireSense_ignitionCovariates`.
#' @param yearTotals The `data.frame` of `outOfFoldYearTotals()`: `year` (calendar year; the fit's is standardised
#'   when `rescaleVars`, so map it with `calendarYearOf()` first), `observed`, `predicted`; or `NULL`.
#' @return `list(variable, center, scale, beta, norm)`, `variable` being the covariate's column name.
fitYearClimate <- function(variable, raw, yearTotals) {
  variable <- fireSenseUtils::yearClimateColumn(variable, names(raw))
  m <- fireSenseUtils::yearClimateMeans(raw, variable)
  scale <- stats::sd(m$mean)
  noTerm <- function(why) {
    message("fireSense_ignitionFit: no year climate term (", variable, "): ", why, ".")
    list(variable = variable, center = if (nrow(m)) mean(m$mean) else 0,
         scale = if (is.finite(scale) && scale > 0) scale else 1, beta = 0, norm = 1)
  }
  if (nrow(m) < 2L || !is.finite(scale) || scale <= 0)
    return(noTerm(paste0(nrow(m), " fit year(s) with a mean, and the yearly means do not vary")))
  center <- mean(m$mean)
  if (is.null(yearTotals)) return(noTerm("no out-of-year totals"))
  d <- data.frame(observed = yearTotals$observed, predicted = yearTotals$predicted,
                  z = (m$mean[match(yearTotals$year, m$year)] - center) / scale)
  d <- d[stats::complete.cases(d) & d$predicted > 0, ]
  if (nrow(d) == 0L || sum(d$observed) <= 0) return(noTerm("nothing observed or expected"))
  b <- tryCatch(suppressWarnings(unname(stats::coef(
    stats::glm(observed ~ 0 + z, offset = log(predicted * sum(observed) / sum(predicted)), family = stats::poisson,
               data = d))["z"])), error = function(e) NA_real_)
  if (!is.finite(b)) return(noTerm("the glm failed"))
  beta <- max(0, b)
  message("fireSense_ignitionFit: year climate term (", variable, "): beta = ", signif(beta, 3), " from ", nrow(d),
          " of ", nrow(m), " fit years.")
  list(variable = variable, center = center, scale = scale, beta = beta, norm = mean(exp(beta * (m$mean - center) / scale)))
}

#' The year climate multiplier of each row
#'
#' @param yearClimate `modelList$yearClimate` of the ignition fit, or `NULL` (no term).
#' @param raw The raw ignition covariates.
#' @param year The calendar year of each element of the result; defaults to the year of each row of `raw`.
#' @return Numeric, `fireSenseUtils::yearClimateMultiplier()` of each year's `fireSenseUtils::yearClimateMeans()`;
#'   1 if `yearClimate` is `NULL`, or for a year with no usable rows.
yearClimateRowMultiplier <- function(yearClimate, raw, year = raw[[fireSenseUtils::yearTxt]]) {
  if (is.null(yearClimate)) return(rep(1, length(year)))
  m <- fireSenseUtils::yearClimateMeans(raw, yearClimate$variable)
  out <- fireSenseUtils::yearClimateMultiplier(yearClimate, m$mean[match(year, m$year)])
  out[is.na(out)] <- 1
  out
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
#' @return The calendar year of each `year`; stops if one is not among `fitYear`.
calendarYearOf <- function(year, fitYear, rawYear) {
  out <- rawYear[match(signif(year, 10), signif(fitYear, 10))]
  if (anyNA(out))
    stop("fireSense_ignitionFit: ", sum(is.na(out)), " of ", length(year), " years of the out-of-year totals are not ",
         "among the years of the covariates, so the year climate term cannot be fitted.", call. = FALSE)
  out
}
