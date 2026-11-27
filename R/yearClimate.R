#' The year climate term of an ignition or an escape fit
#'
#' Fits `E[O_y] = P_y * exp(beta * z_y) / norm`, `beta >= 0`, where `P_y` and `O_y` are the out-of-year predicted and
#' observed ignitions (or escapes) of year `y` (`yearTotals`) and `z_y` is the standardised yearly mean of `variable`
#' (`fireSenseUtils::yearClimateMeans()` of the raw covariates, the unweighted mean over the rows with the variable and with core (`coreFraction` above 0): centre is the mean over the fit years of the yearly
#' means, scale their sd). `beta` is the coefficient of
#' `glm(observed ~ 0 + z, offset = log(scaleToObservedLevel(observed, predicted)), family = poisson)`, floored at 0:
#' the offset removes the overall level gap, as `fireSenseUtils::yearEffectSD()` does, so `beta` is only the
#' year-to-year response. `norm` is the mean of `exp(beta * z_y)` over the fit years, so the multiplier averages 1
#' over them. `beta` is 0 (with a message) if there are fewer than 2 years with a mean, the means do not vary, nothing
#' is observed or expected, or the glm fails; `center` and `scale` are always finite.
#'
#' With `shape = "capped"` the term is `exp(beta * min(z, zCap))` divided by `norm` (`fireSenseUtils::cappedYearCurve()`):
#' log-linear up to `zCap` and flat above it, so it saturates in years drier than `zCap`. `beta` and the level are
#' fitted by the weighted Poisson glm `observed ~ min(z, zCap)` with offset `log(predicted)` for each `zCap`, and `zCap` is
#' profiled over the fit years' `z` (`fitCappedYearClimate()`). Nothing is bounded: a cap below the driest fit year is
#' kept only if it lowers the deviance by more than 2 (one parameter, as AIC counts it); otherwise `zCap` is the
#' driest fit year's `z`, which is the log-linear curve over every fit year, held flat beyond them. `beta >= 0`, as
#' for the log-linear shape. The stored list is `list(variable, center, scale, shape, beta, zCap, norm)`.
#' `shape = "logistic"` is retired for new fits (its steepness and upper level were not identified: with the midpoint
#' fixed at `z = 0` they ran to their bounds); a stored logistic fit is still reused and predicts as before.
#' With `weights = "climate"` each year counts `fireSenseUtils::annualYearWeights()` in the fit. `"loglinear"` and `"equal"` is the
#' term as before, to the last bit.
#'
#' @param variable Covariate name, with or without underscores. The column must be in `raw`.
#' @param raw The raw ignition covariates, `sim$fireSense_ignitionCovariates`.
#' @param yearTotals The `data.frame` of `outOfFoldYearTotals()`: `year` (calendar year; the fit's is standardised
#'   when `rescaleVars`, so map it with `calendarYearOf()` first), `observed`, `predicted`; or `NULL`.
#' @param shape `"loglinear"` (`exp(beta * z)`, the default) or `"capped"`.
#' @param weights `"equal"` (every year counts the same) or `"climate"` (`fireSenseUtils::annualYearWeights()`).
#' @return `list(variable, center, scale, beta, norm)`, or with `shape = "capped"`
#'   `list(variable, center, scale, shape, beta, zCap, norm)`; then `weights` when not `"equal"`. `variable` is the
#'   covariate's column name.
fitYearClimate <- function(variable, raw, yearTotals, shape = "loglinear", weights = "equal") {
  weights <- match.arg(weights, c("equal", "climate"))
  if (identical(shape, "logistic"))
    stop("fireSense_ignitionFit: `yearClimateShape = \"logistic\"` can no longer be fitted: its steepness and upper level ",
         "ran to their bounds (not identified). Use \"capped\" (log-linear up to a fitted cap). Stored logistic fits are ",
         "still reused.", call. = FALSE)
  shape <- match.arg(shape, c("loglinear", "capped"))
  capped <- shape == "capped"
  variable <- fireSenseUtils::yearClimateColumn(variable, names(raw))
  if (is.null(raw[[fireSenseUtils::coreFractionTxt]]))
    warning("fireSense_ignitionFit: the ignition covariates have no `", fireSenseUtils::coreFractionTxt, "` column, so the ",
            "year climate mean is over every pixel, while fireSense_ignitionPredict averages the core pixels only ",
            "(fireSense_dataPrepFit adds the column).", call. = FALSE)
  m <- fireSenseUtils::yearClimateMeans(raw, variable)
  scale <- stats::sd(m$mean)
  noTerm <- function(why) {
    message("fireSense_ignitionFit: no year climate term (", variable, "): ", why, ".")
    yc <- list(variable = variable, center = if (nrow(m)) mean(m$mean) else 0,
               scale = if (is.finite(scale) && scale > 0) scale else 1)
    yc <- if (capped) c(yc, list(shape = shape, beta = 0, zCap = 0, norm = 1)) else c(yc, list(beta = 0, norm = 1))
    if (weights != "equal") yc$weights <- weights
    yc
  }
  if (nrow(m) < 2L || !is.finite(scale) || scale <= 0)
    return(noTerm(paste0(nrow(m), " fit year(s) with a mean, and the yearly means do not vary")))
  center <- mean(m$mean)
  if (is.null(yearTotals)) return(noTerm("no out-of-year totals"))
  d <- data.frame(observed = yearTotals$observed, predicted = yearTotals$predicted,
                  z = fireSenseUtils::standardiseYearClimate(m$mean[match(yearTotals$year, m$year)], center, scale))
  d <- d[stats::complete.cases(d) & d$predicted > 0, ]
  if (nrow(d) == 0L || sum(d$observed) <= 0) return(noTerm("nothing observed or expected"))
  d$w <- fireSenseUtils::annualYearWeights(d$z, weights)
  yc <- list(variable = variable, center = center, scale = scale)
  if (capped) {
    fit <- fitCappedYearClimate(d)
    if (is.null(fit)) return(noTerm("the glm failed"))
    message("fireSense_ignitionFit: year climate term (", variable, "), capped: beta = ", signif(fit$beta, 3), ", zCap = ",
            signif(fit$zCap, 3), " (driest fit year ", signif(max(d$z), 3), "; deviance gain of the cap ",
            signif(fit$gain, 3), ", kept if > 2) from ", nrow(d), " of ", nrow(m), " fit years.")
    yc <- c(yc, list(shape = shape, beta = fit$beta, zCap = fit$zCap))
    keep <- c("variable", "center", "scale", "shape", "beta", "zCap", "norm")
  } else {
    d$offset <- log(fireSenseUtils::scaleToObservedLevel(d$observed, d$predicted))
    f <- tryCatch(suppressWarnings(stats::glm(observed ~ 0 + z, offset = offset, weights = d$w, family = stats::poisson,
                                              data = d)), error = function(e) NULL)
    b <- if (is.null(f)) NA_real_ else unname(stats::coef(f)["z"])
    if (!is.finite(b)) return(noTerm("the glm failed"))
    beta <- max(0, b)
    message("fireSense_ignitionFit: year climate term (", variable, "): beta = ", signif(beta, 3), " from ", nrow(d),
            " of ", nrow(m), " fit years.")
    yc <- c(yc, list(beta = beta))
    keep <- c("variable", "center", "scale", "beta", "norm")
  }
  yc$norm <- mean(fireSenseUtils::yearClimateMultiplier(yc, m$mean))
  if (weights != "equal") yc$weights <- weights   # recorded so a stored fit can be compared with a run's settings
  yc[c(keep, if (weights != "equal") "weights")]
}

#' The bounds of the logistic year climate term
#'
#' The bounds stored logistic fits were made under (`kMax`, `ratioMax`, `capFactor`; the logistic shape is no longer
#' fitted). A stored logistic term is reused only if it records these (`storedYearClimateProblem()`).
logisticYearBounds <- c(capFactor = 3, kMax = 10, ratioMax = 1000)

#' The capped year climate term: `beta` and `zCap`
#'
#' For a cap `zCap`, the weighted Poisson glm `observed ~ x`, `x = min(z, zCap)`, offset `log(predicted)`, gives `beta`
#' (floored at 0, then the level alone is fitted) and a deviance. `zCap` is profiled: `optimize()` within each interval
#' between the sorted distinct `z` from the second smallest up, and at each `z`. The log-linear curve is the cap at the
#' largest `z`; a lower cap is kept only if its deviance is more than 2 below that.
#'
#' @param d `data.frame` with `observed`, `predicted`, `z` and `w` of the years of the fit.
#' @return `list(beta, zCap, gain)` (`gain`: the deviance of the log-linear curve minus the chosen one), or `NULL` if
#'   the glm fails.
fitCappedYearClimate <- function(d) {
  fitAt <- function(zCap) {
    x <- pmin(d$z, zCap)
    f <- tryCatch(suppressWarnings(stats::glm(d$observed ~ x, offset = log(d$predicted), weights = d$w,
                                              family = stats::poisson)), error = function(e) NULL)
    if (is.null(f) || !all(is.finite(stats::coef(f)))) return(list(beta = NA_real_, deviance = Inf))
    if (stats::coef(f)[["x"]] > 0) return(list(beta = stats::coef(f)[["x"]], deviance = stats::deviance(f)))
    f0 <- suppressWarnings(stats::glm(d$observed ~ 1, offset = log(d$predicted), weights = d$w, family = stats::poisson))
    list(beta = 0, deviance = stats::deviance(f0))
  }
  zs <- sort(unique(d$z))
  top <- zs[length(zs)]
  linear <- fitAt(top)
  if (!is.finite(linear$deviance)) return(NULL)
  best <- list(zCap = top, beta = linear$beta, deviance = linear$deviance)
  if (length(zs) >= 3L) {
    cand <- zs[-1]
    for (i in seq_len(length(cand) - 1L)) {
      o <- stats::optimize(function(m) fitAt(m)$deviance, c(cand[i], cand[i + 1L]))
      cand <- c(cand, o$minimum)
    }
    for (m in cand) {
      f <- fitAt(m)
      if (f$deviance < best$deviance) best <- list(zCap = m, beta = f$beta, deviance = f$deviance)
    }
  }
  gain <- linear$deviance - best$deviance
  if (gain <= 2) best <- list(zCap = top, beta = linear$beta)
  if (best$beta == 0) best$zCap <- top
  list(beta = best$beta, zCap = best$zCap, gain = gain)
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

#' Does a stored fit's year climate term differ from this run's settings?
#'
#' A fit stored in the ledger carries the term it was made with (`modelList$yearClimate`; `shape` with
#' `yearClimateShape = "logistic"` (with its `bounds`), `weights` with `yearClimateWeights = "climate"`). The term must be
#' there exactly when `ignitionYearClimate` (`escapeYearClimate`) is set, of the same variable (with or without
#' underscores). Where both the stored fit and this run have a term, `yearClimateShape`, the bounds of a logistic one
#' (`logisticYearBounds`) and `yearClimateWeights` must agree. A fit stored before these parameters existed has
#' neither `shape`, `bounds`, `weights` nor an escape term: it matches only their defaults, and one made before the
#' ignition term existed is refitted when `ignitionYearClimate` is set.
#'
#' @param fitted A stored `list(modelList = ...)`, or `NULL`.
#' @param igOrEsc `"ignition"` or `"escape"`.
#' @param par The module's parameters (`P(sim)`): `ignitionYearClimate`, `escapeYearClimate`, `yearClimateShape`,
#'   `yearClimateWeights`.
#' @return `NULL` if they agree, else a sentence saying how they differ.
storedYearClimateProblem <- function(fitted, igOrEsc, par) {
  yc <- fitted$modelList$yearClimate
  variable <- par[[paste0(igOrEsc, "YearClimate")]]
  want <- length(variable) == 1L && !is.na(variable)
  param <- paste0(igOrEsc, "YearClimate")
  squash <- function(x) gsub("_", "", x, fixed = TRUE)
  if (want != !is.null(yc))
    return(paste0("it has ", if (is.null(yc)) "no" else "an", " year climate term, but `", param, "` is ",
                  if (want) paste0("\"", variable, "\"") else "NA"))
  if (!want) return(NULL)
  if (!identical(squash(yc$variable), squash(variable)))
    return(paste0("its year climate term is of ", yc$variable, ", but `", param, "` is \"", variable, "\""))
  storedShape <- if (is.null(yc$shape)) "loglinear" else yc$shape
  now <- if (is.null(par$yearClimateShape)) "loglinear" else par$yearClimateShape
  if (!identical(storedShape, now))
    return(paste0("it was fitted with `yearClimateShape = \"", storedShape, "\"`, but this run's is \"", now, "\""))
  if (identical(storedShape, "logistic") && !isTRUE(all.equal(yc$bounds, logisticYearBounds)))
    return(paste0("its logistic year climate term was fitted under other bounds (cap, k, ratio) than this version's"))
  stored <- if (is.null(yc$weights)) "equal" else yc$weights
  nowWeights <- if (is.null(par$yearClimateWeights)) "equal" else par$yearClimateWeights
  if (!identical(stored, nowWeights))
    return(paste0("its year climate term was fitted with `yearClimateWeights = \"", stored, "\"`, but this run's is \"", nowWeights, "\""))
  NULL
}
