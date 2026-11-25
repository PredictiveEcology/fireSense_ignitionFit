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
#' With `shape = "logistic"` the term is an S-curve through its midpoint at `z = 0`, the ELF's average year:
#' `g(z) = 1 + (ratio - 1) * plogis(k * z)` divided by `norm`, the mean of `g` over the fit years (so it averages 1,
#' as above). Two parameters are fitted, `k >= 0` (the steepness, in 1 / sd) and `ratio >= 1` (the upper level over the
#' floor); the floor is 1 before `norm`, so the multiplier is always positive, and it saturates beyond the fit range
#' at `ratio / norm` (the dry side, `z > 0`) and `1 / norm` (the wet side). They maximise the
#' (weighted) Poisson likelihood of `O_y ~ c * P_y * g(z_y)`, with the level `c` profiled out (`optim()`, Nelder-Mead,
#' from several starts), and `ratio / norm`, the upper level over the average, may not exceed 3 times the largest
#' observed year total over the mean of the observed year totals, `max(O_y) / mean(O_y)` (a fit with no gain over a
#' flat multiplier stores `k = 0`, `ratio = 1`). The stored list is `list(variable, center, scale, shape, k, ratio, norm)`.
#' With `weights = "climate"` each year counts `yearClimateWeights()` in the fit. `"loglinear"` and `"equal"` is the
#' term as before, to the last bit.
#'
#' @param variable Covariate name, with or without underscores. The column must be in `raw`.
#' @param raw The raw ignition covariates, `sim$fireSense_ignitionCovariates`.
#' @param yearTotals The `data.frame` of `outOfFoldYearTotals()`: `year` (calendar year; the fit's is standardised
#'   when `rescaleVars`, so map it with `calendarYearOf()` first), `observed`, `predicted`; or `NULL`.
#' @param shape `"loglinear"` (`exp(beta * z)`, the default) or `"logistic"`.
#' @param weights `"equal"` (every year counts the same) or `"climate"` (`yearClimateWeights()`).
#' @return `list(variable, center, scale, beta, norm)`, or with `shape = "logistic"`
#'   `list(variable, center, scale, shape, k, ratio, norm)`; then `weights` when not `"equal"`. `variable` is the
#'   covariate's column name.
fitYearClimate <- function(variable, raw, yearTotals, shape = "loglinear", weights = "equal") {
  weights <- match.arg(weights, c("equal", "climate"))
  shape <- match.arg(shape, c("loglinear", "logistic"))
  logistic <- shape == "logistic"
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
    yc <- if (logistic) c(yc, list(shape = shape, k = 0, ratio = 1, norm = 1)) else c(yc, list(beta = 0, norm = 1))
    if (weights != "equal") yc$weights <- weights
    yc
  }
  if (nrow(m) < 2L || !is.finite(scale) || scale <= 0)
    return(noTerm(paste0(nrow(m), " fit year(s) with a mean, and the yearly means do not vary")))
  center <- mean(m$mean)
  if (is.null(yearTotals)) return(noTerm("no out-of-year totals"))
  d <- data.frame(observed = yearTotals$observed, predicted = yearTotals$predicted,
                  z = (m$mean[match(yearTotals$year, m$year)] - center) / scale)
  d <- d[stats::complete.cases(d) & d$predicted > 0, ]
  if (nrow(d) == 0L || sum(d$observed) <= 0) return(noTerm("nothing observed or expected"))
  d$w <- yearClimateWeights(d$z, weights)
  yc <- list(variable = variable, center = center, scale = scale)
  if (logistic) {
    fit <- fitLogisticYearClimate(d, (m$mean - center) / scale)
    if (is.null(fit)) return(noTerm("the fit failed"))
    message("fireSense_ignitionFit: year climate term (", variable, "), logistic: k = ", signif(fit$k, 3), ", ratio = ",
            signif(fit$ratio, 3), " from ", nrow(d), " of ", nrow(m), " fit years.")
    yc <- c(yc, list(shape = shape, k = fit$k, ratio = fit$ratio))
    keep <- c("variable", "center", "scale", "shape", "k", "ratio", "norm")
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

#' The logistic year climate term: `k` and `ratio`
#'
#' Maximises the (weighted) Poisson likelihood of `O_y ~ c * P_y * g(z_y)`, `g = 1 + (ratio - 1) * plogis(k * z)`, over
#' `k` in 0 to 10 and `ratio` in 1 to 1000, with `c` profiled out, from a grid of starting values (`optim()`, Nelder-Mead
#' on an unconstrained scale). The upper level over the average, `ratio / mean(g)` over the fit years, may not exceed 3
#' times `max(O_y) / mean(O_y)` (a point outside is not accepted). A fit no better than the flat multiplier gives
#' `k = 0`, `ratio = 1`.
#'
#' @param d `data.frame` with `observed`, `predicted`, `z` and `w` (`NULL`: equal) of the years of the fit.
#' @param zFit `z` of every fit year (the years `norm` averages over).
#' @return `list(k, ratio)`, or `NULL` if no start gave a finite likelihood.
fitLogisticYearClimate <- function(d, zFit) {
  w <- if (is.null(d$w)) rep(1, nrow(d)) else d$w
  kMax <- 10; ratioMax <- 1000
  cap <- 3 * max(d$observed) / mean(d$observed)
  dev <- function(k, ratio) {
    g <- 1 + (ratio - 1) * stats::plogis(k * d$z)
    mu <- sum(w * d$observed) / sum(w * d$predicted * g) * d$predicted * g
    2 * sum(w * (ifelse(d$observed > 0, d$observed * log(d$observed / mu), 0) - (d$observed - mu)))
  }
  objective <- function(u) {
    k <- kMax * stats::plogis(u[1])
    ratio <- 1 + min(exp(u[2]), ratioMax - 1)
    if (ratio / mean(1 + (ratio - 1) * stats::plogis(k * zFit)) > cap) return(1e12)
    out <- dev(k, ratio)
    if (is.finite(out)) out else 1e12
  }
  best <- list(value = dev(0, 1), k = 0, ratio = 1)
  if (!is.finite(best$value)) return(NULL)
  for (k0 in c(0.5, 1, 2, 4)) for (r0 in c(0.3, 1, 3)) {
    o <- tryCatch(stats::optim(c(stats::qlogis(k0 / kMax), log(r0)), objective, method = "Nelder-Mead"),
                  error = function(e) NULL)
    if (!is.null(o) && o$value < best$value - 1e-9)
      best <- list(value = o$value, k = kMax * stats::plogis(o$par[1]), ratio = 1 + min(exp(o$par[2]), ratioMax - 1))
  }
  list(k = best$k, ratio = best$ratio)
}

#' The weights of the years in the fit of the year climate term
#'
#' `"equal"`: none (`NULL`; every year counts the same). `"climate"`: `|z_y|` floored at 0.25, then rescaled to mean 1,
#' so the extreme years, wet or dry, decide the slope and a year near the mean cannot flatten it.
#'
#' @param z The standardised yearly climate mean of the years in the fit.
#' @param weights `"equal"` or `"climate"`.
#' @return `NULL`, or one weight per `z`, averaging 1.
yearClimateWeights <- function(z, weights = "equal") {
  weights <- match.arg(weights, c("equal", "climate"))
  if (weights == "equal") return(NULL)
  w <- pmax(abs(z), 0.25)
  w / mean(w)
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
#' `yearClimateShape = "logistic"`, `weights` with `yearClimateWeights = "climate"`). The escape term must be there exactly when
#' `escapeYearClimate` is set (same variable, with or without underscores). Where both the stored fit and this run have a
#' term, `yearClimateShape` and `yearClimateWeights` must agree. A fit stored before these parameters existed has
#' neither `shape`, `weights` nor an escape term: it matches only their defaults.
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
  squash <- function(x) gsub("_", "", x, fixed = TRUE)
  if (igOrEsc == "escape" && want != !is.null(yc))
    return(paste0("it has ", if (is.null(yc)) "no" else "an", " year climate term, but `escapeYearClimate` is ",
                  if (want) paste0("\"", variable, "\"") else "NA"))
  if (igOrEsc == "escape" && want && !identical(squash(yc$variable), squash(variable)))
    return(paste0("its year climate term is of ", yc$variable, ", but `escapeYearClimate` is \"", variable, "\""))
  if (!want || is.null(yc)) return(NULL)
  storedShape <- if (is.null(yc$shape)) "loglinear" else yc$shape
  now <- if (is.null(par$yearClimateShape)) "loglinear" else par$yearClimateShape
  if (!identical(storedShape, now))
    return(paste0("it was fitted with `yearClimateShape = \"", storedShape, "\"`, but this run's is \"", now, "\""))
  stored <- if (is.null(yc$weights)) "equal" else yc$weights
  nowWeights <- if (is.null(par$yearClimateWeights)) "equal" else par$yearClimateWeights
  if (!identical(stored, nowWeights))
    return(paste0("its year climate term was fitted with `yearClimateWeights = \"", stored, "\"`, but this run's is \"", nowWeights, "\""))
  NULL
}
