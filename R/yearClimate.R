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
#' With `hinge` the term is `exp(beta * z + beta2 * max(0, z))`, `beta, beta2 >= 0` (`beta2` is stored): the best of the
#' fits with both coefficients, `beta` only and `beta2` only that have no negative coefficient, the level being a free
#' intercept there (the hinge's multiplier does not average 1; `norm` re-centres it). With
#' `weights = "climate"` each year counts `yearClimateWeights()` in the glm. `"equal"` and no hinge is the term
#' as before, to the last bit.
#'
#' @param variable Covariate name, with or without underscores. The column must be in `raw`.
#' @param raw The raw ignition covariates, `sim$fireSense_ignitionCovariates`.
#' @param yearTotals The `data.frame` of `outOfFoldYearTotals()`: `year` (calendar year; the fit's is standardised
#'   when `rescaleVars`, so map it with `calendarYearOf()` first), `observed`, `predicted`; or `NULL`.
#' @param hinge `TRUE`: also fit `beta2`, the extra slope of the dry half (`z > 0`).
#' @param weights `"equal"` (every year counts the same) or `"climate"` (`yearClimateWeights()`).
#' @return `list(variable, center, scale, beta, norm)` (`beta2` after `norm` with `hinge`, then `weights` when not `"equal"`), `variable` being the
#'   covariate's column name.
fitYearClimate <- function(variable, raw, yearTotals, hinge = FALSE, weights = "equal") {
  weights <- match.arg(weights, c("equal", "climate"))
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
               scale = if (is.finite(scale) && scale > 0) scale else 1, beta = 0, norm = 1)
    if (hinge) yc$beta2 <- 0
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
  d$zHigh <- pmax(0, d$z)
  ## Without the hinge the overall level gap is removed by the offset (the term is the one made before). With it the
  ## level is a free intercept: the hinge's multiplier does not average 1 over the years, `norm` re-centres it.
  d$offset <- log(if (hinge) d$predicted else fireSenseUtils::scaleToObservedLevel(d$observed, d$predicted))
  fitTerms <- function(terms) tryCatch(suppressWarnings(stats::glm(
    stats::reformulate(c(if (hinge) "1" else "0", terms), "observed", intercept = TRUE), offset = offset, weights = d$w,
    family = stats::poisson, data = d)), error = function(e) NULL)
  if (hinge) {
    ## best fit with both coefficients >= 0: the unconstrained one if it has them, else one coefficient at 0
    fits <- lapply(list(c("z", "zHigh"), "z", "zHigh"), fitTerms)
    coefs <- lapply(fits, function(f) if (!is.null(f) && all(is.finite(stats::coef(f)))) stats::coef(f)[-1])
    ok <- vapply(coefs, function(cf) !is.null(cf) && all(cf >= 0), logical(1))
    if (all(vapply(fits, is.null, logical(1)))) return(noTerm("the glm failed"))
    b <- b2 <- 0
    if (any(ok)) {
      cf <- coefs[[which(ok)[which.min(vapply(fits[ok], stats::deviance, numeric(1)))]]]
      if ("z" %in% names(cf)) b <- unname(cf["z"])
      if ("zHigh" %in% names(cf)) b2 <- unname(cf["zHigh"])
    }
    beta <- max(0, b)
    beta2 <- max(0, b2)
  } else {
    f <- fitTerms("z")
    b <- if (is.null(f)) NA_real_ else unname(stats::coef(f)["z"])
    if (!is.finite(b)) return(noTerm("the glm failed"))
    beta <- max(0, b)
    beta2 <- NULL
  }
  message("fireSense_ignitionFit: year climate term (", variable, "): beta = ", signif(beta, 3),
          if (hinge) paste0(", beta2 = ", signif(beta2, 3)), " from ", nrow(d), " of ", nrow(m), " fit years.")
  yc <- list(variable = variable, center = center, scale = scale, beta = beta)
  if (hinge) yc$beta2 <- beta2
  yc$norm <- mean(fireSenseUtils::yearClimateMultiplier(yc, m$mean))
  if (weights != "equal") yc$weights <- weights   # recorded so a stored fit can be compared with a run's settings
  yc[c("variable", "center", "scale", "beta", "norm", if (hinge) "beta2", if (weights != "equal") "weights")]
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
#' A fit stored in the ledger carries the term it was made with (`modelList$yearClimate`; `beta2` with
#' `yearClimateHinge`, `weights` with `yearClimateWeights = "climate"`). The escape term must be there exactly when
#' `escapeYearClimate` is set (same variable, with or without underscores). Where both the stored fit and this run have a
#' term, `yearClimateHinge` and `yearClimateWeights` must agree. A fit stored before these parameters existed has
#' neither `beta2`, `weights` nor an escape term: it matches only their defaults.
#'
#' @param fitted A stored `list(modelList = ...)`, or `NULL`.
#' @param igOrEsc `"ignition"` or `"escape"`.
#' @param par The module's parameters (`P(sim)`): `ignitionYearClimate`, `escapeYearClimate`, `yearClimateHinge`,
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
  if (!identical(!is.null(yc$beta2), isTRUE(par$yearClimateHinge)))
    return(paste0("it was fitted ", if (is.null(yc$beta2)) "without" else "with", " the year climate hinge, but ",
                  "`yearClimateHinge` is ", isTRUE(par$yearClimateHinge)))
  stored <- if (is.null(yc$weights)) "equal" else yc$weights
  now <- if (is.null(par$yearClimateWeights)) "equal" else par$yearClimateWeights
  if (!identical(stored, now))
    return(paste0("its year climate term was fitted with `yearClimateWeights = \"", stored, "\"`, but this run's is \"", now, "\""))
  NULL
}
