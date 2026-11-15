#' The fitted year effects of ignitions and escapes, for hindcasts
#'
#' A hindcast of a fit year may use that year's fitted effect instead of a draw of `yearSD`
#' (`fireSense_burn`'s `yearEffects = "fitted"`); forecasts and NRV keep the draw. For each fit year `y` this
#' gives the multiplier of the expected ignitions, `u_I`, and of the expected escapes of an ignited pixel,
#' `u_E`, such that the hindcast's own draw (`fireSense_burn`: Poisson ignitions of mean `u_I * lambda`, each
#' escaping with `fireSenseUtils::escapeProbPerIgnition(u_E * e, u_I * lambda)`, counted in the core) expects
#' the year's target counts. `lambda` and `e` are the fold-mean in-sample predictions
#' (`fireSenseUtils::foldMeanPrediction()`, as `fireSense_ignitionPredict` predicts) of every pixel of the
#' year, so the multipliers absorb at the year level both the models' in-sample bias and where the
#' ignitions land (on pixels with a lower or higher escape mean than the observed ignitions).
#'
#' The target is not the observed count itself: a year with few expected events says little about its
#' effect. With `O_y ~ Poisson(u_y * mu_y)` and a mean-one gamma prior on `u_y` of shape and rate
#' `a = 1 / (exp(s^2) - 1)` (the variance of `fireSenseUtils::yearEffectMultiplier(s, z)`), the posterior
#' mean is `(a + O_y) / (a + mu_y)`, and the target is `mu_y` times it. `s` is
#' `fireSenseUtils::yearEffectSD(O, mu)` of the same in-sample totals (empirical Bayes): the spread of the
#' observed totals around these predictions, not the fit's `yearSD`, which is out of year and larger. A year
#' with nothing expected keeps 1; `s` 0 (no more spread than Poisson) gives 1 everywhere. Escapes are
#' solved after ignitions, with `u_I` applied, and cannot exceed the ignitions (the cap in
#' `escapeProbPerIgnition()`): a target above that gives the multiplier at which every ignition escapes.
#'
#' @param ignitions `data.table` of the ignition fit data, one row per pixel and year: `pixelID`, `year`,
#'   `observed` (ignitions), `mean` (the predicted ignitions of the whole pixel) and `weight` (its share in
#'   the core, `exposureOf()`).
#' @param escapes `data.table` of the escape fit data: `pixelID`, `year`, `observed` (escapes) and `mean`
#'   (the predicted escapes if ignited), for the same pixels and years.
#' @return A list of two `data.frame`s, `ignition` and `escape`, one row per year: `year`, `observed` and
#'   `expected` (the core totals with no effect), `priorSD` (`s`) and `multiplier`.
fittedYearEffects <- function(ignitions, escapes) {
  ign <- yearPosterior(ignitions$year, ignitions$observed, ignitions$weight * ignitions$mean)
  rows <- match(paste(escapes$pixelID, escapes$year), paste(ignitions$pixelID, ignitions$year))
  if (anyNA(rows))
    stop("fittedYearEffects(): ", sum(is.na(rows)), " escape rows have no ignition row of the same pixel and year")
  lambda <- ignitions$mean[rows] * ign$multiplier[match(escapes$year, ign$year)]
  w <- ignitions$weight[rows]
  coreEscapes <- function(u, i) sum(w[i] * lambda[i] * fireSenseUtils::escapeProbPerIgnition(u * escapes$mean[i], lambda[i]))
  byYear <- split(seq_along(rows), escapes$year)
  esc <- yearPosterior(as.numeric(names(byYear)), vapply(byYear, function(i) sum(escapes$observed[i]), numeric(1)),
                       vapply(byYear, function(i) coreEscapes(1, i), numeric(1)), summed = TRUE)
  esc$multiplier <- mapply(function(i, target) solveMultiplier(function(u) coreEscapes(u, i), target),
                           byYear, esc$expected * esc$multiplier, USE.NAMES = FALSE)
  list(ignition = ign, escape = esc)
}

## Per year: observed and expected totals, the prior sd and the posterior-mean multiplier (see fittedYearEffects()).
## `summed`: `observed` and `expected` are already one per year (in `year`), else one per row.
yearPosterior <- function(year, observed, expected, summed = FALSE) {
  if (!summed) {
    tot <- rowsum(cbind(observed = as.numeric(observed), expected = expected), group = year)
    year <- as.numeric(rownames(tot)); observed <- tot[, "observed"]; expected <- tot[, "expected"]
  }
  s <- fireSenseUtils::yearEffectSD(observed, expected)
  a <- 1 / expm1(s^2)
  multiplier <- if (s > 0) (a + observed) / (a + expected) else rep(1, length(observed))
  multiplier[expected <= 0] <- 1
  data.frame(year = year, observed = unname(observed), expected = unname(expected), priorSD = s,
             multiplier = unname(multiplier), row.names = NULL)
}

## The u at which the increasing `f(u)` (expected core escapes) is `target`; f is linear until the cap, so
## the root is between 0 and the u at which every pixel is capped (where f stops rising).
solveMultiplier <- function(f, target) {
  if (f(1) <= 0) return(1)
  hi <- 1
  while (f(2 * hi) > f(hi) * (1 + 1e-12) && f(hi) < target) hi <- 2 * hi
  if (f(hi) <= target) return(hi)
  stats::uniroot(function(u) f(u) - target, c(0, hi), tol = 1e-10)$root
}
