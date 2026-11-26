## The year climate term of the ignitions: E[O_y] = P_y * exp(beta * z_y) / norm, beta >= 0, z_y the standardised
## yearly mean of `CMDsm` (fitYearClimate(), yearClimateRowMultiplier()), and what it does to the fit.

## 5 years x 3 pixels. The yearly mean of `CMDsm` is 1, 2, 3, 4, 5 on the pixels that count (coreFraction 1) and
## is not changed by the pixel outside the core (coreFraction 0, CMDsm 100).
syntheticClimate <- function() {
  data.table::data.table(pixelID = rep(1:3, 5), year = rep(2001:2005, each = 3),
                         CMDsm = c(rbind(1:5, 1:5, 100)),
                         coreFraction = rep(c(1, 1, 0), 5), ignitions = 0L)
}

test_that("the z basis is the unweighted yearly mean, centred and scaled over the fit years", {
  raw <- syntheticClimate()
  m <- fireSenseUtils::yearClimateMeans(raw, "CMD_sm")
  expect_equal(m$year, 2001:2005)
  expect_equal(m$mean, as.numeric(1:5))
  yt <- data.frame(year = 2001:2005, observed = c(10, 20, 30, 40, 50), predicted = 30)
  yc <- fitYearClimate("CMDsm", raw, yt)
  expect_identical(yc$variable, "CMDsm")
  expect_equal(yc$center, 3)
  expect_equal(yc$scale, stats::sd(1:5))
  expect_equal(yearClimateRowMultiplier(yc, raw)[raw$pixelID == 1],
               exp(yc$beta * (1:5 - 3) / stats::sd(1:5)) / yc$norm)
})

test_that("covariates built from the fit table for one year give the fit's multiplier for that year", {
  ## what fireSense_ignitionPredict is given: that year's rows, with an NA row, as the fit table has
  raw <- syntheticClimate()
  data.table::set(raw, 4L, "CMDsm", NA_real_)      # a pixel of 2002 with no CMDsm
  yt <- data.frame(year = 2001:2005, observed = c(10, 20, 30, 40, 50), predicted = 30)
  yc <- fitYearClimate("CMDsm", raw, yt)
  for (y in c(2001, 2003, 2005)) {
    covs <- raw[raw$year == y, ]
    expect_equal(fireSenseUtils::yearClimateMultiplier(yc, fireSenseUtils::yearClimateMeans(covs, "CMDsm")$mean),
                 yearClimateRowMultiplier(yc, raw)[raw$year == y][1])
  }
})

test_that("beta is the Poisson fit on the out-of-year totals, and never below 0", {
  raw <- syntheticClimate()
  z <- (1:5 - 3) / stats::sd(1:5)
  up <- data.frame(year = 2001:2005, predicted = 100, observed = round(100 * exp(0.4 * z)))
  expect_equal(fitYearClimate("CMDsm", raw, up)$beta, 0.4, tolerance = 0.05)
  ## a level gap (observed 2.5 times predicted in every year) is not climate: beta is the same
  expect_equal(fitYearClimate("CMDsm", raw, transform(up, observed = round(250 * exp(0.4 * z))))$beta, 0.4, tolerance = 0.05)
  expect_identical(fitYearClimate("CMDsm", raw, transform(up, observed = 250))$beta, 0)
  down <- data.frame(year = 2001:2005, predicted = 100, observed = round(100 * exp(-0.4 * z)))
  expect_identical(fitYearClimate("CMDsm", raw, down)$beta, 0)
})

test_that("beta is 0 when nothing is observed or expected, or there are no totals", {
  raw <- syntheticClimate()
  expect_identical(fitYearClimate("CMDsm", raw, data.frame(year = 2001:2005, observed = 0, predicted = 5))$beta, 0)
  expect_identical(fitYearClimate("CMDsm", raw, data.frame(year = 2001:2005, observed = 5, predicted = 0))$beta, 0)
  expect_identical(fitYearClimate("CMDsm", raw, NULL)$beta, 0)
})

test_that("the multiplier averages 1 over the fit years (norm)", {
  raw <- syntheticClimate()
  z <- (1:5 - 3) / stats::sd(1:5)
  yc <- fitYearClimate("CMDsm", raw, data.frame(year = 2001:2005, predicted = 100, observed = round(100 * exp(0.4 * z))))
  expect_gt(yc$beta, 0)
  expect_equal(mean(yearClimateRowMultiplier(yc, raw, 2001:2005)), 1)
})

test_that("NA in the covariate, or a year with no rows, never gives a non-finite centre or scale", {
  raw <- syntheticClimate()
  data.table::set(raw, which(raw$year == 2003), "CMDsm", NA_real_)    # no row left in 2003
  data.table::set(raw, 1L, "CMDsm", NA_real_)                          # one row of 2001 left
  yt <- data.frame(year = 2001:2005, observed = c(10, 20, 30, 40, 50), predicted = 30)
  yc <- fitYearClimate("CMDsm", raw, yt)
  expect_true(all(is.finite(unlist(yc[c("center", "scale", "beta", "norm")]))))
  expect_equal(yc$center, mean(c(1, 2, 4, 5)))
  expect_equal(yearClimateRowMultiplier(yc, raw, 2003), 1)             # a year without a mean: no term
  ## nothing usable: no term, with a message, and finite
  allNA <- syntheticClimate(); data.table::set(allNA, NULL, "CMDsm", NA_real_)
  expect_message(yc0 <- fitYearClimate("CMDsm", allNA, yt), "no year climate term")
  expect_identical(c(yc0$beta, yc0$norm, yc0$center, yc0$scale), c(0, 1, 0, 1))
  expect_message(fitYearClimate("CMDsm", syntheticClimate(), NULL), "no year climate term")
})

test_that("raw covariates without coreFraction give a warning: the mean would be over the buffer too", {
  raw <- syntheticClimate(); data.table::set(raw, NULL, "coreFraction", NULL)
  yt <- data.frame(year = 2001:2005, observed = c(10, 20, 30, 40, 50), predicted = 30)
  expect_warning(fitYearClimate("CMDsm", raw, yt), "no `coreFraction` column")
  expect_no_warning(fitYearClimate("CMDsm", syntheticClimate(), yt))
})

test_that("calendarYearOf stops when a year is not among the covariates' years", {
  expect_error(calendarYearOf(c(0, 0.5), c(0, 1), c(2001, 2002)), "1 of 2 years")
})

test_that("a variable that is not a covariate stops", {
  expect_error(fitYearClimate("CMDfoo", syntheticClimate(), NULL), "not among the covariates")
})

test_that("calendarYearOf maps standardised years, also through rownames(), and does not assume every year has rows", {
  raw <- c(2001, 2001, 2003, 2003, 2007)
  fit <- (raw - mean(raw)) / stats::sd(raw)
  through <- as.numeric(rownames(rowsum(rep(1, 5), fit)))   # as outOfFoldYearTotals() makes them
  expect_identical(calendarYearOf(through, fit, raw), c(2001, 2003, 2007))
  expect_identical(calendarYearOf(2003, 2003, 2003), 2003)  # rescaleVars FALSE
})

test_that("no term is a multiplier of 1", {
  expect_equal(yearClimateRowMultiplier(NULL, syntheticClimate()), rep(1, 15))
})

test_that("yearSD shrinks when the climate term explains the totals", {
  raw <- syntheticClimate()
  z <- (1:5 - 3) / stats::sd(1:5)
  yt <- data.frame(year = 2001:2005, predicted = 1000)
  yt$observed <- round(1000 * exp(0.5 * z))
  yc <- fitYearClimate("CMDsm", raw, yt)
  withTerm <- fireSenseUtils::yearEffectSD(yt$observed, yt$predicted * yearClimateRowMultiplier(yc, raw, yt$year))
  expect_lt(withTerm, fireSenseUtils::yearEffectSD(yt$observed, yt$predicted) / 3)
})

## The yearly mean of CMDsm falls with the year (3.5 down to 0.5) and the ignitions fall with it, 8:4:2:1, at the
## same pixel pattern: out of year the model under-predicts the dry first year most, so beta > 0.
climateYears <- function() {
  dt <- makeIgnitionCovariates()
  k <- c(8, 4, 2, 1)[match(dt$year, 2001:2004)]
  data.table::set(dt, NULL, "ignitions", as.integer((dt$youngAge == 1 & dt$pixelID %% 3 == 0) * k))
  data.table::set(dt, NULL, "CMDsm", dt$CMDsm * 0.1 + (2004 - dt$year))
  dt
}

fitWith <- function(par) {
  esc <- climateYears()
  data.table::set(esc, NULL, "escapes", as.integer(esc$ignitions > 0 & esc$lightning >= 6))
  root <- withr::local_tempdir(.local_envir = parent.frame())
  paths <- testPaths
  paths$cachePath <- file.path(root, "cache")
  paths$outputPath <- file.path(root, "outputs")
  sim <- NULL
  utils::capture.output(suppressWarnings(suppressMessages({
    sim <- SpaDES.core::simInit(times = list(start = 1, end = 1), modules = moduleName, paths = paths,
                                objects = list(fireSense_ignitionCovariates = climateYears(), fireSense_escapeCovariates = esc,
                                               ignitionFitRTM = makeIgnitionFitRTM()),
                                params = stats::setNames(list(testParams(c(list(.plots = NA), par))), moduleName))
    sim <- SpaDES.core::spades(sim)
  })))
  sim
}

test_that("the fit stores yearClimate, maps the standardised years, and NA leaves the fit as before", {
  withTerm <- fitWith(list())
  without <- fitWith(list(ignitionYearClimate = NA_character_))
  yc <- withTerm$fireSense_IgnitionFitted$modelList$yearClimate
  expect_named(yc, c("variable", "center", "scale", "beta", "norm"))
  expect_identical(yc$variable, "CMDsm")
  expect_gt(yc$beta, 0)
  m <- fireSenseUtils::yearClimateMeans(climateYears(), "CMDsm")
  expect_equal(c(yc$center, yc$scale), c(mean(m$mean), stats::sd(m$mean)))

  ## NA: no yearClimate, the names are the ones of a fit made before the parameter
  expect_null(without$fireSense_IgnitionFitted$modelList$yearClimate)
  expect_identical(names(without$fireSense_IgnitionFitted$modelList),
                   c("model", "rescales", "fittingRes", "lambdaRescaleFactor", "yearSD", "concentration", "yearEffects"))
  ## the boosters do not depend on the term; yearSD is the spread left after it
  yt <- without$fireSense_IgnitionFitted$modelList$model$yearTotals
  expect_identical(without$fireSense_IgnitionFitted$modelList$yearSD, fireSenseUtils::yearEffectSD(yt$observed, yt$predicted))
  expect_equal(withTerm$fireSense_IgnitionFitted$modelList$model$yearTotals, yt)
  expect_lt(withTerm$fireSense_IgnitionFitted$modelList$yearSD, without$fireSense_IgnitionFitted$modelList$yearSD)

  ## the fitted year effects are relative to the term: the expected totals are those of the fit without it times exp(beta z)
  ei <- withTerm$fireSense_IgnitionFitted$modelList$yearEffects
  ei0 <- without$fireSense_IgnitionFitted$modelList$yearEffects
  expect_equal(ei$expected / ei0$expected,
               fireSenseUtils::yearClimateMultiplier(yc, m$mean[match(ei$year, m$year)]))
})

test_that("rescaleVars = FALSE gives the same term as the standardised fit, and a logical NA is no term", {
  yc <- fitWith(list())$fireSense_IgnitionFitted$modelList$yearClimate
  raw <- fitWith(list(rescaleVars = FALSE))$fireSense_IgnitionFitted$modelList$yearClimate
  expect_equal(raw[c("variable", "center", "scale")], yc[c("variable", "center", "scale")])
  expect_gt(raw$beta, 0)
  expect_equal(raw$beta, yc$beta, tolerance = 0.05)
  expect_null(fitWith(list(ignitionYearClimate = NA))$fireSense_IgnitionFitted$modelList$yearClimate)
})

test_that("a year climate variable that is not a covariate stops the fit", {
  expect_error(fitWith(list(ignitionYearClimate = "CMDfoo")), "not among the covariates")
})

## ---- the year climate term of the escapes, the logistic shape and the year weights ----

test_that("the defaults give the term as before, to the last bit", {
  raw <- syntheticClimate()
  z <- (1:5 - 3) / stats::sd(1:5)
  yt <- data.frame(year = 2001:2005, predicted = c(900, 1100, 1000, 950, 1050))
  yt$observed <- round(yt$predicted * exp(0.4 * z + c(0.1, -0.1, 0, 0.05, -0.05)))
  ## the code of the term before `shape` and `weights` existed
  d <- data.frame(observed = yt$observed, predicted = yt$predicted, z = (1:5 - 3) / stats::sd(1:5))
  b <- unname(stats::coef(stats::glm(observed ~ 0 + z, offset = log(fireSenseUtils::scaleToObservedLevel(observed, predicted)),
                                     family = stats::poisson, data = d))["z"])
  m <- 1:5; beta <- max(0, b)
  old <- list(variable = "CMDsm", center = mean(m), scale = stats::sd(m), beta = beta,
              norm = mean(exp(beta * (m - mean(m)) / stats::sd(m))))
  expect_identical(suppressMessages(fitYearClimate("CMDsm", raw, yt)), old)
  expect_identical(suppressMessages(fitYearClimate("CMDsm", raw, yt, shape = "loglinear", weights = "equal")), old)
  expect_identical(fireSenseUtils::yearClimateMultiplier(old, m), exp(old$beta * (m - old$center) / old$scale) / old$norm)
})

test_that("an escape term is fitted to the escapes' totals with the ignitions' yearly means", {
  raw <- syntheticClimate()
  z <- (1:5 - 3) / stats::sd(1:5)
  yt <- data.frame(year = 2001:2005, predicted = 500)
  yt$observed <- round(500 * exp(0.3 * z))
  yc <- suppressMessages(fitYearClimate("CMDsm", raw, yt))
  expect_equal(yc$beta, 0.3, tolerance = 0.05)
  expect_equal(c(yc$center, yc$scale), c(3, stats::sd(1:5)))     # the ignition covariates' means, not the escapes' rows
  ## escape rows of other years than the ignition rows' own get the multiplier of their year
  esc <- raw[raw$pixelID == 1 & raw$year %in% c(2002, 2005), ]
  expect_equal(yearClimateRowMultiplier(yc, raw, esc$year),
               fireSenseUtils::yearClimateMultiplier(yc, c(2, 5)))
})

test_that("the logistic shape recovers a known S-curve, saturates beyond the fit range, and is flat for a flat response", {
  raw <- syntheticClimate()
  z <- (1:5 - 3) / stats::sd(1:5)
  yt <- data.frame(year = 2001:2005, predicted = 1e5)
  curve <- transform(yt, observed = round(1e5 * (1 + 3 * stats::plogis(3 * z))))      # k = 3, ratio = 4
  yc <- suppressMessages(fitYearClimate("CMDsm", raw, curve, shape = "logistic"))
  expect_named(yc, c("variable", "center", "scale", "shape", "k", "ratio", "norm", "bounds"))
  expect_equal(yc$bounds, logisticYearBounds)
  expect_identical(yc$shape, "logistic")
  expect_equal(c(yc$k, yc$ratio), c(3, 4), tolerance = 0.05)
  expect_equal(mean(fireSenseUtils::yearClimateMultiplier(yc, 1:5)), 1)             # renormalised
  ## the fitted multipliers are the curve of the data, up to its level
  mult <- fireSenseUtils::yearClimateMultiplier(yc, 1:5)
  expect_equal(mult, curve$observed / mean(curve$observed), tolerance = 0.01)
  ## beyond the fit range it saturates at the floor and the upper level, over norm
  far <- fireSenseUtils::yearClimateMultiplier(yc, 3 + yc$scale * c(-50, -20, 20, 50))
  expect_equal(far[1:2], rep(1 / yc$norm, 2), tolerance = 1e-6)
  expect_equal(far[3:4], rep(yc$ratio / yc$norm, 2), tolerance = 1e-6)
  expect_lt(max(far), 3 * max(curve$observed) / mean(curve$observed))               # the upper level is bounded
  ## a flat response gives no term; a falling one too (k >= 0)
  flat <- suppressMessages(fitYearClimate("CMDsm", raw, transform(yt, observed = 1e5), shape = "logistic"))
  expect_equal(c(flat$k, flat$ratio, flat$norm), c(0, 1, 1))
  falling <- suppressMessages(fitYearClimate("CMDsm", raw, transform(yt, observed = round(1e5 * exp(-0.5 * z))), shape = "logistic"))
  expect_lt(falling$k * (falling$ratio - 1), 0.5)
  expect_identical(suppressMessages(fitYearClimate("CMDsm", raw, NULL, shape = "logistic"))$k, 0)
})

test_that("the logistic upper level over the average is bounded at 3 times the largest year over the mean", {
  raw <- syntheticClimate()
  ## one very dry year with a huge total: the unbounded fit would put the upper level far higher
  yt <- data.frame(year = 2001:2005, predicted = 100, observed = c(100, 100, 100, 100, 400))
  yc <- suppressMessages(fitYearClimate("CMDsm", raw, yt, shape = "logistic"))
  expect_lte(yc$ratio / yc$norm, 3 * max(yt$observed) / mean(yt$observed) + 1e-6)
})

test_that("the climate year weights are the shared rule (fireSenseUtils::annualYearWeights), and weight the extremes", {
  z <- c(-1.5, -0.1, 0, 0.5, 1.1)
  w <- fireSenseUtils::annualYearWeights(z, "climate")
  expect_equal(mean(w), 1)
  expect_equal(w / min(w), pmax(abs(z), 0.25) / 0.25)
  ## a response that is steep only in the dry extreme: the weighted slope follows the extreme years
  raw <- syntheticClimate()
  zz <- (1:5 - 3) / stats::sd(1:5)
  yt <- data.frame(year = 2001:2005, predicted = 1e4, observed = round(1e4 * exp(c(0, 0, 0, 0.1, 0.9))))
  eq <- suppressMessages(fitYearClimate("CMDsm", raw, yt))
  wt <- suppressMessages(fitYearClimate("CMDsm", raw, yt, weights = "climate"))
  expect_false(isTRUE(all.equal(eq$beta, wt$beta)))
  expect_gt(wt$beta, 0)
})

test_that("an escape term is stored, shares z with the ignition term, and applies to the fitted escape effects", {
  both <- fitWith(list(escapeYearClimate = "CMDsm"))
  ign <- both$fireSense_IgnitionFitted$modelList
  esc <- both$fireSense_EscapeFitted$modelList
  expect_named(esc$yearClimate, c("variable", "center", "scale", "beta", "norm"))
  expect_equal(esc$yearClimate[c("variable", "center", "scale")], ign$yearClimate[c("variable", "center", "scale")])
  ## default: the escapes have no term, as before
  none <- fitWith(list())$fireSense_EscapeFitted$modelList
  expect_null(none$yearClimate)
  expect_equal(esc$model$yearTotals, none$model$yearTotals)       # the boosters do not depend on it
  ## the fitted escape effects are relative to the term
  m <- fireSenseUtils::yearClimateMeans(climateYears(), "CMDsm")
  ## (not exactly: the expected escapes of an ignited pixel saturate with the mean, `escapeProbPerIgnition()`)
  expect_equal(esc$yearEffects$expected / none$yearEffects$expected,
               fireSenseUtils::yearClimateMultiplier(esc$yearClimate, m$mean[match(esc$yearEffects$year, m$year)]),
               tolerance = 1e-3)
  expect_identical(both$fireSense_IgnitionFitted$modelList$yearClimate, fitWith(list())$fireSense_IgnitionFitted$modelList$yearClimate)
})

test_that("the shape and the weights reach both fits", {
  fit <- fitWith(list(escapeYearClimate = "CMDsm", yearClimateShape = "logistic", yearClimateWeights = "climate"))
  for (f in list(fit$fireSense_IgnitionFitted, fit$fireSense_EscapeFitted)) {
    yc <- f$modelList$yearClimate
    expect_identical(yc$shape, "logistic")
    expect_identical(yc$weights, "climate")
    expect_gte(yc$k, 0)
  }
})

test_that("fitYearClimate records the weights only when they are not equal", {
  raw <- syntheticClimate()
  yt <- data.frame(year = 2001:2005, predicted = 1e4, observed = round(1e4 * exp(c(0, 0, 0, 0.1, 0.9))))
  expect_null(suppressMessages(fitYearClimate("CMDsm", raw, yt))$weights)
  expect_identical(suppressMessages(fitYearClimate("CMDsm", raw, yt, weights = "climate"))$weights, "climate")
  expect_identical(suppressMessages(fitYearClimate("CMDsm", raw, NULL, weights = "climate"))$weights, "climate")
})
