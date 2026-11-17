## fittedYearEffects(): per fit year, the multipliers of the expected ignitions and escapes with which a hindcast's
## own draw (fireSense_burn) expects the year's target counts: the observed, shrunk toward the prediction

## 3 years x 2000 pixels; the predictions are known, the observed totals are the predictions times a year effect
syntheticYears <- function(uIgn = c(2, 0.5, 1.3), uEsc = c(0.6, 1.5, 1), core = TRUE) {
  nPix <- 2000
  i <- seq_len(nPix)
  years <- 2001:2003
  ign <- data.table::data.table(pixelID = rep(i, 3), year = rep(years, each = nPix),
                                mean = rep(0.2 + (i %% 7) / 2, 3),             # 0.2 .. 3.2 per pixel
                                weight = if (core) rep(ifelse(i %% 4 == 0, 0, 0.5 + (i %% 3) / 4), 3) else 1)
  esc <- data.table::data.table(pixelID = ign$pixelID, year = ign$year, mean = rep(0.2 + (i %% 5) / 5, 3))
  ## observed: what the draw expects under the true effects, rounded (large totals: little shrinkage)
  lam <- ign$mean * uIgn[match(ign$year, years)]
  ue <- uEsc[match(esc$year, years)]
  expEsc <- ign$weight * lam * fireSenseUtils::escapeProbPerIgnition(ue * esc$mean, lam)
  totI <- round(tapply(ign$weight * lam, ign$year, sum)); totE <- round(tapply(expEsc, esc$year, sum))
  ## each year's total on its first pixel: only the totals matter
  first <- !duplicated(ign$year)
  data.table::set(ign, NULL, "observed", ifelse(first, totI[as.character(ign$year)], 0))
  data.table::set(esc, NULL, "observed", ifelse(first, totE[as.character(esc$year)], 0))
  list(ign = ign, esc = esc, years = years, totI = as.numeric(totI), totE = as.numeric(totE))
}

## what fireSense_burn's draw expects in the core, with the multipliers
drawExpects <- function(d, eff) {
  uI <- eff$ignition$multiplier[match(d$ign$year, eff$ignition$year)]
  uE <- eff$escape$multiplier[match(d$esc$year, eff$escape$year)]
  lam <- d$ign$mean * uI
  list(ignitions = as.numeric(tapply(d$ign$weight * lam, d$ign$year, sum)),
       escapes = as.numeric(tapply(d$ign$weight * lam * fireSenseUtils::escapeProbPerIgnition(uE * eff$escapeCalibration * d$esc$mean, lam),
                                   d$esc$year, sum)))
}

test_that("with the fitted effects the draw expects each year's target, close to the observed when counts are large", {
  d <- syntheticYears()
  eff <- fittedYearEffects(d$ign, d$esc)
  expect_identical(eff$ignition$year, as.numeric(d$years))
  expect_identical(eff$ignition$observed, d$totI)
  expect_identical(eff$escape$observed, d$totE)
  got <- drawExpects(d, eff)
  ## the target: the expected count times the gamma-Poisson posterior mean of the multiplier
  post <- function(e) { a <- 1 / expm1(e$priorSD[1]^2); e$expected * (a + e$observed) / (a + e$expected) }
  expect_equal(got$ignitions, post(eff$ignition), tolerance = 1e-10)
  expect_equal(got$escapes, post(eff$escape), tolerance = 1e-8)
  ## thousands of events a year: the target is within 1% of the observed
  expect_equal(got$ignitions, d$totI, tolerance = 0.01)
  expect_equal(got$escapes, d$totE, tolerance = 0.01)
  expect_equal(eff$ignition$multiplier, c(2, 0.5, 1.3), tolerance = 0.01)
  ## the multipliers are the true effects relative to the calibration (the in-sample total of escapes is not the models')
  expect_equal(eff$escape$multiplier * eff$escapeCalibration, c(0.6, 1.5, 1), tolerance = 0.01)
  ## the prior sd is the in-sample spread of the totals around the predictions
  expect_identical(eff$ignition$priorSD, rep(fireSenseUtils::yearEffectSD(d$totI, eff$ignition$expected), 3))
})

test_that("a year with little expected is shrunk toward 1, one with nothing expected keeps 1", {
  d <- syntheticYears()
  ## year 2003: almost no expected ignitions, yet 3 observed
  y3 <- which(d$ign$year == 2003)
  data.table::set(d$ign, y3, "mean", 1e-4)
  data.table::set(d$ign, y3[1], "observed", 3)
  eff <- fittedYearEffects(d$ign, d$esc)
  e3 <- eff$ignition[eff$ignition$year == 2003, ]
  expect_lt(e3$multiplier, e3$observed / e3$expected / 10)   # far below the raw ratio
  expect_gt(e3$multiplier, 1)
  ## nothing expected at all
  data.table::set(d$ign, y3, "mean", 0)
  eff <- fittedYearEffects(d$ign, d$esc)
  expect_identical(eff$ignition$multiplier[3], 1)
  expect_identical(eff$escape$multiplier[3], 1)
})

test_that("no more spread than Poisson around the predictions: every multiplier is 1", {
  d <- syntheticYears(uIgn = c(1, 1, 1), uEsc = c(1, 1, 1))
  eff <- fittedYearEffects(d$ign, d$esc)
  expect_identical(eff$ignition$priorSD, rep(0, 3))
  expect_identical(eff$ignition$multiplier, rep(1, 3))
  expect_identical(eff$escape$multiplier, rep(1, 3))
})

test_that("an escape target above the ignitions gives the multiplier at which every ignition escapes", {
  d <- syntheticYears()
  data.table::set(d$esc, which(d$esc$year == 2001)[1], "observed", 1e6)
  eff <- fittedYearEffects(d$ign, d$esc)
  got <- drawExpects(d, eff)
  expect_equal(got$escapes[1], drawExpects(d, eff)$ignitions[1])   # all ignitions escape
})


test_that("the posterior is stored: a rep's multiplier * qgamma(p, shape, shape) is a draw of Gamma(a + O, a + mu)", {
  d <- syntheticYears()
  ## a small year, so the posterior is wide
  data.table::set(d$ign, which(d$ign$year == 2003), "mean", 0.002)
  data.table::set(d$ign, which(d$ign$year == 2003)[1], "observed", 4)
  eff <- fittedYearEffects(d$ign, d$esc)
  e <- eff$ignition
  a <- 1 / expm1(e$priorSD[1]^2)
  expect_equal(e$shape, a + e$observed)
  p <- c(0.1, 0.5, 0.9)
  ## fireSense_burn's draw of the year-2003 multiplier is the posterior's quantile
  expect_equal(e$multiplier[3] * stats::qgamma(p, e$shape[3], e$shape[3]),
               stats::qgamma(p, a + e$observed[3], a + e$expected[3]))
  expect_true(all(is.finite(eff$escape$shape)))
  ## no uncertainty where the multiplier is fixed at 1
  d0 <- syntheticYears(uIgn = c(1, 1, 1), uEsc = c(1, 1, 1))
  expect_true(all(is.na(fittedYearEffects(d0$ign, d0$esc)$ignition$shape)))
})

## The models' ignitions are spread evenly over the pixels, the observed ones sit in the pixels where an ignition escapes
## more often: the escape fit, on observed ignitions only, cannot see that
locationData <- function(obsEscapes = 0.3) {
  nPix <- 2000
  years <- 2001:2003
  i <- rep(seq_len(nPix), 3)
  high <- i %% 2 == 0
  ign <- data.table::data.table(pixelID = i, year = rep(years, each = nPix), mean = 0.5, weight = 1,
                                observed = as.numeric(high))
  esc <- data.table::data.table(pixelID = i, year = ign$year, mean = ifelse(high, 0.3, 0.1),
                                observed = ifelse(high, obsEscapes, 0))
  list(ign = ign, esc = esc)
}

test_that("ignitions that fall where escapes are rare give a calibration above 1, and the calibrated draw expects the observed escapes", {
  d <- locationData()
  eff <- fittedYearEffects(d$ign, d$esc)
  ## uncalibrated, the models' ignitions expect sum(w * lambda * q) of the observed escapes
  uncal <- sum(d$ign$weight * d$ign$mean * fireSenseUtils::escapeProbPerIgnition(d$esc$mean, d$ign$mean))
  expect_lt(uncal, sum(d$esc$observed))
  expect_gt(eff$escapeCalibration, 1)
  expect_equal(eff$escapeCalibration, sum(d$esc$observed) / uncal, tolerance = 1e-8)
  ## the calibrated draw (fireSense_ignitionPredict multiplies the escape mean by it) reproduces the observed escapes
  cal <- sum(d$ign$weight * d$ign$mean *
               fireSenseUtils::escapeProbPerIgnition(eff$escapeCalibration * d$esc$mean, d$ign$mean))
  expect_equal(cal, sum(d$esc$observed), tolerance = 1e-8)
})

test_that("the year multipliers of the escapes centre on 1 with the calibration, and not without it", {
  d <- locationData()
  eff <- fittedYearEffects(d$ign, d$esc)
  ## this year-balanced toy has no year effect: the multipliers are 1 whatever c
  expect_equal(eff$escape$multiplier, rep(1, 3), tolerance = 1e-6)
  ## expected totals of the escapes are the calibrated ones
  expect_equal(eff$escape$expected, eff$escape$observed, tolerance = 1e-8)
  ## a year with twice the escapes is above 1 and the others below it, relative to the calibrated total
  data.table::set(d$esc, which(d$esc$year == 2002 & d$esc$observed > 0), "observed", 0.6)
  eff2 <- fittedYearEffects(d$ign, d$esc)
  expect_gt(eff2$escape$multiplier[2], 1)
  expect_lt(eff2$escape$multiplier[1], 1)
})

test_that("a fit with no escape gives a calibration of 1, and the cap is respected", {
  d <- syntheticYears(uIgn = c(1, 1, 1), uEsc = c(1, 1, 1))
  expect_equal(fittedYearEffects(d$ign, d$esc)$escapeCalibration, 1, tolerance = 0.01)
  ## more escapes observed than the ignitions can give: c is where every ignition escapes, no higher
  d <- locationData(obsEscapes = 5)
  eff <- fittedYearEffects(d$ign, d$esc)
  q <- fireSenseUtils::escapeProbPerIgnition(eff$escapeCalibration * d$esc$mean, d$ign$mean)
  expect_lte(max(q), 1)
  expect_equal(max(q), 1)
  expect_lte(sum(d$ign$weight * d$ign$mean * q), sum(d$ign$weight * d$ign$mean))
})
