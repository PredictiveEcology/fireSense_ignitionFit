## concentrationOf() and chainedEscapeFactor(), on toys with known answers

## 100 pixel-years: the 10 highest predictions are 0.5 each, the other 90 are 5/90 each (10 events predicted);
## 8 observed events fall in the top 10 pixel-years and 2 in the rest
toyConcentration <- function() {
  expected <- c(rep(0.5, 10), rep(5 / 90, 90))
  observed <- c(rep(1, 8), 0, 0, 1, 1, rep(0, 88))
  list(observed = observed, expected = expected)
}

test_that("concentrationOf: shares of events in the top pixel-years, ties shared equally", {
  d <- toyConcentration()
  co <- concentrationOf(d$observed, d$expected, basis = "out-of-year")
  expect_identical(co$basis, "out-of-year")
  expect_identical(co$n, 100L)
  expect_equal(c(co$observed, co$expected), c(10, 10))
  expect_equal(co$lift$top, c(0.01, 0.05, 0.10, 0.25))
  ## the top 10 are tied: 1% and 5% of the pixel-years are 1/10 and 1/2 of them
  expect_equal(co$lift$observed, c(0.08, 0.40, 0.80, (8 + 2 * 15 / 90) / 10))
  expect_equal(co$lift$predicted, c(0.05, 0.25, 0.50, (5 + 15 * 5 / 90) / 10))
  ## ranking by the prediction: moving the rows does not change it
  p <- c(31:100, 1:30)
  expect_equal(concentrationOf(d$observed[p], d$expected[p])$lift, co$lift)
})

test_that("concentrationOf: Gini index of the observed and the predicted counts", {
  d <- toyConcentration()
  co <- concentrationOf(d$observed, d$expected)
  ## mean absolute difference over twice the mean: 10 ones and 90 zeros give 0.9; 10 of 0.5 and 90 of 1/18 give 0.4
  expect_equal(unname(co$gini), c(0.9, 0.4))
  flat <- concentrationOf(d$observed, rep(0.1, 100))
  expect_equal(unname(flat$gini["predicted"]), 0)
  expect_equal(flat$lift$predicted, flat$lift$top)             # nothing is concentrated: the top x% hold x%
  expect_equal(flat$lift$observed, flat$lift$top)
  expect_true(is.na(giniIndex(numeric(5))))
})

test_that("concentrationOf: the Poisson deviance is 0 when the prediction is the observation", {
  d <- toyConcentration()
  expect_equal(concentrationOf(d$observed, d$observed)$deviance, 0, tolerance = 1e-6)
  mu <- d$expected
  y <- d$observed
  manual <- mean(2 * (ifelse(y > 0, y * log(y / mu), 0) - (y - mu)))
  expect_equal(concentrationOf(y, mu)$deviance, manual)
})

## the escape probability of an ignition, from its definition (fireSenseUtils::escapeProbPerIgnition)
qOf <- function(e, lambda) pmin(1, e * (1 - exp(-lambda)) / lambda)

toyChain <- function(observedIgnitions) {
  n <- 6
  ign <- data.table::data.table(pixelID = 1:n, year = 2001L, observed = observedIgnitions,
                                mean = c(0.9, 0.9, 0.3, 0.3, 0.1, 0.1), weight = c(1, 0.5, 1, 1, 0.5, 1))
  esc <- data.table::data.table(pixelID = 1:n, year = 2001L, observed = c(0, 0, 0, 0, 0, 0),
                                mean = c(0.9, 0.2, 0.5, 0.1, 0.8, 0.4))
  list(ign = ign, esc = esc)
}

test_that("chainedEscapeFactor: c is the observed escapes over those of model-placed ignitions, in three factors", {
  d <- toyChain(c(1, 0, 0, 1, 0, 1))
  data.table::set(d$esc, NULL, "observed", c(0, 0, 0, 1, 0, 0))
  r <- chainedEscapeFactor(d$ign, d$esc)
  lambda <- d$ign$mean; q <- qOf(d$esc$mean, lambda); w <- d$ign$weight
  expect_equal(r$expectedIgnitions, sum(w * lambda))
  expect_equal(r$expectedEscapes, sum(w * lambda * q))
  expect_equal(c(r$observedIgnitions, r$observedEscapes), c(3, 1))
  expect_equal(r$c, 1 / sum(w * lambda * q))
  expect_equal(r$qModel, sum(w * lambda * q) / sum(w * lambda))
  expect_equal(r$qObserved, sum(d$ign$observed * q) / 3)
  expect_equal(r$c, r$ignitionRatio * r$escapeCalibration * r$placement)
  expect_equal(r$ignitionRatio, 3 / sum(w * lambda))
  expect_equal(r$placement, r$qObserved / r$qModel)
})

test_that("chainedEscapeFactor: placement is 1 when the observed ignitions fall as the model places them, above 1 when they sit where escaping is unlikely", {
  d <- toyChain(c(1, 0, 0, 0, 0, 0))
  lambda <- d$ign$mean; w <- d$ign$weight
  ## ignitions in proportion to the model's: the model-weighted and the observed-weighted mean q agree
  proportional <- data.table::copy(d$ign)
  data.table::set(proportional, NULL, "observed", w * lambda * 100)
  expect_equal(chainedEscapeFactor(proportional, d$esc)$placement, 1)
  ## all the observed ignitions in the pixel with the lowest q: observed escapes come out above what the model places
  q <- qOf(d$esc$mean, lambda)
  lowest <- data.table::copy(d$ign)
  data.table::set(lowest, NULL, "observed", as.numeric(seq_len(6) == which.min(q)) * 5)
  expect_lt(chainedEscapeFactor(lowest, d$esc)$qObserved, chainedEscapeFactor(lowest, d$esc)$qModel)
  highest <- data.table::copy(d$ign)
  data.table::set(highest, NULL, "observed", as.numeric(seq_len(6) == which.max(q)) * 5)
  expect_gt(chainedEscapeFactor(highest, d$esc)$placement, 1)
})

test_that("chainedEscapeFactor stops for an escape row without an ignition row; the messages carry the numbers", {
  d <- toyChain(c(1, 0, 0, 1, 0, 1))
  data.table::set(d$esc, 3L, "pixelID", 99L)
  expect_error(chainedEscapeFactor(d$ign, d$esc), "no ignition row")
  d <- toyChain(c(1, 0, 0, 1, 0, 1))
  data.table::set(d$esc, NULL, "observed", c(0, 0, 0, 1, 0, 0))
  r <- chainedEscapeFactor(d$ign, d$esc)
  expect_match(chainMessage(r), paste0("c = ", format(round(r$c, 3))), fixed = TRUE)
  dd <- toyConcentration()
  expect_match(concentrationMessage(concentrationOf(dd$observed, dd$expected, "out-of-year"), "ignition"),
               "ignition concentration (out-of-year, 100 pixel-years)", fixed = TRUE)
})

test_that("a fit of both processes stores its concentration, the escape's with the chained-escape factor, and logs it", {
  paths <- localModulePaths()
  msgs <- new.env()
  sim <- runModuleMsgs(paths, msgs)
  ig <- sim$fireSense_IgnitionFitted$modelList$concentration
  es <- sim$fireSense_EscapeFitted$modelList$concentration
  expect_identical(c(ig$basis, es$basis), c("out-of-year", "out-of-year"))
  expect_identical(ig$n, 400L)
  expect_identical(es$n, 61L)                                         # the ignited pixel-years
  expect_equal(ig$observed, 83)
  expect_equal(es$observed, 33)
  expect_null(ig$chain)
  expect_named(es$chain, c("c", "ignitionRatio", "escapeCalibration", "placement", "qObserved", "qModel",
                           "observedIgnitions", "expectedIgnitions", "observedEscapes", "expectedEscapes"))
  expect_equal(es$chain$c, es$chain$ignitionRatio * es$chain$escapeCalibration * es$chain$placement)
  expect_equal(es$chain$observedEscapes, 33)
  ## the booster list holds the boosters, not the diagnostic
  expect_null(sim$fireSense_IgnitionFitted$modelList$model$concentration)
  expect_length(grep("fireSense_ignitionFit ignition concentration", msgs$text), 1L)
  expect_length(grep("fireSense_ignitionFit escape concentration", msgs$text), 1L)
  expect_length(grep("chained-escape factor c = ", msgs$text), 1L)
  ## one process: its own concentration and no chain
  one <- runModuleMsgs(localModulePaths(), new.env(), params = list(whichProcessesToFit = "ignition"),
                       drop = "fireSense_escapeCovariates")
  expect_null(one$fireSense_IgnitionFitted$modelList$concentration$chain)
  expect_false(is.null(one$fireSense_IgnitionFitted$modelList$concentration))
})

