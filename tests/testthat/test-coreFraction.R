## The ignition counts start in the core of the ELF, but the training rows cover every coarse pixel, buffer
## included. A pixel that is half buffer can hold only half of its ignitions, so a fit that ignores that
## learns the mean rate over core and buffer and predicts too few ignitions for the core. `coreFraction`
## (fireSense_dataPrepFit) is the share of each pixel in the core; the fit is weighted by it.

## 2000 pixels x 4 years. Half the pixels have coreFraction 0.5, the rest 1; ignitions are drawn as
## Poisson(rate * coreFraction), where the rate of a whole pixel is 0.1, or 0.4 when x is high.
makeCoreFractionCovariates <- function(nPix = 2000, nYr = 4, seed = 1) {
  set.seed(seed)
  n <- nPix * nYr
  x <- ((seq_len(n) * 37) %% 101) / 100
  coreFraction <- rep(rep(c(0.5, 1), length.out = nPix), times = nYr)
  rate <- ifelse(x > 0.5, 0.4, 0.1)
  data.table::data.table(pixelID = rep(seq_len(nPix), times = nYr), year = rep(2001:(2000 + nYr), each = nPix),
                         x = x, z = ((seq_len(n) * 53) %% 97) / 97,
                         ignitions = as.integer(stats::rpois(n, rate * coreFraction)),
                         coreFraction = coreFraction, rate = rate)
}

fitCoreFraction <- function(dat, type = "ignition") {
  fitCache()
  set.seed(42)
  fitQuietly(data.table::copy(dat[, setdiff(names(dat), "rate"), with = FALSE]),
             dig = paste("coreFraction", type), type = type, nFolds = 5)
}

test_that("the ignition fit predicts the core's ignitions when pixels are partly buffer", {
  dat <- makeCoreFractionCovariates()
  m <- fitCoreFraction(dat)
  pred <- predict(m$Fold1, as.data.frame(dat)[, c("x", "z")])
  ## the prediction is of a whole pixel; the core's share of it is compared with the core's ignitions
  expect_equal(sum(pred * dat$coreFraction), sum(dat$ignitions), tolerance = 0.06)
  ## and it learns the rate of the whole pixel, not the rate diluted by the buffer
  expect_equal(mean(pred[dat$x > 0.5]), 0.4, tolerance = 0.1)
  expect_equal(mean(pred[dat$x <= 0.5]), 0.1, tolerance = 0.15)
  ## the out-of-year totals are of the core's ignitions too
  expect_equal(sum(m$yearTotals$predicted), sum(m$yearTotals$observed), tolerance = 0.08)
  expect_identical(sum(m$yearTotals$observed), sum(as.numeric(dat$ignitions)))
})

test_that("the weighting is what recovers the core's ignitions: the same rows without it are about 25% low", {
  dat <- makeCoreFractionCovariates()
  unweighted <- fitCoreFraction(dat[, setdiff(names(dat), "coreFraction"), with = FALSE])
  pred <- predict(unweighted$Fold1, as.data.frame(dat)[, c("x", "z")])
  expect_lt(sum(pred * dat$coreFraction) / sum(dat$ignitions), 0.85)
})

test_that("coreFraction is never a feature, and rows with no core are dropped", {
  dat <- makeCoreFractionCovariates()
  data.table::set(dat, which(dat$pixelID %% 10 == 0), "coreFraction", 0)
  data.table::set(dat, which(dat$coreFraction == 0), "ignitions", 0L)
  m <- fitCoreFraction(dat)
  for (fold in paste0("Fold", 1:5)) expect_setequal(stats::variable.names(m[[fold]]), c("x", "z"))
  ## the folds hold the integer response of the rows that remain
  nHeld <- vapply(m$rocs, function(r) length(r$response), numeric(1))
  expect_identical(sum(nHeld), as.numeric(sum(dat$coreFraction > 0)))
  expect_true(all(unlist(lapply(m$rocs, function(r) r$response)) %in% 0:1))
  expect_identical(sum(m$yearTotals$observed), sum(as.numeric(dat$ignitions)))

  esc <- data.table::copy(dat)
  data.table::set(esc, NULL, "escapes", as.integer(esc$ignitions > 0 & esc$z > 0.5))
  data.table::set(esc, NULL, "ignitions", pmax(esc$ignitions, esc$escapes))
  me <- fitCoreFraction(esc, "escape")
  for (fold in paste0("Fold", 1:5)) expect_setequal(stats::variable.names(me[[fold]]), c("x", "z"))
})

test_that("without coreFraction the fit is the one it was: every weight is 1", {
  dat <- makeCoreFractionCovariates()
  without <- dat[, setdiff(names(dat), "coreFraction"), with = FALSE]
  ones <- data.table::copy(dat)
  data.table::set(ones, NULL, "coreFraction", 1)
  mWithout <- fitCoreFraction(without)
  mOnes <- fitCoreFraction(ones)
  X <- as.data.frame(dat)[, c("x", "z")]
  expect_equal(predict(mWithout$Fold1, X), predict(mOnes$Fold1, X))
  expect_equal(mWithout$yearTotals, mOnes$yearTotals)
  expect_identical(exposureOf(without, "ignition"), rep(1, nrow(without)))
  expect_identical(exposureOf(dat, "escape"), rep(1, nrow(dat)))
  expect_identical(exposureOf(dat, "ignition"), dat$coreFraction)
})

test_that("coreFraction is not a feature of the stored-fit check or the plots", {
  expect_identical(xgbFeatureNames(c("pixelID", "year", "ignitions", "x", "coreFraction", "escapes", "z")), c("x", "z"))
  dat <- makeCoreFractionCovariates(nPix = 200, nYr = 2)
  data.table::set(dat, NULL, "u", dat$x * dat$z)
  models <- list(Fold1 = stats::lm(ignitions ~ x + z + u, data = dat))
  plots <- setupPlots(models, dat = dat[, c("x", "z", "u", "ignitions", "coreFraction", "pixelID", "year")], igOrEsc = "ignition")
  expect_setequal(c(as.character(plots$Climate$data$variable), as.character(plots$Fuel$data$variable)), c("x", "z", "u"))
})
