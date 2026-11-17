## The xgboost settings are module parameters (`xgb*`); their defaults are the values the fits were hard-coded with.

## A toy where more rounds overfit: a weak signal in x1 under Poisson noise, 4 years x 300 pixels
makeNoisyCovariates <- function(seed = 7, nPix = 300, nYr = 4) {
  withr::local_seed(seed)
  n <- nPix * nYr
  dt <- data.table::data.table(pixelID = rep(seq_len(nPix), nYr), year = rep(2001:(2000 + nYr), each = nPix),
                               x1 = stats::rnorm(n), x2 = stats::rnorm(n), x3 = stats::rnorm(n))
  data.table::set(dt, NULL, "ignitions", as.integer(stats::rpois(n, exp(-1.2 + 0.35 * dt$x1))))
  dt
}

## A toy with a planted hotspot: a four-way interaction of binary-ish covariates holds 1/16 of the pixel-years
makeHotspotCovariates <- function(seed = 11, nPix = 1500, nYr = 4) {
  withr::local_seed(seed)
  n <- nPix * nYr
  dt <- data.table::data.table(pixelID = rep(seq_len(nPix), nYr), year = rep(2001:(2000 + nYr), each = nPix),
                               a = stats::runif(n), b = stats::runif(n), c = stats::runif(n), d = stats::runif(n))
  hot <- dt$a > 0.5 & dt$b > 0.5 & dt$c > 0.5 & dt$d > 0.5
  data.table::set(dt, NULL, "ignitions", as.integer(stats::rpois(n, ifelse(hot, 0.8, 0.01))))
  dt
}

test_that("the xgb* parameter defaults are xgbSettings()'s, which are the old hard-coded values", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  p <- md$parameters
  default <- function(name) p$default[[match(name, p$paramName)]]
  par <- list(xgbMaxDepth = default("xgbMaxDepth"), xgbLearningRate = default("xgbLearningRate"),
              xgbNrounds = default("xgbNrounds"), xgbRegLambda = default("xgbRegLambda"),
              xgbEarlyStopping = default("xgbEarlyStopping"))
  expect_identical(xgbSettingsOfParams(par), xgbSettings())
  expect_identical(xgbSettings(), list(nrounds = 100, max_depth = 2, reg_lambda = 0.5, learning_rate = 0.15))
  expect_null(xgbSettings(earlyStopping = NA)$early_stopping_rounds)
  expect_null(xgbSettings(earlyStopping = 0)$early_stopping_rounds)
  expect_identical(xgbSettings(earlyStopping = 20)$early_stopping_rounds, 20L)
})

test_that("with the default settings a fold's model is the one the hard-coded call made", {
  fitCache()
  dat <- scaleCovariates(makeIgnitionCovariates())
  x <- as.data.frame(dat)[, c("CMDsm", "lightning", "Pice_mar", "youngAge")]
  evalSet <- seq(1, 400, by = 5)
  w <- rep(c(1, 0.5), 200)
  fit <- function(...) fitFoldXgboost(x, dat$ignitions / w, evalSet, "defaultsTest", "defaultsTest", weight = w, ...)
  utils::capture.output(m <- fit())
  utils::capture.output(old <- xgboost(x = x, y = dat$ignitions / w, weights = w, objective = "count:poisson",
                                       nthread = 1, eval_set = evalSet, monitor_training = TRUE, verbosity = 0,
                                       nrounds = 100, max_depth = 2, reg_lambda = 0.5, learning_rate = 0.15))
  expect_identical(predict(m, x), predict(old, x))
  expect_equal(attr(m, "params")[c("max_depth", "learning_rate", "reg_lambda")],
               list(max_depth = 2, learning_rate = 0.15, reg_lambda = 0.5))
  expect_identical(nrow(attr(m, "evaluation_log")), 100L)
})

test_that("the settings reach the model: depth, learning rate, rounds and lambda", {
  fitCache()
  dat <- scaleCovariates(makeIgnitionCovariates())
  set.seed(42)
  m <- fitQuietly(data.table::copy(dat), dig = "test", type = "ignition", nFolds = 3,
                  xgb = xgbSettings(maxDepth = 4, learningRate = 0.1, nrounds = 30, regLambda = 2))
  for (fold in paste0("Fold", 1:3)) {
    expect_equal(attr(m[[fold]], "params")[c("max_depth", "learning_rate", "reg_lambda")],
                 list(max_depth = 4, learning_rate = 0.1, reg_lambda = 2))
    expect_identical(nrow(attr(m[[fold]], "evaluation_log")), 30L)
  }
  expect_identical(m$concentration$rounds, rep(30, 3))
})

test_that("the module parameters reach the fits of both processes", {
  sim <- runModuleMsgs(localModulePaths(), new.env(),
                       params = list(xgbMaxDepth = 3, xgbLearningRate = 0.2, xgbNrounds = 20, xgbRegLambda = 1))
  for (nm in c("fireSense_IgnitionFitted", "fireSense_EscapeFitted")) {
    m <- sim[[nm]]$modelList$model$Fold1
    expect_equal(attr(m, "params")[c("max_depth", "learning_rate", "reg_lambda")],
                 list(max_depth = 3, learning_rate = 0.2, reg_lambda = 1))
    expect_identical(nrow(attr(m, "evaluation_log")), 20L)
  }
})

test_that("early stopping keeps fewer rounds where more rounds overfit, and the best one", {
  fitCache()
  dat <- makeNoisyCovariates()
  set.seed(42)
  cap <- xgbSettings(maxDepth = 4, learningRate = 0.3, nrounds = 200)
  full <- fitQuietly(data.table::copy(dat), dig = "noisy", type = "ignition", nFolds = 4, xgb = cap)
  early <- fitQuietly(data.table::copy(dat), dig = "noisy", type = "ignition", nFolds = 4,
                      xgb = modifyList(cap, list(early_stopping_rounds = 10L)))
  expect_identical(full$concentration$rounds, rep(200, 4))
  expect_true(all(early$concentration$rounds >= 1 & early$concentration$rounds < 100))
  ## the kept round is the one with the lowest held-out loss: later rounds only overfit
  for (fold in paste0("Fold", 1:4)) {
    heldOut <- attr(full[[fold]], "evaluation_log")$eval_poisson_nloglik
    keptLoss <- attr(early[[fold]], "evaluation_log")$eval_poisson_nloglik
    expect_lt(min(keptLoss), heldOut[200])
    expect_equal(min(keptLoss), min(heldOut[seq_along(keptLoss)]))
  }
  ## and it is a better model out of year
  expect_lt(early$concentration$deviance, full$concentration$deviance)
})

test_that("deeper trees concentrate the predictions on a planted hotspot: the lift rises", {
  fitCache()
  dat <- makeHotspotCovariates()
  set.seed(42)
  shallow <- fitQuietly(data.table::copy(dat), dig = "hot", type = "ignition", nFolds = 4)   # the defaults
  deep <- fitQuietly(data.table::copy(dat), dig = "hot", type = "ignition", nFolds = 4,
                     xgb = xgbSettings(maxDepth = 6, learningRate = 0.1, nrounds = 300, earlyStopping = 20))
  top5 <- function(m) m$concentration$lift[m$concentration$lift$top == 0.05, ]
  ## the hotspot holds 1/16 = 6.25% of the pixel-years and about 86% of the events
  expect_gt(top5(deep)$observed, 0.55)
  expect_gt(top5(deep)$predicted, top5(shallow)$predicted + 0.05)
  expect_gt(deep$concentration$gini[["predicted"]], shallow$concentration$gini[["predicted"]])
  expect_lt(deep$concentration$deviance, shallow$concentration$deviance)
})

test_that("the refit with the chosen rounds is the xgb.train() model at its best round", {
  fitCache()
  dat <- makeNoisyCovariates()
  x <- as.data.frame(dat)[, c("x1", "x2", "x3")]
  w <- rep(c(1, 0.5, 0.8), length.out = nrow(x))
  y <- dat$ignitions / w
  evalSet <- seq(1, nrow(x), by = 5)
  xgb <- xgbSettings(maxDepth = 4, learningRate = 0.3, nrounds = 200, earlyStopping = 10)
  utils::capture.output(m <- fitFoldXgboost(x, y, evalSet, "sameModel", "sameModel", weight = w, xgb = xgb))
  ## the same training, in one xgb.train() call with early stopping, and its model at the best round
  dm <- xgb.DMatrix(as.matrix(x), label = y, weight = w)
  trained <- xgb.train(params = list(objective = "count:poisson", max_depth = 4, learning_rate = 0.3, reg_lambda = 0.5,
                                     nthread = 1, verbosity = 0),
                       data = xgb.slice.DMatrix(dm, setdiff(seq_len(nrow(x)), evalSet)), nrounds = 200,
                       evals = list(eval = xgb.slice.DMatrix(dm, as.integer(evalSet))), early_stopping_rounds = 10, verbose = 0)
  best <- which.min(attributes(trained)$evaluation_log$eval_poisson_nloglik)
  expect_identical(nrow(attr(m, "evaluation_log")), as.integer(best))
  expect_lt(best, 200)
  atBest <- predict(trained, as.matrix(x), iterationrange = c(1, best))
  expect_equal(unname(predict(m, x)), unname(atBest), tolerance = 1e-6)
  expect_equal(unname(predict(trained, as.matrix(x))), unname(atBest))             # the default prediction is at the best round
  ## the comparison is not trivial: the rounds after the best one predict differently
  expect_false(isTRUE(all.equal(unname(predict(trained, as.matrix(x), iterationrange = c(1, nrow(attributes(trained)$evaluation_log)))),
                                unname(atBest), tolerance = 1e-6)))
})
