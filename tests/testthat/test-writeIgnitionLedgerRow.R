## The shared, geo-keyed ignition/escape fit ledger (mirrors fireSense_SpreadFit's).
##
## `ignitionFitGoogleDriveFolder = NULL` keeps `CacheGeo()` purely local (see its source: the
## Google Drive branch only runs when `isTRUE(useCloud) || !is.null(cloudFolderID)`), so these
## tests read and write a real ledger file on disk without any Google Drive access at all.

studyAreaToy <- function() terra::vect(terra::ext(0, 2500, 0, 2500), crs = "EPSG:3005")

## Like runModule() in test-events.R, but paths are supplied (not a fresh tempdir each call), so
## a ledger file written by one call can be read back by the next.
runIgnitionModule <- function(paths, params = list(), objects = list(), end = 1, msgs = new.env()) {
  defaults <- list(
    fireSense_ignitionCovariates = makeIgnitionCovariates(),
    fireSense_escapeCovariates   = makeIgnitionCovariates(escapes = TRUE),
    ignitionFitRTM               = makeIgnitionFitRTM()
  )
  defaults[names(objects)] <- objects
  objects <- defaults
  sim <- NULL
  msgs$text <- character()
  utils::capture.output(suppressWarnings(withCallingHandlers({
    sim <- SpaDES.core::simInit(times = list(start = 1, end = end), modules = moduleName,
                                paths = paths, objects = objects,
                                params = stats::setNames(list(testParams(params)), moduleName))
    sim <- SpaDES.core::spades(sim)
  }, message = function(m) { msgs$text <- c(msgs$text, conditionMessage(m)); invokeRestart("muffleMessage") })))
  sim
}

localLedgerPaths <- function() {
  root <- withr::local_tempdir(.local_envir = parent.frame())
  paths <- list(cachePath = file.path(root, "cache"), inputPath = file.path(root, "inputs"),
                modulePath = modulePath, outputPath = file.path(root, "outputs"))
  for (p in paths[c("cachePath", "inputPath", "outputPath")])
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  paths
}

test_that("without studyArea, nothing changes: no ledger is read or written", {
  paths <- localLedgerPaths()
  sim <- runIgnitionModule(paths, params = list(.plots = NA), objects = list(.ELFind = "6.1.1"))
  expect_null(sim$ignitionFitPreRun)
  expect_length(list.files(paths$inputPath, pattern = "fireSenseIgnitionParams"), 0L)
  expect_identical(names(sim$fireSense_IgnitionFitted$modelList$model), c(paste0("Fold", 1:5), "rocs", "yearTotals"))
})

## the model tag is fireSenseUtils' to set; only the years come from this module
ledgerName <- paste0("fireSenseIgnitionParams_2001-2004", fireSenseUtils::ignitionFitFileTag, ".rds")

test_that("a fit writes one ledger row per polygon, with rocs stripped", {
  paths <- localLedgerPaths()
  sa <- studyAreaToy()
  sim <- runIgnitionModule(paths, params = list(.plots = NA, ignitionFitGoogleDriveFolder = NULL),
                           objects = list(.ELFind = "6.1.1", studyArea = sa))

  files <- list.files(paths$inputPath, pattern = paste0("^fireSenseIgnitionParams_.*", fireSenseUtils::ignitionFitFileTag, "\\.rds$"))
  expect_length(files, 1L)
  expect_identical(files, ledgerName)   # from makeIgnitionCovariates()'s years

  expect_identical(as.character(sim$ignitionFitPreRun$polygonID), "6.1.1")
  ledgerFit <- sim$ignitionFitPreRun$fireSense_IgnitionFitted[[1]]
  ledgerEsc <- sim$ignitionFitPreRun$fireSense_EscapeFitted[[1]]
  expect_null(ledgerFit$modelList$model$rocs)                     # stripped before writing
  expect_null(ledgerEsc$modelList$model$rocs)
  expect_identical(names(ledgerFit$modelList$model), c(paste0("Fold", 1:5), "yearTotals"))
  ## the year effect's sd travels with the fit, for fireSense_ignitionPredict to pass to fireSense_burn
  expect_identical(ledgerFit$modelList$yearSD, sim$fireSense_IgnitionFitted$modelList$yearSD)
  expect_identical(ledgerEsc$modelList$yearSD, sim$fireSense_EscapeFitted$modelList$yearSD)

  ## the live sim objects (not round-tripped) still have rocs: only the ledger copy is stripped
  expect_false(is.null(sim$fireSense_IgnitionFitted$modelList$model$rocs))
  ## everything else about the fit is unchanged by stripping
  expect_identical(ledgerFit$modelList$model$Fold1, sim$fireSense_IgnitionFitted$modelList$model$Fold1)
  expect_identical(ledgerFit$scaleData, sim$fireSense_IgnitionFitted$scaleData)
})

test_that("a named ignitionFitFilename (an experiment's) receives the fit, and the row records the escape definition", {
  ## fitFileSuffix is gone: an experiment names its ledger file in full. The escape definition is recorded so
  ## fireSenseUtils::latestIgnitionFits() uses only fits made under the run's.
  paths <- localLedgerPaths()
  runIgnitionModule(paths, params = list(.plots = NA, ignitionFitGoogleDriveFolder = NULL,
                                         ignitionFitFilename = "fireSenseIgnitionParams_myExperiment.rds"),
                    objects = list(.ELFind = "6.1.1", studyArea = studyAreaToy()))
  files <- list.files(paths$inputPath, pattern = "^fireSenseIgnitionParams_.*\\.rds$")
  expect_identical(files, "fireSenseIgnitionParams_myExperiment.rds")
  row <- readRDS(file.path(paths$inputPath, files))
  expect_equal(as.data.frame(row)[, names(fireSenseUtils::escapeDefinition())],
               as.data.frame(fireSenseUtils::escapeDefinition()), ignore_attr = TRUE)
  expect_identical(ignitionLedgerWriteFile("mine.rds", 2001:2004), "mine.rds")  # a named file is used as given
})

test_that("a stored ledger fit is used instead of refitting, and predict() still works on it", {
  paths <- localLedgerPaths()
  sa <- studyAreaToy()
  commonParams <- list(.plots = NA, ignitionFitGoogleDriveFolder = NULL)

  sim1 <- runIgnitionModule(paths, params = commonParams, objects = list(.ELFind = "6.1.1", studyArea = sa))
  expect_true(file.exists(file.path(paths$inputPath, ledgerName)))

  ## second run, same ledger location: the fit must not run again. The converted package's name
  ## has dots where the module name has underscores (R package names cannot have underscores).
  local_mocked_bindings(
    buildModelsFitModels = function(...) stop("must not fit: the ledger already has a fit for this polygon"),
    .package = gsub("_", ".", moduleName, fixed = TRUE))
  sim2 <- runIgnitionModule(paths, params = commonParams, objects = list(.ELFind = "6.1.1", studyArea = sa))

  ## equal, not identical: each is an independent readRDS() of xgboost Booster objects, whose
  ## external pointers differ between deserializations even when the model content is the same.
  expect_equal(sim2$fireSense_IgnitionFitted, sim1$ignitionFitPreRun$fireSense_IgnitionFitted[[1]])
  expect_equal(sim2$fireSense_EscapeFitted, sim1$ignitionFitPreRun$fireSense_EscapeFitted[[1]])
  expect_identical(as.character(sim2$ignitionFitPreRun$polygonID), "6.1.1")

  ## predict() still works on the ledger's stripped, round-tripped object
  sc <- sim2$fireSense_IgnitionFitted$scaleData
  covs <- c("CMDsm", "lightning", "Pice_mar", "youngAge")
  raw <- makeIgnitionCovariates()
  newdata <- as.data.frame(raw)[, covs]
  for (col in covs) newdata[[col]] <- (newdata[[col]] - sc$`scaled:center`[[col]]) / sc$`scaled:scale`[[col]]
  pred <- predict(sim2$fireSense_IgnitionFitted$modelList$model$Fold1, newdata)
  expect_equal(as.numeric(tapply(pred, raw$ignitions, mean)), c(0, 1, 2), tolerance = 0.1)
})

test_that("another polygon's fit is appended to the ledger and the first polygon's row is kept", {
  paths <- localLedgerPaths()
  commonParams <- list(.plots = NA, ignitionFitGoogleDriveFolder = NULL)
  next_to <- terra::vect(terra::ext(10000, 12500, 0, 2500), crs = "EPSG:3005")
  runIgnitionModule(paths, params = commonParams, objects = list(.ELFind = "6.1.1", studyArea = studyAreaToy()))
  sim2 <- runIgnitionModule(paths, params = commonParams, objects = list(.ELFind = "6.1.2", studyArea = next_to))
  ledger <- readRDS(file.path(paths$inputPath, ledgerName))
  expect_setequal(as.character(ledger$polygonID), c("6.1.1", "6.1.2"))
  expect_equal(ledger$fireSense_IgnitionFitted[[which(ledger$polygonID == "6.1.2")]],
               sim2$ignitionFitPreRun$fireSense_IgnitionFitted[[which(sim2$ignitionFitPreRun$polygonID == "6.1.2")]])
})

test_that("doFit = TRUE fits again even when the ledger already has a row", {
  paths <- localLedgerPaths()
  sa <- studyAreaToy()
  commonParams <- list(.plots = NA, ignitionFitGoogleDriveFolder = NULL)

  runIgnitionModule(paths, params = commonParams, objects = list(.ELFind = "6.1.1", studyArea = sa))

  local_mocked_bindings(
    hasPreRunIgnitionFitForThisPolygon = function(sim) stop("must not be consulted when doFit is TRUE"),
    .package = gsub("_", ".", moduleName, fixed = TRUE))
  sim2 <- runIgnitionModule(paths, params = c(commonParams, list(doFit = TRUE)),
                            objects = list(.ELFind = "6.1.1", studyArea = sa))
  expect_identical(names(sim2$fireSense_IgnitionFitted$modelList$model), c(paste0("Fold", 1:5), "rocs", "yearTotals"))
  ## the refit replaces the polygon's row; it does not add a second one
  ledger <- readRDS(file.path(paths$inputPath, ledgerName))
  expect_identical(as.character(ledger$polygonID), "6.1.1")
})

## Does a run fit? A run that reuses the ledger's fit says so; a run that fits writes its row. A ledger is
## seeded first when wanted. (buildModelsFitModels() cannot be counted by mocking: the module calls it from
## its own environment.)
fitsInRun <- function(paths, params, ledger, seed = list()) {
  commonParams <- list(.plots = NA, ignitionFitGoogleDriveFolder = NULL)
  objs <- list(.ELFind = "6.1.1", studyArea = studyAreaToy())
  if (ledger) runIgnitionModule(paths, params = c(commonParams, seed), objects = objs)
  msgs <- new.env()
  out <- try(runIgnitionModule(paths, params = c(commonParams, params), objects = objs, msgs = msgs), silent = TRUE)
  failed <- inherits(out, "try-error")
  list(fits = if (failed) 0L else as.integer(!any(grepl("using those instead of fitting", msgs$text))),
       error = if (failed) as.character(out) else NULL, msgs = msgs$text)
}

test_that("doFit FALSE (the default), NA and TRUE, with and without a ledger fit", {
  paths <- localLedgerPaths()
  ## no ledger fit
  r <- fitsInRun(localLedgerPaths(), list(doFit = NULL), ledger = FALSE)                  # default: error
  expect_identical(r$fits, 0L)
  expect_match(r$error, "polygon 6\\.1\\.1 in the ledger.*the default.*doFit = NA")
  expect_gt(fitsInRun(localLedgerPaths(), list(doFit = NA), ledger = FALSE)$fits, 0L)     # fit
  expect_gt(fitsInRun(localLedgerPaths(), list(doFit = TRUE), ledger = FALSE)$fits, 0L)
  r <- fitsInRun(localLedgerPaths(), list(doFit = FALSE), ledger = FALSE)                  # never fit: error
  expect_identical(r$fits, 0L)
  expect_match(r$error, "polygon 6\\.1\\.1.*doFit")
  ## a ledger fit
  expect_identical(fitsInRun(localLedgerPaths(), list(doFit = NULL), ledger = TRUE)$fits, 0L)  # default: reuse
  expect_identical(fitsInRun(localLedgerPaths(), list(doFit = NA), ledger = TRUE)$fits, 0L)
  expect_identical(fitsInRun(localLedgerPaths(), list(doFit = FALSE), ledger = TRUE)$fits, 0L)
  expect_gt(fitsInRun(localLedgerPaths(), list(doFit = TRUE), ledger = TRUE)$fits, 0L)    # replace it
})

test_that("refitExisting is a deprecated alias of doFit = TRUE, with a message naming doFit", {
  r <- fitsInRun(localLedgerPaths(), list(refitExisting = TRUE), ledger = TRUE)
  expect_gt(r$fits, 0L)
  expect_length(grep("refitExisting.*deprecated.*doFit", r$msgs), 1L)
  ## FALSE is the old default: it does nothing, and says nothing about doFit = TRUE
  expect_identical(fitsInRun(localLedgerPaths(), list(refitExisting = FALSE), ledger = TRUE)$fits, 0L)
  ## doFit FALSE, the default, yields to the alias (it cannot be told apart from an unset doFit); NA does not
  expect_gt(fitsInRun(localLedgerPaths(), list(doFit = NULL, refitExisting = TRUE), ledger = TRUE)$fits, 0L)
  expect_identical(fitsInRun(localLedgerPaths(), list(doFit = NA, refitExisting = FALSE), ledger = TRUE)$fits, 0L)
  ## no deprecated parameter set: no deprecation message
  expect_length(grep("deprecated", fitsInRun(localLedgerPaths(), list(doFit = NA), ledger = TRUE)$msgs), 0L)
})

test_that("the default doFit never fits: no ledger row stops the run before any xgboost fit", {
  ## a user who adds fireSense only to predict must never start a fit by accident.
  ## Each fold's xgboost fit is cached (fitFoldXgboost(), function name "xgboost_<type>_<fold>"), so the
  ## run's cache shows whether any fit ran; with doFit = NA the same count is positive.
  xgboostFits <- function(params) {
    paths <- localLedgerPaths()
    r <- fitsInRun(paths, params, ledger = FALSE)
    sc <- reproducible::showCache(paths$cachePath, verbose = -2)
    r$xgboost <- length(unique(sc$cacheId[grepl("^xgboost", sc$tagValue)]))
    r
  }
  r <- xgboostFits(list(doFit = NULL))                     # the module's default
  expect_match(r$error, "no usable stored fit for polygon 6\\.1\\.1.*the default")
  expect_identical(r$xgboost, 0L)
  expect_gt(xgboostFits(list(doFit = NA))$xgboost, 0L)
})

test_that("writeIgnitionLedgerRow stops on an invalid .ELFind, as fireSense_SpreadFit does", {
  paths <- localLedgerPaths()
  sa <- studyAreaToy()
  expect_error(
    runIgnitionModule(paths, params = list(.plots = NA, ignitionFitGoogleDriveFolder = NULL),
                      objects = list(.ELFind = character(0), studyArea = sa)),
    "sim\\$.ELFind.*must be a single non-empty character"
  )
})

test_that("stripFittedForLedger drops rocs, shrinks the object, and predict() is unaffected", {
  paths <- localLedgerPaths()
  sim <- runIgnitionModule(paths, params = list(.plots = NA), objects = list(.ELFind = "6.1.1"))
  fitted <- sim$fireSense_IgnitionFitted
  expect_false(is.null(fitted$modelList$model$rocs))

  stripped <- stripFittedForLedger(fitted)
  expect_null(stripped$modelList$model$rocs)
  expect_lt(as.numeric(object.size(stripped)), as.numeric(object.size(fitted)))
  expect_identical(stripped$modelList$model$Fold1, fitted$modelList$model$Fold1)

  raw <- makeIgnitionCovariates()
  newdata <- as.data.frame(raw)[, c("CMDsm", "lightning", "Pice_mar", "youngAge")]
  expect_identical(predict(stripped$modelList$model$Fold1, newdata),
                   predict(fitted$modelList$model$Fold1, newdata))

  expect_null(stripFittedForLedger(NULL))
})

## A covariate set whose `Pice_mar` column has been renamed, as when land-cover groups are merged
## under a new name: the stored fit's features are no longer all among the covariates.
renamedCovariates <- function(escapes = FALSE) {
  dt <- makeIgnitionCovariates(escapes = escapes)
  data.table::setnames(dt, "Pice_mar", "Pice_new")
  dt
}

test_that("storedFitIsCompatible: all stored features present -> TRUE; one absent -> FALSE, named", {
  paths <- localLedgerPaths()
  sim <- runIgnitionModule(paths, params = list(.plots = NA), objects = list(.ELFind = "6.1.1"))
  for (type in c("ignition", "escape")) {
    fitted <- if (type == "ignition") sim$fireSense_IgnitionFitted else sim$fireSense_EscapeFitted
    covs <- makeIgnitionCovariates(escapes = type == "escape")
    expect_true(storedFitIsCompatible(fitted, covs, "6.1.1", type))
    ## a current column the model never saw does not stop reuse
    covs2 <- data.table::copy(covs)
    data.table::set(covs2, NULL, "extraCol", 1)
    expect_true(storedFitIsCompatible(fitted, covs2, "6.1.1", type))
    expect_message(
      ok <- storedFitIsCompatible(fitted, renamedCovariates(type == "escape"), "6.1.1", type),
      "6.1.1.*Pice_mar.*Pice_new")
    expect_false(ok)
  }
  expect_false(storedFitIsCompatible(NULL, makeIgnitionCovariates(), "6.1.1", "ignition"))
})

test_that("a stored ledger fit whose features are missing from the covariates is refitted and replaced", {
  paths <- localLedgerPaths()
  sa <- studyAreaToy()
  commonParams <- list(.plots = NA, ignitionFitGoogleDriveFolder = NULL)
  runIgnitionModule(paths, params = commonParams, objects = list(.ELFind = "6.1.1", studyArea = sa))

  sim2 <- runIgnitionModule(paths, params = commonParams,
                            objects = list(.ELFind = "6.1.1", studyArea = sa,
                                           fireSense_ignitionCovariates = renamedCovariates(),
                                           fireSense_escapeCovariates = renamedCovariates(escapes = TRUE)))
  for (fitted in list(sim2$fireSense_IgnitionFitted, sim2$fireSense_EscapeFitted))
    expect_setequal(stats::variable.names(fitted$modelList$model$Fold1),
                    c("CMDsm", "lightning", "Pice_new", "youngAge"))
  ## the new fit replaced the polygon's row
  ledger <- readRDS(file.path(paths$inputPath, ledgerName))
  expect_identical(as.character(ledger$polygonID), "6.1.1")
  expect_true("Pice_new" %in% stats::variable.names(ledger$fireSense_IgnitionFitted[[1]]$modelList$model$Fold1))
})

## Real small boosters of a given objective, in the shape of a stored ledger fit.
storedFitWith <- function(objective, folds = 2) {
  dat <- scaleCovariates(makeIgnitionCovariates())
  x <- dat[, c("CMDsm", "youngAge", "Pice_mar", "lightning")]
  model <- lapply(seq_len(folds), function(i)
    xgboost::xgboost(x = x, y = dat$ignitions, objective = objective, nrounds = 5, max_depth = 2,
                     nthread = 1, verbosity = 0))
  names(model) <- paste0("Fold", seq_len(folds))
  list(modelList = list(model = model))
}

test_that("a stored fit made with another objective (reg:tweedie) is not reused, count:poisson is", {
  covs <- makeIgnitionCovariates()
  expect_true(storedFitIsCompatible(storedFitWith("count:poisson"), covs, "6.1.1", "ignition"))
  expect_message(ok <- storedFitIsCompatible(storedFitWith("reg:tweedie"), covs, "6.1.1", "ignition"),
                 "reg:tweedie.*count:poisson")
  expect_false(ok)
  ## one stray fold is enough
  mixed <- storedFitWith("count:poisson")
  mixed$modelList$model$Fold2 <- storedFitWith("reg:tweedie")$modelList$model$Fold1
  expect_false(suppressMessages(storedFitIsCompatible(mixed, covs, "6.1.1", "ignition")))
})

test_that("writeIgnitionLedgerRow refuses a model fitted with another objective", {
  sim <- list(.ELFind = "6.1.1", fireSense_IgnitionFitted = storedFitWith("count:poisson"),
              fireSense_EscapeFitted = storedFitWith("reg:tweedie"))
  expect_error(writeIgnitionLedgerRow(sim), "reg:tweedie.*count:poisson")
})

test_that("buildModel's cache key changes when the body of a function it calls changes", {
  key <- function(callees) reproducible:::.robustDigest(callees)
  base <- buildModelCallees()
  expect_identical(key(base), key(buildModelCallees()))
  edited <- base
  body(edited$runXGBOOST) <- quote(stop("a different fit"))
  expect_false(identical(key(base), key(edited)))
  edited <- base
  edited$ignitionObjectiveTxt <- "reg:tweedie"
  expect_false(identical(key(base), key(edited)))
})

## A run over several ELFs (fireSense_ELFs makes `.ELFind` the ELFs a studyAreaLarge covers): each needs
## its stored fit; a fit is of one ELF.
test_that("hasPreRunIgnitionFitForThisPolygon(): every ELF of .ELFind needs its row", {
  rows <- data.frame(polygonID = c("6.1.1", "6.1.2"))
  has <- function(id, r = rows) hasPreRunIgnitionFitForThisPolygon(list(ignitionFitPreRun = r, .ELFind = id))
  expect_true(has("6.1.1"))
  expect_true(has(c("6.1.2", "6.1.1")))
  expect_false(has(c("6.1.1", "6.1.3")))
  expect_false(has(NA_character_, data.frame(polygonID = NA)))
  expect_false(has(character(0)))
})

test_that("over several ELFs the default doFit uses their stored fits; a fit of several stops", {
  paths <- localLedgerPaths()
  commonParams <- list(.plots = NA, ignitionFitGoogleDriveFolder = NULL)
  next_to <- terra::vect(terra::ext(2500, 5000, 0, 2500), crs = "EPSG:3005")
  runIgnitionModule(paths, params = commonParams, objects = list(.ELFind = "6.1.1", studyArea = studyAreaToy()))
  runIgnitionModule(paths, params = commonParams, objects = list(.ELFind = "6.1.2", studyArea = next_to))
  both <- terra::vect(terra::ext(0, 5000, 0, 2500), crs = "EPSG:3005")
  objs <- list(.ELFind = c("6.1.1", "6.1.2"), studyArea = both)
  msgs <- new.env()
  sim <- runIgnitionModule(paths, params = c(commonParams, list(doFit = NULL)), objects = objs, msgs = msgs)
  expect_length(grep("using those instead of fitting", msgs$text), 1L)
  ## the pooled covariates of several ELFs are not compared with one ELF's stored model
  dropLightning <- function(d) data.table::set(d, j = "lightning", value = NULL)
  noLightning <- c(objs, list(fireSense_ignitionCovariates = dropLightning(makeIgnitionCovariates()),
                              fireSense_escapeCovariates = dropLightning(makeIgnitionCovariates(escapes = TRUE))))
  msgs <- new.env()
  runIgnitionModule(paths, params = c(commonParams, list(doFit = NULL)), objects = noLightning, msgs = msgs)
  expect_length(grep("using those instead of fitting", msgs$text), 1L)
  expect_error(runIgnitionModule(paths, params = c(commonParams, list(doFit = TRUE)), objects = objs),
               "a fit is of one ELF.*6\\.1\\.1, 6\\.1\\.2")
})

test_that("a ledger fit made before the concentration diagnostic binds with new rows and is still used", {
  paths <- localLedgerPaths()
  commonParams <- list(.plots = NA, ignitionFitGoogleDriveFolder = NULL)
  next_to <- terra::vect(terra::ext(10000, 12500, 0, 2500), crs = "EPSG:3005")
  runIgnitionModule(paths, params = commonParams, objects = list(.ELFind = "6.1.1", studyArea = studyAreaToy()))
  ## make the stored row an old one: no concentration, no chain
  file <- file.path(paths$inputPath, ledgerName)
  ledger <- readRDS(file)
  for (col in c("fireSense_IgnitionFitted", "fireSense_EscapeFitted")) {
    expect_false(is.null(ledger[[col]][[1]]$modelList$concentration))
    ledger[[col]][[1]]$modelList$concentration <- NULL
  }
  saveRDS(ledger, file)
  ## a new fit of another polygon is appended to it, and the old row is read back without a concentration
  sim2 <- runIgnitionModule(paths, params = commonParams, objects = list(.ELFind = "6.1.2", studyArea = next_to))
  ledger <- readRDS(file)
  expect_setequal(as.character(ledger$polygonID), c("6.1.1", "6.1.2"))
  expect_null(ledger$fireSense_EscapeFitted[[which(ledger$polygonID == "6.1.1")]]$modelList$concentration)
  expect_false(is.null(ledger$fireSense_EscapeFitted[[which(ledger$polygonID == "6.1.2")]]$modelList$concentration$chain))
  sim3 <- runIgnitionModule(paths, params = commonParams, objects = list(.ELFind = "6.1.1", studyArea = studyAreaToy()))
  expect_null(sim3$fireSense_IgnitionFitted$modelList$concentration)
  expect_false(is.null(sim3$fireSense_IgnitionFitted$modelList$yearEffects))
})

## fireSense_dataPrepFit leaves the ignition fit data out (and the objects are NULL here, as the module that
## declares them has not made them) when doFit is FALSE and skipUnusedFitPrep is TRUE.
noFitData <- list(fireSense_ignitionCovariates = NULL, fireSense_escapeCovariates = NULL, ignitionFitRTM = NULL)

test_that("a ledger fit is reused without the fit data or ignitionFitRTM when doFit is FALSE", {
  paths <- localLedgerPaths()
  sa <- studyAreaToy()
  commonParams <- list(.plots = NA, ignitionFitGoogleDriveFolder = NULL)
  runIgnitionModule(paths, params = commonParams, objects = list(.ELFind = "6.1.1", studyArea = sa))

  msgs <- new.env()
  ## a named ledger file: with "latest" the file is named for the covariates' years (next test)
  sim <- runIgnitionModule(paths, params = c(commonParams, list(doFit = FALSE, ignitionFitFilename = ledgerName)),
                           msgs = msgs,
                           objects = c(list(.ELFind = "6.1.1", studyArea = sa), noFitData))
  expect_true(any(grepl("using those instead of fitting", msgs$text)))
  expect_true(any(grepl("fireSense_ignitionPredict checks it", msgs$text)))
  expect_false(is.null(sim$fireSense_IgnitionFitted$modelList$model$Fold1))
  expect_false(is.null(sim$fireSense_EscapeFitted$modelList$model$Fold1))
})

test_that("a fit still needs the nonNAs of ignitionFitRTM", {
  rtm <- makeIgnitionFitRTM()
  attr(rtm, "nonNAs") <- NULL
  expect_error(runIgnitionModule(localLedgerPaths(), params = list(.plots = NA, doFit = TRUE),
                                 objects = list(ignitionFitRTM = rtm)),
               "nonNAs must be a non-empty/non-NULL numeric", fixed = TRUE)
  ## with a ledger fit and doFit = TRUE, too
  paths <- localLedgerPaths()
  objs <- list(.ELFind = "6.1.1", studyArea = studyAreaToy())
  commonParams <- list(.plots = NA, ignitionFitGoogleDriveFolder = NULL)
  runIgnitionModule(paths, params = commonParams, objects = objs)
  expect_error(runIgnitionModule(paths, params = c(commonParams, list(doFit = TRUE)),
                                 objects = c(objs, list(ignitionFitRTM = rtm))),
               "nonNAs must be a non-empty/non-NULL numeric", fixed = TRUE)
})

test_that("without fit data, ignitionFitFilename = \"latest\" reads the newest row of the polygon", {
  rows <- data.frame(polygonID = "6.1.1")
  called <- NULL
  local_mocked_bindings(
    latestIgnitionFits = function(cloudFolderID, destinationPath, polygonIDs, escapeDefinition) {
      called <<- list(cloudFolderID = cloudFolderID, polygonIDs = polygonIDs, escapeDefinition = escapeDefinition)
      rows
    }, .package = "fireSenseUtils")
  sim <- SpaDES.core::simInit(times = list(start = 1, end = 1), modules = moduleName, paths = localLedgerPaths(),
                              objects = c(list(.ELFind = "6.1.1", studyArea = studyAreaToy()), noFitData),
                              params = stats::setNames(list(testParams(list(.plots = NA, doFit = FALSE,
                                                              ignitionFitGoogleDriveFolder = "folderID"))), moduleName))
  sim <- suppressMessages(SpaDES.core::spades(sim, events = list(fireSense_ignitionFit = c("init", "checkData"))))
  expect_identical(sim$ignitionFitPreRun, rows)
  ## the run's escape definition, so a fit made under another is skipped
  expect_identical(called, list(cloudFolderID = "folderID", polygonIDs = "6.1.1",
                                escapeDefinition = fireSenseUtils::escapeDefinition()))
})

test_that("a stored fit is reused only by a run with the same year climate settings, in both directions", {
  esc <- list(escapeYearClimate = "CMDsm")
  hinge <- list(yearClimateShape = "logistic")
  wts <- list(yearClimateWeights = "climate")
  reuse <- function(seed, run) fitsInRun(localLedgerPaths(), c(list(doFit = NA), run), ledger = TRUE, seed = seed)
  expect_identical(reuse(list(), list())$fits, 0L)                 # the same settings
  expect_identical(reuse(esc, esc)$fits, 0L)
  expect_identical(reuse(c(esc, hinge, wts), c(esc, hinge, wts))$fits, 0L)
  for (case in list(list(list(), esc), list(esc, list()),                       # an escape term set, or dropped
                    list(list(), hinge), list(hinge, list()),
                    list(list(), wts), list(wts, list()),
                    list(esc, c(esc, hinge)), list(c(esc, wts), esc),
                    list(list(), list(ignitionYearClimate = NA_character_)),     # the ignition term, set or dropped
                    list(list(ignitionYearClimate = NA_character_), list()),
                    list(list(), list(ignitionYearClimate = "lightning")))) {
    r <- reuse(case[[1]], case[[2]])
    expect_gt(r$fits, 0L)
    expect_match(paste(r$msgs, collapse = "\n"), "Fitting again")
  }
})

test_that("storedYearClimateProblem compares the escape term, the shape and the weights", {
  term <- function(...) list(modelList = list(yearClimate = list(variable = "CMDsm", center = 0, scale = 1, beta = 1, norm = 1, ...)))
  none <- list(modelList = list())
  par <- list(ignitionYearClimate = "CMDsm", escapeYearClimate = NA_character_, yearClimateShape = "loglinear", yearClimateWeights = "equal")
  set <- function(...) modifyList(par, list(...))
  expect_null(storedYearClimateProblem(term(), "ignition", par))
  expect_null(storedYearClimateProblem(none, "escape", par))
  expect_match(storedYearClimateProblem(none, "ignition", par), "no year climate term")   # a fit from before the ignition term
  expect_null(storedYearClimateProblem(none, "ignition", set(ignitionYearClimate = NA_character_)))
  expect_match(storedYearClimateProblem(term(), "ignition", set(ignitionYearClimate = NA_character_)), "ignitionYearClimate. is NA")
  expect_match(storedYearClimateProblem(term(), "ignition", set(ignitionYearClimate = "lightning")), "term is of CMDsm")
  expect_null(storedYearClimateProblem(term(), "ignition", set(ignitionYearClimate = "CMD_sm")))
  expect_match(storedYearClimateProblem(none, "escape", set(escapeYearClimate = "CMDsm")), "no year climate term")
  expect_match(storedYearClimateProblem(term(), "escape", par), "escapeYearClimate. is NA")
  expect_null(storedYearClimateProblem(term(), "escape", set(escapeYearClimate = "CMD_sm")))   # underscores do not matter
  expect_match(storedYearClimateProblem(term(), "escape", set(escapeYearClimate = "lightning")), "term is of CMDsm")
  expect_match(storedYearClimateProblem(term(shape = "logistic"), "ignition", par), "logistic.*loglinear")
  expect_match(storedYearClimateProblem(term(), "ignition", set(yearClimateShape = "logistic")), "loglinear.*logistic")
  lg <- term(shape = "logistic", bounds = logisticYearBounds)
  expect_null(storedYearClimateProblem(lg, "ignition", set(yearClimateShape = "logistic")))
  ## a logistic term fitted under other bounds (cap, k, ratio) than this version's, or recording none, is not reused
  other <- term(shape = "logistic", bounds = replace(logisticYearBounds, "capFactor", 5))
  expect_match(storedYearClimateProblem(other, "ignition", set(yearClimateShape = "logistic")), "other bounds")
  expect_match(storedYearClimateProblem(term(shape = "logistic"), "ignition", set(yearClimateShape = "logistic")), "other bounds")
  expect_match(storedYearClimateProblem(term(weights = "climate"), "ignition", par), "climate.*equal")
  expect_match(storedYearClimateProblem(term(), "ignition", set(yearClimateWeights = "climate")), "equal.*climate")
  expect_null(storedYearClimateProblem(term(weights = "climate"), "ignition", set(yearClimateWeights = "climate")))
})
