## The shared, geo-keyed ignition/escape fit ledger (mirrors fireSense_SpreadFit's).
##
## `ignitionFitGoogleDriveFolder = NULL` keeps `CacheGeo()` purely local (see its source: the
## Google Drive branch only runs when `isTRUE(useCloud) || !is.null(cloudFolderID)`), so these
## tests read and write a real ledger file on disk without any Google Drive access at all.

studyAreaToy <- function() terra::vect(terra::ext(0, 2500, 0, 2500), crs = "EPSG:3005")

## Like runModule() in test-events.R, but paths are supplied (not a fresh tempdir each call), so
## a ledger file written by one call can be read back by the next.
runIgnitionModule <- function(paths, params = list(), objects = list(), end = 1) {
  defaults <- list(
    fireSense_ignitionCovariates = makeIgnitionCovariates(),
    fireSense_escapeCovariates   = makeIgnitionCovariates(escapes = TRUE),
    ignitionFitRTM               = makeIgnitionFitRTM()
  )
  defaults[names(objects)] <- objects
  objects <- defaults
  sim <- NULL
  utils::capture.output(suppressWarnings(suppressMessages({
    sim <- SpaDES.core::simInit(times = list(start = 1, end = end), modules = moduleName,
                                paths = paths, objects = objects,
                                params = stats::setNames(list(params), moduleName))
    sim <- SpaDES.core::spades(sim)
  })))
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
  expect_identical(names(sim$fireSense_IgnitionFitted$modelList$model), c(paste0("Fold", 1:5), "rocs"))
})

test_that("a fit writes one ledger row per polygon, with rocs stripped", {
  paths <- localLedgerPaths()
  sa <- studyAreaToy()
  sim <- runIgnitionModule(paths, params = list(.plots = NA, ignitionFitGoogleDriveFolder = NULL),
                           objects = list(.ELFind = "6.1.1", studyArea = sa))

  files <- list.files(paths$inputPath, pattern = "^fireSenseIgnitionParams_.*_xgboost\\.rds$")
  expect_length(files, 1L)
  expect_identical(files, "fireSenseIgnitionParams_2001-2004_xgboost.rds")   # from makeIgnitionCovariates()'s years

  expect_identical(as.character(sim$ignitionFitPreRun$polygonID), "6.1.1")
  ledgerFit <- sim$ignitionFitPreRun$fireSense_IgnitionFitted[[1]]
  ledgerEsc <- sim$ignitionFitPreRun$fireSense_EscapeFitted[[1]]
  expect_null(ledgerFit$modelList$model$rocs)                     # stripped before writing
  expect_null(ledgerEsc$modelList$model$rocs)
  expect_identical(names(ledgerFit$modelList$model), paste0("Fold", 1:5))

  ## the live sim objects (not round-tripped) still have rocs: only the ledger copy is stripped
  expect_false(is.null(sim$fireSense_IgnitionFitted$modelList$model$rocs))
  ## everything else about the fit is unchanged by stripping
  expect_identical(ledgerFit$modelList$model$Fold1, sim$fireSense_IgnitionFitted$modelList$model$Fold1)
  expect_identical(ledgerFit$scaleData, sim$fireSense_IgnitionFitted$scaleData)
})

test_that("a stored ledger fit is used instead of refitting, and predict() still works on it", {
  paths <- localLedgerPaths()
  sa <- studyAreaToy()
  commonParams <- list(.plots = NA, ignitionFitGoogleDriveFolder = NULL)

  sim1 <- runIgnitionModule(paths, params = commonParams, objects = list(.ELFind = "6.1.1", studyArea = sa))
  expect_true(file.exists(file.path(paths$inputPath, "fireSenseIgnitionParams_2001-2004_xgboost.rds")))

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
  ledger <- readRDS(file.path(paths$inputPath, "fireSenseIgnitionParams_2001-2004_xgboost.rds"))
  expect_setequal(as.character(ledger$polygonID), c("6.1.1", "6.1.2"))
  expect_equal(ledger$fireSense_IgnitionFitted[[which(ledger$polygonID == "6.1.2")]],
               sim2$ignitionFitPreRun$fireSense_IgnitionFitted[[which(sim2$ignitionFitPreRun$polygonID == "6.1.2")]])
})

test_that("refitExisting fits again even when the ledger already has a row", {
  paths <- localLedgerPaths()
  sa <- studyAreaToy()
  commonParams <- list(.plots = NA, ignitionFitGoogleDriveFolder = NULL)

  runIgnitionModule(paths, params = commonParams, objects = list(.ELFind = "6.1.1", studyArea = sa))

  local_mocked_bindings(
    hasPreRunIgnitionFitForThisPolygon = function(sim) stop("must not be consulted when refitExisting is TRUE"),
    .package = gsub("_", ".", moduleName, fixed = TRUE))
  sim2 <- runIgnitionModule(paths, params = c(commonParams, list(refitExisting = TRUE)),
                            objects = list(.ELFind = "6.1.1", studyArea = sa))
  expect_identical(names(sim2$fireSense_IgnitionFitted$modelList$model), c(paste0("Fold", 1:5), "rocs"))
  ## the refit replaces the polygon's row; it does not add a second one
  ledger <- readRDS(file.path(paths$inputPath, "fireSenseIgnitionParams_2001-2004_xgboost.rds"))
  expect_identical(as.character(ledger$polygonID), "6.1.1")
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
  ledger <- readRDS(file.path(paths$inputPath, "fireSenseIgnitionParams_2001-2004_xgboost.rds"))
  expect_identical(as.character(ledger$polygonID), "6.1.1")
  expect_true("Pice_new" %in% stats::variable.names(ledger$fireSense_IgnitionFitted[[1]]$modelList$model$Fold1))
})
