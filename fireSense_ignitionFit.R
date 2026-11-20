defineModule(sim, list(
  name = "fireSense_ignitionFit",
  description = paste("Fit statistical models that can be used to parameterize (calibrate)",
                      "the fire ignition component of landscape fire models (e.g. fireSense)."),
  keywords = c("fire frequency", "optimization", "additive property", "poisson",
               "negative binomial", "fireSense"),
  authors = c(
    person("Eliot", "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = "aut"),
    person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = "aut"),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = "ctb")
  ),
  childModules = character(),
  version = list(fireSense_ignitionFit = "1.2.0.9016"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_ignitionFit.Rmd"),
  loadOrder = list(after = c("fireSense_dataPrepFit", "fireSense_ELFs"),
                   before = "fireSense_dataPrepPredict"),
  reqdPkgs = list("data.table", "dplyr", "PredictiveEcology/SpaDES.core@development (>= 3.2.1.9037)",
                  "PredictiveEcology/fireSenseUtils@development (>= 0.2.3.9184)", # 9184: scaleToObservedLevel(), yearClimateMeans(covariates, variable) unweighted over core rows, yearClimateMultiplier() norm; 9182: yearClimateColumn(), yearClimateMeans(), yearClimateMultiplier(); 9179: latestIgnitionFits(escapeDefinition); 9177: storedFitProblem(), xgbFeatureNames(), foreignObjectives(); 9168: foldMeanPrediction(), escapeProbPerIgnition()
                  "ggplot2", "ggpubr", "googledrive", "magrittr", "purrr", "RColorBrewer",
                  "numDeriv", "parallel", "parallelly",
                  "PredictiveEcology/pemisc@development",
                  "PredictiveEcology/reproducible@development (>= 3.2.1.9067)", # CacheGeo: ledger with xgboost models; `le = le` no longer partially matches `ledger`
                  "RhpcBLASctl", "sf",
                  "caret", "pROC",
                  "PredictiveEcology/SHAPforxgboost (>= 0.1.3.9001)", "xgboost (>= 3.2.1.1)", "lightgbm", # 3.2.1.1: xgb.train(early_stopping_rounds), xgb.slice.DMatrix and evaluation_log used by earlyStoppedRounds() were checked on it only; install.packages('xgboost', repos = c('https://dmlc.r-universe.dev', 'https://cloud.r-project.org'))
                  "terra"),
  parameters = bindrows(
    defineParameter("crossValType", "character", c("time-ordered", "crossValidation"), NA, NA,
                    paste("Has no effect on the fit: `runXGBOOST()` always uses k-fold cross validation.",
                          "The first element is only used in the plot filename.")),
    defineParameter("rescaleVars", "logical", default = TRUE,
                    desc = paste("If `TRUE`, covariates are centred and scaled with `scale()` before fitting.",
                                 "The centring and scaling values are returned in the `scaleData` element of the outputs.")),
    defineParameter("whichProcessesToFit", "character", c("ignition", "escape"), NA, NA,
                    "which processes to fit: ignition, escape, or both (the default)"),
    defineParameter("modelAlgorithm", "character", "xgboost", NA, NA,
                    "Model type. Only `xgboost` (any value containing 'xgb') works."),
    defineParameter(".plots", "character", default = "screen",
                    desc = paste("See `?Plots`. If set, plots the predicted response against each covariate,",
                                 "for climate and fuel covariates separately, and saves it as png in `figurePath(sim)`.")),
    defineParameter(".runInitialTime", "numeric", default = start(sim),
                    desc = "when to start this module? By default, the start time of the simulation."),
    defineParameter(".runInterval", "numeric", default = NA,
                    desc = paste("optional. Interval between two runs of this module,",
                                 "expressed in units of simulation time.",
                                 "By default, NA, which means that this module only runs once per simulation.")),
    defineParameter(".seed", "list", NULL, NA, NA,
                    paste("Named list of seeds to use for each event (names).",
                          "E.g., `list('init' = 123)` will `set.seed(123)`",
                          "at the start of the init event and unset it at the end.",
                          "Defaults to `NULL`, meaning that no seeds will be set.")),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    desc = paste("Should this entire module be run with caching activated?",
                                 "This is generally intended for data-type modules,",
                                 "where stochasticity and time are not relevant.")),
    defineParameter(".neverCache", "character", "checkData", NA, NA,
                    paste("Events that are never cached, whatever `.useCache` says (read by SpaDES.core, 3.2.1.9037 or later):",
                          "`checkData` reads the ignition-fit ledger on Google Drive, which is not in the event's cache key, so a cached `checkData` would replay a stale answer to whether this polygon already has a fit.")),
    defineParameter("ignitionFitGoogleDriveFolder", "character",
                    fireSenseUtils::ignitionFitLedgerFolder(),
                    NA, NA, paste("Google Drive folder url holding the shared ignition/escape fit ledger",
                                  "(`ignitionFitFilename`). `NULL` keeps the ledger purely local, in",
                                  "`inputPath(sim)`, with no Google Drive access at all (e.g. for tests).")),
    defineParameter("escapeSizeHa", "numeric", fireSenseUtils::fireSenseEscapeSizeHa, NA, NA,
                    paste("Size (ha) of the whole fire polygon at and above which a fire escaped: the escape definition",
                          "(`fireSenseUtils::escapeDefinition()`) the escape model is fitted under, as in `fireSense_dataPrepFit`",
                          "(every module that has it must agree: an error otherwise). Recorded in the ledger row, so",
                          "`fireSenseUtils::latestIgnitionFits()` uses only fits made under the run's definition.")),
    defineParameter("fireCauses", "character", fireSenseUtils::fireSenseFireCauses, NA, NA,
                    paste("The fire cause groups of the escape definition, as in `fireSense_dataPrepFit` (which must agree:",
                          "an error otherwise). Recorded in the ledger row with `escapeSizeHa`.")),
    defineParameter("ignitionFitFilename", "character", "latest",
                    NA, NA, paste("File name of the shared fit ledger: an `sf` object with one row per polygon, holding",
                                  "its fitted `fireSense_IgnitionFitted` and `fireSense_EscapeFitted`. A name (e.g. an",
                                  "experiment's, starting `fireSenseIgnitionParams_`) is read and written as it is. `\"latest\"`",
                                  "(the default) reads and writes the file named for this fit's fire years,",
                                  "`fireSenseUtils::ignitionFitFilenameFor()`, e.g.",
                                  "`fireSenseIgnitionParams_1986-2024_xgboostPoisson_20261013.rds`; with no fit data, this ELF's",
                                  "newest usable fit in any ledger file (`fireSenseUtils::latestIgnitionFits()`).")),
    defineParameter("xgbMaxDepth", "numeric", 2, NA, NA,
                    paste("xgboost `max_depth`, the depth of each tree, of the ignition and the escape fits. Few, shallow trees",
                          "(the default 2) and shrinkage pull the predictions toward the mean; deeper trees concentrate them.",
                          "The defaults of the five `xgb*` parameters are the values the fits were hard-coded with,",
                          "so a fit with them is the one made before they existed (`xgbSettings()`).")),
    defineParameter("xgbLearningRate", "numeric", 0.15, NA, NA,
                    "xgboost `learning_rate`, the shrinkage of each boosting round."),
    defineParameter("xgbNrounds", "numeric", 100, NA, NA,
                    paste("xgboost `nrounds`, the number of boosting rounds; with `xgbEarlyStopping`, the most a fold",
                          "may use.")),
    defineParameter("xgbRegLambda", "numeric", 0.5, NA, NA, "xgboost `reg_lambda`, the L2 penalty on the leaf weights."),
    defineParameter("xgbEarlyStopping", "numeric", NA_real_, NA, NA,
                    paste("`NA` (default): off, every fold boosts `xgbNrounds` rounds. A number of rounds: each fold's",
                          "boosting stops once the Poisson loss of its held-out rows (`cvFolds()`, or the whole years of",
                          "`yearBlockedFolds()` for the year-blocked models) has not improved for this many rounds,",
                          "and the fold's model keeps the best round (at most `xgbNrounds`). The fold is chosen",
                          "with its own held-out rows, so its predictions on them are slightly optimistic.")),
    defineParameter("ignitionYearClimate", "character", "CMDsm", NA, NA,
                    paste("The climate covariate whose yearly ELF mean scales each year's expected ignitions,",
                          "`exp(beta * z)` with `beta >= 0` fitted to the out-of-year totals and `z` the standardised yearly",
                          "mean (`fitYearClimate()`); fixed for every ELF; revisit when all ELFs are fitted. With or without",
                          "underscore (`CMD_sm` is `CMDsm`); stops if it is not a column of `fireSense_ignitionCovariates`.",
                          "`NA`: no term.")),
    defineParameter("doFit", "logical", FALSE, NA, NA,
                    paste("Whether to fit this polygon. FALSE (default): never fit -- use the ledger's fits, and",
                          "stop with an error if it holds none for this polygon (`ignitionFitFilename` picks the ledger), so",
                          "a run that only predicts can never start a fit by accident. NA: do whatever is needed --",
                          "use the ledger's ignition and escape fits if it holds them for this polygon (and they fit",
                          "the current covariates), otherwise fit. TRUE: always fit, replacing the ledger row (FOR",
                          "DEVELOPERS ONLY: set it when the fit's INPUTS have changed, so the stored row is stale).")),
    defineParameter("refitExisting", "logical", NA, NA, NA,
                    paste("DEPRECATED, use `doFit`. TRUE is `doFit = TRUE`; NA (unset) and FALSE do nothing.",
                          "Overrides `doFit = FALSE` (the default), not `doFit = NA` or TRUE."))
  ),
  inputObjects = bindrows(
    expectsInput(".ELFind", "character",
                 desc = paste("Identifier of the polygon being fit, e.g. '6.1.1'. This becomes the",
                              "`polygonID` of the row this module writes to the shared cloud fit ledger",
                              "(`ignitionFitFilename` in `ignitionFitGoogleDriveFolder`). Only used, and only",
                              "required, when `studyArea` is supplied."),
                 sourceURL = NA),
    expectsInput("studyArea", "SpatVector",
                 desc = paste("Polygon being fit; its geometry and crs go in the ledger row. Optional: without it,",
                              "this module fits every time and never reads or writes the shared ledger, as before",
                              "this parameter existed."),
                 sourceURL = NA),
    expectsInput("fireSense_ignitionCovariates", "data.frame",
                 desc = paste("Table of aggregated ignition covariates with annual `ignitions` counts,",
                              "one row per `pixelID` and `year`. An optional column `coreFraction` (from fireSense_dataPrepFit)",
                              "is the share of the pixel in the ELF's core, where the ignitions are counted: the ignition fit",
                              "weights the pixel by it and is not given it as a covariate. Without it every pixel counts fully.")),
    expectsInput("ignitionFitRTM", "SpatRaster",
                 desc = paste("Template raster at the resolution and extent of `fireSense_ignitionCovariates`.",
                              "Its resolution and its number of non-NA cells, which must be in the attribute `nonNAs`",
                              "(`attributes(ignitionFitRTM)$nonNAs`), are stored in the fitted objects."))
  ),
  outputObjects = bindrows(
    createsOutput("fireSense_EscapeFitted", "fireSense_EscapeFit",
                  desc = paste("List of `modelList` and `scaleData`, as `fireSense_IgnitionFitted`,",
                               "with `modelList` of class `fireSense_EscapeFit` (its `yearEffects` multiply the expected",
                               "escapes of an ignited pixel). Either freshly fitted, or, when",
                               "`studyArea` matches a row of the shared ledger, read from it.")),
    createsOutput("fireSense_IgnitionFitted", "fireSense_IgnitionFit",
                  desc = paste("List of `modelList` (class `fireSense_IgnitionFit`: `model`, the per-fold xgboost models",
                               "and their ROC curves and out-of-fold `yearTotals`, `fittingRes`, `lambdaRescaleFactor`,",
                               "`rescales`, `yearSD`: the sd of the year effect that `fireSense_burn` draws; absent in fits",
                               "stored before it existed, which means 0; and, when both processes are fitted, `yearEffects`:",
                               "the fitted multiplier of each fit year (`fittedYearEffects()`), which a hindcast of that year",
                               "may use instead of a draw, `fireSense_burn`'s `yearEffects = \"fitted\"`; absent in older fits)",
                               "and `scaleData` (centre and scale used to standardise the covariates). Either freshly",
                               "fitted, or, when `studyArea` matches a row of the shared ledger, read from it.")),
    createsOutput("ignitionFitPreRun", "data.frame",
                  desc = paste("Only when `studyArea` is supplied: ledger rows that overlap it, from `CacheGeo`",
                               "(`NULL` if there is no previous fit): a `geometry` column, `polygonID`, and",
                               "`fireSense_IgnitionFitted`/`fireSense_EscapeFitted` list-columns."))
  )
))

## event types
#   - type `init` is required for initialiazation

#' Event dispatcher
#'
#' `init` schedules `checkData` and `run` at `.runInitialTime`; `run` repeats every
#' `.runInterval` if that is not `NA`.
#'
#' @param sim A `simList`.
#' @param eventTime Time of the event.
#' @param eventType One of `init`, `checkData`, `run`.
#' @param debug Unused.
#' @return `sim`, invisibly.
doEvent.fireSense_ignitionFit = function(sim, eventTime, eventType, debug = FALSE) {
  moduleName <- current(sim)$moduleName

  switch(
    eventType,
    init = {
      #the reason this is not scheduling Init is because it is only performing sanity checks on data
      # which is created during other modules non-init events. And since Inits are all scheduled first...
      sim <- scheduleEvent(sim, P(sim)$.runInitialTime, moduleName, "checkData")

      sim <- scheduleEvent(sim, P(sim)$.runInitialTime, moduleName, "run")
    },
    checkData = {

      sim <- Init(sim)
    },
    run = {
      sim <- frequencyFitRun(sim)

      if (!is.na(P(sim)$.runInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, moduleName, "run")
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )

  return(invisible(sim))
}

### template initialization
#' Check inputs before fitting (the `checkData` event)
#'
#' When `sim$studyArea`
#' is supplied, also reads the shared ignition/escape fit ledger with `CacheGeo` into
#' `sim$ignitionFitPreRun`, so `run` can tell whether this polygon already has a fit.
#'
#' @param sim A `simList`.
#' @return `sim`, invisibly, unchanged except for `sim$ignitionFitPreRun`.
Init <- function(sim) {

  if (!any(c("ignition", "escape") %in% P(sim)$whichProcessesToFit)) {
    stop("please review P(sim)$whichProcesesToFit...ensure lower-case")
  }

  if (!is.null(sim$studyArea)) {
    sa <- sim$studyArea
    if (inherits(sa, "SpatVector")) sa <- sf::st_as_sf(sa)
    sim$ignitionFitPreRun <-
      if (identical(P(sim)$ignitionFitFilename, "latest") && !length(ignitionFitYears(sim))) {
        ## no fit data to name the file from: the newest row of this ELF over the ledger files, as fireSense_ELFs reads it
        fireSenseUtils::latestIgnitionFits(P(sim)$ignitionFitGoogleDriveFolder, destinationPath = inputPath(sim),
                                           polygonIDs = as.character(sim$.ELFind),
                                           escapeDefinition = fireSenseUtils::escapeDefinitionOf(sim))
      } else {
        ignitionLedgerCacheGeo(sim, domain = sa, action = "nothing", useCache = FALSE, bufferOK = TRUE)
      }
  }

  return(invisible(sim))
}

#' Fit the requested processes (the `run` event)
#'
#' When `sim$studyArea` is supplied and the ledger already has a row for `sim$.ELFind` (and
#' `doFit` is not `TRUE`) and the stored models' features are all among the current
#' covariate columns (`useExistingFitIsCompatible()`), the stored
#' `fireSense_IgnitionFitted`/`fireSense_EscapeFitted` are used and the fit is skipped. With `doFit` FALSE
#' and no usable ledger fit the run stops. Otherwise calls `buildModelsFitModels()` for each of
#' `P(sim)$whichProcessesToFit`, assigns the results to `sim$fireSense_IgnitionFitted` and/or
#' `sim$fireSense_EscapeFitted`, and, when `sim$studyArea` is supplied, writes them to the ledger.
#'
#' @param sim A `simList`. A fit of ignition needs the `nonNAs` attribute on `sim$ignitionFitRTM`; a run that
#'   uses the ledger's fit does not need `sim$ignitionFitRTM`.
#' @return `sim`, invisibly.
frequencyFitRun <- function(sim) {

  doFit <- fireSenseUtils::resolveDoFit(P(sim)$doFit, P(sim)$refitExisting,
                                        module = "fireSense_ignitionFit", message = TRUE)
  if (!fireSenseUtils::mustFit(doFit, hasPreRunIgnitionFitForThisPolygon(sim) && useExistingFitIsCompatible(sim),
                               paste(sim$.ELFind, collapse = ", "))) {
    return(invisible(useExistingIgnitionFit(sim)))
  }
  if (length(sim$.ELFind) > 1L)
    stop("fireSense_ignitionFit: a fit is of one ELF, but this run covers several (", paste(sim$.ELFind, collapse = ", "),
         "). Fit each on its own (`.ELFind`), or predict from their stored fits (`doFit = FALSE`, the default).",
         call. = FALSE)

  if ("ignition" %in% P(sim)$whichProcessesToFit &&
      (is.null(attributes(sim$ignitionFitRTM)$nonNAs) || length(attributes(sim$ignitionFitRTM)$nonNAs) == 0))
    stop("sim$ignitionFitRTM@data@attributes$nonNAs must be a non-empty/non-NULL numeric")

  #TODO: make cache smart by digest args in advance
  whichProcessesToFit <- P(sim)$whichProcessesToFit
  names(whichProcessesToFit) <- igOrEscNames(whichProcessesToFit, post = "Fitted", case = "sen")

  fireSense_FittedModels <-
    purrr::pmap(.l  = list(igOrEsc = whichProcessesToFit), sim = sim,
                .f = buildModelsFitModels)

  # put to sim sim$fireSense_IgnitionFitted, sim$fireSense_EscapeFitted
  list2env(fireSense_FittedModels, envir = envir(sim))
  ## the escape effect is that of the escapes the fitted ignitions give, so it needs both fits
  if (all(c("ignition", "escape") %in% whichProcessesToFit))
    sim <- addFittedYearEffects(sim)

  if (!is.null(sim$studyArea))
    sim <- writeIgnitionLedgerRow(sim)

  return(invisible(sim))
}

#' Fit the model for one process
#'
#' Calls `runXGBOOST()`. Only xgboost algorithms are supported; the non-xgboost modelling path
#' was removed.
#'
#' @param covariates data.table of response and covariates. The response columns are
#'   converted to integer by reference.
#' @param type `"ignition"` or `"escape"`.
#' @param dig Digest of the unscaled covariates, passed to `runXGBOOST()` as its cache key.
#' @param modelAlgorithm See the `modelAlgorithm` module parameter; must contain `"xgb"`.
#' @param nFolds Number of cross-validation folds.
#' @param xgb Boosting settings, `xgbSettings()`.
#' @return See `runXGBOOST()`.
buildModel <- function(covariates, type = "ignition",
                       dig,
                       modelAlgorithm, nFolds = 5, xgb = xgbSettings()) {

  if (!any(grepl("xgb", modelAlgorithm)))
    stop("modelAlgorithm must be an xgboost algorithm; the non-xgboost path was removed.")

  dat <- covariates

  varToInt <- fireSenseUtils::ignitionsTxt
  if (type == "escape") {
    varToInt <- c(fireSenseUtils::ignitionsTxt, fireSenseUtils::escapesTxt)
  }

  for (i in varToInt)
    set(dat, NULL, i, as.integer(dat[[i]]))

  runXGBOOST(dat, dig, type = type, nFolds = nFolds, xgb = xgb)
}


#' Default inputs
#'
#' Stops if `fireSense_ignitionCovariates` is not supplied.
#'
#' @param sim A `simList`.
#' @return `sim`.
.inputObjects <- function(sim) {
  if (!suppliedElsewhere("fireSense_ignitionCovariates", sim)) {
    stop("this module does not produce data - consider usin the module 'PredictiveEcology/fireSense_dataPrepFit'")
  }

  return(sim)
}

#' Digest of the covariates used for the fit
#'
#' `fireSenseUtils::prepareCovariatesOuter()` returns its digest as
#' `digestOfData = list(unscaledData = <hash>)`. The module used to ask for
#' `digestOfData$fireSense_ignitionCovariates`, which is not an element of that list, so `NULL`
#' was silently passed down as the cache key of `runXGBOOST()`, and the data did not take part
#' in it. This takes the digest that is actually there, and stops if it is missing.
#'
#' @param data Output of `fireSenseUtils::prepareCovariatesOuter()`.
#' @return Character digest of the unscaled covariates.
dataDigest <- function(data) {
  dig <- data$digestOfData
  if (is.list(dig))
    dig <- dig[["unscaledData"]]
  if (is.null(dig) || !length(dig) || !nzchar(paste0(dig, collapse = "")))
    stop("prepareCovariatesOuter() returned no digest of the covariates; ",
         "the fit would be cached without the data in its key.", call. = FALSE)
  dig
}

#' The covariates of one process, prepared as its fit has them
#'
#' `prepareCovariatesOuter()` of `sim$fireSense_ignitionCovariates` or `sim$fireSense_escapeCovariates`, with the
#' module's `modelAlgorithm` and `rescaleVars`. The rows stay in the order of the unprepared table.
#'
#' @inheritParams buildModelsFitModels
#' @param modelAlgorithm,rescaleVars The module's parameters.
#' @return The result of `prepareCovariatesOuter()`: `covariates` (standardised when `rescaleVars`) and more.
fitCovariates <- function(igOrEsc, sim, modelAlgorithm = Par$modelAlgorithm, rescaleVars = Par$rescaleVars) {
  prepareCovariatesOuter(sim[[igOrEscNames(igOrEsc, post = "Covariates")]],
                         algorithm = modelAlgorithm, rescaleVars = rescaleVars)
}

#' Add the fitted year effects and the chained-escape factor to the ignition and escape fits
#'
#' Predicts every pixel and year of the fit data with each fit (`fireSenseUtils::foldMeanPrediction()`; the ignition
#' prediction is multiplied by the year climate term, `modelList$yearClimate`, so the fitted multipliers are relative
#' to it) and stores `fittedYearEffects()` as `modelList$yearEffects` of `sim$fireSense_IgnitionFitted` and
#' `sim$fireSense_EscapeFitted`, for `fireSense_ignitionPredict` to pass to `fireSense_burn`, and
#' `chainedEscapeFactor()` of the same predictions as `modelList$concentration$chain` of the escape fit (logged).
#'
#' @param sim A `simList` with both fits.
#' @return `sim`, invisibly.
addFittedYearEffects <- function(sim) {
  one <- function(igOrEsc, fit, response) {
    raw <- sim[[igOrEscNames(igOrEsc, post = "Covariates")]]
    data.table(pixelID = raw$pixelID, year = raw[[fireSenseUtils::yearTxt]], observed = raw[[response]],
               mean = fireSenseUtils::foldMeanPrediction(fit, fitCovariates(igOrEsc, sim)$covariates) *
                 (if (igOrEsc == "ignition") yearClimateRowMultiplier(fit$modelList$yearClimate, raw) else 1),
               weight = exposureOf(raw, "ignition"))
  }
  ignitions <- one("ignition", sim$fireSense_IgnitionFitted, fireSenseUtils::ignitionsTxt)
  escapes <- one("escape", sim$fireSense_EscapeFitted, fireSenseUtils::escapesTxt)
  eff <- fittedYearEffects(ignitions, escapes)
  sim$fireSense_IgnitionFitted$modelList$yearEffects <- eff$ignition
  sim$fireSense_EscapeFitted$modelList$yearEffects <- eff$escape
  sim$fireSense_EscapeFitted$modelList$concentration$chain <- chainedEscapeFactor(ignitions, escapes)
  message(chainMessage(sim$fireSense_EscapeFitted$modelList$concentration$chain))
  invisible(sim)
}

#' Prepare covariates, fit one process, and plot it
#'
#' @param igOrEsc `"ignition"` or `"escape"`. Covariates are read from
#'   `sim$fireSense_ignitionCovariates` or `sim$fireSense_escapeCovariates`.
#' @param sim A `simList`.
#' @return List of `modelList`, of class `fireSense_IgnitionFit` or `fireSense_EscapeFit`, and
#'   `scaleData`, the centre and scale used to standardise the covariates. `modelList` holds
#'   `model` (from `buildModel()`), `rescales`, `fittingRes` (resolution of `ignitionFitRTM`),
#'   `lambdaRescaleFactor` (rows of covariates / non-NA cells of `ignitionFitRTM`) and `yearSD`, the sd
#'   of the year effect on the predicted mean (`fireSenseUtils::yearEffectSD()` of the out-of-fold
#'   `model$yearTotals`, for ignitions after the year climate term), `yearClimate` (ignitions only, absent
#'   when `ignitionYearClimate` is `NA`: `fitYearClimate()`), which `fireSense_burn` draws each year, with one z shared by ignitions, escapes and spread,
#'   and `concentration` (`concentrationOf()`, logged by `concentrationMessage()`).
buildModelsFitModels <- function(igOrEsc, sim) {
  data <- fitCovariates(igOrEsc, sim)
  ## `rescaleCovariates()` does not scale coreFraction, and it is in the unscaled data `digestOfData` hashes
  nFolds <- 5

  modelHere <- buildModel(covariates = data$covariates,
                          type = igOrEsc, nFolds = nFolds,
                          dig = dataDigest(data),
                          modelAlgorithm = Par$modelAlgorithm, xgb = xgbSettingsOfParams(Par)) |>
    Cache(omitArgs = c("covariates"),
          .functionName = paste0("buildModel for ", igOrEsc),
          .cacheExtra = c(data$digestOfData, list(callees = buildModelCallees())),
          cacheSaveFormat = "rds") # doesn't work with qs
  #ignition specific
  origNoPix <- attributes(sim$ignitionFitRTM)$nonNAs   ## nrow(preSampleData) in eg above
  finalNoPix <- nrow(data$covariates)     ## nrow(postSampleData) in eg above
  lambdaRescaleFactor <- finalNoPix/origNoPix

  concentration <- modelHere$concentration # beside `yearEffects`, not among the boosters; absent in older fits
  modelHere$concentration <- NULL
  message(concentrationMessage(concentration, igOrEsc))

  ## the year climate term of the ignitions; `yearSD` is the spread left after it
  yearClimate <- NULL
  yearTotals <- modelHere$yearTotals
  if (identical(igOrEsc, "ignition") && length(Par$ignitionYearClimate) == 1L && !is.na(Par$ignitionYearClimate)) {
    raw <- sim$fireSense_ignitionCovariates
    if (!is.null(yearTotals))
      yearTotals$year <- calendarYearOf(yearTotals$year, data$covariates[[fireSenseUtils::yearTxt]], raw[[fireSenseUtils::yearTxt]])
    yearClimate <- fitYearClimate(Par$ignitionYearClimate, raw, yearTotals)
    if (!is.null(yearTotals))
      yearTotals$predicted <- yearTotals$predicted * yearClimateRowMultiplier(yearClimate, raw, yearTotals$year)
  }

  modelList <- list(
    model = modelHere,
    rescales = data$ignitionRescalers,
    fittingRes = res(sim$ignitionFitRTM)[1],
    lambdaRescaleFactor = lambdaRescaleFactor,
    yearSD = fireSenseUtils::yearEffectSD(yearTotals$observed, yearTotals$predicted),
    concentration = concentration)
  modelList$yearClimate <- yearClimate   # absent (not NULL) without a term, as in fits made before it

  class(modelList) <- igOrEscNames(igOrEsc, post = "Fit", case = "Title")

  if (anyPlotting(P(sim)$.plots)) {
    message("Plotting ", igOrEsc, "...")

    figPath <- figurePath(sim)
    digModels <- attr(modelHere, "tags")
    modelOnly <- fireSenseUtils::foldModels(list(modelList = list(model = modelHere)))
    aa <- setupPlots(modelOnly, dat = data$covariates, igOrEsc) |>
      Cache(omitArgs = c("dat", "modelOnly"),
            .cacheExtra = list(digestOfData = data$digestOfData, digModels = digModels, plotPredictions = plotPredictions),
            .functionName = paste0(".functionName_", igOrEsc))
    fn1 <- functionNameHelper("FuelClimate", ifelse(igOrEsc == "escape", "", "Lightning"), "predicted", igOrEsc, Par$crossValType[1])
    fn <- functionNameHelper(fn1, format(Sys.time()))
    bb <- quote(ggarrange(plotlist = aa))

    a <- Plots(bb, path = figPath,
               filename = fn,
               type = "png", ggsaveArgs = list(width = 8, height = 5, scale = 1.5), useCache = TRUE)
  }
  list(modelList = modelList, scaleData = attr(data$covariates, "scaleData"))
}


#' Fit cross-validated xgboost models
#'
#' One model per fold, Poisson objective (`fireSenseUtils::ignitionObjectiveTxt`, `count:poisson`), trained on all rows except the fold's, which are the
#' `eval_set`. Folds come from `cvFolds()`. For `type = "escape"` only rows with
#' `ignitions > 0` are used. The RNG seed is fixed, and restored on exit, so folds repeat and
#' caching works.
#'
#' @param dat data.table of response and covariates. `pixelID`, `year`, `yearChar` and `coreFraction`
#'   columns are not features. The ignition fit has the weight `coreFraction` (see `exposureOf()`): rows
#'   without core are dropped, and the model is fitted to `ignitions / coreFraction`. Everything else keeps
#'   the integer response: the folds, `stopIfFoldsLackPositives()` and the ROC curves. The predictions
#'   in `valData` and `yearTotals` are expected counts, the model's prediction times `coreFraction`.
#' @param dig Digest of the unscaled covariates, used to cache the creation of dummy variables.
#' @param type `"ignition"` or `"escape"`; also the pattern that finds the response column.
#' @param nFolds Number of folds.
#' @param xgb Boosting settings, `xgbSettings()` (the module's `xgb*` parameters).
#' @return List of the `nFolds` xgboost models, named `Fold1`, ..., plus `rocs` from `rocPerFold()`,
#'   `yearTotals` from `outOfFoldYearTotals()`: out of year, from extra models fitted on `yearBlockedFolds()`
#'   (not returned), or out of fold from the `cvFolds()` models with fewer than 2 years or when a year-blocked
#'   fold would train without a positive, and `concentration`: `concentrationOf()` of the same out-of-year (or
#'   out-of-fold) predictions and the response, with the rounds each fold model kept. With `xgbEarlyStopping` the
#'   held-out rows of a fold also chose its rounds (`earlyStoppedRounds()`), so these out-of-year (out-of-fold) lift,
#'   Gini index and deviance are slightly optimistic.
runXGBOOST <- function(dat, dig, type = "ignition", nFolds = 5, xgb = xgbSettings()) {
  dat3Forxgboost <- dat

  # Add dummy variables for factor columns -- i.e., the random effects
  if (all(sapply(dat3Forxgboost, is.numeric)) %in% FALSE)
    dat3Forxgboost <-
      model.matrix(~ . + 0, data = dat3Forxgboost) |>
      Cache(omitArgs = c("object", "data", "x"), .cacheExtra = dig)# Creates dummy variables


  savedSeed <- .Random.seed
  on.exit(assign(".Random.seed", savedSeed, envir = .GlobalEnv), add = TRUE)
  set.seed(12345) # so kfolds are same, so Caching works correctly below; if dat3 changes number of rows,
  # it will be a totally different sequence; but it will be the same sequence
  # if number of rows doesn't change

  # Setup k-folds
  indexNames <- c("keepAll", "keepEval")

  if (identical(type, "escape")) {
    dat3Forxgboost <- dat3Forxgboost[ignitions > 0]
  }
  if (identical(type, "ignition") && !is.null(dat3Forxgboost[[fireSenseUtils::coreFractionTxt]]))
    dat3Forxgboost <- dat3Forxgboost[dat3Forxgboost[[fireSenseUtils::coreFractionTxt]] > 0, ]
  weight <- exposureOf(dat3Forxgboost, type)
  trainIndexK <- cvFolds(dat3Forxgboost[[grep(value = TRUE, type, colnames(dat3Forxgboost))]], nFolds)
  yearOfRow <- dat3Forxgboost[[fireSenseUtils::yearTxt]] # dropped from the features below

  colOrder <- setdiff(colnames(dat3Forxgboost), c("pixelID", "year", fireSenseUtils::coreFractionTxt))
  colOrder <- colOrder[colOrder %in% grep("yearChar", colnames(dat3Forxgboost), invert = TRUE, value = TRUE)]
  colOrder <- sample(colOrder)
  dat3Forxgboost <- dat3Forxgboost[, ..colOrder]
  dig <- c(.robustDigest(dat3Forxgboost), if (!all(weight == 1)) .robustDigest(weight))
  colnamesNoIgn <- fireSenseUtils::xgbFeatureNames(colnames(dat3Forxgboost)) # alphabetical
  dat3ForxgboostNoIgn <- dat3Forxgboost[, ..colnamesNoIgn]

  ignOrEscapeColName <- grep(value = TRUE,type, colnames(dat3Forxgboost))
  response <- dat3Forxgboost[[ignOrEscapeColName]]
  stopIfFoldsLackPositives(response, trainIndexK, type)
  st <- system.time(
    mm <- purrr::pmap(
      list(valInd = trainIndexK, kFold = seq(nFolds)),
      function(valInd, kFold)
      {

        wholeDataset <- valInd[[indexNames[[1]]]]
        valInd <- valInd[[indexNames[[2]]]]
        digValInd <- .robustDigest(valInd) # should be eval set

        mPoisson <- fitFoldXgboost(dat3ForxgboostNoIgn[wholeDataset], # should be whole set
                                   response[wholeDataset] / weight[wholeDataset],
                                   valInd, # should be keepEval
                                   functionNameHelper("xgboost", type, kFold), c(dig, digValInd, type),
                                   weight = weight[wholeDataset], xgb = xgb)
        # Predict probabilities
        valData <- dat3Forxgboost[valInd, ]
        pred2 <- predict(mPoisson, valData)
        if (!all(is.finite(pred2))) {
          stop(type, " cross-validation fold ", kFold, ": the model predicted ", sum(!is.finite(pred2)),
               " non-finite values out of ", length(pred2), call. = FALSE)
        }
        valData <- cbind(valData, predPoisson = pred2 * weight[valInd]) # expected count

        shap_values <- shap.values(mPoisson, dat3ForxgboostNoIgn) |>
          Cache(omitArgs = formalArgs(shap.values),
                .functionName = functionNameHelper("shap.values", type, kFold),
                .cacheExtra = c(dig, digValInd, type))
        shapContrib <- shap_values$shap_score
        shapContrib <- shapContrib[, -"(Intercept)"]
        shap_long <- shap.prep(
          shap_contrib = shapContrib, X_train = dat3ForxgboostNoIgn)  |>
          Cache(omitArgs = formalArgs(shap.prep),
                .functionName = functionNameHelper("shap.prep", type, kFold),
                .cacheExtra = c(dig, digValInd, type))
        list(valData = valData, mod = mPoisson) # not shap_long: nFolds long tables exhausted memory
      })
  )
  rocs <- rocPerFold(mm, ignOrEscapeColName)
  aucs <- aucPerFold(rocs)
  print(meanRocMessage(aucs))
  folds <- trainIndexK
  basis <- "out-of-fold"
  heldOut <- lapply(mm, function(m) m$valData$predPoisson)
  yearTotals <- outOfFoldYearTotals(yearOfRow, response, folds, heldOut)
  ## the yearSD comes from out-of-year totals: refit with whole years held out, when every fold can train
  yearFolds <- yearBlockedFolds(yearOfRow, nFolds)
  if (!is.null(yearFolds) && !any(foldsLackPositives(response, yearFolds))) {
    folds <- yearFolds
    basis <- "out-of-year"
    heldOut <- Map(fold = yearFolds, kFold = seq_along(yearFolds), f = function(fold, kFold) {
      m <- fitFoldXgboost(dat3ForxgboostNoIgn, response / weight, fold$keepEval,
                          functionNameHelper("xgboostYearBlocked", type, kFold),
                          c(dig, .robustDigest(fold$keepEval), type), weight = weight, xgb = xgb)
      predict(m, dat3Forxgboost[fold$keepEval, ]) * weight[fold$keepEval]
    })
    yearTotals <- outOfFoldYearTotals(yearOfRow, response, folds, heldOut)
  }
  ## the same out-of-year (else out-of-fold) predictions, row by row, for the concentration diagnostic
  mm2 <- Map(m = mm, function(m) m$mod)
  concentration <- concentrationOf(response, outOfFoldPrediction(length(response), folds, heldOut), basis = basis)
  concentration$rounds <- vapply(mm2, function(m) nrow(attr(m, "evaluation_log")), numeric(1), USE.NAMES = FALSE)
  mm2 <- append(mm2, list(rocs = rocs, yearTotals = yearTotals, concentration = concentration))

  return(mm2)
}


#' ROC curve per validation fold
#'
#' AUC needs both outcomes present. A validation fold can legitimately hold only one --
#' for escape in a small study area, every fire escaped or none did -- and `pROC::roc()`
#' stops there with "'response' must have two levels". That used to abort the whole fit,
#' which is the wrong trade: these curves are a diagnostic, printed and returned alongside
#' the models, and nothing about the fitted model (`m$mod`) depends on them. A fold with
#' one outcome therefore reports no AUC and the fit stands.
#'
#' @param mm List with one element per fold, each holding `valData`: the validation rows with
#'   the response and the prediction, `predPoisson`.
#' @param ignOrEscapeColName Name of the response column.
#' @return List of `pROC::roc` objects; `NULL` for a fold with only one outcome.
rocPerFold <- function(mm, ignOrEscapeColName) {
  lapply(mm, function(d) {
    ignZeroAndOnes <- pmin(1L, d$valData[[ignOrEscapeColName]])
    if (length(unique(ignZeroAndOnes)) < 2L) {
      message("  ... a validation fold of ", ignOrEscapeColName, " holds only the value ",
              unique(ignZeroAndOnes), "; AUC is undefined for it, so it is skipped")
      return(NULL)
    }
    pROC::roc(ignZeroAndOnes, d$valData[["predPoisson"]])
  })
}

#' AUC of each fold
#'
#' @param rocs Output of `rocPerFold()`.
#' @return Numeric vector; `NA` where the fold has no ROC curve.
aucPerFold <- function(rocs) {
  vapply(rocs, function(r) if (is.null(r)) NA_real_ else as.numeric(r$auc), numeric(1))
}

#' Cross-validation folds for runXGBOOST(), stratified on whether each row is positive.
#'
#' Folds built on the raw response put rare positives in one group, so a single fold could hold every
#' positive and train on zeros only (ELF 3.1.2). Stratified on positive / not, the positives are spread
#' across the folds, and every fold trains with some once there are 2 or more.
#'
#' @param response Numeric response.
#' @param nFolds Number of folds.
#' @return List of `nFolds` lists, each with `keepAll` (all row indices) and `keepEval` (the
#'   fold's held-out row indices).
cvFolds <- function(response, nFolds) {
  folds <- caret::createFolds(factor(response > 0), k = nFolds, list = TRUE, returnTrain = FALSE)
  lapply(folds, function(f) list(keepAll = seq_along(response), keepEval = f))
}

#' Stop when a cross-validation fold would train without any positive observation.
#'
#' xgboost trains each fold's model on `keepAll` without the fold's `eval_set` (`keepEval`). If those
#' rows hold no positive, the model predicts NaN for every row (ELF 3.1.2: one ignition in 2002-2022,
#' which landed in one fold's held-out rows), and pROC::roc() then stops with "No control observation".
#'
#' @param response Numeric response.
#' @param folds Output of `cvFolds()`.
#' @param type `"ignition"` or `"escape"`; used in the error message.
#' @return `NULL`, invisibly. Stops if any fold would train without a positive.
stopIfFoldsLackPositives <- function(response, folds, type) {
  noPositive <- foldsLackPositives(response, folds)
  if (any(noPositive)) {
    stop("Too few ", type, "s to fit ", length(folds), " cross-validation folds: ",
         sum(response > 0), " positive observation(s), and ", sum(noPositive),
         " fold(s) would train without any.", call. = FALSE)
  }
}

#' Which folds would train without any positive observation
#'
#' @param response Numeric response.
#' @param folds Output of `cvFolds()` or `yearBlockedFolds()`.
#' @return Logical, one per fold.
foldsLackPositives <- function(response, folds) {
  vapply(folds, function(f) !any(response[setdiff(f$keepAll, f$keepEval)] > 0), logical(1))
}

#' Cross-validation folds of whole years
#'
#' Each fold holds out every row of its years, so a fold model's predictions for those years are out of
#' year. Years go to folds at random (the caller sets the seed); with fewer years than `nFolds`, one fold
#' per year.
#'
#' @param year The year of each row; `NULL` if the data had no year column.
#' @param nFolds Number of folds.
#' @return Like `cvFolds()`; `NULL` without `year` or with fewer than 2 years.
yearBlockedFolds <- function(year, nFolds) {
  years <- sort(unique(year))
  if (length(years) < 2) return(NULL)
  nFolds <- min(nFolds, length(years))
  foldOfYear <- sample(rep_len(seq_len(nFolds), length(years)))
  lapply(seq_len(nFolds), function(k)
    list(keepAll = seq_along(year), keepEval = which(year %in% years[foldOfYear == k])))
}

#' Fit one fold's xgboost model
#'
#' Poisson objective (`fireSenseUtils::ignitionObjectiveTxt`), trained on the rows of `x` not in
#' `evalSet`, which xgboost holds out for evaluation. Cached; xgboost objects do not save with `qs`,
#' so `rds`. With early stopping (`xgb$early_stopping_rounds`) the number of rounds is first chosen on the
#' held-out rows by `earlyStoppedRounds()` (also cached), and the model is fitted with that many.
#'
#' @param x Covariates (no response).
#' @param y Response (label); `ignitions / weight` when the rows are weighted.
#' @param evalSet Row indices held out of training.
#' @param functionName `.functionName` for `Cache()`.
#' @param cacheExtra `.cacheExtra` for `Cache()`: digests of the data and of `evalSet`, and the type.
#' @param weight Weight of each row of `x`, also for the held-out rows (`NULL`: all 1).
#' @param xgb Boosting settings, `xgbSettings()`.
#' @return An xgboost model.
fitFoldXgboost <- function(x, y, evalSet, functionName, cacheExtra, weight = NULL, xgb = xgbSettings()) {
  nrounds <- xgb$nrounds
  if (!is.null(xgb$early_stopping_rounds))
    nrounds <- earlyStoppedRounds(x, y, evalSet, weight, xgb) |>
      Cache(omitArgs = c("x", "y", "weight", "evalSet"),
            .functionName = paste0(functionName, "_rounds"),
            .cacheExtra = c(cacheExtra, list(xgb = xgb)),
            cacheSaveFormat = "rds")
  xgboost(x = x, y = y, weights = weight,
          objective = fireSenseUtils::ignitionObjectiveTxt,
          nthread = fireSenseUtils::nThreads(),
          eval_set = evalSet,
          monitor_training = TRUE,
          verbosity = 0,
          nrounds = nrounds,
          max_depth = xgb$max_depth,
          reg_lambda = xgb$reg_lambda,
          learning_rate = xgb$learning_rate
  ) |> Cache(omitArgs = c("x", "y", "weights", "eval_set", "nthread"),
             .functionName = functionName,
             .cacheExtra = cacheExtra,
             cacheSaveFormat = "rds")
}

#' The number of boosting rounds that early stopping picks for one fold
#'
#' `xgboost()` takes `early_stopping_rounds` but does not pass it on (xgboost 3.2.1.1), so the rounds are chosen
#' here with `xgb.train()`: boosting runs to `xgb$nrounds`, stopping once the Poisson loss of the held-out rows
#' (`evalSet`) has not improved for `xgb$early_stopping_rounds` rounds. Training is deterministic, so
#' `fitFoldXgboost()` refits with the best round as `nrounds` and gets that round's model.
#'
#' @inheritParams fitFoldXgboost
#' @return The best round: the one with the lowest held-out Poisson loss (a number from 1 to `xgb$nrounds`).
earlyStoppedRounds <- function(x, y, evalSet, weight = NULL, xgb = xgbSettings()) {
  evalSet <- as.integer(evalSet)
  dm <- xgb.DMatrix(as.matrix(x), label = y, weight = weight)
  trained <- xgb.train(
    params = list(objective = fireSenseUtils::ignitionObjectiveTxt, max_depth = xgb$max_depth,
                  learning_rate = xgb$learning_rate, reg_lambda = xgb$reg_lambda,
                  nthread = fireSenseUtils::nThreads(), verbosity = 0),
    data = xgb.slice.DMatrix(dm, setdiff(seq_len(nrow(x)), evalSet)),
    nrounds = xgb$nrounds, evals = list(eval = xgb.slice.DMatrix(dm, evalSet)),
    early_stopping_rounds = xgb$early_stopping_rounds, verbose = 0)
  which.min(attributes(trained)$evaluation_log[[2]])
}

#' Message reporting the mean AUC across folds
#'
#' @param aucs Output of `aucPerFold()`.
#' @return Character string.
meanRocMessage <- function(aucs) {
  if (all(is.na(aucs)))
    return("mean roc:  not computed (no validation fold had both outcomes)")
  paste0("mean roc:  ", format(mean(aucs, na.rm = TRUE), digits = 3),
         if (anyNA(aucs)) paste0(" (over ", sum(!is.na(aucs)), " of ", length(aucs),
                                 " folds)") else "")
}

#' Paste with `_`, for cache function names and filenames
#'
#' @param ... Passed to `paste()`.
#' @param sep Separator.
#' @return Character string.
functionNameHelper <- function(..., sep = "_") {
  paste(..., sep = sep)
}

#' Plot the predicted response against each covariate of one group
#'
#' @param df data.table made in `setupPlots()`, with columns `value`, `predictedProb`, `varFac`.
#' @param fuelOrClimate `"Fuel"` or `"Climate"`: the group of covariates to plot.
#' @param labels Named character vector giving the group of each covariate.
#' @param value Unused; `value` in the plot is the column of `df`.
#' @param jitter Jitter width on x, so the folds, which share x values, do not overlap exactly. y is never jittered.
#' @param colors Named vector of colours, one per covariate.
#' @param fuelOrClimateInd Index into `colors` and `labels` of the covariates in this group.
#' @param igOrEsc `"ignition"` or `"escape"`; used in the y-axis label. The models are Poisson
#'   (`count:poisson`) fits to counts, so the y-axis is an expected count: ignitions per pixel-year,
#'   or escapes per ignited pixel-year.
#' @return A ggplot.
plotPredictions <- function(df, fuelOrClimate, labels, value, jitter, colors, fuelOrClimateInd, igOrEsc) {

  ggplot(df[varFac %in% names(labels)[labels %in% fuelOrClimate]],
         aes(x = value, y = predictedProb, group = varFac, col = varFac)) +
    geom_point(position = position_jitter(width = jitter, height = 0)) +
    geom_smooth(span = 1) +
    scale_color_manual(aesthetics = "colour", values = colors[fuelOrClimateInd],
                       labels = names(labels)[fuelOrClimateInd],
                       name = paste0(fuelOrClimate, ":"),
                       guide = guide_legend(reverse = TRUE, title.position = "top", order = 0)) +
    ggplot2::xlab(paste0("Scaled, centred")) +
    ggplot2::ylab(paste0("Predicted ", igOrEsc, "s per ", if (igOrEsc == "escape") "ignited ", "pixel-year")) +
    theme_bw()
}

#' Build the response plots for the fitted xgboost models
#'
#' For each fold and covariate, predicts along the covariate's whole observed range (minimum to
#' maximum of the standardised data, `nGrid` equally spaced points, so the cost does not grow with
#' the length of the covariate's tail) with every other covariate at 0, its mean.
#'
#' @param modelOnly List of the per-fold xgboost models.
#' @param dat data.table of standardised covariates used for the fit.
#' @param igOrEsc `"ignition"` or `"escape"`.
#' @param nGrid Number of points per covariate at which to predict.
#' @return List of two ggplots, `Climate` and `Fuel`. Covariates whose names match
#'   `CMD|light|positiveCG` are the climate ones.
setupPlots <- function(modelOnly, dat, igOrEsc, nGrid = 50) {

  cnNoIgnNoEsc <- colnames(dat) |> setdiff(c(fireSenseUtils::ignitionsTxt, fireSenseUtils::escapesTxt, "year", "pixelID", fireSenseUtils::coreFractionTxt))
  dat <- dat[, ..cnNoIgnNoEsc]

  df <- Map(fold = seq_along(modelOnly), function(fold) {
    df <- Map(nam = cnNoIgnNoEsc, function(nam) {
      rr <- range(dat[[nam]], na.rm = TRUE)
      df1 <- data.frame(seq(rr[1], rr[2], length.out = nGrid)) |> setNames(nam)
      cnHere <- setdiff(cnNoIgnNoEsc, nam)
      df0 <- lapply(cnHere, function(x) list(0)) |> data.frame() |> setNames(cnHere)
      df <- data.frame(df1, df0)
      # Predict along the regularl sequence of x-axis values
      predictedProb <- predict(modelOnly[[fold]], newdata = df)
      data.table(df1, predictedProb) |>
        melt(measure.vars = nam, stringsAsFactor = FALSE)
    })
    rbindlist(df, use.names = TRUE)
  })
  df <- rbindlist(df, use.names = TRUE)
  set(df, NULL, "variable", as.character(df$variable))
  setorderv(df, "variable")

  set(df, NULL, "varFac", factor(df$variable, levels = cnNoIgnNoEsc)) # keeps colours constant
  set(df, NULL, "varInt", as.integer(df$varFac))

  levels(df$varFac) # gives order that ggplot2 will use
  setorderv(df, "varFac")

  vals <- cnNoIgnNoEsc # stays constant colour regardless of importance
  colors <- NULL
  colOptions <- c("Paired", "Gr")
  colors <- RColorBrewer::brewer.pal(length(vals), "Paired")
  if (length(vals) != length(colors))
    colors <- colorRampPalette(colors)(length(vals))
  names(colors) <- vals
  labels <- rep("Fuel", length(vals))
  names(labels) <- vals
  climateGrep <- "CMD|light|positiveCG"
  climateInd <- grep(climateGrep, vals) # this is how I identify climate vars: not robust!!!!!
  set(df, NULL, "FuelOrClimate", "Fuel")
  set(df, which(df$variable %in% names(labels[climateInd])), "FuelOrClimate", "Climate")

  labels[climateInd] <- "Climate"
  FuelInd <- -climateInd

  jitter <- 0.05

  aa <- Map(foc = c("Climate", "Fuel"), focInd = list(climateInd, FuelInd), function(foc, focInd) {
    plotPredictions(df, fuelOrClimate = foc, labels, value,
                    jitter, colors, fuelOrClimateInd = focInd, igOrEsc)
  })
  aa
}

