#' Strip a fitted model down to what `predict()` needs before it goes in the shared ledger
#'
#' `modelList$model` holds one xgboost `Booster` per fold plus `rocs`, the `pROC::roc()` curves
#' `runXGBOOST()` builds from each fold's held-out response and prediction, purely to print and
#' return the mean AUC (`meanRocMessage()`). Nothing downstream -- `fireSense_IgnitionPredict`'s
#' `predict()` calls, or `useExistingIgnitionFit()` reading this same object back -- reads `rocs`;
#' it only inflates the ledger row. `test-writeIgnitionLedgerRow.R` checks that `predict()` still
#' works on the stripped object.
#'
#' @param fitted One `list(modelList = list(model = ...), scaleData = ...)`, the shape of
#'   `sim$fireSense_IgnitionFitted`/`sim$fireSense_EscapeFitted`. `NULL` is returned as is (a
#'   process that was not fit, e.g. `whichProcessesToFit = "ignition"`, has no `EscapeFitted`).
#' @return The same shape, with `modelList$model$rocs` removed.
stripFittedForLedger <- function(fitted) {
  if (is.null(fitted)) return(fitted)
  fitted$modelList$model$rocs <- NULL
  fitted
}

#' Write this run's ignition and escape fits to the shared ledger
#'
#' One row per polygon: `sim$studyArea`'s geometry and crs, `polygonID` (`sim$.ELFind`, validated
#' exactly as fireSense_SpreadFit validates it -- the ledger is shared cloud state every other
#' project reads, keyed on the polygon's identity, not a run label), and the stripped
#' (`stripFittedForLedger()`) `fireSense_IgnitionFitted`/`fireSense_EscapeFitted`. Parallel to
#' fireSense_SpreadFit's ledger write in its `run` event.
#'
#' @param sim A `simList`.
#' @return `sim`, with `sim$ignitionFitPreRun` updated to what `CacheGeo()` returns.
writeIgnitionLedgerRow <- function(sim) {
  polygonID <- sim$.ELFind
  if (!is.character(polygonID) || length(polygonID) != 1L ||
      is.na(polygonID) || !nzchar(polygonID))
    stop("fireSense_ignitionFit: `sim$.ELFind` must be a single non-empty character ",
         "identifying the polygon being fit; got: ",
         paste(format(polygonID), collapse = ", "))

  df <- data.frame(I(list(stripFittedForLedger(sim$fireSense_IgnitionFitted))),
                    I(list(stripFittedForLedger(sim$fireSense_EscapeFitted)))) |>
    setNames(fireSenseUtils::ignitionFitAdditionalColNamesTxt)
  df <- data.frame(df, polygonID = polygonID)

  crses <- terra::crs(sim$studyArea)
  df <- dplyr::mutate(df, crs = I(crses))

  saHere <- if (is(sim$studyArea, "SpatVector")) sf::st_as_sf(sim$studyArea) else sim$studyArea
  saHere <- sf::st_as_sf(sf::st_geometry(saHere))
  sf::st_geometry(saHere) <- "geometry"
  ledgerRow <- dplyr::mutate(saHere, df)

  le <- function(x) {x}
  sim$ignitionFitPreRun <- ignitionLedgerCacheGeo(sim, domain = saHere, FUN = le(ledgerRow),
                                                  le = le, ledgerRow = ledgerRow, action = "update")
  invisible(sim)
}

#' The one place the ignition/escape ledger is read or written through `CacheGeo()`
#'
#' Both the read in `Init()` and the write in `writeIgnitionLedgerRow()` come here, so moving the
#' ledger to another reproducible API changes this function only.
#'
#' @param sim A `simList`.
#' @param ... Passed on to `CacheGeo()`: `domain`, `action`, and, for a write, `FUN`.
#' @return What `CacheGeo()` returns.
ignitionLedgerCacheGeo <- function(sim, ...) {
  CacheGeo(
    cloudFolderID = Par$ignitionFitGoogleDriveFolder,
    targetFile = ignitionLedgerWriteFile(Par$ignitionFitFilename, ignitionFitYears(sim)),
    destinationPath = inputPath(sim),
    ## `purge` re-downloads from the cloud copy; with no cloud folder there is nothing to
    ## re-download, and `prepInputs()` then cannot find the local file it just set aside.
    purge = if (is.null(Par$ignitionFitGoogleDriveFolder)) FALSE else 7,
    ...)
}
