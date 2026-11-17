#' The ledger file this run reads and writes
#'
#' `"latest"` is resolved deterministically from this run's own fire years, so the file this
#' module reads for `sim$ignitionFitPreRun` (in `Init()`) is exactly the file it would write to
#' (in `writeIgnitionLedgerRow()`) -- unlike `fireSense_dataPrepFit`'s `spreadFitFilename`, which
#' reads a fixed name because that module has no fit of its own to name the file after. Any other
#' `ignitionFitFilename` is used as it is, and `fitFileSuffix` does not apply. Parallel to fireSense_SpreadFit's `ledgerWriteFile()`.
#'
#' @param ignitionFitFilename The module's `ignitionFitFilename` parameter.
#' @param fitFileSuffix The `fitFileSuffix` parameter, agreed with the other fit modules
#'   (`SpaDES.core::paramCheckOtherMods()`).
#' @param fireYears The years of the covariates this run fit on (see [ignitionFitYears()]).
#' @return A file name.
ignitionLedgerWriteFile <- function(ignitionFitFilename, fireYears, fitFileSuffix = "") {
  if (!identical(ignitionFitFilename, "latest"))
    return(ignitionFitFilename)
  fireSenseUtils::ignitionFitFilenameFor(fireYears, fitFileSuffix)
}

#' The fire years this run's covariates cover
#'
#' Used to resolve `ignitionFitFilename = "latest"` to a concrete file name. Reads whichever of
#' `sim$fireSense_ignitionCovariates`/`sim$fireSense_escapeCovariates` is present -- one may be
#' absent depending on `whichProcessesToFit`.
#'
#' @param sim A `simList`.
#' @return Integer vector of years.
ignitionFitYears <- function(sim) {
  yrs <- c(sim$fireSense_ignitionCovariates[[fireSenseUtils::yearTxt]],
           sim$fireSense_escapeCovariates[[fireSenseUtils::yearTxt]])
  as.integer(yrs)
}

#' The ledger rows of this polygon, when the run has no fit data to name the file from
#'
#' With `ignitionFitFilename = "latest"` the file is named for the fit's fire years, which come from the
#' covariates ([ignitionFitYears()]). A run that only reuses the ledger's fit has none (`fireSense_dataPrepFit`
#' does not prepare them when `doFit` is `FALSE`), so the newest row of `sim$.ELFind` over the ledger files is
#' read, as `fireSense_ELFs` reads it (`fireSenseUtils::latestIgnitionFits()`).
#'
#' @param sim A `simList`.
#' @return The ledger rows of `sim$.ELFind`, or `NULL` when no file has one.
ignitionLedgerLatestRows <- function(sim) {
  if (is.null(Par$ignitionFitGoogleDriveFolder))
    stop("fireSense_ignitionFit: ignitionFitFilename = \"latest\" without fit data needs ignitionFitGoogleDriveFolder; ",
         "name the ledger file to read a local one", call. = FALSE)
  fireSenseUtils::latestIgnitionFits(Par$ignitionFitGoogleDriveFolder, destinationPath = inputPath(sim),
                                     polygonIDs = as.character(sim$.ELFind),
                                     suffix = SpaDES.core::paramCheckOtherMods(sim, "fitFileSuffix"))
}
