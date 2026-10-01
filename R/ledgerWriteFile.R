#' The ledger file this run reads and writes
#'
#' `"latest"` is resolved deterministically from this run's own fire years, so the file this
#' module reads for `sim$ignitionFitPreRun` (in `Init()`) is exactly the file it would write to
#' (in `writeIgnitionLedgerRow()`) -- unlike `fireSense_dataPrepFit`'s `spreadFitFilename`, which
#' reads a fixed name because that module has no fit of its own to name the file after. Any other
#' `ignitionFitFilename` is used as it is. Parallel to fireSense_SpreadFit's `ledgerWriteFile()`.
#'
#' @param ignitionFitFilename The module's `ignitionFitFilename` parameter.
#' @param fireYears The years of the covariates this run fit on (see [ignitionFitYears()]).
#' @return A file name.
ignitionLedgerWriteFile <- function(ignitionFitFilename, fireYears) {
  if (!identical(ignitionFitFilename, "latest"))
    return(ignitionFitFilename)
  fireSenseUtils::ignitionFitFilenameFor(fireYears)
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
