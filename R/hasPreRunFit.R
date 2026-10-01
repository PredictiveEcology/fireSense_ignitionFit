#' Does the shared ledger already hold a fit for the polygon this run fits?
#'
#' `sim$ignitionFitPreRun` is whatever ledger rows intersect `sim$studyArea`: neighbours' fits,
#' this polygon's fit, or nothing at all. Its class therefore says nothing about *this* polygon.
#' Parallel to fireSense_SpreadFit's `hasPreRunFitForThisPolygon()`.
#'
#' @param sim A `simList`. Uses `sim$.ELFind` as the polygon identity (the ledger key).
#' @return Logical, length one. Always `FALSE` when `sim$studyArea` was not supplied (no
#'   `sim$ignitionFitPreRun` was read).
hasPreRunIgnitionFitForThisPolygon <- function(sim) {
  sa <- sim$ignitionFitPreRun
  if (!(is(sa, "sf") || is(sa, "data.frame")) || NROW(sa) == 0L) return(FALSE)
  if (!"polygonID" %in% names(sa)) return(FALSE)
  id <- sim$.ELFind
  if (is.null(id) || !length(id)) return(FALSE)
  isTRUE(any(as.character(sa$polygonID) == as.character(id[1])))
}

#' Take this polygon's ignition and escape fits from the ledger instead of fitting
#'
#' @param sim A `simList`; `sim$ignitionFitPreRun` must already have a row for `sim$.ELFind`
#'   (see [hasPreRunIgnitionFitForThisPolygon()]).
#' @return `sim`, with `sim$fireSense_IgnitionFitted` and `sim$fireSense_EscapeFitted` set from
#'   the matching ledger row.
useExistingIgnitionFit <- function(sim) {
  row <- sim$ignitionFitPreRun[which(as.character(sim$ignitionFitPreRun$polygonID) ==
                                      as.character(sim$.ELFind))[1], ]
  message("fireSense_IgnitionFit: ledger already holds ignition/escape fits for polygon ",
          sim$.ELFind, "; using those instead of fitting.")
  sim$fireSense_IgnitionFitted <- row[[fireSenseUtils::ignitionFitAdditionalColNamesTxt[1]]][[1]]
  sim$fireSense_EscapeFitted <- row[[fireSenseUtils::ignitionFitAdditionalColNamesTxt[2]]][[1]]
  invisible(sim)
}
