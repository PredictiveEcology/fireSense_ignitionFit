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
  row <- ledgerRowForThisPolygon(sim)
  message("fireSense_ignitionFit: ledger already holds ignition/escape fits for polygon ",
          sim$.ELFind, "; using those instead of fitting.")
  sim$fireSense_IgnitionFitted <- row[[fireSenseUtils::ignitionFitAdditionalColNamesTxt[1]]][[1]]
  sim$fireSense_EscapeFitted <- row[[fireSenseUtils::ignitionFitAdditionalColNamesTxt[2]]][[1]]
  invisible(sim)
}

#' The ledger row for the polygon this run fits
#'
#' @param sim A `simList` for which [hasPreRunIgnitionFitForThisPolygon()] is `TRUE`.
#' @return The one row of `sim$ignitionFitPreRun` whose `polygonID` is `sim$.ELFind`.
ledgerRowForThisPolygon <- function(sim) {
  sim$ignitionFitPreRun[which(as.character(sim$ignitionFitPreRun$polygonID) ==
                                as.character(sim$.ELFind))[1], ]
}

#' Names of the features `runXGBOOST()` gives the model
#'
#' Every column except `pixelID`, `year`, `yearChar*` and the responses, in alphabetical order.
#' `runXGBOOST()` and [storedFitIsCompatible()] both use it, so they cannot disagree.
#'
#' @param cols Character vector of the column names of the fit data.
#' @return Character vector.
xgbFeatureNames <- function(cols) {
  cols <- setdiff(cols, c("pixelID", "year"))
  cols <- grep("yearChar", cols, invert = TRUE, value = TRUE)
  sort(grep(paste0(fireSenseUtils::ignitionsTxt, "|", fireSenseUtils::escapesTxt), cols,
            value = TRUE, invert = TRUE))
}

#' Can a stored fit predict with the covariates this run would fit with?
#'
#' A fit stored in the ledger by an earlier run records the features it was trained on (its
#' xgboost boosters' `variable.names()`). If the covariates have since changed (e.g. non-forest
#' land-cover groups merged into one), `predict()` fails on the missing columns. The stored fit is
#' compatible when all its features are among the current covariate columns. Only xgboost fits
#' exist (`buildModel()` stops for any other algorithm). Messages name what is missing and extra.
#'
#' @param fitted A stored `list(modelList = list(model = ...), ...)`, or `NULL`.
#' @param covariates The data.table `buildModelsFitModels()` would fit with.
#' @param polygonID,type Used in the message only.
#' @return Logical, length one; `FALSE` if `fitted` is `NULL` or has no boosters.
storedFitIsCompatible <- function(fitted, covariates, polygonID = "", type = "ignition") {
  boosters <- fitted$modelList$model
  boosters <- boosters[grep("^Fold", names(boosters))]
  if (!length(boosters)) {
    message("fireSense_ignitionFit: no stored ", type, " model for polygon ", polygonID, "; fitting.")
    return(FALSE)
  }
  stored <- unique(unlist(lapply(boosters, stats::variable.names)))
  current <- xgbFeatureNames(colnames(covariates))
  missingNow <- setdiff(stored, current)
  if (length(missingNow)) {
    message("fireSense_ignitionFit: the stored ", type, " model for polygon ", polygonID,
            " uses features absent from the current covariates (", paste(missingNow, collapse = ", "),
            "); current covariates not in the model: ",
            paste(setdiff(current, stored), collapse = ", "), ". Fitting again.")
    return(FALSE)
  }
  TRUE
}

#' Is every stored fit this run needs compatible with the current covariates?
#'
#' @param sim A `simList` for which [hasPreRunIgnitionFitForThisPolygon()] is `TRUE`.
#' @return Logical, length one: `storedFitIsCompatible()` for each of `P(sim)$whichProcessesToFit`.
useExistingFitIsCompatible <- function(sim) {
  row <- ledgerRowForThisPolygon(sim)
  procs <- P(sim)$whichProcessesToFit
  ok <- vapply(procs, function(igOrEsc) {
    col <- fireSenseUtils::ignitionFitAdditionalColNamesTxt[match(igOrEsc, c("ignition", "escape"))]
    storedFitIsCompatible(row[[col]][[1]],
                          sim[[igOrEscNames(igOrEsc, post = "Covariates")]],
                          polygonID = sim$.ELFind, type = igOrEsc)
  }, logical(1))
  all(ok)
}
