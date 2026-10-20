#' Observed and out-of-fold predicted totals per year
#'
#' Each row is predicted by the one fold model that did not train on it (its fold's `keepEval`), so the
#' predicted totals are out of fold. The folds are random rows (`cvFolds()`), so every year is still in
#' every model's training rows: the totals are less optimistic than in-sample ones, not out of year.
#' `fireSenseUtils::yearEffectSD()` of these is the fit's `yearSD`.
#'
#' @param year The year of each row (as fitted, so standardised when `rescaleVars`); `NULL` if the data
#'   had no year column.
#' @param observed The response of each row.
#' @param folds Output of `cvFolds()`.
#' @param heldOutPreds List, one element per fold: the fold model's predictions on its `keepEval` rows.
#' @return A `data.frame` with `year`, `observed` and `predicted`, one row per year, or `NULL` without `year`.
outOfFoldYearTotals <- function(year, observed, folds, heldOutPreds) {
  if (is.null(year)) return(NULL)
  pred <- rep(NA_real_, length(observed))
  for (k in seq_along(folds)) pred[folds[[k]]$keepEval] <- heldOutPreds[[k]]
  tot <- rowsum(cbind(observed = as.numeric(observed), predicted = pred), group = year)
  data.frame(year = as.numeric(rownames(tot)), observed = tot[, "observed"], predicted = tot[, "predicted"],
             row.names = NULL)
}
