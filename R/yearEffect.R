#' Observed and out-of-fold predicted totals per year
#'
#' Each row is predicted by the one fold model that did not train on it (its fold's `keepEval`), so the
#' predicted totals are out of fold. With random-row folds (`cvFolds()`) every year is still in every
#' model's training rows, so the year totals track the observed ones too closely and understate the year
#' effect. `runXGBOOST()` therefore gives this year-blocked folds (`yearBlockedFolds()`): the totals are
#' out of year, and `fireSenseUtils::yearEffectSD()` of them is the fit's `yearSD`.
#'
#' @param year The year of each row (as fitted, so standardised when `rescaleVars`); `NULL` if the data
#'   had no year column.
#' @param observed The response of each row.
#' @param folds Output of `cvFolds()` or `yearBlockedFolds()`.
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
