#' The settings of the xgboost fits
#'
#' The one place the boosting settings of the ignition and the escape fits are put together, from the module's
#' `xgb*` parameters; the defaults are the values the fits were hard-coded with before (`nrounds = 100`,
#' `max_depth = 2`, `learning_rate = 0.15`, `reg_lambda = 0.5`, no early stopping), so a fit with the default
#' parameters is the one made before they existed. `runXGBOOST()` and `fitFoldXgboost()` take this list as it
#' is, so a setting is never listed twice.
#'
#' @param maxDepth,learningRate,regLambda xgboost's `max_depth`, `learning_rate` and `reg_lambda`.
#' @param nrounds xgboost's `nrounds`: the number of boosting rounds, or, with `earlyStopping`, their cap.
#' @param earlyStopping `NA` (or `NULL`, `0`): off. Otherwise xgboost's `early_stopping_rounds`: each fold's
#'   boosting stops once the held-out fold's Poisson loss (`poisson-nloglik`, the deviance up to a constant) has
#'   not improved for this many rounds, and the fold's model keeps the best round.
#' @return A named list: `max_depth`, `learning_rate`, `reg_lambda`, `nrounds` and, only with early
#'   stopping, `early_stopping_rounds` (so an unset one adds nothing to the fit's cache key).
xgbSettings <- function(maxDepth = 2, learningRate = 0.15, nrounds = 100, regLambda = 0.5, earlyStopping = NA) {
  settings <- list(nrounds = nrounds, max_depth = maxDepth, reg_lambda = regLambda, learning_rate = learningRate)
  if (length(earlyStopping) && !is.na(earlyStopping) && earlyStopping > 0)
    settings$early_stopping_rounds <- as.integer(earlyStopping)
  settings
}

#' The `xgbSettings()` of the module's parameters
#'
#' @param par The module's parameters, `P(sim)` or `Par`.
#' @return `xgbSettings()` of `xgbMaxDepth`, `xgbLearningRate`, `xgbNrounds`, `xgbRegLambda` and `xgbEarlyStopping`.
xgbSettingsOfParams <- function(par) {
  xgbSettings(maxDepth = par$xgbMaxDepth, learningRate = par$xgbLearningRate, nrounds = par$xgbNrounds,
              regLambda = par$xgbRegLambda, earlyStopping = par$xgbEarlyStopping)
}
