# fireSense_ignitionFit (development version)

- The response-curve figures are saved next to the fit's ledger file, in `<inputPath>/fits/<.ELFind>_<ledger file name>` (`fireSenseUtils::fitOutputPath()`), not in `figurePath(sim)`, which in FireSense runs is the output folder of whichever scenario and replicate fitted first. A fit depends only on its polygon, fire years and model. New parameter `fitOutputPath` (default `NULL`) overrides the folder. Floor: `fireSenseUtils` (>= 0.2.3.9087, for `fitOutputPath()`).
- The response-curve figures now cover each covariate's whole observed range (50 points from its minimum to its maximum), not just -2 to 2 standardised units. The fuel panel's legend says "Fuel:" (it said "Climate:"), and the y axis reads "Predicted ignitions per pixel-year" or "Predicted escapes per ignited pixel-year" rather than "probability", because the Tweedie models predict a count. Each fold's prediction is drawn once, with x jitter only (the folds share x values); the y jitter, which was large relative to ignition rates near 3e-5, is gone.
- A fit stored in the ledger is reused only if its features (the xgboost boosters' `variable.names()`) are all among the current covariate columns, for each process in `whichProcessesToFit`. Otherwise the module says which features are missing and fits again, replacing the ledger row, as with `refitExisting = TRUE`. Before, a stored fit made with older covariates (e.g. non-forest groups since merged) was reused and `fireSense_ignitionPredict` failed on the missing columns. `runXGBOOST()` and the check share `xgbFeatureNames()`.

# fireSense_ignitionFit 1.1.2

- Ignition and escape fits are shared through a cloud, geo-keyed ledger, as `fireSense_spreadFit` does its fits. When `studyArea` and `.ELFind` are supplied, the module reads the ledger (new parameters `ignitionFitGoogleDriveFolder` and `ignitionFitFilename`) for the study area. If a row exists for this polygon, its `fireSense_IgnitionFitted` and `fireSense_EscapeFitted` are used and the fit is skipped. Otherwise the fit runs and its result is written to the ledger. New parameter `refitExisting` fits even when the ledger has a row. Without `studyArea` nothing changes. New output `ignitionFitPreRun`, the ledger rows read for `studyArea`.
- The ledger row drops `modelList$model$rocs` (the per-fold `pROC::roc()` curves, unused by `predict()`), so rows stay small.
- Floors: `fireSenseUtils` (>= 0.2.3.9078, for `ignitionFitFilenameFor()` and `ignitionFitAdditionalColNamesTxt`), `reproducible` (>= 3.2.1.9058, whose `CacheGeo()` can append to a ledger that holds xgboost models).

# fireSense_ignitionFit 1.1.1

- reqdPkgs now lists `purrr` and `RColorBrewer`, which the module calls with `::` but did not list. Version 1.1.1.

# fireSense_ignitionFit 1.1.0

- Renamed from `fireSense_IgnitionFit` to `fireSense_ignitionFit` (module naming convention `<model>_<camelCaseComponent>`); projects must rename the module and its `params` key. The class `fireSense_IgnitionFit` and the objects `fireSense_IgnitionFitted`, `fireSense_IgnitionFittedList` and `fireSense_IgnitionPredicted` keep their names. Version 1.1.0.

- New parameter `.studyAreaName` (default `NA`), the name PredictiveEcology modules use for the study area. This module does not use it yet.

# fireSense_IgnitionFit 1.0.2

## Breaking changes

- Removed the non-xgboost modelling path entirely (`runGLMMAdaptiveWithSimplifications()`,
  `runGLM.NB()`, `checkData()`, `messageFormulaFn()`, and the formula/AIC code in `buildModel()`),
  together with the disabled time-ordered branch of `runXGBOOST()`.
  `buildModel()` now stops unless `modelAlgorithm` contains `"xgb"`.
- Removed defunct parameters: `escapeFamily`, `ignitionFamily`, `plot_fuelBiomassPerPrediction`,
  `.studyAreaName`, `.plotInitialTime`, `.saveInitialTime`, `.saveInterval`.
- Removed defunct inputs: `fireSense_ignitionFormula`, `climateVariablesForFire`.
- `modelList` no longer carries a `family` element: it was always the `stats::family` function
  itself, never a family object, and nothing read it.

## Bug fixes

- The digest of the covariates is now actually passed to `runXGBOOST()`. The module asked
  `prepareCovariatesOuter()` for `digestOfData$fireSense_ignitionCovariates`, which is not an
  element of what that function returns, so `NULL` was passed and the data played no part in the
  cache key of the per-fold models. **This changes the cache keys: existing caches of the fit are
  invalidated and the models will be refitted once.**

# fireSense_IgnitionFit 1.0.1

First release from `development` since `master` was last updated (2020-04-02). Full history: https://github.com/PredictiveEcology/fireSense_IgnitionFit/compare/8ea03bc...v1.0.1

## Breaking changes

- Removed input `dataFireSense_IgnitionFit` (data.frame).
- Removed parameters: `cores`, `data`, `family`, `formula`, `iterDEoptim`, `iterNlminb`, `lb`, `nlminb.control`, `plot`, `start`, `trace`, `ub`.

## New features

- New inputs: `climateVariablesForFire`, `fireSense_ignitionCovariates`, `fireSense_ignitionFormula`, `ignitionFitRTM`.
- New output: `fireSense_EscapeFitted`.
- New parameters: `.plotInitialTime`, `.plots`, `.seed`, `.studyAreaName`, `crossValType`, `escapeFamily`, `ignitionFamily`, `modelAlgorithm`, `plot_fuelBiomassPerPrediction`, `rescaleVars`, `whichProcessesToFit`.

## Dependencies

- No longer depends on `DEoptim`.
- Now depends on `RhpcBLASctl`, `SHAPforxgboost`, `SpaDES.core`, `caret`, `data.table`, `fireSenseUtils`, `ggplot2`, `ggpubr`, `glmmTMB`, `lightgbm`, `pROC`, `parallelly`, `pemisc`, `reproducible`, `terra`, `xgboost`.

## Testing

- testthat suite and CI (`testthat-module`), including a snapshot of the module's inputs, outputs and parameters in `tests/testthat/test-metadata.R`.
