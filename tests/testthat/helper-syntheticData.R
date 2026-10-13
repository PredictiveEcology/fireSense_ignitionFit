## Synthetic inputs for the tests. No random numbers: every value is arithmetic on the row
## index, so the data are identical on every machine and R version.

## 100 pixels x 4 years = 400 rows. `ignitions` is a function of `CMDsm` and `youngAge` only, so a
## model that sees those two columns can separate the outcomes perfectly (AUC = 1). `escapes` is a
## function of `lightning` among the rows with an ignition, for the same reason.
makeIgnitionCovariates <- function(escapes = FALSE) {
  i <- 1:400
  CMDsm <- ((i * 37) %% 101) / 10                  # 0 .. 10
  dt <- data.table::data.table(
    pixelID   = rep(1:100, times = 4),
    year      = rep(2001:2004, each = 100),
    CMDsm     = CMDsm,
    youngAge  = as.numeric(i %% 2),
    Pice_mar  = as.numeric((i * 53) %% 97),
    lightning = as.numeric((i * 17) %% 13),
    ignitions = as.integer(CMDsm > 8 & i %% 2 == 0) + as.integer(CMDsm > 8.9)
  )
  if (escapes)
    data.table::set(dt, NULL, "escapes", as.integer(dt$ignitions > 0 & dt$lightning >= 6))
  dt
}

## The covariates as `fireSenseUtils::rescaleCovariates()` standardises them for xgboost: every
## column except the responses is centred and scaled.
scaleCovariates <- function(dt) {
  out <- data.table::copy(dt)
  for (col in setdiff(names(out), c("ignitions", "escapes")))
    data.table::set(out, NULL, col, (out[[col]] - mean(out[[col]])) / stats::sd(out[[col]]))
  out
}

## 10 x 10 template raster, 250 m cells. `nonNAs` is set to 1000, not 100, so that
## `lambdaRescaleFactor` (= rows of covariates / nonNAs) is not 1.
makeIgnitionFitRTM <- function(nonNAs = 1000) {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 2500, ymin = 0, ymax = 2500,
                   crs = "EPSG:3005", vals = 1)
  attr(r, "nonNAs") <- nonNAs
  r
}

## 2000 pixels x 10 years = 20000 rows with 33 ignitions, all where `x1` is high (0.17% of rows).
## Sparse enough that a log-link Tweedie objective loses predicted mass every boosting round.
makeSparseCovariates <- function(nPix = 2000, nYr = 10) {
  i <- seq_len(nPix * nYr)
  x1 <- ((i * 37) %% 1009) / 100.9
  dt <- data.table::data.table(
    pixelID = rep(seq_len(nPix), times = nYr),
    year    = rep(2001:(2000 + nYr), each = nPix),
    x1      = x1,
    x2      = as.numeric((i * 53) %% 97),
    x3      = as.numeric((i * 17) %% 13),
    ignitions = as.integer(x1 > 9.5 & i %% 25 == 0)
  )
  dt
}
