## outOfFoldYearTotals(): the yearly totals fireSenseUtils::yearEffectSD() estimates the year effect from

test_that("outOfFoldYearTotals sums each row's held-out prediction by year", {
  folds <- list(list(keepAll = 1:6, keepEval = c(1L, 4L)), list(keepAll = 1:6, keepEval = c(2L, 3L, 5L, 6L)))
  yt <- outOfFoldYearTotals(year = c(1, 1, 1, 2, 2, 2), observed = c(1, 0, 2, 0, 0, 3), folds = folds,
                            heldOutPreds = list(c(0.5, 0.25), c(1, 2, 4, 8)))
  expect_identical(yt, data.frame(year = c(1, 2), observed = c(3, 3), predicted = c(0.5 + 1 + 2, 0.25 + 4 + 8)))
  expect_null(outOfFoldYearTotals(NULL, 1, folds, list()))
})
