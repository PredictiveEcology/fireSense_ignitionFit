## A project that sets no `.useCache` gets the module's default. It must cache what a project script
## used to switch on with `fireSense_ignitionFit = list(.useCache = TRUE)`: every event and
## `.inputObjects`, except what `.neverCache` removes (`checkData`).

md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
paramDefault <- function(name) md$parameters$default[[which(md$parameters$paramName == name)]]
code <- readLines(file.path(moduleRoot, paste0(moduleName, ".R")))
eventIsCached <- utils::getFromNamespace(".eventIsCached", "SpaDES.core")
defaultParams <- list(.useCache = paramDefault(".useCache"), .neverCache = paramDefault(".neverCache"))

test_that("the default .useCache caches every event and .inputObjects except those in .neverCache", {
  eventsScheduled <- unique(sub('.*moduleName, "([[:alnum:]]+)".*', "\\1", grep('scheduleEvent\\(.*moduleName, "', code, value = TRUE)))
  expect_true(any(grepl("^\\.inputObjects <- function", code)))
  steps <- c(".inputObjects", "init", eventsScheduled)
  cached <- steps[vapply(steps, function(s) eventIsCached(defaultParams, s)$cache, logical(1))]
  expect_setequal(cached, setdiff(steps, paramDefault(".neverCache")))
  expect_false("checkData" %in% cached)
})

test_that(".useCache accepts logical or character", {
  expect_setequal(unlist(md$parameters$paramClass[md$parameters$paramName == ".useCache"]), c("logical", "character"))
})
