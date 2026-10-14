## `.neverCache` (read by SpaDES.core >= 3.2.1.9037) wins over every form of `.useCache`, including a
## `.globals$.useCache` that lists the event. `checkData` reads the ignition-fit ledger on Google Drive, and `run` decides from it whether to fit and then writes the ledger. Neither is in the event's cache key, and a cache hit would skip the write.
test_that(".neverCache is declared and SpaDES.core refuses to cache those events whatever .useCache says", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  declared <- md$parameters$default[[which(md$parameters$paramName == ".neverCache")]]
  expect_setequal(declared, c("checkData", "run"))
  eventIsCached <- utils::getFromNamespace(".eventIsCached", "SpaDES.core")
  for (useCache in list(c(".inputObjects", declared), TRUE)) {
    for (ev in declared)
      expect_false(eventIsCached(list(.useCache = useCache, .neverCache = declared), ev)$cache)
  }
  ## and the helper still caches what `.useCache` names that `.neverCache` does not list
  expect_true(eventIsCached(list(.useCache = ".inputObjects", .neverCache = declared), ".inputObjects")$cache)
})
