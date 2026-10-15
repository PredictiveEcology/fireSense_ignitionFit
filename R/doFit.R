#' Should this run fit, use the ledger's fit, or stop?
#'
#' `doFit` is `NA` (fit only when the ledger has no usable fit for this polygon), `TRUE` (always
#' fit, replacing the ledger row) or `FALSE` (never fit: use the ledger fit, or stop if there is
#' none). The deprecated `refitExisting = TRUE` is `doFit = TRUE`; an explicit `doFit` is not
#' overridden by it.
#'
#' @param par the module's parameters, `P(sim)`.
#' @param message logical; send the deprecation message.
#' @return `NA`, `TRUE` or `FALSE`.
resolveDoFit <- function(par, message = FALSE) {
  if (!is.na(par$refitExisting) && message)
    message("fireSense_ignitionFit: parameter `refitExisting` is deprecated; ",
            "use `doFit` (NA: fit only if the ledger has no fit; TRUE: always fit; FALSE: never fit).")
  if (!is.na(par$doFit)) return(isTRUE(par$doFit))
  if (isTRUE(par$refitExisting)) return(TRUE)
  NA
}
