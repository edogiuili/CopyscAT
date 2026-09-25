#' Package-internal state for a CopyscAT session
#'
#' CopyscAT keeps genome references, binning parameters and output paths in a
#' dedicated environment rather than passing them through every call. Populate
#' it with [initialiseEnvironment()] and [setOutputFile()] before running any
#' analysis step.
#'
#' The environment is exported because the original tutorial workflow reads
#' `scCNVCaller$locPrefix` and `scCNVCaller$outPrefix` directly when writing
#' user-side output files.
#'
#' @format An [environment].
#' @export
scCNVCaller <- new.env(parent = emptyenv())

#' Report whether the session environment has been initialised
#'
#' @param fields Character vector of names that must be present in
#'   [scCNVCaller].
#' @param call Function name to report in the error message.
#' @return Invisibly `TRUE`; throws an error listing the missing fields
#'   otherwise.
#' @noRd
assert_initialised <- function(fields, call = NULL) {
  missing <- fields[!vapply(fields, exists, logical(1), envir = scCNVCaller,
                            inherits = FALSE)]
  if (length(missing) > 0L) {
    hint <- if (any(c("locPrefix", "outPrefix") %in% missing)) {
      "Call setOutputFile() first."
    } else {
      "Call initialiseEnvironment() first."
    }
    stop(
      sprintf(
        "%sCopyscAT session is not initialised: missing %s. %s",
        if (is.null(call)) "" else paste0(call, "(): "),
        paste(missing, collapse = ", "),
        hint
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Is a numeric vector sorted in non-decreasing order?
#'
#' Replaces `S4Vectors::isSorted()`, which the package previously reached only
#' because a Bioconductor dependency happened to be attached.
#'
#' @param x Numeric vector.
#' @return `TRUE` when `x` is non-decreasing, ignoring `NA`s.
#' @noRd
is_sorted_ascending <- function(x) {
  x <- x[!is.na(x)]
  length(x) < 2L || all(diff(x) >= 0)
}
