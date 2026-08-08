## ============================================================================
## Data-contract schema query functions
## ============================================================================
##
## Public helpers that let users (and tests) inspect the frozen contract.

#' Return the required metadata columns
#'
#' Lists the column names that must be present in the sample metadata passed
#' to \code{\link{tps_import}}. These are frozen by the data contract:
#' \code{sample_id}, \code{condition}, \code{time}, \code{is_baseline}.
#'
#' @return a character vector of required column names.
#'
#' @examples
#' tps_metadata_schema()
#'
#' @export
tps_metadata_schema <- function() {
  structure(
    list(
      required = tps_required_metadata_cols(),
      optional = tps_known_optional_cols()
    ),
    class = "tps_metadata_schema"
  )
}

#' @export
print.tps_metadata_schema <- function(x, ...) {
  cat("<timeProtState metadata schema>\n")
  cat("  required:", paste(x$required, collapse = ", "), "\n")
  cat("  optional:", paste(x$optional, collapse = ", "), "\n")
  invisible(x)
}
