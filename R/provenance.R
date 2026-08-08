## ============================================================================
## Provenance: track function, params, versions, checksum for every stage
## ============================================================================
##
## Every public function appends to metadata(se)$timeProtState$provenance.
## This is mandatory (STATISTICAL_CONTRACT 4.3) — never print-and-forget.

.make_provenance <- function(arg_capture, se, stage) {
  ## The input object may be a SummarizedExperiment (matrix fingerprint),
  ## a tps_* result list (summary fingerprint), or anything else.
  if (inherits(se, "SummarizedExperiment")) {
    assay_name <- S4Vectors::metadata(se)$timeProtState$assay_name %||% "log2_abundance"
    expr <- SummarizedExperiment::assay(se, assay_name)
    cs <- tps_checksum(expr)
    input_checksum <- cs$checksum
    input_dims <- cs$dims
    input_n <- cs$n
  } else if (is.list(se)) {
    ## fingerprint a result object by its class + a digest of top-level names
    ## and key scalar fields; cheap and stable.
    input_checksum <- sprintf("tps%s", substr(digest_list(se), 1, 24))
    input_dims <- length(se)
    input_n <- NA_integer_
  } else {
    input_checksum <- sprintf("tps%s", substr(digest_list(list(obj = se)), 1, 24))
    input_dims <- length(se)
    input_n <- NA_integer_
  }
  pv <- as.character(getNamespaceVersion("timeProtState"))
  list(
    stage = stage,
    function_call = arg_capture$call,
    timestamp = arg_capture$timestamp,
    assay_name = arg_capture$assay_name,
    input_checksum = input_checksum,
    input_dims = input_dims,
    input_n = input_n,
    package_version = pv,
    r_version = R.version.string,
    params = arg_capture[setdiff(names(arg_capture),
                                  c("call","timestamp","assay_name"))]
  )
}

## lightweight list digest without a digest dependency: hash a stable
## representation of names + numeric summaries of any matrices/vectors.
digest_list <- function(x) {
  nm <- paste(names(x), collapse = ",")
  parts <- nm
  for (nn in names(x)) {
    v <- x[[nn]]
    if (is.matrix(v) || is.numeric(v)) {
      vv <- v[is.finite(v)]
      if (length(vv)) parts <- paste0(parts, "|", round(sum(vv^2), 4),
                                       ":", length(vv))
    } else if (is.character(v) && length(v) < 50L) {
      parts <- paste0(parts, "|c:", paste(head(v, 5L), collapse = ";"))
    } else if (is.logical(v)) {
      parts <- paste0(parts, "|l:", sum(v, na.rm = TRUE))
    }
  }
  sprintf("%016.0f", sum(utf8ToInt(parts) ^ 2) %% 1e16)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Retrieve the provenance record from an object
#'
#' Extracts the recorded provenance list from a
#' \code{SummarizedExperiment} (or a result object carrying one). Each
#' entry records the stage, function call, timestamp, input checksum, and
#' package version.
#'
#' @param object a \code{SummarizedExperiment} produced by \code{tps_import}
#'   or a downstream result object.
#' @return a list of provenance entries (one per stage applied), or
#'   \code{NULL} if none recorded.
#'
#' @examples
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' pv <- tps_provenance(se)
#' if (!is.null(pv)) pv[[1]]$stage
#'
#' @export
tps_provenance <- function(object) {
  if (inherits(object, "SummarizedExperiment")) {
    md <- S4Vectors::metadata(object)$timeProtState
  } else if (is.list(object) && "provenance" %in% names(object)) {
    md <- object
  } else if (is.list(object) && "timeProtState" %in% names(object)) {
    md <- object$timeProtState
  } else {
    return(NULL)
  }
  md$provenance
}
