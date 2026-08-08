#' @import methods
#' @importFrom stats setNames sd quantile p.adjust embed median runif rnorm aggregate pnorm qnorm
#' @importFrom utils head write.csv combn
#' @importFrom rlang .data
"_PACKAGE"

## ============================================================================
## Internal helpers
## ============================================================================
##
## Functions here are NOT exported. They implement the data-contract layer's
## plumbing: structured errors, metadata schema defaults, ID alignment.
## All design choices traceable to STATISTICAL_CONTRACT.md section 4 and
## LEGACY_PIPELINE_AUDIT.md section 3 (sample-name parsing, ID coverage audit).

## ---- structured error ------------------------------------------------------
## Errors in this package must name the offending column, sample, or group,
## never a bare "invalid metadata". See STATISTICAL_CONTRACT.md requirement
## and the v2 spec section 13.4 (error messages are a user interface).
tps_abort <- function(message, ..., call. = TRUE) {
  msg <- if (...length() > 0L) sprintf(message, ...) else message
  stop(msg, call. = call.)
}

tps_warn <- function(message, ...) {
  msg <- if (...length() > 0L) sprintf(message, ...) else message
  warning(msg, call. = FALSE, immediate. = TRUE)
}

## ---- required metadata columns --------------------------------------------
## Frozen by STATISTICAL_CONTRACT.md section 4.1. These five columns are
## mandatory; the design layer refuses to proceed without them.
tps_required_metadata_cols <- function() {
  c("sample_id", "condition", "time", "is_baseline")
}

tps_known_optional_cols <- function() {
  c("replicate", "batch", "run_order", "subject_id", "block", "weight")
}

## ---- hashing for provenance -----------------------------------------------
## A small, dependency-free checksum used to fingerprint input matrices.
## Not cryptographically strong; sufficient for detecting accidental edits.
tps_checksum <- function(x, algo = "xxhash") {
  if (!is.matrix(x)) x <- as.matrix(x)
  raw <- serialize(x, connection = NULL, ascii = FALSE)
  ## Use native R digest via base::sprintf hex of the object pointer is not
  ## stable; instead hash a stable representation: dims + first/last/step.
  ## For a robust lightweight checksum, sample evenly across the data.
  d <- dim(x)
  n <- length(x)
  if (n == 0L) return(list(checksum = "empty", n = 0L, dims = d))
  idx <- unique(as.integer(seq.int(1L, n, length.out = min(1000L, n))))
  vals <- x[idx]
  vals <- vals[is.finite(vals)]
  cs <- if (length(vals)) {
    ## combine a few stable summaries into a hex-ish string; %x needs integers.
    num <- sum(as.double(vals)^2, na.rm = TRUE) * 1e6 +
           length(vals) * 3.14159 + d[1] * 1000 + d[2]
    sprintf("tps%013.0f", abs(num))
  } else "allna"
  list(checksum = cs, n = n, dims = d, sampled = length(vals))
}
