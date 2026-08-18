## ============================================================================
## Import: build a SummarizedExperiment from a protein x sample matrix + metadata
## ============================================================================
##
## Design contract: STATISTICAL_CONTRACT.md section 4. The import function is
## the ONLY entry point that constructs the working object. It enforces:
##   - numeric abundance matrix
##   - unique protein ids
##   - sample_id <-> colnames alignment
##   - required metadata columns present
##   - time is finite numeric or coercible to numeric
##   - at least one baseline and >=2 non-baseline timepoints per design
##   - shared-baseline detection (all baselines share one timepoint)
##   - estimability per condition x time group (>=1 sample)
##   - no all-NA proteins (warning, with names)
## Importantly, this stage does NOT impute (DECISIONS_NEEDED D3: imputation is
## opt-in via a separate step, never the default).

## ---- default data-contract validator --------------------------------------
## A closure holding every check, so tests can drive each independently.

## ---- tps_import ------------------------------------------------------------
#' Construct a timeProtState experiment object
#'
#' Builds a \code{\link[SummarizedExperiment]{SummarizedExperiment}} from a
#' protein-by-sample abundance matrix and sample metadata, checking the
#' data contract required by the rest of the package: numeric time,
#' estimable condition-by-time groups, and aligned ids. A shared baseline
#' (one baseline time and one baseline condition label) is reported by
#' \code{\link{tps_validate}}; non-shared baselines warn at import and are
#' rejected later by \code{\link{tps_make_contrasts}} and
#' \code{\link{tps_trajectory_matrix}}.
#'
#' Imputation is NOT performed at import (missing values are preserved).
#' Missing values are preserved and surface naturally in downstream limma
#' fitting.
#'
#' @param abundance numeric matrix (or data.frame coercible to one), proteins
#'   x samples. Row names must be unique protein identifiers.
#' @param col_data data.frame of sample metadata. Must contain columns
#'   \code{sample_id}, \code{condition}, \code{time}, and \code{is_baseline}.
#'   \code{sample_id} must match \code{colnames(abundance)} one-to-one.
#' @param row_data optional data.frame of protein annotation, row-aligned to
#'   \code{abundance}. Useful for locus-tag / gene-name mapping.
#' @param assay_name character scalar, name of the assay stored in the object.
#' @param validate logical; if \code{TRUE} (default) run \code{\link{tps_validate}}
#'   on the result and abort on contract violations.
#' @return a \code{SummarizedExperiment} with provenance recorded in
#'   \code{metadata(se)$timeProtState}.
#'
#' @examples
#' data(tps_example, package = "timeProtState")
#' se <- tps_import(
#'   abundance = tps_example_abundance,
#'   col_data  = tps_example_metadata,
#'   row_data  = tps_example_annotation
#' )
#' se
#'
#' @export
tps_import <- function(abundance, col_data, row_data = NULL,
                        assay_name = "log2_abundance", validate = TRUE) {
  arg_capture <- list(call = "tps_import", assay_name = assay_name,
                       timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  se <- .build_se(abundance, col_data, row_data, assay_name)

  if (isTRUE(validate)) {
    report <- tps_validate(se, action = "abort")
    S4Vectors::metadata(se)$tps_validation <- report
  }
  S4Vectors::metadata(se)$timeProtState <- c(
    S4Vectors::metadata(se)$timeProtState,
    list(provenance = list(.make_provenance(arg_capture, se, stage = "import")))
  )
  se
}

## ---- internal builder ------------------------------------------------------
.build_se <- function(abundance, col_data, row_data, assay_name) {
  ## Coerce abundance to matrix.
  if (!is.matrix(abundance)) {
    if (is.data.frame(abundance)) {
      rns <- rownames(abundance)
      abundance <- as.matrix(abundance)
      rownames(abundance) <- rns
    } else {
      tps_abort("abundance must be a matrix or data.frame, got %s",
                 class(abundance)[1L])
    }
  }
  if (!is.numeric(abundance)) {
    original_na <- is.na(abundance)
    converted <- suppressWarnings(as.numeric(abundance))
    converted <- matrix(converted, nrow = nrow(abundance), ncol = ncol(abundance),
                        dimnames = dimnames(abundance))
    introduced_na <- is.na(converted) & !original_na
    if (any(introduced_na)) {
      tps_abort("abundance contains non-numeric values that cannot be converted safely.")
    }
    abundance <- converted
  } else {
    storage.mode(abundance) <- "double"
  }
  if (nrow(abundance) == 0L || ncol(abundance) == 0L) {
    tps_abort("abundance must have at least one feature and one sample.")
  }
  if (length(assay_name) != 1L || is.na(assay_name) || !nzchar(assay_name)) {
    tps_abort("assay_name must be one non-empty string.")
  }

  ## Protein ids.
  if (is.null(rownames(abundance))) {
    tps_abort("abundance must have unique row names (protein identifiers).")
  }
  if (anyNA(rownames(abundance)) || any(!nzchar(rownames(abundance)))) {
    tps_abort("abundance row names must be non-missing, non-empty protein ids.")
  }
  dup <- duplicated(rownames(abundance))
  if (any(dup)) {
    tps_abort("abundance has %d duplicated protein ids, e.g. %s.",
               sum(dup), head(rownames(abundance)[dup], 3L))
  }

  ## col_data must be a data.frame with required columns.
  if (!is.data.frame(col_data)) {
    tps_abort("col_data must be a data.frame, got %s.", class(col_data)[1L])
  }
  missing_cols <- setdiff(tps_required_metadata_cols(), names(col_data))
  if (length(missing_cols)) {
    tps_abort("col_data is missing required column(s): %s.",
               paste(missing_cols, collapse = ", "))
  }
  if (anyNA(col_data$sample_id) || any(!nzchar(as.character(col_data$sample_id)))) {
    tps_abort("col_data$sample_id must contain non-missing, non-empty ids.")
  }
  if (anyDuplicated(col_data$sample_id)) {
    tps_abort("col_data$sample_id must be unique.")
  }
  ## Ordinary data.frames have automatic row names (1, 2, ...), not NULL.
  ## Promote sample_id in that common input form so alignment works as documented.
  if (.has_automatic_rownames(col_data)) {
    rownames(col_data) <- col_data$sample_id
  }
  if (anyDuplicated(rownames(col_data))) {
    tps_abort("col_data row names must be unique sample ids.")
  }
  ## Enforce the documented invariant: rownames must equal sample_id values.
  ## Without this, manual rownames matching abundance colnames but with
  ## mismatched sample_id would silently produce a corrupt SE.
  if (!identical(rownames(col_data), as.character(col_data$sample_id))) {
    tps_abort("col_data row names must equal col_data$sample_id; got mismatched values.")
  }

  ## Align samples: colnames(abundance) <-> col_data row order.
  csamp <- colnames(abundance)
  if (is.null(csamp)) {
    tps_abort("abundance must have column names matching sample_id.")
  }
  if (anyNA(csamp) || any(!nzchar(csamp)) || anyDuplicated(csamp)) {
    tps_abort("abundance column names must be unique, non-missing sample ids.")
  }
  unknown_samples <- setdiff(csamp, rownames(col_data))
  if (length(unknown_samples)) {
    tps_abort("abundance column(s) have no metadata row: %s.",
               head(unknown_samples, 5L))
  }
  extra_meta <- setdiff(rownames(col_data), csamp)
  if (length(extra_meta)) {
    tps_abort("col_data has row(s) with no abundance column: %s.",
               head(extra_meta, 5L))
  }
  col_data <- col_data[csamp, , drop = FALSE]

  ## Coerce the baseline flag before any downstream grouping logic uses it.
  col_data$is_baseline <- .as_logical_strict(col_data$is_baseline,
                                              col_name = "is_baseline")
  ## Factor hygiene: condition -> factor (baseline kept as its own level).
  col_data$condition <- .as_condition_factor(col_data$condition,
                                              col_data$is_baseline)
  ## Downstream contrast names and dynamic ordering require numeric time.
  col_data$time <- .as_time(col_data$time)

  ## row_data alignment.
  if (is.null(row_data)) {
    row_data <- data.frame(rowname = rownames(abundance),
                            stringsAsFactors = FALSE)
    rownames(row_data) <- rownames(abundance)
  } else {
    if (!is.data.frame(row_data)) {
      tps_abort("row_data must be a data.frame or NULL.")
    }
    if (.has_automatic_rownames(row_data)) {
      if (nrow(row_data) == nrow(abundance)) {
        rownames(row_data) <- rownames(abundance)
      } else {
        tps_abort("row_data has no row names and length mismatches abundance.")
      }
    }
    if (anyDuplicated(rownames(row_data))) {
      tps_abort("row_data row names must be unique feature ids.")
    }
    if (anyNA(rownames(row_data)) || any(!nzchar(rownames(row_data)))) {
      tps_abort("row_data row names must be non-missing, non-empty feature ids.")
    }
    missing_prots <- setdiff(rownames(abundance), rownames(row_data))
    if (length(missing_prots)) {
      tps_abort("abundance row(s) have no row_data annotation: %s.",
                 head(missing_prots, 5L))
    }
    row_data <- row_data[rownames(abundance), , drop = FALSE]
  }

  ## Assemble. rowData uses a DataFrame for SummarizedExperiment interop.
  assays <- structure(list(abundance), names = assay_name)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    rowData = S4Vectors::DataFrame(row_data, row.names = rownames(abundance)),
    colData = S4Vectors::DataFrame(col_data, row.names = csamp),
    metadata = list(timeProtState = list(assay_name = assay_name,
                                          package = "timeProtState"))
  )
  se
}

.has_automatic_rownames <- function(x) {
  .row_names_info(x, type = 1L) < 0L
}

## ---- coercion helpers ------------------------------------------------------
.as_condition_factor <- function(x, is_baseline) {
  if (is.factor(x)) return(droplevels(x))
  if (!is.character(x)) x <- as.character(x)
  factor(x)
}

.as_time <- function(x) {
  out <- if (is.numeric(x)) x else {
    suppressWarnings(as.numeric(as.character(x)))
  }
  if (length(out) != length(x) || any(!is.finite(out))) {
    tps_abort("time must be fully observed, finite numeric values or coercible to numeric.")
  }
  out
}

.as_logical_strict <- function(x, col_name) {
  if (is.logical(x)) {
    if (anyNA(x)) tps_abort("column %s must not contain NA.", col_name)
    return(x)
  }
  if (is.character(x)) {
    ok <- tolower(x) %in% c("true","false","t","f","yes","no","1","0")
    if (!all(ok)) {
      tps_abort("column %s has non-logical values: %s.",
                 col_name, paste(unique(x[!ok]), collapse = ", "))
    }
    lx <- tolower(x) %in% c("true","t","yes","1")
    return(lx)
  }
  if (is.numeric(x)) {
    if (anyNA(x) || any(!x %in% c(0, 1))) {
      tps_abort("column %s numeric values must be only 0 or 1.", col_name)
    }
    return(as.logical(x))
  }
  tps_abort("column %s must be logical.", col_name)
}
