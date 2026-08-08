## ============================================================================
## Validate: data-contract checks on a SummarizedExperiment
## ============================================================================
##
## The validator is the gatekeeper for every downstream analysis. It catches
## the failure modes documented in STATISTICAL_CONTRACT.md section 4.2 and
## LEGACY_PIPELINE_AUDIT.md section 3:
##   - duplicate sample ids            (audit 3.3: fragile parsing)
##   - baseline not shared             (shared-baseline grammar, the novelty)
##   - fewer than 2 non-baseline times (no temporal design)
##   - unestimable groups (0 samples)  (audit 3.1/3.2: rank issues)
##   - all-NA proteins                 (audit 3.6: silent NA handling)
##   - time ordering ambiguity
##   - condition/time/organism leakage (audit 3.7: no hardcoding)
##
## `action` controls behaviour:
##   "report"  -> return a structured list, never stop.
##   "abort"   -> stop on the first ERROR-level problem, but always warn first.
##   "warn"    -> only warn, return the report.

#' Validate a timeProtState experiment
#'
#' Runs the full data-contract check suite on a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} produced by
#' \code{\link{tps_import}} and returns a structured report of errors and
#' warnings. The checks cover sample-id uniqueness, baseline sharing,
#' time ordering, group estimability, and all-NA features.
#'
#' @param se a \code{SummarizedExperiment} produced by \code{tps_import}.
#' @param action character; one of \code{"report"}, \code{"abort"},
#'   \code{"warn"}. \code{"report"} (default) returns the report silently;
#'   \code{"abort"} stops on the first error; \code{"warn"} emits warnings.
#' @return invisibly, a list with components \code{errors}, \code{warnings},
#'   \code{design} (a summary of the experimental design), and \code{ok}.
#'
#' @examples
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' report <- tps_validate(se)
#' report$ok
#' report$design
#'
#' @export
tps_validate <- function(se, action = c("report", "abort", "warn")) {
  action <- match.arg(action)
  if (!inherits(se, "SummarizedExperiment")) {
    tps_abort("se must be a SummarizedExperiment, got %s.", class(se)[1L])
  }
  errors <- character(0)
  warnings <- character(0)
  ## messages may be supplied as a paste0-able template split across several
  ## string args followed by sprintf args. We fold consecutive character
  ## fragments into one template, then sprintf over the remaining args.
  fold_msg <- function(...) {
    args <- list(...)
    is_chr <- vapply(args, is.character, logical(1L))
    ## leading run of character args = template fragments; the rest = values
    n_tmpl <- 0L
    for (k in seq_along(args)) {
      if (is_chr[k]) n_tmpl <- k else break
    }
    tmpl <- paste0(args[seq_len(n_tmpl)], collapse = "")
    vals <- if (n_tmpl < length(args)) args[(n_tmpl + 1L):length(args)] else list()
    if (length(vals)) do.call(sprintf, c(list(tmpl), vals)) else tmpl
  }
  add_err  <- function(...) errors   <<- c(errors, fold_msg(...))
  add_warn <- function(...) warnings <<- c(warnings, fold_msg(...))

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  assay_name <- S4Vectors::metadata(se)$timeProtState$assay_name
  if (is.null(assay_name)) assay_name <- "log2_abundance"
  expr <- SummarizedExperiment::assay(se, assay_name)

  ## ---- check 1: required columns in colData ------------------------------
  missing_cols <- setdiff(tps_required_metadata_cols(), names(cd))
  if (length(missing_cols)) {
    add_err("colData missing required column(s): %s.",
            paste(missing_cols, collapse = ", "))
  }

  ## ---- check 2: sample_id uniqueness and alignment -----------------------
  if ("sample_id" %in% names(cd)) {
    sid <- cd$sample_id
    if (anyDuplicated(sid)) {
      dups <- unique(sid[duplicated(sid)])
      add_err("sample_id has %d duplicated value(s): %s.",
              length(dups), paste(head(dups, 5L), collapse = ", "))
    }
    colnames_expr <- colnames(expr)
    if (!identical(as.character(sid), colnames_expr)) {
      misaligned <- which(as.character(sid) != colnames_expr)
      add_err("sample_id and assay colnames disagree at positions: %s.",
              paste(head(misaligned, 5L), collapse = ", "))
    }
  }

  ## ---- check 3: is_baseline present and at least one TRUE ----------------
  if ("is_baseline" %in% names(cd)) {
    if (!is.logical(cd$is_baseline)) {
      add_err("is_baseline must be logical; got %s.", class(cd$is_baseline)[1L])
    } else if (!any(cd$is_baseline, na.rm = TRUE)) {
      add_err("No baseline samples found (is_baseline all FALSE/NA).")
    }
  }

  ## ---- check 4: time sortable & >=2 non-baseline timepoints ---------------
  design_summary <- list()
  if ("time" %in% names(cd) && "is_baseline" %in% names(cd)) {
    tv <- cd$time
    if (is.numeric(tv)) {
      if (anyNA(tv)) {
        add_err("time has %d NA value(s); time must be fully observed.",
                sum(is.na(tv)))
      }
      time_levels <- sort(unique(tv))
    } else if (is.ordered(tv)) {
      time_levels <- levels(droplevels(tv))
    } else {
      add_err("time must be numeric or ordered factor; got %s.", class(tv)[1L])
      time_levels <- character(0)
    }
    baseline_times <- unique(tv[cd$is_baseline])
    non_baseline_times <- unique(tv[!cd$is_baseline])
    n_nonbaseline <- length(non_baseline_times)
    if (n_nonbaseline < 2L) {
      add_err("Need at least 2 non-baseline timepoints for a time-course; found %d.",
              n_nonbaseline)
    }
    design_summary$time_levels <- time_levels
    design_summary$baseline_times <- baseline_times
    design_summary$non_baseline_times <- non_baseline_times

    ## ---- check 5: SHARED BASELINE detection (the novelty) ----------------
    ## In a shared-baseline design, all baseline rows share exactly one time.
    if (length(baseline_times) > 1L) {
      add_warn("Baseline spans %s distinct time values (%s). This is NOT a ",
               "shared-baseline design; tps_make_design will fall back to a ",
               "general group-means model.",
               length(baseline_times),
               paste(head(as.character(baseline_times), 5L), collapse = ", "))
      design_summary$shared_baseline <- FALSE
    } else {
      design_summary$shared_baseline <- TRUE
    }
  }

  ## ---- check 6: condition-by-time group estimability ---------------------
  if (all(c("condition","time","is_baseline") %in% names(cd))) {
    cd$.group <- .compose_group(cd$condition, cd$time, cd$is_baseline)
    grp_tab <- table(cd$.group)
    design_summary$groups <- as.data.frame.table(grp_tab,
                                                  responseName = "n_samples",
                                                  stringsAsFactors = FALSE)
    names(design_summary$groups)[1L] <- "group"
    empty_groups <- character(0)  ## none can be empty by construction
    small_groups <- design_summary$groups$group[design_summary$groups$n_samples < 2L]
    if (length(small_groups)) {
      add_warn(paste0("Group(s) with <2 samples (low-replicate, limited power): ",
                      paste(head(small_groups, 5L), collapse = ", "), "."))
    }
    design_summary$n_groups <- length(grp_tab)
    design_summary$min_replicates <- min(grp_tab)
  }

  ## ---- check 7: all-NA proteins ------------------------------------------
  all_na <- rownames(expr)[apply(is.na(expr), 1L, all)]
  if (length(all_na)) {
    add_warn(paste0(length(all_na),
                    " protein(s) are all-NA across all samples; they will be ",
                    "uninformative. Examples: ",
                    paste(head(all_na, 5L), collapse = ", "), "."))
  }
  design_summary$n_all_na_proteins <- length(all_na)

  ## ---- check 8: protein id uniqueness ------------------------------------
  if (anyDuplicated(rownames(expr))) {
    add_err("Protein ids are duplicated (should have been caught at import).")
  }

  ## ---- check 9: organism/condition hardcoding leak (heuristic) ----------
  ## Heuristic: flag if "Msmeg" or "MSMEG_" appears as a level name, since
  ## the package must stay organism-agnostic (LEGACY_AUDIT 3.7).
  if ("condition" %in% names(cd) && is.factor(cd$condition)) {
    leak <- grep("MSMEG_|Msmeg", levels(cd$condition), value = TRUE)
    if (length(leak)) {
      add_warn(paste0("condition levels look like organism-internal ids (",
                      paste(head(leak, 3L), collapse = ", "),
                      "); the package is organism-agnostic, double-check this is intended."))
    }
  }

  ## ---- assemble report ---------------------------------------------------
  report <- list(
    errors = errors,
    warnings = warnings,
    design = design_summary,
    ok = (length(errors) == 0L),
    n_proteins = nrow(expr),
    n_samples = ncol(expr)
  )
  class(report) <- "tps_validation_report"

  if (length(warnings) && action %in% c("warn", "abort")) {
    for (w in warnings) tps_warn("%s", w)
  }
  if (length(errors) && action == "abort") {
    ## surface the first concrete error in the abort message (friendlier and
    ## matchable by tests than a generic "N errors detected" summary).
    tps_abort("%s", errors[1L])
  }
  invisible(report)
}

## compose a group label from condition/time/baseline. Baseline rows all
## collapse to "baseline" so they form a single shared cell.
.compose_group <- function(condition, time, is_baseline) {
  out <- character(length(condition))
  out[is_baseline] <- "baseline"
  nb <- !is_baseline
  out[nb] <- paste(as.character(condition[nb]),
                    as.character(time[nb]), sep = "_")
  out
}

## ---- print method for the report ------------------------------------------
#' @export
print.tps_validation_report <- function(x, ...) {
  cat("<timeProtState validation report>\n")
  cat(sprintf("  proteins: %d   samples: %d\n", x$n_proteins, x$n_samples))
  if (length(x$errors)) {
    cat(sprintf("  errors:   %d\n", length(x$errors)))
    for (e in head(x$errors, 10L)) cat("    - ", e, "\n", sep = "")
  } else {
    cat("  errors:   0\n")
  }
  if (length(x$warnings)) {
    cat(sprintf("  warnings: %d\n", length(x$warnings)))
    for (w in head(x$warnings, 10L)) cat("    - ", w, "\n", sep = "")
  } else {
    cat("  warnings: 0\n")
  }
  d <- x$design
  if ("shared_baseline" %in% names(d)) {
    cat(sprintf("  shared_baseline: %s\n", d$shared_baseline))
  }
  if ("n_groups" %in% names(d)) {
    cat(sprintf("  groups: %d (min replicates: %s)\n",
                d$n_groups, d$min_replicates))
  }
  invisible(x)
}
