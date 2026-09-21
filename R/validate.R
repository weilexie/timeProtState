## ============================================================================
## Validate: data-contract checks on a SummarizedExperiment
## ============================================================================
##
## The validator is the gatekeeper for every downstream analysis. It catches
## the failure modes documented in STATISTICAL_CONTRACT.md and
## ID / sample-name checks:
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
#' data(tps_example, package = "timeProtState")
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
  issues <- new.env(parent = emptyenv())
  issues$errors <- character(0)
  issues$warnings <- character(0)
  ## The first argument is always the sprintf template; remaining args
  ## (character or numeric) are sprintf values. This avoids the old
  ## is.character() heuristic that mangled %s when all args were character.
  fold_msg <- function(...) {
    args <- list(...)
    if (length(args) == 0L) return(character(0L))
    tmpl <- args[[1L]]
    vals <- if (length(args) > 1L) args[-1L] else list()
    if (length(vals)) do.call(sprintf, c(list(tmpl), vals)) else tmpl
  }
  add_err <- function(...) {
    issues$errors <- c(issues$errors, fold_msg(...))
  }
  add_warn <- function(...) {
    issues$warnings <- c(issues$warnings, fold_msg(...))
  }

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  assay_name <- S4Vectors::metadata(se)$timeProtState$assay_name
  if (is.null(assay_name)) assay_name <- "log2_abundance"
  if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
    add_err("assay_name '%s' not found in assays (%s).",
            assay_name,
            paste(SummarizedExperiment::assayNames(se), collapse = ", "))
    expr <- SummarizedExperiment::assay(se, 1)
  } else {
    expr <- SummarizedExperiment::assay(se, assay_name)
  }

  ## ---- check 1: required columns in colData ------------------------------
  missing_cols <- setdiff(tps_required_metadata_cols(), names(cd))
  if (length(missing_cols)) {
    add_err("colData missing required column(s): %s.",
            paste(missing_cols, collapse = ", "))
  }

  ## ---- check 2: sample_id uniqueness and alignment -----------------------
  if ("sample_id" %in% names(cd)) {
    sid <- cd$sample_id
    if (anyNA(sid) || any(!nzchar(as.character(sid)))) {
      add_err("sample_id must contain non-missing, non-empty values.")
    }
    if (anyDuplicated(sid)) {
      dups <- unique(sid[duplicated(sid)])
      add_err("sample_id has %d duplicated value(s): %s.",
              length(dups), paste(head(dups, 5L), collapse = ", "))
    }
    colnames_expr <- colnames(expr)
    if (!identical(as.character(sid), colnames_expr)) {
      misaligned <- which(is.na(sid) | as.character(sid) != colnames_expr)
      add_err("sample_id and assay colnames disagree at positions: %s.",
              paste(head(misaligned, 5L), collapse = ", "))
    }
  }

  ## ---- check 3: is_baseline present and at least one TRUE ----------------
  baseline_flag_ok <- FALSE
  if ("is_baseline" %in% names(cd)) {
    if (!is.logical(cd$is_baseline)) {
      add_err("is_baseline must be logical; got %s.", class(cd$is_baseline)[1L])
    } else if (anyNA(cd$is_baseline)) {
      add_err("is_baseline must not contain NA.")
    } else if (!any(cd$is_baseline, na.rm = TRUE)) {
      add_err("No baseline samples found (is_baseline all FALSE/NA).")
    } else {
      baseline_flag_ok <- TRUE
    }
  }

  if ("condition" %in% names(cd) &&
      (anyNA(cd$condition) || any(!nzchar(as.character(cd$condition))))) {
    add_err("condition must contain non-missing, non-empty labels.")
  }

  ## ---- check 4: finite numeric time & >=2 non-baseline timepoints ----------
  design_summary <- list()
  if ("time" %in% names(cd) && baseline_flag_ok) {
    tv <- cd$time
    if (is.numeric(tv) && all(is.finite(tv))) {
      if (any(tv < 0)) add_warn("time contains %d negative value(s); the package assumes forward time progression.",
                                sum(tv < 0))
      time_levels <- sort(unique(tv))
    } else {
      add_err("time must contain fully observed, finite numeric values.")
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
    ## A shared baseline is one physical baseline group: one time and one
    ## condition label. Equal time alone is insufficient because condition-
    ## specific controls measured at t0 must not be pooled.
    if ("condition" %in% names(cd)) {
      binfo <- .baseline_info(cd$condition, tv, cd$is_baseline)
      design_summary$shared_baseline <- binfo$shared_baseline
      design_summary$baseline_conditions <- binfo$baseline_conditions
      if (!binfo$shared_baseline) {
        add_warn(paste0(
          "Baseline rows are NOT a valid shared-baseline group (times: ",
          paste(as.character(binfo$baseline_time), collapse = ", "),
          "; conditions: ",
          paste(binfo$baseline_conditions, collapse = ", "),
          "). Shared-baseline contrasts are unavailable."
        ))
      }
    }

    ## All non-baseline conditions must share the same time grid because the
    ## contrast families and pathway dynamics compare condition x time cells.
    if ("condition" %in% names(cd)) {
      nb <- !cd$is_baseline
      conditions <- sort(unique(as.character(cd$condition[nb])))
      expected <- expand.grid(
        condition = conditions,
        time = sort(unique(tv[nb])),
        stringsAsFactors = FALSE
      )
      expected_key <- paste(expected$condition, expected$time, sep = "\r")
      observed_key <- unique(paste(as.character(cd$condition[nb]),
                                   as.character(tv[nb]), sep = "\r"))
      missing_key <- expected_key[!expected_key %in% observed_key]
      if (length(missing_key)) {
        add_err("Incomplete condition-by-time grid; missing cell(s): %s.",
                paste(head(gsub("\r", "@", missing_key, fixed = TRUE), 8L),
                      collapse = ", "))
      }
    }
  }

  ## ---- check 6: condition-by-time group estimability ---------------------
  if (all(c("condition", "time") %in% names(cd)) && baseline_flag_ok &&
      is.numeric(cd$time) && all(is.finite(cd$time)) &&
      !anyNA(cd$condition)) {
    cd$.group <- .compose_group(cd$condition, cd$time, cd$is_baseline)
    grp_tab <- table(cd$.group)
    design_summary$groups <- as.data.frame.table(grp_tab,
                                                  responseName = "n_samples",
                                                  stringsAsFactors = FALSE)
    names(design_summary$groups)[1L] <- "group"
    small_groups <- design_summary$groups$group[design_summary$groups$n_samples < 2L]
    if (length(small_groups)) {
      add_warn(paste0("Group(s) with <2 samples (low-replicate, limited power): ",
                      paste(head(small_groups, 5L), collapse = ", "), "."))
    }
    design_summary$n_groups <- length(grp_tab)
    design_summary$min_replicates <- min(grp_tab)
  }

  ## ---- check 7: all-NA proteins and samples -------------------------------
  all_na <- rownames(expr)[apply(is.na(expr), 1L, all)]
  if (length(all_na)) {
    add_warn(paste0(length(all_na),
                    " protein(s) are all-NA across all samples; they will be ",
                    "uninformative. Examples: ",
                    paste(head(all_na, 5L), collapse = ", "), "."))
  }
  design_summary$n_all_na_proteins <- length(all_na)
  all_na_samples <- colnames(expr)[apply(is.na(expr), 2L, all)]
  if (length(all_na_samples)) {
    add_err("Sample(s) are all-NA across all proteins: %s. Remove or fix them before fitting.",
            paste(head(all_na_samples, 5L), collapse = ", "))
  }
  if (any(is.infinite(expr))) {
    add_err("assay contains %d infinite value(s); replace them with NA or finite values.",
            sum(is.infinite(expr)))
  }

  ## ---- check 8: protein id uniqueness ------------------------------------
  if (anyDuplicated(rownames(expr))) {
    add_err("Protein ids are duplicated (should have been caught at import).")
  }

  ## ---- check 9: id-like condition levels (generic heuristic) --------------
  ## Flag locus-tag-like condition levels (e.g. "MSMEG_1234", "YAL001W"),
  ## since condition names should describe experimental conditions, not
  ## feature ids. Pattern-based and organism-agnostic: matches short
  ## alphabetic prefixes with a numeric suffix typical of systematic ids.
  if ("condition" %in% names(cd) && is.factor(cd$condition)) {
    leak <- grep("^[A-Za-z]{2,8}_[0-9]{3,}$", levels(cd$condition),
                 value = TRUE)
    if (length(leak)) {
      add_warn(paste0("condition levels look like systematic feature ids (",
                      paste(head(leak, 3L), collapse = ", "),
                      "); conditions should name experimental groups, ",
                      "double-check this is intended."))
    }
  }

  ## ---- assemble report ---------------------------------------------------
  errors <- issues$errors
  warnings <- issues$warnings
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
  binfo <- .baseline_info(condition, time, is_baseline)
  out <- paste(as.character(condition), as.character(time), sep = "_")
  if (binfo$shared_baseline) out[is_baseline] <- "baseline"
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
