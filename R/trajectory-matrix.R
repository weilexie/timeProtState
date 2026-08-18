## ============================================================================
## Trajectory matrix: shared-baseline-aware per-protein trajectories
## ============================================================================
##
## Implements STATISTICAL_CONTRACT.md section 7.1 and DECISIONS D6.
## The legacy Mfuzz analysis copied day0 baseline as t0=0 for each condition
## (3x inflation, artificially weighting baseline). This module fixes that:
##   - unique-cell matrix: baseline appears ONCE (the correct statistical input)
##   - display matrix: baseline optionally copied per condition FOR PLOTTING ONLY
##   - delta-from-baseline: log2FC vs baseline per condition x time
##
## All three views are derived from the SAME underlying group means, so they
## are mutually consistent and auditable.

## ---- group means ----------------------------------------------------------
.group_means <- function(se, condition_col, time_col, is_baseline_col,
                         assay_name, summary) {
  expr <- SummarizedExperiment::assay(se, assay_name)
  cd <- SummarizedExperiment::colData(se)
  cond <- as.character(cd[[condition_col]])
  tms <- cd[[time_col]]
  bl <- as.logical(cd[[is_baseline_col]])
  ## build cell labels: baseline rows -> "baseline"; others -> cond_time
  cell <- character(ncol(se))
  cell[bl] <- "baseline"
  cell[!bl] <- paste(cond[!bl], tms[!bl], sep = "_")
  cells <- unique(cell)
  ## summarize each protein within each cell.
  gm <- vapply(cells, function(cc) {
    cols <- which(cell == cc)
    if (length(cols) == 0L) rep(NA_real_, nrow(expr))
    else if (summary == "mean") {
      rowMeans(expr[, cols, drop = FALSE], na.rm = TRUE)
    } else {
      apply(expr[, cols, drop = FALSE], 1L, stats::median, na.rm = TRUE)
    }
  }, numeric(nrow(expr)))
  if (!is.matrix(gm)) gm <- matrix(gm, nrow = nrow(expr))
  rownames(gm) <- rownames(expr)
  colnames(gm) <- cells
  gm[is.nan(gm)] <- NA_real_
  list(means = gm, cells = cells)
}

## ---- tps_trajectory_matrix ------------------------------------------------
#' Build trajectory matrices for shared-baseline time-course data
#'
#' Constructs three mutually-consistent views of per-protein dynamics from a
#' shared-baseline design:
#' \itemize{
#'   \item \code{unique_cell}: group means with baseline appearing ONCE (the
#'     statistically correct input for clustering -- fixes the legacy 3x
#'     baseline-copy inflation).
#'   \item \code{display}: baseline copied per condition, for plotting only.
#'   \item \code{delta}: log2FC vs baseline per condition x time.
#' }
#'
#' @param se a \code{SummarizedExperiment} from \code{\link{tps_import}}.
#' @param features optional character vector of feature ids to subset; default
#'   all features.
#' @param condition,time,is_baseline column names in \code{colData}.
#' @param assay_name assay to use; default the one recorded at import.
#' @param summary function to summarize replicates within a cell; default
#'   \code{"mean"} (also \code{"median"}).
#' @return a list of class \code{tps_trajectory} with the three matrices plus
#'   cell ordering, missingness summary, and provenance.
#'
#' @examples
#' data(tps_example, package = "timeProtState")
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' traj <- tps_trajectory_matrix(se)
#' names(traj)
#' dim(traj$unique_cell)
#'
#' @export
tps_trajectory_matrix <- function(se, features = NULL,
                                   condition = "condition", time = "time",
                                   is_baseline = "is_baseline",
                                   assay_name = NULL,
                                   summary = c("mean","median")) {
  summary <- match.arg(summary)
  arg_capture <- list(call = "tps_trajectory_matrix",
                       condition = condition, time = time,
                       is_baseline = is_baseline, summary = summary,
                       timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  if (!inherits(se, "SummarizedExperiment")) {
    tps_abort("se must be a SummarizedExperiment.")
  }
  if (is.null(assay_name)) {
    assay_name <- S4Vectors::metadata(se)$timeProtState$assay_name
    if (is.null(assay_name)) assay_name <- "log2_abundance"
  }
  if (!all(c(condition, time, is_baseline) %in% names(SummarizedExperiment::colData(se)))) {
    tps_abort("colData must contain condition, time, is_baseline columns.")
  }

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  binfo <- .baseline_info(cd[[condition]], cd[[time]], cd[[is_baseline]])
  if (!binfo$shared_baseline) {
    tps_abort("tps_trajectory_matrix requires one shared baseline group.")
  }
  gm <- .group_means(se, condition, time, is_baseline, assay_name, summary)
  means <- gm$means
  if (!is.null(features)) {
    unknown <- setdiff(features, rownames(means))
    if (length(unknown)) {
      tps_abort("features not found in assay: %s.",
                paste(head(unknown, 5L), collapse = ", "))
    }
    means <- means[features, , drop = FALSE]
  }
  cells <- gm$cells
  bl <- as.logical(cd[[is_baseline]])
  cond_vec <- as.character(cd[[condition]])
  tms <- cd[[time]]
  conditions <- sort(unique(cond_vec[!bl]))
  times <- sort(unique(tms[!bl]))

  ## unique-cell matrix: baseline once + cond_time cells, ordered
  baseline_cell <- "baseline"
  cond_cells <- paste(rep(conditions, each = length(times)),
                       rep(times, length(conditions)), sep = "_")
  missing_cells <- setdiff(cond_cells, colnames(means))
  if (length(missing_cells)) {
    tps_abort("condition-by-time cells missing from trajectory input: %s.",
              paste(head(missing_cells, 8L), collapse = ", "))
  }
  uc_cols <- c(if (baseline_cell %in% cells) baseline_cell, cond_cells)
  uc_cols <- intersect(uc_cols, colnames(means))
  unique_cell <- means[, uc_cols, drop = FALSE]

  ## display matrix: baseline copied per condition (t0 prepended)
  display_cols <- unlist(lapply(conditions, function(cc) {
    c(baseline_cell, paste(cc, times, sep = "_"))
  }))
  display_cols <- display_cols[display_cols %in% colnames(means)]
  display_mat <- means[, display_cols, drop = FALSE]
  baseline_label <- format(binfo$baseline_time, trim = TRUE,
                           scientific = FALSE)
  display_names <- unlist(lapply(conditions, function(cc) {
    c(paste0(cc, "_baseline_t", baseline_label),
      paste(cc, times, sep = "_"))
  }))
  colnames(display_mat) <- display_names

  ## delta-from-baseline: log2FC vs baseline, per cond_time
  if (baseline_cell %in% colnames(means)) {
    bl_vec <- means[, baseline_cell]
    delta_mat <- sweep(means[, cond_cells, drop = FALSE], 1L, bl_vec, "-")
  } else {
    delta_mat <- means[, cond_cells, drop = FALSE]
    delta_mat[] <- NA_real_
  }

  ## missingness summary per feature
  if (nrow(means) == 0L) {
    miss_summary <- data.frame(
      feature_id = character(0),
      n_cells_observed = integer(0),
      n_cells_total = integer(0),
      pct_missing = numeric(0),
      stringsAsFactors = FALSE)
  } else {
    n_obs <- rowSums(!is.na(means))
    miss_summary <- data.frame(
      feature_id = rownames(means),
      n_cells_observed = n_obs,
      n_cells_total = ncol(means),
      pct_missing = round(100 * (1 - n_obs / ncol(means)), 2),
      stringsAsFactors = FALSE)
  }

  out <- structure(
    list(
      unique_cell = unique_cell,
      display = display_mat,
      delta = delta_mat,
      cells = uc_cols,
      conditions = conditions,
      times = times,
      summary = summary,
      missingness = miss_summary,
      provenance = list(.make_provenance(arg_capture, se, stage = "trajectory_matrix"))
    ),
    class = "tps_trajectory"
  )
  out
}

## ---- methods ---------------------------------------------------------------
#' @export
print.tps_trajectory <- function(x, ...) {
  cat("<timeProtState trajectory>\n")
  cat(sprintf("  features: %d   unique cells: %d\n",
              nrow(x$unique_cell), ncol(x$unique_cell)))
  cat(sprintf("  conditions: %s   times: %s\n",
              paste(x$conditions, collapse = ", "),
              paste(x$times, collapse = ",")))
  cat(sprintf("  display cells: %d (baseline copied per condition)\n",
              ncol(x$display)))
  cat(sprintf("  delta cells: %d (log2FC vs baseline)\n", ncol(x$delta)))
  invisible(x)
}
