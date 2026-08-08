## ============================================================================
## Design: saturated group-means model for shared-baseline time-course
## ============================================================================
##
## Implements STATISTICAL_CONTRACT.md section 5.1. The design is a no-intercept
## cell-means model: X = model.matrix(~ 0 + group). In a shared-baseline design
## the baseline (day0) collapses to a single "baseline" cell shared across all
## conditions; otherwise each condition x time combination forms its own cell.
##
## The function also returns:
##   - the design matrix (full column rank expected)
##   - a group table (cell -> samples, n, estimable?)
##   - rank audit (rank, null space if any)
##   - provenance
##
## This layer does NOT fit anything; it only encodes the experimental design.

## ---- group composition ----------------------------------------------------
.compose_design_groups <- function(cd) {
  ## cd = colData as data.frame with condition, time, is_baseline.
  ## In shared-baseline design, all is_baseline rows -> one "baseline" group
  ## at the shared baseline time. Non-baseline rows -> "<condition>_<time>".
  bl <- cd$is_baseline
  group <- character(nrow(cd))
  if (any(bl)) {
    bl_time <- unique(cd$time[bl])
    group[bl] <- "baseline"
  } else {
    bl_time <- NA
  }
  nb <- !bl
  group[nb] <- paste(as.character(cd$condition[nb]),
                      as.character(cd$time[nb]), sep = "_")
  list(group = factor(group), baseline_time = bl_time,
       shared_baseline = length(bl_time) == 1L && !anyNA(bl_time))
}

## ---- tps_make_design -------------------------------------------------------
#' Construct a saturated group-means design
#'
#' Builds a no-intercept cell-means design matrix from the experimental design
#' encoded in a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' produced by \code{\link{tps_import}}. For a shared-baseline design (all
#' baseline rows at one timepoint), the baseline collapses to a single cell
#' shared by all conditions; otherwise a general group-means model is used.
#'
#' @param se a \code{SummarizedExperiment} produced by \code{tps_import}.
#' @param condition character; name of the condition column in \code{colData}.
#' @param time character; name of the time column.
#' @param is_baseline character; name of the baseline flag column.
#' @param batch optional character; name of a batch covariate column. If
#'   supplied, batch indicator columns are appended to the design.
#' @return a list of class \code{tps_design} with components:
#'   \item{design_matrix}{the numeric design matrix (samples x cells).}
#'   \item{group_table}{data.frame of cells, sample ids, and counts.}
#'   \item{estimability}{data.frame per-cell estimability (n_samples, estimable).}
#'   \item{rank_audit}{list with rank, ncol, full_rank, null_space dimension.}
#'   \item{baseline}{list with baseline_time, shared_baseline flag.}
#'   \item{provenance}{provenance entry.}
#'
#' @examples
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' design <- tps_make_design(se, condition = "condition", time = "time",
#'                            is_baseline = "is_baseline")
#' design$rank_audit$full_rank
#' head(design$group_table)
#'
#' @export
tps_make_design <- function(se, condition = "condition", time = "time",
                             is_baseline = "is_baseline", batch = NULL) {
  arg_capture <- list(call = "tps_make_design",
                       condition = condition, time = time,
                       is_baseline = is_baseline, batch = batch,
                       timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  if (!inherits(se, "SummarizedExperiment")) {
    tps_abort("se must be a SummarizedExperiment, got %s.", class(se)[1L])
  }
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  for (col in c(condition, time, is_baseline)) {
    if (!col %in% names(cd)) {
      tps_abort("colData has no column '%s'.", col)
    }
  }
  cd$.condition <- cd[[condition]]
  cd$.time <- cd[[time]]
  cd$.is_baseline <- as.logical(cd[[is_baseline]])

  ginfo <- .compose_design_groups(cd)
  group <- ginfo$group
  group_levels <- levels(group)
  ## order: baseline first, then condition x time sorted naturally
  is_baseline_level <- group_levels == "baseline"
  if (any(is_baseline_level)) {
    non_bl <- group_levels[!is_baseline_level]
  } else {
    non_bl <- group_levels
  }
  non_bl_ord <- non_bl[order(.natural_sort_key(non_bl))]
  group_levels <- c(if (any(is_baseline_level)) "baseline", non_bl_ord)
  group <- factor(as.character(group), levels = group_levels)

  ## design matrix ~0 + group
  design_matrix <- stats::model.matrix(~ 0 + group)
  colnames(design_matrix) <- sub("^group", "", colnames(design_matrix))

  ## optional batch covariate
  if (!is.null(batch)) {
    if (!batch %in% names(cd)) tps_abort("colData has no batch column '%s'.", batch)
    bcol <- cd[[batch]]
    if (!is.factor(bcol)) bcol <- factor(bcol)
    bmat <- stats::model.matrix(~ 0 + bcol)[, -1L, drop = FALSE]
    colnames(bmat) <- paste0("batch", levels(droplevels(bcol))[-1L])
    design_matrix <- cbind(design_matrix, bmat)
  }

  ## group table + estimability
  gt <- as.data.frame(table(group), responseName = "n_samples",
                       stringsAsFactors = FALSE)
  names(gt)[1L] <- "group"
  gt$estimable <- gt$n_samples >= 1L
  ## attach sample ids per group
  gt$samples <- vapply(gt$group, function(g) {
    paste(rownames(cd)[group == g], collapse = ",")
  }, character(1L))

  ## rank audit
  qr_d <- qr(design_matrix)
  rank <- qr_d$rank
  ncol_d <- ncol(design_matrix)
  rank_audit <- list(rank = rank, ncol = ncol_d,
                     full_rank = (rank == ncol_d),
                     null_space_dim = ncol_d - rank)

  out <- structure(
    list(
      design_matrix = design_matrix,
      group_table = gt,
      estimability = gt[, c("group","n_samples","estimable")],
      rank_audit = rank_audit,
      baseline = list(baseline_time = ginfo$baseline_time,
                       shared_baseline = ginfo$shared_baseline),
      condition = condition, time = time, is_baseline = is_baseline,
      group_levels = group_levels,
      provenance = list(.make_provenance(arg_capture, se, stage = "design"))
    ),
    class = "tps_design"
  )
  out
}

## natural sort key: split "<condition>_<time>" into (condition, numeric time)
.natural_sort_key <- function(x) {
  parts <- strsplit(x, "_", fixed = TRUE)
  cond <- vapply(parts, `[`, character(1L), 1L)
  tt <- vapply(parts, function(p) {
    if (length(p) >= 2L) suppressWarnings(as.numeric(p[2L])) else NA_real_
  }, numeric(1L))
  order(cond, tt)
}

## ---- methods for tps_design ------------------------------------------------
#' @export
print.tps_design <- function(x, ...) {
  cat("<timeProtState design>\n")
  ra <- x$rank_audit
  cat(sprintf("  cells: %d   rank: %d   full_rank: %s\n",
              ra$ncol, ra$rank, ra$full_rank))
  cat(sprintf("  shared_baseline: %s (baseline_time = %s)\n",
              x$baseline$shared_baseline,
              format(x$baseline$baseline_time)))
  bl <- x$group_table
  cat(sprintf("  groups: min=%d  median=%d  max=%d replicates\n",
              min(bl$n_samples), stats::median(bl$n_samples), max(bl$n_samples)))
  if (!ra$full_rank) {
    cat(sprintf("  WARNING: design is rank-deficient (null space dim = %d)\n",
                ra$null_space_dim))
  }
  invisible(x)
}

#' @export
summary.tps_design <- function(object, ...) {
  print(object)
  cat("\n  group table:\n")
  print(object$group_table[, c("group","n_samples","estimable")], row.names = FALSE)
  invisible(object)
}
