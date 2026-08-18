## ============================================================================
## Contrasts: estimable contrast families for shared-baseline time-course
## ============================================================================
##
## Implements STATISTICAL_CONTRACT.md section 5.2. Six contrast families,
## each a sub-matrix of the full contrast matrix, with per-family rank audit.
##
##   A. vs_baseline      : mu_{cond,time} - mu_baseline          (single t)
##   B. within_condition : consecutive or all-pairs within a condition (single t)
##   C. average_time     : mean_c(mu_{c,t}) - mu_baseline        (omnibus F)
##   D. heterogeneity    : mu_{c,t} - mu_{ref_c,t}, ref per time (omnibus F, full rank)
##   E. within_omnibus   : all vs_baseline of a condition jointly (omnibus F)
##   F. endpoint         : mu_{cond,endpoint} - mu_baseline      (single t)
##
## Each family records its rank so downstream F-tests are auditable.

## ---- helpers to parse condition x time cells from design levels ----------
## A design level is "baseline" or "<condition>_<time>". Returns a data.frame
## with columns cell, condition, time, is_baseline.
.parse_cells <- function(group_levels) {
  is_bl <- group_levels == "baseline"
  cond <- time <- rep(NA_character_, length(group_levels))
  if (any(is_bl)) {
    cond[is_bl] <- "baseline"
    time[is_bl] <- "baseline"
  }
  nb <- !is_bl
  base <- group_levels[nb]
  cond[nb] <- sub("_[^_]+$", "", base)
  time[nb] <- sub("^.*_", "", base)
  data.frame(cell = group_levels, condition = cond, time = time,
             is_baseline = is_bl, stringsAsFactors = FALSE)
}

## ---- tps_make_contrasts ----------------------------------------------------
#' Construct contrast families for a shared-baseline design
#'
#' Generates estimable, rank-audited contrast families from a
#' \code{tps_design} object (produced by \code{\link{tps_make_design}}):
#' vs-baseline, within-condition temporal, average-time, condition
#' heterogeneity, within-condition omnibus, and endpoint contrasts. Each family
#' is returned with its effective rank.
#'
#' @param design a \code{tps_design} object from \code{\link{tps_make_design}}.
#' @param endpoint the timepoint value to treat as the endpoint (e.g. the last
#'   time). Used for the endpoint family. If \code{NULL}, the largest numeric
#'   non-baseline time is used.
#' @param within character; temporal contrast style within a condition: one of
#'   \code{"consecutive"} (default) or \code{"all_pairs"}.
#' @param reference_condition character; the condition used as reference for
#'   the heterogeneity family. If \code{NULL}, the first non-baseline
#'   condition (alphabetical) is used. The heterogeneity contrast subspace is
#'   reference-invariant. With feature-level missing values, numerical limma
#'   F-statistics can differ slightly because \code{contrasts.fit} uses an
#'   approximate covariance transformation.
#' @param include_average_time logical; include the average-time family (C).
#' @param include_heterogeneity logical; include the heterogeneity family (D).
#' @param include_within_omnibus logical; include within-condition omnibus (E).
#' @return a list of class \code{tps_contrasts} with \code{families} (named
#'   list of contrast matrices), \code{rank_audit} (per-family rank), a
#'   flattened \code{contrast_table}, and \code{provenance}.
#'
#' @examples
#' data(tps_example, package = "timeProtState")
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' design <- tps_make_design(se)
#' contrasts <- tps_make_contrasts(design, endpoint = 7)
#' names(contrasts$families)
#' contrasts$rank_audit
#'
#' @export
tps_make_contrasts <- function(design, endpoint = NULL,
                                within = c("consecutive", "all_pairs"),
                                reference_condition = NULL,
                                include_average_time = TRUE,
                                include_heterogeneity = TRUE,
                                include_within_omnibus = TRUE) {
  within <- match.arg(within)
  arg_capture <- list(call = "tps_make_contrasts", endpoint = endpoint,
                       within = within, reference_condition = reference_condition,
                       include_average_time = include_average_time,
                       include_heterogeneity = include_heterogeneity,
                       include_within_omnibus = include_within_omnibus,
                       timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  if (!inherits(design, "tps_design")) {
    tps_abort("design must be a tps_design object from tps_make_design.")
  }
  flags <- c(include_average_time, include_heterogeneity,
             include_within_omnibus)
  if (!is.logical(flags) || length(flags) != 3L || anyNA(flags)) {
    tps_abort("contrast-family include flags must be TRUE or FALSE.")
  }
  levels_d <- colnames(design$design_matrix)
  cells <- .parse_cells(design$group_levels)
  if (!design$baseline$shared_baseline) {
    tps_abort(paste0(
      "tps_make_contrasts requires a shared-baseline design; this ",
      "design is not shared-baseline. Use a general group-means ",
      "approach or fix the input."
    ))
  }
  if (!is.null(design$rank_audit) && !isTRUE(design$rank_audit$full_rank)) {
    tps_abort(paste0(
      "design matrix is rank-deficient (null space dimension = ",
      design$rank_audit$null_space_dim,
      "); resolve confounded cells or batch covariates before fitting."
    ))
  }
  conditions <- sort(unique(cells$condition[!cells$is_baseline]))
  if (length(conditions) < 1L) {
    tps_abort("No non-baseline conditions found; cannot build contrasts.")
  }
  times <- sort(unique(as.numeric(cells$time[!cells$is_baseline])))
  times <- times[is.finite(times)]
  if (length(times) < 2L) {
    tps_abort("Need at least 2 non-baseline timepoints; found %s.",
              length(times))
  }
  if (is.null(endpoint)) endpoint <- max(times)
  if (!is.numeric(endpoint) || length(endpoint) != 1L ||
      is.na(endpoint) || !is.finite(endpoint)) {
    tps_abort("endpoint must be NULL or one finite numeric timepoint.")
  }
  if (!endpoint %in% times) {
    tps_abort("endpoint = %s is not among non-baseline timepoints (%s).",
              endpoint, paste(times, collapse = ", "))
  }
  if (is.null(reference_condition)) {
    reference_condition <- conditions[1L]
  }
  if (!is.character(reference_condition) || length(reference_condition) != 1L ||
      is.na(reference_condition) || !nzchar(reference_condition)) {
    tps_abort("reference_condition must be one non-empty condition name.")
  }
  if (!reference_condition %in% conditions) {
    tps_abort("reference_condition '%s' not among conditions (%s).",
              reference_condition, paste(conditions, collapse = ", "))
  }

  families <- list()
  family_meta <- list()

  ## ---- A. vs_baseline -----------------------------------------------------
  fam_A <- .build_vs_baseline(cells, conditions, times, levels_d)
  families[["vs_baseline"]] <- fam_A$cmat
  family_meta[["vs_baseline"]] <- fam_A$meta

  ## ---- B. within_condition ------------------------------------------------
  fam_B <- .build_within_condition(cells, conditions, times, within, levels_d)
  families[["within_condition"]] <- fam_B$cmat
  family_meta[["within_condition"]] <- fam_B$meta

  ## ---- C. average_time ----------------------------------------------------
  if (include_average_time) {
    fam_C <- .build_average_time(cells, conditions, times, levels_d)
    families[["average_time"]] <- fam_C$cmat
    family_meta[["average_time"]] <- fam_C$meta
  }

  ## ---- D. heterogeneity ---------------------------------------------------
  if (include_heterogeneity && length(conditions) >= 2L) {
    fam_D <- .build_heterogeneity(
      cells, conditions, times, reference_condition, levels_d
    )
    families[["heterogeneity"]] <- fam_D$cmat
    family_meta[["heterogeneity"]] <- fam_D$meta
  } else if (include_heterogeneity) {
    tps_warn("heterogeneity family omitted because the design has one condition.")
  }

  ## ---- E. within_omnibus --------------------------------------------------
  if (include_within_omnibus) {
    fam_E <- .build_within_omnibus(cells, conditions, times, levels_d)
    families[["within_omnibus"]] <- fam_E$cmat
    family_meta[["within_omnibus"]] <- fam_E$meta
  }

  ## ---- F. endpoint --------------------------------------------------------
  fam_F <- .build_endpoint(cells, conditions, endpoint, levels_d)
  families[["endpoint"]] <- fam_F$cmat
  family_meta[["endpoint"]] <- fam_F$meta

  ## ---- rank audit per family ---------------------------------------------
  rank_audit <- lapply(families, function(cm) {
    if (is.null(cm) || length(cm) == 0L) return(list(rank = 0L, ncol = 0L))
    r <- qr(cm)$rank
    list(rank = r, ncol = ncol(cm), full_rank = (r == ncol(cm)))
  })

  ## ---- flattened contrast table ------------------------------------------
  contrast_table <- .flatten_contrast_table(families, family_meta)

  out <- structure(
    list(
      families = families,
      family_meta = family_meta,
      rank_audit = rank_audit,
      contrast_table = contrast_table,
      design_levels = levels_d,
      conditions = conditions,
      times = times,
      endpoint = endpoint,
      reference_condition = reference_condition,
      provenance = list(.make_provenance(arg_capture, design, stage = "contrasts"))
    ),
    class = "tps_contrasts"
  )
  out
}

## ---- family builders -------------------------------------------------------
## Each returns list(cmat = matrix (levels x contrasts), meta = data.frame).
## A contrast vector over levels; 0 for non-involved cells.
.mk_contrast <- function(levels_d, coeffs) {
  v <- setNames(numeric(length(levels_d)), levels_d)
  v[names(coeffs)] <- coeffs
  v
}

.build_vs_baseline <- function(cells, conditions, times, levels_d) {
  bl <- "baseline"
  rows <- list()
  for (cond in conditions) {
    for (tt in times) {
      cell <- paste(cond, tt, sep = "_")
      name <- sprintf("%s_t%s_vs_baseline", cond, tt)
      v <- .mk_contrast(levels_d, setNames(c(1, -1), c(cell, bl)))
      rows[[name]] <- v
    }
  }
  cmat <- do.call(cbind, rows)
  meta <- data.frame(contrast = colnames(cmat),
                      family = "vs_baseline", test = "t",
                      stringsAsFactors = FALSE)
  list(cmat = cmat, meta = meta)
}

.build_within_condition <- function(cells, conditions, times, within, levels_d) {
  rows <- list()
  tt_sorted <- sort(times)
  pairs <- if (within == "consecutive") {
    as.list(as.data.frame(t(embed(tt_sorted, 2))))  # consecutive pairs reversed
  } else {
    ## all pairs (t_i, t_j), i < j
    idx <- utils::combn(seq_along(tt_sorted), 2L)
    lapply(seq_len(ncol(idx)), function(k) tt_sorted[idx[, k]])
  }
  for (cond in conditions) {
    for (pr in pairs) {
      t_hi <- max(pr); t_lo <- min(pr)
      cell_hi <- paste(cond, t_hi, sep = "_")
      cell_lo <- paste(cond, t_lo, sep = "_")
      name <- sprintf("%s_t%s_vs_t%s", cond, t_hi, t_lo)
      v <- .mk_contrast(levels_d, setNames(c(1, -1), c(cell_hi, cell_lo)))
      rows[[name]] <- v
    }
  }
  cmat <- do.call(cbind, rows)
  meta <- data.frame(contrast = colnames(cmat),
                      family = "within_condition", test = "t",
                      stringsAsFactors = FALSE)
  list(cmat = cmat, meta = meta)
}

.build_average_time <- function(cells, conditions, times, levels_d) {
  C <- length(conditions)
  rows <- list()
  for (tt in times) {
    cells_t <- paste(conditions, tt, sep = "_")
    ## A_t = mean_c(mu_{c,t}) - mu_baseline
    coeffs <- c(setNames(rep(1 / C, C), cells_t), baseline = -1)
    name <- sprintf("average_t%s_vs_baseline", tt)
    rows[[name]] <- .mk_contrast(levels_d, coeffs)
  }
  cmat <- do.call(cbind, rows)
  meta <- data.frame(contrast = colnames(cmat),
                      family = "average_time", test = "F_omnibus",
                      stringsAsFactors = FALSE)
  list(cmat = cmat, meta = meta)
}

.build_heterogeneity <- function(cells, conditions, times, ref, levels_d) {
  ## Full-rank: for each non-ref condition, for each time, mu_{c,t} - mu_{ref,t}.
  ## (C-1) * length(times) contrasts, rank = (C-1)*T.
  non_ref <- setdiff(conditions, ref)
  rows <- list()
  for (cond in non_ref) {
    for (tt in times) {
      cell_c <- paste(cond, tt, sep = "_")
      cell_r <- paste(ref, tt, sep = "_")
      name <- sprintf("%s_minus_%s_t%s", cond, ref, tt)
      v <- .mk_contrast(levels_d, setNames(c(1, -1), c(cell_c, cell_r)))
      rows[[name]] <- v
    }
  }
  cmat <- do.call(cbind, rows)
  meta <- data.frame(contrast = colnames(cmat),
                      family = "heterogeneity", test = "F_omnibus",
                      reference = ref, stringsAsFactors = FALSE)
  list(cmat = cmat, meta = meta)
}

.build_within_omnibus <- function(cells, conditions, times, levels_d) {
  ## For each condition, all vs-baseline contrasts jointly -> omnibus F.
  ## We split the vs-baseline family by condition into per-condition blocks.
  rows <- list()
  meta_rows <- list()
  for (cond in conditions) {
    cells_c <- paste(cond, times, sep = "_")
    for (cell in cells_c) {
      name <- sprintf("within_omni_%s", cell)
      v <- .mk_contrast(levels_d, setNames(c(1, -1), c(cell, "baseline")))
      rows[[name]] <- v
      meta_rows[[name]] <- data.frame(contrast = name, family = "within_omnibus",
                                       test = "F_omnibus", condition = cond,
                                       stringsAsFactors = FALSE)
    }
  }
  cmat <- do.call(cbind, rows)
  meta <- do.call(rbind, meta_rows)
  list(cmat = cmat, meta = meta)
}

.build_endpoint <- function(cells, conditions, endpoint, levels_d) {
  rows <- list()
  for (cond in conditions) {
    cell <- paste(cond, endpoint, sep = "_")
    name <- sprintf("endpoint_%s_t%s_vs_baseline", cond, endpoint)
    v <- .mk_contrast(levels_d, setNames(c(1, -1), c(cell, "baseline")))
    rows[[name]] <- v
  }
  cmat <- do.call(cbind, rows)
  meta <- data.frame(contrast = colnames(cmat),
                      family = "endpoint", test = "t",
                      stringsAsFactors = FALSE)
  list(cmat = cmat, meta = meta)
}

## ---- flatten --------------------------------------------------------------
.flatten_contrast_table <- function(families, family_meta) {
  rows <- list()
  for (fn in names(families)) {
    cm <- families[[fn]]
    mm <- family_meta[[fn]]
    if (is.null(cm) || length(cm) == 0L) next
    for (k in seq_len(ncol(cm))) {
      rows[[length(rows) + 1L]] <- data.frame(
        family = fn,
        contrast = colnames(cm)[k],
        test = mm$test[k],
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

## ---- methods ---------------------------------------------------------------
#' @export
print.tps_contrasts <- function(x, ...) {
  cat("<timeProtState contrasts>\n")
  cat(sprintf("  conditions: %s\n", paste(x$conditions, collapse = ", ")))
  cat(sprintf("  times: %s   endpoint: %s\n",
              paste(x$times, collapse = ","), x$endpoint))
  cat(sprintf("  reference_condition: %s\n", x$reference_condition))
  cat("  families:\n")
  for (fn in names(x$families)) {
    ra <- x$rank_audit[[fn]]
    cat(sprintf("    %-18s %d contrast(s), rank %d (full_rank: %s)\n",
                fn, ra$ncol, ra$rank, ra$full_rank))
  }
  invisible(x)
}
