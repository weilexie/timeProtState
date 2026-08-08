## ============================================================================
## Pathway scores: NET / REM / COH (deterministic, no permutation)
## ============================================================================
##
## Implements STATISTICAL_CONTRACT.md section 8.1-8.3 and DECISIONS D5.
## This is the deterministic scoring layer (Goal 11). Permutation calibration
## lives in pathway-null.R (Goal 12).
##
## For pathway g and contrast j, with protein-level statistic s_{i,j}:
##   NET[g,j] = sum_i M[g,i] * s[i,j] / size[g]      (net direction)
##   REM[g,j] = sum_i M[g,i] * |s[i,j]| / size[g]    (total remodeling magnitude)
##   COH[g,j] = |NET| / (REM + eps)                   (directional coherence, [0,1])
##
## Default statistic is the moderated t (carries direction + uncertainty).
## Faithful to the validated spike scripts (pathway_remodeling_spike*.R).

## ---- membership matrix construction ---------------------------------------
.build_membership <- function(term2gene, universe, min_size) {
  ## universe: character vector of feature ids that appear in the statistic
  ## matrix (rownames). term2gene: data.frame with columns term, gene.
  if (!all(c("term","gene") %in% names(term2gene))) {
    tps_abort("term2gene must have columns 'term' and 'gene'.")
  }
  pathways <- unique(term2gene$term)
  n_path <- length(pathways)
  n_prot <- length(universe)
  M <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                             dims = c(n_path, n_prot),
                             dimnames = list(pathways, universe))
  observed_size <- integer(n_path)
  library_size <- integer(n_path)
  for (p in seq_len(n_path)) {
    members_lib <- term2gene$gene[term2gene$term == pathways[p]]
    library_size[p] <- length(unique(members_lib))
    members_obs <- intersect(unique(members_lib), universe)
    observed_size[p] <- length(members_obs)
    if (length(members_obs) > 0L) {
      j <- match(members_obs, universe)
      M[p, j] <- 1
    }
  }
  keep <- observed_size >= min_size
  M <- M[keep, , drop = FALSE]
  attr(M, "library_size") <- library_size[keep]
  attr(M, "observed_size") <- observed_size[keep]
  attr(M, "dropped") <- pathways[!keep]
  M
}

## ---- tps_pathway_scores ----------------------------------------------------
#' Compute deterministic pathway NET/REM/COH scores
#'
#' Computes the directional net (NET), remodeling magnitude (REM), and
#' directional coherence (COH) for each pathway x contrast, from a moderated
#' statistic matrix (default: moderated t from \code{\link{tps_fit}}). This is
#' the deterministic layer; permutation calibration is applied separately via
#' \code{\link{tps_calibrate_pathways}}.
#'
#' Pathways with fewer than \code{min_size} observed members (intersected with
#' the analysis universe) are dropped.
#'
#' @param fit a \code{tps_fit} object, OR a numeric statistic matrix
#'   (features x contrasts).
#' @param term2gene data.frame with columns \code{term} and \code{gene}. The
#'   \code{gene} values must be matchable to the fit feature ids (see
#'   \code{id_col}).
#' @param type character; statistic to aggregate: \code{"moderated_t"} (default),
#'   \code{"logFC"}, or \code{"signed_log10p"}.
#' @param family character; if \code{fit} is a \code{tps_fit}, restrict to this
#'   contrast family (e.g. \code{"endpoint"}). Ignored if a matrix is supplied.
#' @param id_col character; if feature ids in the fit come from a rowData column
#'   (e.g. locus tags) rather than row names, name that column. If \code{NULL},
#'   row names of the statistic matrix are used.
#' @param min_size integer; minimum observed pathway size to retain.
#' @param min_coverage numeric in [0,1]; minimum library coverage
#'   (observed/library) to retain. Default 0 (no coverage filter beyond size).
#' @return a list of class \code{tps_pathway_scores} with:
#'   \item{scores}{list with NET, REM, COH matrices (pathways x contrasts).}
#'   \item{membership}{the sparse membership matrix (with size attributes).}
#'   \item{pathway_table}{per-pathway size/coverage audit.}
#'   \item{statistic_type}{the statistic used.}
#'   \item{provenance}{provenance entry.}
#'
#' @examples
#' se <- tps_import(tps_example_abundance, tps_example_metadata,
#'                  tps_example_annotation)
#' design <- tps_make_design(se)
#' contrasts <- tps_make_contrasts(design, endpoint = 7)
#' fit <- tps_fit(se, design, contrasts)
#' pws <- tps_pathway_scores(fit, tps_example_term2gene,
#'                            id_col = "locus_tag", family = "endpoint")
#' dim(pws$scores$NET)
#'
#' @export
tps_pathway_scores <- function(fit, term2gene, type = "moderated_t",
                                family = NULL, id_col = NULL,
                                min_size = 5L, min_coverage = 0) {
  arg_capture <- list(call = "tps_pathway_scores", type = type, family = family,
                       id_col = id_col, min_size = min_size,
                       min_coverage = min_coverage,
                       timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  eps <- 1e-12
  if (min_size < 1L) {
    tps_abort("min_size must be >= 1 (got %d); pathways with 0 members would ",
              "cause division-by-zero in NET/REM scores.", min_size)
  }

  ## obtain statistic matrix
  if (inherits(fit, "tps_fit")) {
    stat_mat <- tps_statistic_matrix(fit, type = type, family = family)
    ## if id_col given, remap rownames via the annotation carried by the fit
    if (!is.null(id_col)) {
      ann <- fit$annotation
      ## annotation rownames align with stat_mat rownames (both = feature ids)
      if (is.null(ann)) {
        tps_abort(paste0("fit carries no annotation; cannot remap via id_col='",
                          id_col, "'. Pass a pre-mapped statistic matrix or omit id_col."))
      }
      if (!id_col %in% names(ann)) {
        tps_abort("annotation has no column '%s'.", id_col)
      }
      ## ann rownames == stat_mat rownames; remap
      new_ids <- ann[[id_col]][match(rownames(stat_mat), rownames(ann))]
      unmapped <- is.na(new_ids)
      new_ids[unmapped] <- rownames(stat_mat)[unmapped]
      rownames(stat_mat) <- new_ids
    }
  } else if (is.matrix(fit) && is.numeric(fit)) {
    stat_mat <- fit
    type <- "user_matrix"
  } else {
    tps_abort("fit must be a tps_fit object or a numeric matrix.")
  }

  universe <- rownames(stat_mat)
  if (is.null(universe)) {
    tps_abort("statistic matrix must have row names (feature ids).")
  }

  M <- .build_membership(term2gene, universe, min_size = 0L)
  lib_size <- attr(M, "library_size")
  obs_size <- attr(M, "observed_size")

  ## coverage filter
  if (min_coverage > 0) {
    cov <- obs_size / pmax(lib_size, 1L)
    keep_cov <- cov >= min_coverage
    M <- M[keep_cov, , drop = FALSE]
    lib_size <- lib_size[keep_cov]
    obs_size <- obs_size[keep_cov]
  }
  ## min_size filter (applied after coverage for correct counts)
  keep_size <- obs_size >= min_size
  M <- M[keep_size, , drop = FALSE]
  lib_size <- lib_size[keep_size]
  obs_size <- obs_size[keep_size]

  if (nrow(M) == 0L) {
    tps_abort(paste0("No pathways retained after min_size/min_coverage filtering; ",
                     "check term2gene mapping and id_col. universe size = ",
                     length(universe), ", universe ids e.g. ",
                     paste(head(universe, 3), collapse = ", "), "; term2gene genes e.g. ",
                     paste(head(unique(term2gene$gene), 3), collapse = ", "), "."))
  }

  sizes <- obs_size
  Mt <- as(M, "CsparseMatrix")
  NET <- as.matrix(Mt %*% stat_mat) / sizes
  REM <- as.matrix(Mt %*% abs(stat_mat)) / sizes
  COH <- abs(NET) / (REM + eps)

  pathway_table <- data.frame(
    term = rownames(M),
    library_size = lib_size,
    observed_size = obs_size,
    coverage = obs_size / pmax(lib_size, 1L),
    stringsAsFactors = FALSE)

  out <- structure(
    list(
      scores = list(NET = NET, REM = REM, COH = COH),
      membership = M,
      pathway_table = pathway_table,
      statistic_matrix = stat_mat,
      statistic_type = type,
      universe = universe,
      min_size = min_size,
      min_coverage = min_coverage,
      dropped_pathways = attr(.build_membership(term2gene, universe, 0L), "dropped"),
      provenance = list(.make_provenance(arg_capture, fit, stage = "pathway_scores"))
    ),
    class = "tps_pathway_scores"
  )
  out
}

## ---- methods ---------------------------------------------------------------
#' @export
print.tps_pathway_scores <- function(x, ...) {
  cat("<timeProtState pathway scores>\n")
  cat(sprintf("  pathways: %d   contrasts: %d   statistic: %s\n",
              nrow(x$scores$NET), ncol(x$scores$NET), x$statistic_type))
  pt <- x$pathway_table
  os <- as.integer(pt$observed_size)
  cat(sprintf("  observed size: min=%s median=%s max=%s\n",
              min(os), as.integer(stats::median(os)), max(os)))
  cat(sprintf("  coverage: min=%.2f median=%.2f\n",
              min(pt$coverage), stats::median(pt$coverage)))
  if (length(x$dropped_pathways)) {
    cat(sprintf("  dropped (< min_size %d): %d\n",
                x$min_size, length(x$dropped_pathways)))
  }
  invisible(x)
}
