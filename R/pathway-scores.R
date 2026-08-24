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
  if (!is.data.frame(term2gene) ||
      !all(c("term", "gene") %in% names(term2gene))) {
    tps_abort("term2gene must have columns 'term' and 'gene'.")
  }
  term2gene <- unique(data.frame(
    term = as.character(term2gene$term),
    gene = as.character(term2gene$gene),
    stringsAsFactors = FALSE
  ))
  valid <- !is.na(term2gene$term) & nzchar(term2gene$term) &
    !is.na(term2gene$gene) & nzchar(term2gene$gene)
  term2gene <- term2gene[valid, , drop = FALSE]
  if (nrow(term2gene) == 0L) {
    tps_abort("term2gene has no non-missing term-gene mappings.")
  }
  pathways <- unique(term2gene$term)
  n_path <- length(pathways)
  n_prot <- length(universe)
  M <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                             x = numeric(0),
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
#' The scores are: \code{NET = mean(t)}, \code{REM = mean(|t|)},
#' \code{COH = |NET| / (REM + eps)}. COH is a magnitude-weighted coherence
#' measure in [0,1]: it equals 1 when all members share the same direction
#' and approaches 0 when opposing directions cancel. This is not a simple
#' count-based fraction of same-sign members; it weights by magnitude.
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
#'   contrast family (e.g. \code{"vs_baseline"} or \code{"endpoint"}). If
#'   \code{NULL} and exactly one single-contrast (t) family is present, that
#'   family is used. If \code{NULL} and multiple t-families are present, the
#'   call aborts because pooling can double-count algebraically identical
#'   contrasts (e.g. endpoint vs vs_baseline at the endpoint timepoint).
#'   Ignored if a matrix is supplied.
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
#' data(tps_example, package = "timeProtState")
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
  if (!is.numeric(min_size) || length(min_size) != 1L || is.na(min_size) ||
      !is.finite(min_size) || min_size < 1L ||
      min_size != as.integer(min_size)) {
    tps_abort("min_size must be one positive integer.")
  }
  if (!is.numeric(min_coverage) || length(min_coverage) != 1L ||
      is.na(min_coverage) || !is.finite(min_coverage) ||
      min_coverage < 0 || min_coverage > 1) {
    tps_abort("min_coverage must be in [0,1].")
  }
  if (!is.null(id_col) &&
      (!is.character(id_col) || length(id_col) != 1L ||
       is.na(id_col) || !nzchar(id_col))) {
    tps_abort("id_col must be NULL or one non-empty column name.")
  }

  ## obtain statistic matrix
  pathway_ids <- NULL
  residual_mat <- NULL
  if (inherits(fit, "tps_fit")) {
    ## T01/B5: family=NULL pools every single-contrast (t) family. In a
    ## shared-baseline design, endpoint is algebraically a subset of
    ## vs_baseline at the endpoint timepoint, so pooling double-counts the
    ## same contrast in NET/REM. Require an explicit family when more than
    ## one t-family is present.
    if (is.null(family)) {
      ct <- fit$contrasts$contrast_table
      t_families <- if (is.data.frame(ct) && "test" %in% names(ct) &&
                        "family" %in% names(ct)) {
        unique(as.character(ct$family[ct$test == "t"]))
      } else character(0)
      if (length(t_families) > 1L) {
        tps_abort(paste0(
          "family=NULL would pool ", length(t_families),
          " single-contrast families (", paste(t_families, collapse = ", "),
          "). This can double-count algebraically identical contrasts ",
          "(e.g. endpoint vs vs_baseline at the endpoint timepoint). ",
          "Specify family= explicitly (commonly 'vs_baseline' or 'endpoint')."
        ))
      }
      if (length(t_families) == 0L) {
        tps_abort("fit has no single-contrast (t) families to score.")
      }
    }
    stat_mat <- tps_statistic_matrix(fit, type = type, family = family)
    ## T01: detect algebraically duplicate contrasts (identical statistic
    ## columns). Within a single well-formed family this should not occur;
    ## a hit signals a family-construction or pooling bug.
    if (ncol(stat_mat) >= 2L) {
      dup_cols <- anyDuplicated(as.data.frame(stat_mat))
      if (dup_cols > 0L) {
        tps_warn(paste0(
          "statistic matrix has duplicate columns at position ", dup_cols,
          " and ", dup_cols + 1L, "; two contrasts are algebraically ",
          "identical. Check family selection and contrast construction."
        ))
      }
    }
    if (!is.null(fit$residuals)) {
      residual_mat <- fit$residuals
    }
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
      new_ids <- as.character(new_ids)
      unmapped <- is.na(new_ids)
      new_ids[unmapped] <- rownames(stat_mat)[unmapped]
      if (anyDuplicated(new_ids)) {
        dup <- unique(new_ids[duplicated(new_ids)])
        tps_abort(paste0(
          "id_col='", id_col, "' maps multiple features to the same id (e.g. ",
          paste(head(dup, 3L), collapse = ", "),
          "). Resolve or aggregate duplicate mappings before pathway scoring."
        ))
      }
      pathway_ids <- new_ids
    }
  } else if (is.matrix(fit) && is.numeric(fit)) {
    stat_mat <- fit
    type <- "user_matrix"
  } else {
    tps_abort("fit must be a tps_fit object or a numeric matrix.")
  }

  feature_ids <- rownames(stat_mat)
  if (is.null(feature_ids)) {
    tps_abort("statistic matrix must have row names (feature ids).")
  }
  if (anyDuplicated(feature_ids)) {
    tps_abort("statistic matrix row names must be unique feature ids.")
  }
  if (is.null(pathway_ids)) pathway_ids <- feature_ids

  ## Keep one common, finite feature universe across all contrasts. Treating an
  ## unestimable statistic as zero would silently turn missing evidence into a
  ## null effect and dilute pathway scores.
  finite_rows <- apply(is.finite(stat_mat), 1L, all)
  dropped_nonfinite <- feature_ids[!finite_rows]
  stat_mat <- stat_mat[finite_rows, , drop = FALSE]
  pathway_ids <- pathway_ids[finite_rows]
  if (!is.null(residual_mat)) {
    rid <- feature_ids[finite_rows]
    if (!is.matrix(residual_mat) || !all(rid %in% rownames(residual_mat))) {
      residual_mat <- NULL
    } else {
      residual_mat <- residual_mat[rid, , drop = FALSE]
    }
  }
  universe <- pathway_ids
  if (length(universe) == 0L) {
    tps_abort("statistic matrix has no features finite across all selected contrasts.")
  }

  ## T03/M2: carry AveExpr from the fit so tps_calibrate_pathways can derive
  ## abundance_stratified bins without the user manually rebuilding it.
  ave_expr <- NULL
  if (inherits(fit, "tps_fit") && "AveExpr" %in% names(fit$statistics)) {
    st <- fit$statistics[fit$statistics$test == "t", , drop = FALSE]
    if (nrow(st) > 0L) {
      ## one AveExpr per feature; AveExpr is feature-level (same across contrasts)
      ave_df <- unique(st[, c("feature_id", "AveExpr")])
      ave_df <- ave_df[match(rownames(stat_mat), ave_df$feature_id), ,
                       drop = FALSE]
      if (nrow(ave_df) == nrow(stat_mat) && all(!is.na(ave_df$AveExpr))) {
        ave_expr <- setNames(ave_df$AveExpr, ave_df$feature_id)
      }
    }
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

  ## T03/M2: attrition/coverage audit — counts from input features down to
  ## the final pathway universe so reviewers can trace why N differs per step.
  attrition <- data.frame(
    step = c("input_features", "nonfinite_dropped", "in_universe",
             "pathways_in_term2gene", "pathways_after_filters"),
    count = c(length(feature_ids),
             length(dropped_nonfinite),
             length(universe),
             length(unique(as.character(term2gene$term))),
             nrow(M)),
    stringsAsFactors = FALSE)

  ## T02/M1: carry contrast metadata from the fit's contrast_table so
  ## downstream consumers (pathway-dynamics, plotting) read condition/time
  ## from the object instead of re-parsing contrast name strings. NULL for
  ## bare matrices (no upstream contrast construction).
  contrast_meta <- NULL
  if (inherits(fit, "tps_fit") && !is.null(fit$contrasts) &&
      is.data.frame(fit$contrasts$contrast_table)) {
    ct <- fit$contrasts$contrast_table
    keep_t <- ct$test == "t"
    if (is.null(family)) {
      ## T01 guard already ensured exactly one t-family; use it.
      cm_sub <- ct[keep_t, , drop = FALSE]
    } else {
      cm_sub <- ct[keep_t & ct$family == family, , drop = FALSE]
    }
    ## align row order to stat_mat columns
    cm_sub <- cm_sub[match(colnames(stat_mat), cm_sub$contrast), ,
                     drop = FALSE]
    if (nrow(cm_sub) == ncol(stat_mat)) {
      contrast_meta <- cm_sub
    }
  }

  all_pathways <- unique(as.character(term2gene$term))
  all_pathways <- all_pathways[!is.na(all_pathways) & nzchar(all_pathways)]
  out <- structure(
    list(
      scores = list(NET = NET, REM = REM, COH = COH),
      membership = M,
      pathway_table = pathway_table,
      statistic_matrix = stat_mat,
      residuals = residual_mat,
      statistic_type = type,
      contrast_meta = contrast_meta,
      universe = universe,
      feature_ids = rownames(stat_mat),
      id_mapping = data.frame(feature_id = rownames(stat_mat),
                              pathway_id = pathway_ids,
                              stringsAsFactors = FALSE),
      min_size = min_size,
      min_coverage = min_coverage,
      abundance = ave_expr,
      attrition = attrition,
      ## T04/M3: carry the fit's input fingerprint so calibrate/dynamics can
      ## detect stale reuse (e.g. a calibrated object from an old fit reused
      ## with a new design). NULL for bare matrices.
      input_fingerprint = if (inherits(fit, "tps_fit")) fit$input_fingerprint else NULL,
      dropped_pathways = setdiff(all_pathways, rownames(M)),
      dropped_features_nonfinite = dropped_nonfinite,
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
    cat(sprintf("  pathways dropped by size/coverage filters: %d\n",
                length(x$dropped_pathways)))
  }
  if (length(x$dropped_features_nonfinite)) {
    cat(sprintf("  features dropped (non-finite statistic): %d\n",
                length(x$dropped_features_nonfinite)))
  }
  invisible(x)
}
