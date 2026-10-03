## ============================================================================
## Pathway calibration: CAMERA / ROMER significance + classification
## ============================================================================
##
## Implements STATISTICAL_CONTRACT.md, DECISIONS D5, and ADR 0011.
## Default significance is limma competitive gene-set tests (CAMERA, or
## ROMER when a repeated-measures block is present). NET/REM/COH are
## summaries computed in tps_pathway_scores and are not recomputed here.
## Permutation+VIF (ADR 0010) remains an explicit diagnostic engine.
##
## Significance columns: p_directional / q_directional and p_mixed / q_mixed.
## COH has no p-value.

## ---- stratification helpers ------------------------------------------------
.make_abundance_bins <- function(abundance, n_bins) {
  ## equal-frequency bins over abundance (AveExpr); returns integer bin per protein
  if (anyNA(abundance)) {
    tps_abort("abundance must not contain NA for bin assignment.")
  }
  if (length(unique(abundance)) < n_bins) {
    ## too few distinct values; fall back to as many bins as distinct values
    n_bins <- length(unique(abundance))
  }
  br <- stats::quantile(abundance, probs = seq(0, 1, length.out = n_bins + 1),
                        na.rm = TRUE)
  br <- unique(br)  # collapse ties
  if (length(br) < 2L) return(rep(1L, length(abundance)))
  as.integer(cut(abundance, breaks = br, include.lowest = TRUE, labels = FALSE))
}

.adjust_p_matrix <- function(p, correction) {
  out <- if (correction == "per_contrast") {
    vapply(seq_len(ncol(p)), function(j) {
      stats::p.adjust(p[, j], method = "BH")
    }, numeric(nrow(p)))
  } else {
    matrix(stats::p.adjust(as.vector(p), method = "BH"),
           nrow = nrow(p), ncol = ncol(p))
  }
  if (!is.matrix(out)) out <- matrix(out, nrow = nrow(p))
  dimnames(out) <- dimnames(p)
  out
}

## generate a single permutation index under a given mode
## - unstratified: global sample.int(n)
## - abundance_stratified: shuffle within each bin, concatenate in bin-template order
.perm_index <- function(n_prot, mode, bin_list, bin_template) {
  if (mode == "unstratified") {
    return(sample.int(n_prot))
  }
  if (mode == "abundance_stratified") {
    perm_per_bin <- lapply(bin_list, function(ii) ii[sample.int(length(ii))])
    perm <- unlist(perm_per_bin, use.names = FALSE)
    idx <- integer(n_prot)
    idx[bin_template] <- perm
    return(idx)
  }
  tps_abort("unknown permutation mode: %s", mode)
}

## Mean pairwise correlation of member *residuals across samples* after
## the group-means fit. Residual correlation is the competitive extra
## covariance that label permutation destroys. Do not estimate rho from
## moderated-t profiles: a shared pathway effect would be mistaken for
## noise correlation and the real NET would be over-corrected.
##
## Clipping-bias debiasing (ADR 0010, 2026-08-25): the clipped estimator
## max(0, mean pairwise rho) is positively biased at true rho = 0 — half
## the estimation-noise mass is floored at 0 while the other half stays
## positive, so E[rho_hat | rho = 0] > 0. With k members the bias is
## amplified by (k-1) in the VIF, inflating every pathway's VIF above 1
## even on independent data and systematically shrinking real signal.
## The bias floor is estimated directly by drawing size-matched random
## gene sets from the same residual matrix (permutation of membership is
## exactly the null in question) and subtracted from each pathway's rho.
.estimate_pathway_rho <- function(resid_mat, membership, n_floor_draws = 50L,
                                   seed = 1L) {
  n_path <- nrow(membership)
  rho <- rep(0, n_path)
  if (is.null(resid_mat) || !is.matrix(resid_mat) ||
      ncol(resid_mat) < 4L || n_path < 1L) {
    return(rho)
  }
  if (nrow(resid_mat) != ncol(membership)) {
    return(rho)
  }
  .mean_pairwise_rho <- function(idx) {
    S <- resid_mat[idx, , drop = FALSE]
    keep <- apply(is.finite(S), 1L, function(r) sum(r) >= 4L)
    S <- S[keep, , drop = FALSE]
    if (nrow(S) < 3L) return(NA_real_)
    C <- suppressWarnings(stats::cor(t(S), use = "pairwise.complete.obs"))
    off <- C[upper.tri(C)]
    off <- off[is.finite(off)]
    if (length(off)) mean(off) else NA_real_
  }
  for (p in seq_len(n_path)) {
    idx <- which(as.numeric(membership[p, ]) != 0)
    if (length(idx) < 3L) next
    r <- .mean_pairwise_rho(idx)
    rho[p] <- if (is.finite(r)) r else 0
  }
  ## clipping-bias floor: mean rho over size-matched random sets
  n_prot <- nrow(resid_mat)
  sizes_raw <- as.integer(rowSums(as.matrix(membership != 0)))
  floor_by_size <- list()
  withr::with_seed(seed, {
    for (p in seq_len(n_path)) {
      k <- sizes_raw[p]
      if (k < 3L || k >= n_prot) next
      if (is.null(floor_by_size[[as.character(k)]])) {
        draws <- vapply(seq_len(n_floor_draws), function(i) {
          r <- .mean_pairwise_rho(sample.int(n_prot, k))
          if (is.finite(r)) r else NA_real_
        }, numeric(1L))
        draws <- draws[is.finite(draws)]
        floor_by_size[[as.character(k)]] <-
          if (length(draws)) mean(draws) else 0
      }
    }
  })
  rho_deb <- vapply(seq_len(n_path), function(p) {
    fl <- floor_by_size[[as.character(sizes_raw[p])]]
    if (is.null(fl)) fl <- 0
    max(0, rho[p] - fl)
  }, numeric(1L))
  pmax(0, pmin(rho_deb, 0.99))
}

.pathway_vif <- function(rho, sizes) {
  1 + pmax(as.numeric(sizes) - 1, 0) * pmax(as.numeric(rho), 0)
}

## ---- CAMERA / ROMER / FRY helpers -----------------------------------------
.limma_default_formals <- function(generic) {
  ns <- asNamespace("limma")
  meth <- paste0(generic, ".default")
  fun <- if (exists(meth, ns, inherits = FALSE)) {
    get(meth, ns, inherits = FALSE)
  } else {
    getExportedValue("limma", generic)
  }
  names(formals(fun))
}

.camera_has_directional <- function() {
  fmls <- unique(c(names(formals(limma::camera)),
                   .limma_default_formals("camera")))
  "directional" %in% fmls
}

.empty_p_matrix <- function(terms, contrasts) {
  matrix(NA_real_, length(terms), length(contrasts),
         dimnames = list(terms, contrasts))
}

.match_set_stat <- function(tab, terms, column) {
  out <- setNames(rep(NA_real_, length(terms)), terms)
  if (is.null(tab) || is.null(dim(tab)) || is.null(rownames(tab))) {
    return(out)
  }
  if (!column %in% colnames(tab)) return(out)
  ii <- match(terms, rownames(tab))
  ok <- !is.na(ii)
  out[ok] <- as.numeric(tab[ii[ok], column])
  out
}

.match_set_chr <- function(tab, terms, column) {
  out <- setNames(rep(NA_character_, length(terms)), terms)
  if (is.null(tab) || is.null(dim(tab)) || is.null(rownames(tab))) {
    return(out)
  }
  if (!column %in% colnames(tab)) return(out)
  ii <- match(terms, rownames(tab))
  ok <- !is.na(ii)
  out[ok] <- as.character(tab[ii[ok], column])
  out
}

.run_camera <- function(y, index, design, contrast, use.ranks, directional,
                         trend) {
  args <- list(
    y = y, index = index, design = design, contrast = contrast,
    sort = FALSE, inter.gene.cor = NA, trend.var = isTRUE(trend),
    use.ranks = isTRUE(use.ranks)
  )
  if (.camera_has_directional()) {
    args$directional <- isTRUE(directional)
  }
  do.call(limma::camera, args)
}

.derive_test_direction <- function(p_up, p_down, hint = NULL) {
  td <- matrix(NA_character_, nrow(p_up), ncol(p_up), dimnames = dimnames(p_up))
  both <- is.finite(p_up) & is.finite(p_down)
  td[both & p_up < p_down] <- "Up"
  td[both & p_down < p_up] <- "Down"
  if (!is.null(hint)) {
    fill <- is.na(td) & !is.na(hint) & hint %in% c("Up", "Down")
    td[fill] <- hint[fill]
  }
  td
}

.net_direction_matrix <- function(NET) {
  nd <- matrix(NA_character_, nrow(NET), ncol(NET), dimnames = dimnames(NET))
  ok <- is.finite(NET)
  nd[ok & NET > 0] <- "Up"
  nd[ok & NET < 0] <- "Down"
  nd
}

.direction_concordant <- function(test_direction, net_direction) {
  dc <- matrix(NA, nrow(net_direction), ncol(net_direction),
               dimnames = dimnames(net_direction))
  if (!is.matrix(test_direction) ||
      !identical(dim(test_direction), dim(net_direction))) {
    return(dc)
  }
  ok <- !is.na(test_direction) & !is.na(net_direction)
  dc[ok] <- test_direction[ok] == net_direction[ok]
  dc
}

.pathway_index <- function(scores, terms) {
  index <- scores$index
  if (is.null(index)) {
    index <- .index_from_membership(scores$membership)
  }
  if (is.null(names(index)) || !all(terms %in% names(index))) {
    index <- .index_from_membership(scores$membership)
  }
  index <- index[terms]
  if (any(vapply(index, function(z) length(z) < 1L, logical(1L)))) {
    tps_abort("CAMERA/ROMER index has an empty pathway; check membership.")
  }
  n_feat <- nrow(scores$expression)
  bad <- vapply(index, function(z) {
    any(!is.finite(z) | z < 1L | z > n_feat)
  }, logical(1L))
  if (any(bad)) {
    tps_abort("pathway index is not aligned with expression rows.")
  }
  index
}

.require_competitive_inputs <- function(scores, engine_used) {
  if (is.null(scores$expression) || !is.matrix(scores$expression) ||
      is.null(scores$design_matrix) || !is.matrix(scores$design_matrix)) {
    tps_abort(paste0(
      "engine='", engine_used, "' requires scores$expression and ",
      "scores$design_matrix. Bare statistic-matrix scores only support ",
      "engine='permutation'."
    ))
  }
  if (is.null(scores$contrast_matrix) || !is.matrix(scores$contrast_matrix)) {
    tps_abort(paste0(
      "engine='", engine_used, "' requires scores$contrast_matrix ",
      "(build scores from a tps_fit, not a bare statistic matrix)."
    ))
  }
  if (nrow(scores$design_matrix) != ncol(scores$expression)) {
    tps_abort("design_matrix rows (%d) != expression columns (%d).",
              nrow(scores$design_matrix), ncol(scores$expression))
  }
  if (nrow(scores$expression) != ncol(scores$membership)) {
    tps_abort("expression rows (%d) != membership columns (%d).",
              nrow(scores$expression), ncol(scores$membership))
  }
  NET <- scores$scores$NET
  need <- colnames(NET)
  cm <- scores$contrast_matrix
  if (is.null(need) || is.null(colnames(cm)) || !all(need %in% colnames(cm))) {
    tps_abort(paste0(
      "engine='", engine_used,
      "' requires contrast_matrix columns to match colnames of NET."
    ))
  }
  if (identical(engine_used, "romer") && !is.null(scores$block_value)) {
    if (length(scores$block_value) != ncol(scores$expression)) {
      tps_abort("block_value length (%d) != expression columns (%d).",
                length(scores$block_value), ncol(scores$expression))
    }
    if (is.null(scores$correlation) || length(scores$correlation) != 1L ||
        !is.finite(scores$correlation)) {
      tps_abort(paste0(
        "engine='romer' with a block requires a finite consensus ",
        "correlation from tps_fit."
      ))
    }
  }
  invisible(TRUE)
}

.contrast_vector <- function(contrast_matrix, contrast_name) {
  as.numeric(contrast_matrix[, contrast_name, drop = TRUE])
}

.finalize_calibrated <- function(scores, arg_capture, fields) {
  new_prov <- .make_provenance(arg_capture, scores, stage = "calibrate_pathways")
  old_prov <- if (is.list(scores$provenance)) scores$provenance else list()
  out <- scores
  for (nm in names(fields)) {
    out[[nm]] <- fields[[nm]]
  }
  n_old <- length(old_prov)
  combined_prov <- old_prov
  combined_prov[[n_old + 1L]] <- new_prov
  out$provenance <- combined_prov
  class(out) <- c("tps_pathway_calibrated", class(scores))
  out
}

.calibrate_competitive <- function(scores, engine_used, correction, nrot, seed,
                                    arg_capture) {
  .require_competitive_inputs(scores, engine_used)
  y <- as.matrix(scores$expression)
  design <- as.matrix(scores$design_matrix)
  NET <- scores$scores$NET
  terms <- rownames(NET)
  contrasts <- colnames(NET)
  cm <- scores$contrast_matrix[, contrasts, drop = FALSE]
  index <- .pathway_index(scores, terms)
  trend <- scores$trend

  p_directional <- .empty_p_matrix(terms, contrasts)
  p_mixed <- .empty_p_matrix(terms, contrasts)
  p_up <- .empty_p_matrix(terms, contrasts)
  p_down <- .empty_p_matrix(terms, contrasts)
  test_hint <- matrix(NA_character_, length(terms), length(contrasts),
                      dimnames = list(terms, contrasts))

  for (j in seq_along(contrasts)) {
    cn <- contrasts[[j]]
    cj <- .contrast_vector(cm, cn)
    if (identical(engine_used, "camera")) {
      cam_dir <- tryCatch(
        .run_camera(y, index, design, cj, use.ranks = FALSE,
                    directional = TRUE, trend = trend),
        error = function(e) e
      )
      if (inherits(cam_dir, "error")) {
        tps_abort(paste0("camera() failed for contrast '", cn, "': ",
                         conditionMessage(cam_dir)))
      }
      cam_mix <- tryCatch(
        .run_camera(y, index, design, cj, use.ranks = TRUE,
                    directional = FALSE, trend = trend),
        error = function(e) e
      )
      if (inherits(cam_mix, "error")) {
        tps_abort(paste0("camera() failed for contrast '", cn, "': ",
                         conditionMessage(cam_mix)))
      }
      p_directional[, j] <- .match_set_stat(cam_dir, terms, "PValue")
      p_mixed[, j] <- .match_set_stat(cam_mix, terms, "PValue")
      test_hint[, j] <- .match_set_chr(cam_dir, terms, "Direction")
    } else if (identical(engine_used, "romer")) {
      rom_args <- list(y = y, index = index, design = design,
                       contrast = cj, nrot = as.integer(nrot))
      if (!is.null(scores$block_value)) {
        rom_args$block <- scores$block_value
        rom_args$correlation <- scores$correlation
      }
      rom <- tryCatch(
        withr::with_seed(seed, do.call(limma::romer, rom_args)),
        error = function(e) e
      )
      if (inherits(rom, "error")) {
        tps_abort(paste0("romer() failed for contrast '", cn, "': ",
                         conditionMessage(rom)))
      }
      p_up[, j] <- .match_set_stat(rom, terms, "Up")
      p_down[, j] <- .match_set_stat(rom, terms, "Down")
      p_mixed[, j] <- .match_set_stat(rom, terms, "Mixed")
      p_directional[, j] <- pmin(1, 2 * pmin(p_up[, j], p_down[, j]))
    } else if (identical(engine_used, "fry")) {
      fr <- tryCatch(
        limma::fry(y, index = index, design = design, contrast = cj,
                   sort = "none"),
        error = function(e) e
      )
      if (inherits(fr, "error")) {
        tps_abort(paste0("fry() failed for contrast '", cn, "': ",
                         conditionMessage(fr)))
      }
      p_directional[, j] <- .match_set_stat(fr, terms, "PValue")
      p_mixed[, j] <- .match_set_stat(fr, terms, "PValue.Mixed")
      test_hint[, j] <- .match_set_chr(fr, terms, "Direction")
    }
  }

  q_directional <- .adjust_p_matrix(p_directional, correction)
  q_mixed <- .adjust_p_matrix(p_mixed, correction)
  test_direction <- .derive_test_direction(p_up, p_down, test_hint)
  net_direction <- .net_direction_matrix(NET)
  null_type <- switch(engine_used,
                      camera = "camera_competitive",
                      romer = "romer_rotation",
                      fry = "fry_rotation",
                      engine_used)

  .finalize_calibrated(scores, arg_capture, list(
    p_directional = p_directional,
    p_mixed = p_mixed,
    q_directional = q_directional,
    q_mixed = q_mixed,
    p_up = p_up,
    p_down = p_down,
    p_net = p_directional,
    p_rem = p_mixed,
    q_net = q_directional,
    q_rem = q_mixed,
    test_direction = test_direction,
    net_direction = net_direction,
    direction_concordant = .direction_concordant(test_direction, net_direction),
    engine_used = engine_used,
    null_type = null_type,
    nrot = as.integer(nrot),
    seed = seed,
    correction = correction
  ))
}

## ---- tps_calibrate_pathways ------------------------------------------------
#' Calibrate pathway scores with CAMERA, ROMER, or a diagnostic permutation
#'
#' Applies a gene-set significance engine to a \code{tps_pathway_scores}
#' object. The default \code{engine = "auto"} uses \code{\link[limma]{camera}}
#' when no repeated-measures block is present and \code{\link[limma]{romer}}
#' when \code{scores$block_value} is not \code{NULL} (ADR 0011). Explicit
#' \code{engine = "camera"} aborts if a block is present. \code{engine = "romer"}
#' is allowed with or without a block. \code{engine = "fry"} is a sensitivity
#' engine. \code{engine = "permutation"} is the ADR 0010 diagnostic
#' permutation+VIF path.
#'
#' NET, REM, and COH are summaries from \code{\link{tps_pathway_scores}} and
#' are not recomputed here. Significance columns are \code{p_directional} /
#' \code{q_directional} (net direction) and \code{p_mixed} / \code{q_mixed}
#' (mixed remodeling). COH has no p-value.
#'
#' CAMERA is called with \code{sort = FALSE}, \code{inter.gene.cor = NA}, and
#' \code{trend.var} from the fit. Directional tests are parametric
#' (\code{use.ranks = FALSE}); mixed tests use ranks. This function never
#' calls \code{cameraPR}. ROMER uses \code{nrot} rotations with a local RNG
#' stream (\code{seed}).
#'
#' Auto never switches engine because of missing values: the feature universe
#' is formed once in \code{\link{tps_pathway_scores}}.
#'
#' The permutation engine remains available for diagnostics. Default
#' \code{mode} is \code{"abundance_stratified"}: proteins are binned into
#' \code{n_bins} equal-frequency abundance quantiles and permuted within each
#' bin. \code{intergene = "vif"} estimates each pathway's mean pairwise
#' residual correlation and scales the observed NET / REM excess by
#' \code{1/sqrt(VIF)} before the permutation comparison, with
#' \code{VIF = 1 + (k - 1) * rho}. Do not retune \code{vif_cap} for CAMERA
#' or ROMER; it applies only to the permutation engine.
#'
#' @param scores a \code{tps_pathway_scores} object from
#'   \code{\link{tps_pathway_scores}}.
#' @param abundance numeric vector of per-feature abundance (e.g. AveExpr from
#'   limma), aligned to \code{scores$feature_ids}. A named vector is reordered
#'   when its names match either feature ids or mapped pathway ids. Required for
#'   permutation \code{abundance_stratified} mode; ignored otherwise.
#' @param n_perm integer; number of permutations for
#'   \code{engine = "permutation"} (default 1000; use 5000+ for publication).
#' @param mode character; \code{"abundance_stratified"} (default) or
#'   \code{"unstratified"}. Used only by the permutation engine.
#' @param n_bins integer; abundance bins for stratified permutation (default 4).
#' @param seed integer; RNG seed (default 1). Uses a local RNG stream so the
#'   global state is not perturbed. Used by ROMER and the permutation engine.
#' @param correction character; \code{"per_contrast"} (default) applies BH
#'   within each contrast, \code{"global"} applies BH across all
#'   pathway x contrast tests (more conservative for multi-timepoint designs).
#' @param intergene character; \code{"vif"} (default) applies the inter-gene
#'   variance inflation, or \code{"none"} for the unadjusted permutation.
#'   Used only by the permutation engine.
#' @param vif_cap numeric; upper bound on the VIF for the permutation engine
#'   (default \code{3}). See ADR 0010. Not used by CAMERA/ROMER/FRY.
#' @param engine character; \code{"auto"} (default), \code{"camera"},
#'   \code{"romer"}, \code{"fry"}, or \code{"permutation"}.
#' @param nrot integer; rotations for ROMER (default 9999).
#' @return a list of class \code{tps_pathway_calibrated} extending the input
#'   with \code{p_directional}, \code{p_mixed}, \code{q_directional},
#'   \code{q_mixed}, \code{engine_used}, and \code{null_type}. Permutation
#'   also stores \code{p_net}/\code{p_rem} (aliased to the new names) and
#'   null summaries. CAMERA/ROMER store the same aliases so downstream
#'   code can read either naming scheme.
#' @export
#' @examples
#' data(tps_example, package = "timeProtState")
#' se <- tps_import(tps_example_abundance, tps_example_metadata,
#'                  tps_example_annotation)
#' design <- tps_make_design(se)
#' contrasts <- tps_make_contrasts(design, endpoint = 7)
#' fit <- tps_fit(se, design, contrasts)
#' pws <- tps_pathway_scores(fit, tps_example_term2gene,
#'                            id_col = "locus_tag", family = "vs_baseline")
#' pws_c <- tps_calibrate_pathways(pws, nrot = 99)
tps_calibrate_pathways <- function(scores, abundance = NULL, n_perm = 1000L,
                                    mode = c("abundance_stratified","unstratified"),
                                    n_bins = 4L, seed = 1L,
                                    correction = c("per_contrast","global"),
                                    intergene = c("vif", "none"),
                                    vif_cap = 3,
                                    engine = c("auto", "camera", "romer",
                                               "fry", "permutation"),
                                    nrot = 9999L) {
  engine <- match.arg(engine)
  correction <- match.arg(correction)
  if (!is.numeric(vif_cap) || length(vif_cap) != 1L || is.na(vif_cap) ||
      vif_cap < 1) {
    tps_abort("vif_cap must be one numeric value >= 1 (use Inf for no cap).")
  }
  if (!is.numeric(nrot) || length(nrot) != 1L || is.na(nrot) ||
      !is.finite(nrot) || nrot < 1L || nrot != as.integer(nrot)) {
    tps_abort("nrot must be one positive integer.")
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed)) {
    tps_abort("seed must be one finite value.")
  }
  arg_capture <- list(call = "tps_calibrate_pathways", n_perm = n_perm,
                       mode = mode, n_bins = n_bins, seed = seed,
                       correction = correction, intergene = intergene,
                       vif_cap = vif_cap, engine = engine, nrot = nrot,
                       timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  if (!inherits(scores, "tps_pathway_scores")) {
    tps_abort("scores must be a tps_pathway_scores object.")
  }

  ## Missingness must not change the engine: auto keys only off block_value.
  engine_used <- if (identical(engine, "auto")) {
    if (!is.null(scores$block_value)) "romer" else "camera"
  } else {
    engine
  }
  if (identical(engine_used, "camera") && !is.null(scores$block_value)) {
    tps_abort(paste0(
      "engine='camera' cannot be used when a repeated-measures block is ",
      "present. Use engine='romer' (the auto default when block is set) ",
      "or engine='permutation'."
    ))
  }
  if (identical(engine_used, "fry") && !is.null(scores$block_value)) {
    tps_abort(paste0(
      "engine='fry' cannot be used when a repeated-measures block is ",
      "present. Use engine='romer'."
    ))
  }
  if (engine_used %in% c("camera", "romer", "fry")) {
    return(.calibrate_competitive(scores, engine_used, correction, nrot, seed,
                                  arg_capture))
  }

  mode <- match.arg(mode)
  intergene <- match.arg(intergene)
  arg_capture$mode <- mode
  arg_capture$intergene <- intergene
  eps <- 1e-12
  M <- scores$membership
  stat_mat <- scores$statistic_matrix
  sizes <- scores$pathway_table$observed_size
  n_prot <- nrow(stat_mat)
  n_path <- nrow(M)
  n_contr <- ncol(stat_mat)
  if (!is.numeric(n_perm) || length(n_perm) != 1L || is.na(n_perm) ||
      !is.finite(n_perm) || n_perm < 1L ||
      n_perm != as.integer(n_perm)) {
    tps_abort("n_perm must be one positive integer.")
  }
  if (!is.numeric(n_bins) || length(n_bins) != 1L || is.na(n_bins) ||
      !is.finite(n_bins) || n_bins < 1L ||
      n_bins != as.integer(n_bins)) {
    tps_abort("n_bins must be one positive integer.")
  }

  ## abundance alignment for stratified mode
  bin_list <- NULL; bin_template <- NULL
  ave_bin <- NULL
  if (mode == "abundance_stratified") {
    ## T03/M2: if the user did not pass abundance but the scores object
    ## carries AveExpr from the fit, use it automatically. This makes
    ## score -> calibration a linear path that does not require the user to
    ## manually rebuild a hidden field.
    if (is.null(abundance) && !is.null(scores$abundance) &&
        is.numeric(scores$abundance)) {
      abundance <- scores$abundance
    }
    if (is.null(abundance)) {
      tps_abort(paste0(
        "mode='abundance_stratified' requires an 'abundance' vector. ",
        "Either pass one explicitly, or construct the scores object from ",
        "a tps_fit (which carries AveExpr automatically) via tps_pathway_scores."
      ))
    }
    if (length(abundance) != n_prot) {
      tps_abort("abundance length (%d) != universe size (%d).",
                length(abundance), n_prot)
    }
    if (!is.numeric(abundance) || any(!is.finite(abundance))) {
      tps_abort("abundance must be numeric and finite for every feature.")
    }
    if (!is.null(names(abundance))) {
      if (anyDuplicated(names(abundance))) {
        tps_abort("named abundance values must have unique names.")
      }
      target <- if (all(scores$feature_ids %in% names(abundance))) {
        scores$feature_ids
      } else if (all(scores$universe %in% names(abundance))) {
        scores$universe
      } else {
        NULL
      }
      if (is.null(target)) {
        tps_abort("abundance names do not match score feature or pathway ids.")
      }
      abundance <- abundance[target]
    }
    ave_bin <- .make_abundance_bins(abundance, n_bins)
    bin_list <- split(seq_len(n_prot), ave_bin)
    bin_template <- unlist(bin_list, use.names = FALSE)
  }

  obs_net <- scores$scores$NET
  obs_rem <- scores$scores$REM
  Mt <- as(M, "CsparseMatrix")
  abs_stat <- abs(stat_mat)

  ## accumulators
  rem_null_ge <- net_null_sum <- net_null_sqsum <- rem_null_sum <-
    rem_null_sqsum <-
    matrix(0.0, n_path, n_contr)

  withr::with_seed(seed, {
    for (b in seq_len(n_perm)) {
      idx <- .perm_index(n_prot, mode, bin_list, bin_template)
      Tp <- stat_mat[idx, , drop = FALSE]
      net_b <- as.matrix(Mt %*% Tp) / sizes
      rem_b <- as.matrix(Mt %*% abs_stat[idx, , drop = FALSE]) / sizes
      rem_null_ge     <- rem_null_ge + (rem_b >= obs_rem - eps)
      net_null_sum    <- net_null_sum + net_b
      net_null_sqsum  <- net_null_sqsum + net_b^2
      rem_null_sum    <- rem_null_sum + rem_b
      rem_null_sqsum  <- rem_null_sqsum + rem_b^2
    }
  })

  net_null_mean <- net_null_sum / n_perm
  net_null_sd <- sqrt(pmax(net_null_sqsum / n_perm - net_null_mean^2, 0))
  rem_null_mean <- rem_null_sum / n_perm
  rem_null_sd <- sqrt(pmax(rem_null_sqsum / n_perm - rem_null_mean^2, 0))
  dimnames(net_null_mean) <- dimnames(obs_net)
  dimnames(net_null_sd) <- dimnames(obs_net)
  dimnames(rem_null_mean) <- dimnames(obs_rem)
  dimnames(rem_null_sd) <- dimnames(obs_rem)

  ## Inter-gene VIF: label permutation destroys member correlation, so the
  ## permutation scale is too small when rho > 0. Scale the observed
  ## departure (ADR 0010). intergene='none' keeps the ADR 0005 p-values.
  ## vif_cap limits the VIF to prevent over-correction in large pathways
  ## with high estimated rho (the compound-symmetric VIF grows linearly
  ## with k, which can over-penalize real signal).
  intergene_cor <- .estimate_pathway_rho(scores$residuals, M,
                                         n_floor_draws = 50L, seed = seed)
  vif <- .pathway_vif(intergene_cor, sizes)
  if (!identical(intergene, "vif")) {
    intergene_cor[] <- 0
    vif[] <- 1
  }
  vif <- pmin(vif, vif_cap)
  vif_mat <- matrix(vif, nrow = n_path, ncol = n_contr)
  vif_scale <- sqrt(pmax(vif_mat, 1))

  ## A two-sided permutation test must be centered on the permutation-null
  ## mean. Comparing |NET| directly with |NET_null| is biased when the whole
  ## protein universe has a directional shift. Replay the deterministic
  ## permutation stream now that the null centers are known.
  net_null_ge <- matrix(0.0, n_path, n_contr)
  observed_net_departure <- abs(obs_net - net_null_mean) / vif_scale
  ## REM must also use a centered second-pass comparison, consistent with
  ## NET. Previously (ADR 0005), REM retained the first-pass uncentered
  ## comparison, which was inconsistent with the centered NET test and
  ## biased when the protein universe has a directional shift.
  rem_null_ge <- matrix(0.0, n_path, n_contr)
  observed_rem_excess <- pmax(obs_rem - rem_null_mean, 0) / vif_scale
  withr::with_seed(seed, {
    for (b in seq_len(n_perm)) {
      idx <- .perm_index(n_prot, mode, bin_list, bin_template)
      Tp <- stat_mat[idx, , drop = FALSE]
      net_b <- as.matrix(Mt %*% Tp) / sizes
      rem_b <- as.matrix(Mt %*% abs_stat[idx, , drop = FALSE]) / sizes
      net_null_ge <- net_null_ge +
        (abs(net_b - net_null_mean) >= observed_net_departure - eps)
      rem_null_ge <- rem_null_ge +
        (pmax(rem_b - rem_null_mean, 0) >= observed_rem_excess - eps)
    }
  })

  p_net <- (net_null_ge + 1) / (n_perm + 1)
  p_rem <- (rem_null_ge + 1) / (n_perm + 1)
  dimnames(p_net) <- dimnames(obs_net)
  dimnames(p_rem) <- dimnames(obs_rem)
  q_net <- .adjust_p_matrix(p_net, correction)
  q_rem <- .adjust_p_matrix(p_rem, correction)
  p_up <- matrix(NA_real_, nrow(obs_net), ncol(obs_net),
                 dimnames = dimnames(obs_net))
  p_down <- matrix(NA_real_, nrow(obs_net), ncol(obs_net),
                   dimnames = dimnames(obs_net))
  test_direction <- .derive_test_direction(p_up, p_down)
  net_direction <- .net_direction_matrix(obs_net)

  .finalize_calibrated(scores, arg_capture, list(
    p_net = p_net, p_rem = p_rem,
    q_net = q_net, q_rem = q_rem,
    p_directional = p_net, p_mixed = p_rem,
    q_directional = q_net, q_mixed = q_rem,
    p_up = p_up, p_down = p_down,
    test_direction = test_direction,
    net_direction = net_direction,
    direction_concordant = .direction_concordant(test_direction, net_direction),
    net_null_mean = net_null_mean, net_null_sd = net_null_sd,
    rem_null_mean = rem_null_mean, rem_null_sd = rem_null_sd,
    n_perm = n_perm, mode = mode, n_bins = n_bins, seed = seed,
    abundance_bins = if (mode == "abundance_stratified") ave_bin else NULL,
    correction = correction,
    intergene = intergene,
    intergene_cor = setNames(intergene_cor, rownames(obs_net)),
    vif = setNames(vif, rownames(obs_net)),
    vif_cap = vif_cap,
    engine_used = "permutation",
    null_type = "permutation",
    nrot = as.integer(nrot)
  ))
}

## ---- classification --------------------------------------------------------
.classify_one_pathway <- function(net, rem, coh, q_directional, q_mixed,
                                   alpha = 0.05, coh_cut = 0.5) {
  if (isTRUE(q_mixed >= alpha) && isTRUE(q_directional >= alpha)) {
    return("no_strong_state")
  }
  if (isTRUE(q_directional < alpha) && isTRUE(net > 0) && isTRUE(coh >= coh_cut)) {
    return("directional_up")
  }
  if (isTRUE(q_directional < alpha) && isTRUE(net < 0) && isTRUE(coh >= coh_cut)) {
    return("directional_down")
  }
  if (isTRUE(q_mixed < alpha) &&
      (isTRUE(q_directional >= alpha) || isTRUE(coh < coh_cut))) {
    return("mixed_remodeling")
  }
  if (isTRUE(q_directional < alpha) && isTRUE(coh < coh_cut)) {
    return("directional_with_internal_opposition")
  }
  "no_strong_state"
}

#' Classify pathways into directional / remodeling states
#'
#' Assigns each pathway x contrast a state label from a calibrated
#' \code{tps_pathway_calibrated} object: \code{directional_up},
#' \code{directional_down}, \code{mixed_remodeling},
#' \code{directional_with_internal_opposition}, or \code{no_strong_state}.
#' Labels use CAMERA/ROMER (or permutation) q-values
#' \code{q_directional} and \code{q_mixed}.
#'
#' @param calibrated a \code{tps_pathway_calibrated} object from
#'   \code{\link{tps_calibrate_pathways}}.
#' @param alpha numeric; q-value threshold (default 0.05).
#' @param coh_cut numeric; COH threshold separating directional from
#'   remodeling (default 0.5).
#' @return a data.frame (long format) with one row per pathway x contrast:
#'   term, contrast, size, NET, REM, COH, q_directional, q_mixed, state.
#'   \code{q_net} and \code{q_rem} are aliases of \code{q_directional} and
#'   \code{q_mixed} for one release.
#'
#' @examples
#' data(tps_example, package = "timeProtState")
#' se <- tps_import(tps_example_abundance, tps_example_metadata,
#'                  tps_example_annotation)
#' design <- tps_make_design(se)
#' contrasts <- tps_make_contrasts(design, endpoint = 7)
#' fit <- tps_fit(se, design, contrasts)
#' pws <- tps_pathway_scores(fit, tps_example_term2gene,
#'                            id_col = "locus_tag", family = "vs_baseline")
#' pws_c <- tps_calibrate_pathways(pws, nrot = 99)
#' tps_classify_pathways(pws_c)
#' @export
tps_classify_pathways <- function(calibrated, alpha = 0.05, coh_cut = 0.5) {
  if (!inherits(calibrated, "tps_pathway_calibrated")) {
    tps_abort(paste0("calibrated must be a tps_pathway_calibrated object; run ",
                     "tps_calibrate_pathways first."))
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    tps_abort("alpha must be in (0,1).")
  }
  if (!is.numeric(coh_cut) || length(coh_cut) != 1L || is.na(coh_cut) ||
      !is.finite(coh_cut) || coh_cut < 0 || coh_cut > 1) {
    tps_abort("coh_cut must be in [0,1].")
  }
  NET <- calibrated$scores$NET; REM <- calibrated$scores$REM
  COH <- calibrated$scores$COH
  q_directional <- calibrated$q_directional
  q_mixed <- calibrated$q_mixed
  if (is.null(q_directional)) q_directional <- calibrated$q_net
  if (is.null(q_mixed)) q_mixed <- calibrated$q_rem
  if (is.null(q_directional) || is.null(q_mixed)) {
    tps_abort("calibrated object is missing q_directional/q_mixed.")
  }
  sizes <- calibrated$pathway_table$observed_size
  terms <- rownames(NET); contrasts <- colnames(NET)
  rows <- list()
  for (j in seq_along(contrasts)) {
    for (i in seq_along(terms)) {
      rows[[length(rows) + 1L]] <- data.frame(
        term = terms[i], contrast = contrasts[j], size = sizes[i],
        NET = NET[i, j], REM = REM[i, j], COH = COH[i, j],
        q_directional = q_directional[i, j], q_mixed = q_mixed[i, j],
        q_net = q_directional[i, j], q_rem = q_mixed[i, j],
        state = .classify_one_pathway(NET[i, j], REM[i, j], COH[i, j],
                                       q_directional[i, j], q_mixed[i, j],
                                       alpha = alpha, coh_cut = coh_cut),
        stringsAsFactors = FALSE)
    }
  }
  result <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame(
    term = character(0), contrast = character(0), size = integer(0),
    NET = numeric(0), REM = numeric(0), COH = numeric(0),
    q_directional = numeric(0), q_mixed = numeric(0),
    q_net = numeric(0), q_rem = numeric(0),
    state = character(0), stringsAsFactors = FALSE)
  attr(result, "provenance") <- list(
    stage = "classify_pathways", alpha = alpha, coh_cut = coh_cut,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    package_version = as.character(getNamespaceVersion("timeProtState")))
  result
}

## ---- methods ---------------------------------------------------------------
#' @export
print.tps_pathway_calibrated <- function(x, ...) {
  cat("<timeProtState pathway calibrated>\n")
  engine <- if (!is.null(x$engine_used)) x$engine_used else "unspecified"
  cat(sprintf("  pathways: %d   contrasts: %d   engine: %s\n",
              nrow(x$scores$NET), ncol(x$scores$NET), engine))
  if (identical(engine, "permutation")) {
    n_perm <- if (!is.null(x$n_perm)) x$n_perm else NA_integer_
    mode <- if (!is.null(x$mode)) x$mode else NA_character_
    cat(sprintf("  n_perm: %d   mode: %s\n", n_perm, mode))
    if (identical(x$mode, "abundance_stratified")) {
      cat(sprintf("  abundance bins: %d\n", x$n_bins))
    }
    ig <- if (!is.null(x$intergene)) x$intergene else "none"
    if (identical(ig, "vif") && !is.null(x$vif)) {
      cat(sprintf("  intergene: vif   median VIF: %.2f   median rho: %.2f\n",
                  stats::median(x$vif), stats::median(x$intergene_cor)))
    } else {
      cat(sprintf("  intergene: %s\n", ig))
    }
  } else if (identical(engine, "romer") && !is.null(x$nrot)) {
    cat(sprintf("  nrot: %d   seed: %s\n", x$nrot, x$seed))
  }
  invisible(x)
}
