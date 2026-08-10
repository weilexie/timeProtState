## ============================================================================
## Pathway calibration: permutation null + classification
## ============================================================================
##
## Implements STATISTICAL_CONTRACT.md section 8.4-8.5 and DECISIONS D5.
## Applies a size-matched permutation null (default: abundance-stratified)
## to compute empirical p_net / p_rem, BH q-values, and pathway state labels.
##
## Uses (+1) smoothing, per-contrast BH by default, and four-bin abundance
## stratification with one shared permutation index across pathways.

## ---- stratification helpers ------------------------------------------------
.make_abundance_bins <- function(abundance, n_bins) {
  ## equal-frequency bins over abundance (AveExpr); returns integer bin per protein
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

## ---- tps_calibrate_pathways ------------------------------------------------
#' Calibrate pathway scores with a permutation null
#'
#' Applies a size-matched permutation null to a \code{tps_pathway_scores}
#' object and computes empirical p-values for NET (two-sided around its
#' permutation-null mean) and REM (upper-tail), Benjamini-Hochberg q-values,
#' and null mean/SD per pathway.
#'
#' The default mode is \code{"abundance_stratified"}: proteins are binned into
#' \code{n_bins} equal-frequency abundance quantiles and permuted within each
#' bin. One shared permutation index retains overlap among pathway definitions;
#' stratification reduces differences in abundance distributions but does not
#' preserve within-pathway protein covariance.
#'
#' @param scores a \code{tps_pathway_scores} object from
#'   \code{\link{tps_pathway_scores}}.
#' @param abundance numeric vector of per-feature abundance (e.g. AveExpr from
#'   limma), aligned to \code{scores$feature_ids}. A named vector is reordered
#'   when its names match either feature ids or mapped pathway ids. Required for
#'   \code{abundance_stratified} mode; ignored for \code{unstratified}.
#' @param n_perm integer; number of permutations (default 1000; use 5000+ for
#'   publication).
#' @param mode character; \code{"abundance_stratified"} (default) or
#'   \code{"unstratified"}.
#' @param n_bins integer; abundance bins for stratified mode (default 4).
#' @param seed integer; RNG seed (default 1). Uses a local RNG stream so the
#'   global state is not perturbed.
#' @param correction character; \code{"per_contrast"} (default) applies BH
#'   within each contrast, \code{"global"} applies BH across all
#'   pathway x contrast tests (more conservative for multi-timepoint designs).
#' @return a list of class \code{tps_pathway_calibrated} extending the input
#'   with \code{p_net}, \code{p_rem}, \code{q_net}, \code{q_rem}, null summaries.
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
#' ave <- rep(1, nrow(pws$statistic_matrix))
#' pws_c <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 10, seed = 1)
tps_calibrate_pathways <- function(scores, abundance = NULL, n_perm = 1000L,
                                    mode = c("abundance_stratified","unstratified"),
                                    n_bins = 4L, seed = 1L,
                                    correction = c("per_contrast","global")) {
  mode <- match.arg(mode)
  correction <- match.arg(correction)
  arg_capture <- list(call = "tps_calibrate_pathways", n_perm = n_perm,
                       mode = mode, n_bins = n_bins, seed = seed,
                       correction = correction,
                       timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  if (!inherits(scores, "tps_pathway_scores")) {
    tps_abort("scores must be a tps_pathway_scores object.")
  }
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
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed)) {
    tps_abort("seed must be one finite value.")
  }

  ## abundance alignment for stratified mode
  bin_list <- NULL; bin_template <- NULL
  if (mode == "abundance_stratified") {
    if (is.null(abundance)) {
      tps_abort("mode='abundance_stratified' requires an 'abundance' vector.")
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

  ## A two-sided permutation test must be centered on the permutation-null
  ## mean. Comparing |NET| directly with |NET_null| is biased when the whole
  ## protein universe has a directional shift. Replay the deterministic
  ## permutation stream now that the null centers are known.
  net_null_ge <- matrix(0.0, n_path, n_contr)
  observed_net_departure <- abs(obs_net - net_null_mean)
  withr::with_seed(seed, {
    for (b in seq_len(n_perm)) {
      idx <- .perm_index(n_prot, mode, bin_list, bin_template)
      net_b <- as.matrix(Mt %*% stat_mat[idx, , drop = FALSE]) / sizes
      net_null_ge <- net_null_ge +
        (abs(net_b - net_null_mean) >= observed_net_departure - eps)
    }
  })

  p_net <- (net_null_ge + 1) / (n_perm + 1)
  p_rem <- (rem_null_ge + 1) / (n_perm + 1)
  dimnames(p_net) <- dimnames(obs_net)
  dimnames(p_rem) <- dimnames(obs_rem)
  q_net <- .adjust_p_matrix(p_net, correction)
  q_rem <- .adjust_p_matrix(p_rem, correction)

  ## append calibrate provenance as a new element in the provenance list.
  new_prov <- .make_provenance(arg_capture, scores, stage = "calibrate_pathways")
  old_prov <- if (is.list(scores$provenance)) scores$provenance else list()
  n_old <- length(old_prov)
  combined_prov <- old_prov
  combined_prov[[n_old + 1L]] <- new_prov

  ## build output by copying scores then overwriting/adding fields
  ## (c() would append rather than overwrite same-named elements)
  out <- scores
  out$p_net <- p_net; out$p_rem <- p_rem
  out$q_net <- q_net; out$q_rem <- q_rem
  out$net_null_mean <- net_null_mean; out$net_null_sd <- net_null_sd
  out$rem_null_mean <- rem_null_mean; out$rem_null_sd <- rem_null_sd
  out$n_perm <- n_perm; out$mode <- mode; out$n_bins <- n_bins; out$seed <- seed
  out$abundance_bins <- if (mode == "abundance_stratified") ave_bin else NULL
  out$provenance <- combined_prov
  out$correction <- correction
  class(out) <- c("tps_pathway_calibrated", class(scores))
  out
}

## ---- classification --------------------------------------------------------
.classify_one_pathway <- function(net, rem, coh, q_net, q_rem,
                                   alpha = 0.05, coh_cut = 0.5) {
  if (isTRUE(q_rem >= alpha) && isTRUE(q_net >= alpha)) return("no_strong_state")
  if (isTRUE(q_net < alpha) && isTRUE(net > 0) && isTRUE(coh >= coh_cut)) {
    return("directional_up")
  }
  if (isTRUE(q_net < alpha) && isTRUE(net < 0) && isTRUE(coh >= coh_cut)) {
    return("directional_down")
  }
  if (isTRUE(q_rem < alpha) &&
      (isTRUE(q_net >= alpha) || isTRUE(coh < coh_cut))) {
    return("mixed_remodeling")
  }
  if (isTRUE(q_net < alpha) && isTRUE(coh < coh_cut)) {
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
#'
#' @param calibrated a \code{tps_pathway_calibrated} object from
#'   \code{\link{tps_calibrate_pathways}}.
#' @param alpha numeric; q-value threshold (default 0.05).
#' @param coh_cut numeric; COH threshold separating directional from
#'   remodeling (default 0.5).
#' @return a data.frame (long format) with one row per pathway x contrast:
#'   term, contrast, size, NET, REM, COH, q_net, q_rem, state.
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
#' ave <- rep(1, nrow(pws$statistic_matrix))
#' pws_c <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 10, seed = 1)
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
  q_net <- calibrated$q_net; q_rem <- calibrated$q_rem
  sizes <- calibrated$pathway_table$observed_size
  terms <- rownames(NET); contrasts <- colnames(NET)
  rows <- list()
  for (j in seq_along(contrasts)) {
    for (i in seq_along(terms)) {
      rows[[length(rows) + 1L]] <- data.frame(
        term = terms[i], contrast = contrasts[j], size = sizes[i],
        NET = NET[i, j], REM = REM[i, j], COH = COH[i, j],
        q_net = q_net[i, j], q_rem = q_rem[i, j],
        state = .classify_one_pathway(NET[i, j], REM[i, j], COH[i, j],
                                       q_net[i, j], q_rem[i, j],
                                       alpha = alpha, coh_cut = coh_cut),
        stringsAsFactors = FALSE)
    }
  }
  result <- do.call(rbind, rows)
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
  cat(sprintf("  pathways: %d   contrasts: %d   n_perm: %d   mode: %s\n",
              nrow(x$scores$NET), ncol(x$scores$NET), x$n_perm, x$mode))
  if (x$mode == "abundance_stratified") {
    cat(sprintf("  abundance bins: %d\n", x$n_bins))
  }
  invisible(x)
}
