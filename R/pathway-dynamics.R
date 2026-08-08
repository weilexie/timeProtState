## ============================================================================
## Pathway dynamics: contrast-driven remodeling trajectories
## ============================================================================
##
## Implements the core methodological upgrade: from endpoint pathway snapshots
## to contrast-driven pathway remodeling DYNAMICS across the time course.
##
## The insight (confirmed by code audit): tps_pathway_scores and
## tps_calibrate_pathways are already contrast-agnostic. Passing
## family="vs_baseline" yields NET/REM/COH for EVERY (condition, timepoint)
## pair. This module classifies those trajectories into biologically
## interpretable dynamic patterns and links them to protein-level states.
##
## This is the differentiator vs TiDEomics: TiDEomics does protein-level DE
## then standard enrichment (two-step, loses pathway direction structure).
## We use the full rank-audited contrast family's moderated statistics to
## drive pathway-level inference directly -- one integrated step that
## preserves direction and uncertainty, and is calibrated per-timepoint.

## ---- contrast-time parser --------------------------------------------------
## Parse contrast names like "Acid_t7_vs_baseline" into (condition, time).
## Returns a data.frame with columns: contrast, condition, time.
.parse_contrast_time <- function(contrast_names) {
  ## strip "_vs_baseline" suffix
  base <- sub("_vs_baseline$", "", contrast_names)
  ## split on "_t" -- condition may contain underscores, time is last segment
  ## pattern: <condition>_t<time>  e.g. "Acid_t7", "Starvation_t3"
  cond <- sub("_t[0-9.]+$", "", base)
  time_str <- sub("^.*_t", "", base)
  time <- suppressWarnings(as.numeric(time_str))
  data.frame(contrast = contrast_names, condition = cond, time = time,
             stringsAsFactors = FALSE)
}

## ---- classify one pathway-condition trajectory -----------------------------
## Given NET/REM/COH/q_net vectors ordered by time, classify the dynamic pattern.
.classify_one_trajectory <- function(net, rem, coh, q_net, q_rem, alpha, coh_cut) {
  n <- length(net)
  if (n == 0L) return("no_dynamic_signal")

  sig <- !is.na(q_net) & q_net < alpha
  sig_rem <- !is.na(q_rem) & q_rem < alpha

  if (!any(sig) && !any(sig_rem)) return("no_dynamic_signal")

  ## NA-safe helpers: comparisons that may involve NA use isTRUE()
  rem_n <- rem[n]; rem_1 <- rem[1L]; coh_1 <- coh[1L]

  ## Priority: specific patterns first (transient, gradual, delayed),
  ## then the general "early_directional" catch-all.

  ## transient_burst: REM peaks significantly in the middle, drops at the end.
  if (n >= 3L) {
    rem_mid <- rem[seq(2L, n - 1L)]
    sig_rem_mid <- sig_rem[seq(2L, n - 1L)]
    mid_max <- max(rem_mid, na.rm = TRUE)
    if (any(sig_rem_mid) &&
        isTRUE(mid_max > rem_n * 1.5) &&
        isTRUE(mid_max > rem_1 * 1.5)) {
      return("transient_burst")
    }
  }

  ## delayed_opposition: only the LAST timepoint is NET-significant
  if (n >= 2L && isTRUE(sig[n]) && !any(sig[-n], na.rm = TRUE)) {
    return("delayed_opposition")
  }

  ## gradual_remodeling: REM increases over time, signal builds progressively.
  if (n >= 2L && isTRUE(rem_n > rem_1) && isTRUE(sig_rem[n]) &&
      isTRUE(rem_n > rem_1 * 1.2)) {
    return("gradual_remodeling")
  }

  ## early_directional: first timepoint NET-significant with high COH AND
  ## signal persists (at least 2 timepoints significant).
  if (isTRUE(sig[1L]) && isTRUE(coh_1 >= coh_cut) && sum(sig, na.rm = TRUE) >= 2L) {
    return("early_directional")
  }

  ## any significant NET -- directional but not early
  if (any(sig, na.rm = TRUE)) return("directional_delayed")

  ## significant REM but never directional NET -- pure remodeling
  if (any(sig_rem, na.rm = TRUE) && !any(sig, na.rm = TRUE)) return("sustained_remodeling")

  "no_dynamic_signal"
}

## ---- detect cross-condition divergence -------------------------------------
## For a pathway, check if conditions diverge in NET direction over time.
## Returns the COUNT of divergent timepoints (0 = no divergence).
.count_divergent_timepoints <- function(traj_long, alpha) {
  conditions <- unique(traj_long$condition)
  if (length(conditions) < 2L) return(0L)
  n_div <- 0L
  for (tt in unique(traj_long$time)) {
    if (is.na(tt)) next
    sub <- traj_long[traj_long$time == tt & !is.na(traj_long$time), , drop = FALSE]
    sig <- sub[!is.na(sub$q_net) & sub$q_net < alpha, , drop = FALSE]
    if (nrow(sig) >= 2L) {
      signs <- sign(sig$NET)
      if (any(signs > 0) && any(signs < 0)) n_div <- n_div + 1L
    }
  }
  n_div
}

## boolean wrapper for backward compatibility in tps_classify_pathway_dynamics
.detect_divergence <- function(traj_long, alpha) {
  .count_divergent_timepoints(traj_long, alpha) > 0L
}

## ---- tps_classify_pathway_dynamics -----------------------------------------
#' Classify pathway remodeling dynamics across the time course
#'
#' Converts a calibrated pathway scores object (from
#' \code{\link{tps_calibrate_pathways}} with multiple timepoints) into
#' per-pathway-per-condition dynamic state labels describing HOW the pathway
#' remodels over time: \code{early_directional}, \code{gradual_remodeling},
#' \code{delayed_opposition}, \code{transient_burst},
#' \code{directional_delayed}, \code{sustained_remodeling}, or
#' \code{no_dynamic_signal}. Also detects \code{divergent} pathways where
#' conditions separate in NET direction.
#'
#' This is the time-course-aware upgrade of \code{\link{tps_classify_pathways}}
#' (which gives a single endpoint snapshot). Use with
#' \code{family = "vs_baseline"} in \code{\link{tps_pathway_scores}} to get all
#' timepoints.
#'
#' @param calibrated a \code{tps_pathway_calibrated} object from
#'   \code{\link{tps_calibrate_pathways}}.
#' @param alpha numeric; q-value threshold (default 0.05).
#' @param coh_cut numeric; COH threshold for "directional" (default 0.5).
#' @return a data.frame with columns: \code{term, condition, time, NET, REM,
#'   COH, q_net, q_rem, dynamic_state} (long format, one row per
#'   pathway-condition-time), plus a summary column \code{is_divergent} marking
#'   pathways where conditions diverge.
#'
#' @examples
#' \dontrun{
#' se <- tps_import(tps_example_abundance, tps_example_metadata,
#'                  tps_example_annotation)
#' design <- tps_make_design(se)
#' contrasts <- tps_make_contrasts(design, endpoint = 7)
#' fit <- tps_fit(se, design, contrasts)
#' pws <- tps_pathway_scores(fit, tps_example_term2gene,
#'                            id_col = "locus_tag", family = "vs_baseline")
#' ave <- as.numeric(tapply(
#'   tps_results(fit, "vs_baseline")$AveExpr,
#'   factor(tps_results(fit, "vs_baseline")$feature_id,
#'          levels = rownames(pws$statistic_matrix)), function(x) x[1]))
#' pws_c <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 200, seed = 1)
#' dyn <- tps_classify_pathway_dynamics(pws_c)
#' head(dyn)
#' }
#'
#' @export
tps_classify_pathway_dynamics <- function(calibrated, alpha = 0.05, coh_cut = 0.5) {
  if (!inherits(calibrated, "tps_pathway_calibrated")) {
    tps_abort(paste0("calibrated must be a tps_pathway_calibrated object; ",
                     "run tps_calibrate_pathways first."))
  }
  NET <- calibrated$scores$NET
  REM <- calibrated$scores$REM
  COH <- calibrated$scores$COH
  q_net <- calibrated$q_net
  q_rem <- calibrated$q_rem
  terms <- rownames(NET)
  contrasts <- colnames(NET)

  ## parse time + condition from contrast names
  meta <- .parse_contrast_time(contrasts)
  if (all(is.na(meta$time))) {
    tps_abort("Could not parse timepoints from contrast names. Expected ",
              "format '<condition>_t<time>_vs_baseline' (from family='vs_baseline').")
  }

  ## build long table: one row per term x condition x time
  rows <- list()
  for (i in seq_along(terms)) {
    for (cc in unique(meta$condition)) {
      idx <- which(meta$condition == cc)
      idx <- idx[order(meta$time[idx])]  # time-ordered
      if (length(idx) == 0L) next
      dyn <- .classify_one_trajectory(
        net = NET[i, idx], rem = REM[i, idx], coh = COH[i, idx],
        q_net = q_net[i, idx], q_rem = q_rem[i, idx],
        alpha = alpha, coh_cut = coh_cut)
      for (k in seq_along(idx)) {
        j <- idx[k]
        rows[[length(rows) + 1L]] <- data.frame(
          term = terms[i], condition = cc, time = meta$time[j],
          NET = NET[i, j], REM = REM[i, j], COH = COH[i, j],
          q_net = q_net[i, j], q_rem = q_rem[i, j],
          dynamic_state = dyn, stringsAsFactors = FALSE)
      }
    }
  }
  long <- do.call(rbind, rows)

  ## detect divergence per pathway
  divergent_terms <- character(0)
  for (tm in terms) {
    sub <- long[long$term == tm, , drop = FALSE]
    if (.detect_divergence(sub, alpha)) divergent_terms <- c(divergent_terms, tm)
  }
  long$is_divergent <- long$term %in% divergent_terms

  ## attach summary (one row per term)
  attr(long, "summary") <- .dynamics_summary(long)
  attr(long, "alpha") <- alpha
  attr(long, "coh_cut") <- coh_cut
  long
}

## ---- dynamics summary (one row per pathway) --------------------------------
.dynamics_summary <- function(long) {
  do.call(rbind, lapply(unique(long$term), function(tm) {
    sub <- long[long$term == tm, , drop = FALSE]
    data.frame(
      term = tm,
      n_conditions = length(unique(sub$condition)),
      n_timepoints = length(unique(sub$time)),
      any_divergent = any(sub$is_divergent),
      dynamic_states = paste(unique(sub$dynamic_state), collapse = ";"),
      stringsAsFactors = FALSE)
  }))
}

## ---- tps_link_states_pathways ----------------------------------------------
#' Link protein-level states with pathway-level dynamics
#'
#' Produces a consistency diagnosis linking protein evidence states (from
#' \code{\link{tps_classify_states}}) with pathway remodeling dynamics (from
#' \code{\link{tps_classify_pathway_dynamics}}). For each protein state x
#' pathway dynamic combination, reports whether proteins in that state are
#' enriched in pathways showing that dynamic pattern.
#'
#' This is not standard ORA -- it connects protein-level evidence axes to
#' pathway-level remodeling axes bidirectionally, enabling diagnoses like:
#' "do temporal-reversal proteins cluster in transient-burst pathways?"
#'
#' @param states a \code{tps_states} object.
#' @param pathway_dynamics a data.frame from \code{tps_classify_pathway_dynamics}.
#' @param term2gene pathway membership table (columns: term, gene).
#' @param protein_ids optional character vector of protein ids aligned to
#'   \code{states$state_table$feature_id}. If \code{NULL}, uses feature_id directly.
#' @return a data.frame with columns: \code{protein_state, pathway_dynamic,
#'   n_proteins, n_pathways, observed_overlap, expected_overlap, enrichment}.
#'
#' @examples
#' \dontrun{
#' dyn <- tps_classify_pathway_dynamics(pws_c)
#' link <- tps_link_states_pathways(states, dyn, tps_example_term2gene)
#' head(link)
#' }
#'
#' @export
tps_link_states_pathways <- function(states, pathway_dynamics, term2gene,
                                      protein_ids = NULL) {
  if (!inherits(states, "tps_states")) {
    tps_abort("states must be a tps_states object.")
  }
  if (!is.data.frame(pathway_dynamics) || !"dynamic_state" %in% names(pathway_dynamics)) {
    tps_abort("pathway_dynamics must be from tps_classify_pathway_dynamics.")
  }
  st <- states$state_table
  if (is.null(protein_ids)) protein_ids <- st$feature_id

  ## analysis universe: proteins observed in the state table (NOT library)
  universe <- unique(protein_ids)
  n_universe <- length(universe)

  ## protein states
  prot_states <- unique(st$state)

  ## modal dynamic state per pathway (handles multi-condition)
  dyn_agg <- stats::aggregate(dynamic_state ~ term, data = pathway_dynamics,
                        FUN = function(x) {
                          tab <- table(x)
                          names(tab)[which.max(tab)]
                        })
  names(dyn_agg)[2] <- "dynamic_state"

  rows <- list()
  for (ps in prot_states) {
    ps_proteins <- st$feature_id[st$state == ps]
    n_ps <- length(ps_proteins)
    for (dyn in unique(dyn_agg$dynamic_state)) {
      dyn_pathways <- dyn_agg$term[dyn_agg$dynamic_state == dyn]
      ## members of these pathways, restricted to analysis universe
      members <- intersect(unique(term2gene$gene[term2gene$term %in% dyn_pathways]),
                            universe)
      n_members <- length(members)
      overlap <- length(intersect(ps_proteins, members))
      n_dyn_pathways <- length(dyn_pathways)

      ## hypergeometric test: is the overlap enriched?
      ## P(X >= overlap) where X ~ Hypergeometric(N, K, n)
      ## N = universe, K = members, n = ps_proteins
      p_hyper <- if (n_members > 0 && n_ps > 0) {
        stats::phyper(overlap - 1, n_members, n_universe - n_members, n_ps,
                       lower.tail = FALSE)
      } else 1

      expected <- n_ps * n_members / max(n_universe, 1L)
      enrich <- if (expected > 0) overlap / expected else 0
      rows[[length(rows) + 1L]] <- data.frame(
        protein_state = ps, pathway_dynamic = dyn,
        n_proteins_in_state = n_ps, n_pathways_in_dynamic = n_dyn_pathways,
        observed_overlap = overlap, expected_overlap = round(expected, 2),
        enrichment = round(enrich, 3),
        p_hypergeometric = round(p_hyper, 4),
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  ## BH correction over all tests
  out$q_hypergeometric <- stats::p.adjust(out$p_hypergeometric, method = "BH")
  out
}

## ============================================================================
## Permutation-calibrated dynamic states + divergence scoring
## ============================================================================
##
## NOVELTY: TiDEomics and all competitors calibrate protein-level DE or do
## standard enrichment. NO tool calibrates the pathway-level DYNAMIC STATE
## CLASSIFICATION itself. This function generates a null distribution of
## dynamic states by permuting protein labels, then computes the empirical
## probability that a pathway-condition pair receives a non-null dynamic
## state under the null. This is the strongest methodological differentiator:
## the dynamic classification decision is itself calibrated, not just the
## underlying NET/REM statistics.

## ---- divergence score (continuous, not just boolean) -----------------------
.compute_divergence_score <- function(traj_long, alpha) {
  ## For a pathway, compute a continuous divergence score based on the
  ## maximum absolute NET difference between any pair of conditions at any
  ## timepoint, weighted by significance.
  conditions <- unique(traj_long$condition)
  if (length(conditions) < 2L) return(0)
  max_div <- 0
  for (tt in unique(traj_long$time)) {
    if (is.na(tt)) next
    sub <- traj_long[traj_long$time == tt & !is.na(traj_long$time), , drop = FALSE]
    if (nrow(sub) < 2L) next
    for (i in seq_len(nrow(sub) - 1L)) {
      for (j in seq(i + 1L, nrow(sub))) {
        ## only count if at least one is significant
        if (!is.na(sub$q_net[i]) && !is.na(sub$q_net[j])) {
          if (sub$q_net[i] < alpha || sub$q_net[j] < alpha) {
            div <- abs(sub$NET[i] - sub$NET[j])
            if (div > max_div) max_div <- div
          }
        }
      }
    }
  }
  max_div
}

#' Compute cross-condition pathway divergence scores
#'
#' Computes a continuous divergence score for each pathway: the maximum
#' absolute NET difference between any two conditions at any timepoint,
#' gated on at least one condition being significant. Also identifies
#' which condition pair and timepoint contribute the maximum divergence.
#'
#' This quantifies HOW MUCH a pathway's remodeling direction differs across
#' conditions, complementing the boolean \code{is_divergent} flag in
#' \code{\link{tps_classify_pathway_dynamics}}. The score is not calibrated
#' (no permutation null); use it as a descriptive ranking aid.
#'
#' @param calibrated a \code{tps_pathway_calibrated} object.
#' @param dyn optional pre-computed dynamics table from
#'   \code{\link{tps_classify_pathway_dynamics}}. If \code{NULL}, computed.
#' @param alpha numeric; q-value threshold (default 0.05).
#' @return a data.frame with columns: \code{term, divergence_score,
#'   n_divergent_timepoints, max_divergent_pair, max_divergent_time}.
#'
#' @examples
#' \dontrun{
#' dyn <- tps_classify_pathway_dynamics(pws_c)
#' div <- tps_divergence_score(pws_c, dyn)
#' head(div[order(-div$divergence_score), ])
#' }
#'
#' @export
tps_divergence_score <- function(calibrated, dyn = NULL, alpha = 0.05) {
  if (!inherits(calibrated, "tps_pathway_calibrated")) {
    tps_abort("calibrated must be a tps_pathway_calibrated object.")
  }
  if (is.null(dyn)) dyn <- tps_classify_pathway_dynamics(calibrated, alpha = alpha)
  terms <- unique(dyn$term)
  rows <- list()
  for (tm in terms) {
    sub <- dyn[dyn$term == tm, , drop = FALSE]
    score <- .compute_divergence_score(sub, alpha)
    ## find the pair and timepoint of max divergence
    conditions <- unique(sub$condition)
    max_pair <- NA_character_
    max_time <- NA_real_
    if (length(conditions) >= 2L && score > 0) {
      for (tt in unique(sub$time)) {
        if (is.na(tt)) next
        st <- sub[sub$time == tt & !is.na(sub$time), , drop = FALSE]
        if (nrow(st) >= 2L) {
          for (i in seq_len(nrow(st) - 1L)) {
            for (j in seq(i + 1L, nrow(st))) {
              if (abs(st$NET[i] - st$NET[j]) >= score - 1e-10) {
                max_pair <- paste(st$condition[i], st$condition[j], sep = " vs ")
                max_time <- tt
              }
            }
          }
        }
      }
    }
    n_div_tp <- .count_divergent_timepoints(sub, alpha)
    rows[[length(rows) + 1L]] <- data.frame(
      term = tm, divergence_score = round(score, 4),
      n_divergent_timepoints = n_div_tp,
      max_divergent_pair = max_pair, max_divergent_time = max_time,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

## ---- permutation-calibrated dynamic states ---------------------------------
#' Calibrate pathway dynamic states with a permutation null
#'
#' Generates a null distribution of dynamic states by permuting protein
#' labels within the calibrated pathway scores, then computes the empirical
#' probability that each pathway-condition pair receives a non-null dynamic
#' state under the null. This calibrates the CLASSIFICATION DECISION itself,
#' not just the underlying NET/REM statistics.
#'
#' This is the key methodological differentiator: no existing tool calibrates
#' pathway-level dynamic state classification. The result tells you whether
#' a pathway classified as 'early_directional' is genuinely non-null or could
#' arise by chance.
#'
#' @param calibrated a \code{tps_pathway_calibrated} object.
#' @param n_perm integer; number of permutations for the null (default 500).
#' @param alpha numeric; q-value threshold for dynamics (default 0.05).
#' @param coh_cut numeric; COH threshold (default 0.5).
#' @param seed integer; RNG seed.
#' @return a data.frame with the original dynamics plus a column
#'   \code{p_dynamic} (empirical probability of non-null state under null)
#'   and \code{q_dynamic} (BH-adjusted).
#'
#' @examples
#' \dontrun{
#' dyn <- tps_classify_pathway_dynamics(pws_c)
#' dyn_cal <- tps_calibrate_dynamics(pws_c, n_perm = 200, seed = 1)
#' head(dyn_cal[dyn_cal$q_dynamic < 0.1, ])
#' }
#'
#' @export
tps_calibrate_dynamics <- function(calibrated, n_perm = 500L, alpha = 0.05,
                                    coh_cut = 0.5, seed = 1L) {
  if (!inherits(calibrated, "tps_pathway_calibrated")) {
    tps_abort("calibrated must be a tps_pathway_calibrated object.")
  }
  if (n_perm < 1L) tps_abort("n_perm must be >= 1.")

  ## get observed dynamics
  dyn_obs <- tps_classify_pathway_dynamics(calibrated, alpha = alpha,
                                            coh_cut = coh_cut)

  ## local RNG
  rng_old <- .GlobalEnv$.Random.seed
  on.exit({
    if (is.null(rng_old)) rm(".Random.seed", envir = .GlobalEnv)
    else .GlobalEnv$.Random.seed <- rng_old
  }, add = TRUE)
  set.seed(seed)

  ## For each permutation: shuffle protein labels in the statistic matrix,
  ## recompute NET/REM/COH, re-classify, count non-null.
  stat_mat <- calibrated$statistic_matrix
  M <- calibrated$membership
  sizes <- calibrated$pathway_table$observed_size
  n_prot <- nrow(stat_mat)
  Mt <- as(M, "CsparseMatrix")
  abs_stat <- abs(stat_mat)
  eps <- 1e-12

  ## Use the parent null's mean/SD (from tps_calibrate_pathways) to construct
  ## a z-score gate for the null classification. Under permutation, NET_null
  ## is a draw from the same null distribution, so its z-score is ~N(0,1).
  ## Gate on |z| > z_alpha (two-sided) to approximate the per-timepoint
  ## significance test under the null. This is the key fix: the previous
  ## implementation used q=0.5 (never significant), making the null inert.
  net_null_mean <- calibrated$net_null_mean
  net_null_sd <- calibrated$net_null_sd
  rem_null_mean <- calibrated$rem_null_mean
  rem_null_sd <- calibrated$rem_null_sd
  ## Guard: when net_null_sd == 0 (degenerate null, e.g. singleton pathway
  ## or all-identical statistics), z-score is undefined. Set q=1 (non-sig)
  ## for those cells to avoid false positives from division by eps.
  sd_zero_net <- net_null_sd < eps
  sd_zero_rem <- rem_null_sd < eps

  ## accumulate: for each (term, condition), count #{null classified non-null}
  meta <- .parse_contrast_time(colnames(stat_mat))
  term_cond_pairs <- unique(dyn_obs[, c("term", "condition")])
  null_count <- rep(0L, nrow(term_cond_pairs))
  names(null_count) <- paste(term_cond_pairs$term, term_cond_pairs$condition, sep = "##")

  for (b in seq_len(n_perm)) {
    idx <- sample.int(n_prot)
    Tp <- stat_mat[idx, , drop = FALSE]
    net_null <- as.matrix(Mt %*% Tp) / sizes
    rem_null <- as.matrix(Mt %*% abs_stat[idx, , drop = FALSE]) / sizes
    coh_null <- abs(net_null) / (rem_null + eps)

    ## z-score gate: under null, NET_null ~ N(net_null_mean, net_null_sd)
    ## so q = 2*pnorm(-|z|) approximates the per-timepoint p-value.
    ## When net_null_sd == 0 (degenerate), set q=1 (non-significant).
    z_net <- (net_null - net_null_mean) / pmax(net_null_sd, eps)
    z_rem <- (rem_null - rem_null_mean) / pmax(rem_null_sd, eps)
    q_null <- 2 * stats::pnorm(-abs(z_net))
    q_rem_null <- 2 * stats::pnorm(-abs(z_rem))
    q_null[sd_zero_net] <- 1
    q_rem_null[sd_zero_rem] <- 1

    ## classify each term-condition under null
    terms <- rownames(net_null)
    for (i in seq_along(terms)) {
      for (cc in unique(meta$condition)) {
        idx_cc <- which(meta$condition == cc)
        idx_cc <- idx_cc[order(meta$time[idx_cc])]
        if (length(idx_cc) == 0L) next
        dyn_null <- .classify_one_trajectory(
          net = net_null[i, idx_cc], rem = rem_null[i, idx_cc],
          coh = coh_null[i, idx_cc], q_net = q_null[i, idx_cc],
          q_rem = q_rem_null[i, idx_cc], alpha = alpha, coh_cut = coh_cut)
        if (dyn_null != "no_dynamic_signal") {
          key <- paste(terms[i], cc, sep = "##")
          if (key %in% names(null_count)) {
            null_count[key] <- null_count[key] + 1L
          }
        }
      }
    }
  }

  ## empirical p(dynamic != null | permuted) per term-condition
  p_dynamic <- (null_count + 1) / (n_perm + 1)
  q_dynamic <- stats::p.adjust(p_dynamic, method = "BH")

  ## attach to observed dynamics (one p_dynamic per term-condition)
  dyn_obs$key <- paste(dyn_obs$term, dyn_obs$condition, sep = "##")
  dyn_obs$p_dynamic <- p_dynamic[dyn_obs$key]
  dyn_obs$q_dynamic <- q_dynamic[dyn_obs$key]
  dyn_obs$key <- NULL

  dyn_obs
}
