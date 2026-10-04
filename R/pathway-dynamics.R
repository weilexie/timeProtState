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
## The full rank-audited contrast family's moderated statistics drive the
## pathway summaries directly, retaining direction and model uncertainty at
## every timepoint. Significance is CAMERA/ROMER (ADR 0010);
## tps_calibrate_dynamics uses Bonferroni-minP then global BH (q_dynamic),
## not a permutation test of the shape label.

## ---- contrast-time parser --------------------------------------------------
## Resolve per-contrast (condition, time) metadata. Prefers the contrast_meta
## data.frame carried on the calibrated object (authoritative, set at
## construction by tps_make_contrasts). Falls back to string parsing of
## '<condition>_t<time>_vs_baseline' names for objects without metadata
## (bare matrices, legacy objects). Returns a data.frame with columns:
## contrast, condition, time.
.contrast_meta <- function(calibrated, contrast_names) {
  cm <- calibrated$contrast_meta
  if (is.data.frame(cm) &&
      all(c("contrast", "condition", "time") %in% names(cm)) &&
      length(contrast_names) == nrow(cm) &&
      all(cm$contrast %in% contrast_names) &&
      any(!is.na(cm$time))) {
    ## use authoritative metadata, aligned to contrast order
    cm <- cm[match(contrast_names, cm$contrast), , drop = FALSE]
    return(cm[, c("contrast", "condition", "time")])
  }
  ## fallback: parse from names (legacy / bare matrix path)
  .parse_contrast_time_string(contrast_names)
}

.parse_contrast_time_string <- function(contrast_names) {
  base <- sub("_vs_baseline$", "", contrast_names)
  cond <- sub("_t[^_]+$", "", base)
  time_str <- sub("^.*_t", "", base)
  time <- suppressWarnings(as.numeric(time_str))
  data.frame(contrast = contrast_names, condition = cond, time = time,
             stringsAsFactors = FALSE)
}

## Column contract (ADR 0010): calibrated objects carry q_directional /
## q_mixed. Permutation-era aliases q_net / q_rem are accepted as fallbacks.
.q_dir <- function(obj) {
  if (!is.null(obj$q_directional)) obj$q_directional else obj$q_net
}
.q_mix <- function(obj) {
  if (!is.null(obj$q_mixed)) obj$q_mixed else obj$q_rem
}
.p_dir <- function(obj) {
  if (!is.null(obj$p_directional)) obj$p_directional else obj$p_net
}
.p_mix <- function(obj) {
  if (!is.null(obj$p_mixed)) obj$p_mixed else obj$p_rem
}
.q_dir_long <- function(df) {
  if ("q_directional" %in% names(df)) df$q_directional else df$q_net
}

## ---- classify one pathway-condition trajectory -----------------------------
## Given NET/REM/COH/q_directional/q_mixed vectors ordered by time, classify
## the dynamic pattern. Shape rules are descriptive; they read per-time q
## values, not q_dynamic.
.classify_one_trajectory <- function(net, rem, coh, q_dir, q_mix, alpha, coh_cut) {
  n <- length(net)
  if (n == 0L) return("no_dynamic_signal")

  sig <- !is.na(q_dir) & q_dir < alpha
  sig_rem <- !is.na(q_mix) & q_mix < alpha

  if (!any(sig) && !any(sig_rem)) return("no_dynamic_signal")

  ## NA-safe helpers: comparisons that may involve NA use isTRUE()
  rem_n <- rem[n]; rem_1 <- rem[1L]; coh_1 <- coh[1L]

  ## Priority: specific patterns first (transient, gradual, delayed),
  ## then the general "early_directional" catch-all.

  ## transient_burst: REM peaks significantly in the middle, drops at the end.
  if (n >= 3L) {
    rem_mid <- rem[seq(2L, n - 1L)]
    sig_rem_mid <- sig_rem[seq(2L, n - 1L)]
    mid_max <- if (any(is.finite(rem_mid))) max(rem_mid, na.rm = TRUE) else -Inf
    if (any(sig_rem_mid) &&
        isTRUE(mid_max > rem_n * 1.5) &&
        isTRUE(mid_max > rem_1 * 1.5)) {
      return("transient_burst")
    }
  }

  ## delayed_opposition: the final supported NET direction reverses an earlier
  ## supported direction. A last-only signal is delayed, not opposition.
  if (n >= 2L && isTRUE(sig[n]) && any(sig[-n], na.rm = TRUE)) {
    earlier <- net[-n][sig[-n] & !is.na(net[-n])]
    if (length(earlier) && isTRUE(sign(net[n]) != 0) &&
        any(sign(earlier) == -sign(net[n]))) {
      return("delayed_opposition")
    }
  }

  ## last-only directional signal.
  if (n >= 2L && isTRUE(sig[n]) && !any(sig[-n], na.rm = TRUE)) {
    return("directional_delayed")
  }

  ## early_directional: first timepoint NET-significant with high COH and a
  ## persistent, directionally consistent signal. Checked BEFORE gradual_remodeling
  ## so that a growing coherent directional signal is labeled early_directional,
  ## not gradual_remodeling (which is for non-directional REM growth).
  sig_net <- net[sig & !is.na(net)]
  same_direction <- length(sig_net) >= 2L &&
    length(unique(sign(sig_net[sign(sig_net) != 0]))) == 1L
  if (isTRUE(sig[1L]) && isTRUE(coh_1 >= coh_cut) && same_direction) {
    return("early_directional")
  }

  ## gradual_remodeling: REM increases over time, signal builds progressively.
  ## Only reached if early_directional did not fire (i.e., NOT a coherent
  ## directional signal from the start).
  if (n >= 2L && isTRUE(rem_n > rem_1) && isTRUE(sig_rem[n]) &&
      isTRUE(rem_n > rem_1 * 1.2)) {
    return("gradual_remodeling")
  }

  ## any significant NET -- directional but not early
  if (any(sig, na.rm = TRUE)) return("directional_delayed")

  ## significant REM but never directional NET -- pure remodeling
  if (any(sig_rem, na.rm = TRUE) && !any(sig, na.rm = TRUE)) {
    return("sustained_remodeling")
  }

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
    q_dir <- .q_dir_long(sub)
    sig <- sub[!is.na(q_dir) & q_dir < alpha, , drop = FALSE]
    if (nrow(sig) >= 1L) {
      sig_signs <- sign(sig$NET)
      div_found <- FALSE
      ## Case 1: two or more significant conditions with opposite signs
      ## (strongest divergence signal — all conditions significant)
      if (nrow(sig) >= 2L && any(sig_signs > 0) && any(sig_signs < 0)) {
        div_found <- TRUE
      }
      ## Case 2: a significant condition has opposite sign from a
      ## non-significant condition with a non-trivial NET (magnitude floor
      ## prevents false positives from near-zero non-significant NETs)
      if (!div_found) {
        nonsig <- sub[is.na(q_dir) | q_dir >= alpha, , drop = FALSE]
        if (nrow(nonsig) > 0L) {
          if (any(sig_signs > 0) && any(nonsig$NET < -0.1, na.rm = TRUE))
            div_found <- TRUE
          if (!div_found && any(sig_signs < 0) &&
              any(nonsig$NET > 0.1, na.rm = TRUE))
            div_found <- TRUE
        }
      }
      if (div_found) n_div <- n_div + 1L
    }
  }
  n_div
}

## boolean wrapper for backward compatibility in tps_classify_pathway_dynamics
.detect_divergence <- function(traj_long, alpha) {
  .count_divergent_timepoints(traj_long, alpha) > 0L
}

## ---- tps_classify_pathway_dynamics -----------------------------------------
#' Assign rule-based trajectory shape labels to pathway dynamics
#'
#' Converts a calibrated pathway scores object (from
#' \code{\link{tps_calibrate_pathways}} with multiple timepoints) into
#' per-pathway-per-condition \strong{rule-based trajectory shape labels}
#' describing HOW the pathway remodels over time:
#' \code{early_directional}, \code{gradual_remodeling},
#' \code{delayed_opposition}, \code{transient_burst},
#' \code{directional_delayed}, \code{sustained_remodeling}, or
#' \code{no_dynamic_signal}.
#'
#' These labels are produced by a transparent, prioritized rule engine
#' (not a trained classifier). Each label is a descriptive phenotype:
#' \code{gradual_remodeling} is not proven to be strictly monotonic;
#' \code{transient_burst} is not proven to have a complete rise-then-fall;
#' \code{sustained_remodeling} is not proven significant at every timepoint;
#' \code{directional_delayed} does not prove a delayed mechanism. The rules
#' are auditable (see the source of \code{.classify_one_trajectory}) and
#' intended for exploratory trajectory interpretation, not as class-specific
#' statistical inferences. Use \code{\link{tps_calibrate_dynamics}} for a
#' trajectory-level Bonferroni-minP + global BH combination
#' (\code{p_dynamic}/\code{q_dynamic}); that combination is not a shape
#' test and is independent of the shape label.
#'
#' Cross-condition NET opposition is reported as the boolean column
#' \code{is_divergent}, not as a \code{dynamic_state} label.
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
#'   COH, q_directional, q_mixed, dynamic_state} (long format, one row per
#'   pathway-condition-time), plus a summary column \code{is_divergent} marking
#'   pathways where conditions diverge. Legacy aliases \code{q_net} and
#'   \code{q_rem} duplicate \code{q_directional} and \code{q_mixed}.
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
#' tps_classify_pathway_dynamics(pws_c)
#' @export
tps_classify_pathway_dynamics <- function(calibrated, alpha = 0.05, coh_cut = 0.5) {
  if (!inherits(calibrated, "tps_pathway_calibrated")) {
    tps_abort(paste0("calibrated must be a tps_pathway_calibrated object; ",
                     "run tps_calibrate_pathways first."))
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    tps_abort("alpha must be in (0,1).")
  }
  if (!is.numeric(coh_cut) || length(coh_cut) != 1L || is.na(coh_cut) ||
      !is.finite(coh_cut) || coh_cut < 0 || coh_cut > 1) {
    tps_abort("coh_cut must be in [0,1].")
  }
  NET <- calibrated$scores$NET
  REM <- calibrated$scores$REM
  COH <- calibrated$scores$COH
  q_dir <- .q_dir(calibrated)
  q_mix <- .q_mix(calibrated)
  if (is.null(q_dir) || is.null(q_mix)) {
    tps_abort(paste0(
      "calibrated object is missing q_directional; ",
      "run tps_calibrate_pathways first."
    ))
  }
  terms <- rownames(NET)
  contrasts <- colnames(NET)
  if (is.null(colnames(q_dir)) || !all(contrasts %in% colnames(q_dir))) {
    if (!is.null(calibrated$q_net) &&
        all(contrasts %in% colnames(calibrated$q_net))) {
      q_dir <- calibrated$q_net
    }
  }
  if (is.null(colnames(q_mix)) || !all(contrasts %in% colnames(q_mix))) {
    if (!is.null(calibrated$q_rem) &&
        all(contrasts %in% colnames(calibrated$q_rem))) {
      q_mix <- calibrated$q_rem
    }
  }
  if (is.null(colnames(q_dir)) || !all(contrasts %in% colnames(q_dir)) ||
      is.null(colnames(q_mix)) || !all(contrasts %in% colnames(q_mix))) {
    tps_abort(paste0(
      "Could not resolve timepoints from contrast metadata. Expected ",
      "family='vs_baseline' (format '<condition>_t<time>_vs_baseline')."
    ))
  }
  q_dir <- q_dir[terms, contrasts, drop = FALSE]
  q_mix <- q_mix[terms, contrasts, drop = FALSE]

  ## parse time + condition from contrast names
  meta <- .contrast_meta(calibrated, contrasts)
  if (anyNA(meta$time)) {
    tps_abort(paste0(
      "Could not resolve timepoints from contrast metadata. Expected ",
      "family='vs_baseline' (format '<condition>_t<time>_vs_baseline')."
    ))
  }
  if (any(grepl("^endpoint_", meta$condition))) {
    tps_abort(paste0(
      "endpoint-family contrasts detected. tps_classify_pathway_dynamics ",
      "requires family='vs_baseline' in tps_pathway_scores; endpoint family ",
      "produces only one timepoint per condition."
    ))
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
        q_dir = q_dir[i, idx], q_mix = q_mix[i, idx],
        alpha = alpha, coh_cut = coh_cut)
      for (k in seq_along(idx)) {
        j <- idx[k]
        rows[[length(rows) + 1L]] <- data.frame(
          term = terms[i], condition = cc, time = meta$time[j],
          NET = NET[i, j], REM = REM[i, j], COH = COH[i, j],
          q_directional = q_dir[i, j], q_mixed = q_mix[i, j],
          q_net = q_dir[i, j], q_rem = q_mix[i, j],
          dynamic_state = dyn, stringsAsFactors = FALSE)
      }
    }
  }
  long <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame(
    term = character(0), condition = character(0), time = numeric(0),
    NET = numeric(0), REM = numeric(0), COH = numeric(0),
    q_directional = numeric(0), q_mixed = numeric(0),
    q_net = numeric(0), q_rem = numeric(0),
    dynamic_state = character(0), stringsAsFactors = FALSE)

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
  attr(long, "provenance") <- list(
    stage = "classify_pathway_dynamics", alpha = alpha, coh_cut = coh_cut,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    package_version = as.character(getNamespaceVersion("timeProtState")))
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
#' over-represented in pathways showing that dynamic pattern.
#'
#' \strong{Circularity caveat}: the hypergeometric p-value assumes the
#' protein state labels and pathway dynamic labels are independent. In
#' practice both are derived from the same protein set and the same
#' underlying fit, so the independence assumption is violated. The results
#' should be interpreted as \strong{descriptive consistency scores}, not as
#' validated inferential enrichment. To obtain a valid inferential test, an
#' independent partition (e.g., sample split or pathway permutation null)
#' would be required; this function does not provide one.
#'
#' @param states a \code{tps_states} object.
#' @param pathway_dynamics a data.frame from \code{tps_classify_pathway_dynamics}.
#' @param term2gene pathway membership table (columns: term, gene).
#' @param protein_ids optional character vector of pathway-space ids. If
#'   \code{NULL}, \code{feature_id} is used. A \emph{named} vector is
#'   reordered to \code{states$state_table$feature_id} (safer than a
#'   positional vector, because the state table is sorted). An unnamed
#'   vector must already be aligned to the state-table rows.
#' @return a data.frame with columns: \code{protein_state},
#'   \code{pathway_dynamic}, \code{n_proteins_in_state},
#'   \code{n_pathways_in_dynamic}, \code{observed_overlap},
#'   \code{expected_overlap}, \code{enrichment},
#'   \code{p_hypergeometric}, \code{q_hypergeometric}.
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
#' dyn <- tps_classify_pathway_dynamics(pws_c)
#' ids <- stats::setNames(tps_example_annotation$locus_tag,
#'                        tps_example_annotation$protein_id)
#' link <- tps_link_states_pathways(tps_classify_states(fit), dyn,
#'                                  tps_example_term2gene, protein_ids = ids)
#' @export
tps_link_states_pathways <- function(states, pathway_dynamics, term2gene,
                                      protein_ids = NULL) {
  if (!inherits(states, "tps_states")) {
    tps_abort("states must be a tps_states object.")
  }
  if (!is.data.frame(pathway_dynamics) ||
      !all(c("term", "dynamic_state") %in% names(pathway_dynamics))) {
    tps_abort("pathway_dynamics must be from tps_classify_pathway_dynamics.")
  }
  if (!is.data.frame(term2gene) ||
      !all(c("term", "gene") %in% names(term2gene))) {
    tps_abort("term2gene must have columns 'term' and 'gene'.")
  }
  st <- states$state_table
  if (!is.data.frame(st) || !all(c("feature_id", "state") %in% names(st)) ||
      nrow(st) == 0L) {
    tps_abort("states$state_table must contain feature_id and state rows.")
  }
  if (is.null(protein_ids)) {
    protein_ids <- as.character(st$feature_id)
  } else {
    id_names <- names(protein_ids)
    protein_ids <- as.character(unname(protein_ids))
    if (!is.null(id_names) && any(nzchar(id_names))) {
      names(protein_ids) <- id_names
      mapped <- protein_ids[as.character(st$feature_id)]
      if (anyNA(mapped)) {
        tps_abort("named protein_ids must cover every feature_id in states.")
      }
      protein_ids <- unname(as.character(mapped))
    }
  }
  if (length(protein_ids) != nrow(st) || anyNA(protein_ids) ||
      any(!nzchar(protein_ids)) || anyDuplicated(protein_ids)) {
    tps_abort(paste0(
      "protein_ids must be unique, non-empty, non-missing, and aligned ",
      "with state rows (or named by feature_id)."
    ))
  }
  st$.pathway_id <- protein_ids

  pathway_dynamics <- unique(data.frame(
    term = as.character(pathway_dynamics$term),
    dynamic_state = as.character(pathway_dynamics$dynamic_state),
    stringsAsFactors = FALSE
  ))
  valid_dyn <- !is.na(pathway_dynamics$term) & nzchar(pathway_dynamics$term) &
    !is.na(pathway_dynamics$dynamic_state) &
    nzchar(pathway_dynamics$dynamic_state)
  pathway_dynamics <- pathway_dynamics[valid_dyn, , drop = FALSE]
  term2gene <- unique(data.frame(
    term = as.character(term2gene$term),
    gene = as.character(term2gene$gene),
    stringsAsFactors = FALSE
  ))
  valid_map <- !is.na(term2gene$term) & nzchar(term2gene$term) &
    !is.na(term2gene$gene) & nzchar(term2gene$gene)
  term2gene <- term2gene[valid_map, , drop = FALSE]
  if (nrow(pathway_dynamics) == 0L || nrow(term2gene) == 0L) {
    tps_abort("pathway_dynamics and term2gene must contain valid rows.")
  }

  ## analysis universe: proteins observed in the state table (NOT library)
  universe <- unique(st$.pathway_id)
  n_universe <- length(universe)
  if (n_universe > 0L && length(unique(term2gene$gene)) > 0L) {
    n_match <- length(intersect(universe, unique(term2gene$gene)))
    if (n_match == 0L) {
      tps_abort(paste0(
        "protein_ids do not match any term2gene$gene; check namespace ",
        "(e.g., locus_tag vs protein_id) or pass a named vector. ",
        "Universe: ", paste(head(universe, 3), collapse = ", "),
        "; term2gene genes: ", paste(head(unique(term2gene$gene), 3), collapse = ", "), "."))
    }
  }

  ## protein states
  prot_states <- unique(st$state)

  rows <- list()
  for (ps in prot_states) {
    ps_proteins <- unique(st$.pathway_id[st$state == ps])
    n_ps <- length(ps_proteins)
    for (dyn in unique(pathway_dynamics$dynamic_state)) {
      ## A pathway may carry different condition-specific dynamic labels. Keep
      ## every observed term-state association instead of selecting an
      ## arbitrary modal label.
      dyn_pathways <- unique(
        pathway_dynamics$term[pathway_dynamics$dynamic_state == dyn]
      )
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
        observed_overlap = overlap, expected_overlap = expected,
        enrichment = enrich,
        p_hypergeometric = p_hyper,
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  out$q_hypergeometric <- stats::p.adjust(out$p_hypergeometric, method = "BH")
  attr(out, "provenance") <- list(
    stage = "link_states_pathways",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    package_version = as.character(getNamespaceVersion("timeProtState")))
  out
}

## ============================================================================
## Trajectory p-values (Bonferroni-minP) + divergence scoring
## ============================================================================
##
## p_dynamic is a Bonferroni-minP combination of per-timepoint p_directional
## and p_mixed; q_dynamic is one global BH across pathway x condition.
## Shape labels stay the descriptive rules on per-time q_directional/q_mixed.

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
    q_dir <- .q_dir_long(sub)
    for (i in seq_len(nrow(sub) - 1L)) {
      for (j in seq(i + 1L, nrow(sub))) {
        ## only count if at least one is significant
        if (!is.na(q_dir[i]) && !is.na(q_dir[j])) {
          if (q_dir[i] < alpha || q_dir[j] < alpha) {
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
#' dyn <- tps_classify_pathway_dynamics(pws_c)
#' div <- tps_divergence_score(pws_c, dyn)
#' @export
tps_divergence_score <- function(calibrated, dyn = NULL, alpha = 0.05) {
  if (!inherits(calibrated, "tps_pathway_calibrated")) {
    tps_abort("calibrated must be a tps_pathway_calibrated object.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    tps_abort("alpha must be one finite value in (0,1).")
  }
  if (is.null(dyn)) dyn <- tps_classify_pathway_dynamics(calibrated, alpha = alpha)
  has_q <- is.data.frame(dyn) &&
    (("q_directional" %in% names(dyn)) || ("q_net" %in% names(dyn)))
  if (!is.data.frame(dyn) ||
      !all(c("term", "condition", "time", "NET") %in% names(dyn)) ||
      !has_q) {
    tps_abort("dyn must be a pathway dynamics table.")
  }
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
  result <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame(
    term = character(0), divergence_score = numeric(0),
    n_divergent_timepoints = integer(0),
    max_divergent_pair = character(0),
    max_divergent_time = numeric(0),
    stringsAsFactors = FALSE)
  attr(result, "provenance") <- list(
    stage = "divergence_score", alpha = alpha,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    package_version = as.character(getNamespaceVersion("timeProtState")))
  result
}

## ---- Bonferroni-minP trajectory p-values -----------------------------------
#' Combine per-timepoint directional and mixed p-values into a trajectory p-value
#'
#' For each pathway-condition trajectory, collects the per-timepoint
#' \code{p_directional} and \code{p_mixed} values and applies a Bonferroni-minP
#' combination: \code{p_traj = min(1, m * min(p_j))}, where \code{m} is the
#' length of that p-value vector. \code{q_dynamic} is one global BH
#' adjustment across all pathway-condition \code{p_traj} values (not a
#' per-condition BH).
#'
#' This is \strong{not} a shape test. \code{p_dynamic}/\code{q_dynamic} flag
#' whether any timepoint shows calibrated directional or mixed evidence;
#' they do not certify a particular \code{dynamic_state} label. Shape labels
#' remain the descriptive rules in \code{\link{tps_classify_pathway_dynamics}},
#' which read per-time \code{q_directional}/\code{q_mixed} (not
#' \code{q_dynamic}). Overlapping pathways share proteins, so
#' \code{q_dynamic} should not be interpreted as a strict FDR among
#' overlapping gene sets.
#'
#' @param calibrated a \code{tps_pathway_calibrated} object from
#'   \code{\link{tps_calibrate_pathways}}. Must contain \code{p_directional}
#'   (or permutation alias \code{p_net}) and \code{p_mixed} (or \code{p_rem}).
#' @param n_perm ignored; retained for backward compatibility. Trajectory
#'   p-values are no longer permutation-based.
#' @param alpha numeric; q-value threshold for the descriptive shape rules
#'   (default 0.05). Does not enter \code{p_dynamic}.
#' @param coh_cut numeric; COH threshold for shape rules (default 0.5).
#' @param seed ignored; retained for backward compatibility.
#' @return a data.frame with the original dynamics plus \code{p_dynamic}
#'   (Bonferroni-minP of per-time directional and mixed p-values) and
#'   \code{q_dynamic} (global BH across pathway x condition).
#'
#' @examples
#' sim <- tps_simulate_timecourse(
#'   n_features = 12, conditions = c("A", "B"), times = c(1, 3),
#'   n_replicates = 2, effect_mix = c(null = 1), seed = 1)
#' se <- tps_import(sim$abundance, sim$col_data)
#' design <- tps_make_design(se)
#' contrasts <- tps_make_contrasts(design, endpoint = 3)
#' fit <- tps_fit(se, design, contrasts)
#' t2g <- data.frame(term = "pw1", gene = rownames(sim$abundance)[1:4],
#'                   stringsAsFactors = FALSE)
#' pws <- tps_pathway_scores(fit, t2g, family = "vs_baseline", min_size = 3)
#' pws_c <- tps_calibrate_pathways(pws, mode = "unstratified",
#'                                 n_perm = 2, seed = 1)
#' dyn_cal <- tps_calibrate_dynamics(pws_c)
#' nrow(dyn_cal)
#' @export
tps_calibrate_dynamics <- function(calibrated, n_perm = 500L, alpha = 0.05,
                                    coh_cut = 0.5, seed = 1L) {
  if (!inherits(calibrated, "tps_pathway_calibrated")) {
    tps_abort("calibrated must be a tps_pathway_calibrated object.")
  }
  ## n_perm and seed are ignored (no permutation loop). Kept in the signature
  ## so existing callers do not error on unused arguments.
  p_dir <- .p_dir(calibrated)
  p_mix <- .p_mix(calibrated)
  if (is.null(p_dir)) {
    tps_abort(paste0(
      "calibrated object is missing p_directional; ",
      "run tps_calibrate_pathways first."
    ))
  }
  if (is.null(p_mix)) {
    tps_abort(paste0(
      "calibrated object is missing p_mixed; ",
      "run tps_calibrate_pathways first."
    ))
  }

  ## Shape labels stay the descriptive per-time q_directional/q_mixed rules.
  dyn_obs <- tps_classify_pathway_dynamics(calibrated, alpha = alpha,
                                            coh_cut = coh_cut)

  contrasts <- colnames(calibrated$scores$NET)
  terms <- rownames(calibrated$scores$NET)
  p_dir <- p_dir[terms, contrasts, drop = FALSE]
  p_mix <- p_mix[terms, contrasts, drop = FALSE]
  meta <- .contrast_meta(calibrated, contrasts)
  if (anyNA(meta$time)) {
    tps_abort("Could not resolve every contrast as '<condition>_t<time>_vs_baseline' from metadata.")
  }
  if (any(grepl("^endpoint_", meta$condition))) {
    tps_abort("endpoint-family contrasts detected; use family='vs_baseline' in tps_pathway_scores.")
  }

  ## Bonferroni-minP over p_directional and p_mixed at this condition's
  ## timepoints; one p_traj per pathway x condition.
  keys <- character(0)
  p_traj <- numeric(0)
  for (i in seq_along(terms)) {
    for (cc in unique(meta$condition)) {
      idx <- which(meta$condition == cc)
      if (length(idx) == 0L) next
      p_j <- c(as.numeric(p_dir[i, idx]), as.numeric(p_mix[i, idx]))
      m <- length(p_j)
      val <- if (!any(is.finite(p_j))) {
        NA_real_
      } else {
        min(1, m * min(p_j, na.rm = TRUE))
      }
      keys <- c(keys, paste(terms[i], cc, sep = "##"))
      p_traj <- c(p_traj, val)
    }
  }
  names(p_traj) <- keys
  ## One global BH across all pathway x condition p_traj (not per-condition).
  q_dyn <- stats::p.adjust(p_traj, method = "BH")

  dyn_obs$key <- paste(dyn_obs$term, dyn_obs$condition, sep = "##")
  dyn_obs$p_dynamic <- unname(p_traj[dyn_obs$key])
  dyn_obs$q_dynamic <- unname(q_dyn[dyn_obs$key])
  dyn_obs$key <- NULL
  ## preserve the classify_pathway_dynamics provenance from the internal call,
  ## then append the calibrate_dynamics entry (do not overwrite).
  classify_prov <- attr(dyn_obs, "provenance")
  attr(dyn_obs, "provenance") <- list(
    classify_pathway_dynamics = classify_prov,
    calibrate_dynamics = list(
      stage = "calibrate_dynamics",
      method = "bonferroni_minp_global_bh",
      n_perm = n_perm, seed = seed, ignored = c("n_perm", "seed"),
      alpha = alpha, coh_cut = coh_cut,
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      package_version = as.character(getNamespaceVersion("timeProtState"))))
  dyn_obs
}
