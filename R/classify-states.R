## ============================================================================
## Classify states: evidence-decomposed protein dynamic states
## ============================================================================
##
## Implements STATISTICAL_CONTRACT.md section 6. The classification NEVER
## collapses to a black box first. It computes three transparent evidence axes:
##   1. temporal_evidence    (from average_time / within_omnibus families)
##   2. condition_heterogeneity (from heterogeneity family)
##   3. direction_pattern    (from endpoint family, cross-condition signs)
## Then combines them into a final `state` label, but always retains the three
## raw axes plus the rule manifest for auditability.
##
## Significance default: adj.P.Val < 0.05 with lfc_threshold = 0 (DECISIONS D2).
## The legacy_msm profile (raw P + 0.74) is achievable via tps_decision_rule().

## ---- decision rule constructor --------------------------------------------
#' Define a decision rule for state classification
#'
#' Constructs a rule list controlling how evidence is thresholded. The package
#' default uses Benjamini-Hochberg FDR < 0.05 with no logFC gate; the
#' \code{legacy_msm} profile (raw P < 0.05 and |logFC| >= 0.74) reproduces the
#' historical dormancy analysis.
#'
#' @param significance character; \code{"adj.P.Val"} (default) or \code{"P.Value"}.
#' @param alpha numeric; significance threshold (default 0.05).
#' @param lfc_threshold numeric; minimum |logFC| (default 0).
#' @param endpoint the endpoint timepoint used for direction; if \code{NULL},
#'   taken from the contrasts object.
#' @return a list of class \code{tps_decision_rule}.
#'
#' @examples
#' tps_decision_rule()
#' tps_decision_rule(significance = "P.Value", alpha = 0.05, lfc_threshold = 0.74)
#'
#' @export
tps_decision_rule <- function(significance = c("adj.P.Val", "P.Value"),
                                alpha = 0.05, lfc_threshold = 0, endpoint = NULL) {
  significance <- match.arg(significance)
  if (length(alpha) != 1L || !is.numeric(alpha) || is.na(alpha) ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    tps_abort("alpha must be in (0,1) and supplied as one finite value.")
  }
  if (length(lfc_threshold) != 1L || !is.numeric(lfc_threshold) ||
      is.na(lfc_threshold) || !is.finite(lfc_threshold) || lfc_threshold < 0) {
    tps_abort("lfc_threshold must be one non-negative finite value.")
  }
  if (!is.null(endpoint) &&
      (length(endpoint) != 1L || is.na(endpoint) || !is.finite(endpoint))) {
    tps_abort("endpoint must be NULL or one finite value.")
  }
  structure(
    list(significance = significance, alpha = alpha,
         lfc_threshold = lfc_threshold, endpoint = endpoint),
    class = "tps_decision_rule"
  )
}

## ---- per-feature hit helper ----------------------------------------------
.is_hit <- function(stat_row, rule) {
  sig_col <- rule$significance
  if (!sig_col %in% names(stat_row)) return(FALSE)
  sig_ok <- !is.na(stat_row[[sig_col]]) && stat_row[[sig_col]] < rule$alpha
  lfc_ok <- is.na(stat_row$logFC) || abs(stat_row$logFC) >= rule$lfc_threshold
  sig_ok && lfc_ok
}

## ---- tps_classify_states ---------------------------------------------------
#' Classify proteins into dynamic evidence states
#'
#' Decomposes per-protein statistical evidence into three transparent axes
#' (temporal, condition heterogeneity, direction pattern) and combines them
#' into a final state label. All three raw axes and the rule manifest are
#' retained for auditability.
#'
#' @param fit a \code{tps_fit} object.
#' @param rule a \code{tps_decision_rule} from \code{\link{tps_decision_rule}}.
#'   If \code{NULL}, the package default (FDR < 0.05, lfc 0) is used.
#' @return a list of class \code{tps_states} with \code{state_table} (per-feature
#'   axes + final state), \code{rule_manifest}, and \code{provenance}.
#'
#' @examples
#' data(tps_example, package = "timeProtState")
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' design <- tps_make_design(se)
#' contrasts <- tps_make_contrasts(design, endpoint = 7)
#' fit <- tps_fit(se, design, contrasts)
#' states <- tps_classify_states(fit)
#' head(states$state_table)
#'
#' @export
tps_classify_states <- function(fit, rule = NULL) {
  if (!inherits(fit, "tps_fit")) tps_abort("fit must be a tps_fit object.")
  if (is.null(rule)) rule <- tps_decision_rule()
  if (!inherits(rule, "tps_decision_rule")) {
    tps_abort("rule must be a tps_decision_rule object.")
  }
  if (!is.null(rule$endpoint) &&
      !identical(as.character(rule$endpoint),
                 as.character(fit$contrasts$endpoint))) {
    tps_abort("rule endpoint does not match the fitted endpoint; rebuild contrasts and refit.")
  }
  arg_capture <- list(call = "tps_classify_states",
                       rule = rule,
                       timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  stats <- fit$statistics
  conditions <- fit$contrasts$conditions

  ## ---- axis 1: temporal_evidence (from average_time + within_omnibus) ----
  ## A feature has temporal evidence if average_time OR any within_omnibus
  ## condition-family is significant.
  avg_rows <- stats[stats$family == "average_time" & stats$test == "F_omnibus", ]
  within_rows <- stats[stats$family == "within_omnibus" & stats$test == "F_omnibus", ]
  tmp_avg <- .sig_per_feature(avg_rows, rule)
  tmp_within <- .sig_per_feature(within_rows, rule)

  ## ---- axis 2: condition_heterogeneity ----------------------------------
  het_rows <- stats[stats$family == "heterogeneity" & stats$test == "F_omnibus", ]
  het <- .sig_per_feature(het_rows, rule)

  ## ---- axis 3: direction_pattern (from endpoint family) -----------------
  end_rows <- stats[stats$family == "endpoint" & stats$test == "t", ]
  end_dir <- .endpoint_directions(end_rows, rule)

  ## ---- reversal detection (from vs_baseline family) ---------------------
  ## A temporal reversal = within a condition, significant vs-baseline logFC
  ## changes sign across timepoints (up-then-down or vice versa), with
  ## effect-size gate (|logFC| >= lfc_threshold at the reversal points).
  vsbl_rows <- stats[stats$family == "vs_baseline" & stats$test == "t", ]
  reversal <- .detect_reversal(vsbl_rows, conditions, rule)

  ## ---- assemble per-feature table ---------------------------------------
  all_features <- sort(unique(stats$feature_id))
  st <- data.frame(
    feature_id = all_features,
    temporal_evidence = .axis_temporal(all_features, tmp_avg, tmp_within),
    condition_heterogeneity = .axis_heterogeneity(all_features, het),
    stringsAsFactors = FALSE)
  ## join endpoint directions per condition
  for (cond in conditions) {
    by_condition <- end_dir[end_dir$condition == cond, , drop = FALSE]
    idx <- match(st$feature_id, by_condition$feature_id)
    direction <- by_condition$dir[idx]
    supported <- by_condition$sig[idx]
    direction[is.na(direction)] <- "ns"
    supported[is.na(supported)] <- FALSE
    st[[paste0("endpoint_dir_", cond)]] <- direction
    st[[paste0("endpoint_sig_", cond)]] <- supported
  }
  st$direction_pattern <- .axis_direction(st, conditions)
  st$is_reversal <- reversal[match(st$feature_id, names(reversal))]
  st$is_reversal[is.na(st$is_reversal)] <- FALSE
  st$state <- .combine_state(st)

  ## rule manifest for audit
  rule_manifest <- data.frame(
    axis = c("temporal","heterogeneity","direction","reversal","state"),
    source_family = c("average_time | within_omnibus","heterogeneity",
                       "endpoint","vs_baseline (sign-change + effect-size gate)",
                       "combination"),
    significance = rule$significance, alpha = rule$alpha,
    lfc_threshold = rule$lfc_threshold, stringsAsFactors = FALSE)

  out <- structure(
    list(
      state_table = st,
      rule_manifest = rule_manifest,
      rule = rule,
      conditions = conditions,
      provenance = list(.make_provenance(arg_capture, fit, stage = "classify_states"))
    ),
    class = "tps_states"
  )
  out
}

## ---- axis helpers ---------------------------------------------------------
.sig_per_feature <- function(rows, rule) {
  if (nrow(rows) == 0L) {
    return(setNames(logical(0), character(0)))
  }
  sig_col <- rule$significance
  hit <- !is.na(rows[[sig_col]]) & rows[[sig_col]] < rule$alpha
  if (rule$lfc_threshold > 0 && "logFC" %in% names(rows)) {
    ## Omnibus F rows have no single logFC; the effect-size gate applies only
    ## where a directional coefficient is defined.
    hit <- hit & (is.na(rows$logFC) | abs(rows$logFC) >= rule$lfc_threshold)
  }
  by_feature <- split(hit, rows$feature_id)
  vapply(by_feature, any, logical(1L), na.rm = TRUE)
}

.axis_temporal <- function(features, tmp_avg, tmp_within) {
  a <- tmp_avg[features]; w <- tmp_within[features]
  a[is.na(a)] <- FALSE; w[is.na(w)] <- FALSE
  out <- ifelse(a & w, "average_and_within",
          ifelse(a, "average_time",
          ifelse(w, "within_condition", "none")))
  out[is.na(out)] <- "none"
  unname(out)
}

.axis_heterogeneity <- function(features, het) {
  h <- het[features]; h[is.na(h)] <- FALSE
  unname(ifelse(h, "condition_dependent", "homogeneous_or_unresolved"))
}

.endpoint_directions <- function(end_rows, rule) {
  if (nrow(end_rows) == 0L) {
    return(data.frame(feature_id = character(0),
                      condition = character(0), dir = character(0),
                      sig = logical(0)))
  }
  dir <- sig <- rep(NA, nrow(end_rows))
  for (i in seq_len(nrow(end_rows))) {
    s <- .is_hit(end_rows[i, ], rule)
    sig[i] <- s
    if (s) {
      dir[i] <- ifelse(end_rows$logFC[i] > 0, "up", "down")
    } else {
      dir[i] <- "ns"
    }
  }
  parsed <- .parse_endpoint_contrast(end_rows$contrast)
  data.frame(feature_id = end_rows$feature_id, condition = parsed$condition,
             dir = dir, sig = sig,
             stringsAsFactors = FALSE)
}

.parse_endpoint_contrast <- function(contrast_names) {
  base <- sub("^endpoint_", "", contrast_names)
  base <- sub("_vs_baseline$", "", base)
  condition <- sub("_t[^_]+$", "", base)
  time <- suppressWarnings(as.numeric(sub("^.*_t", "", base)))
  list(condition = condition, time = time)
}

.axis_direction <- function(st, conditions) {
  cols <- paste0("endpoint_dir_", conditions)
  dirs <- as.data.frame(st[, cols, drop = FALSE])
  apply(dirs, 1L, function(d) {
    ds <- d[d != "ns"]
    if (length(ds) == 0L) return("no_supported_change")
    if (length(unique(ds)) > 1L) return("cross_condition_opposite")
    if (unique(ds) == "up") return("concordant_up")
    "concordant_down"
  })
}

## ---- reversal detection ---------------------------------------------------
## For each feature, check if any condition has significant vs-baseline logFC
## that changes sign across timepoints (with effect-size gate).
.detect_reversal <- function(vsbl_rows, conditions, rule) {
  if (nrow(vsbl_rows) == 0L) return(setNames(logical(0), character(0)))
  feats <- unique(vsbl_rows$feature_id)
  result <- setNames(rep(FALSE, length(feats)), feats)
  ## parse condition + time from contrast name "Cond_tX_vs_baseline"
  parsed <- .parse_vsbl_contrast(vsbl_rows$contrast)
  vsbl_rows$.condition <- parsed$condition
  vsbl_rows$.time <- parsed$time
  for (f in feats) {
    fr <- vsbl_rows[vsbl_rows$feature_id == f, , drop = FALSE]
    if (nrow(fr) < 2L) next
    for (cc in conditions) {
      fr_c <- fr[fr$.condition == cc, , drop = FALSE]
      if (nrow(fr_c) < 2L) next
      fr_c <- fr_c[order(fr_c$.time), ]
      ## determine hit + sign per timepoint (with effect-size gate)
      hits <- vapply(seq_len(nrow(fr_c)), function(i) .is_hit(fr_c[i, ], rule), logical(1))
      if (sum(hits) < 2L) next
      signs <- sign(fr_c$logFC[hits])
      ## reversal if both + and - signs present among significant timepoints
      if (any(signs > 0) && any(signs < 0)) {
        result[f] <- TRUE
        break
      }
    }
  }
  result
}

.parse_vsbl_contrast <- function(contrast_names) {
  ## contrast name format: "Condition_tX_vs_baseline"
  ## Use regex (consistent with .parse_contrast_time in pathway-dynamics.R)
  ## to handle condition names containing "_t<digits>" patterns.
  base <- sub("_vs_baseline$", "", contrast_names)
  condition <- sub("_t[^_]+$", "", base)
  time_str <- sub("^.*_t", "", base)
  time <- suppressWarnings(as.numeric(time_str))
  list(condition = condition, time = time)
}

## ---- final state combination ---------------------------------------------
.combine_state <- function(st) {
  tmp <- st$temporal_evidence
  het <- st$condition_heterogeneity
  dir <- st$direction_pattern
  is_rev <- if ("is_reversal" %in% names(st)) st$is_reversal else FALSE

  out <- character(nrow(st))
  has_temporal <- tmp != "none"
  is_het <- het == "condition_dependent"
  is_opposite <- dir == "cross_condition_opposite"

  out[!has_temporal & !is_het] <- "stable_or_unresolved"
  out[has_temporal & !is_het & dir %in% c("concordant_up","concordant_down")] <-
    ifelse(dir[has_temporal & !is_het & dir %in% c("concordant_up","concordant_down")] == "concordant_up",
           "shared_dynamic_up", "shared_dynamic_down")
  out[has_temporal & is_het] <- "condition_dependent_dynamic"
  out[!has_temporal & is_het] <- "endpoint_only"
  out[is_opposite] <- "cross_condition_opposite"
  ## temporal_reversal: significant sign-change across timepoints within a condition
  out[is_rev] <- "temporal_reversal"
  out[out == ""] <- "stable_or_unresolved"
  out
}

## ---- methods ---------------------------------------------------------------
#' @export
print.tps_states <- function(x, ...) {
  cat("<timeProtState states>\n")
  cat(sprintf("  features: %d   conditions: %s\n",
              nrow(x$state_table), paste(x$conditions, collapse = ", ")))
  tab <- table(x$state_table$state)
  cat("  state distribution:\n")
  for (nm in names(sort(tab, decreasing = TRUE))) {
    cat(sprintf("    %-30s %d\n", nm, tab[nm]))
  }
  invisible(x)
}

#' State table accessor
#'
#' Returns the per-feature state table (evidence axes + final state) from a
#' \code{tps_states} object.
#'
#' @param states a \code{tps_states} object.
#' @return a data.frame with one row per feature, containing the three
#'   evidence axes (\code{temporal_evidence}, \code{condition_heterogeneity},
#'   \code{direction_pattern}), per-condition endpoint direction columns, and
#'   the combined \code{state} label.
#'
#' @examples
#' data(tps_example, package = "timeProtState")
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' design <- tps_make_design(se)
#' contrasts <- tps_make_contrasts(design, endpoint = 7)
#' fit <- tps_fit(se, design, contrasts)
#' st <- tps_classify_states(fit)
#' head(tps_state_table(st))
#'
#' @export
tps_state_table <- function(states) {
  if (!inherits(states, "tps_states")) tps_abort("states must be a tps_states object.")
  states$state_table
}
