## ============================================================================
## Plotting: diagnostic and result visualizations
## ============================================================================
##
## Implements STATISTICAL_CONTRACT.md section 7 and Goal 13.
## Principles (v2 spec section 13):
##   - every plot function returns a ggplot object (composable, themeable)
##   - plot functions NEVER re-fit; they consume already-computed result objects
##   - color-blind friendly palettes
##   - data underlyng each plot is extractable (returned as attribute)
##
## ggplot2 is in Suggests; ComplexHeatmap (for heatmaps) also Suggests.

## ---- helper: attach result data as attribute for extraction ----------------
.with_data <- function(p, data) {
  attr(p, "data") <- data
  p
}

## ---- color palette ---------------------------------------------------------
.tps_colors <- function() {
  ## Okabe-Ito color-blind friendly palette
  c(up = "#D55E00", down = "#0072B2", ns = "grey70",
    sig = "#CC79A7", directional_up = "#D55E00",
    directional_down = "#0072B2", mixed_remodeling = "#009E73",
    no_strong_state = "grey70",
    directional_with_internal_opposition = "#CC79A7")
}

## ===========================================================================
## tps_plot_volcano
## ===========================================================================
#' Volcano plot from a fit
#'
#' Produces a ggplot volcano plot of -log10(P) vs logFC for a single contrast
#' from a \code{tps_fit}. Consumes already-computed statistics; does not re-fit.
#'
#' @param fit a \code{tps_fit} object.
#' @param contrast the contrast name to plot (must be a single-contrast t-test).
#' @param significance character; \code{"adj.P.Val"} (default) or \code{"P.Value"}.
#' @param alpha numeric; significance threshold (default 0.05).
#' @param lfc_threshold numeric; |logFC| threshold (default 0).
#' @return a ggplot object with the underlying data attached as attribute
#'   \code{"data"}.
#'
#' @examples
#' data(tps_example, package = "timeProtState")
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' design <- tps_make_design(se)
#' contrasts <- tps_make_contrasts(design, endpoint = 7)
#' fit <- tps_fit(se, design, contrasts)
#' p <- tps_plot_volcano(fit, contrast = "endpoint_Acid_t7_vs_baseline")
#'
#' @export
tps_plot_volcano <- function(fit, contrast, significance = "adj.P.Val",
                              alpha = 0.05, lfc_threshold = 0) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    tps_abort("ggplot2 is required for plotting; install it or it is in Suggests.")
  }
  if (!inherits(fit, "tps_fit")) tps_abort("fit must be a tps_fit object.")
  if (length(significance) != 1L ||
      !significance %in% c("adj.P.Val", "P.Value")) {
    tps_abort("significance must be 'adj.P.Val' or 'P.Value'.")
  }
  if (length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    tps_abort("alpha must be in (0,1).")
  }
  if (length(lfc_threshold) != 1L || is.na(lfc_threshold) ||
      !is.finite(lfc_threshold) || lfc_threshold < 0) {
    tps_abort("lfc_threshold must be one non-negative finite value.")
  }
  d <- tps_results(fit, contrast = contrast)
  d <- d[d$test == "t", , drop = FALSE]
  if (nrow(d) == 0L) tps_abort("contrast '%s' has no t-test rows.", contrast)
  ## significance column may be P.Value or adj.P.Val; pull by name
  sig_vals <- d[[significance]]
  d$direction <- ifelse(!is.na(sig_vals) & sig_vals < alpha &
                          abs(d$logFC) >= lfc_threshold,
                        ifelse(d$logFC > 0, "up", "down"), "ns")
  d$neglog10p <- -log10(pmax(sig_vals, 1e-300))
  cols <- .tps_colors()
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$logFC, .data$neglog10p,
                                        color = .data$direction)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed",
                        color = "grey50") +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                        linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = cols, name = NULL) +
    ggplot2::labs(title = contrast,
                  x = "log2 fold change", y = sprintf("-log10(%s)", significance)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right")
  .with_data(p, d)
}

## ===========================================================================
## tps_plot_trajectories
## ===========================================================================
#' Plot trajectory centroids or individual trajectories
#'
#' Produces a ggplot of cluster centroids (from
#' \code{\link{tps_cluster_trajectories}}) or individual feature trajectories.
#'
#' @param clusters a \code{tps_clusters} object (plots centroids), OR a
#'   \code{tps_trajectory} object (plots individual trajectories as thin lines).
#' @param n_features integer; if plotting individual trajectories, subsample
#'   to this many (default 100) for readability.
#' @return a ggplot object.
#'
#' @examples
#' data(tps_example, package = "timeProtState")
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' traj <- tps_trajectory_matrix(se)
#' cl <- tps_cluster_trajectories(traj, method = "kmeans", k = 4, seed = 1)
#' tps_plot_trajectories(cl)
#'
#' @export
tps_plot_trajectories <- function(clusters, n_features = 100L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    tps_abort("ggplot2 is required for plotting.")
  }
  if (!is.numeric(n_features) || length(n_features) != 1L ||
      is.na(n_features) || !is.finite(n_features) || n_features < 1L ||
      n_features != as.integer(n_features)) {
    tps_abort("n_features must be one positive integer.")
  }
  if (inherits(clusters, "tps_clusters")) {
    ## plot centroids: cells x clusters matrix; rownames = cells
    cen <- clusters$centroids  # matrix cells x clusters
    long <- reshape2_or_gather(cen, key = "cluster", value = "value",
                               id = "cell")
    long$cell <- factor(long$cell, levels = rownames(cen))
    long$cluster <- factor(long$cluster)
    p <- ggplot2::ggplot(long, ggplot2::aes(.data$cell, .data$value,
                                             group = .data$cluster,
                                             color = .data$cluster)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(title = "Cluster centroids",
                    x = NULL, y = "standardized value", color = "cluster") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    return(.with_data(p, long))
  }
  if (inherits(clusters, "tps_trajectory")) {
    m <- clusters$unique_cell
    if (nrow(m) > n_features) {
      ## deterministic subsample to avoid RNG perturbation; use first n_features
      m <- m[seq_len(n_features), , drop = FALSE]
    }
    long <- reshape2_or_gather(m, key = "cell", value = "value",
                               id = "feature_id")
    long$cell <- factor(long$cell, levels = colnames(m))
    p <- ggplot2::ggplot(long, ggplot2::aes(.data$cell, .data$value,
                                             group = .data$feature_id)) +
      ggplot2::geom_line(alpha = 0.15, color = "grey40") +
      ggplot2::labs(title = sprintf("Individual trajectories (n=%d)", nrow(m)),
                    x = NULL, y = "value") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    return(.with_data(p, long))
  }
  tps_abort("clusters must be a tps_clusters or tps_trajectory object.")
}

## lightweight wide->long without reshape2 dependency.
## `df` may be a matrix (rownames become the id) or a data.frame with an id col.
reshape2_or_gather <- function(df, key, value, id) {
  if (is.matrix(df)) {
    rn <- rownames(df)
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    df[[id]] <- rn
  }
  if (!id %in% names(df)) {
    df[[id]] <- seq_len(nrow(df))
  }
  cols <- setdiff(names(df), id)
  long <- do.call(rbind, lapply(cols, function(cc) {
    data.frame(id = df[[id]], key = cc, val = df[[cc]],
               stringsAsFactors = FALSE)
  }))
  names(long) <- c(id, key, value)
  long
}

## ===========================================================================
## tps_plot_pathway_states
## ===========================================================================
#' Plot pathway state distribution
#'
#' Produces a ggplot bar chart of pathway state counts (from
#' \code{\link{tps_classify_pathways}}) across contrasts.
#'
#' @param pathway_states a data.frame from \code{tps_classify_pathways}.
#' @return a ggplot object.
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
#' tps_plot_pathway_states(tps_classify_pathways(pws_c))
#' @export
tps_plot_pathway_states <- function(pathway_states) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    tps_abort("ggplot2 is required for plotting.")
  }
  if (!is.data.frame(pathway_states) || !"state" %in% names(pathway_states)) {
    tps_abort("pathway_states must be a data.frame from tps_classify_pathways.")
  }
  tab <- as.data.frame(table(pathway_states$contrast, pathway_states$state))
  names(tab) <- c("contrast","state","n")
  cols <- .tps_colors()
  p <- ggplot2::ggplot(tab, ggplot2::aes(.data$contrast, .data$n,
                                          fill = .data$state)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = cols, name = "state") +
    ggplot2::labs(title = "Pathway state distribution",
                  x = "contrast", y = "number of pathways") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  .with_data(p, tab)
}

## ===========================================================================
## tps_plot_state_distribution (protein-level)
## ===========================================================================
#' Plot protein state distribution
#'
#' Produces a ggplot bar chart of protein-level state counts (from
#' \code{\link{tps_classify_states}}).
#'
#' @param states a \code{tps_states} object.
#' @return a ggplot object.
#'
#' @examples
#' data(tps_example, package = "timeProtState")
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' design <- tps_make_design(se)
#' contrasts <- tps_make_contrasts(design, endpoint = 7)
#' fit <- tps_fit(se, design, contrasts)
#' st <- tps_classify_states(fit)
#' tps_plot_state_distribution(st)
#'
#' @export
tps_plot_state_distribution <- function(states) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    tps_abort("ggplot2 is required for plotting.")
  }
  if (!inherits(states, "tps_states")) {
    tps_abort("states must be a tps_states object.")
  }
  tab <- as.data.frame(table(states$state_table$state))
  names(tab) <- c("state","n")
  tab <- tab[order(-tab$n), ]
  p <- ggplot2::ggplot(tab, ggplot2::aes(stats::reorder(.data$state, -.data$n),
                                           .data$n)) +
    ggplot2::geom_col(fill = "#0072B2") +
    ggplot2::labs(title = "Protein dynamic state distribution",
                  x = NULL, y = "number of proteins") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  .with_data(p, tab)
}

## ===========================================================================
## tps_plot_pathway_trajectory
## ===========================================================================
#' Plot pathway remodeling trajectory over time
#'
#' Produces a ggplot showing NET, REM, or COH for a single pathway across all
#' timepoints and conditions, with significant timepoints marked using
#' q_rem for REM and q_net for NET/COH. Requires the pathway scores to have been computed with
#' \code{family = "vs_baseline"} (all timepoints), not just endpoint.
#'
#' @param calibrated a \code{tps_pathway_calibrated} object.
#' @param term character; the pathway term to plot.
#' @param metric character; \code{"NET"} (default), \code{"REM"}, or \code{"COH"}.
#' @param alpha numeric; q-value threshold for marking significant timepoints.
#' @return a ggplot object with underlying data attached as attribute \code{"data"}.
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
#' tps_plot_pathway_trajectory(pws_c, term = rownames(pws_c$scores$NET)[1])
#' @export
tps_plot_pathway_trajectory <- function(calibrated, term,
                                         metric = c("NET", "REM", "COH"),
                                         alpha = 0.05) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    tps_abort("ggplot2 is required for plotting.")
  }
  if (!inherits(calibrated, "tps_pathway_calibrated")) {
    tps_abort("calibrated must be a tps_pathway_calibrated object.")
  }
  if (length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    tps_abort("alpha must be in (0,1).")
  }
  metric <- match.arg(metric)
  if (!term %in% rownames(calibrated$scores[[metric]])) {
    tps_abort("term '%s' not found in pathway scores.", term)
  }
  ## build long data for this pathway
  meta <- .parse_contrast_time_pub(colnames(calibrated$scores[[metric]]))
  vals <- calibrated$scores[[metric]][term, ]
  q_name <- if (identical(metric, "REM")) "q_rem" else "q_net"
  qvals <- calibrated[[q_name]][term, ]
  df <- data.frame(
    time = meta$time, condition = meta$condition,
    value = as.numeric(vals),
    significant = !is.na(qvals) & qvals < alpha,
    stringsAsFactors = FALSE)
  df <- df[!is.na(df$time), , drop = FALSE]
  df <- df[order(df$condition, df$time), ]

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$time, .data$value,
                                         color = .data$condition,
                                         group = .data$condition)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$significant), size = 2.5) +
    ggplot2::scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
                                 name = sprintf("%s < %s", q_name, alpha)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
    ggplot2::labs(title = sprintf("%s -- %s trajectory", term, metric),
                  x = "time", y = metric, color = "condition") +
    ggplot2::theme_bw()
  .with_data(p, df)
}

## public wrapper for the contrast-time parser (used by plotting)
.parse_contrast_time_pub <- function(contrast_names) {
  base <- sub("_vs_baseline$", "", contrast_names)
  cond <- sub("_t[^_]+$", "", base)
  time_str <- sub("^.*_t", "", base)
  time <- suppressWarnings(as.numeric(time_str))
  data.frame(contrast = contrast_names, condition = cond, time = time,
             stringsAsFactors = FALSE)
}
