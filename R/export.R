## ============================================================================
## Export: write results to CSV/RDS/TSV with manifests
## ============================================================================
##
## Implements Goal 14 (partial -- export layer; report rendering is separate
## and optional). Provides structured export of fit/state/pathway/cluster
## result objects to a directory tree (spec section 7.2).

## ---- tps_export ------------------------------------------------------------
#' Export analysis results to a directory
#'
#' Writes the standard result tables from a completed analysis (fit, states,
#' trajectory, clusters, pathway scores/states) to a structured directory,
#' along with a manifest recording provenance and software versions.
#'
#' @param output_dir character; target directory (created if missing).
#' @param fit optional \code{tps_fit}.
#' @param states optional \code{tps_states}.
#' @param trajectory optional \code{tps_trajectory}.
#' @param clusters optional \code{tps_clusters}.
#' @param pathway_scores optional \code{tps_pathway_scores} or calibrated.
#' @param pathway_states optional pathway state data.frame.
#' @return invisibly, a data.frame of written files.
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
#' tps_export(tempdir(), fit = fit, states = tps_classify_states(fit))
#' @export
tps_export <- function(output_dir, fit = NULL, states = NULL,
                        trajectory = NULL, clusters = NULL,
                        pathway_scores = NULL, pathway_states = NULL) {
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    tps_abort("output_dir must be one non-empty path.")
  }
  ## Validate every supplied object before writing anything, so a late type
  ## error cannot leave a misleading partial export.
  if (!is.null(fit) && !inherits(fit, "tps_fit")) {
    tps_abort("fit must be a tps_fit.")
  }
  if (!is.null(states) && !inherits(states, "tps_states")) {
    tps_abort("states must be a tps_states.")
  }
  if (!is.null(trajectory) && !inherits(trajectory, "tps_trajectory")) {
    tps_abort("trajectory must be a tps_trajectory.")
  }
  if (!is.null(clusters) && !inherits(clusters, "tps_clusters")) {
    tps_abort("clusters must be a tps_clusters.")
  }
  if (!is.null(pathway_scores) &&
      !inherits(pathway_scores, "tps_pathway_scores") &&
      !inherits(pathway_scores, "tps_pathway_calibrated")) {
    tps_abort("pathway_scores must be a pathway score or calibrated object.")
  }
  if (!is.null(pathway_states) && !is.data.frame(pathway_states)) {
    tps_abort("pathway_states must be a data.frame.")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) tps_abort("Cannot create output_dir '%s'.", output_dir)
  export_state <- new.env(parent = emptyenv())
  export_state$written <- character(0)
  write_tbl <- function(obj, name) {
    if (is.null(obj)) return()
    path <- file.path(output_dir, name)
    if (inherits(obj, "data.frame")) {
      utils::write.csv(obj, path, row.names = FALSE)
    } else {
      saveRDS(obj, path)
    }
    export_state$written <- c(export_state$written, path)
  }
  matrix_table <- function(x, id_name) {
    if (is.null(x)) return(NULL)
    if (is.null(dim(x)) || length(dim(x)) != 2L) {
      tps_abort("export matrix '%s' is not two-dimensional.", id_name)
    }
    out <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
    ids <- rownames(x)
    if (is.null(ids)) ids <- as.character(seq_len(nrow(x)))
    out <- cbind(setNames(data.frame(ids, stringsAsFactors = FALSE), id_name),
                 out)
    rownames(out) <- NULL
    out
  }

  if (!is.null(fit)) {
    write_tbl(fit$statistics, "protein_statistics.csv")
    write_tbl(fit, "fit.rds")
  }
  if (!is.null(states)) {
    write_tbl(states$state_table, "protein_states.csv")
    write_tbl(states$rule_manifest, "state_rule_manifest.csv")
    write_tbl(states, "states.rds")
  }
  if (!is.null(trajectory)) {
    write_tbl(matrix_table(trajectory$unique_cell, "feature_id"),
              "trajectories_unique_cell.csv")
    write_tbl(matrix_table(trajectory$display, "feature_id"),
              "trajectories_display.csv")
    write_tbl(matrix_table(trajectory$delta, "feature_id"),
              "trajectories_delta.csv")
    write_tbl(trajectory$missingness, "trajectory_missingness.csv")
    write_tbl(trajectory, "trajectory.rds")
  }
  if (!is.null(clusters)) {
    assign_df <- data.frame(feature_id = names(clusters$assignment),
                             cluster = as.integer(clusters$assignment),
                             stringsAsFactors = FALSE)
    write_tbl(assign_df, "cluster_assignment.csv")
    write_tbl(clusters$qc, "cluster_qc.csv")
    write_tbl(matrix_table(clusters$centroids, "cell"),
              "cluster_centroids.csv")
    if (!is.null(clusters$membership)) {
      write_tbl(matrix_table(clusters$membership, "feature_id"),
                "cluster_membership.csv")
    }
    write_tbl(clusters, "clusters.rds")
  }
  if (!is.null(pathway_scores)) {
    write_tbl(pathway_scores$pathway_table, "pathway_audit.csv")
    if (!is.null(pathway_scores$scores$NET)) {
      for (nm in names(pathway_scores$scores)) {
        m <- pathway_scores$scores[[nm]]
        write_tbl(matrix_table(m, "term"), sprintf("pathway_%s.csv", nm))
      }
    }
    for (nm in c("p_net", "p_rem", "q_net", "q_rem",
                 "net_null_mean", "net_null_sd",
                 "rem_null_mean", "rem_null_sd",
                 "intergene_cor", "vif")) {
      m <- pathway_scores[[nm]]
      if (!is.null(m) && is.numeric(m) && length(m) > 0L) {
        write_tbl(matrix_table(matrix(m, ncol = 1, dimnames = list(names(m), nm)),
                               "term"), sprintf("pathway_%s.csv", nm))
      }
    }
    write_tbl(pathway_scores, "pathway_scores.rds")
  }
  if (!is.null(pathway_states)) {
    write_tbl(pathway_states, "pathway_states.csv")
  }

  ## manifest (always written, even with no result inputs)
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  pv <- as.character(getNamespaceVersion("timeProtState"))
  written <- export_state$written
  if (length(written) == 0L) {
    manifest <- data.frame(file = character(0), path = character(0),
                            written_at = character(0),
                            package_version = character(0),
                            stringsAsFactors = FALSE)
  } else {
    manifest <- data.frame(
      file = basename(written),
      path = written,
      written_at = ts,
      package_version = pv,
      stringsAsFactors = FALSE)
  }
  utils::write.csv(manifest, file.path(output_dir, "export_manifest.csv"),
                   row.names = FALSE)
  invisible(manifest)
}
