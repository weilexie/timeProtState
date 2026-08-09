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
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) tps_abort("Cannot create output_dir '%s'.", output_dir)
  written <- character(0)
  write_tbl <- function(obj, name) {
    if (is.null(obj)) return()
    path <- file.path(output_dir, name)
    if (inherits(obj, "data.frame")) {
      utils::write.csv(obj, path, row.names = FALSE)
    } else {
      saveRDS(obj, path)
    }
    written <<- c(written, path)
  }

  if (!is.null(fit)) {
    if (!inherits(fit, "tps_fit")) tps_abort("fit must be a tps_fit.")
    write_tbl(fit$statistics, "protein_statistics.csv")
    saveRDS(fit, file.path(output_dir, "fit.rds"))
    written <<- c(written, file.path(output_dir, "fit.rds"))
  }
  if (!is.null(states)) {
    if (!inherits(states, "tps_states")) tps_abort("states must be a tps_states.")
    write_tbl(states$state_table, "protein_states.csv")
    write_tbl(states$rule_manifest, "state_rule_manifest.csv")
  }
  if (!is.null(trajectory)) {
    if (!inherits(trajectory, "tps_trajectory")) {
      tps_abort("trajectory must be a tps_trajectory.")
    }
    write_tbl(as.data.frame(trajectory$delta), "trajectories_delta.csv")
    write_tbl(trajectory$missingness, "trajectory_missingness.csv")
  }
  if (!is.null(clusters)) {
    if (!inherits(clusters, "tps_clusters")) {
      tps_abort("clusters must be a tps_clusters.")
    }
    assign_df <- data.frame(feature_id = names(clusters$assignment),
                             cluster = as.integer(clusters$assignment),
                             stringsAsFactors = FALSE)
    write_tbl(assign_df, "cluster_assignment.csv")
    write_tbl(clusters$qc, "cluster_qc.csv")
  }
  if (!is.null(pathway_scores)) {
    if (inherits(pathway_scores, "tps_pathway_scores") ||
        inherits(pathway_scores, "tps_pathway_calibrated")) {
      write_tbl(pathway_scores$pathway_table, "pathway_audit.csv")
      if (!is.null(pathway_scores$scores$NET)) {
        for (nm in names(pathway_scores$scores)) {
          m <- pathway_scores$scores[[nm]]
          write_tbl(as.data.frame(m), sprintf("pathway_%s.csv", nm))
        }
      }
    }
  }
  if (!is.null(pathway_states)) {
    write_tbl(pathway_states, "pathway_states.csv")
  }

  ## manifest (always written, even with no result inputs)
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  pv <- as.character(getNamespaceVersion("timeProtState"))
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
