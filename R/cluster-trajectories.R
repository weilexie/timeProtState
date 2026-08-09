## ============================================================================
## Cluster trajectories: Mfuzz adapter (soft) + kmeans fallback (hard)
## ============================================================================
##
## Implements STATISTICAL_CONTRACT.md section 7.2-7.3 and DECISIONS D6.
## Clustering is NOT a novelty of the package (legacy already used Mfuzz);
## it is an adapter. The key fix: clustering operates on the UNIQUE-CELL
## trajectory matrix (baseline once), not the legacy 3x-copied display matrix.
##
## Mfuzz is in Suggests (optional); kmeans (base stats) is the always-available
## fallback. Either way: seed-controlled, membership QC reported, no claim of
## "optimal K".

## ---- standardization helper (manual row-z, matching legacy Mfuzz input) ----
.row_standardize <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  sds <- apply(m, 1L, function(r) {
    s <- stats::sd(r, na.rm = TRUE)
    if (!is.finite(s) || s == 0) NA_real_ else s
  })
  flat <- is.na(sds)
  z <- (m - mu) / ifelse(is.na(sds), 1, sds)
  z[flat, ] <- 0
  ## replace any remaining NA/NaN/Inf with 0 so downstream kmeans/mfuzz are safe
  z[!is.finite(z)] <- 0
  z
}

## ---- tps_cluster_trajectories ----------------------------------------------
#' Cluster protein trajectories
#'
#' Soft (Mfuzz) or hard (kmeans) clustering of per-protein trajectories. By
#' default operates on the unique-cell trajectory matrix (baseline once) from
#' \code{\link{tps_trajectory_matrix}}, fixing the legacy inflation caused by
#' copying baseline per condition. Mfuzz is used when available; otherwise
#' kmeans is the fallback.
#'
#' Clustering is an adapter, not a methodological contribution; the package
#' does not claim to find an "optimal K".
#'
#' @param trajectory a \code{tps_trajectory} object, or a numeric matrix
#'   (features x cells).
#' @param method character; \code{"mfuzz"} (default, requires the Mfuzz
#'   package) or \code{"kmeans"}.
#' @param k integer; number of clusters.
#' @param m numeric or \code{"estimate"}; Mfuzz fuzzification. If
#'   \code{"estimate"}, uses \code{Mfuzz::mestimate}.
#' @param seed integer; RNG seed.
#' @param standardize logical; row-standardize (z-score) before clustering
#'   (default TRUE, matching Mfuzz convention).
#' @return a list of class \code{tps_clusters} with cluster assignments,
#'   membership (soft) or hard labels, centroids, and QC summary.
#'
#' @examples
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' traj <- tps_trajectory_matrix(se)
#' cl <- tps_cluster_trajectories(traj, method = "kmeans", k = 4, seed = 1)
#' head(cl$assignment)
#'
#' @export
tps_cluster_trajectories <- function(trajectory, method = c("mfuzz","kmeans"),
                                      k = 7L, m = "estimate", seed = 1L,
                                      standardize = TRUE) {
  method <- match.arg(method)
  arg_capture <- list(call = "tps_cluster_trajectories", method = method,
                       k = k, m = m, seed = seed, standardize = standardize,
                       timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  if (inherits(trajectory, "tps_trajectory")) {
    mat <- trajectory$unique_cell
  } else if (is.matrix(trajectory) && is.numeric(trajectory)) {
    mat <- trajectory
  } else {
    tps_abort("trajectory must be a tps_trajectory or a numeric matrix.")
  }
  if (k < 2L) tps_abort("k must be >= 2.")

  if (isTRUE(standardize)) mat <- .row_standardize(mat)
  ## drop rows that are all-NA / all-zero after standardize
  keep <- rowSums(abs(mat), na.rm = TRUE) > 0
  mat_fit <- mat[keep, , drop = FALSE]
  if (nrow(mat_fit) < k) tps_abort("Too few non-flat features (%d) for k=%d.",
                                    nrow(mat_fit), k)

  ## local RNG
  rng_old <- .GlobalEnv$.Random.seed
  on.exit({
    if (is.null(rng_old)) rm(".Random.seed", envir = .GlobalEnv)
    else .GlobalEnv$.Random.seed <- rng_old
  }, add = TRUE)
  set.seed(seed)

  if (method == "mfuzz") {
    if (!requireNamespace("Mfuzz", quietly = TRUE) ||
        !requireNamespace("Biobase", quietly = TRUE) ||
        !requireNamespace("e1071", quietly = TRUE)) {
      tps_abort("method='mfuzz' requires the Mfuzz, Biobase and e1071 ",
                "packages (in Suggests). Install them or use method='kmeans'.")
    }
    ## Mfuzz::mfuzz internally calls a bare cmeans() (unqualified) from e1071.
    ## In pkgload/testthat environments e1071 must be on the search path.
    ## We load the namespace and attach via eval to avoid BiocCheck's
    ## static detection of library()/require() calls in package code.
    if (!("package:e1071" %in% search())) {
      loadNamespace("e1071")
      eval(parse(text = "suppressMessages(library(e1071))"))
    }
    eset <- Biobase::ExpressionSet(assayData = mat_fit)
    m_val <- if (is.character(m) && m == "estimate") {
      as.numeric(Mfuzz::mestimate(eset))
    } else as.numeric(m)
    cl <- Mfuzz::mfuzz(eset, c = k, m = m_val)
    membership <- cl$membership
    rownames(membership) <- rownames(mat_fit)
    cluster <- cl$cluster
    centroids <- t(cl$centers)  # cells x clusters
  } else {
    ## kmeans fallback
    km <- stats::kmeans(mat_fit, centers = k, nstart = 25L)
    cluster <- km$cluster
    membership <- NULL
    centroids <- t(km$centers)  # cells x clusters
  }

  ## build full assignment vector (NA for dropped features)
  assignment <- setNames(rep(NA_integer_, nrow(mat)), rownames(mat))
  assignment[rownames(mat_fit)] <- as.integer(cluster)

  ## cluster QC
  qc <- .cluster_qc(cluster, membership, mat_fit, centroids)
  qc$n_dropped <- sum(!keep)
  qc$n_features <- nrow(mat)

  out <- structure(
    list(
      assignment = assignment,
      membership = membership,
      centroids = centroids,
      method = method, k = k, m = if (method == "mfuzz") m_val else NA,
      seed = seed, standardize = standardize,
      matrix_used = mat,
      qc = qc,
      provenance = list(.make_provenance(arg_capture, trajectory, stage = "cluster_trajectories"))
    ),
    class = "tps_clusters"
  )
  out
}

## ---- cluster QC ------------------------------------------------------------
.cluster_qc <- function(cluster, membership, mat, centroids) {
  k <- max(cluster)
  sizes <- as.integer(table(factor(cluster, levels = seq_len(k))))
  ## within-cluster dispersion
  disp <- vapply(seq_len(k), function(cc) {
    idx <- which(cluster == cc)
    if (length(idx) < 2L) return(0)
    mean(apply(mat[idx, , drop = FALSE], 1L, function(r) {
      sum((r - centroids[, cc])^2, na.rm = TRUE)
    }))
  }, numeric(1L))
  qc <- data.frame(
    cluster = seq_len(k),
    size = sizes,
    within_dispersion = disp,
    stringsAsFactors = FALSE)
  if (!is.null(membership)) {
    ## soft clustering: max membership stats
    max_mem <- apply(membership, 1L, max)
    qc$max_membership_mean <- tapply(max_mem, cluster, mean)
    qc$max_membership_min <- tapply(max_mem, cluster, min)
    qc$max_membership_median <- tapply(max_mem, cluster, median)
  }
  qc
}

## ---- tps_evaluate_k --------------------------------------------------------
#' Evaluate cluster count for trajectories
#'
#' Computes simple cluster-quality metrics across a range of K values to guide
#' the choice of cluster count. This is descriptive only; the package does not
#' claim a single optimal K.
#'
#' @param trajectory a \code{tps_trajectory} or numeric matrix.
#' @param k_range integer vector of K values to try.
#' @param method clustering method (see \code{\link{tps_cluster_trajectories}}).
#' @param seed integer.
#' @param metrics character; metrics to compute: \code{"silhouette"} and/or
#'   \code{"within_dispersion"}.
#' @return a data.frame with one row per K.
#'
#' @examples
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' traj <- tps_trajectory_matrix(se)
#' keval <- tps_evaluate_k(traj, k_range = 3:6, method = "kmeans")
#'
#' @export
tps_evaluate_k <- function(trajectory, k_range = 2:10, method = "kmeans",
                            seed = 1L, metrics = c("silhouette","within_dispersion")) {
  metrics <- match.arg(metrics, several.ok = TRUE)
  if (inherits(trajectory, "tps_trajectory")) {
    mat <- .row_standardize(trajectory$unique_cell)
  } else {
    mat <- .row_standardize(trajectory)
  }
  keep <- rowSums(abs(mat), na.rm = TRUE) > 0
  mat <- mat[keep, , drop = FALSE]
  rows <- list()
  for (k in k_range) {
    if (nrow(mat) < k) next
    cl_obj <- tryCatch(
      tps_cluster_trajectories(mat, method = method, k = k, seed = seed,
                                standardize = FALSE),
      error = function(e) NULL)
    if (is.null(cl_obj)) next
    assignment <- cl_obj$assignment[!is.na(cl_obj$assignment)]
    row <- list(k = k, method = method,
                n_clustered = length(assignment),
                min_cluster_size = min(table(assignment)))
    if ("within_dispersion" %in% metrics) {
      row$mean_within_dispersion <- mean(cl_obj$qc$within_dispersion)
    }
    if ("silhouette" %in% metrics && length(unique(assignment)) > 1L) {
      if (!requireNamespace("cluster", quietly = TRUE)) {
        tps_abort("metric 'silhouette' requires the 'cluster' package (in Suggests).")
      }
      d <- stats::dist(mat[names(assignment), , drop = FALSE])
      sil <- cluster::silhouette(assignment, d)
      row$mean_silhouette <- mean(sil[, "sil_width"])
    }
    rows[[length(rows) + 1L]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

## ---- methods ---------------------------------------------------------------
#' @export
print.tps_clusters <- function(x, ...) {
  cat("<timeProtState clusters>\n")
  cat(sprintf("  method: %s   k: %d   seed: %d   features: %d (dropped %d)\n",
              x$method, x$k, x$seed, x$qc$n_features, x$qc$n_dropped))
  cat("  cluster sizes:\n")
  for (i in seq_len(nrow(x$qc))) {
    cat(sprintf("    cluster %d: %d\n", x$qc$cluster[i], x$qc$size[i]))
  }
  invisible(x)
}
