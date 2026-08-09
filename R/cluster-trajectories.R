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
#' data(tps_example, package = "timeProtState")
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
  if (nrow(mat) == 0L || ncol(mat) < 2L) {
    tps_abort("trajectory matrix must have features and at least two cells.")
  }
  if (is.null(rownames(mat))) {
    rownames(mat) <- sprintf("feature_%d", seq_len(nrow(mat)))
  }
  if (anyNA(rownames(mat)) || any(!nzchar(rownames(mat))) ||
      anyDuplicated(rownames(mat))) {
    tps_abort("trajectory matrix row names must be unique, non-empty feature ids.")
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- sprintf("cell_%d", seq_len(ncol(mat)))
  }
  if (any(is.infinite(mat))) {
    tps_abort("trajectory matrix must not contain infinite values.")
  }
  if (length(k) != 1L || is.na(k) || k < 2L || k != as.integer(k)) {
    tps_abort("k must be one integer >= 2.")
  }
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
    tps_abort("seed must be one finite value.")
  }
  if (length(standardize) != 1L || is.na(standardize)) {
    tps_abort("standardize must be TRUE or FALSE.")
  }

  if (isTRUE(standardize)) {
    mat <- .row_standardize(mat)
  } else if (any(!is.finite(mat))) {
    tps_abort("standardize=FALSE requires a fully finite trajectory matrix.")
  }
  ## drop rows that are all-NA / all-zero after standardize
  keep <- rowSums(abs(mat), na.rm = TRUE) > 0
  mat_fit <- mat[keep, , drop = FALSE]
  if (nrow(mat_fit) < k) tps_abort("Too few non-flat features (%d) for k=%d.",
                                    nrow(mat_fit), k)

  withr::with_seed(seed, {
    if (method == "mfuzz") {
    if (!requireNamespace("Mfuzz", quietly = TRUE) ||
        !requireNamespace("Biobase", quietly = TRUE) ||
        !requireNamespace("e1071", quietly = TRUE)) {
      tps_abort(paste0(
        "method='mfuzz' requires the Mfuzz, Biobase and e1071 ",
        "packages (in Suggests). Install them or use method='kmeans'."
      ))
    }
    eset <- Biobase::ExpressionSet(assayData = mat_fit)
    if (!(is.character(m) && length(m) == 1L && identical(m, "estimate")) &&
        !(is.numeric(m) && length(m) == 1L && !is.na(m) &&
          is.finite(m) && m > 1)) {
      tps_abort("m must be 'estimate' or one finite value greater than 1.")
    }
    m_val <- if (is.character(m)) {
      as.numeric(Mfuzz::mestimate(eset))
    } else {
      as.numeric(m)
    }
    if (length(m_val) != 1L || is.na(m_val) || !is.finite(m_val) || m_val <= 1) {
      tps_abort("m must be 'estimate' or one finite value greater than 1.")
    }
    ## Mfuzz::mfuzz is a thin wrapper around e1071::cmeans but relies on
    ## e1071 being attached. Call the underlying routine explicitly to avoid
    ## changing the caller's search path.
    cl <- e1071::cmeans(Biobase::exprs(eset), centers = k,
                        method = "cmeans", m = m_val)
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
  })
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
#' data(tps_example, package = "timeProtState")
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' traj <- tps_trajectory_matrix(se)
#' keval <- tps_evaluate_k(traj, k_range = 3:6, method = "kmeans")
#'
#' @export
tps_evaluate_k <- function(trajectory, k_range = 2:10, method = "kmeans",
                            seed = 1L, metrics = c("silhouette","within_dispersion")) {
  method <- match.arg(method, c("mfuzz", "kmeans"))
  metrics <- match.arg(metrics, several.ok = TRUE)
  if (inherits(trajectory, "tps_trajectory")) {
    mat <- .row_standardize(trajectory$unique_cell)
  } else if (is.matrix(trajectory) && is.numeric(trajectory)) {
    mat <- .row_standardize(trajectory)
  } else {
    tps_abort("trajectory must be a tps_trajectory or a numeric matrix.")
  }
  if (!is.numeric(k_range) || length(k_range) == 0L || anyNA(k_range) ||
      any(!is.finite(k_range)) || any(k_range < 2) ||
      any(k_range != as.integer(k_range))) {
    tps_abort("k_range must contain integers >= 2.")
  }
  k_range <- unique(as.integer(k_range))
  keep <- rowSums(abs(mat), na.rm = TRUE) > 0
  mat <- mat[keep, , drop = FALSE]
  rows <- list()
  failures <- character(0)
  for (k in k_range) {
    if (nrow(mat) < k) next
    attempt <- tryCatch(
      list(
        value = tps_cluster_trajectories(
          mat, method = method, k = k, seed = seed, standardize = FALSE
        ),
        error = NULL
      ),
      error = function(e) list(value = NULL, error = conditionMessage(e))
    )
    if (!is.null(attempt$error)) {
      failures <- c(failures, sprintf("k=%d: %s", k, attempt$error))
      next
    }
    cl_obj <- attempt$value
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
  if (length(rows) == 0L) {
    detail <- if (length(failures)) utils::tail(failures, 1L) else
      "all K exceed feature count"
    tps_abort("No K evaluation succeeded (%s).", detail)
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
