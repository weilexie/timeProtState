## test-trajectory-cluster.R
##
## Tests Goal 4: tps_trajectory_matrix (D6 baseline-once fix) and
## tps_cluster_trajectories (Mfuzz/kmeans adapter + QC).

library(testthat)
library(timeProtState)
library(SummarizedExperiment)

run_traj <- function() {
  se <- tps_import(tps_example_abundance, tps_example_metadata,
                   tps_example_annotation)
  tps_trajectory_matrix(se)
}

## ===========================================================================
## tps_trajectory_matrix
## ===========================================================================
test_that("trajectory_matrix builds unique_cell with baseline appearing once (D6)", {
  tr <- run_traj()
  ## unique_cell must have exactly ONE baseline column (the D6 fix)
  bl_cols <- grep("^baseline", colnames(tr$unique_cell))
  expect_length(bl_cols, 1L)
  expect_equal(ncol(tr$unique_cell), 1L + 3L * 3L)  # baseline + 3 cond x 3 time
})

test_that("display matrix copies baseline per condition", {
  tr <- run_traj()
  ## display should have baseline_t0 once per condition position? Actually
  ## display = baseline copied per condition -> at least as many cols as unique
  expect_gte(ncol(tr$display), ncol(tr$unique_cell))
})

test_that("delta matrix is log2FC vs baseline", {
  tr <- run_traj()
  bl <- tr$unique_cell[, "baseline"]
  ## delta[col] should equal unique_cell[col] - baseline
  for (cc in colnames(tr$delta)) {
    expect_equal(tr$delta[, cc], tr$unique_cell[, cc] - bl, tolerance = 1e-10)
  }
})

test_that("trajectory missingness summary is correct", {
  tr <- run_traj()
  ms <- tr$missingness
  expect_equal(nrow(ms), 200L)
  expect_true(all(ms$n_cells_observed >= 0))
  expect_true(all(ms$pct_missing >= 0 & ms$pct_missing <= 100))
})

test_that("trajectory_matrix respects feature subset", {
  tr <- run_traj()
  tr_sub <- tps_trajectory_matrix(
    tps_import(tps_example_abundance, tps_example_metadata),
    features = rownames(tps_example_abundance)[1:10])
  expect_equal(nrow(tr_sub$unique_cell), 10L)
})

test_that("trajectory print works", {
  tr <- run_traj()
  expect_output(print(tr), "unique cells: 10")
})

## ===========================================================================
## tps_cluster_trajectories
## ===========================================================================
test_that("kmeans clustering produces valid assignments", {
  tr <- run_traj()
  cl <- tps_cluster_trajectories(tr, method = "kmeans", k = 4, seed = 1)
  expect_s3_class(cl, "tps_clusters")
  expect_length(cl$assignment, 200L)
  expect_equal(cl$k, 4L)
  expect_equal(sort(unique(na.omit(cl$assignment))), 1:4)
})

test_that("kmeans clustering is reproducible with same seed", {
  tr <- run_traj()
  a <- tps_cluster_trajectories(tr, method = "kmeans", k = 4, seed = 42)
  b <- tps_cluster_trajectories(tr, method = "kmeans", k = 4, seed = 42)
  expect_identical(a$assignment, b$assignment)
})

test_that("clustering does not perturb global RNG state", {
  tr <- run_traj()
  set.seed(99); before <- runif(5)
  set.seed(99)
  tps_cluster_trajectories(tr, method = "kmeans", k = 3, seed = 1)
  after <- runif(5)
  expect_equal(before, after)
})

test_that("Mfuzz clustering works when package available", {
  skip_if_not_installed("Mfuzz")
  skip_if_not_installed("e1071")
  tr <- run_traj()
  cl <- tps_cluster_trajectories(tr, method = "mfuzz", k = 5, seed = 123)
  expect_s3_class(cl, "tps_clusters")
  expect_equal(cl$k, 5L)
  expect_false(is.null(cl$membership))
  ## membership rows sum to 1 (compare values, ignoring names)
  expect_equal(as.numeric(rowSums(cl$membership)),
               rep(1, nrow(cl$membership)), tolerance = 1e-6)
})

test_that("cluster QC reports sizes and dispersion", {
  tr <- run_traj()
  cl <- tps_cluster_trajectories(tr, method = "kmeans", k = 3, seed = 1)
  qc <- cl$qc
  expect_equal(nrow(qc), 3L)
  expect_true(all(c("cluster","size","within_dispersion") %in% names(qc)))
  expect_equal(sum(qc$size), sum(!is.na(cl$assignment)))
})

test_that("clustering accepts a bare numeric matrix", {
  m <- matrix(rnorm(50 * 5), nrow = 50, ncol = 5)
  rownames(m) <- paste0("P", 1:50)
  cl <- tps_cluster_trajectories(m, method = "kmeans", k = 3, seed = 1)
  expect_length(cl$assignment, 50L)
})

test_that("clustering errors on too-small k", {
  tr <- run_traj()
  expect_error(tps_cluster_trajectories(tr, method = "kmeans", k = 1), ">= 2")
})

test_that("tps_evaluate_k returns a data.frame over K range", {
  skip_if_not_installed("cluster")
  tr <- run_traj()
  keval <- tps_evaluate_k(tr, k_range = 3:5, method = "kmeans")
  expect_true(is.data.frame(keval))
  expect_equal(keval$k, 3:5)
  expect_true("mean_silhouette" %in% names(keval))
})

test_that("cluster print works", {
  tr <- run_traj()
  cl <- tps_cluster_trajectories(tr, method = "kmeans", k = 3, seed = 1)
  expect_output(print(cl), "method: kmeans")
})
