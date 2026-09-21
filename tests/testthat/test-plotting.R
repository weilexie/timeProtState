## test-plotting.R
##
## Plotting functions return a ggplot object
## (composable), carries extractable data, and does not re-fit.

library(testthat)
library(timeProtState)
library(SummarizedExperiment)

skip_if_no_ggplot2 <- function() {
  skip_if_not_installed("ggplot2")
}

run_pipeline <- function() {
  se <- tps_import(tps_example_abundance, tps_example_metadata,
                   tps_example_annotation)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)
  states <- tps_classify_states(fit)
  traj <- tps_trajectory_matrix(se)
  cl <- tps_cluster_trajectories(traj, method = "kmeans", k = 4, seed = 1)
  list(se = se, design = design, contrasts = contrasts, fit = fit,
       states = states, traj = traj, cl = cl)
}

## ===========================================================================
## tps_plot_volcano
## ===========================================================================
test_that("tps_plot_volcano returns a ggplot with data", {
  skip_if_no_ggplot2()
  p <- run_pipeline()
  plt <- tps_plot_volcano(p$fit, contrast = "endpoint_Acid_t7_vs_baseline")
  expect_true(ggplot2::is_ggplot(plt))
  expect_true(is.data.frame(attr(plt, "data")))
  expect_true("direction" %in% names(attr(plt, "data")))
})

test_that("tps_plot_volcano errors on unknown contrast", {
  skip_if_no_ggplot2()
  p <- run_pipeline()
  expect_error(tps_plot_volcano(p$fit, contrast = "nope"), "no t-test rows|contrast")
})

test_that("tps_plot_volcano respects significance and lfc thresholds", {
  skip_if_no_ggplot2()
  p <- run_pipeline()
  plt <- tps_plot_volcano(p$fit, contrast = "endpoint_Acid_t7_vs_baseline",
                          significance = "P.Value", alpha = 0.01,
                          lfc_threshold = 1)
  d <- attr(plt, "data")
  ## with lfc_threshold=1, 'up'/'down' must satisfy |logFC|>=1
  ups <- d[d$direction == "up", ]
  downs <- d[d$direction == "down", ]
  ## assert at least some points pass the threshold (non-vacuous test)
  expect_gt(nrow(ups) + nrow(downs), 0L)
  if (nrow(ups) > 0) expect_true(all(abs(ups$logFC) >= 1 - 1e-10))
  if (nrow(downs) > 0) expect_true(all(abs(downs$logFC) >= 1 - 1e-10))
})

## ===========================================================================
## tps_plot_trajectories
## ===========================================================================
test_that("tps_plot_trajectories plots centroids from tps_clusters", {
  skip_if_no_ggplot2()
  p <- run_pipeline()
  plt <- tps_plot_trajectories(p$cl)
  expect_true(ggplot2::is_ggplot(plt))
  d <- attr(plt, "data")
  expect_true("cluster" %in% names(d))
})

test_that("tps_plot_trajectories plots individual from tps_trajectory", {
  skip_if_no_ggplot2()
  p <- run_pipeline()
  plt <- tps_plot_trajectories(p$traj, n_features = 20)
  expect_true(ggplot2::is_ggplot(plt))
  d <- attr(plt, "data")
  expect_true("feature_id" %in% names(d))
  expect_lte(length(unique(d$feature_id)), 20L)
})

test_that("tps_plot_trajectories errors on bad input", {
  skip_if_no_ggplot2()
  expect_error(tps_plot_trajectories(123), "tps_clusters or tps_trajectory")
})

## ===========================================================================
## tps_plot_state_distribution
## ===========================================================================
test_that("tps_plot_state_distribution returns a ggplot", {
  skip_if_no_ggplot2()
  p <- run_pipeline()
  plt <- tps_plot_state_distribution(p$states)
  expect_true(ggplot2::is_ggplot(plt))
  d <- attr(plt, "data")
  expect_true(all(c("state","n") %in% names(d)))
})

test_that("tps_plot_state_distribution errors on bad input", {
  skip_if_no_ggplot2()
  expect_error(tps_plot_state_distribution(123), "tps_states")
})

## ===========================================================================
## tps_plot_pathway_states
## ===========================================================================
test_that("tps_plot_pathway_states returns a ggplot", {
  skip_if_no_ggplot2()
  p <- run_pipeline()
  pws <- tps_pathway_scores(p$fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  ave <- as.numeric(tapply(tps_results(p$fit,"endpoint")$AveExpr,
                            factor(tps_results(p$fit,"endpoint")$feature_id,
                                   levels = rownames(pws$statistic_matrix)),
                            function(x) x[1]))
  pws_c <- tps_calibrate_pathways(engine = "permutation", pws, abundance = ave, n_perm = 50L, seed = 1L)
  pwstates <- tps_classify_pathways(pws_c)
  plt <- tps_plot_pathway_states(pwstates)
  expect_true(ggplot2::is_ggplot(plt))
  d <- attr(plt, "data")
  expect_true(all(c("contrast","state","n") %in% names(d)))
})

test_that("tps_plot_pathway_states errors on bad input", {
  skip_if_no_ggplot2()
  expect_error(tps_plot_pathway_states(123), "tps_classify_pathways")
})

test_that("tps_plot_pathway_trajectory renders vs_baseline calibrated object", {
  skip_if_no_ggplot2()
  data(tps_example, package = "timeProtState")
  se <- tps_import(tps_example_abundance, tps_example_metadata,
                   tps_example_annotation)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                            id_col = "locus_tag", family = "vs_baseline")
  ave <- as.numeric(tapply(
    tps_results(fit, "vs_baseline")$AveExpr,
    factor(tps_results(fit, "vs_baseline")$feature_id,
           levels = rownames(pws$statistic_matrix)),
    function(x) x[1]))
  pws_c <- tps_calibrate_pathways(engine = "permutation", pws, abundance = ave, n_perm = 50, seed = 1)
  p <- tps_plot_pathway_trajectory(pws_c, term = rownames(pws_c$scores$NET)[1])
  expect_true(ggplot2::is_ggplot(p))
  expect_true(nrow(p$data) > 0)
  expect_true("condition" %in% names(p$data))
  expect_true("time" %in% names(p$data))
})

test_that("tps_plot_pathway_trajectory handles endpoint family prefix", {
  skip_if_no_ggplot2()
  data(tps_example, package = "timeProtState")
  se <- tps_import(tps_example_abundance, tps_example_metadata,
                   tps_example_annotation)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                            id_col = "locus_tag", family = "endpoint")
  ave <- as.numeric(tapply(
    tps_results(fit, "endpoint")$AveExpr,
    factor(tps_results(fit, "endpoint")$feature_id,
           levels = rownames(pws$statistic_matrix)),
    function(x) x[1]))
  pws_c <- tps_calibrate_pathways(engine = "permutation", pws, abundance = ave, n_perm = 50, seed = 1)
  p <- tps_plot_pathway_trajectory(pws_c, term = rownames(pws_c$scores$NET)[1])
  expect_true(ggplot2::is_ggplot(p))
  ## endpoint_ prefix must be stripped — condition should NOT start with "endpoint_"
  expect_false(any(grepl("^endpoint_", p$data$condition)))
})
