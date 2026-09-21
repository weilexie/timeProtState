## Tests for tps_simulate_correlated_pathways: within-pathway correlation
## generation, signal injection, determinism, and pipeline integration.
##
## Validates the simulation tool that exposes the core calibration risk
## flagged in STATUS section 7: the abundance-stratified permutation null
## breaks within-pathway protein correlation.

## ---- structure -------------------------------------------------------------
test_that("tps_simulate_correlated_pathways returns expected structure", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 5, pathway_size = 8, n_background = 40,
    conditions = c("A", "B"), times = c(1, 3, 7), n_replicates = 3,
    within_pathway_cor = 0.5, pathway_signal = "codirectional", seed = 1)
  expect_s3_class(sim, "tps_simulation_correlated")
  expect_s3_class(sim, "tps_simulation")
  ## 5*8 + 40 = 80 features; baseline(3) + 2*3*3 = 21 samples
  expect_equal(dim(sim$abundance), c(80, 21))
  expect_equal(nrow(sim$col_data), 21)
  expect_true(all(c("term", "gene") %in% names(sim$term2gene)))
  expect_equal(nrow(sim$term2gene), 40)  # 5 pathways * 8 members
  ## baseline appears once (shared), not per-condition
  expect_equal(sum(sim$col_data$is_baseline), 3)
  ## feature ids are P%04d, matchable to term2gene$gene
  expect_true(all(sim$term2gene$gene %in% rownames(sim$abundance)))
  ## truth has one row per feature per condition
  expect_equal(nrow(sim$truth), 80 * 2)
})

test_that("pathway membership is disjoint and background has no pathway", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 3, pathway_size = 5, n_background = 10,
    conditions = c("A"), times = c(1, 3), n_replicates = 2, seed = 1)
  ## each protein belongs to at most one pathway
  expect_equal(anyDuplicated(sim$term2gene$gene), 0)
  ## background proteins (P0016..P0025) not in term2gene
  bg_ids <- sprintf("P%04d", 16:25)
  expect_false(any(bg_ids %in% sim$term2gene$gene))
})

## ---- correlation behavior --------------------------------------------------
test_that("within_pathway_cor=0.5 induces positive within-pathway correlation", {
  ## Large enough for stable estimate; center per-protein to remove mu0 offset.
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 1, pathway_size = 30, n_background = 30,
    conditions = c("A", "B"), times = c(1, 3, 7, 10, 15), n_replicates = 10,
    within_pathway_cor = 0.5, pathway_signal = "none", seed = 7)
  pw <- sim$abundance[1:30, , drop = FALSE]
  bg <- sim$abundance[31:60, , drop = FALSE]
  pw_c <- sweep(pw, 1, rowMeans(pw), "-")
  bg_c <- sweep(bg, 1, rowMeans(bg), "-")
  pw_cor <- mean(cor(t(pw_c))[lower.tri(diag(30))])
  bg_cor <- mean(cor(t(bg_c))[lower.tri(diag(30))])
  ## within-pathway correlation should be near 0.5 (tolerance for MC noise)
  expect_gt(pw_cor, 0.35)
  expect_lt(pw_cor, 0.65)
  ## background proteins remain ~independent
  expect_lt(abs(bg_cor), 0.15)
})

test_that("within_pathway_cor=0 gives independent noise (backwards compatible)", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 1, pathway_size = 30, n_background = 30,
    conditions = c("A", "B"), times = c(1, 3, 7, 10, 15), n_replicates = 10,
    within_pathway_cor = 0, pathway_signal = "none", seed = 7)
  pw <- sim$abundance[1:30, , drop = FALSE]
  pw_c <- sweep(pw, 1, rowMeans(pw), "-")
  pw_cor <- mean(cor(t(pw_c))[lower.tri(diag(30))])
  expect_lt(abs(pw_cor), 0.15)
})

test_that("background proteins are independent under correlation", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 2, pathway_size = 10, n_background = 30,
    conditions = c("A", "B"), times = c(1, 3, 7, 10, 15), n_replicates = 10,
    within_pathway_cor = 0.5, pathway_signal = "none", seed = 11)
  bg <- sim$abundance[21:50, , drop = FALSE]
  bg_c <- sweep(bg, 1, rowMeans(bg), "-")
  bg_cor <- mean(cor(t(bg_c))[lower.tri(diag(30))])
  expect_lt(abs(bg_cor), 0.15)
})

## ---- signal injection ------------------------------------------------------
test_that("pathway_signal='none' marks no signal members", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 4, pathway_size = 5, n_background = 10,
    conditions = c("A", "B"), times = c(1, 3, 7), n_replicates = 2,
    within_pathway_cor = 0.3, pathway_signal = "none", seed = 1)
  expect_false(any(sim$truth$is_signal_member))
  expect_equal(length(sim$config$signal_pathways), 0)
})

test_that("pathway_signal='codirectional' marks half of pathways as signal", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 4, pathway_size = 5, n_background = 10,
    conditions = c("A", "B"), times = c(1, 3, 7), n_replicates = 2,
    within_pathway_cor = 0.3, pathway_signal = "codirectional", seed = 1)
  ## signal_fraction default 0.5 -> 2 of 4 pathways
  expect_equal(length(sim$config$signal_pathways), 2)
  ## signal members exist and are in signal pathways
  sig_truth <- sim$truth[sim$truth$is_signal_member, ]
  expect_true(all(sig_truth$pathway %in% sim$config$signal_pathways))
  ## signal members span both conditions
  expect_setequal(unique(sig_truth$condition), c("A", "B"))
})

test_that("pathway_signal='mixed' creates opposing signs within pathway", {
  ## With noise_sd=0 and a single condition/timepoint, the mixed signal makes
  ## half the pathway members +amplitude and half -amplitude at the endpoint.
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 1, pathway_size = 10, n_background = 0,
    conditions = c("A"), times = c(1, 7), n_replicates = 1,
    within_pathway_cor = 0, pathway_signal = "mixed",
    signal_fraction = 1, amplitude = 3, noise_sd = 0, seed = 1)
  ## At endpoint (t=7, last timepoint), ramp = amplitude (full).
  ## Members 1-5 get +3, members 6-10 get -3.
  endpoint_col <- grep("t7", colnames(sim$abundance), value = TRUE)
  vals <- sim$abundance[, endpoint_col]
  ## member_index stored in truth; align by feature_id
  tr <- sim$truth[sim$truth$condition == "A", ]
  tr <- tr[match(rownames(sim$abundance), tr$feature_id), ]
  ## subtract per-protein baseline mean to isolate the effect
  base_col <- grep("^baseline", colnames(sim$abundance))
  eff <- vals - sim$abundance[, base_col]
  up_idx <- which(tr$member_index <= 5)
  down_idx <- which(tr$member_index > 5)
  expect_true(all(eff[up_idx] > 2.5))
  expect_true(all(eff[down_idx] < -2.5))
})

test_that("signal_fraction=0 disables signal even with codirectional type", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 4, pathway_size = 5, n_background = 10,
    conditions = c("A"), times = c(1, 3), n_replicates = 2,
    within_pathway_cor = 0.3, pathway_signal = "codirectional",
    signal_fraction = 0, seed = 1)
  expect_false(any(sim$truth$is_signal_member))
})

## ---- determinism -----------------------------------------------------------
test_that("same seed reproduces abundance and truth exactly", {
  a <- tps_simulate_correlated_pathways(
    n_pathways = 3, pathway_size = 6, n_background = 15,
    conditions = c("A", "B"), times = c(1, 3, 7), n_replicates = 2,
    within_pathway_cor = 0.4, pathway_signal = "mixed", seed = 99)
  b <- tps_simulate_correlated_pathways(
    n_pathways = 3, pathway_size = 6, n_background = 15,
    conditions = c("A", "B"), times = c(1, 3, 7), n_replicates = 2,
    within_pathway_cor = 0.4, pathway_signal = "mixed", seed = 99)
  expect_identical(a$abundance, b$abundance)
  expect_identical(a$truth, b$truth)
  expect_identical(a$term2gene, b$term2gene)
})

test_that("different seeds give different abundance", {
  a <- tps_simulate_correlated_pathways(
    n_pathways = 3, pathway_size = 6, n_background = 15,
    conditions = c("A"), times = c(1, 3), n_replicates = 2,
    within_pathway_cor = 0.4, pathway_signal = "none", seed = 1)
  b <- tps_simulate_correlated_pathways(
    n_pathways = 3, pathway_size = 6, n_background = 15,
    conditions = c("A"), times = c(1, 3), n_replicates = 2,
    within_pathway_cor = 0.4, pathway_signal = "none", seed = 2)
  expect_false(identical(a$abundance, b$abundance))
})

## ---- argument validation ---------------------------------------------------
test_that("invalid within_pathway_cor errors", {
  expect_error(tps_simulate_correlated_pathways(within_pathway_cor = -0.1),
               "within_pathway_cor")
  expect_error(tps_simulate_correlated_pathways(within_pathway_cor = 1),
               "within_pathway_cor")
  expect_error(tps_simulate_correlated_pathways(within_pathway_cor = 1.5),
               "within_pathway_cor")
})

test_that("invalid pathway_signal errors", {
  expect_error(tps_simulate_correlated_pathways(pathway_signal = "bogus"))
})

test_that("pathway_size >= 2 required", {
  expect_error(tps_simulate_correlated_pathways(pathway_size = 1),
               "pathway_size")
})

test_that("signal_fraction out of range errors", {
  expect_error(tps_simulate_correlated_pathways(signal_fraction = -0.1),
               "signal_fraction")
  expect_error(tps_simulate_correlated_pathways(signal_fraction = 1.5),
               "signal_fraction")
})

## ---- pipeline integration smoke --------------------------------------------
test_that("simulation feeds tps_pathway_scores end-to-end", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 5, pathway_size = 8, n_background = 40,
    conditions = c("A", "B"), times = c(1, 3, 7), n_replicates = 3,
    within_pathway_cor = 0.3, pathway_signal = "codirectional", seed = 1)
  se <- tps_import(sim$abundance, sim$col_data)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)
  pws <- tps_pathway_scores(fit, sim$term2gene, family = "vs_baseline")
  ave <- rep(1, nrow(pws$statistic_matrix))
  pws_c <- tps_calibrate_pathways(engine = "permutation", pws, abundance = ave, n_perm = 50, seed = 1)
  expect_s3_class(pws_c, "tps_pathway_calibrated")
  expect_equal(nrow(pws_c$q_net), 5)
})

## ---- missingness -----------------------------------------------------------
test_that("missingness injects NA at approximately the right fraction", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 2, pathway_size = 6, n_background = 12,
    conditions = c("A"), times = c(1, 3), n_replicates = 2,
    within_pathway_cor = 0.3, pathway_signal = "none",
    missingness = 0.2, seed = 1)
  frac <- mean(is.na(sim$abundance))
  expect_gt(frac, 0.1)
  expect_lt(frac, 0.3)
})
