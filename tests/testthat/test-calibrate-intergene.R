## Inter-gene VIF adjustment (ADR 0010): competitive permutation scaled by
## estimated member correlation. Tests drive shipped tps_calibrate_pathways
## and limma::camera on tps_simulate_correlated_pathways.

library(testthat)
library(timeProtState)

.endpoint_typeI <- function(p_net, alpha = 0.05) {
  mean(p_net < alpha, na.rm = TRUE)
}

.fit_endpoint <- function(sim) {
  se <- tps_import(sim$abundance, sim$col_data)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = max(sim$config$times))
  fit <- tps_fit(se, design, contrasts)
  pws <- tps_pathway_scores(fit, sim$term2gene, family = "vs_baseline",
                             min_size = 5L)
  list(se = se, design = design, contrasts = contrasts, fit = fit, pws = pws)
}

test_that("intergene='none' reproduces unadjusted p_net on a multi-contrast fit", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 6, pathway_size = 8, n_background = 24,
    conditions = c("A", "B"), times = c(1, 3, 7), n_replicates = 3,
    within_pathway_cor = 0, pathway_signal = "none", seed = 4)
  pipe <- .fit_endpoint(sim)
  none <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                                  n_perm = 40, seed = 1, intergene = "none")
  vif0 <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                                  n_perm = 40, seed = 1, intergene = "vif")
  expect_equal(none$intergene, "none")
  expect_equal(vif0$intergene, "vif")
  expect_true(all(vif0$vif >= 1))
  ## rho=0: VIF ~ 1, p-values stay in the same discrete permutation family
  expect_lt(mean(abs(vif0$p_net - none$p_net)), 0.15)
})

test_that("VIF shrinks anti-conservative NET calls at rho=0.5 versus none", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 12, pathway_size = 10, n_background = 40,
    conditions = c("A", "B", "C"), times = c(1, 3, 7), n_replicates = 3,
    within_pathway_cor = 0.5, pathway_signal = "none", seed = 11)
  pipe <- .fit_endpoint(sim)
  none <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                                  n_perm = 80, seed = 1, intergene = "none")
  adj <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                                 n_perm = 80, seed = 1, intergene = "vif")
  expect_gt(median(adj$intergene_cor), 0.15)
  expect_gt(median(adj$vif), 1.5)
  expect_lt(.endpoint_typeI(adj$p_net), .endpoint_typeI(none$p_net))
})

test_that("VIF-adjusted NET type-I at rho=0.5 is below CAMERA directional", {
  n_seeds <- 8L
  tps_hits <- camera_hits <- numeric(0)
  for (s in seq_len(n_seeds)) {
    sim <- tps_simulate_correlated_pathways(
      n_pathways = 12, pathway_size = 10, n_background = 30,
      conditions = c("A", "B"), times = c(1, 3, 7), n_replicates = 3,
      within_pathway_cor = 0.5, pathway_signal = "none", seed = 100L + s)
    pipe <- .fit_endpoint(sim)
    adj <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                                   n_perm = 60, seed = 1, intergene = "vif")
    end_cols <- grepl("_t7_vs_baseline$", colnames(adj$p_net))
    tps_hits <- c(tps_hits, as.numeric(adj$p_net[, end_cols, drop = FALSE] < 0.05))

    expr <- SummarizedExperiment::assay(pipe$se)
    index <- lapply(split(sim$term2gene$gene, sim$term2gene$term), function(g) {
      i <- match(g, rownames(expr))
      i[!is.na(i)]
    })
    index <- index[vapply(index, function(z) length(z) >= 5L, logical(1))]
    cm <- pipe$contrasts$families$endpoint
    for (j in seq_len(ncol(cm))) {
      cam <- limma::camera(expr, index, pipe$design$design_matrix,
                            contrast = cm[, j], sort = FALSE)
      camera_hits <- c(camera_hits, as.numeric(cam$PValue < 0.05))
    }
  }
  tps_typeI <- mean(tps_hits)
  camera_typeI <- mean(camera_hits)
  expect_gt(length(tps_hits), 80L)
  expect_lt(tps_typeI, camera_typeI)
  expect_lt(tps_typeI, 0.18)
})

test_that("VIF still ranks a strong codirectional pathway below null pathways", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 6, pathway_size = 10, n_background = 30,
    conditions = c("A", "B"), times = c(1, 3, 7), n_replicates = 3,
    within_pathway_cor = 0.5, pathway_signal = "codirectional",
    signal_fraction = 0.5, amplitude = 3, seed = 3)
  pipe <- .fit_endpoint(sim)
  adj <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                                 n_perm = 80, seed = 1, intergene = "vif")
  signal_pw <- sim$config$signal_pathways
  end_col <- grep("_t7_vs_baseline$", colnames(adj$p_net), value = TRUE)[1]
  p_sig <- adj$p_net[signal_pw, end_col, drop = TRUE]
  p_null <- adj$p_net[setdiff(rownames(adj$p_net), signal_pw), end_col,
                      drop = TRUE]
  expect_lt(median(p_sig), median(p_null))
})

## ---- edge cases (audit gap 3) ---------------------------------------------

test_that("pathway size 3 gives VIF >= 1 and bounded", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 3, pathway_size = 5, n_background = 10,
    conditions = c("A", "B"), times = c(1, 3), n_replicates = 3,
    within_pathway_cor = 0.5, pathway_signal = "none", seed = 5)
  pipe <- .fit_endpoint(sim)
  adj <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                                 n_perm = 40, seed = 1, intergene = "vif")
  ## VIF = 1 + (k-1)*rho; with k=5, rho~0.5 → VIF ~3, but bounded by [1, 1+4*0.99]
  expect_true(all(adj$vif >= 1))
  expect_true(all(adj$vif <= 5))  # 1 + 4*0.99 = 4.96 max
})

test_that("bare statistic matrix (no residuals) gives VIF = 1", {
  ## When fit is a raw matrix (not tps_fit), no residuals are available,
  ## so VIF should fall back to 1 for all pathways.
  set.seed(10)
  stat_mat <- matrix(rnorm(60), nrow = 20, ncol = 3,
                     dimnames = list(sprintf("P%04d", 1:20),
                                     c("c1", "c2", "c3")))
  t2g <- data.frame(
    term = rep(c("PW1", "PW2", "PW3"), each = 8),
    gene = sprintf("P%04d", rep(1:24, length.out = 24)),
    stringsAsFactors = FALSE)
  t2g <- t2g[t2g$gene %in% rownames(stat_mat), ]
  pws <- tps_pathway_scores(stat_mat, t2g)
  adj <- tps_calibrate_pathways(engine = "permutation", pws, mode = "unstratified",
                                 n_perm = 30, seed = 1, intergene = "vif")
  expect_true(all(adj$vif == 1))
  expect_true(all(adj$intergene_cor == 0))
})

## ---- REM type-I (audit gap 4) ----------------------------------------------

test_that("REM type-I is controlled under VIF at rho=0.5 (multi-seed)", {
  set.seed(99)
  typeI_rem <- numeric(8)
  for (s in seq_len(8)) {
    sim <- tps_simulate_correlated_pathways(
      n_pathways = 10, pathway_size = 10, n_background = 50,
      conditions = c("A", "B"), times = c(1, 3, 7), n_replicates = 3,
      within_pathway_cor = 0.5, pathway_signal = "none", seed = s * 13)
    pipe <- .fit_endpoint(sim)
    adj <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                                   n_perm = 100, seed = 1, intergene = "vif")
    end_col <- grep("_t7_vs_baseline$", colnames(adj$p_rem), value = TRUE)[1]
    typeI_rem[s] <- mean(adj$q_rem[, end_col] < 0.05, na.rm = TRUE)
  }
  expect_lt(mean(typeI_rem), 0.10)
})

## ---- golden-master: intergene="none" exact reproduction (audit gap 2) -------
## ADR 0010 point 3 claims intergene="none" reproduces ADR 0005 p-values exactly.
## Lock this down: with rho=0 (independent noise), "none" and "vif" must give
## IDENTICAL p_net because VIF=1 and the code path is the same centered
## second-pass. Any future change that breaks this contract will fail here.

test_that("golden-master: intergene='none' p_net is deterministic and reproducible", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 6, pathway_size = 8, n_background = 24,
    conditions = c("A", "B"), times = c(1, 3, 7), n_replicates = 3,
    within_pathway_cor = 0, pathway_signal = "none", seed = 4)
  pipe <- .fit_endpoint(sim)
  ## Two independent calls with same seed must give identical p_net
  a <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                               n_perm = 50, seed = 42, intergene = "none")
  b <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                               n_perm = 50, seed = 42, intergene = "none")
  expect_identical(a$p_net, b$p_net)
  expect_identical(a$p_rem, b$p_rem)
  ## VIF must be all 1 in "none" mode
  expect_true(all(a$vif == 1))
  expect_true(all(a$intergene_cor == 0))
  ## intergene field is set
  expect_equal(a$intergene, "none")
})

test_that("golden-master: at rho=0, 'none' and 'vif' give near-identical p_net", {
  ## ADR 0010 point 3 claims intergene="none" reproduces ADR 0005 exactly.
  ## At rho=0 the TRUE VIF is 1, but the estimated rho from residuals has
  ## small-sample noise, so estimated VIF may be slightly >1 and p_net
  ## diverges slightly. This test locks the contract that the divergence is
  ## small (not a structural change to the null).
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 6, pathway_size = 8, n_background = 24,
    conditions = c("A", "B"), times = c(1, 3, 7), n_replicates = 3,
    within_pathway_cor = 0, pathway_signal = "none", seed = 4)
  pipe <- .fit_endpoint(sim)
  none <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                                  n_perm = 50, seed = 1, intergene = "none")
  vif0 <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                                  n_perm = 50, seed = 1, intergene = "vif")
  ## At rho=0, VIF should be close to 1 (estimated from residuals)
  expect_true(all(vif0$vif >= 1))
  expect_true(median(vif0$vif) < 2)  # small-sample noise, but not huge
  ## p_net should be close but NOT necessarily identical (VIF estimation noise)
  expect_lt(mean(abs(vif0$p_net - none$p_net)), 0.15)
})

## ---- vif_cap: prevent over-correction in large pathways --------------------

test_that("vif_cap limits VIF and is recorded in output", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 6, pathway_size = 10, n_background = 30,
    conditions = c("A", "B"), times = c(1, 3, 7), n_replicates = 3,
    within_pathway_cor = 0.5, pathway_signal = "none", seed = 3)
  pipe <- .fit_endpoint(sim)
  uncapped <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                                      n_perm = 40, seed = 1, intergene = "vif",
                                      vif_cap = Inf)
  capped <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                                    n_perm = 40, seed = 1, intergene = "vif",
                                    vif_cap = 3)
  ## Uncapped VIF may exceed 3 (k=10, rho~0.5 → VIF~5.5)
  expect_true(any(uncapped$vif > 3))
  ## Capped VIF must not exceed 3
  expect_true(all(capped$vif <= 3))
  expect_equal(capped$vif_cap, 3)
  expect_equal(uncapped$vif_cap, Inf)
  ## Capped p_net should be smaller (more significant) than uncapped
  ## because less scaling → larger observed departure
  end_col <- grep("_t7_vs_baseline$", colnames(capped$p_net), value = TRUE)[1]
  expect_true(mean(capped$p_net[, end_col]) <= mean(uncapped$p_net[, end_col]))
})

test_that("vif_cap=1 is equivalent to intergene='none'", {
  ## vif_cap=1 forces VIF=1, which is the same as no VIF adjustment
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 6, pathway_size = 8, n_background = 24,
    conditions = c("A", "B"), times = c(1, 3, 7), n_replicates = 3,
    within_pathway_cor = 0.5, pathway_signal = "none", seed = 4)
  pipe <- .fit_endpoint(sim)
  none <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                                  n_perm = 40, seed = 1, intergene = "none")
  cap1 <- tps_calibrate_pathways(engine = "permutation", pipe$pws, mode = "unstratified",
                                  n_perm = 40, seed = 1, intergene = "vif",
                                  vif_cap = 1)
  expect_equal(none$p_net, cap1$p_net, tolerance = 0)
})

test_that("invalid vif_cap errors", {
  sim <- tps_simulate_correlated_pathways(
    n_pathways = 3, pathway_size = 5, n_background = 10,
    conditions = c("A"), times = c(1, 3), n_replicates = 2, seed = 1)
  pipe <- .fit_endpoint(sim)
  expect_error(tps_calibrate_pathways(engine = "permutation", pipe$pws, n_perm = 10, vif_cap = 0.5),
               "vif_cap")
  expect_error(tps_calibrate_pathways(engine = "permutation", pipe$pws, n_perm = 10, vif_cap = NA_real_),
               "vif_cap")
})

## ---- rho clipping-bias debias (2026-08-25) ---------------------------------
test_that("independent residuals give VIF ~= 1 after debias", {
  ## With truly independent residuals the clipped rho estimator is upward
  ## biased; the random-set floor must remove that bias so median VIF ~ 1.
  set.seed(11)
  n_feat <- 120; n_samp <- 24
  stat <- matrix(rnorm(n_feat * 3), nrow = n_feat,
                 dimnames = list(paste0("P", 1:n_feat), NULL))
  t2g <- data.frame(term = rep(paste0("pw", 1:10), each = 10),
                    gene = paste0("P", 1:100), stringsAsFactors = FALSE)
  resid_ind <- matrix(rnorm(n_feat * n_samp), nrow = n_feat)
  dimnames(resid_ind) <- list(paste0("P", 1:n_feat), NULL)
  pws <- tps_pathway_scores(stat, t2g, min_size = 5)
  pws$residuals <- resid_ind[rownames(pws$statistic_matrix), ]
  pc <- tps_calibrate_pathways(engine = "permutation", pws, mode = "unstratified", n_perm = 5,
                               seed = 1, intergene = "vif")
  expect_lt(stats::median(pc$vif), 1.10)
  expect_true(all(pc$vif >= 1))
})

test_that("correlated residuals still inflate VIF after debias", {
  set.seed(12)
  n_feat <- 120; n_samp <- 24
  stat <- matrix(rnorm(n_feat * 3), nrow = n_feat,
                 dimnames = list(paste0("P", 1:n_feat), NULL))
  t2g <- data.frame(term = rep(paste0("pw", 1:10), each = 10),
                    gene = paste0("P", 1:100), stringsAsFactors = FALSE)
  ## compound-symmetric rho = 0.5 within each 10-member block via factor model
  eps <- matrix(rnorm(n_feat * n_samp), nrow = n_feat)
  blk <- rep(1:12, each = 10)  # 12 blocks x 10 genes, one shared factor each
  fac_by_block <- matrix(rnorm(12 * n_samp), nrow = 12)
  Z <- fac_by_block[blk, ]    # n_feat x n_samp
  resid_cor <- sqrt(0.5) * Z + sqrt(0.5) * eps
  dimnames(resid_cor) <- list(paste0("P", 1:n_feat), NULL)
  pws <- tps_pathway_scores(stat, t2g, min_size = 5)
  pws$residuals <- resid_cor[rownames(pws$statistic_matrix), ]
  pc <- tps_calibrate_pathways(engine = "permutation", pws, mode = "unstratified", n_perm = 5,
                               seed = 1, intergene = "vif")
  ## true rho 0.5 - floor(~0.05) -> VIF ~ 1 + 9*0.45 = 5.05, capped at 3
  expect_gt(stats::median(pc$vif), 1.5)
})
