## test-pathway-scores.R
##
## Tests Goal 11 (tps_pathway_scores) and Goal 12 (tps_calibrate_pathways,
## tps_classify_pathways). Covers: membership construction, NET/REM/COH formula,
## size/coverage filtering, id_col mapping, permutation reproducibility,
## null mode switching, and golden-master equivalence on real dormancy data.

library(testthat)
library(timeProtState)
library(SummarizedExperiment)

## ---- shared pipeline fixture ----------------------------------------------
run_fit <- function() {
  se <- tps_import(tps_example_abundance, tps_example_metadata,
                   tps_example_annotation)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  tps_fit(se, design, contrasts)
}

## ===========================================================================
## Goal 11: tps_pathway_scores
## ===========================================================================
test_that("tps_pathway_scores builds scores from a fit + term2gene", {
  fit <- run_fit()
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  expect_s3_class(pws, "tps_pathway_scores")
  expect_equal(dim(pws$scores$NET), c(12L, 3L))
  expect_equal(dim(pws$scores$REM), c(12L, 3L))
  expect_equal(dim(pws$scores$COH), c(12L, 3L))
})

test_that("NET/REM/COH satisfy the definitional relationships", {
  fit <- run_fit()
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  NET <- pws$scores$NET; REM <- pws$scores$REM; COH <- pws$scores$COH
  ## REM >= |NET| always (since |.| sum >= signed sum)
  expect_true(all(REM >= abs(NET) - 1e-10))
  ## COH = |NET|/(REM+eps) in [0,1]
  expect_true(all(COH >= -1e-10) && all(COH <= 1 + 1e-10))
  ## COH matches the formula directly
  expect_equal(COH, abs(NET) / (REM + 1e-12), tolerance = 1e-10)
})

test_that("min_size filtering drops small pathways", {
  fit <- run_fit()
  pws_strict <- tps_pathway_scores(fit, tps_example_term2gene,
                                    id_col = "locus_tag", family = "endpoint",
                                    min_size = 15L)
  pws_loose <- tps_pathway_scores(fit, tps_example_term2gene,
                                   id_col = "locus_tag", family = "endpoint",
                                   min_size = 1L)
  expect_lt(nrow(pws_strict$scores$NET), nrow(pws_loose$scores$NET))
  expect_true(all(pws_strict$pathway_table$observed_size >= 15L))
})

test_that("id_col locus_tag mapping is used for membership", {
  fit <- run_fit()
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  ## universe reported by pws is the statistic_matrix rownames (feature ids),
  ## but membership is built on locus tags. At least one pathway should have
  ## observed members.
  expect_gt(sum(pws$pathway_table$observed_size), 0L)
})

test_that("tps_pathway_scores accepts a bare numeric matrix", {
  fit <- run_fit()
  m <- tps_statistic_matrix(fit, family = "endpoint")
  ## annotation rownames align with statistic_matrix rownames (= feature ids);
  ## remap rownames to locus so the example term2gene (gene=locus) matches.
  ann <- fit$annotation
  rownames(m) <- ann$locus_tag[match(rownames(m), rownames(ann))]
  pws <- tps_pathway_scores(m, tps_example_term2gene, type = "user_matrix")
  expect_s3_class(pws, "tps_pathway_scores")
  expect_equal(pws$statistic_type, "user_matrix")
  expect_gt(nrow(pws$scores$NET), 0L)
})

test_that("tps_pathway_scores errors on bad term2gene columns", {
  fit <- run_fit()
  bad <- tps_example_term2gene; names(bad)[1] <- "wrong"
  expect_error(tps_pathway_scores(fit, bad, id_col = "locus_tag"),
               "term.*term.*gene")
})

test_that("pathway_table records library vs observed size and coverage", {
  fit <- run_fit()
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  pt <- pws$pathway_table
  expect_true(all(c("term","library_size","observed_size","coverage") %in%
                  names(pt)))
  expect_true(all(pt$observed_size <= pt$library_size))
  expect_true(all(pt$coverage <= 1 + 1e-10))
})

test_that("print method works", {
  fit <- run_fit()
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  expect_output(print(pws), "pathway scores")
})

## ===========================================================================
## Goal 12: tps_calibrate_pathways
## ===========================================================================
abundance_for <- function(fit, pws) {
  end_stats <- tps_results(fit, family = "endpoint")
  ave <- tapply(end_stats$AveExpr,
                factor(end_stats$feature_id, levels = rownames(pws$statistic_matrix)),
                function(x) x[1L])
  as.numeric(ave)
}

test_that("tps_calibrate_pathways returns empirical p-values in [0,1]", {
  fit <- run_fit()
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  ave <- abundance_for(fit, pws)
  pws_c <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 100L,
                                   mode = "abundance_stratified", seed = 1L)
  expect_s3_class(pws_c, "tps_pathway_calibrated")
  expect_true(all(pws_c$p_net >= 0) && all(pws_c$p_net <= 1))
  expect_true(all(pws_c$p_rem >= 0) && all(pws_c$p_rem <= 1))
  expect_equal(dim(pws_c$q_net), dim(pws_c$p_net))
})

test_that("same seed reproduces p_net exactly (permutation determinism)", {
  fit <- run_fit()
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  ave <- abundance_for(fit, pws)
  c1 <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 100L, seed = 42L)
  c2 <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 100L, seed = 42L)
  expect_identical(c1$p_net, c2$p_net)
  expect_identical(c1$p_rem, c2$p_rem)
})

test_that("different seeds give different (but valid) p_net", {
  fit <- run_fit()
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  ave <- abundance_for(fit, pws)
  c1 <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 100L, seed = 1L)
  c2 <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 100L, seed = 2L)
  ## not necessarily all different, but the matrices should not be identical
  expect_false(identical(c1$p_net, c2$p_net))
})

test_that("unstratified mode works without abundance", {
  fit <- run_fit()
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  pws_c <- tps_calibrate_pathways(pws, n_perm = 100L,
                                   mode = "unstratified", seed = 1L)
  expect_true(all(pws_c$p_net >= 0) && all(pws_c$p_net <= 1))
})

test_that("abundance_stratified errors without abundance vector", {
  fit <- run_fit()
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  expect_error(tps_calibrate_pathways(pws, n_perm = 50L,
                                       mode = "abundance_stratified"),
               "abundance_stratified.*abundance")
})

test_that("calibration does not perturb the global RNG state", {
  fit <- run_fit()
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  ave <- abundance_for(fit, pws)
  set.seed(123); before <- runif(5)
  set.seed(123)
  tps_calibrate_pathways(pws, abundance = ave, n_perm = 50L, seed = 1L)
  after <- runif(5)
  expect_equal(before, after)
})

## ===========================================================================
## Goal 12: tps_classify_pathways
## ===========================================================================
test_that("tps_classify_pathways returns a long state table", {
  fit <- run_fit()
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  ave <- abundance_for(fit, pws)
  pws_c <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 100L, seed = 1L)
  st <- tps_classify_pathways(pws_c)
  expect_true(is.data.frame(st))
  expect_equal(nrow(st), nrow(pws$scores$NET) * ncol(pws$scores$NET))
  expect_true(all(c("term","contrast","NET","REM","COH","q_net","q_rem","state") %in%
                  names(st)))
})

test_that("state labels are from the allowed vocabulary", {
  fit <- run_fit()
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  ave <- abundance_for(fit, pws)
  pws_c <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 100L, seed = 1L)
  st <- tps_classify_pathways(pws_c)
  allowed <- c("directional_up","directional_down","mixed_remodeling",
               "directional_with_internal_opposition","no_strong_state")
  expect_true(all(st$state %in% allowed))
})

test_that("tps_classify_pathways errors on uncalibrated scores", {
  fit <- run_fit()
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  expect_error(tps_classify_pathways(pws), "calibrate_pathways first")
})

## ===========================================================================
## Golden-master pathway equivalence (real dormancy data)
## ===========================================================================
test_that("golden_pathway: NET/REM/COH match oracle (tol 1e-10)", {
  g <- readRDS(testthat::test_path("fixtures", "golden_pathway.rds"))
  ## statistic_matrix rownames are accessions; term2gene genes are locus tags.
  ## Remap rownames to locus so membership builds correctly.
  stat_m <- g$statistic_matrix
  rownames(stat_m) <- g$annotation$locus_tag[match(rownames(stat_m),
                                                     g$annotation$feature_id)]
  pws <- tps_pathway_scores(stat_m, g$term2gene,
                             type = "user_matrix", min_size = g$min_size)
  ## align pathway order
  o <- match(g$pathway_terms, rownames(pws$scores$NET))
  expect_equal(pws$scores$NET[o, ], g$NET, tolerance = 1e-10)
  expect_equal(pws$scores$REM[o, ], g$REM, tolerance = 1e-10)
  expect_equal(pws$scores$COH[o, ], g$COH, tolerance = 1e-10)
})

test_that("golden_pathway: unstratified p_net matches oracle exactly (seed=1)", {
  g <- readRDS(testthat::test_path("fixtures", "golden_pathway.rds"))
  stat_m <- g$statistic_matrix
  rownames(stat_m) <- g$annotation$locus_tag[match(rownames(stat_m),
                                                     g$annotation$feature_id)]
  pws <- tps_pathway_scores(stat_m, g$term2gene,
                             type = "user_matrix", min_size = g$min_size)
  pws_c <- tps_calibrate_pathways(pws, n_perm = g$n_perm,
                                   mode = "unstratified", seed = g$seed,
                                   intergene = "none")
  o <- match(g$pathway_terms, rownames(pws_c$p_net))
  expect_equal(pws_c$p_net[o, ], g$p_net_unstratified, tolerance = 0)
})

test_that("golden_pathway: p_rem is deterministic and in [0,1] (centered REM oracle)", {
  ## Round 2 changed REM to centered second-pass comparison. This test locks
  ## down the p_rem determinism (same seed → same p_rem) and range, since the
  ## golden fixture does not store p_rem (it predates the Round 2 fix).
  g <- readRDS(testthat::test_path("fixtures", "golden_pathway.rds"))
  stat_m <- g$statistic_matrix
  rownames(stat_m) <- g$annotation$locus_tag[match(rownames(stat_m),
                                                     g$annotation$feature_id)]
  pws <- tps_pathway_scores(stat_m, g$term2gene,
                             type = "user_matrix", min_size = g$min_size)
  pws_c <- tps_calibrate_pathways(pws, n_perm = g$n_perm,
                                   mode = "unstratified", seed = g$seed,
                                   intergene = "none")
  pr <- pws_c$p_rem
  expect_true(all(pr >= 0 & pr <= 1))
  ## determinism: re-run with same seed, must match exactly
  pws_c2 <- tps_calibrate_pathways(pws, n_perm = g$n_perm,
                                    mode = "unstratified", seed = g$seed,
                                    intergene = "none")
  expect_equal(pr, pws_c2$p_rem, tolerance = 0)
})

## ===========================================================================
## Edge cases (auditor flagged missing)
## ===========================================================================
test_that("tps_pathway_scores errors when no genes match universe", {
  fit <- run_fit()
  bad_t2g <- data.frame(term = "nonexistent", gene = "FAKE_GENE_999",
                         stringsAsFactors = FALSE)
  expect_error(tps_pathway_scores(fit, bad_t2g, id_col = "locus_tag",
                                   family = "endpoint"),
               "No pathways retained")
})

test_that("tps_calibrate_pathways works with n_perm = 1", {
  fit <- run_fit()
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "endpoint")
  ave <- abundance_for(fit, pws)
  pws_c <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 1L, seed = 1)
  expect_true(all(pws_c$p_net >= 0 & pws_c$p_net <= 1))
})
