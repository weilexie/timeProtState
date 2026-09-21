## Default CAMERA/ROMER engine (ADR 0011).

library(testthat)
library(timeProtState)

.test_pipe <- function() {
  se <- tps_import(tps_example_abundance, tps_example_metadata,
                   tps_example_annotation)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "vs_baseline",
                             min_size = 5)
  list(se = se, design = design, contrasts = contrasts, fit = fit, pws = pws)
}

test_that("auto without block uses camera and fills directional/mixed p", {
  pipe <- .test_pipe()
  cal <- tps_calibrate_pathways(pipe$pws, nrot = 99, seed = 1)
  expect_equal(cal$engine_used, "camera")
  expect_true(is.matrix(cal$p_directional))
  expect_true(is.matrix(cal$p_mixed))
  expect_equal(dim(cal$p_directional), dim(pipe$pws$scores$NET))
  expect_true(all(cal$p_directional >= 0 & cal$p_directional <= 1, na.rm = TRUE))
  expect_true(all(cal$p_mixed >= 0 & cal$p_mixed <= 1, na.rm = TRUE))
  expect_equal(cal$q_directional, cal$q_net)
  expect_false(is.null(cal$engine_used))
})

test_that("explicit camera with block aborts", {
  se <- tps_import(tps_example_abundance, tps_example_metadata,
                   tps_example_annotation)
  SummarizedExperiment::colData(se)$subject_id <- rep(seq_len(10L), each = 3L)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts, block = "subject_id")
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "vs_baseline",
                             min_size = 5)
  expect_equal(
    tps_calibrate_pathways(pws, nrot = 99, seed = 1)$engine_used,
    "romer"
  )
  expect_error(
    tps_calibrate_pathways(pws, engine = "camera", nrot = 99),
    "camera"
  )
})

test_that("bare matrix scores abort auto/camera and allow permutation", {
  pipe <- .test_pipe()
  mat <- pipe$pws$statistic_matrix
  t2g <- tps_example_term2gene
  names(t2g) <- c("term", "gene")
  ## statistic matrix rownames are feature ids, not locus tags
  scores <- tps_pathway_scores(mat, data.frame(
    term = t2g$term, gene = rownames(mat)[match(t2g$gene, tps_example_annotation$locus_tag)],
    stringsAsFactors = FALSE
  ), min_size = 5)
  expect_error(tps_calibrate_pathways(scores, engine = "camera"), "permutation")
  perm <- tps_calibrate_pathways(scores, engine = "permutation",
                                  mode = "unstratified", n_perm = 20, seed = 1)
  expect_equal(perm$engine_used, "permutation")
  expect_true(all(perm$p_directional >= 0 & perm$p_directional <= 1))
})

test_that("shuffled term2gene order still matches CAMERA p-values by pathway name", {
  pipe <- .test_pipe()
  cal1 <- tps_calibrate_pathways(pipe$pws, nrot = 99, seed = 1)
  t2g <- tps_example_term2gene
  set.seed(1)
  t2g2 <- t2g[sample.int(nrow(t2g)), , drop = FALSE]
  pws2 <- tps_pathway_scores(pipe$fit, t2g2, id_col = "locus_tag",
                             family = "vs_baseline", min_size = 5)
  cal2 <- tps_calibrate_pathways(pws2, nrot = 99, seed = 1)
  common <- intersect(rownames(cal1$p_directional), rownames(cal2$p_directional))
  expect_gt(length(common), 1L)
  expect_equal(
    cal1$p_directional[common, , drop = FALSE],
    cal2$p_directional[common, , drop = FALSE],
    tolerance = 1e-10
  )
  expect_equal(
    cal1$p_mixed[common, , drop = FALSE],
    cal2$p_mixed[common, , drop = FALSE],
    tolerance = 1e-10
  )
})

test_that("example abundance has no missing values", {
  expect_equal(sum(is.na(tps_example_abundance)), 0L)
  expect_equal(dim(tps_example_abundance), c(200L, 30L))
})

test_that("protein lmFit still allows NA while CAMERA universe is complete-case", {
  ab <- tps_example_abundance
  ab[1:5, 1] <- NA_real_
  se <- tps_import(ab, tps_example_metadata, tps_example_annotation)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)
  expect_true(anyNA(fit$expression))
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "vs_baseline",
                             min_size = 5)
  expect_gt(pws$n_dropped_expr_na, 0L)
  expect_equal(sum(is.na(pws$expression)), 0L)
  cal <- tps_calibrate_pathways(pws, nrot = 99, seed = 1)
  expect_equal(cal$engine_used, "camera")
})
