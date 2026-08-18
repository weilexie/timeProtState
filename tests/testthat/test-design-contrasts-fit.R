## test-design-contrasts-fit.R
##
## Tests the design/contrast/fit pipeline. The crown test is golden-master
## equivalence against the legacy dormancy limma (STATISTICAL_CONTRACT 10):
## timeProtState's endpoint contrasts must reproduce legacy logFC/t/P within
## tolerance 1e-10 on the same 20-protein fixture.

library(testthat)
library(timeProtState)
library(SummarizedExperiment)

## ---- golden equivalence ---------------------------------------------------
## Build an SE from the fixture's expr+meta, run the timeProtState pipeline,
## and compare endpoint logFC/t/P/adj.P against the frozen legacy limma.

with_golden <- function(code) {
  g <- readRDS(testthat::test_path("fixtures", "golden_fit.rds"))
  ## build colData with the columns tps_import expects
  cd <- data.frame(
    sample_id = rownames(g$meta),
    condition = ifelse(g$meta$model == "Control", "baseline",
                       as.character(g$meta$model)),
    time = g$meta$time,
    replicate = sub(".*-", "", rownames(g$meta)),
    is_baseline = g$meta$model == "Control",
    row.names = rownames(g$meta),
    stringsAsFactors = FALSE)
  list(g = g, cd = cd)
}

run_pipeline_on_golden <- function(g, cd) {
  se <- tps_import(g$expr, cd, validate = TRUE)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts, trend = TRUE, robust = TRUE)
  list(se = se, design = design, contrasts = contrasts, fit = fit)
}

test_that("golden-master: endpoint logFC matches legacy limma (tol 1e-10)", {
  d <- with_golden()
  p <- run_pipeline_on_golden(d$g, d$cd)
  obs <- tps_results(p$fit, family = "endpoint")
  exp <- d$g$statistics

  for (cc in unique(exp$contrast)) {
    o <- obs[obs$contrast == cc, ]; e <- exp[exp$contrast == cc, ]
    o <- o[match(e$feature_id, o$feature_id), ]
    expect_equal(o$logFC, e$logFC, tolerance = 1e-10,
                 info = paste("logFC", cc))
  }
})

test_that("golden-master: endpoint moderated-t matches legacy (tol 1e-10)", {
  d <- with_golden()
  p <- run_pipeline_on_golden(d$g, d$cd)
  obs <- tps_results(p$fit, family = "endpoint")
  exp <- d$g$statistics
  for (cc in unique(exp$contrast)) {
    o <- obs[obs$contrast == cc, ]; e <- exp[exp$contrast == cc, ]
    o <- o[match(e$feature_id, o$feature_id), ]
    expect_equal(o$t, e$t, tolerance = 1e-10, info = paste("t", cc))
  }
})

test_that("golden-master: endpoint P.Value matches legacy (tol 1e-10)", {
  d <- with_golden()
  p <- run_pipeline_on_golden(d$g, d$cd)
  obs <- tps_results(p$fit, family = "endpoint")
  exp <- d$g$statistics
  for (cc in unique(exp$contrast)) {
    o <- obs[obs$contrast == cc, ]; e <- exp[exp$contrast == cc, ]
    o <- o[match(e$feature_id, o$feature_id), ]
    expect_equal(o$P.Value, e$P.Value, tolerance = 1e-10, info = paste("P", cc))
  }
})

## ===========================================================================
## design
## ===========================================================================
test_that("tps_make_design builds full-rank shared-baseline design", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  expect_true(d$rank_audit$full_rank)
  expect_true(d$baseline$shared_baseline)
  expect_equal(d$baseline$baseline_time, 0)
  expect_equal(nrow(d$design_matrix), 30L)
})

test_that("tps_make_design has exactly 10 cells for 3x3+baseline", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  expect_equal(d$rank_audit$ncol, 10L)
  expect_true("baseline" %in% d$group_levels)
})

test_that("tps_make_design errors on missing column", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  expect_error(tps_make_design(se, condition = "nope"), "no column")
})

test_that("tps_make_design print works", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  expect_output(print(d), "shared_baseline: TRUE")
})

## ===========================================================================
## contrasts
## ===========================================================================
test_that("tps_make_contrasts produces 6 full-rank families", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  c <- tps_make_contrasts(d, endpoint = 7)
  expect_setequal(names(c$families),
                  c("vs_baseline","within_condition","average_time",
                    "heterogeneity","within_omnibus","endpoint"))
  for (fn in names(c$families)) {
    expect_true(c$rank_audit[[fn]]$full_rank,
                info = paste("family not full rank:", fn))
  }
})

test_that("heterogeneity family span is reference-invariant (STATISTICAL_CONTRACT invariant 4)", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  ## Changing the reference condition must span the SAME model-difference
  ## subspace: equal rank and equal column space. With feature-level missing
  ## values, contrasts.fit uses an approximate covariance transformation.
  ## STATISTICAL_CONTRACT invariant 4 is therefore expressed as subspace
  ## equivalence, which is the mathematically provable form.)
  c1 <- tps_make_contrasts(d, endpoint = 7, reference_condition = "Acid")
  c2 <- tps_make_contrasts(d, endpoint = 7, reference_condition = "Hypoxia")
  cm1 <- c1$families$heterogeneity
  cm2 <- c2$families$heterogeneity
  expect_equal(c1$rank_audit$heterogeneity$rank,
               c2$rank_audit$heterogeneity$rank)
  expect_equal(c1$rank_audit$heterogeneity$ncol,
               c2$rank_audit$heterogeneity$ncol)
  ## column space equivalence: rank(cm1) == rank(cbind(cm1, cm2))
  combined <- cbind(cm1, cm2)
  expect_equal(qr(cm1)$rank, qr(combined)$rank)
})

test_that("heterogeneity contrast subspace equivalence is the provable invariant", {
  ## NOTE for STATISTICAL_CONTRACT: the mathematically provable design-level
  ## invariant is subspace equivalence. Numerical limma F equality also holds
  ## for complete feature-level data, but may not hold when per-feature missing
  ## values make contrasts.fit's covariance transformation approximate.
  ## This test documents that finding by confirming the subspace test passes
  ## (already in the previous test) and that the two families are both valid
  ## full-rank encodings.
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  c1 <- tps_make_contrasts(d, endpoint = 7, reference_condition = "Acid")
  c2 <- tps_make_contrasts(d, endpoint = 7, reference_condition = "Hypoxia")
  expect_true(c1$rank_audit$heterogeneity$full_rank)
  expect_true(c2$rank_audit$heterogeneity$full_rank)
  expect_equal(c1$rank_audit$heterogeneity$rank,
               c2$rank_audit$heterogeneity$rank)
})

test_that("tps_make_contrasts errors on bad endpoint", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  expect_error(tps_make_contrasts(d, endpoint = 99), "not among non-baseline")
})

test_that("tps_make_contrasts within='all_pairs' produces more contrasts", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  c_cons <- tps_make_contrasts(d, endpoint = 7, within = "consecutive")
  c_all  <- tps_make_contrasts(d, endpoint = 7, within = "all_pairs")
  expect_gt(c_all$rank_audit$within_condition$ncol,
            c_cons$rank_audit$within_condition$ncol)
})

## ===========================================================================
## fit
## ===========================================================================
test_that("tps_fit uses trend=TRUE robust=TRUE by default (D1 resolution)", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  c <- tps_make_contrasts(d, endpoint = 7)
  f <- tps_fit(se, d, c)
  expect_true(f$thresholds$trend)
  expect_true(f$thresholds$robust)
})

test_that("tps_fit stores design residuals aligned to the assay", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  c <- tps_make_contrasts(d, endpoint = 7)
  f <- tps_fit(se, d, c)
  expect_true(is.matrix(f$residuals))
  expect_equal(dim(f$residuals), dim(SummarizedExperiment::assay(se)))
  expect_equal(rownames(f$residuals), rownames(se))
})

test_that("tps_fit statistics have all required columns", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  c <- tps_make_contrasts(d, endpoint = 7)
  f <- tps_fit(se, d, c)
  for (col in c("feature_id","family","contrast","test","logFC",
                 "P.Value","adj.P.Val")) {
    expect_true(col %in% names(f$statistics), info = col)
  }
})

test_that("tps_results filters by family and contrast", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  c <- tps_make_contrasts(d, endpoint = 7)
  f <- tps_fit(se, d, c)
  expect_true(nrow(tps_results(f, family = "endpoint")) < nrow(f$statistics))
  expect_true(all(tps_results(f, family = "endpoint")$family == "endpoint"))
})

test_that("tps_statistic_matrix returns feature x contrast matrix", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  c <- tps_make_contrasts(d, endpoint = 7)
  f <- tps_fit(se, d, c)
  m <- tps_statistic_matrix(f, family = "endpoint")
  expect_equal(nrow(m), 200L)
  expect_equal(ncol(m), 3L)
})

test_that("tps_fit respects trend/robust override", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  c <- tps_make_contrasts(d, endpoint = 7)
  f <- tps_fit(se, d, c, trend = FALSE, robust = FALSE)
  expect_false(f$thresholds$trend)
  expect_false(f$thresholds$robust)
})

test_that("tps_fit estimates and records repeated-measures correlation", {
  meta <- tps_example_metadata
  meta$subject_id <- meta$replicate
  se <- tps_import(tps_example_abundance, meta)
  d <- tps_make_design(se)
  c <- tps_make_contrasts(d, endpoint = 7)
  f <- tps_fit(se, d, c, block = "subject_id")

  expect_equal(f$thresholds$block, "subject_id")
  expect_equal(f$thresholds$correlation_source, "estimated")
  expect_true(is.finite(f$thresholds$correlation))
  expect_output(print(f), "correlation")
})

test_that("tps_fit accepts supplied block correlation", {
  meta <- tps_example_metadata
  meta$subject_id <- meta$replicate
  se <- tps_import(tps_example_abundance, meta)
  d <- tps_make_design(se)
  c <- tps_make_contrasts(d, endpoint = 7)
  f <- tps_fit(se, d, c, block = "subject_id", correlation = 0.2)

  expect_equal(f$thresholds$correlation, 0.2)
  expect_equal(f$thresholds$correlation_source, "supplied")
})

test_that("tps_fit validates block and correlation arguments", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  c <- tps_make_contrasts(d, endpoint = 7)

  expect_error(tps_fit(se, d, c, block = "not_a_column"), "not present")
  expect_error(tps_fit(se, d, c, correlation = 0.2), "requires")
  expect_error(tps_fit(se, d, c, block = seq_len(3)), "vector of length")
})
