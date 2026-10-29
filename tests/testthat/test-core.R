library(testthat)
library(timeProtState)

example_fit <- function() {
  se <- tps_import(tps_example_abundance, tps_example_metadata,
                   tps_example_annotation)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)
  list(se = se, design = design, contrasts = contrasts, fit = fit)
}

test_that("tps_import builds a SummarizedExperiment from the example", {
  se <- tps_import(tps_example_abundance, tps_example_metadata,
                   tps_example_annotation)
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(dim(se), c(200L, 30L))
  expect_equal(colnames(se), tps_example_metadata$sample_id)
})

test_that("tps_validate accepts the shared-baseline example", {
  se <- tps_import(tps_example_abundance, tps_example_metadata,
                   tps_example_annotation)
  report <- tps_validate(se)
  expect_true(report$ok)
  expect_true(report$design$shared_baseline)
})

test_that("tps_make_design is full rank with a shared baseline", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  design <- tps_make_design(se)
  expect_true(design$rank_audit$full_rank)
  expect_true(design$baseline$shared_baseline)
  expect_equal(nrow(design$design_matrix), 30L)
})

test_that("tps_make_contrasts and tps_fit return endpoint results", {
  pipe <- example_fit()
  expect_true(pipe$contrasts$rank_audit$endpoint$full_rank)
  res <- tps_results(pipe$fit, family = "endpoint")
  expect_true(is.data.frame(res))
  expect_gt(nrow(res), 0L)
  expect_true(all(c("feature_id", "contrast", "logFC", "t", "P.Value") %in% names(res)))
})

test_that("tps_pathway_scores returns NET, REM, and COH", {
  pipe <- example_fit()
  pws <- tps_pathway_scores(pipe$fit, tps_example_term2gene,
                            id_col = "locus_tag", family = "endpoint")
  expect_s3_class(pws, "tps_pathway_scores")
  expect_equal(dim(pws$scores$NET), dim(pws$scores$REM))
  expect_equal(dim(pws$scores$NET), dim(pws$scores$COH))
  expect_gt(nrow(pws$scores$NET), 0L)
})

test_that("CAMERA calibration and pathway time-path labels run", {
  pipe <- example_fit()
  pws <- tps_pathway_scores(pipe$fit, tps_example_term2gene,
                            id_col = "locus_tag", family = "vs_baseline")
  pws_c <- tps_calibrate_pathways(pws)
  expect_s3_class(pws_c, "tps_pathway_calibrated")
  expect_true(all(c("p_directional", "q_directional") %in% names(pws_c)))
  dyn <- tps_classify_pathway_dynamics(pws_c)
  expect_true(is.data.frame(dyn))
  expect_true(all(c("term", "condition", "time", "dynamic_state") %in% names(dyn)))
  dyn_cal <- tps_calibrate_dynamics(pws_c)
  expect_true(all(c("p_dynamic", "q_dynamic") %in% names(dyn_cal)))
})
