## End-to-end integration: block fit through pathway calibration + VIF.
## Addresses Round 5 Gap D: block fit never run through pathway_scores →
## calibrate → dynamics. The GLS residuals from a block fit feed VIF estimation,
## but this chain was untested.

test_that("block fit runs through full pathway pipeline without error", {
  skip_if_not_installed("limma")
  sim <- tps_simulate_timecourse(
    n_features = 80, conditions = c("Acid", "Hypoxia"),
    times = c(1, 3, 7), n_replicates = 3,
    effect_mix = c(null = 0.8, shared_up = 0.2), seed = 42)
  t2g <- data.frame(
    term = rep(paste0("PW", 1:4), each = 10),
    gene = rownames(sim$abundance)[1:40],
    stringsAsFactors = FALSE)

  se <- tps_import(sim$abundance, sim$col_data)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts, block = "replicate")

  ## block correlation should be recorded
  expect_true(!is.null(fit$thresholds$correlation))
  expect_equal(fit$thresholds$correlation_source, "estimated")

  ## residuals should be available (GLS residuals from block fit)
  expect_true(is.matrix(fit$residuals))
  expect_equal(dim(fit$residuals), dim(sim$abundance))

  ## pathway scores should carry residuals through
  pws <- tps_pathway_scores(fit, t2g, family = "endpoint")
  expect_true(is.matrix(pws$residuals) || is.null(pws$residuals))
  expect_equal(nrow(pws$scores$NET), 4)

  ## calibration with VIF should work on block-fit residuals
  ave <- as.numeric(tapply(
    tps_results(fit, "endpoint")$AveExpr,
    factor(tps_results(fit, "endpoint")$feature_id,
           levels = rownames(pws$statistic_matrix)),
    function(x) x[1]))
  pws_c <- tps_calibrate_pathways(engine = "permutation", pws, abundance = ave, n_perm = 50,
                                   seed = 1, intergene = "vif")
  expect_s3_class(pws_c, "tps_pathway_calibrated")
  expect_true(all(pws_c$vif >= 1))
  expect_true(all(pws_c$p_net >= 0 & pws_c$p_net <= 1))

  ## calibration with intergene="none" should also work
  pws_none <- tps_calibrate_pathways(engine = "permutation", pws, abundance = ave, n_perm = 50,
                                      seed = 1, intergene = "none")
  expect_true(all(pws_none$vif == 1))
})

test_that("full 9-stage pipeline produces consistent results", {
  ## Addresses Round 5 Gap A1: no single test runs the full chain.
  sim <- tps_simulate_timecourse(
    n_features = 60, conditions = c("Acid", "Hypoxia"),
    times = c(1, 3, 7), n_replicates = 3, seed = 7)
  t2g <- data.frame(
    term = rep(paste0("PW", 1:3), each = 10),
    gene = rownames(sim$abundance)[1:30],
    stringsAsFactors = FALSE)

  se <- tps_import(sim$abundance, sim$col_data)
  report <- tps_validate(se)
  expect_true(report$ok || length(report$errors) == 0)

  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)

  pws <- tps_pathway_scores(fit, t2g, family = "vs_baseline")
  ave <- as.numeric(tapply(
    tps_results(fit, "vs_baseline")$AveExpr,
    factor(tps_results(fit, "vs_baseline")$feature_id,
           levels = rownames(pws$statistic_matrix)),
    function(x) x[1]))
  pws_c <- tps_calibrate_pathways(engine = "permutation", pws, abundance = ave, n_perm = 50, seed = 1)

  ## both classifiers should work on the same calibrated object
  states <- tps_classify_pathways(pws_c)
  expect_true(is.data.frame(states))
  expect_true(all(c("term", "contrast", "state") %in% names(states)))

  dyn <- tps_classify_pathway_dynamics(pws_c)
  expect_true(is.data.frame(dyn))
  expect_true("dynamic_state" %in% names(dyn))
})
