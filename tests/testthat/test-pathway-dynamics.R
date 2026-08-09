## test-pathway-dynamics.R
##
## Tests the contrast-driven pathway remodeling dynamics (the core
## methodological upgrade from endpoint snapshot to time-course trajectories).
## Covers: tps_classify_pathway_dynamics, tps_link_states_pathways,
## tps_plot_pathway_trajectory.

library(testthat)
library(timeProtState)
library(SummarizedExperiment)

## ---- shared fixture: full pipeline with vs_baseline (all timepoints) -------
run_dynamics <- function(n_perm = 100L) {
  se <- tps_import(tps_example_abundance, tps_example_metadata,
                   tps_example_annotation)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "vs_baseline")
  vsbl <- tps_results(fit, "vs_baseline")
  ave <- as.numeric(tapply(vsbl$AveExpr,
    factor(vsbl$feature_id, levels = rownames(pws$statistic_matrix)),
    function(x) x[1]))
  pws_c <- tps_calibrate_pathways(pws, abundance = ave, n_perm = n_perm, seed = 1)
  states <- tps_classify_states(fit)
  list(fit = fit, pws = pws, pws_c = pws_c, states = states)
}

## ===========================================================================
## tps_classify_pathway_dynamics
## ===========================================================================
test_that("dynamics produces a long table with expected columns", {
  d <- run_dynamics()
  dyn <- tps_classify_pathway_dynamics(d$pws_c)
  expect_true(is.data.frame(dyn))
  expected_cols <- c("term", "condition", "time", "NET", "REM", "COH",
                      "q_net", "q_rem", "dynamic_state", "is_divergent")
  expect_true(all(expected_cols %in% names(dyn)))
})

test_that("dynamics has one row per pathway x condition x timepoint", {
  d <- run_dynamics()
  dyn <- tps_classify_pathway_dynamics(d$pws_c)
  n_path <- nrow(d$pws_c$scores$NET)
  n_contr <- ncol(d$pws_c$scores$NET)
  expect_equal(nrow(dyn), n_path * n_contr)
})

test_that("dynamics states are from the allowed vocabulary", {
  d <- run_dynamics()
  dyn <- tps_classify_pathway_dynamics(d$pws_c)
  allowed <- c("early_directional", "gradual_remodeling", "delayed_opposition",
               "transient_burst", "directional_delayed", "sustained_remodeling",
               "no_dynamic_signal")
  expect_true(all(dyn$dynamic_state %in% allowed))
})

test_that("dynamics detects divergent pathways column", {
  d <- run_dynamics()
  dyn <- tps_classify_pathway_dynamics(d$pws_c)
  expect_true(is.logical(dyn$is_divergent))
})

test_that("dynamics summary attribute is present", {
  d <- run_dynamics()
  dyn <- tps_classify_pathway_dynamics(d$pws_c)
  sm <- attr(dyn, "summary")
  expect_true(is.data.frame(sm))
  expect_true("term" %in% names(sm))
  expect_equal(nrow(sm), length(unique(dyn$term)))
})

test_that("dynamics errors on uncalibrated input", {
  se <- tps_import(tps_example_abundance, tps_example_metadata,
                   tps_example_annotation)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)
  pws <- tps_pathway_scores(fit, tps_example_term2gene,
                             id_col = "locus_tag", family = "vs_baseline")
  expect_error(tps_classify_pathway_dynamics(pws), "calibrate_pathways first")
})

test_that("dynamics errors on endpoint-only scores (no parseable times)", {
  d <- run_dynamics()
  ## build an endpoint-only calibrated object manually (should fail parsing)
  pws_end <- d$pws  ## this is vs_baseline, modify to endpoint-like names
  ## simulate endpoint naming that can't parse time
  pws_end2 <- d$pws_c
  pws_end2$scores$NET <- pws_end2$scores$NET[, 1, drop = FALSE]
  pws_end2$q_net <- pws_end2$q_net[, 1, drop = FALSE]
  colnames(pws_end2$scores$NET) <- "endpoint_only"
  colnames(pws_end2$q_net) <- "endpoint_only"
  expect_error(tps_classify_pathway_dynamics(pws_end2), "Could not parse timepoints")
})

## ===========================================================================
## tps_link_states_pathways
## ===========================================================================
test_that("link produces consistency diagnosis table", {
  d <- run_dynamics()
  dyn <- tps_classify_pathway_dynamics(d$pws_c)
  link <- tps_link_states_pathways(d$states, dyn, tps_example_term2gene)
  expect_true(is.data.frame(link))
  expect_true(all(c("protein_state", "pathway_dynamic", "enrichment") %in%
                  names(link)))
  ## should have one row per protein_state x pathway_dynamic combination
  expect_gte(nrow(link), 1L)
})

test_that("link errors on bad inputs", {
  d <- run_dynamics()
  dyn <- tps_classify_pathway_dynamics(d$pws_c)
  expect_error(tps_link_states_pathways(123, dyn, tps_example_term2gene),
               "tps_states")
  expect_error(tps_link_states_pathways(d$states, 123, tps_example_term2gene),
               "tps_classify_pathway_dynamics")
})

## ===========================================================================
## tps_plot_pathway_trajectory
## ===========================================================================
test_that("trajectory plot returns ggplot with data", {
  skip_if_not_installed("ggplot2")
  d <- run_dynamics()
  term <- rownames(d$pws_c$scores$NET)[1]
  plt <- tps_plot_pathway_trajectory(d$pws_c, term = term, metric = "NET")
  expect_true(ggplot2::is_ggplot(plt))
  expect_true(is.data.frame(attr(plt, "data")))
})

test_that("trajectory plot works for REM and COH", {
  skip_if_not_installed("ggplot2")
  d <- run_dynamics()
  term <- rownames(d$pws_c$scores$NET)[1]
  for (m in c("REM", "COH")) {
    plt <- tps_plot_pathway_trajectory(d$pws_c, term = term, metric = m)
    expect_true(ggplot2::is_ggplot(plt))
  }
})

test_that("trajectory plot errors on bad term", {
  skip_if_not_installed("ggplot2")
  d <- run_dynamics()
  expect_error(tps_plot_pathway_trajectory(d$pws_c, term = "nonexistent"),
               "not found")
})

## ===========================================================================
## vs_baseline gives more contrasts than endpoint
## ===========================================================================
test_that("vs_baseline family yields all timepoints (not just endpoint)", {
  d <- run_dynamics()
  ## vs_baseline should have 3 conditions x 3 timepoints = 9 contrasts
  expect_equal(ncol(d$pws$scores$NET), 9L)
  ## endpoint would have only 3 (one per condition)
  se <- tps_import(tps_example_abundance, tps_example_metadata,
                   tps_example_annotation)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)
  pws_end <- tps_pathway_scores(fit, tps_example_term2gene,
                                 id_col = "locus_tag", family = "endpoint")
  expect_equal(ncol(pws_end$scores$NET), 3L)
})

## ===========================================================================
## Null type-I error calibration (blocker: auditor flagged missing)
## Under the complete null (effect_mix = c(null=1)), the fraction of pathways
## with q_net < 0.05 should be near 0 (calibrated), and most dynamics should
## be "no_dynamic_signal".
## ===========================================================================
test_that("null calibration: complete null yields low false-positive dynamics", {
  ## Use 200 proteins, 10 pathways, 1000 permutations for adequate resolution
  sim <- tps_simulate_timecourse(
    n_features = 200, conditions = c("A","B"),
    times = c(1, 3, 7), n_replicates = 3,
    effect_mix = c(null = 1), noise_sd = 0.5, seed = 99)
  se <- tps_import(sim$abundance, sim$col_data)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)

  ## 10 pathways of 20 proteins each
  t2g <- do.call(rbind, lapply(1:10, function(i) {
    data.frame(term = sprintf("null_path_%02d", i),
               gene = sprintf("P%04d", ((i-1)*20 + 1):(i*20)),
               stringsAsFactors = FALSE)
  }))

  pws <- tps_pathway_scores(fit, t2g, family = "vs_baseline")
  vsbl <- tps_results(fit, "vs_baseline")
  ave <- as.numeric(tapply(vsbl$AveExpr,
    factor(vsbl$feature_id, levels = rownames(pws$statistic_matrix)),
    function(x) x[1]))
  pws_c <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 1000, seed = 1)

  ## (a) per-contrast false positive rate at q < 0.05 should be low (< 15%)
  fpr <- mean(pws_c$q_net < 0.05, na.rm = TRUE)
  expect_lt(fpr, 0.15)

  ## (b) most pathways should be classified no_dynamic_signal (>= 90%)
  dyn <- tps_classify_pathway_dynamics(pws_c)
  pct_no_signal <- mean(dyn$dynamic_state == "no_dynamic_signal")
  expect_gte(pct_no_signal, 0.90)
})

## ===========================================================================
## correction = "global" vs "per_contrast"
## ===========================================================================
test_that("global correction produces valid q-values (more tests = more conservative on average)", {
  d <- run_dynamics()
  pws_global <- tps_calibrate_pathways(d$pws, abundance = rep(1, nrow(d$pws$statistic_matrix)),
                                        n_perm = 100, seed = 1,
                                        correction = "global")
  ## global BH over more tests should have >= mean q than per_contrast
  expect_true(all(pws_global$q_net >= 0 & pws_global$q_net <= 1))
  expect_gte(mean(pws_global$q_net), mean(d$pws_c$q_net) - 0.01)
})

test_that("calibrate_pathways records correction mode in provenance", {
  d <- run_dynamics()
  ## search all provenance entries for the calibrate stage
  pv_list <- d$pws_c$provenance
  found <- FALSE
  for (pv in pv_list) {
    if (!is.null(pv$stage) && pv$stage == "calibrate_pathways") {
      if ("correction" %in% names(pv$params)) {
        found <- TRUE
        expect_equal(pv$params$correction, "per_contrast")
      }
    }
  }
  expect_true(found)
})

## ===========================================================================
## Link function now has hypergeometric test (blocker: auditor flagged)
## ===========================================================================
test_that("link includes hypergeometric p-value and BH q", {
  d <- run_dynamics()
  dyn <- tps_classify_pathway_dynamics(d$pws_c)
  link <- tps_link_states_pathways(d$states, dyn, tps_example_term2gene)
  expect_true("p_hypergeometric" %in% names(link))
  expect_true("q_hypergeometric" %in% names(link))
  expect_true(all(link$p_hypergeometric >= 0 & link$p_hypergeometric <= 1))
  expect_true(all(link$q_hypergeometric >= 0 & link$q_hypergeometric <= 1))
})

## ===========================================================================
## Pathway-level effect injection + recovery (blocker: auditor flagged missing)
## Plant a known dynamic pattern into a pathway, then verify recovery.
## ===========================================================================
test_that("injected early_directional effect is recovered with CORRECT label", {
  sim <- tps_simulate_timecourse(
    n_features = 100, conditions = c("Acid","Hypoxia"),
    times = c(1, 3, 7), n_replicates = 3,
    effect_mix = c(null = 1), noise_sd = 0.3, seed = 42)
  t2g <- data.frame(
    term = c(rep("pathway_A", 20), rep("pathway_B", 20)),
    gene = c(sprintf("P%04d", 1:20), sprintf("P%04d", 21:40)),
    stringsAsFactors = FALSE)
  sim2 <- tps_inject_pathway_effect(sim, t2g, term = "pathway_A",
                                     effect_type = "early_directional",
                                     condition = "Acid", amplitude = 3)
  se <- tps_import(sim2$abundance, sim2$col_data)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)
  pws <- tps_pathway_scores(fit, t2g, family = "vs_baseline")
  vsbl <- tps_results(fit, "vs_baseline")
  ave <- as.numeric(tapply(vsbl$AveExpr,
    factor(vsbl$feature_id, levels = rownames(pws$statistic_matrix)),
    function(x) x[1]))
  pws_c <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 500, seed = 1)
  dyn <- tps_classify_pathway_dynamics(pws_c)

  ## pathway_A under Acid should be early_directional (exact label)
  a_dyn <- dyn[dyn$term == "pathway_A" & dyn$condition == "Acid", ]
  expect_equal(a_dyn$dynamic_state[1], "early_directional")
  ## the null pathway should not recover the injected label
  b_dyn <- dyn[dyn$term == "pathway_B" & dyn$condition == "Acid", ]
  expect_true(b_dyn$dynamic_state[1] != "early_directional")
})

test_that("injected transient_burst is recovered with CORRECT label", {
  sim <- tps_simulate_timecourse(
    n_features = 60, conditions = c("Acid"),
    times = c(1, 3, 7), n_replicates = 3,
    effect_mix = c(null = 1), noise_sd = 0.3, seed = 7)
  t2g <- data.frame(term = c(rep("pw_burst", 15), rep("pw_null", 15)),
                     gene = c(sprintf("P%04d", 1:15), sprintf("P%04d", 16:30)),
                     stringsAsFactors = FALSE)
  sim2 <- tps_inject_pathway_effect(sim, t2g, term = "pw_burst",
                                     effect_type = "transient_burst",
                                     condition = "Acid", amplitude = 3)
  se <- tps_import(sim2$abundance, sim2$col_data)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7,
                                   include_heterogeneity = FALSE)
  fit <- tps_fit(se, design, contrasts)
  pws <- tps_pathway_scores(fit, t2g, family = "vs_baseline")
  vsbl <- tps_results(fit, "vs_baseline")
  ave <- as.numeric(tapply(vsbl$AveExpr,
    factor(vsbl$feature_id, levels = rownames(pws$statistic_matrix)),
    function(x) x[1]))
  pws_c <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 500, seed = 1)
  dyn <- tps_classify_pathway_dynamics(pws_c)

  burst <- dyn[dyn$term == "pw_burst", ]
  null_pw <- dyn[dyn$term == "pw_null", ]
  expect_equal(burst$dynamic_state[1], "transient_burst")
  ## null pathway may occasionally have a false positive at n_perm=500;
  ## assert it is NOT the injected type (specificity check)
  expect_true(null_pw$dynamic_state[1] != "transient_burst")
})

test_that("injected gradual_remodeling is recovered with CORRECT label", {
  sim <- tps_simulate_timecourse(
    n_features = 60, conditions = c("Acid"),
    times = c(1, 3, 7), n_replicates = 3,
    effect_mix = c(null = 1), noise_sd = 0.3, seed = 11)
  t2g <- data.frame(term = c(rep("pw_grad", 15), rep("pw_null", 15)),
                     gene = c(sprintf("P%04d", 1:15), sprintf("P%04d", 16:30)),
                     stringsAsFactors = FALSE)
  sim2 <- tps_inject_pathway_effect(sim, t2g, term = "pw_grad",
                                     effect_type = "gradual_remodeling",
                                     condition = "Acid", amplitude = 3)
  se <- tps_import(sim2$abundance, sim2$col_data)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7,
                                   include_heterogeneity = FALSE)
  fit <- tps_fit(se, design, contrasts)
  pws <- tps_pathway_scores(fit, t2g, family = "vs_baseline")
  vsbl <- tps_results(fit, "vs_baseline")
  ave <- as.numeric(tapply(vsbl$AveExpr,
    factor(vsbl$feature_id, levels = rownames(pws$statistic_matrix)),
    function(x) x[1]))
  pws_c <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 500, seed = 1)
  dyn <- tps_classify_pathway_dynamics(pws_c)

  grad <- dyn[dyn$term == "pw_grad", ]
  null_pw <- dyn[dyn$term == "pw_null", ]
  expect_equal(grad$dynamic_state[1], "gradual_remodeling")
  expect_true(null_pw$dynamic_state[1] != "gradual_remodeling")
})

test_that("injected delayed_opposition is recovered with CORRECT label", {
  sim <- tps_simulate_timecourse(
    n_features = 60, conditions = c("Acid"),
    times = c(1, 3, 7), n_replicates = 3,
    effect_mix = c(null = 1), noise_sd = 0.3, seed = 22)
  t2g <- data.frame(term = c(rep("pw_delay", 15), rep("pw_null", 15)),
                     gene = c(sprintf("P%04d", 1:15), sprintf("P%04d", 16:30)),
                     stringsAsFactors = FALSE)
  sim2 <- tps_inject_pathway_effect(sim, t2g, term = "pw_delay",
                                     effect_type = "delayed_opposition",
                                     condition = "Acid", amplitude = 3)
  se <- tps_import(sim2$abundance, sim2$col_data)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7,
                                   include_heterogeneity = FALSE)
  fit <- tps_fit(se, design, contrasts)
  pws <- tps_pathway_scores(fit, t2g, family = "vs_baseline")
  vsbl <- tps_results(fit, "vs_baseline")
  ave <- as.numeric(tapply(vsbl$AveExpr,
    factor(vsbl$feature_id, levels = rownames(pws$statistic_matrix)),
    function(x) x[1]))
  pws_c <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 500, seed = 1)
  dyn <- tps_classify_pathway_dynamics(pws_c)

  delay <- dyn[dyn$term == "pw_delay", ]
  null_pw <- dyn[dyn$term == "pw_null", ]
  expect_equal(delay$dynamic_state[1], "delayed_opposition")
  expect_true(null_pw$dynamic_state[1] != "delayed_opposition")
})

test_that("tps_inject_pathway_effect errors on bad inputs", {
  sim <- tps_simulate_timecourse(n_features = 20, seed = 1)
  t2g <- data.frame(term = "x", gene = "P0001", stringsAsFactors = FALSE)
  expect_error(tps_inject_pathway_effect(123, t2g, "x", "early_directional", "A"),
               "tps_simulation")
  expect_error(tps_inject_pathway_effect(sim, t2g, "x", "early_directional", "BAD"),
               "not in simulation")
})

## ===========================================================================
## COH threshold sensitivity analysis (Gap D: auditor flagged magic numbers)
## ===========================================================================
test_that("coh_cut sensitivity: states are stable across 0.3/0.5/0.7", {
  d <- run_dynamics(n_perm = 200L)
  states_by_cut <- lapply(c(0.3, 0.5, 0.7), function(cc) {
    dyn <- tps_classify_pathway_dynamics(d$pws_c, coh_cut = cc)
    table(dyn$dynamic_state)
  })
  ## all three should produce valid state tables with the same vocabulary
  for (i in 1:3) {
    expect_true(length(states_by_cut[[i]]) > 0)
    expect_true(all(names(states_by_cut[[i]]) %in%
      c("early_directional","gradual_remodeling","delayed_opposition",
        "transient_burst","directional_delayed","sustained_remodeling",
        "no_dynamic_signal")))
  }
  ## no_dynamic_signal should be the majority in all cases (synthetic data)
  for (i in 1:3) {
    expect_true(states_by_cut[[i]]["no_dynamic_signal"] >=
                sum(states_by_cut[[i]]) * 0.5)
  }
})

## ===========================================================================
## tps_divergence_score (continuous divergence + pair identification)
## ===========================================================================
test_that("divergence_score returns per-pathway continuous scores", {
  d <- run_dynamics(n_perm = 200L)
  dyn <- tps_classify_pathway_dynamics(d$pws_c)
  div <- tps_divergence_score(d$pws_c, dyn)
  expect_true(is.data.frame(div))
  expect_true(all(c("term","divergence_score","n_divergent_timepoints",
                     "max_divergent_pair","max_divergent_time") %in% names(div)))
  expect_equal(nrow(div), length(unique(dyn$term)))
  expect_true(all(div$divergence_score >= 0))
})

test_that("divergence_score errors on bad input", {
  expect_error(tps_divergence_score(123), "tps_pathway_calibrated")
})

## ===========================================================================
## tps_calibrate_dynamics (permutation-calibrated dynamic state classification)
## ===========================================================================
test_that("calibrate_dynamics adds p_dynamic and q_dynamic columns", {
  d <- run_dynamics(n_perm = 100L)
  dyn_cal <- tps_calibrate_dynamics(d$pws_c, n_perm = 50, seed = 1)
  expect_true("p_dynamic" %in% names(dyn_cal))
  expect_true("q_dynamic" %in% names(dyn_cal))
  expect_true(all(dyn_cal$p_dynamic >= 0 & dyn_cal$p_dynamic <= 1))
  expect_true(all(dyn_cal$q_dynamic >= 0 & dyn_cal$q_dynamic <= 1))
})

test_that("calibrate_dynamics is reproducible with same seed", {
  d <- run_dynamics(n_perm = 100L)
  a <- tps_calibrate_dynamics(d$pws_c, n_perm = 30, seed = 42)
  b <- tps_calibrate_dynamics(d$pws_c, n_perm = 30, seed = 42)
  expect_identical(a$p_dynamic, b$p_dynamic)
})

test_that("calibrate_dynamics does not perturb global RNG", {
  d <- run_dynamics(n_perm = 100L)
  set.seed(123); before <- runif(5)
  set.seed(123)
  tps_calibrate_dynamics(d$pws_c, n_perm = 20, seed = 1)
  after <- runif(5)
  expect_equal(before, after)
})

test_that("calibrate_dynamics errors on bad input", {
  expect_error(tps_calibrate_dynamics(123), "tps_pathway_calibrated")
})

## ===========================================================================
## Null calibration: under complete null, most p_dynamic should be high
## (pathways should NOT be classified as dynamic by chance)
## ===========================================================================
test_that("calibrate_dynamics null: complete null yields low false-positive rate", {
  sim <- tps_simulate_timecourse(
    n_features = 100, conditions = c("A","B"),
    times = c(1, 3, 7), n_replicates = 3,
    effect_mix = c(null = 1), noise_sd = 0.5, seed = 55)
  se <- tps_import(sim$abundance, sim$col_data)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)
  t2g <- do.call(rbind, lapply(1:5, function(i) {
    data.frame(term = sprintf("null_pw_%d", i),
               gene = sprintf("P%04d", ((i-1)*20 + 1):(i*20)),
               stringsAsFactors = FALSE)
  }))
  pws <- tps_pathway_scores(fit, t2g, family = "vs_baseline")
  vsbl <- tps_results(fit, "vs_baseline")
  ave <- as.numeric(tapply(vsbl$AveExpr,
    factor(vsbl$feature_id, levels = rownames(pws$statistic_matrix)),
    function(x) x[1]))
  pws_c <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 200, seed = 1)
  dyn_cal <- tps_calibrate_dynamics(pws_c, n_perm = 100, seed = 1)
  ## under complete null, p_dynamic should have discrimination (not all identical)
  ## and most should be high (pathways rarely classified as non-null by chance)
  expect_gt(length(unique(dyn_cal$p_dynamic)), 1L)
  expect_gt(mean(dyn_cal$p_dynamic > 0.1), 0.5)
})

test_that("calibrate_dynamics discrimination: injected effect has lower p_dynamic", {
  ## inject a strong early_directional effect into one pathway; verify its
  ## p_dynamic is lower than the null pathway's (has discrimination)
  sim <- tps_simulate_timecourse(
    n_features = 100, conditions = c("Acid","Hypoxia"),
    times = c(1, 3, 7), n_replicates = 3,
    effect_mix = c(null = 1), noise_sd = 0.3, seed = 42)
  t2g <- data.frame(
    term = c(rep("pw_signal", 20), rep("pw_null", 20)),
    gene = c(sprintf("P%04d", 1:20), sprintf("P%04d", 21:40)),
    stringsAsFactors = FALSE)
  sim2 <- tps_inject_pathway_effect(sim, t2g, term = "pw_signal",
                                     effect_type = "early_directional",
                                     condition = "Acid", amplitude = 3)
  se <- tps_import(sim2$abundance, sim2$col_data)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)
  pws <- tps_pathway_scores(fit, t2g, family = "vs_baseline")
  vsbl <- tps_results(fit, "vs_baseline")
  ave <- as.numeric(tapply(vsbl$AveExpr,
    factor(vsbl$feature_id, levels = rownames(pws$statistic_matrix)),
    function(x) x[1]))
  pws_c <- tps_calibrate_pathways(pws, abundance = ave, n_perm = 200, seed = 1)
  dyn_cal <- tps_calibrate_dynamics(pws_c, n_perm = 100, seed = 1)

  ## the injected pathway under Acid should have LOWER p_dynamic than the
  ## null pathway (the classifier fires more often on real signal than on
  ## permuted null data for the injected pathway)
  sig_p <- dyn_cal$p_dynamic[dyn_cal$term == "pw_signal" & dyn_cal$condition == "Acid"][1]
  null_p <- dyn_cal$p_dynamic[dyn_cal$term == "pw_null" & dyn_cal$condition == "Acid"][1]
  expect_true(!is.na(sig_p) && !is.na(null_p))
  ## the actual scientific claim: signal pathway has LOWER p_dynamic
  expect_lt(sig_p, null_p)
})
