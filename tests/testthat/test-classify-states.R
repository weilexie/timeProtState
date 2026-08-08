## test-classify-states.R
##
## Tests the evidence-decomposed state classification (STATISTICAL_CONTRACT 6).
## Verifies the three axes are computed, the default rule uses FDR (D2),
## and the legacy_msm profile reproduces raw-P behavior.

library(testthat)
library(timeProtState)

run_states <- function(rule = NULL) {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  d <- tps_make_design(se)
  c <- tps_make_contrasts(d, endpoint = 7)
  f <- tps_fit(se, d, c)
  tps_classify_states(f, rule = rule)
}

test_that("tps_decision_rule defaults to FDR < 0.05 with lfc 0 (D2)", {
  r <- tps_decision_rule()
  expect_equal(r$significance, "adj.P.Val")
  expect_equal(r$alpha, 0.05)
  expect_equal(r$lfc_threshold, 0)
})

test_that("tps_decision_rule legacy_msm profile uses raw P + 0.74", {
  r <- tps_decision_rule(significance = "P.Value", alpha = 0.05, lfc_threshold = 0.74)
  expect_equal(r$significance, "P.Value")
  expect_equal(r$lfc_threshold, 0.74)
})

test_that("tps_decision_rule rejects bad alpha", {
  expect_error(tps_decision_rule(alpha = 0), "alpha must be in")
  expect_error(tps_decision_rule(alpha = 1), "alpha must be in")
})

test_that("classify_states produces three evidence axes + final state", {
  s <- run_states()
  for (col in c("temporal_evidence","condition_heterogeneity",
                 "direction_pattern","state")) {
    expect_true(col %in% names(s$state_table), info = col)
  }
  expect_true(nrow(s$state_table) == 200L)
})

test_that("classify_states default uses adj.P.Val (more conservative than raw P)", {
  s_default <- run_states()
  s_rawp <- run_states(tps_decision_rule(significance = "P.Value"))
  n_default <- sum(s_default$state_table$state != "stable_or_unresolved")
  n_rawp <- sum(s_rawp$state_table$state != "stable_or_unresolved")
  expect_gt(n_rawp, 0L)
  expect_lte(n_default, n_rawp)
})

test_that("state distribution is non-degenerate (has at least 2 states)", {
  s <- run_states()
  expect_gte(length(unique(s$state_table$state)), 2L)
})

test_that("state rule manifest is recorded for audit", {
  s <- run_states()
  expect_true(is.data.frame(s$rule_manifest))
  expect_true("axis" %in% names(s$rule_manifest))
  expect_true("significance" %in% names(s$rule_manifest))
})

test_that("tps_state_table accessor works", {
  s <- run_states()
  st <- tps_state_table(s)
  expect_equal(nrow(st), 200L)
})

test_that("classify_states print works", {
  s <- run_states()
  expect_output(print(s), "state distribution")
})

test_that("endpoint direction columns are per-condition", {
  s <- run_states()
  expect_true("endpoint_dir_Acid" %in% names(s$state_table))
  expect_true("endpoint_dir_Hypoxia" %in% names(s$state_table))
  expect_true("endpoint_dir_Starvation" %in% names(s$state_table))
})

test_that("states provenance records classify_states stage", {
  s <- run_states()
  pv <- s$provenance[[1]]
  expect_equal(pv$stage, "classify_states")
})

test_that("temporal_reversal state is detected for sign-changing proteins", {
  ## Use simulated data with reversal effects + raw P rule (higher sensitivity)
  ## to verify the classifier can label sign-changing trajectories.
  sim <- tps_simulate_timecourse(n_features = 200,
                                  conditions = c("Acid","Hypoxia"),
                                  times = c(1, 3, 7), n_replicates = 4,
                                  effect_mix = c(null = 0.2, reversal = 0.8),
                                  noise_sd = 0.2, seed = 42)
  se <- tps_import(sim$abundance, sim$col_data)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)
  ## use raw P rule for higher sensitivity in detecting reversal
  states <- tps_classify_states(fit,
    rule = tps_decision_rule(significance = "P.Value", alpha = 0.05,
                              lfc_threshold = 0))
  st <- tps_state_table(states)
  ## the is_reversal column must exist and some proteins flagged
  expect_true("is_reversal" %in% names(st))
  expect_gt(sum(st$is_reversal), 0L)
})

test_that("rule manifest includes reversal axis", {
  s <- run_states()
  expect_true("reversal" %in% s$rule_manifest$axis)
})
