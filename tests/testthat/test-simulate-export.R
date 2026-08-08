## test-simulate-export.R
##
## Tests Goal 15 (simulate) and Goal 14/export. The simulator must be
## reproducible, respect the shared-baseline structure, inject missingness,
## and produce a usable SE via tps_import for the full pipeline.

library(testthat)
library(timeProtState)
library(SummarizedExperiment)

## ===========================================================================
## tps_simulate_timecourse
## ===========================================================================
test_that("simulation produces correctly-sized abundance + metadata", {
  sim <- tps_simulate_timecourse(n_features = 50, conditions = c("A","B"),
                                  times = c(1,3,7), n_replicates = 3, seed = 1)
  expect_s3_class(sim, "tps_simulation")
  ## baseline (3) + 2 cond x 3 time x 3 rep = 3 + 18 = 21
  expect_equal(ncol(sim$abundance), 3L + 2L * 3L * 3L)
  expect_equal(nrow(sim$abundance), 50L)
  expect_true(sum(sim$col_data$is_baseline) == 3L)  # shared baseline once
})

test_that("simulation is reproducible with the same seed", {
  a <- tps_simulate_timecourse(n_features = 30, seed = 42)
  b <- tps_simulate_timecourse(n_features = 30, seed = 42)
  expect_identical(a$abundance, b$abundance)
  expect_identical(a$truth, b$truth)
})

test_that("simulation injects the requested missingness fraction", {
  sim <- tps_simulate_timecourse(n_features = 100, missingness = 0.1, seed = 1)
  frac <- mean(is.na(sim$abundance))
  expect_gt(frac, 0.05)
  expect_lt(frac, 0.15)
})

test_that("simulation truth covers all feature-condition combinations", {
  sim <- tps_simulate_timecourse(n_features = 20,
                                  conditions = c("A","B","C"), seed = 1)
  expect_equal(nrow(sim$truth), 20L * 3L)
  expect_setequal(unique(sim$truth$condition), c("A","B","C"))
})

test_that("simulation effect_mix must sum to 1", {
  expect_error(tps_simulate_timecourse(effect_mix = c(null = 0.5, up = 0.6)),
               "sum to 1")
})

test_that("simulation produces a null-dominant distribution by default", {
  sim <- tps_simulate_timecourse(n_features = 100, seed = 1)
  tab <- table(sim$truth$effect_type)
  expect_true(tab["null"] > tab["shared_up"])
})

test_that("simulated data runs through the full pipeline", {
  sim <- tps_simulate_timecourse(n_features = 80, conditions = c("Acid","Hypoxia"),
                                  times = c(1,3,7), n_replicates = 3, seed = 7)
  se <- tps_import(sim$abundance, sim$col_data)
  expect_true(tps_validate(se)$ok)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- tps_fit(se, design, contrasts)
  states <- tps_classify_states(fit)
  expect_gte(nrow(states$state_table), 80L)
})

test_that("simulation print works", {
  sim <- tps_simulate_timecourse(n_features = 30, seed = 1)
  expect_output(print(sim), "simulation")
})

test_that("simulation requires >= 2 non-baseline timepoints", {
  expect_error(tps_simulate_timecourse(times = 1), "at least 2")
})

## ===========================================================================
## tps_export
## ===========================================================================
test_that("tps_export writes result tables + manifest", {
  sim <- tps_simulate_timecourse(n_features = 40, conditions = c("A","B"),
                                  times = c(1,3), n_replicates = 2, seed = 1)
  se <- tps_import(sim$abundance, sim$col_data)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 3)
  fit <- tps_fit(se, design, contrasts)
  states <- tps_classify_states(fit)
  out <- file.path(tempdir(), "tps_export_test")
  unlink(out, recursive = TRUE)
  manifest <- tps_export(out, fit = fit, states = states)
  expect_true(dir.exists(out))
  expect_true(file.exists(file.path(out, "protein_statistics.csv")))
  expect_true(file.exists(file.path(out, "protein_states.csv")))
  expect_true(file.exists(file.path(out, "fit.rds")))
  expect_true(file.exists(file.path(out, "export_manifest.csv")))
  expect_gte(nrow(manifest), 3L)
})

test_that("tps_export handles partial inputs gracefully", {
  out <- file.path(tempdir(), "tps_export_partial")
  unlink(out, recursive = TRUE)
  manifest <- tps_export(out)
  expect_true(file.exists(file.path(out, "export_manifest.csv")))
  expect_equal(nrow(manifest), 0L)
})

test_that("tps_export errors on bad types", {
  expect_error(tps_export(tempdir(), fit = 123), "tps_fit")
  expect_error(tps_export(tempdir(), states = 123), "tps_states")
})
