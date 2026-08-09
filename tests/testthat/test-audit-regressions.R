library(testthat)
library(timeProtState)

test_that("ordinary automatic data.frame row names are aligned by documented ids", {
  md <- tps_example_metadata
  rownames(md) <- NULL
  rd <- tps_example_annotation
  rownames(rd) <- NULL
  se <- tps_import(tps_example_abundance, md, rd)
  expect_identical(colnames(se), md$sample_id)
  expect_identical(rownames(se), rownames(tps_example_abundance))
})

test_that("baseline and time metadata are validated strictly", {
  md <- tps_example_metadata
  md$is_baseline <- as.numeric(md$is_baseline)
  md$is_baseline[1] <- 2
  expect_error(tps_import(tps_example_abundance, md), "only 0 or 1")

  md <- tps_example_metadata
  md$time <- paste0("day", md$time)
  expect_error(tps_import(tps_example_abundance, md), "finite numeric")

  bad <- matrix("not-a-number", nrow = 2, ncol = ncol(tps_example_abundance),
                dimnames = list(c("p1", "p2"), colnames(tps_example_abundance)))
  expect_error(tps_import(bad, tps_example_metadata), "non-numeric values")
})

test_that("condition-specific t0 controls are not pooled as a shared baseline", {
  md <- tps_example_metadata
  md$condition[md$is_baseline] <- c("Acid", "Hypoxia", "Starvation")
  se <- tps_import(tps_example_abundance, md, validate = FALSE)
  report <- tps_validate(se)
  expect_false(report$design$shared_baseline)
  expect_match(paste(report$warnings, collapse = " "), "NOT a valid shared-baseline")
  expect_error(tps_make_contrasts(tps_make_design(se)), "not shared-baseline")
})

test_that("incomplete condition by time grids are rejected", {
  keep <- !(tps_example_metadata$condition == "Acid" &
              tps_example_metadata$time == 1)
  md <- tps_example_metadata[keep, , drop = FALSE]
  mat <- tps_example_abundance[, md$sample_id, drop = FALSE]
  expect_error(tps_import(mat, md), "Incomplete condition-by-time grid")
})

test_that("batch columns are excluded from biological contrasts and fit safely", {
  md <- tps_example_metadata
  md$batch <- md$replicate
  se <- tps_import(tps_example_abundance, md)
  design <- tps_make_design(se, batch = "batch")
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  batch_rows <- grepl("^batch", rownames(contrasts$families$endpoint))
  expect_true(any(batch_rows))
  for (cm in contrasts$families) {
    expect_equal(nrow(cm), ncol(design$design_matrix))
    expect_true(all(cm[batch_rows, , drop = FALSE] == 0))
  }
  expect_s3_class(suppressWarnings(tps_fit(se, design, contrasts)), "tps_fit")
})

test_that("within omnibus is reported separately for every condition", {
  se <- tps_import(tps_example_abundance, tps_example_metadata)
  design <- tps_make_design(se)
  contrasts <- tps_make_contrasts(design, endpoint = 7)
  fit <- suppressWarnings(tps_fit(se, design, contrasts))
  out <- tps_results(fit, family = "within_omnibus")
  expect_equal(nrow(out), nrow(se) * length(contrasts$conditions))
  expect_setequal(unique(out$contrast),
                  paste0("within_omnibus_", contrasts$conditions))
})

test_that("endpoint directions are extracted per condition", {
  stats <- data.frame(
    feature_id = c("f1", "f1"), family = "endpoint",
    contrast = c("endpoint_A_t7_vs_baseline", "endpoint_B_t7_vs_baseline"),
    test = "t", logFC = c(2, -2), P.Value = c(0.001, 0.001),
    adj.P.Val = c(0.001, 0.001), stringsAsFactors = FALSE)
  fit <- structure(list(
    statistics = stats,
    contrasts = list(conditions = c("A", "B"), endpoint = 7)),
    class = "tps_fit")
  states <- tps_classify_states(fit)
  expect_equal(states$state_table$endpoint_dir_A, "up")
  expect_equal(states$state_table$endpoint_dir_B, "down")
  expect_equal(states$state_table$state, "cross_condition_opposite")
})

test_that("logFC gates do not discard omnibus evidence with undefined logFC", {
  stats <- data.frame(
    feature_id = "f1", family = "average_time",
    contrast = "average_time_omnibus", test = "F_omnibus",
    logFC = NA_real_, P.Value = 0.001, adj.P.Val = 0.001,
    stringsAsFactors = FALSE)
  fit <- structure(list(
    statistics = stats,
    contrasts = list(conditions = "A", endpoint = 7)),
    class = "tps_fit")
  rule <- tps_decision_rule(lfc_threshold = 0.74)
  states <- tps_classify_states(fit, rule)
  expect_equal(states$state_table$temporal_evidence, "average_time")
})

test_that("trajectory summaries honor median and repeat baseline only for display", {
  md <- data.frame(
    sample_id = c(paste0("b", 1:3),
                  paste0(rep(c("A1", "A2", "B1", "B2"), each = 3), "r", 1:3)),
    condition = c(rep("baseline", 3), rep(c("A", "A", "B", "B"), each = 3)),
    time = c(rep(0, 3), rep(c(1, 2, 1, 2), each = 3)),
    is_baseline = c(rep(TRUE, 3), rep(FALSE, 12)),
    stringsAsFactors = FALSE
  )
  mat <- matrix(0, nrow = 1, ncol = nrow(md),
                dimnames = list("P1", md$sample_id))
  mat[1, 4:6] <- c(1, 1, 100)
  se <- tps_import(mat, md)
  mean_tr <- tps_trajectory_matrix(se, summary = "mean")
  median_tr <- tps_trajectory_matrix(se, summary = "median")
  expect_equal(mean_tr$unique_cell["P1", "A_1"], 34)
  expect_equal(median_tr$unique_cell["P1", "A_1"], 1)
  expect_equal(sum(grepl("baseline", colnames(median_tr$display))), 2L)
  expect_equal(ncol(median_tr$unique_cell), 5L)
})

test_that("pathway membership is numeric and single-contrast calibration stays matrix-shaped", {
  stat <- matrix(seq(-3, 3, length.out = 8), ncol = 1,
                 dimnames = list(paste0("g", 1:8), "A_t1_vs_baseline"))
  t2g <- data.frame(term = rep(c("p1", "p2"), each = 4),
                    gene = paste0("g", 1:8))
  expect_no_warning(scores <- tps_pathway_scores(stat, t2g, min_size = 2))
  expect_s4_class(scores$membership, "dgCMatrix")
  cal <- tps_calibrate_pathways(scores, mode = "unstratified",
                                n_perm = 10, seed = 1)
  expect_equal(dim(cal$q_net), c(2L, 1L))
  expect_equal(dim(cal$q_rem), c(2L, 1L))
})

test_that("two-sided NET calibration is centered on the permutation null", {
  stat <- matrix(c(rep(10, 4), rep(8, 4)), ncol = 1,
                 dimnames = list(paste0("g", 1:8), "A_t1_vs_baseline"))
  t2g <- data.frame(term = rep(c("high", "low"), each = 4),
                    gene = paste0("g", 1:8))
  scores <- tps_pathway_scores(stat, t2g, min_size = 2)
  cal <- tps_calibrate_pathways(scores, mode = "unstratified",
                                n_perm = 100, seed = 9)
  expect_equal(unname(cal$p_net["high", ]), unname(cal$p_net["low", ]))
})

test_that("non-finite pathway statistics are dropped, never recoded as zero", {
  stat <- matrix(c(1, 2, NA, 4, 5, 6), ncol = 2,
                 dimnames = list(c("g1", "g2", "g3"), c("c1", "c2")))
  t2g <- data.frame(term = "p1", gene = c("g1", "g2", "g3"))
  scores <- tps_pathway_scores(stat, t2g, min_size = 2)
  expect_equal(scores$dropped_features_nonfinite, "g3")
  expect_false("g3" %in% rownames(scores$statistic_matrix))
})

test_that("protein id mapping is honored by state-pathway linking", {
  states <- structure(list(state_table = data.frame(
    feature_id = c("f1", "f2"), state = c("state1", "state2"))),
    class = "tps_states")
  dynamics <- data.frame(term = c("p1", "p2"),
                         dynamic_state = c("dyn1", "dyn2"))
  t2g <- data.frame(term = c("p1", "p2"), gene = c("g1", "g2"))
  linked <- tps_link_states_pathways(states, dynamics, t2g,
                                     protein_ids = c("g1", "g2"))
  hit <- linked$protein_state == "state1" & linked$pathway_dynamic == "dyn1"
  expect_equal(linked$observed_overlap[hit], 1L)
  expect_equal(linked$expected_overlap[hit], 0.5)
})

test_that("export preserves matrix identifiers and includes fit RDS in manifest", {
  fit <- structure(list(statistics = data.frame(feature_id = "f1")),
                   class = "tps_fit")
  trajectory <- structure(list(
    unique_cell = matrix(1:4, nrow = 2,
                         dimnames = list(c("f1", "f2"), c("A_1", "A_2"))),
    display = matrix(1:4, nrow = 2,
                     dimnames = list(c("f1", "f2"), c("A_1", "A_2"))),
    delta = matrix(1:4, nrow = 2,
                   dimnames = list(c("f1", "f2"), c("A_1", "A_2"))),
    missingness = data.frame(feature_id = c("f1", "f2"))),
    class = "tps_trajectory")
  stat <- matrix(1:8, nrow = 4,
                 dimnames = list(paste0("g", 1:4), c("c1", "c2")))
  scores <- tps_pathway_scores(
    stat, data.frame(term = "p1", gene = paste0("g", 1:4)), min_size = 2)
  out <- tempfile("tps-export-")
  manifest <- tps_export(out, fit = fit, trajectory = trajectory,
                         pathway_scores = scores)
  expect_true("fit.rds" %in% manifest$file)
  expect_identical(utils::read.csv(file.path(out, "trajectories_delta.csv"))$feature_id,
                   c("f1", "f2"))
  expect_identical(utils::read.csv(file.path(out, "pathway_NET.csv"))$term, "p1")
  expect_false(exists("written", envir = .GlobalEnv, inherits = FALSE))
})

test_that("simulated shared and opposite effects are coherent across conditions", {
  shared <- tps_simulate_timecourse(
    n_features = 5, conditions = c("A", "B"), times = c(1, 2),
    n_replicates = 1, effect_mix = c(shared_up = 1), noise_sd = 0, seed = 1)
  by_feature <- split(shared$truth, shared$truth$feature_id)
  expect_true(all(vapply(by_feature, function(x) length(unique(x$amplitude)) == 1L,
                         logical(1))))
  expect_true(all(vapply(by_feature, function(x) identical(x$sign, c(1, 1)),
                         logical(1))))

  opposite <- tps_simulate_timecourse(
    n_features = 5, conditions = c("A", "B"), times = c(1, 2),
    n_replicates = 1, effect_mix = c(opposite = 1), noise_sd = 0, seed = 1)
  by_feature <- split(opposite$truth, opposite$truth$feature_id)
  expect_true(all(vapply(by_feature, function(x) sum(x$sign) == 0, logical(1))))
})

test_that("delayed opposition injection reverses direction at the endpoint", {
  sim <- tps_simulate_timecourse(
    n_features = 10, conditions = "A", times = c(1, 3, 7),
    n_replicates = 1, effect_mix = c(null = 1), noise_sd = 0, seed = 1)
  t2g <- data.frame(term = "p1", gene = paste0("P", sprintf("%04d", 1:5)))
  injected <- tps_inject_pathway_effect(
    sim, t2g, term = "p1", effect_type = "delayed_opposition",
    condition = "A", amplitude = 3)
  delta <- injected$abundance["P0001", ] - injected$abundance["P0001", 1]
  expect_equal(as.numeric(delta[grepl("A_t", names(delta))]), c(-3, -3, 3))
})

test_that("REM trajectory significance uses q_rem and centroid order is stable", {
  skip_if_not_installed("ggplot2")
  nm <- c("A_t1_vs_baseline", "A_t3_vs_baseline", "A_t7_vs_baseline")
  score <- matrix(c(1, 2, 3), nrow = 1, dimnames = list("p1", nm))
  calibrated <- structure(list(
    scores = list(NET = score, REM = score, COH = score),
    q_net = matrix(1, 1, 3, dimnames = dimnames(score)),
    q_rem = matrix(c(0.01, 1, 1), 1, 3, dimnames = dimnames(score))),
    class = "tps_pathway_calibrated")
  plot <- tps_plot_pathway_trajectory(calibrated, "p1", metric = "REM")
  expect_identical(attr(plot, "data")$significant, c(TRUE, FALSE, FALSE))

  clusters <- structure(list(
    centroids = matrix(1:6, nrow = 3,
                       dimnames = list(c("t1", "t3", "t10"), c("1", "2")))),
    class = "tps_clusters")
  plot2 <- tps_plot_trajectories(clusters)
  expect_identical(levels(attr(plot2, "data")$cell), c("t1", "t3", "t10"))
})

test_that("display baseline labels use the actual baseline time", {
  sim <- tps_simulate_timecourse(
    n_features = 5, conditions = c("A", "B"), times = c(1, 2),
    n_replicates = 1, effect_mix = c(null = 1), seed = 1)
  sim$col_data$time[sim$col_data$is_baseline] <- -1
  se <- suppressWarnings(tps_import(sim$abundance, sim$col_data))
  trajectory <- tps_trajectory_matrix(se)
  expect_true(all(grepl("baseline_t-1$",
                         colnames(trajectory$display)[c(1, 4)])))
})
