## ============================================================================
## Simulate time-course proteomics for benchmarking
## ============================================================================
##
## Implements STATISTICAL_CONTRACT.md.
## Generates synthetic shared-baseline time-course proteomics with known truth,
## for type-I error, power, FDR, and state-recovery benchmarks.
##
## The simulator follows the data-generating model (spec 9.1):
##   y[i,c,t,r] = mu_i + alpha_ic + f_ic(t) + b_batch + epsilon
## with shared baseline generated ONCE (not tripled).
##
## Effect types (spec 9.2): null, monotonic, delayed, transient, reversal,
## condition-specific, cross-condition opposite, etc.

## ---- dynamic function library (spec 9.2) ----------------------------------
## Each returns a numeric vector of effect values at the given times.
.f_null <- function(t, amplitude = 0) rep(0, length(t))
.f_monotonic <- function(t, amplitude = 1) amplitude * (t - min(t)) /
  (max(t) - min(t))
.f_transient <- function(t, amplitude = 1, peak_frac = 0.4, width_frac = 0.15) {
  span <- max(t) - min(t)
  peak <- min(t) + peak_frac * span
  width <- width_frac * span
  amplitude * exp(-(t - peak)^2 / (2 * width^2))
}
.f_reversal <- function(t, amplitude = 1) {
  ## up then down: peaks mid-course, returns below baseline by end
  span <- max(t) - min(t)
  x <- (t - min(t)) / span
  amplitude * (4 * x * (1 - x) - 0.5) * 2  # parabola, range [-1, +1]
}

## ---- effect assignment plans ----------------------------------------------
## returns one coherent feature-level effect expanded across conditions.
.make_plan <- function(n_features, conditions, effect_mix) {
  ## Do not reset the RNG here: tps_simulate_timecourse already opens a local
  ## seeded stream. Resetting it would couple the plan and downstream noise in
  ## an undocumented way.
  rows <- list()
  ## every feature has an entry per condition (null by default)
  eff_types <- names(effect_mix)
  eff_probs <- as.numeric(effect_mix)
  for (i in seq_len(n_features)) {
    et <- sample(eff_types, size = 1L, prob = eff_probs)
    amp <- switch(et,
      null = 0,
      shared_up = runif(1, 0.8, 2),
      shared_down = runif(1, 0.8, 2),
      condition_specific_up = runif(1, 0.8, 2),
      condition_specific_down = runif(1, 0.8, 2),
      opposite = runif(1, 0.8, 2),
      transient = runif(1, 0.8, 1.5),
      reversal = runif(1, 0.8, 1.5))
    targeted <- et %in% c("condition_specific_up", "condition_specific_down",
                           "opposite")
    target <- if (targeted) sample(conditions, 1L) else NA_character_
    anchor_sign <- if (et == "opposite") sample(c(-1, 1), 1L) else 1
    for (cc in conditions) {
      sgn <- switch(et,
        null = 0,
        shared_up = 1,
        shared_down = -1,
        condition_specific_up = if (cc == target) 1 else 0,
        condition_specific_down = if (cc == target) -1 else 0,
        opposite = if (cc == target) anchor_sign else -anchor_sign,
        transient = 1,
        reversal = 1)
      rows[[length(rows) + 1L]] <- data.frame(
        feature_id = sprintf("P%04d", i),
        condition = cc, effect_type = et,
        amplitude = amp, sign = sgn, target_condition = target,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

## resolve an effect type + condition context to a per-time function value
.resolve_effect <- function(effect_type, t, amp, sgn) {
  base <- switch(effect_type,
    null = .f_null(t),
    shared_up = sgn * .f_monotonic(t, amp),
    shared_down = sgn * .f_monotonic(t, amp),
    condition_specific_up = sgn * .f_monotonic(t, amp),
    condition_specific_down = sgn * .f_monotonic(t, amp),
    opposite = sgn * .f_monotonic(t, amp),
    transient = sgn * .f_transient(t, amp),
    reversal = sgn * .f_reversal(t, amp),
    .f_null(t))
  base
}

## ---- tps_simulate_timecourse -----------------------------------------------
#' Simulate a shared-baseline time-course proteomics experiment
#'
#' Generates synthetic protein-level abundance data with a known ground-truth
#' effect structure, suitable for benchmarking type-I error, power, FDR, and
#' state recovery. The design is shared-baseline: a single baseline timepoint
#' shared across all conditions (spec 9.1).
#'
#' @param n_features integer; number of proteins.
#' @param conditions character vector of condition names (excluding baseline).
#' @param times numeric vector of non-baseline timepoints (baseline = min, added
#'   automatically as time 0 relative).
#' @param n_replicates integer; biological replicates per cell.
#' @param effect_mix named list of effect-type probabilities (must sum to 1).
#'   See details for available types.
#' @param noise_sd numeric; residual SD on log2 scale.
#' @param missingness numeric in [0,1]; fraction of MCAR missing values.
#' @param seed integer.
#' @return a list of class \code{tps_simulation} with:
#'   \code{abundance} (matrix), \code{col_data} (metadata),
#'   \code{truth} (per-feature-condition effect plan), \code{config}.
#'
#' @details Available effect types in \code{effect_mix}: \code{null},
#'   \code{shared_up}, \code{shared_down}, \code{condition_specific_up},
#'   \code{condition_specific_down}, \code{opposite}, \code{transient},
#'   \code{reversal}.
#'
#' @examples
#' sim <- tps_simulate_timecourse(n_features = 100, conditions = c("A","B"),
#'                                 times = c(1, 3, 7), n_replicates = 3,
#'                                 seed = 1)
#' dim(sim$abundance)
#' head(sim$truth)
#'
#' @export
tps_simulate_timecourse <- function(n_features = 200L,
                                     conditions = c("Acid","Hypoxia","Starvation"),
                                     times = c(1, 3, 7),
                                     n_replicates = 3L,
                                     effect_mix = c(
                                       null = 0.5,
                                       shared_up = 0.1,
                                       shared_down = 0.1,
                                       condition_specific_up = 0.08,
                                       condition_specific_down = 0.07,
                                       opposite = 0.05,
                                       transient = 0.05,
                                       reversal = 0.05),
                                     noise_sd = 0.5,
                                     missingness = 0,
                                     seed = 1L) {
  arg_capture <- list(call = "tps_simulate_timecourse",
                       n_features = n_features, conditions = conditions,
                       times = times, n_replicates = n_replicates,
                       noise_sd = noise_sd, missingness = missingness,
                       seed = seed,
                       timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  if (length(n_features) != 1L || is.na(n_features) || n_features < 1L ||
      n_features != as.integer(n_features)) {
    tps_abort("n_features must be one positive integer.")
  }
  if (!is.character(conditions) || length(conditions) < 1L ||
      anyNA(conditions) || any(!nzchar(conditions)) || anyDuplicated(conditions)) {
    tps_abort("conditions must be unique, non-empty character values.")
  }
  if (!is.numeric(times) || length(times) < 2L || any(!is.finite(times)) ||
      any(times <= 0) || anyDuplicated(times)) {
    tps_abort("times must contain at least 2 unique, positive finite values.")
  }
  times <- sort(times)
  if (length(n_replicates) != 1L || is.na(n_replicates) ||
      n_replicates < 1L || n_replicates != as.integer(n_replicates)) {
    tps_abort("n_replicates must be one positive integer.")
  }
  allowed_effects <- c("null", "shared_up", "shared_down",
                       "condition_specific_up", "condition_specific_down",
                       "opposite", "transient", "reversal")
  if (!is.numeric(effect_mix) || is.null(names(effect_mix)) ||
      any(!nzchar(names(effect_mix))) || anyDuplicated(names(effect_mix)) ||
      any(!is.finite(effect_mix)) || any(effect_mix < 0)) {
    tps_abort("effect_mix must be a named, non-negative finite numeric vector.")
  }
  if (abs(sum(effect_mix) - 1) > 1e-6) {
    tps_abort("effect_mix probabilities must sum to 1 (got %.4f).", sum(effect_mix))
  }
  if (!all(names(effect_mix) %in% allowed_effects)) {
    tps_abort("effect_mix contains unsupported effect type(s): %s.",
              paste(setdiff(names(effect_mix), allowed_effects), collapse = ", "))
  }
  if (length(conditions) < 2L && "opposite" %in% names(effect_mix) &&
      effect_mix[["opposite"]] > 0) {
    tps_abort("effect type 'opposite' requires at least two conditions.")
  }
  if (length(noise_sd) != 1L || is.na(noise_sd) || !is.finite(noise_sd) ||
      noise_sd < 0) {
    tps_abort("noise_sd must be one non-negative finite value.")
  }
  if (length(missingness) != 1L || is.na(missingness) ||
      !is.finite(missingness) || missingness < 0 || missingness > 1) {
    tps_abort("missingness must be in [0,1].")
  }
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
    tps_abort("seed must be one finite value.")
  }

  withr::with_seed(seed, {
    ## baseline time = 0 (relative); the shared baseline
    baseline_time <- 0
  all_times <- c(baseline_time, times)

  ## truth plan: per feature per condition
  plan <- .make_plan(n_features, conditions, effect_mix)

  ## baseline abundance per feature
  mu0 <- stats::rnorm(n_features, mean = 18, sd = 4)

  ## build sample metadata + abundance matrix
  sample_rows <- list()
  expr_cols <- list()
  sid <- 0L
  ## baseline cells (one shared baseline group)
  for (r in seq_len(n_replicates)) {
    sid <- sid + 1L
    sname <- sprintf("baseline_R%d", r)
    sample_rows[[sid]] <- data.frame(
      sample_id = sname, condition = "baseline", time = baseline_time,
      replicate = sprintf("R%d", r), is_baseline = TRUE,
      stringsAsFactors = FALSE)
  }
  ## condition x time cells
  for (cc in conditions) {
    for (tt in times) {
      for (r in seq_len(n_replicates)) {
        sid <- sid + 1L
        sname <- sprintf("%s_t%s_R%d", cc, tt, r)
        sample_rows[[sid]] <- data.frame(
          sample_id = sname, condition = cc, time = tt,
          replicate = sprintf("R%d", r), is_baseline = FALSE,
          stringsAsFactors = FALSE)
      }
    }
  }
  col_data <- do.call(rbind, sample_rows)
  rownames(col_data) <- col_data$sample_id

  ## abundance matrix: features x samples
  abundance <- matrix(NA_real_, nrow = n_features, ncol = nrow(col_data),
                       dimnames = list(sprintf("P%04d", seq_len(n_features)),
                                       col_data$sample_id))

  for (i in seq_len(n_features)) {
    feat_id <- sprintf("P%04d", i)
    for (j in seq_len(nrow(col_data))) {
      cd <- col_data[j, ]
      if (cd$is_baseline) {
        abundance[i, j] <- mu0[i] + stats::rnorm(1, sd = noise_sd)
      } else {
        cc <- cd$condition; tt <- cd$time
        pr <- plan[plan$feature_id == feat_id & plan$condition == cc, ]
        ## pass the FULL time grid so shape functions (esp. reversal/transient)
        ## normalize correctly; then extract the current timepoint's value.
        eff <- .resolve_effect(pr$effect_type, all_times,
                               pr$amplitude, pr$sign)
        eff_val <- eff[match(tt, all_times)]
        abundance[i, j] <- mu0[i] + eff_val + stats::rnorm(1, sd = noise_sd)
      }
    }
  }

  ## inject MCAR missingness
  if (missingness > 0) {
    n_miss <- floor(missingness * length(abundance))
    miss_idx <- sample(seq_along(abundance), n_miss)
    abundance[miss_idx] <- NA
  }

    structure(
    list(
      abundance = abundance,
      col_data = col_data,
      truth = plan,
      config = list(n_features = n_features, conditions = conditions,
                    times = times, n_replicates = n_replicates,
                    effect_mix = effect_mix, noise_sd = noise_sd,
                    missingness = missingness, seed = seed,
                    baseline_time = baseline_time),
      provenance = list(.make_provenance(arg_capture, list(), stage = "simulate"))
    ),
    class = "tps_simulation"
    )
  })
}

## ---- methods ---------------------------------------------------------------
#' @export
print.tps_simulation <- function(x, ...) {
  cat("<timeProtState simulation>\n")
  cat(sprintf("  features: %d   samples: %d   conditions: %s\n",
              x$config$n_features, ncol(x$abundance),
              paste(x$config$conditions, collapse = ", ")))
  cat(sprintf("  times: %s   replicates: %d   noise_sd: %s\n",
              paste(x$config$times, collapse = ","),
              x$config$n_replicates, x$config$noise_sd))
  if (x$config$missingness > 0) {
    cat(sprintf("  missingness: %.0f%%\n", 100 * x$config$missingness))
  }
  tab <- table(x$truth$effect_type)
  cat("  effect distribution (per feature-condition):\n")
  for (nm in names(tab)) cat(sprintf("    %-28s %d\n", nm, tab[nm]))
  invisible(x)
}

## ---- pathway-level effect injection (for validating dynamics) --------------
## Injects a coherent effect shape into all members of a specified pathway,
## so tps_classify_pathway_dynamics can be validated by recovery.
##
## Returns a MODIFIED simulation object (abundance + truth updated).

#' Inject a pathway-level effect into a simulation
#'
#' Modifies a \code{tps_simulation} object by injecting a coherent temporal
#' effect into all proteins that are members of a specified pathway. This
#' enables validation of \code{\link{tps_classify_pathway_dynamics}} by
#' planting a known dynamic pattern and checking recovery.
#'
#' @param sim a \code{tps_simulation} object from
#'   \code{\link{tps_simulate_timecourse}}.
#' @param term2gene pathway membership (columns: term, gene).
#' @param term the pathway term to inject into.
#' @param effect_type one of \code{"early_directional"},
#'   \code{"gradual_remodeling"}, \code{"transient_burst"},
#'   \code{"delayed_opposition"}.
#' @param condition the condition to inject into (must be in
#'   \code{sim$config$conditions}).
#' @param amplitude numeric; effect magnitude on log2 scale.
#' @return a modified \code{tps_simulation} object.
#'
#' @examples
#' sim <- tps_simulate_timecourse(n_features = 50, seed = 1)
#' t2g_sim <- data.frame(term = "pw1", gene = paste0("P", sprintf("%04d", 1:10)),
#'                       stringsAsFactors = FALSE)
#' sim2 <- tps_inject_pathway_effect(sim, t2g_sim, term = "pw1",
#'   effect_type = "early_directional", condition = "Acid")
#' @export
tps_inject_pathway_effect <- function(sim, term2gene, term,
                                       effect_type = c("early_directional",
                                                       "gradual_remodeling",
                                                       "transient_burst",
                                                       "delayed_opposition"),
                                       condition, amplitude = 2) {
  effect_type <- match.arg(effect_type)
  if (!inherits(sim, "tps_simulation") ||
      inherits(sim, "tps_simulation_correlated")) {
    tps_abort("sim must be a tps_simulation object from tps_simulate_timecourse (not tps_simulate_correlated_pathways, which has a different truth schema).")
  }
  if (!condition %in% sim$config$conditions) {
    tps_abort("condition '%s' not in simulation conditions (%s).",
              condition, paste(sim$config$conditions, collapse = ", "))
  }
  if (!is.data.frame(term2gene) ||
      !all(c("term", "gene") %in% names(term2gene))) {
    tps_abort("term2gene must have columns 'term' and 'gene'.")
  }
  if (length(term) != 1L || is.na(term) || !nzchar(term)) {
    tps_abort("term must be one non-empty value.")
  }
  if (length(amplitude) != 1L || is.na(amplitude) || !is.finite(amplitude) ||
      amplitude <= 0) {
    tps_abort("amplitude must be one positive finite value.")
  }
  members <- unique(as.character(term2gene$gene[term2gene$term == term]))
  ## map to feature ids present in the simulation
  feat_ids <- rownames(sim$abundance)
  members_in_sim <- intersect(members, feat_ids)
  if (length(members_in_sim) == 0L) {
    tps_abort("No members of term '%s' found in simulation features.", term)
  }

  times <- sim$config$times
  baseline_time <- sim$config$baseline_time
  all_times <- c(baseline_time, times)

  ## define the effect shape over all_times
  shape <- switch(effect_type,
    early_directional = rep(amplitude, length(all_times)),
    gradual_remodeling = {
      ## REM increases over time; t1 has NO signal (distinguishes from early_directional)
      v <- rep(0, length(all_times))
      for (k in 2:length(all_times)) {
        v[k] <- amplitude * (k - 1) / (length(all_times) - 1)
      }
      v[1] <- 0  # baseline: no effect
      v[2] <- 0  # t1: no effect (signal starts from t2)
      v
    },
    transient_burst = {
      v <- rep(0, length(all_times))
      ## peak at the true middle non-baseline timepoint
      ## all_times = c(baseline, t1, ..., tN); peak at ceiling of midpoint
      mid <- ceiling((length(all_times) + 1) / 2)
      mid <- min(mid, length(all_times) - 1L)  # never the last (endpoint)
      mid <- max(mid, 2L)  # never the first (baseline)
      v[mid] <- amplitude
      v
    },
    delayed_opposition = {
      v <- rep(0, length(all_times))
      if (length(all_times) > 2L) v[2:(length(all_times) - 1L)] <- -amplitude
      v[length(all_times)] <- amplitude
      v
    },
    rep(0, length(all_times))
  )
  ## shape[1] is baseline_time (effect = 0 at baseline)
  shape[1] <- 0

  ## inject into abundance for the target condition's samples
  abund <- sim$abundance
  cd <- sim$col_data
  for (feat in members_in_sim) {
    for (j in seq_len(ncol(abund))) {
      if (cd$condition[j] == condition && !cd$is_baseline[j]) {
        tt <- cd$time[j]
        idx <- match(tt, all_times)
        if (!is.na(idx)) {
          abund[feat, j] <- abund[feat, j] + shape[idx]
        }
      }
    }
  }
  sim$abundance <- abund

  ## record injection in truth
  inj <- data.frame(
    feature_id = members_in_sim, condition = condition,
    effect_type = paste0("injected_", effect_type),
    amplitude = amplitude, sign = 1, target_condition = condition,
    stringsAsFactors = FALSE)
  sim$truth <- rbind(sim$truth, inj)
  sim$provenance[[length(sim$provenance) + 1L]] <- .make_provenance(
    list(call = "tps_inject_pathway_effect", term = term,
         effect_type = effect_type, condition = condition,
         amplitude = amplitude,
         timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    sim, stage = "inject_pathway_effect")

  sim
}

## ---- correlated-pathway simulator ------------------------------------------
## Generates protein-level data with explicit pathway structure AND within-
## pathway protein correlation (compound symmetric via a one-factor model).
##
## This function only generates data; it does not run calibration.
## The package default significance engine is CAMERA/ROMER (ADR 0010).
## Within-pathway correlation inflates NET/REM relative to an
## abundance-stratified permutation null; that risk applies to the
## diagnostic permutation engine (`engine = "permutation"`), not to the
## default CAMERA/ROMER path.
##
## Noise model (compound symmetric within each pathway):
##   epsilon[i in pathway k] = sqrt(rho) * Z_k + sqrt(1-rho) * E_i
## where Z_k ~ N(0, noise_sd^2) is a pathway-shared factor and
## E_i ~ N(0, noise_sd^2) is protein-specific. This gives pairwise
## correlation rho within a pathway and 0 across pathways/background.
## rho = 0 reduces to the independent-noise case (backwards compatible).

## resolve a signal shape for a given pathway over the time grid.
## Returns a per-time effect value (length = length(all_times)).
.signal_shape <- function(signal_type, all_times, amplitude, n_members,
                           member_index) {
  span <- max(all_times) - min(all_times)
  t_norm <- (all_times - min(all_times)) / if (span > 0) span else 1
  ## monotonic ramp: 0 at baseline, amplitude at last timepoint
  ramp <- amplitude * t_norm
  ramp[1] <- 0  # baseline unaffected
  switch(signal_type,
    none = rep(0, length(all_times)),
    codirectional = ramp,
    mixed = {
      ## half the members go up, half go down; member_index is 1-based within
      ## the pathway, so the sign flips at the midpoint.
      sgn <- ifelse(member_index <= ceiling(n_members / 2), 1, -1)
      sgn * ramp
    },
    rep(0, length(all_times))
  )
}

#' Simulate a time-course experiment with within-pathway protein correlation
#'
#' Generates synthetic protein-level abundance with an explicit pathway
#' membership structure and compound-symmetric within-pathway protein
#' correlation, for validating the type-I error and power of
#' \code{\link{tps_calibrate_pathways}} under correlated pathway members.
#'
#' The within-pathway correlation is induced by a one-factor model:
#' \code{epsilon[i in k] = sqrt(rho)*Z_k + sqrt(1-rho)*E_i}, giving pairwise
#' correlation \code{rho} within a pathway and 0 across pathways/background.
#' \code{within_pathway_cor = 0} reduces to independent noise.
#'
#' @param n_pathways integer; number of structured pathways.
#' @param pathway_size integer; proteins per pathway (pathways are disjoint).
#' @param n_background integer; additional unstructured (no-signal) proteins.
#' @param conditions character vector of non-baseline conditions.
#' @param times numeric vector of non-baseline timepoints (baseline = 0 auto).
#' @param n_replicates integer; biological replicates per cell.
#' @param within_pathway_cor numeric in [0,1); pairwise correlation of proteins
#'   within each pathway. 0 = independent (backwards compatible).
#' @param pathway_signal character; \code{"none"} (type-I/FDR benchmark),
#'   \code{"codirectional"} (all members same-direction monotonic; tests NET),
#'   or \code{"mixed"} (half up / half down; tests REM, NET should be null).
#' @param signal_fraction numeric in [0,1]; fraction of pathways that carry the
#'   signal (default 0.5; the rest are null pathways).
#' @param amplitude numeric; signal magnitude on log2 scale.
#' @param noise_sd numeric; residual SD on log2 scale.
#' @param missingness numeric in [0,1]; fraction of MCAR missing values.
#' @param seed integer.
#' @return a list of class \code{tps_simulation_correlated} with:
#'   \code{abundance} (matrix), \code{col_data} (data.frame),
#'   \code{truth} (per-feature-condition signal plan), \code{term2gene}
#'   (data.frame with columns term/gene, ready for \code{tps_pathway_scores}),
#'   \code{config}.
#' @export
#' @examples
#' sim <- tps_simulate_correlated_pathways(n_pathways = 5, pathway_size = 8,
#'   n_background = 40, within_pathway_cor = 0.5,
#'   pathway_signal = "codirectional", seed = 1)
#' dim(sim$abundance)
#' head(sim$term2gene)
#' head(sim$truth)
tps_simulate_correlated_pathways <- function(n_pathways = 20L,
                                              pathway_size = 10L,
                                              n_background = 100L,
                                              conditions = c("A", "B", "C"),
                                              times = c(1, 3, 7),
                                              n_replicates = 3L,
                                              within_pathway_cor = 0,
                                              pathway_signal = c("none",
                                                                  "codirectional",
                                                                  "mixed"),
                                              signal_fraction = 0.5,
                                              amplitude = 2,
                                              noise_sd = 0.5,
                                              missingness = 0,
                                              seed = 1L) {
  pathway_signal <- match.arg(pathway_signal)
  arg_capture <- list(
    call = "tps_simulate_correlated_pathways", n_pathways = n_pathways,
    pathway_size = pathway_size, n_background = n_background,
    conditions = conditions, times = times, n_replicates = n_replicates,
    within_pathway_cor = within_pathway_cor, pathway_signal = pathway_signal,
    signal_fraction = signal_fraction, amplitude = amplitude,
    noise_sd = noise_sd, missingness = missingness, seed = seed,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  ## ---- validate scalar args ------------------------------------------------
  if (length(n_pathways) != 1L || is.na(n_pathways) || n_pathways < 1L ||
      n_pathways != as.integer(n_pathways)) {
    tps_abort("n_pathways must be one positive integer.")
  }
  if (length(pathway_size) != 1L || is.na(pathway_size) || pathway_size < 2L ||
      pathway_size != as.integer(pathway_size)) {
    tps_abort("pathway_size must be one integer >= 2.")
  }
  if (pathway_signal == "mixed" && pathway_size %% 2L != 0L) {
    tps_abort("pathway_signal='mixed' requires even pathway_size so that half the members go up and half down (got %d).", pathway_size)
  }
  if (length(n_background) != 1L || is.na(n_background) || n_background < 0L ||
      n_background != as.integer(n_background)) {
    tps_abort("n_background must be one non-negative integer.")
  }
  if (!is.character(conditions) || length(conditions) < 1L ||
      anyNA(conditions) || any(!nzchar(conditions)) ||
      anyDuplicated(conditions)) {
    tps_abort("conditions must be unique, non-empty character values.")
  }
  if (!is.numeric(times) || length(times) < 2L || any(!is.finite(times)) ||
      any(times <= 0) || anyDuplicated(times)) {
    tps_abort("times must contain at least 2 unique, positive finite values.")
  }
  times <- sort(times)
  if (length(n_replicates) != 1L || is.na(n_replicates) ||
      n_replicates < 1L || n_replicates != as.integer(n_replicates)) {
    tps_abort("n_replicates must be one positive integer.")
  }
  if (length(within_pathway_cor) != 1L || is.na(within_pathway_cor) ||
      !is.finite(within_pathway_cor) || within_pathway_cor < 0 ||
      within_pathway_cor >= 1) {
    tps_abort("within_pathway_cor must be in [0,1).")
  }
  if (length(signal_fraction) != 1L || is.na(signal_fraction) ||
      !is.finite(signal_fraction) || signal_fraction < 0 ||
      signal_fraction > 1) {
    tps_abort("signal_fraction must be in [0,1].")
  }
  if (length(amplitude) != 1L || is.na(amplitude) || !is.finite(amplitude) ||
      amplitude < 0) {
    tps_abort("amplitude must be one non-negative finite value.")
  }
  if (length(noise_sd) != 1L || is.na(noise_sd) || !is.finite(noise_sd) ||
      noise_sd < 0) {
    tps_abort("noise_sd must be one non-negative finite value.")
  }
  if (length(missingness) != 1L || is.na(missingness) ||
      !is.finite(missingness) || missingness < 0 || missingness > 1) {
    tps_abort("missingness must be in [0,1].")
  }
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
    tps_abort("seed must be one finite value.")
  }

  withr::with_seed(seed, {
    baseline_time <- 0
    all_times <- c(baseline_time, times)

    n_path_members <- n_pathways * pathway_size
    n_features <- n_path_members + n_background
    feat_ids <- sprintf("P%04d", seq_len(n_features))

    ## pathway membership: pathways are disjoint blocks of consecutive ids.
    ## Background proteins (tail) belong to no pathway.
    term2gene <- data.frame(
      term = rep(sprintf("PW%02d", seq_len(n_pathways)), each = pathway_size),
      gene = feat_ids[seq_len(n_path_members)],
      stringsAsFactors = FALSE)

    ## which pathways carry the signal (empty when signal is "none" or
    ## signal_fraction is 0)
    n_signal <- if (pathway_signal == "none") 0L else floor(signal_fraction * n_pathways)
    signal_pathways <- if (n_signal > 0L) sprintf("PW%02d", seq_len(n_signal)) else character(0)

    ## per-feature pathway assignment (NA for background)
    feat_pathway <- rep(NA_character_, n_features)
    feat_pathway[seq_len(n_path_members)] <- term2gene$term

    ## per-feature within-pathway index (for mixed-signal sign assignment)
    feat_member_idx <- rep(NA_integer_, n_features)
    feat_member_idx[seq_len(n_path_members)] <-
      rep(seq_len(pathway_size), n_pathways)

    ## truth plan: one row per feature per condition
    truth_rows <- list()
    for (i in seq_len(n_features)) {
      pw <- feat_pathway[i]
      is_signal <- !is.na(pw) && pw %in% signal_pathways &&
        pathway_signal != "none"
      for (cc in conditions) {
        truth_rows[[length(truth_rows) + 1L]] <- data.frame(
          feature_id = feat_ids[i], condition = cc,
          pathway = pw, is_signal_member = is_signal,
          signal_type = if (is_signal) pathway_signal else "none",
          member_index = feat_member_idx[i],
          stringsAsFactors = FALSE)
      }
    }
    truth <- do.call(rbind, truth_rows)

    ## baseline abundance
    mu0 <- stats::rnorm(n_features, mean = 18, sd = 4)

    ## sample metadata
    sample_rows <- list(); sid <- 0L
    for (r in seq_len(n_replicates)) {
      sid <- sid + 1L
      sample_rows[[sid]] <- data.frame(
        sample_id = sprintf("baseline_R%d", r), condition = "baseline",
        time = baseline_time, replicate = sprintf("R%d", r),
        is_baseline = TRUE, stringsAsFactors = FALSE)
    }
    for (cc in conditions) {
      for (tt in times) {
        for (r in seq_len(n_replicates)) {
          sid <- sid + 1L
          sample_rows[[sid]] <- data.frame(
            sample_id = sprintf("%s_t%s_R%d", cc, tt, r), condition = cc,
            time = tt, replicate = sprintf("R%d", r), is_baseline = FALSE,
            stringsAsFactors = FALSE)
        }
      }
    }
    col_data <- do.call(rbind, sample_rows)
    rownames(col_data) <- col_data$sample_id
    n_samples <- nrow(col_data)

    ## ---- correlated noise (compound symmetric within pathway) --------------
    rho <- within_pathway_cor
    sd_factor <- sqrt(rho)
    sd_specific <- sqrt(1 - rho)

    ## one pathway-shared factor draw per pathway per sample (matrix: pathways x samples)
    ## Background proteins and rho==0 use pure E_i noise.
    pathway_factors <- NULL
    if (rho > 0 && n_pathways >= 1L) {
      pathway_factors <- matrix(
        stats::rnorm(n_pathways * n_samples, sd = noise_sd),
        nrow = n_pathways, ncol = n_samples)
    }

    abundance <- matrix(NA_real_, nrow = n_features, ncol = n_samples,
                         dimnames = list(feat_ids, col_data$sample_id))

    for (i in seq_len(n_features)) {
      pw <- feat_pathway[i]
      pw_idx <- match(pw, sprintf("PW%02d", seq_len(n_pathways)))
      has_factor <- rho > 0 && !is.na(pw_idx)
      for (j in seq_len(n_samples)) {
        cd <- col_data[j, ]
        ## residual noise: factor + specific
        specific <- stats::rnorm(1, sd = noise_sd)
        factor_part <- if (has_factor) sd_factor * pathway_factors[pw_idx, j] else 0
        noise <- factor_part + sd_specific * specific
        if (cd$is_baseline) {
          abundance[i, j] <- mu0[i] + noise
        } else {
          cc <- cd$condition; tt <- cd$time
          tr <- truth[truth$feature_id == feat_ids[i] &
                      truth$condition == cc, ]
          eff_val <- 0
          if (isTRUE(tr$is_signal_member)) {
            shape <- .signal_shape(tr$signal_type, all_times, amplitude,
                                    pathway_size, tr$member_index)
            eff_val <- shape[match(tt, all_times)]
          }
          abundance[i, j] <- mu0[i] + eff_val + noise
        }
      }
    }

    ## inject MCAR missingness
    if (missingness > 0) {
      n_miss <- floor(missingness * length(abundance))
      miss_idx <- sample(seq_along(abundance), n_miss)
      abundance[miss_idx] <- NA
    }

    structure(
      list(
        abundance = abundance,
        col_data = col_data,
        truth = truth,
        term2gene = term2gene,
        config = list(n_pathways = n_pathways, pathway_size = pathway_size,
                      n_background = n_background, conditions = conditions,
                      times = times, n_replicates = n_replicates,
                      within_pathway_cor = within_pathway_cor,
                      pathway_signal = pathway_signal,
                      signal_fraction = signal_fraction,
                      signal_pathways = signal_pathways,
                      amplitude = amplitude, noise_sd = noise_sd,
                      missingness = missingness, seed = seed,
                      baseline_time = baseline_time),
        provenance = list(.make_provenance(arg_capture, list(),
                                            stage = "simulate_correlated"))
      ),
      class = c("tps_simulation_correlated", "tps_simulation")
    )
  })
}

## ---- methods ---------------------------------------------------------------
#' @export
print.tps_simulation_correlated <- function(x, ...) {
  cat("<timeProtState correlated-pathway simulation>\n")
  cat(sprintf("  features: %d   pathways: %d (size %d)   background: %d\n",
              nrow(x$abundance), x$config$n_pathways,
              x$config$pathway_size, x$config$n_background))
  cat(sprintf("  conditions: %s   times: %s   replicates: %d\n",
              paste(x$config$conditions, collapse = ", "),
              paste(x$config$times, collapse = ","),
              x$config$n_replicates))
  cat(sprintf("  within_pathway_cor: %.2f   signal: %s   amplitude: %.1f\n",
              x$config$within_pathway_cor, x$config$pathway_signal,
              x$config$amplitude))
  n_sig <- length(x$config$signal_pathways)
  cat(sprintf("  signal pathways: %d / %d\n", n_sig, x$config$n_pathways))
  invisible(x)
}
