## ============================================================================
## Simulate time-course proteomics for benchmarking
## ============================================================================
##
## Implements STATISTICAL_CONTRACT.md section 9 and Goal 15.
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
.f_delayed <- function(t, amplitude = 1, onset_frac = 0.5) {
  onset <- min(t) + onset_frac * (max(t) - min(t))
  amplitude * as.numeric(t >= onset)
}
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
## returns a data.frame: feature_id, condition, effect_type, amplitude, sign
.make_plan <- function(n_features, conditions, effect_mix, seed) {
  ## NOTE: do NOT set.seed here — the caller (tps_simulate_timecourse) already
  ## set the seed in a local RNG stream. Re-seeding here would couple the plan
  ## RNG to the downstream noise RNG in an undocumented way.
  rows <- list()
  ## every feature has an entry per condition (null by default)
  eff_types <- names(effect_mix)
  eff_probs <- as.numeric(effect_mix)
  for (i in seq_len(n_features)) {
    for (cc in conditions) {
      et <- sample(eff_types, size = 1L, prob = eff_probs)
      amp <- switch(et,
        null = 0,
        shared_up = runif(1, 0.8, 2),
        shared_down = runif(1, 0.8, 2),
        condition_specific_up = runif(1, 0.8, 2),
        condition_specific_down = runif(1, 0.8, 2),
        opposite = runif(1, 0.8, 2),
        transient = runif(1, 0.8, 1.5),
        reversal = runif(1, 0.8, 1.5),
        runif(1, 0.5, 1.5))
      sgn <- switch(et,
        shared_up = , condition_specific_up = 1,
        shared_down = , condition_specific_down = -1,
        opposite = sample(c(-1, 1), 1),
        1)
      rows[[length(rows) + 1L]] <- data.frame(
        feature_id = sprintf("P%04d", i),
        condition = cc, effect_type = et,
        amplitude = amp, sign = sgn,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

## resolve an effect type + condition context to a per-time function value
.resolve_effect <- function(effect_type, condition, reference_condition, t, amp, sgn) {
  base <- switch(effect_type,
    null = .f_null(t),
    shared_up = .f_monotonic(t, amp),
    shared_down = -.f_monotonic(t, amp),
    transient = .f_transient(t, amp),
    reversal = .f_reversal(t, amp),
    condition_specific_up = if (condition == reference_condition) .f_monotonic(t, amp) else .f_null(t),
    condition_specific_down = if (condition == reference_condition) -.f_monotonic(t, amp) else .f_null(t),
    opposite = if (condition == reference_condition) sgn * .f_monotonic(t, amp) else -sgn * .f_monotonic(t, amp),
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
  if (n_features < 1L) tps_abort("n_features must be >= 1.")
  if (length(conditions) < 1L) tps_abort("Need at least 1 condition.")
  if (length(times) < 2L) tps_abort("Need at least 2 non-baseline timepoints.")
  if (n_replicates < 1L) tps_abort("n_replicates must be >= 1.")
  if (abs(sum(effect_mix) - 1) > 1e-6) {
    tps_abort("effect_mix probabilities must sum to 1 (got %.4f).", sum(effect_mix))
  }
  if (missingness < 0 || missingness > 1) tps_abort("missingness must be in [0,1].")

  ## local RNG stream (does not perturb global .Random.seed)
  rng_old <- .GlobalEnv$.Random.seed
  on.exit({
    if (is.null(rng_old)) rm(".Random.seed", envir = .GlobalEnv)
    else .GlobalEnv$.Random.seed <- rng_old
  }, add = TRUE)
  set.seed(seed)

  ## baseline time = 0 (relative); the shared baseline
  baseline_time <- 0
  all_times <- c(baseline_time, times)

  ## truth plan: per feature per condition
  ref_cond <- conditions[1L]
  plan <- .make_plan(n_features, conditions, effect_mix, seed)

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
        eff <- .resolve_effect(pr$effect_type, cc, ref_cond, all_times,
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
#' \dontrun{
#' sim <- tps_simulate_timecourse(n_features = 100, seed = 1)
#' sim2 <- tps_inject_pathway_effect(sim, tps_example_term2gene,
#'                                    term = "org20001",
#'                                    effect_type = "transient_burst",
#'                                    condition = "Acid")
#' }
#'
#' @export
tps_inject_pathway_effect <- function(sim, term2gene, term,
                                       effect_type = c("early_directional",
                                                       "gradual_remodeling",
                                                       "transient_burst",
                                                       "delayed_opposition"),
                                       condition, amplitude = 2) {
  effect_type <- match.arg(effect_type)
  if (!inherits(sim, "tps_simulation")) {
    tps_abort("sim must be a tps_simulation object.")
  }
  if (!condition %in% sim$config$conditions) {
    tps_abort("condition '%s' not in simulation conditions (%s).",
              condition, paste(sim$config$conditions, collapse = ", "))
  }
  members <- term2gene$gene[term2gene$term == term]
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
    gradual_remodeling = amplitude * seq(0, 1, length.out = length(all_times)),
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
      v[length(all_times)] <- amplitude  # only last timepoint, same direction
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
    amplitude = amplitude, sign = 1,
    stringsAsFactors = FALSE)
  sim$truth <- rbind(sim$truth, inj)

  sim
}
