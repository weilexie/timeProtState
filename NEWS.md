# timeProtState 0.99.0

## Behavior changes (2026-08-25)

- `tps_calibrate_pathways()`: default `vif_cap` changed from `Inf` to
  `3` (ADR 0010). The uncapped default erased nearly all real-data signal
  (0 significant NET pathways on PXD006262; 0% dynamic trajectories on
  PXD065583). Validated operating point (30 seeds, 200 perm): rho=0.2
  type-I 0.023 / power 0.713; rho=0.5 type-I 0.073 / power 0.500
  (vs 0.010/0.040 uncapped). Use `vif_cap = Inf` for maximal type-I
  safety or `vif_cap = 1` to disable VIF scaling.
- `tps_calibrate_pathways()`: the per-pathway residual rho used by
  `intergene = "vif"` is now debiased. The clipped estimator
  `max(0, mean pairwise rho)` is upward-biased at true rho = 0, which
  inflated every pathway's VIF above 1 even on independent data (ADR 0010,
  2026-08-25). The bias floor is estimated from 50 size-matched random
  gene sets drawn from the same residuals and subtracted. Simulation
  (rho=0): power restored to 1.000 at type-I 0.010. PXD065583: NET q<0.05
  44-51 per ASO (was 14-22), Jaccard vs CAMERA 0.55-0.68.
- `tps_calibrate_dynamics()`: `q_dynamic` BH adjustment changed from
  joint-across-conditions to per-condition families. The joint penalty
  scaled with the number of conditions in the object, so a single-condition
  negative-control arm appeared systematically more dynamic than a
  many-condition arm purely from BH aggregation size.

## Implemented API

- Data contract: `tps_import()`, `tps_validate()`, `tps_metadata_schema()`,
  `tps_provenance()` with shared-baseline detection and ID-mapping audit.
- Design: `tps_make_design()` saturated group-means model with rank audit.
- Contrasts: `tps_make_contrasts()` six estimable families (vs-baseline,
  within-condition, average-time, heterogeneity, within-omnibus, endpoint).
- Fit: `tps_fit()` unified limma (`trend=TRUE`, `robust=TRUE`) with optional
  repeated-measures `block` / `correlation` (`limma::duplicateCorrelation`
  when `block` is supplied and `correlation` is left `NULL`). Accessors:
  `tps_results()`, `tps_statistic_matrix()`.
- States: `tps_classify_states()` three evidence axes + rule manifest;
  `tps_decision_rule()` and `tps_state_table()`. The historical raw-P +
  `|logFC| >= 0.74` rule is reproduced by explicit arguments, not a named
  profile object.
- Trajectory: `tps_trajectory_matrix()` unique-cell / display / delta views
  (baseline modeled once).
- Clustering: `tps_cluster_trajectories()` Mfuzz/kmeans adapter with QC;
  `tps_evaluate_k()`.
- Pathway snapshots: `tps_pathway_scores()` NET/REM/COH;
  `tps_calibrate_pathways()` abundance-stratified permutation with a
  default inter-gene VIF adjustment (`intergene = "vif"`) so competitive
  NET p-values account for member correlation and REM uses the same correction
  as a documented first-order approximation (ADR 0010;
  `intergene = "none"` is the historical permutation);
  `tps_classify_pathways()`.
- Pathway dynamics: `tps_classify_pathway_dynamics()`,
  `tps_calibrate_dynamics()`, `tps_divergence_score()`,
  `tps_link_states_pathways()`, `tps_plot_pathway_trajectory()`.
- Visualization: `tps_plot_volcano()`, `tps_plot_trajectories()`,
  `tps_plot_state_distribution()`, `tps_plot_pathway_states()`.
- Simulation: `tps_simulate_timecourse()`,
  `tps_simulate_correlated_pathways()` (within-pathway factor-model
  correlation), `tps_inject_pathway_effect()`.
- Export: `tps_export()` structured directory output with manifest.

## Changes (2026-08-17, Rounds 2-5 deep audit)

### Bug fixes

- **REM p-value centering** (`tps_calibrate_pathways`): REM now uses the same
  centered second-pass permutation comparison as NET (was uncentered first-pass).
  Both `pmax(..., 0)` clamps applied symmetrically for the one-sided upper-tail test.
- **State classification** (`tps_classify_states`): new labels `endpoint_up`,
  `endpoint_down`, `temporal_only`, `heterogeneity_only` replace the ambiguous
  `endpoint_only`. Proteins with temporal F but no endpoint direction are no
  longer mislabeled as `stable_or_unresolved`.
- **Rank-deficiency check** (`tps_make_contrasts`): now aborts if the design
  matrix is rank-deficient (e.g., batch confounding). Guarded with `!is.null()`
  for backward compatibility with designs lacking `rank_audit`.
- **Mixed signal parity** (`tps_simulate_correlated_pathways`): rejects odd
  `pathway_size` when `pathway_signal = "mixed"` (requires balanced +/− split).
- **All-NA sample detection** (`tps_validate`): now flags all-NA sample columns
  as an ERROR (previously only checked row-wise all-NA proteins).
- **Negative time warning** (`tps_validate`): warns when time values are negative.
- **`fold_msg` formatting** (`tps_validate`): fixes mangled error messages when
  all arguments are character (the `is.character()` heuristic is replaced with
  first-arg-as-template).
- **Volcano plot logFC=0** (`tps_plot_volcano`): zero-logFC significant features
  now labeled `ns` instead of `down`.
- **Pathway trajectory parser** (`tps_plot_pathway_trajectory`): strips
  `endpoint_` prefix from contrast names before parsing condition/time.
- **Kmeans distinct-row check** (`tps_cluster_trajectories`): aborts with a
  clear message when distinct rows < k (was an opaque `stats::kmeans` crash).
- **Row standardization NA handling** (`tps_cluster_trajectories`): rows with
  any NA are now excluded from clustering (was silently imputed to 0, violating
  the package's no-imputation contract).
- **`table()` NA handling** (plotting): `useNA = "ifany"` added to prevent
  silent NA drops in state distribution plots.
- **Provenance preservation** (`tps_calibrate_dynamics`): classify step
  provenance is now preserved as a list entry instead of being overwritten.
- **Dead code** (`tps_checksum`): removed unused `serialize()` call.

### New features

- **`vif_cap` parameter** (`tps_calibrate_pathways`): upper bound on VIF to
  prevent over-correction in large pathways with high estimated rho. Default
  `Inf` (no cap); empirical analysis shows cap=3 controls type-I to ~0.055
  while retaining ~59% power at rho=0.5, but the optimal cap is data-dependent.
  (Superseded 2026-08-25: the default is now `3` — see "Behavior changes"
  above.)

### Test improvements

- Block fit end-to-end integration test (fit → pathway_scores → calibrate → VIF).
- `p_rem` golden-master determinism test (locks down Round 2 REM centering fix).
- Golden-master tests for `intergene = "none"` backward compatibility.
- VIF edge-case tests (pathway size 3, bare matrix, REM type-I, vif_cap bounds).
- Trajectory recovery benchmark (50 seeds, 4 effect types, 98-100% recovery).
