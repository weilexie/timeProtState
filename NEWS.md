# timeProtState 0.99.0

## Initial Bioconductor submission candidate

- Data contract: `tps_import()`, `tps_validate()`, `tps_metadata_schema()`,
  `tps_provenance()` with shared-baseline detection and ID-mapping audit.
- Design: `tps_make_design()` saturated group-means model with rank audit.
- Contrasts: `tps_make_contrasts()` six estimable families (vs-baseline,
  within-condition, average-time, heterogeneity, within-omnibus, endpoint).
- Fit: `tps_fit()` unified limma (trend=TRUE, robust=TRUE) with
  `tps_results()` and `tps_statistic_matrix()` accessors.
- States: `tps_classify_states()` three evidence axes + rule manifest;
  `tps_decision_rule()` with legacy profile support.
- Trajectory: `tps_trajectory_matrix()` unique-cell / display / delta views
  (baseline once — fixes legacy inflation).
- Clustering: `tps_cluster_trajectories()` Mfuzz/kmeans adapter with QC;
  `tps_evaluate_k()`.
- Pathway: `tps_pathway_scores()` NET/REM/COH; `tps_calibrate_pathways()`
  abundance-stratified permutation; `tps_classify_pathways()`.
- Visualization: `tps_plot_volcano()`, `tps_plot_trajectories()`,
  `tps_plot_state_distribution()`, `tps_plot_pathway_states()`.
- Simulation: `tps_simulate_timecourse()` for benchmarks.
- Export: `tps_export()` structured directory output with manifest.
