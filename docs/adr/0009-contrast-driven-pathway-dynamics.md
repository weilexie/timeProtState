# ADR 0009: Contrast-driven pathway remodeling dynamics

- **Date**: 2026-08-06
- **Decision ID**: Reframe (paths 1+2+3 integration)
- **Status**: Approved

## Context

The original pathway module (ADR 0005) computed NET/REM/COH for a single
endpoint contrast — a static snapshot. This was the weakest methodological
position vs TiDEomics (Bioconductor 3.24), which does time-course DE +
WGCNA/Trendy trajectories + standard enrichment. An endpoint-only pathway
score cannot answer "when does this pathway start remodeling?" or "do
conditions diverge in pathway direction over time?"

## Key insight (code audit)

`tps_pathway_scores` and `tps_calibrate_pathways` are already
contrast-agnostic — NET/REM/COH, permutation null, and BH q are all computed
per-contrast (per-column of the statistic matrix). Passing
`family = "vs_baseline"` instead of `family = "endpoint"` immediately yields
scores for all (condition, timepoint) pairs with **zero changes to the
pathway computation code**.

## Decision

1. **Reframe** the three contributions from "parallel" to "integrated":
   - Contrast grammar (contribution 1) = infrastructure layer
   - Evidence states (contribution 2) = protein-level intermediate
   - Pathway remodeling dynamics (contribution 3) = time-course output
   - Core novelty = the integration + time-course-aware pathway inference

2. **New functions**:
   - `tps_classify_pathway_dynamics()`: classifies NET/REM/COH trajectories
     into 7 dynamic patterns (early/gradual/delayed/transient/divergent/...)
   - `tps_link_states_pathways()`: connects protein states to pathway dynamics
   - `tps_plot_pathway_trajectory()`: visualizes pathway NET/REM/COH over time

3. **Vignettes** updated to demonstrate `family = "vs_baseline"` (time-course)
   as the primary pathway analysis mode; endpoint noted as a special case.

## Differentiation vs TiDEomics

TiDEomics does protein-level pairwise DE then standard enrichment (two-step,
loses pathway direction structure). timeProtState uses rank-audited contrast
families' moderated statistics to drive pathway-level inference directly —
one integrated step that preserves direction/uncertainty and is calibrated
per-timepoint. TiDEomics has no contrast grammar, no NET/REM/COH, no
permutation calibration.

## Evidence

- `R/pathway-dynamics.R`: tps_classify_pathway_dynamics + tps_link_states_pathways
- `R/plotting.R`: tps_plot_pathway_trajectory
- `tests/testthat/test-pathway-dynamics.R`: 22 tests
- Existing golden-master tests (endpoint, user_matrix) unaffected
