# ADR 0010: CAMERA/ROMER as default pathway significance

- **Date**: 2026-09-15
- **Status**: Approved

## Decision

Pathway **significance** uses limma competitive gene-set tests:

- `engine = "auto"` (default): no `block` → `limma::camera`; `block` present → `limma::romer`.
- Explicit `engine = "camera"` aborts if the fit has a block.
- `engine = "romer"` is allowed with or without block; with block the original aligned block vector and consensus correlation must be passed through.
- `engine = "fry"` and `engine = "permutation"` are explicit only. Auto never switches engine because of missing values.

NET = mean(t), REM = mean(|t|), COH = |NET|/(REM+ε) are **summaries**. Significance columns are `p_directional` / `q_directional` and `p_mixed` / `q_mixed`. COH has no p-value.

CAMERA: `sort = FALSE`, `inter.gene.cor = NA`, `trend.var` from the fit. Directional uses `directional = TRUE` (parametric); mixed uses `directional = FALSE`, `use.ranks = TRUE`. Do not call `cameraPR`. Match results by pathway name.

ROMER: `p_mixed = P.Mixed`; `p_directional = min(1, 2 * min(P.Up, P.Down))` (package convention). Store `p_up`, `p_down`. Default `nrot = 9999`.

Universe is formed once in `tps_pathway_scores(tps_fit)`: features with finite statistics on the selected family **and** complete expression on inference samples. Calibration does not change the universe or recompute NET/REM/COH.

`tps_calibrate_dynamics()` no longer permutes. `p_traj` is Bonferroni-minP of the time-point directional and mixed p-values; `q_dynamic` is one global BH over pathway × condition. Shape labels remain descriptive rules using per-time `q_directional` / `q_mixed`.

Permutation+VIF (ADR 0009) remains a diagnostic engine only. Do not retune `vif_cap`.

## Consequences

- This is not a claim of a new competitive test, and not a claim of superiority versus CAMERA.
- PXD065583 (no block, CAMERA) and PXD005703 (block, ROMER) are different engines; record `engine_used`.
