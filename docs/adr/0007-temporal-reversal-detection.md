# ADR 0007: Temporal reversal detection with effect-size gate

- **Date**: 2026-08-06
- **Decision ID**: D8
- **Status**: Approved

## Context

Spec §6.3/§6.4 requires `temporal_reversal` as a state label. Legacy used
discrete Up/Down/NS co-occurrence (no effect-size gate → noise judged as
opposite). D8 initially marked resolved but classifier lacked the path.

## Decision

- Detect reversal from vs_baseline family: within a condition, if >=2
  significant timepoints have opposite logFC signs (with effect-size gate
  via `tps_decision_rule$lfc_threshold`), label as `temporal_reversal`.
- The `is_reversal` column is retained for audit.

## Evidence

- `R/classify-states.R`: `.detect_reversal()` + `.combine_state()` includes
  `temporal_reversal` branch.
- Test: `test-classify-states.R` verifies reversal detection on simulated data.
