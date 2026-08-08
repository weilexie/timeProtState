# ADR 0002: FDR default significance rule

- **Date**: 2026-08-05
- **Decision ID**: D2
- **Status**: Approved

## Context

Legacy used raw P.Value < 0.05 + |log2FC| >= 0.74 (Q95 noise cutoff).
This is project-specific and not appropriate as a software default.

## Decision

Default: `adj.P.Val < 0.05`, `lfc_threshold = 0`.
Legacy profile via `tps_decision_rule(significance="P.Value", lfc_threshold=0.74)`.

## Consequences

- Default is more conservative (fewer hits) than legacy raw P.
- Users can reproduce legacy exactly with explicit rule.

## Evidence

- `R/classify-states.R`: `tps_decision_rule()` defaults.
- Test: `test-classify-states.R` asserts FDR default + legacy profile.
