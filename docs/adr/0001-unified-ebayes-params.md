# ADR 0001: Unified eBayes parameters (trend=TRUE, robust=TRUE)

- **Date**: 2026-08-05
- **Decision ID**: D1
- **Status**: Approved

## Context

The legacy M. smegmatis analysis used three different eBayes parameter sets:
- Main limma: trend=TRUE, robust=TRUE
- TimeModel: trend=FALSE, robust=FALSE (defaults)
- Pathway: trend=FALSE, robust=TRUE

This caused the same protein's moderated t/P to differ across scripts.

## Decision

Adopt `trend=TRUE, robust=TRUE` uniformly for all limma fits in the package.

## Alternatives considered

- (b) Use defaults everywhere (trend=FALSE) — less appropriate for proteomics
  where intensity-variance dependence is common.
- (c) Data-driven trend decision — added complexity, deferred.

## Consequences

- All contrast families use consistent moderation.
- `tps_fit(trend=, robust=)` exposes the parameters for override.
- Legacy profile reproduces via explicit parameter setting.

## Evidence

- `R/fit.R`: `tps_fit()` default `trend=TRUE, robust=TRUE`.
- Test: `test-design-contrasts-fit.R` asserts defaults.
