# ADR 0007: trend=TRUE as fixed default (not "auto")

- **Date**: 2026-08-06
- **Decision ID**: D7
- **Status**: Approved (deviation from original D7 proposal noted)

## Context

Original D7 proposed `trend="auto"` (data-driven loess). Implementation chose
fixed `trend=TRUE` for simplicity and because proteomics intensity-variance
dependence is near-universal.

## Decision

Fixed `trend=TRUE` default. The `trend` parameter is exposed for override.

## Rationale

- "auto" adds implementation complexity (loess significance test) for marginal
  benefit when trend=TRUE is almost always correct for proteomics.
- Users who need data-driven trend can call `tps_fit(trend=FALSE)` and inspect
  plotSA manually.

## Evidence

- `R/fit.R`: `tps_fit(trend=TRUE)`.
