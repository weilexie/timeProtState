# ADR 0006: Trajectory baseline modeled once

- **Date**: 2026-08-06
- **Decision ID**: D6
- **Status**: Approved

## Context

Legacy Mfuzz copied day0 baseline as t0=0 for each condition (3x inflation).

## Decision

- `tps_trajectory_matrix()` returns unique_cell (baseline once) for clustering.
- display matrix (baseline copied per condition) for plotting only.
- delta matrix: log2FC vs baseline.

## Evidence

- `R/trajectory-matrix.R`: three-view output.
- Test: `test-trajectory-cluster.R` asserts exactly 1 baseline column.
