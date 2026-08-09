# ADR 0004: Full-rank heterogeneity contrasts + subspace invariance

- **Date**: 2026-08-05
- **Decision ID**: D4
- **Status**: Approved (with empirical finding)

## Context

Legacy pathway script used 9 pairwise contrasts (rank 6, 3 redundant).
TimeModel used 6 (full rank). Spec required "reference-invariance".

## Empirical finding

With complete feature-level observations, limma F-tests agree across reference
encodings up to floating-point tolerance. With missing values,
`limma::contrasts.fit()` applies an approximate covariance transformation and
the numerical F-statistics can differ even with `trend=FALSE, robust=FALSE`.
Column-subspace equivalence remains the design-level invariant.

## Decision

- Use 6-contrast full-rank heterogeneity family.
- Test subspace equivalence (`rank(cm1) == rank(cbind(cm1, cm2))`). Treat
  numerical F equality as expected only for complete feature-level data.

## Evidence

- `R/contrasts.R`: `.build_heterogeneity()` builds (C-1)*T contrasts.
- Test: `test-design-contrasts-fit.R` subspace invariance test.
