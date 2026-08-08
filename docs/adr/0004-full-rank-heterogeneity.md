# ADR 0004: Full-rank heterogeneity contrasts + subspace invariance

- **Date**: 2026-08-05
- **Decision ID**: D4
- **Status**: Approved (with empirical finding)

## Context

Legacy pathway script used 9 pairwise contrasts (rank 6, 3 redundant).
TimeModel used 6 (full rank). Spec required "reference-invariance".

## Empirical finding

limma moderated F/t are NOT numerically invariant under re-parameterization
(even with trend=FALSE, robust=FALSE), because eBayes variance moderation
depends on residual structure. Only column-subspace equivalence is provable.

## Decision

- Use 6-contrast full-rank heterogeneity family.
- Test subspace equivalence (`rank(cm1) == rank(cbind(cm1, cm2))`), not
  numerical F equality.

## Evidence

- `R/contrasts.R`: `.build_heterogeneity()` builds (C-1)*T contrasts.
- Test: `test-design-contrasts-fit.R` subspace invariance test.
