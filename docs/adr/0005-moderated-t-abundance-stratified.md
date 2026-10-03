# ADR 0005: Moderated-t input + abundance-stratified permutation null

- **Date**: 2026-08-05
- **Decision ID**: D5
- **Status**: Historical as package default (superseded by ADR 0011).
  Permutation is diagnostic only (`engine="permutation"`). The decision
  text below is the 2026-08-05 record.

## Context

Legacy pathway score = median of row-z (no permutation). No calibration.
Spike showed unstratified permutation is anti-conservative (p<0.05 in 25-44%).

## Decision

- Pathway input statistic: moderated t (carries direction + uncertainty).
- Default null *(historical)*: abundance_stratified (4-bin equal-frequency).
  ADR 0011 superseded this as the package default; permutation is
  diagnostic only (`engine="permutation"`).
- Documented caveat: does not correct within-pathway correlation.

## Evidence

- `spike/pathway_remodeling_spike_stratified.R`: abundance_stratified best.
- `R/pathway-scores.R`: default type="moderated_t".
- `R/pathway-null.R`: default mode="abundance_stratified".
