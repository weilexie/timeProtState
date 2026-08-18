# ADR 0010: Inter-gene VIF adjustment of the competitive permutation null

- **Date**: 2026-08-15
- **Decision ID**: D11
- **Status**: Approved

## Context

ADR 0005 and STATUS §15 documented that abundance-stratified label
permutation does not preserve within-pathway protein correlation. At
compound-symmetric `rho = 0.5`, timeProtState NET type-I rose to ~0.32
and CAMERA directional rose to ~0.31. STATUS §15.4 left the null
unchanged because the inflation was shared with the same-null comparator.

That decision blocked a claim of better *calibration*. The remaining
reviewer objection is exactly this: a competitive directional test that
assumes independent members is anti-conservative when members are
correlated. CAMERA's own inter-gene correlation / VIF estimate does not
fully correct the same design (few proteomics replicates).

NET for pathway `g` and contrast `j` is the mean of member moderated t.
If members have pairwise correlation `rho`,

```
Var(NET) = (sigma^2 / k) * (1 + (k - 1) * rho)
```

Label permutation yields `rho ≈ 0` in the null, so the permutation
scale is too small. The missing factor is the variance inflation
`VIF = 1 + (k - 1) * rho`.

`rho` is estimated from **model residuals across samples** (expression
minus the fitted group-means stored by `tps_fit()`). That residual
correlation is the extra member covariance not explained by the design.
Estimating `rho` from moderated-t profiles would treat a shared pathway
effect as noise correlation and over-correct true NET.

## Decision

1. `tps_calibrate_pathways()` gains `intergene = c("vif", "none")`.
2. **Default is `vif`**: scale the observed NET (and REM excess) by
   `1/sqrt(VIF)` before comparing with the existing permutation stream
   (`+1` smoothing, centered NET and REM second-pass; see "REM centering"
   section below for the Round 2 change from ADR 0005's uncentered REM).
   **Note on REM**: the VIF derivation above is exact for NET (mean of
   correlated normals). REM is the mean of |t| (folded normals), whose
   variance under correlated members has a different structure. We apply
   the NET-derived VIF to REM as a first-order approximation; this is
   reasonable because the dominant effect of correlation is the same
   variance inflation, but it is not a derived identity.
3. `intergene = "none"` reproduces the ADR 0005 NET p-values exactly
   (golden-master path). REM p-values differ from ADR 0005 because the
   Round 2 fix applied centered REM in both modes.
4. `rho` is truncated to `[0, 0.99]`. Negative estimates are treated as
   0. VIF is 1 when residuals are missing, have fewer than 4 samples,
   or a pathway has fewer than 3 finite members.
5. NET/REM/COH *formulas* are unchanged. Only the competitive p-value
   scale changes.
6. `tps_calibrate_dynamics()` uses the stored per-pathway VIF when
   standardizing NET/REM against the parent null SD.

## Consequences

- Independent-member simulations (`rho = 0`) stay type-I controlled;
  estimated `rho` is near 0 and VIF ≈ 1.
- Correlated-member nulls become substantially less anti-conservative
  than unadjusted permutation and than CAMERA on the same design.
- Power at large amplitude remains, but weaker correlated signals are
  harder to call — that is the cost of an honest competitive null.
- Publication benchmarks that need the historical permutation must pass
  `intergene = "none"`.

## vif_cap parameter (added 2026-08-17)

The compound-symmetric VIF `1 + (k-1)*rho` grows linearly with pathway size,
which can over-penalize real signal in large pathways. The `vif_cap` parameter
(default `Inf`, no cap) limits the VIF to prevent over-correction.

Empirical analysis on real proteomics data shows highly variable VIF:
- PXD006262 (yeast): median VIF 5.43, 82% of pathways exceed 3 (rho median 0.371)
- PXD065583 (human): median VIF 1.66, 12% exceed 3 (rho median 0.022)

Simulation (rho=0.5, k=10): cap=3 gives type-I 0.055 + power 0.590;
uncapped gives type-I 0.000 + power 0.060. The optimal cap is data-dependent;
users should inspect the `$vif` field and choose based on their data's rho distribution.

## REM centering (added 2026-08-17, Round 2 fix)

REM previously used the first-pass uncentered comparison (`rem_b >= obs_rem`),
inconsistent with NET's centered second-pass (`|net_b - mean| >= |obs - mean|`).
REM now uses the same centered second-pass with `pmax(..., 0)` on both observed
and permuted excess, ensuring one-sided upper-tail behavior is correct:
pathways with REM below null mean get p ≈ 1.0 (clearly null), not p ≈ 0.5.

## Evidence

- `R/pathway-null.R`: `.estimate_pathway_rho()`, VIF scaling, `vif_cap`.
- `tests/testthat/test-calibrate-intergene.R`.
- `tests/testthat/test-block-e2e.R` (block fit → VIF chain).
