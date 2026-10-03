# ADR 0010: Inter-gene VIF adjustment of the competitive permutation null

> **ADR 0011 (2026-09-15) makes permutation+VIF diagnostic only
> (`engine="permutation"`).** `vif_cap=3` remains the diagnostic cap, not
> the default significance engine. Method description below is unchanged.

- **Date**: 2026-08-15
- **Decision ID**: D11
- **Status**: Approved (diagnostic engine; package default is ADR 0011)

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

## vif_cap parameter (added 2026-08-17; default changed 2026-08-25)

The compound-symmetric VIF `1 + (k-1)*rho` grows linearly with pathway size,
which can over-penalize real signal in large pathways. The `vif_cap` parameter
limits the VIF to prevent over-correction.

Empirical analysis on real proteomics data shows highly variable VIF:
- PXD006262 (yeast): median VIF 5.43, 82% of pathways exceed 3 (rho median 0.371)
- PXD065583 (human): median VIF 1.66, 12% exceed 3 (rho median 0.022)

Simulation (rho=0.5, k=10, 20 seeds, 200 perm): cap=3 gives type-I 0.055 +
power 0.590; uncapped gives type-I 0.000 + power 0.060.

### rho clipping-bias debias (added 2026-08-25)

**Problem.** The per-pathway rho estimator clips at zero
(`max(0, mean pairwise residual correlation)`). At true rho = 0 half the
estimation noise mass is floored while the other half stays positive, so
`E[rho_hat | rho = 0] > 0` (approximately +SE/2, ~0.03-0.05 at proteomics
sample counts). The VIF multiplies this bias by (k-1), so independent data
acquired a spurious median VIF of ~1.5 and EVERY pathway's observed
departure was systematically shrunk. This was diagnosed as the dominant
cause of the real-data discovery deficit on low-correlation data
(PXD065583: 14-22 NET at q<0.05 vs CAMERA's 68-78, despite both methods
having power 1.0 and controlled type-I at rho=0 in simulation).

**Fix.** The bias floor is estimated directly from the data: draw 50
size-matched random gene sets from the same residual matrix (membership
permutation is exactly the null in question), compute each draw's mean
pairwise rho, average per pathway size, and subtract:
`rho_debiased = max(0, rho_pathway - floor(size))`. Fixed seed for
reproducibility.

**Validated outcome (30 seeds, 200 perm, cap=3):**

| rho | type-I | power (before debias) | power (after) |
|---|---|---|---|
| 0 | 0.010 | (deficit existed only in the estimator bias) | 1.000 |
| 0.2 | 0.023 | 0.713 | 0.747 |
| 0.5 | 0.073 | 0.500 | 0.500 (floor negligible at rho 0.5) |

**Real data:** PXD065583 NET q<0.05 rose from 14-22 to 44-51 per ASO
(CAMERA 68-78; Jaccard 0.55-0.68, Spearman 0.92-0.94). PXD006262
q<0.1 went 0 -> 2 with NET-vs-CAMERA Spearman 0.26 -> 0.83.

### Historical 2026-08-25 decision (superseded 2026-09-15)

**Decision.** On 2026-08-25 `vif_cap = 3` became the permutation-engine
cap; ADR 0011 (2026-09-15) made permutation+VIF diagnostic only
(`engine = "permutation"`). It is not the default significance engine.
CAMERA/ROMER do not use `vif_cap`. `vif_cap = Inf` remains available
for maximal type-I safety; `vif_cap = 1` disables VIF scaling.

**Why.** The uncapped default made the tool unusable on real data: on
PXD006262 the default calibration found 0 significant NET pathways at the
endpoint (vs 6 under `intergene="none"`), and the trajectory-level
`q_dynamic` was 0% everywhere on PXD065583 (max-statistic x joint BH x
uncapped VIF stack up). A default that returns nothing on real data is not
an acceptable operating point for a shipped tool, regardless of how safe it
is in type-I terms. The 20-seed simulation identified cap=3 as the balanced
operating point (type-I 0.055 at nominal 0.05, power 0.590 vs 0.060
uncapped); real-data VIF distributions show the cap binds exactly where the
noise-amplified residual rho estimate is most likely to be overstated
(large pathways, high replicate-noise).

**Mechanism note.** The residual-based rho estimate is noisy at proteomics
replicate counts (3-4 per cell), and the VIF multiplies it by (k-1), so
estimation noise is amplified linearly in pathway size. The cap bounds this
amplification. On low-correlation data (PXD065583) the cap is rarely
binding; on high-correlation data (PXD006262) it truncates most pathways,
which is precisely where the uncapped estimate is least trustworthy.

**Validation obligations attached to this change:**
1. Focused simulation re-run confirming the cap=3 operating point (type-I /
   power) beyond the original 20 seeds.
2. Real-data benchmark re-runs under the historical permutation default
   (scripts 05, 06, 10; later superseded by ADR 0011)
   with all CSVs refreshed and parameters recorded.
3. Documentation states the trade-off and the escape hatches
   (`vif_cap = Inf`, `intergene = "none"`).

**Predecessor decision.** The 2026-08-15 decision to default to
`vif_cap = Inf` ("let users choose based on their rho distribution") is
superseded: shipping a default that erases real-data signal pushes the
calibration burden onto users who have not yet looked at their VIF
distribution, which is the population that most needs a sane default.

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
