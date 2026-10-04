timeProtState
================

# timeProtState

> Contrast-aware dynamic state analysis for time-course proteomics

`timeProtState` is an R/Bioconductor-style package for analyzing
**shared-baseline, multi-condition time-course proteomics** experiments.
It provides:

- **Shared-baseline contrast grammar**: automatic construction of
  estimable, full-rank contrast families from a saturated group-means
  design (baseline modeled once, not copied per condition).
- **Evidence-decomposed protein states**: temporal,
  condition-heterogeneity, and direction axes combined into auditable
  state labels.
- **Pathway summaries and significance**: NET (net direction), REM
  (remodeling magnitude), and COH (coherence) are summaries. Default
  `engine = "auto"` tests pathways with limma CAMERA (independent
  arrays) or ROMER (when a repeat-measures block is set). Reported
  columns are `p_directional` / `p_mixed`. Permutation+VIF is a
  diagnostic engine; `vif_cap` is not retuned.
- **Optional repeated measures**: `tps_fit(block, correlation)` passes a
  subject block to `limma::lmFit`; the consensus correlation is estimated
  with `limma::duplicateCorrelation()` when `correlation` is omitted.
  A block also switches pathway significance from CAMERA to ROMER.
- **Correlated-pathway simulation**: `tps_simulate_correlated_pathways()`
  generates known-truth data with compound-symmetric within-pathway
  correlation for type-I / power checks.

## Installation

``` r
# From a local clone/source directory:
# remotes::install_local(".")
# remotes::install_github("ORG/timeProtState")  # fill after the public repo exists
```

## Quick start

``` r
library(timeProtState)
data("tps_example", package = "timeProtState")

se <- tps_import(tps_example_abundance, tps_example_metadata,
                 tps_example_annotation)
design <- tps_make_design(se)
contrasts <- tps_make_contrasts(design, endpoint = 7)
fit <- tps_fit(se, design, contrasts)
states <- tps_classify_states(fit)

pws <- tps_pathway_scores(fit, tps_example_term2gene,
                          id_col = "locus_tag", family = "vs_baseline")
pws_c <- tps_calibrate_pathways(pws)
tps_classify_pathways(pws_c)
```

When the same biological unit is followed over time, pass `block` (a
`colData` column or a sample-aligned vector). Pathway significance then
uses ROMER instead of CAMERA (`engine = "auto"`). The package does not
treat a column named `replicate` as a block by default.

``` r
## Only when subject_id identifies the SAME biological unit over time.
# fit_rm <- tps_fit(se, design, contrasts, block = "subject_id")
```

See the vignettes for pathway analysis (`family = "vs_baseline"`),
visualization, and `tps_simulate_correlated_pathways()`.

## Public validation

Package examples and tests use the synthetic `tps_example` object (200
features × 30 samples). The publication benchmarks use only public data and
known-truth simulation:

- complete-null CAMERA/ROMER calibration —
  `scripts/01_engine_null_calibration.R`
- CAMERA/ROMER power —
  `scripts/02_engine_power.R`
- PXD065583 (multi-condition shared baseline; CAMERA, no block) —
  `scripts/03_PXD065583_three_contrasts.R`
- PXD005703 (donor block → ROMER) —
  `scripts/04_PXD005703_romer.R`
- PXD065583 universe thresholds —
  `scripts/05_PXD065583_universe_thresholds.R`
- permutation+VIF diagnostic —
  `scripts/06_permutation_vif_diagnostic.R`
- baseline counted once versus copied —
  `scripts/07_baseline_once_vs_copied.R`
- same-data tool scope —
  `scripts/08_same_data_tool_scope.R`
- piano, endpoint statistics only —
  `scripts/09_piano_same_endpoint.R`

These benchmark scripts belong to the manuscript reproducibility bundle and
are not distributed inside this standalone package repository. The accessions
above remain the stable public-data identifiers; use the released manuscript
bundle to reproduce the benchmark sequence.

## Scope

The package starts from a **protein-level abundance matrix** (already
quantified and log2-transformed). It does **not** process raw
mass-spectrometry files, perform peptide-to-protein summarization, or
impute by default.

## License

Artistic-2.0
