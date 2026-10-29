timeProtState
================

# timeProtState

> Contrast-aware dynamic state analysis for time-course proteomics

`timeProtState` fits one shared baseline for multi-condition time-course
proteomics and keeps a pathway time-path label with the endpoint call.

## Installation

``` r
# remotes::install_local(".")
remotes::install_github("weilexie/timeProtState@v0.99.0")
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
