# Release process

This checklist covers the `timeProtState` 0.99.x development quality gate.
Do not publish placeholder repository URLs, issue trackers, ORCIDs, or
DOIs. Do not describe the package as submission-ready.

## Local quality gate

Run with R >= 4.5 and a working Pandoc on `PATH`:

```r
devtools::document("04_package/timeProtState")
devtools::test("04_package/timeProtState")
devtools::build("04_package/timeProtState", vignettes = TRUE)
```

Then check the source tarball:

```powershell
R CMD check --as-cran timeProtState_0.99.0.tar.gz
Rscript -e "BiocCheck::BiocCheck('timeProtState_0.99.0.tar.gz')"
```

Verified on 2026-08-13 (see `02_audit_and_contract/05_STATUS.md` §3.3):

- tests: 0 failures, 0 errors (expected warnings from partial-NA coefficients);
- `R CMD check --as-cran`: 0 ERROR, 0 WARNING, 1 NOTE (`New submission`);
- `BiocCheck`: 1 ERROR, 0 WARNING, 8 NOTES. The error is external: the
  maintainer email is not registered on the Bioconductor Support Site.

The Bioconductor Support Site registration of `weile.xie@outlook.com` is
a human-only step. Add an ORCID or `fnd` role only when the real metadata
are supplied.

## Validation gate

Run the **public** sequence documented in
`05_validation/02_timeProtState-paper/README.md`:

- `scripts/04_simulation_correlation_benchmark.R`
- `scripts/05_PXD006262_closest_longitudinal.R`
- `scripts/06_PXD065583_ASO_benchmark.R`
- `scripts/07_trajectory_recovery_benchmark.R`

Confirm that:

- the named public endpoints remain PXD006262 WT 240 min and PXD065583
  ASO 72 h versus the shared mock baseline;
- comparison methods use the same pathway mappings and correction contract;
- NET is compared with CAMERA directional and REM with CAMERA mixed;
  self-contained methods are not ranked by discovery count;
- the generated public CSVs (`results/04–07`, `08–13`, `14–18`) match
  the numbers reported in STATUS;
- no discovery-count comparison is described as sensitivity, specificity,
  or superiority without pathway-level ground truth.

Do **not** run `04_reproduce_dormancy_dynamics.R`,
`07_comprehensive_benchmark.R`, or `12_render_figure1.R` as part of this
gate. Those scripts and the Acid day-7 / Figure 1 artefacts are frozen
historical leftovers.

## Public release and archive

Leave these steps until a real public tracker and manuscript rewrite exist.
Do not invent repository URLs or a Zenodo DOI beforehand.

1. Create the public repository and issue tracker.
2. Add the real repository and issue URLs to `DESCRIPTION`.
3. Commit a source tarball that has passed the local quality gate.
4. Create the version tag and public release.
5. Archive that release in Zenodo and obtain the DOI.
6. Synchronize the DOI and release metadata across `inst/CITATION`,
   `CITATION.cff`, README, and the manuscript.
7. Rebuild and repeat both checks after metadata changes.

Follow the current Bioconductor submission instructions at submission time;
do not assume that a historical Git remote or review workflow is still valid.

## Version scheme

```text
0.99.x  Bioconductor review candidates
1.0.0   First accepted release
1.0.x   Release-branch fixes
```
