# Release process

This checklist covers the `timeProtState` 0.99.x submission candidate. Do not
publish placeholder repository URLs, issue trackers, ORCIDs, or DOIs.

## Local quality gate

Run with R >= 4.5 and a working Pandoc on `PATH`:

```r
devtools::document("04_package/timeProtState")
devtools::test("04_package/timeProtState")
devtools::build("04_package/timeProtState", vignettes = TRUE)
```

Then check the source tarball:

```powershell
R CMD check --as-cran --no-manual timeProtState_0.99.0.tar.gz
Rscript -e "BiocCheck::BiocCheck('timeProtState_0.99.0.tar.gz')"
```

Verified on 2026-08-09:

- tests: 359 success expectations, 64 expected warnings, no failures/errors;
- `R CMD check`: 0 ERROR, 0 WARNING, 1 NOTE (`New submission`);
- `BiocCheck 1.44.2`: 1 ERROR, 0 WARNING, 8 NOTES. The error is external:
  the maintainer email is not registered on the Bioconductor Support Site.

The check is not submission-clean until the maintainer registers
`weile.xie@outlook.com` on the Support Site. The author should also subscribe to
Bioc-devel. Add an ORCID or `fnd` role only when the real metadata are supplied.

## Validation gate

Run the authoritative sequence documented in
`05_validation/02_timeProtState-paper/README.md`. Confirm that:

- the named endpoints remain Acid day 7 and PXD006262 WT 240 min;
- comparison methods use the same pathway mappings and correction contract;
- the generated CSV files match the numbers reported in the manuscript;
- Figure 1 is regenerated from current CSV files and visually inspected;
- no discovery-count comparison is described as sensitivity, specificity, or
  false-positive performance without pathway-level ground truth.

## Public release and archive

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
