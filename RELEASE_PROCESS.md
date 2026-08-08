# Release Process

This document describes the release workflow for `timeProtState` to support
Zenodo DOI minting and Bioconductor submission.

## Version scheme

```
0.99.0   Bioconductor submission candidate (devel)
1.0.0    First Bioconductor release (after review acceptance)
1.0.x    Patch releases (review fixes)
```

## Pre-release checklist

1. **Tests**: `devtools::test()` → 0 failures (currently 231/231 PASS).
2. **R CMD check**: `R CMD check --as-cran` → 0 ERROR, 0 WARNING (vignette
   compilation requires pandoc; CI provides it). Currently 2 WARNINGs are
   vignette-HTML-only (pandoc environment) and 1 NOTE (e1071 library call,
   documented).
3. **BiocCheck**: `BiocCheck::BiocCheck(".")` — run on Linux CI.
4. **Vignettes**: both build with pandoc installed.
5. **Documentation**: all exported functions have runnable examples.
6. **LICENSE**: Artistic-2.0, no commercial restrictions.

## Git release + Zenodo DOI

1. Tag the release:
   ```bash
   git tag -a v0.99.0 -m "Bioconductor submission candidate"
   git push origin v0.99.0
   ```
2. On GitHub: create a Release from the tag. Include:
   - source tarball (`timeProtState_0.99.0.tar.gz`)
   - NEWS.md excerpt
   - checksums
3. Zenodo: enable GitHub integration (Settings → Webhooks → Zenodo). Pushing
   the tag triggers Zenodo to mint a DOI automatically. The `.zenodo.json`
   file controls metadata.
4. Record the DOI in:
   - `inst/CITATION`
   - `CITATION.cff`
   - README.md
   - manuscript Data Availability section

## Bioconductor submission

1. Ensure package is on Bioconductor devel (R >= 4.5).
2. `git clone git@git.bioconductor.org:packages/timeProtState`
3. Push to `devel` branch.
4. Track build status on Bioconductor single-package builder.
5. Respond to reviewer comments in the Bioc-devel tracker.

## Reproducibility

The companion repository `timeProtState-paper` pins the exact versions used
for manuscript figures via `renv.lock` and provides `make reproduce` or
`targets::tar_make()` entry points.
