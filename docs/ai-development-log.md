# AI Development Log

> Per Bioinformatics Application Note requirements (DECISIONS D10), this log
> discloses AI-assisted code development. AI is NOT listed as an author.

## Summary

The `timeProtState` R package was developed with assistance from **ZCode + GLM-5.2**.
The AI was used as an implementation, testing, and auditing tool. All statistical
decisions (design, contrasts, thresholds, permutation schemes, state rules) were
made and approved by the human authors. The AI did not autonomously change the
statistical contract.

## Log entries

### 2026-08-05 — Goal 0: Legacy audit
- **model**: GLM-5.2
- **interface**: ZCode
- **branch**: goal-0-audit
- **prompt_summary**: Audited M. smegmatis dormancy proteomics legacy scripts; produced STATISTICAL_CONTRACT.md, LEGACY_PIPELINE_AUDIT.md, REFERENCE_OUTPUT_MANIFEST.tsv, DECISIONS_NEEDED.md.
- **files_changed**: docs/*.md
- **human_review**: approved
- **statistical_changes**: none (read-only audit)

### 2026-08-05 — Goal 1: Package skeleton + data contract
- **model**: GLM-5.2
- **prompt_summary**: Implemented tps_import(), tps_validate(), tps_provenance(), tps_metadata_schema() with shared-baseline detection and ID-mapping audit.
- **files_changed**: R/{import,validate,provenance,utils,data-contract}.R, tests/testthat/test-import-validate.R
- **human_review**: approved
- **statistical_changes**: none

### 2026-08-05 — Goal 2/5/6: Design + contrasts + fit
- **model**: GLM-5.2
- **prompt_summary**: Implemented saturated group-means design, 6 contrast families (full-rank), unified limma fit (trend=TRUE, robust=TRUE), golden-master equivalence tests.
- **files_changed**: R/{design,contrasts,fit}.R, tests/testthat/test-design-contrasts-fit.R, data-raw/build_golden_fixture.R
- **human_review**: approved with changes (reference-invariance reframed as subspace equivalence after empirical finding)
- **statistical_changes**: eBayes unified to trend=TRUE/robust=TRUE (D1); significance default to FDR (D2); interaction to 6-contrast full-rank (D4)

### 2026-08-05 — Goal 3: Evidence-decomposed states
- **model**: GLM-5.2
- **prompt_summary**: Implemented three evidence axes (temporal, heterogeneity, direction) + combined state labels + rule manifest.
- **files_changed**: R/classify-states.R, tests/testthat/test-classify-states.R
- **human_review**: approved
- **statistical_changes**: none (rule structure approved by authors)

### 2026-08-05 — Goal 11/12: Pathway NET/REM/COH + permutation calibration
- **model**: GLM-5.2
- **prompt_summary**: Implemented deterministic NET/REM/COH scores + abundance-stratified permutation null + 6-class pathway state classification. Validated against golden fixture.
- **files_changed**: R/{pathway-scores,pathway-null}.R, tests/testthat/test-pathway-scores.R, data-raw/build_pathway_fixture.R, spike/*
- **human_review**: approved
- **statistical_changes**: abundance_stratified adopted as default null (D5, validated by spike comparison)

### 2026-08-06 — Goal 4/13/14/15: Trajectory, plotting, export, simulation
- **model**: GLM-5.2
- **prompt_summary**: Implemented trajectory matrix (baseline-once fix D6), Mfuzz/kmeans adapter, 5 plot functions, simulation framework, structured export.
- **files_changed**: R/{trajectory-matrix,cluster-trajectories,plotting,simulate,export}.R, tests/testthat/test-*.R
- **human_review**: approved
- **statistical_changes**: trajectory baseline modeled once (D6); reversal detection uses effect-size gate (D8)

### 2026-08-06 — Goal 16/17/18: External data, vignettes, release prep
- **model**: GLM-5.2
- **prompt_summary**: Reproduced PXD006262 (yeast) and verified PXD018814 (Arabidopsis) designs; wrote 2 vignettes; created paper repo with _targets/Makefile/renv/Zenodo config.
- **files_changed**: timeProtState-paper/*, vignettes/*, docs/STATUS.md, RELEASE_PROCESS.md
- **human_review**: approved
- **statistical_changes**: none

## What the AI did NOT do

- Did not choose the statistical methods or thresholds autonomously.
- Did not claim tests passed without running them.
- Did not delete failing tests.
- Did not write the manuscript text as author-original.
- Did not alter the golden-master oracle values.
