## test-import-validate.R
##
## Tests the data-contract layer (tps_import, tps_validate, tps_provenance,
## tps_metadata_schema). Every error path documented in STATISTICAL_CONTRACT.md
## section 4.2 and LEGACY_PIPELINE_AUDIT.md section 3 must have a test.

library(testthat)
library(timeProtState)
library(SummarizedExperiment)

## ---- shared fixtures ------------------------------------------------------
good_se <- function() {
  tps_import(tps_example_abundance, tps_example_metadata, tps_example_annotation)
}

## ===========================================================================
## tps_metadata_schema
## ===========================================================================
test_that("tps_metadata_schema returns required columns", {
  s <- tps_metadata_schema()
  expect_s3_class(s, "tps_metadata_schema")
  expect_equal(s$required, c("sample_id","condition","time","is_baseline"))
  expect_true(all(c("replicate","batch") %in% s$optional))
  expect_output(print(s), "required:")
})

## ===========================================================================
## tps_import — happy path
## ===========================================================================
test_that("tps_import builds a valid SummarizedExperiment", {
  se <- good_se()
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(dim(se), c(200L, 30L))
  expect_equal(rownames(se), rownames(tps_example_abundance))
  expect_equal(colnames(se), tps_example_metadata$sample_id)
})

test_that("tps_import records provenance", {
  se <- good_se()
  pv <- tps_provenance(se)
  expect_type(pv, "list")
  expect_equal(pv[[1]]$stage, "import")
  expect_match(pv[[1]]$r_version, "^R version")
  expect_true(nzchar(pv[[1]]$input_checksum))
})

test_that("tps_import coerces data.frame abundance to matrix", {
  df <- as.data.frame(tps_example_abundance)
  se <- tps_import(df, tps_example_metadata)
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(dim(se)[2], 30L)
})

test_that("tps_import accepts character is_baseline", {
  md <- tps_example_metadata
  md$is_baseline <- ifelse(md$is_baseline, "true", "false")
  se <- tps_import(tps_example_abundance, md)
  expect_true(all(colData(se)$is_baseline[1:3]))
})

test_that("tps_import aligns col_data to abundance colnames regardless of order", {
  md <- tps_example_metadata[sample(seq_len(30L)), ]
  se <- tps_import(tps_example_abundance, md)
  expect_equal(colnames(se), colnames(tps_example_abundance))
})

## ===========================================================================
## tps_import — error paths (must give specific messages)
## ===========================================================================
test_that("tps_import rejects non-matrix non-data.frame abundance", {
  expect_error(tps_import(1:10, tps_example_metadata),
               "must be a matrix or data.frame")
})

test_that("tps_import rejects abundance without row names", {
  m <- tps_example_abundance; rownames(m) <- NULL
  expect_error(tps_import(m, tps_example_metadata),
               "must have unique row names")
})

test_that("tps_import rejects duplicated protein ids", {
  m <- tps_example_abundance
  rownames(m)[2] <- rownames(m)[1]
  expect_error(tps_import(m, tps_example_metadata), "duplicated protein ids")
})

test_that("tps_import rejects col_data missing required columns", {
  md <- tps_example_metadata; md$condition <- NULL
  expect_error(tps_import(tps_example_abundance, md), "missing required column")
})

test_that("tps_import rejects abundance columns without metadata", {
  m <- cbind(tps_example_abundance, ghost = 1)
  expect_error(tps_import(m, tps_example_metadata), "no metadata row")
})

test_that("tps_import rejects metadata rows without abundance columns", {
  md <- rbind(tps_example_metadata,
              data.frame(sample_id="ghost", condition="Acid", time=1,
                         replicate="R1", is_baseline=FALSE, row.names="ghost"))
  expect_error(tps_import(tps_example_abundance, md), "no abundance column")
})

test_that("tps_import rejects row_data misaligned to abundance", {
  rd <- tps_example_annotation[1:10, ]
  expect_error(tps_import(tps_example_abundance, tps_example_metadata, rd),
               "no row_data annotation")
})

test_that("tps_import rejects non-logical is_baseline values", {
  md <- tps_example_metadata; md$is_baseline <- c(rep("yes", 3), rep("maybe", 27))
  expect_error(tps_import(tps_example_abundance, md), "non-logical values")
})

## ===========================================================================
## tps_validate — design detection
## ===========================================================================
test_that("tps_validate detects shared-baseline design", {
  se <- good_se()
  rep <- tps_validate(se)
  expect_true(rep$ok)
  expect_true(rep$design$shared_baseline)
  expect_gt(length(rep$design$non_baseline_times), 1L)
})

test_that("tps_validate warns when baseline spans multiple times", {
  md <- tps_example_metadata
  ## make some baseline rows have time 1 instead of 0
  md$time[1] <- 1
  ## need >=1 row still baseline at original time to keep a baseline
  se <- tps_import(tps_example_abundance, md, validate = FALSE)
  rep <- tps_validate(se)
  expect_false(rep$design$shared_baseline)
  expect_true(any(grepl("NOT a.*shared-baseline", rep$warnings)))
})

test_that("tps_validate warns on all-NA proteins", {
  m <- tps_example_abundance
  m[1, ] <- NA
  se <- tps_import(m, tps_example_metadata, validate = FALSE)
  rep <- tps_validate(se)
  expect_true(any(grepl("all-NA", rep$warnings)))
  expect_gt(rep$design$n_all_na_proteins, 0L)
})

test_that("tps_validate errors when fewer than 2 non-baseline timepoints", {
  md <- tps_example_metadata
  md <- md[md$time %in% c(0, 1), ]          # only day0 + day1
  m <- tps_example_abundance[, md$sample_id]
  expect_error(tps_validate(tps_import(m, md, validate = FALSE), action = "abort"),
               "at least 2 non-baseline")
})

test_that("tps_validate warns on low-replicate groups", {
  md <- tps_example_metadata
  ## drop 2 of 3 Acid day1 replicates
  md <- md[!(md$condition == "Acid" & md$time == 1 & md$replicate %in% c("R2","R3")), ]
  m <- tps_example_abundance[, md$sample_id]
  se <- tps_import(m, md)
  rep <- tps_validate(se)
  expect_true(any(grepl("<2 samples", rep$warnings)))
})

test_that("tps_validate action='abort' stops on duplicate sample_id", {
  se <- good_se()
  ## post-hoc corruption: duplicate a sample_id value in colData
  cd <- SummarizedExperiment::colData(se)
  cd$sample_id[2] <- cd$sample_id[1]
  SummarizedExperiment::colData(se) <- cd
  expect_error(tps_validate(se, action = "abort"), "duplicated")
})

test_that("tps_validate action='report' never stops on a corrupted object", {
  se <- good_se()
  ## post-hoc corruption: blank out the condition column to simulate data loss
  SummarizedExperiment::colData(se)$condition <- NULL
  rep <- tps_validate(se, action = "report")
  expect_false(rep$ok)
  expect_gt(length(rep$errors), 0L)
})

test_that("tps_validate print method works", {
  se <- good_se()
  rep <- tps_validate(se)
  expect_output(print(rep), "validation report")
})

## ===========================================================================
## tps_provenance
## ===========================================================================
test_that("tps_provenance returns NULL for foreign objects", {
  expect_null(tps_provenance(list()))
  expect_null(tps_provenance(1:5))
})

test_that("tps_provenance extracts from result list with provenance key", {
  fake <- list(provenance = list(list(stage = "test")))
  pv <- tps_provenance(fake)
  expect_equal(pv[[1]]$stage, "test")
})

## ===========================================================================
## invariants (STATISTICAL_CONTRACT section 2, principles)
## ===========================================================================
test_that("column order of abundance does not affect the imported object", {
  idx <- sample(seq_len(30L))
  m <- tps_example_abundance[, idx]
  se1 <- tps_import(m, tps_example_metadata, validate = TRUE)
  ## colData is re-ordered to follow abundance colnames; sample_id must match
  expect_equal(colnames(se1), colnames(m))
  expect_equal(as.character(colData(se1)$sample_id), colnames(m))
  ## and the assay content must equal the (reordered) input
  expect_equal(SummarizedExperiment::assay(se1), m)
})

test_that("is_baseline truth preserved through import", {
  se <- good_se()
  expect_equal(sum(colData(se)$is_baseline), 3L)
  expect_true(all(colData(se)$is_baseline[colData(se)$time == 0]))
})
