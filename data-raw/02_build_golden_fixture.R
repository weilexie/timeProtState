## data-raw/build_golden_fixture.R
##
## Builds a small golden-master fixture from the real M. smegmatis dormancy
## data, to test that timeProtState reproduces the legacy limma analysis.
##
## Per STATISTICAL_CONTRACT.md section 10 and REFERENCE_OUTPUT_MANIFEST.tsv:
## the fixture must reproduce legacy logFC/t/P/adj.P for a fixed protein subset
## under the legacy_msm profile (imputed input, trend=TRUE, robust=TRUE,
## raw P significance — but the numeric fit is what we freeze).
##
## Output: tests/testthat/fixtures/golden_fit.rds
##         tests/testthat/fixtures/golden_meta.csv (the 20-protein subset metadata)

suppressPackageStartupMessages({ library(limma); library(SummarizedExperiment) })

workspace_root <- normalizePath(
  Sys.getenv("TPS_WORKSPACE", unset = getwd()), mustWork = TRUE
)
PROJ <- normalizePath(Sys.getenv(
  "TPS_DORMANCY_PROJECT",
  unset = file.path(dirname(workspace_root), "Msm_dormancy",
                    "Msm_dormancy_proteomics")
), mustWork = TRUE)
step5 <- read.csv(file.path(PROJ, "01_data/02_processed_data/step5_imputed_Perseus_MNAR.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
meta <- read.csv(file.path(PROJ, "01_data/02_processed_data/Metadata.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)

## abundance matrix: PG.ProteinAccessions as rows, sample columns
## The step5 CSV has columns like "[1] ms251803-C1.raw.PG.Quantity" — strip the
## "[N] " prefix and the ".raw.PG.Quantity" suffix to recover clean sample ids.
feat <- step5$PG.ProteinAccessions
clean_names <- sub("^\\[\\d+\\] ", "", names(step5))
expr_col_idx <- grep("^ms.*\\.raw\\.PG\\.Quantity$", clean_names)
expr <- as.matrix(step5[, expr_col_idx])
colnames(expr) <- sub("\\.raw\\.PG\\.Quantity$", "", clean_names[expr_col_idx])
rownames(expr) <- make.unique(feat)

## align metadata to expr columns
## metadata$sample uses short ids (C1, A1-1...); expr cols are ms*-<short>
## build sample_id matching: strip the msNNNNNN- prefix
col_short <- sub("^ms\\d+-", "", colnames(expr))
meta_map <- match(col_short, meta$sample)
stopifnot(!anyNA(meta_map))
meta_aligned <- meta[meta_map, ]
rownames(meta_aligned) <- colnames(expr)

## take a 20-protein subset (deterministic: first 20 with locus tags)
set.seed(42)
has_locus <- grep("MSMEG_", step5$PG.Genes)
keep_idx <- head(has_locus, 20)
expr_sub <- expr[keep_idx, , drop = FALSE]

## build group factor (legacy: Control + 9 condition_time cells)
group <- with(meta_aligned, paste(ifelse(model == "Control", "baseline",
                                          as.character(model)),
                                   time, sep = "_"))
group <- sub("_0$", "", group)  # Control_0 -> baseline
group_levels <- c("baseline",
                   paste(rep(c("Starvation","Acid","Hypoxia"), each = 3),
                         rep(c(1,3,7), 3), sep = "_"))
group <- factor(group, levels = group_levels)

## legacy design: ~0 + group
design <- model.matrix(~ 0 + group)
colnames(design) <- sub("^group", "", colnames(design))

## 3 vs-baseline contrasts (Acid_day7, Hypoxia_day7, Starvation_day7)
## Names match timeProtState's endpoint family naming convention.
cm <- limma::makeContrasts(
  endpoint_Acid_t7_vs_baseline = Acid_7 - baseline,
  endpoint_Hypoxia_t7_vs_baseline = Hypoxia_7 - baseline,
  endpoint_Starvation_t7_vs_baseline = Starvation_7 - baseline,
  levels = design)

fit0 <- limma::lmFit(expr_sub, design)
fitc <- limma::contrasts.fit(fit0, cm)
fitc <- limma::eBayes(fitc, trend = TRUE, robust = TRUE)

## freeze: per-contrast logFC, t, P.Value, adj.P.Val for the 20 proteins
results <- list()
for (cc in colnames(cm)) {
  tt <- limma::topTable(fitc, coef = cc, number = Inf, sort.by = "none")
  results[[cc]] <- data.frame(
    feature_id = rownames(expr_sub),
    contrast = cc,
    logFC = tt$logFC,
    t = tt$t,
    P.Value = tt$P.Value,
    adj.P.Val = tt$adj.P.Val,
    stringsAsFactors = FALSE)
}
golden <- list(
  statistics = do.call(rbind, results),
  design_matrix = design,
  contrast_matrix = cm,
  group_levels = group_levels,
  conditions = c("Acid","Hypoxia","Starvation"),
  times = c(1, 3, 7),
  endpoint = 7,
  n_proteins = 20L,
  protein_ids = rownames(expr_sub),
  seed = 42L
)

fixdir <- "tests/testthat/fixtures"
dir.create(fixdir, showWarnings = FALSE, recursive = TRUE)
saveRDS(golden, file.path(fixdir, "golden_fit.rds"))

## also save the 20x30 expression + aligned metadata for the equivalence test
golden$expr <- expr_sub
golden$meta <- meta_aligned[, c("sample","model","time")]
saveRDS(golden, file.path(fixdir, "golden_fit.rds"))

cat("golden fixture written:", file.path(fixdir, "golden_fit.rds"), "\n")
cat("  proteins:", golden$n_proteins, " samples:", ncol(expr_sub), "\n")
cat("  contrasts:", paste(colnames(cm), collapse=", "), "\n")
cat("  design rank:", qr(design)$rank, "/", ncol(design), "\n")
