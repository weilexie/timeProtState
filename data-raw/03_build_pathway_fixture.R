## data-raw/build_pathway_fixture.R
##
## Builds a small golden-master fixture for the pathway module, from the real
## M. smegmatis dormancy data. Freezes:
##   - the 3-contrast moderated-t statistic matrix (subset of 50 proteins)
##   - the AveExpr abundance vector (for stratification)
##   - the locus-tag annotation (accession -> MSMEG_)
##   - the KEGG term2gene subset restricted to those 50 proteins' pathways
##   - the expected NET/REM/COH computed directly (oracle)
##   - the expected p_net under unstratified permutation, seed=1, n_perm=500
##     (deterministic given seed)
##
## This lets the package tests assert:
##   (a) tps_pathway_scores reproduces the oracle NET/REM/COH (tol 1e-10)
##   (b) tps_calibrate_pathways reproduces the frozen p_net (tol 0, exact)

suppressPackageStartupMessages({ library(limma); library(Matrix) })

PROJ <- "C:/Users/11377/Desktop/Msm_dormancy/Msm_dormancy_proteomics"
LIMMA_ALL_DIR <- file.path(PROJ, "02_results/04_limma/01_tables_all")
KEGG_T2G <- file.path(PROJ, "01_data/03_mapping_tables/KEGG_mapping/kegg_pathway2gene.tsv")
CONTRASTS <- c(endpoint_Acid_t7_vs_baseline = "Acid_day7_vs_Control_ALL.csv",
               endpoint_Hypoxia_t7_vs_baseline = "Hypoxia_day7_vs_Control_ALL.csv",
               endpoint_Starvation_t7_vs_baseline = "Starvation_day7_vs_Control_ALL.csv")

## read the 3 limma tables; keep proteins with MSMEG_ locus tags
read_tbl <- function(f) {
  d <- read.csv(file.path(LIMMA_ALL_DIR, f), stringsAsFactors = FALSE, check.names = FALSE)
  locus <- vapply(strsplit(d$PG.Genes, ","), function(tk) {
    tk <- trimws(tk); loc <- grep("^MSMEG_", tk, value = TRUE)
    if (length(loc) >= 1L) loc[1L] else NA_character_
  }, character(1L))
  data.frame(id = d$ProteinAccessions, t = as.numeric(d$t),
             AveExpr = as.numeric(d$AveExpr), locus = locus,
             stringsAsFactors = FALSE)
}
L <- lapply(CONTRASTS, read_tbl)
## common ids with a locus tag
ids <- L[[1]]$id
for (i in seq_along(L)[-1]) ids <- intersect(ids, L[[i]]$id)
has_locus <- L[[1]]$locus[match(ids, L[[1]]$id)]
ids <- ids[!is.na(has_locus)]
## take first 50 (deterministic)
ids <- head(ids, 50L)
locus <- L[[1]]$locus[match(ids, L[[1]]$id)]

## statistic matrix
T_MAT <- do.call(cbind, lapply(L, function(x) x$t[match(ids, x$id)]))
rownames(T_MAT) <- ids; colnames(T_MAT) <- names(CONTRASTS)
AVE <- L[[1]]$AveExpr[match(ids, L[[1]]$id)]; names(AVE) <- ids

## annotation: accession -> locus
ann <- data.frame(feature_id = ids, locus_tag = locus, stringsAsFactors = FALSE)

## KEGG term2gene: keep only pathways that contain >=1 of our 50 locus tags
t2g <- read.delim(KEGG_T2G, header = TRUE, stringsAsFactors = FALSE,
                  colClasses = c("character","character"))
keep_terms <- unique(t2g$term[t2g$gene %in% locus])
t2g_sub <- t2g[t2g$term %in% keep_terms, ]

## ---- oracle NET/REM/COH (direct computation, min_size=2 for the small set) ----
MIN_SIZE <- 2L
EPS <- 1e-12
pathways <- unique(t2g_sub$term)
n_path <- length(pathways)
M <- sparseMatrix(i = integer(0), j = integer(0), dims = c(n_path, length(ids)),
                  dimnames = list(pathways, ids))
sizes <- integer(n_path)
for (p in seq_len(n_path)) {
  mem_locus <- t2g_sub$gene[t2g_sub$term == pathways[p]]
  mem_acc <- ids[locus %in% mem_locus]
  sizes[p] <- length(mem_acc)
  if (length(mem_acc) > 0L) M[p, match(mem_acc, ids)] <- 1
}
keep <- sizes >= MIN_SIZE
M <- M[keep, ]; sizes <- sizes[keep]
Mt <- as(M, "CsparseMatrix")
NET <- as.matrix(Mt %*% T_MAT) / sizes
REM <- as.matrix(Mt %*% abs(T_MAT)) / sizes
COH <- abs(NET) / (REM + EPS)

## ---- oracle p_net (unstratified, seed=1, n_perm=500) ----
set.seed(1L)
n_perm <- 500L
n_prot <- length(ids)
obs_net <- NET
net_null_absmax <- matrix(0.0, nrow(M), ncol(T_MAT))
abs_T <- abs(T_MAT)
for (b in seq_len(n_perm)) {
  idx <- sample.int(n_prot)
  net_b <- as.matrix(Mt %*% T_MAT[idx, , drop = FALSE]) / sizes
  net_null_absmax <- net_null_absmax + (abs(net_b) >= abs(obs_net) - EPS)
}
p_net_oracle <- (net_null_absmax + 1) / (n_perm + 1)

golden_pathway <- list(
  statistic_matrix = T_MAT,
  abundance = AVE,
  annotation = ann,
  term2gene = t2g_sub,
  min_size = MIN_SIZE,
  NET = NET, REM = REM, COH = COH,
  pathway_terms = rownames(M),
  pathway_sizes = sizes,
  p_net_unstratified = p_net_oracle,
  n_perm = n_perm, seed = 1L,
  n_proteins = 50L
)

fixdir <- "tests/testthat/fixtures"
dir.create(fixdir, showWarnings = FALSE, recursive = TRUE)
saveRDS(golden_pathway, file.path(fixdir, "golden_pathway.rds"))
cat("pathway fixture written:", file.path(fixdir, "golden_pathway.rds"), "\n")
cat("  proteins:", golden_pathway$n_proteins,
    " pathways:", length(golden_pathway$pathway_terms),
    " contrasts:", ncol(T_MAT), "\n")
cat("  pathway sizes:", paste(golden_pathway$pathway_sizes, collapse = ","), "\n")
cat("  NET range:", paste(round(range(NET), 3), collapse = ","), "\n")
cat("  COH range:", paste(round(range(COH), 3), collapse = ","), "\n")
cat("  p_net range:", paste(round(range(p_net_oracle), 3), collapse = ","), "\n")
