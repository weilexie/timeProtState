## data-raw/01_build_example_data.R
##
## Generates the example datasets shipped with timeProtState. The design is
## SYNTHETIC shared-baseline design: 3 conditions x times 1/3/7 plus shared
## day-0, 3 replicates = 30 samples, ~200 proteins. No real data is included.
##
## This script is NOT run at check time; its outputs are committed under
## data/ (lazy-loaded) and inst/extdata/ (file-based for vignettes).
## Run manually: Rscript data-raw/01_build_example_data.R

suppressPackageStartupMessages({
  library(Matrix)
})

set.seed(2026)

## ---- experimental design --------------------------------------------------
conditions <- c("Control","Acid","Hypoxia","Starvation")
## day0 = Control (shared baseline); Acid/Hypoxia/Starvation at day1/3/7
times <- c(0, 1, 3, 7)
n_rep <- 3L

design_rows <- list()
sid_seq <- 0L
for (cond in conditions) {
  if (cond == "Control") {
    tps_cond <- 0
  } else {
    tps_cond <- times[times != 0]
  }
  for (tt in tps_cond) {
    for (r in seq_len(n_rep)) {
      sid_seq <- sid_seq + 1L
      design_rows[[length(design_rows) + 1L]] <- data.frame(
        sample_id = sprintf("S%02d", sid_seq),
        condition = if (cond == "Control") "baseline" else cond,
        time = tt,
        replicate = sprintf("R%d", r),
        is_baseline = (cond == "Control"),
        stringsAsFactors = FALSE
      )
    }
  }
}
metadata <- do.call(rbind, design_rows)
rownames(metadata) <- metadata$sample_id
stopifnot(nrow(metadata) == 30L)
stopifnot(sum(metadata$is_baseline) == 3L)        # shared baseline = 3 samples

## ---- synthetic proteins ---------------------------------------------------
n_prot <- 200L
proteins <- sprintf("P%04d", seq_len(n_prot))
ann <- data.frame(
  protein_id = proteins,
  locus_tag = sprintf("ORG_%04d", seq_len(n_prot)),
  gene_name = sprintf("gene%04d", seq_len(n_prot)),
  stringsAsFactors = FALSE
)
rownames(ann) <- ann$protein_id

## per-protein baseline abundance on log2 scale, broadly ~ N(18, 4)
mu0 <- rnorm(n_prot, mean = 18, sd = 4)
## per-protein noise sd
sigma <- abs(rnorm(n_prot, mean = 0.5, sd = 0.15)) + 0.1

## ---- build expression -----------------------------------------------------
## Each protein gets an effect profile; some are null, some shared-up,
## some condition-specific, some opposite.
assign_effect <- function(i) {
  e <- i %% 8L
  switch(as.character(e),
    "0" = list(type = "null"),
    "1" = list(type = "shared_up"),
    "2" = list(type = "shared_down"),
    "3" = list(type = "acid_specific_up"),
    "4" = list(type = "hypoxia_specific_up"),
    "5" = list(type = "starvation_specific_down"),
    "6" = list(type = "opposite_ah"),   # Acid up, Hypoxia down
    "7" = list(type = "transient"),
    list(type = "null")
  )
}
effect_trend <- function(t, type) {
  t <- t / 7                              # normalize to [0,1] over the course
  switch(type,
    null = rep(0, length(t)),
    shared_up = 1.5 * t,
    shared_down = -1.5 * t,
    acid_specific_up = 1.2 * t,
    hypoxia_specific_up = 1.2 * t,
    starvation_specific_down = -1.2 * t,
    opposite_ah = 0,                      # handled per-condition below
    transient = 1.0 * exp(-(t - 0.3)^2 / 0.05),
    rep(0, length(t))
  )
}

expr <- matrix(NA_real_, nrow = n_prot, ncol = nrow(metadata),
               dimnames = list(proteins, metadata$sample_id))

for (i in seq_len(n_prot)) {
  eff <- assign_effect(i)$type
  for (j in seq_len(nrow(metadata))) {
    rr <- metadata[j, ]
    if (rr$is_baseline) {
      ## baseline: only noise around mu0
      expr[i, j] <- mu0[i] + rnorm(1, sd = sigma[i])
    } else {
      t <- rr$time
      base_effect <- effect_trend(t, eff)
      ## opposite_ah: Acid up, Hypoxia down
      if (eff == "opposite_ah") {
        if (rr$condition == "Acid")     base_effect <- 1.2 * (t / 7)
        if (rr$condition == "Hypoxia")  base_effect <- -1.2 * (t / 7)
        if (rr$condition == "Starvation") base_effect <- 0
      }
      ## condition-specific only applies to its condition
      if (eff == "acid_specific_up" && rr$condition != "Acid") base_effect <- 0
      if (eff == "hypoxia_specific_up" && rr$condition != "Hypoxia") base_effect <- 0
      if (eff == "starvation_specific_down" && rr$condition != "Starvation") base_effect <- 0
      expr[i, j] <- mu0[i] + base_effect + rnorm(1, sd = sigma[i])
    }
  }
}

## miss_frac is 0 so the shipped example has ZERO missing values.
## Missingness tests inject NA themselves.
miss_frac <- 0
miss_idx <- sample.int(length(expr), size = floor(miss_frac * length(expr)))
expr[miss_idx] <- NA

## ---- minimal pathway mapping ---------------------------------------------
## Twelve deterministic toy pathways. The first eight follow the known
## feature-level effect blocks above, so the shipped example demonstrates
## directional, condition-dependent, transient, and null pathway behavior
## instead of only proving that the API runs. The remaining four overlap
## effect blocks to exercise mixed membership. No real pathway data are used.
effect_type <- vapply(seq_len(n_prot), function(i) assign_effect(i)$type,
                      character(1L))
ids_by_effect <- split(ann$locus_tag, effect_type)
pathway_members <- list(
  toy_shared_up = ids_by_effect$shared_up,
  toy_shared_down = ids_by_effect$shared_down,
  toy_acid_specific_up = ids_by_effect$acid_specific_up,
  toy_hypoxia_specific_up = ids_by_effect$hypoxia_specific_up,
  toy_starvation_specific_down = ids_by_effect$starvation_specific_down,
  toy_opposite_acid_hypoxia = ids_by_effect$opposite_ah,
  toy_transient = ids_by_effect$transient,
  toy_null = ids_by_effect$null,
  toy_mixed_shared = c(head(ids_by_effect$shared_up, 12),
                       head(ids_by_effect$shared_down, 12)),
  toy_mixed_condition = c(head(ids_by_effect$acid_specific_up, 12),
                          head(ids_by_effect$hypoxia_specific_up, 12)),
  toy_condition_response = c(head(ids_by_effect$acid_specific_up, 8),
                             head(ids_by_effect$hypoxia_specific_up, 8),
                             head(ids_by_effect$starvation_specific_down, 8)),
  toy_heterogeneous = unlist(lapply(ids_by_effect[names(ids_by_effect) != "null"],
                                    head, 2), use.names = FALSE)
)
term2gene <- do.call(rbind, lapply(names(pathway_members), function(term) {
  data.frame(term = term, gene = pathway_members[[term]],
             stringsAsFactors = FALSE)
}))
rownames(term2gene) <- NULL

## ---- save: lazy-loaded data/ + file-based inst/extdata -------------------
pkg_root <- normalizePath(".", winslash = "/")
while (!file.exists(file.path(pkg_root, "DESCRIPTION")) &&
       !identical(pkg_root, dirname(pkg_root))) {
  pkg_root <- dirname(pkg_root)
}

## File-based copies under inst/extdata (for vignettes + tests that read files)
extdata <- file.path(pkg_root, "inst", "extdata")
dir.create(extdata, showWarnings = FALSE, recursive = TRUE)
write.csv(expr, file.path(extdata, "example_abundance.csv"), row.names = TRUE)
write.csv(metadata, file.path(extdata, "example_metadata.csv"), row.names = TRUE)
write.csv(ann, file.path(extdata, "example_annotation.csv"), row.names = TRUE)
write.csv(term2gene, file.path(extdata, "example_term2gene.csv"), row.names = FALSE)

## Lazy-loaded copies under data/ (named objects, used in examples)
data_dir <- file.path(pkg_root, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
tps_example_abundance <- expr
tps_example_metadata <- metadata
tps_example_annotation <- ann
tps_example_term2gene <- term2gene
save(tps_example_abundance, tps_example_metadata,
     tps_example_annotation, tps_example_term2gene,
     file = file.path(data_dir, "tps_example.rda"), version = 2)

cat("example data written under data/ and inst/extdata/\n")
cat(sprintf("  abundance: %d proteins x %d samples\n", nrow(expr), ncol(expr)))
cat(sprintf("  pathways: %d\n", length(unique(term2gene$term))))
