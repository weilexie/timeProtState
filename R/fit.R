## ============================================================================
## Fit: limma moderated linear model + contrast testing
## ============================================================================
##
## Implements STATISTICAL_CONTRACT.md section 2.4-2.6 and DECISIONS D1.
## The package default is trend = TRUE, robust = TRUE, applied UNIFORMLY to
## every limma fit (resolving the legacy 3-way inconsistency). The significance
## default is adj.P.Val < 0.05 with lfc_threshold = 0; the legacy_msm profile
## (raw P + 0.74) is handled in tps_decision_rule(), not here.

## ---- tps_fit ---------------------------------------------------------------
#' Fit a limma model to a time-course proteomics experiment
#'
#' Fits a moderated linear model with \code{\link[limma]{limma}} and tests all
#' contrast families from a \code{tps_contrasts} object (produced by
#' \code{\link{tps_make_contrasts}}). Returns a \code{tps_fit} carrying
#' per-contrast moderated t/F statistics, logFC, and P-values, plus omnibus
#' F-tests for the average-time, heterogeneity, and within-condition families.
#'
#' The default \code{trend = TRUE, robust = TRUE} is applied uniformly across
#' all families so every contrast family uses the same eBayes settings.
#'
#' @param se a \code{SummarizedExperiment} from \code{\link{tps_import}}.
#' @param design a \code{tps_design} from \code{\link{tps_make_design}}.
#' @param contrasts a \code{tps_contrasts} from \code{\link{tps_make_contrasts}}.
#' @param assay_name character; assay to fit. If \code{NULL}, uses the assay
#'   recorded at import.
#' @param trend logical; passed to \code{\link[limma]{eBayes}}. Default TRUE.
#' @param robust logical; passed to \code{eBayes}. Default TRUE.
#' @param block optional repeated-measures block. Supply one column name from
#'   \code{colData(se)} (for example, \code{"subject_id"}) or a vector aligned
#'   to assay columns. When supplied without \code{correlation}, the consensus
#'   correlation is estimated with \code{\link[limma]{duplicateCorrelation}}.
#' @param correlation optional within-block correlation supplied to
#'   \code{\link[limma]{lmFit}}. Requires \code{block}; when \code{NULL}, the
#'   correlation is estimated from the data.
#' @return a list of class \code{tps_fit} with:
#'   \item{statistics}{per-contrast (and per-family) result table.}
#'   \item{fit}{the limma MArrayLM object.}
#'   \item{thresholds}{the trend/robust settings plus optional
#'     \code{block}, \code{correlation}, and \code{correlation_source}.}
#'   \item{provenance}{provenance entry.}
#'
#' @examples
#' data(tps_example, package = "timeProtState")
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' design <- tps_make_design(se)
#' contrasts <- tps_make_contrasts(design, endpoint = 7)
#' fit <- tps_fit(se, design, contrasts)
#' head(fit$statistics[fit$statistics$family == "endpoint", ])
#'
#' @export
tps_fit <- function(se, design, contrasts, assay_name = NULL,
                     trend = TRUE, robust = TRUE, block = NULL,
                     correlation = NULL) {
  block_capture <- if (is.null(block)) {
    NULL
  } else if (is.character(block) && length(block) == 1L) {
    block
  } else {
    "<aligned vector>"
  }
  arg_capture <- list(call = "tps_fit", trend = trend, robust = robust,
                       block = block_capture, correlation = correlation,
                       timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  if (!inherits(se, "SummarizedExperiment")) {
    tps_abort("se must be a SummarizedExperiment.")
  }
  if (!inherits(design, "tps_design")) {
    tps_abort("design must be a tps_design object.")
  }
  if (!inherits(contrasts, "tps_contrasts")) {
    tps_abort("contrasts must be a tps_contrasts object.")
  }
  if (!is.null(design$rank_audit) && !isTRUE(design$rank_audit$full_rank)) {
    tps_abort("design matrix is rank-deficient; resolve confounded cells/batches before fitting.")
  }
  if (length(trend) != 1L || is.na(trend) || !is.logical(trend) ||
      length(robust) != 1L || is.na(robust) || !is.logical(robust)) {
    tps_abort("trend and robust must each be one non-missing logical value.")
  }
  if (is.null(assay_name)) {
    assay_name <- S4Vectors::metadata(se)$timeProtState$assay_name
    if (is.null(assay_name)) assay_name <- "log2_abundance"
  }
  if (!is.character(assay_name) || length(assay_name) != 1L ||
      is.na(assay_name) ||
      !assay_name %in% SummarizedExperiment::assayNames(se)) {
    tps_abort("assay_name must identify one assay present in se.")
  }
  if (!identical(contrasts$design_levels, colnames(design$design_matrix))) {
    tps_abort("contrasts were not built from the supplied design matrix.")
  }
  expr <- SummarizedExperiment::assay(se, assay_name)

  ## Align expr columns to design rows
  if (ncol(expr) != nrow(design$design_matrix)) {
    tps_abort("expr columns (%d) != design rows (%d).",
              ncol(expr), nrow(design$design_matrix))
  }
  if (!identical(colnames(expr), rownames(design$design_matrix))) {
    tps_abort("assay columns and design rows are not in identical sample order.")
  }

  ## Optional repeated-measures correlation. This is opt-in because a column
  ## called "replicate" can denote either independent replicate labels or the
  ## same biological unit followed over time; the package must not guess.
  block_value <- NULL
  block_name <- NULL
  correlation_used <- NULL
  correlation_source <- "none"
  if (!is.null(block)) {
    if (is.character(block) && length(block) == 1L && !is.na(block) &&
        nzchar(block)) {
      cd <- as.data.frame(SummarizedExperiment::colData(se))
      if (!block %in% names(cd)) {
        tps_abort("block column '%s' is not present in colData(se).", block)
      }
      block_value <- cd[[block]]
      block_name <- block
    } else if (length(block) == ncol(expr)) {
      block_value <- block
      block_name <- "<aligned vector>"
    } else {
      tps_abort(paste0(
        "block must be NULL, one colData column name, or a vector of length ",
        ncol(expr), "."
      ))
    }
    if (anyNA(block_value)) {
      tps_abort("block contains missing values.")
    }
    block_value <- factor(block_value)
    if (nlevels(block_value) < 2L || !anyDuplicated(block_value)) {
      tps_abort("block must contain at least two repeated biological units.")
    }

    if (is.null(correlation)) {
      corfit <- limma::duplicateCorrelation(
        expr, design$design_matrix, block = block_value
      )
      correlation_used <- unname(corfit$consensus.correlation)
      correlation_source <- "estimated"
      if (length(correlation_used) != 1L || !is.finite(correlation_used)) {
        tps_abort("duplicateCorrelation did not return one finite correlation.")
      }
    } else {
      if (!is.numeric(correlation) || length(correlation) != 1L ||
          is.na(correlation) || !is.finite(correlation) ||
          correlation <= -1 || correlation >= 1) {
        tps_abort("correlation must be one finite numeric value in (-1,1).")
      }
      correlation_used <- unname(correlation)
      correlation_source <- "supplied"
    }
  } else if (!is.null(correlation)) {
    tps_abort("correlation requires a non-NULL block.")
  }

  ## Fit with limma. lmFit handles NAs natively (DECISIONS D3: no imputation).
  fit0 <- limma::lmFit(
    expr, design$design_matrix,
    block = block_value, correlation = correlation_used
  )
  fitted0 <- fit0$coefficients %*% t(design$design_matrix)
  residuals0 <- expr - fitted0
  dimnames(residuals0) <- dimnames(expr)

  ## For each family, contrasts.fit + eBayes (uniform params), then extract.
  stats_list <- list()
  fit_objs <- list()
  for (fn in names(contrasts$families)) {
    cm <- contrasts$families[[fn]]
    if (is.null(cm) || length(cm) == 0L) next
    if (nrow(cm) != ncol(design$design_matrix)) {
      tps_abort("contrast family '%s' has %d rows but design has %d columns.",
                fn, nrow(cm), ncol(design$design_matrix))
    }
    test_kind <- unique(contrasts$family_meta[[fn]]$test)
    fc <- limma::contrasts.fit(fit0, cm)
    fc <- limma::eBayes(fc, trend = trend, robust = robust)
    fit_objs[[fn]] <- fc
    stats_list[[fn]] <- .extract_family(fc, fn, test_kind, contrasts)
  }
  statistics <- do.call(rbind, stats_list)
  rownames(statistics) <- NULL

  out <- structure(
    list(
      statistics = statistics,
      fit = fit_objs,
      thresholds = list(
        trend = trend, robust = robust, assay_name = assay_name,
        block = block_name, correlation = correlation_used,
        correlation_source = correlation_source
      ),
      contrasts = contrasts,
      design = design,
      residuals = residuals0,
      annotation = as.data.frame(SummarizedExperiment::rowData(se)),
      provenance = list(.make_provenance(arg_capture, se, stage = "fit"))
    ),
    class = "tps_fit"
  )
  out
}

## ---- extract per family ---------------------------------------------------
.extract_family <- function(fc, family_name, test_kind, contrasts) {
  if (test_kind == "t") {
    ## Single-contrast t-tests: topTable per coef.
    coefs <- colnames(fc$coefficients)
    rows <- list()
    for (cc in coefs) {
      tt <- limma::topTable(fc, coef = cc, number = Inf, sort.by = "none")
      rows[[cc]] <- data.frame(
        feature_id = rownames(tt),
        family = family_name,
        contrast = cc,
        test = "t",
        logFC = tt$logFC,
        AveExpr = tt$AveExpr,
        t = tt$t,
        F = NA_real_,
        P.Value = tt$P.Value,
        adj.P.Val = tt$adj.P.Val,
        B = tt$B,
        stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
  } else if (test_kind == "F_omnibus") {
    meta <- contrasts$family_meta[[family_name]]
    blocks <- if ("condition" %in% names(meta)) {
      split(seq_len(ncol(fc$coefficients)), meta$condition)
    } else {
      list(omnibus = seq_len(ncol(fc$coefficients)))
    }
    rows <- lapply(names(blocks), function(block) {
      idx <- blocks[[block]]
      tt <- limma::topTable(fc, coef = idx, number = Inf, sort.by = "none")
      label <- if (identical(block, "omnibus")) {
        sprintf("%s_omnibus", family_name)
      } else {
        sprintf("%s_%s", family_name, block)
      }
      data.frame(
        feature_id = rownames(tt),
        family = family_name,
        contrast = label,
        test = "F_omnibus",
        logFC = NA_real_,
        AveExpr = tt$AveExpr,
        t = NA_real_,
        F = tt$F,
        P.Value = tt$P.Value,
        adj.P.Val = tt$adj.P.Val,
        B = NA_real_,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  } else {
    tps_abort("unknown test_kind '%s'; expected 't' or 'F_omnibus'.", test_kind)
  }
}

## ---- tps_results: accessor ------------------------------------------------
#' Extract a tidy results table from a fit
#'
#' Returns the per-feature, per-contrast statistics table produced by
#' \code{\link{tps_fit}}, optionally filtered to a family or contrast.
#'
#' @param fit a \code{tps_fit} object.
#' @param family optional character; restrict to one contrast family.
#' @param contrast optional character; restrict to one contrast name.
#' @return a data.frame.
#'
#' @examples
#' data(tps_example, package = "timeProtState")
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' design <- tps_make_design(se)
#' contrasts <- tps_make_contrasts(design, endpoint = 7)
#' fit <- tps_fit(se, design, contrasts)
#' head(tps_results(fit, family = "endpoint"))
#'
#' @export
tps_results <- function(fit, family = NULL, contrast = NULL) {
  if (!inherits(fit, "tps_fit")) tps_abort("fit must be a tps_fit object.")
  for (arg in list(family = family, contrast = contrast)) {
    value <- arg[[1L]]
    if (!is.null(value) &&
        (!is.character(value) || length(value) != 1L ||
         is.na(value) || !nzchar(value))) {
      tps_abort("family and contrast filters must be NULL or one non-empty name.")
    }
  }
  d <- fit$statistics
  if (!is.null(family)) d <- d[d$family == family, , drop = FALSE]
  if (!is.null(contrast)) d <- d[d$contrast == contrast, , drop = FALSE]
  d
}

## ---- statistic matrix accessor (for pathway scoring) ----------------------
#' Build a feature x contrast statistic matrix
#'
#' Pivots a chosen statistic (moderated t by default) from a \code{tps_fit}
#' into a matrix suitable for pathway-level aggregation (NET/REM/COH).
#'
#' @param fit a \code{tps_fit} object.
#' @param type character; one of \code{"moderated_t"}, \code{"logFC"},
#'   \code{"signed_log10p"}.
#' @param family optional; restrict to a family.
#' @return a numeric matrix (features x contrasts).
#'
#' @examples
#' data(tps_example, package = "timeProtState")
#' se <- tps_import(tps_example_abundance, tps_example_metadata)
#' design <- tps_make_design(se)
#' contrasts <- tps_make_contrasts(design, endpoint = 7)
#' fit <- tps_fit(se, design, contrasts)
#' tmat <- tps_statistic_matrix(fit, family = "endpoint")
#' dim(tmat)
#'
#' @export
tps_statistic_matrix <- function(fit, type = c("moderated_t","logFC","signed_log10p"),
                                  family = NULL) {
  type <- match.arg(type)
  d <- tps_results(fit, family = family)
  d <- d[d$test == "t", , drop = FALSE]  # only single contrasts have direction
  if (nrow(d) == 0L) tps_abort("No single-contrast (t) statistics to pivot.")
  val_col <- switch(type,
                    moderated_t = "t",
                    logFC = "logFC",
                    signed_log10p = {
                      s <- sign(d$logFC) * -log10(pmax(d$P.Value, .Machine$double.xmin))
                      over <- abs(s) > 30 & !is.na(s)
                      s[over] <- sign(s[over]) * 30
                      s
                    })
  if (type != "signed_log10p") val_col <- d[[val_col]]
  ## tapply drops dimnames via unname; build a proper matrix preserving them.
  feat <- unique(d$feature_id)
  ctr  <- unique(d$contrast)
  m <- matrix(NA_real_, nrow = length(feat), ncol = length(ctr),
              dimnames = list(feat, ctr))
  m[cbind(match(d$feature_id, feat), match(d$contrast, ctr))] <- val_col
  m
}

## ---- methods ---------------------------------------------------------------
#' @export
print.tps_fit <- function(x, ...) {
  cat("<timeProtState fit>\n")
  cat(sprintf("  trend: %s   robust: %s\n",
              x$thresholds$trend, x$thresholds$robust))
  if (!is.null(x$thresholds$block)) {
    cat(sprintf("  block: %s   correlation: %.3f (%s)\n",
                x$thresholds$block, x$thresholds$correlation,
                x$thresholds$correlation_source))
  }
  fams <- unique(x$statistics$family)
  cat(sprintf("  families fitted: %s\n", paste(fams, collapse = ", ")))
  cat(sprintf("  statistics rows: %d\n", nrow(x$statistics)))
  invisible(x)
}
