options(testthat.use_colours = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a
suppressMessages(library(testthat))

workspace_root <- normalizePath(
  Sys.getenv("TPS_WORKSPACE", unset = getwd()), mustWork = TRUE
)
pkg <- file.path(workspace_root, "04_package", "timeProtState")
if (!file.exists(file.path(pkg, "DESCRIPTION")) &&
    file.exists(file.path(workspace_root, "DESCRIPTION"))) {
  pkg <- workspace_root
}
if (!file.exists(file.path(pkg, "DESCRIPTION"))) {
  stop("Run from the workspace/package root or set TPS_WORKSPACE.")
}
pkgload::load_all(pkg, quiet = TRUE)
data("tps_example", package = "timeProtState", envir = globalenv())
files <- list.files(file.path(pkg, "tests", "testthat"),
                    pattern = "^test-.*\\.R$", full.names = TRUE)

grand <- list()
for (f in files) {
  r <- tryCatch(test_file(f, package = "timeProtState", reporter = "silent"),
                error = function(e) { cat("ERR", basename(f), conditionMessage(e), "\n"); list() })
  grand <- c(grand, r)
}

all_exp <- list()
for (b in grand) for (er in b$results) all_exp <- c(all_exp, list(er))

fails <- all_exp[sapply(all_exp, inherits, "expectation_failure")]
errs  <- all_exp[sapply(all_exp, inherits, "expectation_error")]
skps  <- all_exp[sapply(all_exp, inherits, "expectation_skip")]
succ  <- all_exp[sapply(all_exp, inherits, "expectation_success")]
warns <- all_exp[sapply(all_exp, inherits, "expectation_warning")]
other <- length(all_exp) - length(succ) - length(fails) - length(errs) - length(skps) - length(warns)

cat("\n=== AUTHORITATIVE STATUS (S3 class) ===\n")
cat(sprintf("test_that blocks : %d\n", length(grand)))
cat(sprintf("expectations     : %d\n", length(all_exp)))
cat(sprintf("SUCCESSES        : %d\n", length(succ)))
cat(sprintf("FAILURES         : %d\n", length(fails)))
cat(sprintf("ERRORS           : %d\n", length(errs)))
cat(sprintf("SKIPS            : %d\n", length(skps)))
cat(sprintf("WARNINGS         : %d\n", length(warns)))
cat(sprintf("OTHER            : %d\n", other))

if (length(fails) > 0) {
  cat("\n=== FAILURE DETAILS ===\n")
  for (i in seq_along(fails)) cat(sprintf("[%d] %s\n", i, fails[[i]]$test))
}
if (length(errs) > 0) {
  cat("\n=== ERROR DETAILS ===\n")
  for (i in seq_along(errs)) cat(sprintf("[%d] %s\n", i, errs[[i]]$test))
}
if (length(skps) > 0) {
  cat("\n=== SKIP DETAILS ===\n")
  for (i in seq_along(skps)) cat(sprintf("[%d] %s\n", i, skps[[i]]$test))
}
