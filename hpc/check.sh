#!/usr/bin/env bash
# Document, test, build and check scdRs inside the Bioconductor devel
# image. Run from the repository root:  bash hpc/check.sh [test|check]
#   test  : roxygen + unit tests only (fast)
#   check : full R CMD build / check + BiocCheck (default)
# Logs go to hpc/logs/; paste the summary back to Claude.
set -uo pipefail

MODE="${1:-check}"
IMG="${IMG:-$HOME/containers/bioc-devel.sif}"
RLIB="${SCDRS_RLIB:-$HOME/R/bioc-devel}"
CT="$(command -v apptainer || command -v singularity)"
LOG_DIR="hpc/logs"
mkdir -p "$LOG_DIR"
STAMP="$(date +%Y%m%d-%H%M%S)"
LOG="$LOG_DIR/${MODE}-${STAMP}.log"

run() {
    "$CT" exec --env SCDRS_RLIB="$RLIB" --env R_LIBS="$RLIB" "$IMG" "$@"
}

{
    echo "== git: $(git rev-parse --short HEAD) =="
    run Rscript -e '
        .libPaths(c(Sys.getenv("SCDRS_RLIB"), .libPaths()))
        roxygen2::roxygenise()
        res <- as.data.frame(devtools::test(stop_on_failure = FALSE))
        cat("\n==== tests: failed =", sum(res$failed),
            " skipped =", sum(res$skipped),
            " passed =", sum(res$passed), "====\n")
    '
    if [ "$MODE" = "check" ]; then
        TMP="$(mktemp -d)"
        run R CMD build --no-build-vignettes . && mv scdRs_*.tar.gz "$TMP"/
        TAR="$(ls "$TMP"/scdRs_*.tar.gz)"
        (cd "$TMP" && run R CMD check --no-manual "$TAR")
        echo "==== R CMD check summary ===="
        tail -n 20 "$TMP"/scdRs.Rcheck/00check.log
        run Rscript -e '
            .libPaths(c(Sys.getenv("SCDRS_RLIB"), .libPaths()))
            BiocCheck::BiocCheck(commandArgs(TRUE)[1])
        ' "$TAR"
    fi
} 2>&1 | tee "$LOG"

echo "Log written to $LOG"
