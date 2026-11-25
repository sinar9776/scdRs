#!/usr/bin/env bash
# Document, test, build and check scdRs inside the Bioconductor devel
# image. Run from the repository root:
#   bash hpc/check.sh [test|vignette|check]
#   test     : roxygen + unit tests only (fast)
#   vignette : roxygen + unit tests + devtools::build_vignettes()
#   check    : unit tests + full R CMD build / check + BiocCheck (default)
# Missing dependencies (Imports and Suggests) are installed first.
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
        d <- read.dcf("DESCRIPTION", c("Depends", "Imports", "Suggests"))
        need <- trimws(sub("[ (].*", "", unlist(strsplit(d, ",[[:space:]]*"))))
        need <- setdiff(need[!is.na(need) & need != ""], c("R",
            rownames(installed.packages(priority = "base"))))
        miss <- need[!vapply(need, requireNamespace, logical(1),
            quietly = TRUE)]
        if (length(miss) > 0L) {
            cat("Installing:", miss, "\n")
            BiocManager::install(miss, lib = Sys.getenv("SCDRS_RLIB"),
                update = FALSE, ask = FALSE)
        }
        roxygen2::roxygenise()
        res <- as.data.frame(devtools::test(stop_on_failure = FALSE))
        cat("\n==== tests: failed =", sum(res$failed),
            " skipped =", sum(res$skipped),
            " passed =", sum(res$passed), "====\n")
    '
    if [ "$MODE" = "vignette" ]; then
        run Rscript -e '
            .libPaths(c(Sys.getenv("SCDRS_RLIB"), .libPaths()))
            devtools::build_vignettes()
            ok <- file.exists("doc/scdRs.html")
            cat("\n==== vignette:", if (ok) "built doc/scdRs.html" else
                "FAILED", "====\n")
        '
    fi
    if [ "$MODE" = "check" ]; then
        TMP="$(mktemp -d)"
        run R CMD build . && mv scdRs_*.tar.gz "$TMP"/
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
