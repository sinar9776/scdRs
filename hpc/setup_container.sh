#!/usr/bin/env bash
# Pull the Bioconductor devel image (R 4.6.x) and install the packages
# needed to develop RscDRS into a personal library.
# Usage: bash hpc/setup_container.sh
set -euo pipefail

IMG_DIR="${IMG_DIR:-$HOME/containers}"
IMG="$IMG_DIR/bioc-devel.sif"
RLIB="${SCDRS_RLIB:-$HOME/R/bioc-devel}"

CT="$(command -v apptainer || command -v singularity || true)"
if [ -z "$CT" ]; then
    echo "apptainer/singularity not found. Try: module load apptainer"
    exit 1
fi

mkdir -p "$IMG_DIR" "$RLIB"
"$CT" pull --force "$IMG" docker://bioconductor/bioconductor_docker:devel

"$CT" exec --env SCDRS_RLIB="$RLIB" "$IMG" Rscript -e '
    .libPaths(c(Sys.getenv("SCDRS_RLIB"), .libPaths()))
    pkgs <- c("SingleCellExperiment", "BiocCheck", "biocthis",
        "testthat", "devtools", "roxygen2", "BiocStyle", "knitr",
        "rmarkdown", "ggplot2", "irlba", "BiocNeighbors", "remotes")
    BiocManager::install(pkgs, lib = Sys.getenv("SCDRS_RLIB"),
        update = FALSE, ask = FALSE)
    cat("\n==== summary ====\n")
    print(R.version.string)
    print(BiocManager::version())
    print(.libPaths())
    print(sapply(pkgs, requireNamespace, quietly = TRUE))
'
echo "Done. Image: $IMG  Library: $RLIB"
