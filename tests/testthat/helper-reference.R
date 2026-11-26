## Helpers for comparing RscDRS output with reference output from the
## Python scDRS (v1.0.4, commit 644a5c7). Reference files are made by
## inst/scripts/python/make_reference_outputs.py (seed 0) and stored,
## gzipped, in tests/testthat/fixtures/.

referencePath <- function(...) {
    testthat::test_path("fixtures", ...)
}

## Read reference table `name` (file fixtures/<name>.csv.gz). With
## rowNames = TRUE the first column holds row names. `mean_var` bin labels
## are read as character ("3.10" must not become 3.1), Python's
## "True"/"False" become logical, and empty fields are NA. Skips the test
## if the file is missing.
readReference <- function(name, rowNames = TRUE) {
    path <- referencePath(paste0(name, ".csv.gz"))
    if (!file.exists(path)) {
        testthat::skip(paste("reference file missing:", name))
    }
    header <- names(utils::read.csv(path, nrows = 1, check.names = FALSE))
    cc <- if ("mean_var" %in% header) c(mean_var = "character") else NA
    df <- utils::read.csv(path, row.names = if (rowNames) 1 else NULL,
        check.names = FALSE, colClasses = cc, na.strings = c("NA", ""))
    for (j in seq_along(df)) {
        v <- df[[j]]
        if (is.character(v) && all(v %in% c("True", "False", NA))) {
            df[[j]] <- v == "True"
        }
    }
    df
}

## Fraction of equal elements.
agreement <- function(x, y) {
    mean(x == y)
}

## Control gene sets exported by Python (long format ctrl_id, gene,
## weight), as a list of named numeric vectors.
readControlSets <- function(name) {
    df <- readReference(name, rowNames = FALSE)
    lapply(split(df, df$ctrl_id), function(s) {
        stats::setNames(s$weight, s$gene)
    })
}

## The toy data preprocessed in R, with the gene statistics replaced by
## Python's (fixture pp.<x>.<cov>.none.gene_stats). The covariate
## correction itself matches Python to 1e-6, but `var_tech` can differ for
## genes near the 10 constant genes (see test-scdrsPreprocess.R); using
## Python's statistics makes the scoring tests exact.
toyWithPythonStats <- function(cov = TRUE, dense = FALSE) {
    sce <- scdrsLoadToyData()
    if (dense) {
        assay(sce, "logcounts") <- as.matrix(assay(sce, "logcounts"))
    }
    sce <- suppressMessages(scdrsPreprocess(sce,
        cov = if (cov) "covariate"))
    tag <- paste("pp", if (dense) "dense" else "sparse",
        if (cov) "cov" else "nocov", "none", sep = ".")
    ref <- readReference(paste0(tag, ".gene_stats"))
    rd <- rowData(sce)
    for (col in names(ref)) {
        rd[[paste0("scdrs_", col)]] <- ref[rownames(sce), col]
    }
    rowData(sce) <- rd
    sce
}

## The toy mouse gene set (unweighted).
toyGeneSet <- function() {
    scdrsReadGeneSets(system.file("extdata", "toydata_mouse.gs",
        package = "RscDRS"))[["toydata_gs_mouse"]]
}

## A Python score table (fixture <prefix>.score) as a matrix-free
## data.frame, ordered like `cells`.
readScores <- function(prefix, cells) {
    readReference(paste0(prefix, ".score"))[cells, ]
}
