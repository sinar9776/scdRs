## scdrsPreprocess() and scdrsComputeGeneStats() compared with Python
## scdrs.preprocess() on the toy data, in all 8 modes
## (sparse/dense assay x covariates or not x adj_prop or not). Reference:
## fixtures/pp.<x>.<cov>.<adj>.*.csv.gz. No random numbers involved.
##
## Tolerances: Python works in float32 when there are no covariates, so
## those statistics are compared at 1e-5 (per-cell variances, summed over
## 2,500 genes in float32 with E[x^2] - E[x]^2, at 1e-4); with covariates
## Python uses double precision and the tolerance is 1e-6. Bin labels can
## change for genes lying on a quantile cut point, so they must agree for
## >= 99.5% of genes.
##
## Technical variance with covariates: 10 toy genes
## (fixtures/pp.constant_genes) are an exact linear function of the
## covariate, so after correction their variance is floating-point noise
## (~1e-28), different in R and Python, which still enters the LOESS fit
## (as in Python). This moves the fit for genes with similar means (up to
## ~20% for ~1% of genes). Whether that noise is tiny or exactly 0 depends
## on the BLAS: with an exact 0 (e.g. OpenBLAS in the Bioconductor devel
## image) the 10 genes are left out of the fit, as Python would do too,
## and the mean relative difference to Python's file rises from ~1e-4 to
## ~4e-3. Therefore:
## - with covariates on the full toy data, var_tech and ct_var_tech are
##   compared loosely (correlation >= 0.9999, mean relative difference
##   <= 1e-2); all other statistics strictly;
## - the same scenarios without those 10 genes (Python run
##   "nonconst") are compared strictly, var_tech included;
## - the LOESS itself is tested exactly on Python's inputs (below).

numCols <- c("mean", "var", "var_tech", "ct_mean", "ct_var", "ct_var_tech")

checkScenario <- function(xMode, covMode, adj, nonconst = FALSE) {
    tag <- paste("pp", xMode, covMode, adj, sep = ".")
    sce <- scdrsLoadToyData()
    if (nonconst) {
        tag <- paste0(tag, ".nonconst")
        const <- readReference("pp.constant_genes", rowNames = FALSE)$gene
        sce <- sce[!rownames(sce) %in% const, ]
    }
    loose <- covMode == "cov" && !nonconst
    if (xMode == "dense") {
        assay(sce, "logcounts") <- as.matrix(assay(sce, "logcounts"))
    }
    res <- suppressMessages(scdrsPreprocess(sce,
        cov = if (covMode == "cov") "covariate",
        adjProp = if (adj == "celltype") "cell_type"))
    tol <- if (covMode == "cov") 1e-6 else 1e-5

    ref <- readReference(paste0(tag, ".gene_stats"))
    rd <- rowData(res)[rownames(ref), ]
    tech <- c("var_tech", "ct_var_tech")
    for (col in if (loose) setdiff(numCols, tech) else numCols) {
        expect_equal(rd[[paste0("scdrs_", col)]], ref[[col]],
            tolerance = tol, label = paste(tag, col))
    }
    if (loose) {
        keep <- !rownames(ref) %in%
            readReference("pp.constant_genes", rowNames = FALSE)$gene
        for (col in tech) {
            r <- rd[[paste0("scdrs_", col)]][keep]
            expect_gte(cor(r, ref[[col]][keep]), 0.9999)
            expect_equal(r, ref[[col]][keep], tolerance = 1e-2,
                label = paste(tag, col))
        }
    }
    expect_gte(agreement(rd$scdrs_mean_var, ref$mean_var), 0.995)

    refCell <- readReference(paste0(tag, ".cell_stats"))
    cd <- colData(res)[rownames(refCell), ]
    expect_equal(cd$scdrs_mean, refCell$mean, tolerance = tol)
    expect_equal(cd$scdrs_var, refCell$var,
        tolerance = if (covMode == "cov") tol else 1e-4)

    flags <- readReference(paste0(tag, ".flags"), rowNames = FALSE)
    md <- metadata(res)$scdrs
    expect_identical(md$flagSparse, flags$FLAG_SPARSE)
    expect_identical(md$flagCov, flags$FLAG_COV)
    expect_identical(md$flagAdjProp, flags$FLAG_ADJ_PROP)
    expect_identical(md$implicit, xMode == "sparse" && covMode == "cov")

    if (xMode == "sparse" && covMode == "cov") {
        beta <- readReference(paste0(tag, ".cov_beta"))
        expect_equal(rowData(res)[rownames(beta), "scdrs_cov_beta"][,
            colnames(beta)], as.matrix(beta), tolerance = 1e-6,
            ignore_attr = TRUE)
        mu <- readReference(paste0(tag, ".cov_gene_mean"))
        expect_equal(rowData(res)[rownames(mu), "scdrs_cov_gene_mean"],
            mu$cov_gene_mean, tolerance = 1e-6)
        cmat <- readReference(paste0(tag, ".cov_mat"))
        expect_equal(reducedDim(res, "scdrs_cov")[rownames(cmat),
            colnames(cmat)], as.matrix(cmat), ignore_attr = TRUE)
    }
    invisible(res)
}

for (xMode in c("sparse", "dense")) {
    for (covMode in c("cov", "nocov")) {
        for (adj in c("none", "celltype")) {
            test_that(paste("scdrsPreprocess matches Python:", xMode,
                covMode, adj), {
                checkScenario(xMode, covMode, adj)
            })
        }
    }
    for (adj in c("none", "celltype")) {
        test_that(paste("scdrsPreprocess matches Python without constant",
            "genes:", xMode, "cov", adj), {
            checkScenario(xMode, "cov", adj, nonconst = TRUE)
        })
    }
}

test_that("dense + cov stores the corrected data like Python's X", {
    sce <- scdrsLoadToyData()
    assay(sce, "logcounts") <- as.matrix(assay(sce, "logcounts"))
    res <- suppressMessages(scdrsPreprocess(sce, cov = "covariate"))
    ref <- as.matrix(readReference("pp.dense.cov.none.corrected_head100"))
    expect_identical(metadata(res)$scdrs$assay, "scdrs_corrected")
    expect_equal(assay(res, "scdrs_corrected")[rownames(ref),
        colnames(ref)], ref, tolerance = 1e-6)
    ## the input assay is never modified
    expect_identical(assay(res, "logcounts"), assay(sce, "logcounts"))
})

test_that("implicit = FALSE on sparse input gives the dense-mode results", {
    sce <- scdrsLoadToyData()
    res <- suppressMessages(scdrsPreprocess(sce, cov = "covariate",
        implicit = FALSE))
    ref <- readReference("pp.dense.cov.none.gene_stats")
    expect_equal(rowData(res)[rownames(ref), "scdrs_var"], ref$var,
        tolerance = 1e-6)
    expect_true("scdrs_corrected" %in% assayNames(res))
})

test_that("cov can be a data.frame like the .cov file", {
    sce <- scdrsLoadToyData()
    covFile <- system.file("extdata", "toydata_covariates.csv",
        package = "scdRs")
    covDf <- read.csv(covFile, row.names = 1)
    a <- suppressMessages(scdrsPreprocess(sce, cov = covDf))
    b <- suppressMessages(scdrsPreprocess(sce, cov = "covariate"))
    expect_identical(rowData(a), rowData(b))
    expect_error(scdrsPreprocess(sce, cov = covDf[1:10, , drop = FALSE]),
        "does not match")
})

test_that("re-running replaces earlier results", {
    sce <- scdrsLoadToyData()
    a <- suppressMessages(scdrsPreprocess(sce, cov = "covariate",
        implicit = FALSE))
    b <- suppressMessages(scdrsPreprocess(a))
    expect_false("scdrs_corrected" %in% assayNames(b))
    expect_false("scdrs_cov" %in% reducedDimNames(b))
    expect_null(rowData(b)$scdrs_cov_beta)
    ref <- readReference("pp.sparse.nocov.none.gene_stats")
    expect_equal(rowData(b)[rownames(ref), "scdrs_mean"], ref$mean,
        tolerance = 1e-5)
})

test_that("input errors are helpful", {
    sce <- scdrsLoadToyData()
    expect_error(scdrsPreprocess(sce, assay = "counts"), "not found")
    expect_error(scdrsPreprocess(sce, cov = "nope"), "not in colData")
    expect_error(scdrsPreprocess(sce, implicit = TRUE), "needs covariates")
    expect_error(scdrsPreprocess(sce, adjProp = "causal_variable"),
        "not categorical")
    expect_error(scdrsPreprocess(as.matrix(assay(sce))),
        "SingleCellExperiment")
})

test_that("scdrsComputeGeneStats matches compute_stats (no covariates)", {
    sce <- scdrsLoadToyData()
    st <- suppressMessages(scdrsComputeGeneStats(sce))
    ref <- readReference("pp.sparse.nocov.none.gene_stats")
    for (col in numCols) {
        expect_equal(st$gene[rownames(ref), col], ref[[col]],
            tolerance = 1e-5, label = col)
    }
    expect_gte(agreement(st$gene[rownames(ref), "mean_var"], ref$mean_var),
        0.995)
    expect_error(scdrsComputeGeneStats(sce, implicit = TRUE),
        "run scdrsPreprocess")
})

test_that("LOESS technical variance matches skmisc.loess exactly", {
    ## Same inputs as Python (its ct_mean / ct_var), so no floating-point
    ## noise from the covariate correction is involved. Without covariates
    ## Python's inputs were float32 and the CSV holds them rounded to ~7
    ## digits, hence 1e-6 there.
    tags <- c(pp.sparse.cov.none = 1e-10, pp.dense.cov.celltype = 1e-10,
        pp.sparse.nocov.none = 1e-6)
    for (tag in names(tags)) {
        ref <- readReference(paste0(tag, ".gene_stats"))
        expect_equal(scdRs:::.techVar(ref$ct_mean, ref$ct_var),
            ref$ct_var_tech, tolerance = tags[[tag]], label = tag)
    }
})

test_that("mean_var bins match Python exactly given Python's mean/var", {
    for (tag in c("pp.sparse.nocov.none", "pp.sparse.cov.none")) {
        ref <- readReference(paste0(tag, ".gene_stats"))
        bins <- suppressMessages(scdRs:::.meanVarBins(ref$mean, ref$var,
            20L, 20L))
        expect_identical(bins, ref$mean_var, label = tag)
    }
})
