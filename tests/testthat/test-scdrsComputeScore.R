## scdrsNormalize() compared with scdrs.util.load_h5ad() (filtering,
## normalize_per_cell, log1p; scanpy 1.10.3, float32, hence 1e-6) on a
## simulated count matrix (Python seed 3, fixtures util.normalize.*);
## scdrsComputeScore() compared with the `scdrs compute-score` steps
## (score_cell per trait, with Python's control sets); score files written
## by R compared with the file Python writes. No R random numbers in the
## Python comparisons.

countsSce <- function() {
    counts <- as.matrix(readReference("util.normalize.counts"))
    SingleCellExperiment(list(counts = Matrix::Matrix(counts, sparse = TRUE)))
}

test_that("scdrsNormalize matches load_h5ad with and without filtering", {
    for (flt in c("filter", "nofilter")) {
        ref <- as.matrix(readReference(paste0("util.normalize.", flt)))
        res <- scdrsNormalize(countsSce(), filter = flt == "filter",
            minGenes = 20, minCells = 6)
        expect_identical(dimnames(res), dimnames(ref))
        expect_equal(as.matrix(assay(res, "logcounts")), ref,
            tolerance = 1e-6, label = flt)
        expect_true(is(assay(res, "logcounts"), "sparseMatrix"))
    }
    ## dense counts give the same values
    sce <- countsSce()
    assay(sce, "counts") <- as.matrix(assay(sce, "counts"))
    res <- scdrsNormalize(sce, minGenes = 20, minCells = 6)
    ref <- as.matrix(readReference("util.normalize.filter"))
    expect_equal(assay(res, "logcounts"), ref, tolerance = 1e-6)
    expect_error(scdrsNormalize(sce, name = "counts"), "never overwritten")
    assay(sce, "counts")[1, 1] <- -1
    expect_error(scdrsNormalize(sce), "negative")
})

test_that("scdrsComputeScore runs compute-score's steps", {
    sce <- scdrsLoadToyData()
    gs <- scdrsReadGeneSets(system.file("extdata", "toydata_mouse.gs",
        package = "scdRs"))
    gs$tiny <- c(Pip4k2a = 1, Chd7 = 1)
    set.seed(7)
    expect_message(res <- scdrsComputeScore(sce, gs, cov = "covariate",
        nCtrl = 1000), "tiny: skipped")
    expect_s4_class(res, "SimpleList")
    expect_named(res, "toydata_gs_mouse")
    r <- res$toydata_gs_mouse
    expect_identical(dim(r$ctrl_norm_score), c(30L, 1000L))
    expect_null(r$ctrl_raw_score)
    ## R's own control sets: statistical agreement with Python (1000 sets)
    ref <- readReference("cov.toydata_gs_mouse.score1000")[rownames(r), ]
    expect_gte(cor(r$norm_score, ref$norm_score), 0.99)
    ## the same as scdrsPreprocess() + scdrsScoreCells() with the same sets
    pp <- suppressMessages(scdrsPreprocess(sce, cov = "covariate"))
    r2 <- scdrsScoreCells(pp, gs$toydata_gs_mouse,
        controlSets = metadata(r)$scdrs$controlSets)
    expect_equal(r2$norm_score, r$norm_score)
    ## a preprocessed object is used as it is; a .gs path is accepted
    f <- system.file("extdata", "toydata_mouse.gs", package = "scdRs")
    set.seed(7)
    res3 <- suppressMessages(scdrsComputeScore(pp, f, nCtrl = 1000))
    expect_equal(res3[[1]]$norm_score, r$norm_score)
})

test_that("score files round-trip and match Python's format", {
    s <- scdrsReadScores(referencePath(
        "downstream.toydata_gs_mouse.full_score.gz"))[[1]]
    dir <- tempfile("scores")
    files <- scdrsWriteScores(list(BMI = s), dir)
    expect_identical(basename(files), c("BMI.score.gz", "BMI.full_score.gz"))
    ## same columns, in the same order, as Python's file
    pyHeader <- strsplit(readLines(gzfile(referencePath(
        "downstream.toydata_gs_mouse.full_score.gz")), n = 1), "\t")[[1]]
    rHeader <- strsplit(readLines(gzfile(files[2]), n = 1), "\t")[[1]]
    expect_identical(rHeader[-1], pyHeader[-1])
    back <- scdrsReadScores(file.path(dir, "@.full_score.gz"))
    expect_named(back, "BMI")
    expect_equal(back$BMI, s, ignore_attr = TRUE)
    main <- scdrsReadScores(file.path(dir, "BMI.score.gz"))$BMI
    expect_identical(names(main), c("raw_score", "norm_score", "mc_pval",
        "pval", "nlog10_pval", "zscore"))
    expect_message(skip <- scdrsReadScores(file.path(dir, "@.full_score.gz"),
        cells = paste0("other", 1:10)), "skipped")
    expect_length(skip, 0L)
    expect_error(scdrsReadScores(file.path(dir, "x.csv")), "expect scDRS")
})
