## The toy data is the Python scDRS toy data (toydata_mouse.h5ad) as
## exported by make_reference_outputs.py; compare with Python's
## preprocessed means (pp.sparse.nocov.none: unweighted gene means of X).
test_that("scdrsLoadToyData returns the Python toy data", {
    sce <- scdrsLoadToyData()
    expect_s4_class(sce, "SingleCellExperiment")
    expect_identical(dim(sce), c(2500L, 30L))
    expect_s4_class(assay(sce, "logcounts"), "dgCMatrix")
    expect_identical(as.vector(table(sce$cell_type)), c(10L, 20L))
    ref <- readReference("pp.sparse.nocov.none.gene_stats")
    expect_identical(rownames(sce), rownames(ref))
    expect_equal(unname(Matrix::rowMeans(assay(sce))), ref$mean,
        tolerance = 1e-6)
})
