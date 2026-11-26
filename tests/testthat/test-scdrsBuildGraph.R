## scdrsBuildGraph() compared with scanpy's sc.pp.pca(n_comps = 20) +
## sc.pp.neighbors(n_neighbors = 15) as run by `scdrs perform-downstream`
## on the toy data (fixtures downstream.pca, downstream.connectivities;
## scanpy 1.10.3, umap-learn 0.5.12). No random numbers: for 30 cells
## scanpy uses an exact (brute-force) neighbour search. scanpy works in
## float32, hence the tolerances.

pyGraph <- function(cells) {
    e <- readReference("downstream.connectivities", rowNames = FALSE)
    Matrix::sparseMatrix(i = match(e$cell_i, cells),
        j = match(e$cell_j, cells), x = e$weight,
        dims = c(length(cells), length(cells)))
}

test_that("UMAP connectivities from scanpy's PCA match scanpy's graph", {
    pcs <- as.matrix(readReference("downstream.pca"))
    nn <- RscDRS:::.exactKnn(pcs, 15L)
    expect_identical(nn$index[, 1], seq_len(nrow(pcs)))
    g <- RscDRS:::.umapConnectivities(nn$index, nn$distance)
    ref <- pyGraph(rownames(pcs))
    expect_identical(Matrix::isSymmetric(g), TRUE)
    expect_identical(as.matrix(g) > 0, as.matrix(ref) > 0)
    expect_equal(as.matrix(g), as.matrix(ref), tolerance = 1e-5)
})

test_that("scdrsBuildGraph matches scanpy end to end", {
    sce <- scdrsLoadToyData()
    res <- suppressMessages(scdrsBuildGraph(sce))
    ## PCA scores equal up to the sign of each component
    pcs <- reducedDim(res, "scdrs_pca")
    ref <- as.matrix(readReference("downstream.pca"))[colnames(sce), ]
    expect_identical(dim(pcs), c(30L, 20L))
    expect_equal(abs(unname(pcs)), abs(unname(ref)), tolerance = 1e-4)
    g <- colPair(res, "connectivities", asSparse = TRUE)
    refG <- pyGraph(colnames(sce))
    expect_identical(as.matrix(g) > 0, as.matrix(refG) > 0)
    expect_equal(as.matrix(g), as.matrix(refG), tolerance = 1e-4)
    ## nPCs is capped as in the CLI
    expect_message(scdrsBuildGraph(sce, nPCs = 50), "nPCs = 29")
})

test_that(".smoothKnnDist solves the UMAP equation", {
    set.seed(6)
    d <- cbind(0, t(apply(matrix(runif(40 * 9), 40), 1, sort)))
    sm <- RscDRS:::.smoothKnnDist(d, 10)
    psum <- rowSums(pmin(exp(-(d[, -1] - sm$rho) / sm$sigma), 1))
    expect_equal(psum, rep(log2(10), 40), tolerance = 1e-4)
    expect_identical(sm$rho, d[, 2])
})
