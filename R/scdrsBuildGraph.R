## kNN graph for the heterogeneity test, as `scdrs perform-downstream`
## builds it with scanpy: sc.pp.pca(n_comps) + sc.pp.neighbors(n_neighbors)
## (exact Euclidean neighbours, UMAP connectivities).

## PCA scores (cells x nPCs) of the gene-centered data, as sc.pp.pca with
## zero_center = TRUE, via BiocSingular: exact SVD for small data, IRLBA
## (deferred centering, keeps sparse data sparse) otherwise. Signs of the
## components are arbitrary; distances between cells do not depend on them.
.pcaScores <- function(x, nPCs) {
    xt <- Matrix::t(x)
    bs <- if (as.numeric(nrow(xt)) * ncol(xt) <= 2e7) {
        BiocSingular::ExactParam()
    } else {
        BiocSingular::IrlbaParam(deferred = TRUE)
    }
    pca <- BiocSingular::runPCA(xt, rank = nPCs, center = TRUE,
        scale = FALSE, get.rotation = FALSE, BSPARAM = bs)
    unname(pca$x[, seq_len(nPCs), drop = FALSE])
}

## Exact k nearest neighbours (Euclidean) with BiocNeighbors, the cell
## itself first, as scanpy's neighbour table. Returns n x k matrices of
## indices and distances.
.exactKnn <- function(pcs, k) {
    nn <- BiocNeighbors::findKNN(pcs, k = k - 1L,
        BNPARAM = BiocNeighbors::KmknnParam())
    list(index = cbind(seq_len(nrow(pcs)), nn$index),
        distance = cbind(0, nn$distance))
}

## Port of umap.umap_.smooth_knn_dist() with local_connectivity = 1 and
## bandwidth = 1, vectorized over cells: rho = distance to the nearest
## neighbour at a positive distance; sigma found by binary search (64
## steps) so that sum_j exp(-(d_j - rho) / sigma) = log2(k), over the
## neighbours other than the cell itself.
.smoothKnnDist <- function(distance, k, nIter = 64L) {
    n <- nrow(distance)
    target <- log2(k)
    pos <- distance
    pos[pos <= 0] <- Inf
    rho <- apply(pos, 1L, min)
    rho[!is.finite(rho)] <- 0
    d <- distance[, -1L, drop = FALSE] - rho
    lo <- numeric(n)
    hi <- rep(Inf, n)
    mid <- rep(1, n)
    active <- rep(TRUE, n)
    for (it in seq_len(nIter)) {
        if (!any(active)) {
            break
        }
        a <- which(active)
        e <- ifelse(d[a, , drop = FALSE] > 0,
            exp(-(d[a, , drop = FALSE] / mid[a])), 1)
        psum <- rowSums(e)
        done <- abs(psum - target) < 1e-5
        up <- !done & psum > target
        down <- !done & !up
        hi[a[up]] <- mid[a[up]]
        mid[a[up]] <- (lo[a[up]] + hi[a[up]]) / 2
        lo[a[down]] <- mid[a[down]]
        inf <- a[down][!is.finite(hi[a[down]])]
        fin <- a[down][is.finite(hi[a[down]])]
        mid[inf] <- mid[inf] * 2
        mid[fin] <- (lo[fin] + hi[fin]) / 2
        active[a[done]] <- FALSE
    }
    sigma <- mid
    rowMean <- rowMeans(distance)
    floorVal <- ifelse(rho > 0, 1e-3 * rowMean, 1e-3 * mean(distance))
    sigma <- pmax(sigma, floorVal)
    list(sigma = sigma, rho = rho)
}

## Port of umap's fuzzy_simplicial_set() from a kNN result (as scanpy's
## _connectivity.umap): membership strengths exp(-(d - rho) / sigma)
## (1 within rho, 0 for the cell itself), then the fuzzy union
## P + P' - P * P'.
.umapConnectivities <- function(index, distance) {
    n <- nrow(index)
    k <- ncol(index)
    sm <- .smoothKnnDist(distance, k)
    dr <- distance - sm$rho
    val <- exp(-(dr / sm$sigma))
    val[dr <= 0 | sm$sigma == 0] <- 1
    self <- index == seq_len(n)
    val[self] <- 0
    keep <- val != 0
    p <- Matrix::sparseMatrix(i = row(index)[keep], j = index[keep],
        x = val[keep], dims = c(n, n))
    pt <- Matrix::t(p)
    g <- p + pt - p * pt
    Matrix::drop0(g)
}

#' Build a cell-cell kNN graph
#'
#' Builds the nearest-neighbour graph that the scDRS heterogeneity test
#' ([scdrsTestGroups()]) uses, in the same way as `scdrs
#' perform-downstream` does with scanpy (`sc.pp.pca(n_comps = 20)` and
#' `sc.pp.neighbors(n_neighbors = 15)`): a PCA of the gene-centered
#' expression, exact Euclidean nearest neighbours in PCA space, and UMAP
#' connectivities (fuzzy simplicial set) as edge weights.
#'
#' @param sce A [SingleCellExperiment::SingleCellExperiment].
#' @param assay Expression assay to use (log-normalized).
#' @param nNeighbors Number of neighbours per cell, the cell itself
#'     included (as in scanpy).
#' @param nPCs Number of principal components; at most
#'     `min(dim(sce)) - 1`.
#' @param name Name under which the graph is stored in `colPair(sce)`.
#'
#' @return `sce` with the graph in `colPair(sce, name)` (a symmetric
#'     weighted graph) and the PCA in `reducedDim(sce, "scdrs_pca")`.
#'
#' @details The PCA uses \pkg{BiocSingular} (exact SVD for up to 2e7
#' matrix entries, IRLBA otherwise) and the exact neighbour search uses
#' \pkg{BiocNeighbors}; the edge weights are UMAP connectivities, as in
#' scanpy and scDRS. Any other graph, e.g. from \pkg{bluster} or
#' \pkg{scran}, can be used instead by storing it in `colPair(sce)` or
#' passing it to [scdrsTestGroups()]. The graph is close to scanpy's but not always
#' identical, because PCA and neighbour search differ slightly between
#' implementations.
#'
#' @examples
#' sce <- scdrsBuildGraph(scdrsLoadToyData())
#' sce
#' g <- SingleCellExperiment::colPair(sce, "connectivities", asSparse = TRUE)
#' summary(Matrix::rowSums(g > 0))
#'
#' @export
scdrsBuildGraph <- function(sce, assay = "logcounts", nNeighbors = 15L,
    nPCs = 20L, name = "connectivities") {
    x <- .getAssay(sce, assay)
    nMax <- min(dim(x)) - 1L
    if (nPCs > nMax) {
        message("'nPCs' = ", nPCs, " > min(dim(sce)) - 1 = ", nMax,
            "; using nPCs = ", nMax)
        nPCs <- nMax
    }
    nNeighbors <- min(as.integer(nNeighbors), ncol(x))
    if (nNeighbors < 2L) {
        stop("'nNeighbors' must be at least 2")
    }
    pcs <- .pcaScores(x, nPCs)
    rownames(pcs) <- colnames(x)
    colnames(pcs) <- paste0("PC", seq_len(nPCs))
    nn <- .exactKnn(pcs, nNeighbors)
    g <- .umapConnectivities(nn$index, nn$distance)
    SingleCellExperiment::reducedDim(sce, "scdrs_pca") <- pcs
    SingleCellExperiment::colPair(sce, name) <- unname(g)
    sce
}
