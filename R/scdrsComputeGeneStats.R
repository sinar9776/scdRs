#' Compute gene- and cell-level statistics for scDRS
#'
#' Computes the gene-level statistics used to match control genes (mean,
#' variance, technical variance and mean-variance bins) and the cell-level
#' mean and variance. This is the R port of `scdrs.pp.compute_stats()`.
#' Most users call [scdrsPreprocess()], which runs this function.
#'
#' @param sce A [SingleCellExperiment::SingleCellExperiment] with
#'     log-normalized expression (genes in rows, cells in columns).
#' @param assay Name of the assay to use.
#' @param implicit Logical. If `TRUE`, the statistics are computed for the
#'     covariate-corrected data `X + COV_BETA %*% t(COV_MAT) + COV_GENE_MEAN`
#'     without building it (implicit covariate correction). The correction
#'     is read from `reducedDim(sce, "scdrs_cov")`,
#'     `rowData(sce)$scdrs_cov_beta` and `rowData(sce)$scdrs_cov_gene_mean`,
#'     as stored by [scdrsPreprocess()]. If `FALSE`, the statistics are
#'     computed for the assay as it is.
#' @param cellWeight `NULL` or a numeric vector with one weight per cell,
#'     used for the gene-level statistics.
#' @param nMeanBin,nVarBin Number of mean-expression bins and of
#'     variance bins within each mean bin. If either is larger than
#'     `floor(sqrt(nrow(sce) / 10))`, both are set to that value.
#' @param nChunk Number of chunks used in implicit mode (controls memory
#'     use only; the results do not depend on it).
#'
#' @return A list of two [S4Vectors::DataFrame]s:
#' \describe{
#'   \item{`gene`}{One row per gene: `mean`, `var` (log scale), `var_tech`
#'     (technical variance, log scale), `ct_mean`, `ct_var`, `ct_var_tech`
#'     (the same on the `expm1()` scale) and `mean_var` (mean-variance bin,
#'     e.g. `"3.12"`).}
#'   \item{`cell`}{One row per cell: `mean` and `var` (log scale,
#'     unweighted).}
#' }
#'
#' @details Variances use the population formula (divisor n). Technical
#' variance follows the Seurat v3 recipe: a LOESS fit (span 0.3, degree 2)
#' of `log10(ct_var)` on `log10(ct_mean)`, and
#' `var_tech = var * ct_var_tech / ct_var`. Bins are quantile bins as
#' `pandas.qcut()`, labelled from 0.
#'
#' @seealso [scdrsPreprocess()]
#' @references Zhang MJ, Hou K, et al. (2022). Polygenic enrichment
#'     distinguishes disease associations of individual cells in single-cell
#'     RNA-seq data. *Nature Genetics* 54, 1572-1580.
#'
#' @examples
#' sce <- scdrsLoadToyData()
#' stats <- scdrsComputeGeneStats(sce)
#' head(stats$gene)
#' head(stats$cell)
#'
#' @export
scdrsComputeGeneStats <- function(sce, assay = "logcounts", implicit = FALSE,
    cellWeight = NULL, nMeanBin = 20L, nVarBin = 20L, nChunk = 20L) {
    x <- .getAssay(sce, assay)
    if (!is.null(cellWeight) && length(cellWeight) != ncol(x)) {
        stop("'cellWeight' must have one value per cell")
    }
    if (isTRUE(implicit)) {
        covMat <- if ("scdrs_cov" %in%
            SingleCellExperiment::reducedDimNames(sce)) {
            SingleCellExperiment::reducedDim(sce, "scdrs_cov")
        }
        rd <- SummarizedExperiment::rowData(sce)
        covBeta <- rd$scdrs_cov_beta
        geneMean <- rd$scdrs_cov_gene_mean
        if (is.null(covMat) || is.null(covBeta) || is.null(geneMean)) {
            stop("covariate correction not found in 'sce'; ",
                "run scdrsPreprocess() with 'cov' first")
        }
        covBeta <- as.matrix(covBeta)[, colnames(covMat), drop = FALSE]
        gs <- .meanVarImplicit(x, covMat, covBeta, geneMean, "gene",
            cellWeight, nChunk = nChunk)
        ct <- .meanVarImplicit(x, covMat, covBeta, geneMean, "gene",
            cellWeight, fun = expm1, nChunk = nChunk)
        cs <- .meanVarImplicit(x, covMat, covBeta, geneMean, "cell",
            nChunk = nChunk)
    } else {
        gs <- .meanVar(x, "gene", cellWeight)
        ct <- .meanVar(.expm1(x), "gene", cellWeight)
        cs <- .meanVar(x, "cell")
    }
    ctVarTech <- .techVar(ct$mean, ct$var)
    varTech <- gs$var * ctVarTech / ct$var
    varTech[is.na(varTech)] <- 0
    gene <- S4Vectors::DataFrame(
        mean = gs$mean, var = gs$var, var_tech = varTech,
        ct_mean = ct$mean, ct_var = ct$var, ct_var_tech = ctVarTech,
        mean_var = .meanVarBins(gs$mean, gs$var, nMeanBin, nVarBin),
        row.names = rownames(x))
    cell <- S4Vectors::DataFrame(mean = cs$mean, var = cs$var,
        row.names = colnames(x))
    list(gene = gene, cell = cell)
}
