## Names of everything scdrsPreprocess() stores, so a re-run starts clean.
.ppGeneCols <- paste0("scdrs_", c("mean", "var", "var_tech", "ct_mean",
    "ct_var", "ct_var_tech", "mean_var", "cov_beta", "cov_gene_mean"))
.ppCellCols <- c("scdrs_mean", "scdrs_var")

.clearPreprocess <- function(sce) {
    rd <- SummarizedExperiment::rowData(sce)
    rd <- rd[, setdiff(names(rd), .ppGeneCols), drop = FALSE]
    SummarizedExperiment::rowData(sce) <- rd
    cd <- SummarizedExperiment::colData(sce)
    cd <- cd[, setdiff(names(cd), .ppCellCols), drop = FALSE]
    SummarizedExperiment::colData(sce) <- cd
    if ("scdrs_cov" %in% SingleCellExperiment::reducedDimNames(sce)) {
        SingleCellExperiment::reducedDim(sce, "scdrs_cov") <- NULL
    }
    if ("scdrs_corrected" %in% SummarizedExperiment::assayNames(sce)) {
        SummarizedExperiment::assay(sce, "scdrs_corrected") <- NULL
    }
    S4Vectors::metadata(sce)$scdrs <- NULL
    sce
}

#' Preprocess single-cell data for scDRS
#'
#' Regresses covariates out of the expression data (adding back each gene's
#' mean), then computes the gene- and cell-level statistics that scDRS uses
#' to choose control genes. This is the R port of `scdrs.preprocess()`.
#'
#' @param sce A [SingleCellExperiment::SingleCellExperiment]. The assay
#'     must be size-factor normalized and `log1p`-transformed (natural log),
#'     with genes in rows and cells in columns. Gene and cell names must be
#'     unique.
#' @param assay Name of the assay to use. It is never modified.
#' @param cov Covariates to regress out, or `NULL` for none. Either a
#'     character vector of `colData(sce)` column names, or a `data.frame`
#'     / `DataFrame` with cell names as row names (like the scDRS `.cov`
#'     file; at least 75% of the cells must be found). Non-numeric columns
#'     are turned into dummy variables, missing values are replaced by the
#'     column mean, and a constant term `SCDRS_CONST` is added unless the
#'     covariates already contain one.
#' @param adjProp `NULL`, or the name of a `colData(sce)` column with cell
#'     groups (e.g. cell types). Cells are then weighted inversely to the
#'     size of their group when computing gene statistics (groups smaller
#'     than 1% of the largest one count as that size).
#' @param nMeanBin,nVarBin Number of mean-expression bins and of variance
#'     bins within each mean bin, used to match control genes.
#' @param nChunk Number of chunks for implicit covariate correction
#'     (memory use only). `NULL` uses `5 * nGene * nCell / nnz + 1` in
#'     implicit mode and 20 otherwise, as in Python.
#' @param implicit Use implicit covariate correction? `NULL` (default)
#'     follows scDRS: `TRUE` when the assay is sparse and `cov` is given.
#'     In implicit mode the corrected data is never stored, which keeps the
#'     memory use low; otherwise it is stored as a new dense assay
#'     `"scdrs_corrected"`.
#'
#' @return `sce`, with:
#' \describe{
#'   \item{`rowData(sce)`}{`scdrs_mean`, `scdrs_var`, `scdrs_var_tech`,
#'     `scdrs_ct_mean`, `scdrs_ct_var`, `scdrs_ct_var_tech`,
#'     `scdrs_mean_var` (see [scdrsComputeGeneStats()]); with covariates
#'     also `scdrs_cov_beta` (a genes x covariates matrix, the negated
#'     regression coefficients, Python's `COV_BETA`) and
#'     `scdrs_cov_gene_mean`.}
#'   \item{`colData(sce)`}{`scdrs_mean`, `scdrs_var`.}
#'   \item{`reducedDim(sce, "scdrs_cov")`}{The covariate matrix (cells x
#'     covariates), with covariates only.}
#'   \item{`assay(sce, "scdrs_corrected")`}{The covariate-corrected data,
#'     with covariates and `implicit = FALSE` only.}
#'   \item{`metadata(sce)$scdrs`}{A list of flags and settings (`assay` is
#'     the assay that later steps should score).}
#' }
#' Results of an earlier run of `scdrsPreprocess()` are replaced.
#'
#' @details The corrected data is
#' \deqn{\tilde X = X - C\beta + \mu,}{X~ = X - C beta + mu,}
#' where \eqn{C} is the covariate matrix, \eqn{\beta} the least-squares
#' coefficients and \eqn{\mu} the gene means. Python's `COV_BETA` is
#' \eqn{-\beta^\top}, so
#' `X~ = X + scdrs_cov_beta %*% t(scdrs_cov) + scdrs_cov_gene_mean`.
#'
#' @seealso [scdrsComputeGeneStats()]
#' @references Zhang MJ, Hou K, et al. (2022). Polygenic enrichment
#'     distinguishes disease associations of individual cells in single-cell
#'     RNA-seq data. *Nature Genetics* 54, 1572-1580.
#'
#' @examples
#' sce <- scdrsLoadToyData()
#' sce <- scdrsPreprocess(sce, cov = "covariate")
#' rowData(sce)[1:3, c("scdrs_mean", "scdrs_var", "scdrs_mean_var")]
#' head(reducedDim(sce, "scdrs_cov"))
#' str(metadata(sce)$scdrs)
#'
#' @importFrom methods as is
#' @export
scdrsPreprocess <- function(sce, assay = "logcounts", cov = NULL,
    adjProp = NULL, nMeanBin = 20L, nVarBin = 20L, nChunk = NULL,
    implicit = NULL) {
    x <- .getAssay(sce, assay)
    nGene <- nrow(x)
    nCell <- ncol(x)
    flagSparse <- is(x, "sparseMatrix")
    flagCov <- !is.null(cov)
    flagAdjProp <- !is.null(adjProp)
    if (is.null(implicit)) {
        implicit <- flagSparse && flagCov
    }
    if (isTRUE(implicit) && !flagCov) {
        stop("'implicit = TRUE' needs covariates ('cov')")
    }
    sce <- .clearPreprocess(sce)
    statsAssay <- assay

    if (flagCov) {
        covMat <- .covariateMatrix(sce, cov)
        geneMean <- unname(Matrix::rowMeans(x))
        beta <- solve(crossprod(covMat) / nCell,
            t(as.matrix(x %*% covMat)) / nCell)
        covBeta <- -t(beta)
        dimnames(covBeta) <- list(rownames(x), colnames(covMat))
        SingleCellExperiment::reducedDim(sce, "scdrs_cov") <- covMat
        rd <- SummarizedExperiment::rowData(sce)
        rd$scdrs_cov_beta <- covBeta
        rd$scdrs_cov_gene_mean <- geneMean
        SummarizedExperiment::rowData(sce) <- rd
        if (!implicit) {
            corrected <- as.matrix(x) + covBeta %*% t(covMat) + geneMean
            dimnames(corrected) <- dimnames(x)
            SummarizedExperiment::assay(sce, "scdrs_corrected") <- corrected
            statsAssay <- "scdrs_corrected"
        }
    }

    if (is.null(nChunk)) {
        nnz <- if (flagSparse) length(x@x) else sum(x != 0)
        nChunk <- if (implicit) {
            floor(5 * nCell * nGene / max(nnz, 1)) + 1
        } else {
            20L
        }
    }
    cellWeight <- if (flagAdjProp) .cellWeights(sce, adjProp)

    st <- scdrsComputeGeneStats(sce, assay = statsAssay,
        implicit = implicit, cellWeight = cellWeight, nMeanBin = nMeanBin,
        nVarBin = nVarBin, nChunk = nChunk)
    rd <- SummarizedExperiment::rowData(sce)
    for (nm in names(st$gene)) {
        rd[[paste0("scdrs_", nm)]] <- st$gene[[nm]]
    }
    SummarizedExperiment::rowData(sce) <- rd
    cd <- SummarizedExperiment::colData(sce)
    for (nm in names(st$cell)) {
        cd[[paste0("scdrs_", nm)]] <- st$cell[[nm]]
    }
    SummarizedExperiment::colData(sce) <- cd

    S4Vectors::metadata(sce)$scdrs <- list(
        flagSparse = flagSparse, flagCov = flagCov,
        flagAdjProp = flagAdjProp, implicit = implicit,
        assay = statsAssay, inputAssay = assay, adjProp = adjProp,
        nMeanBin = nMeanBin, nVarBin = nVarBin, nChunk = nChunk)
    sce
}
