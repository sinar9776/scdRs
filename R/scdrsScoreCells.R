#' Select matched control gene sets
#'
#' Draws `nCtrl` control gene sets that match the disease gene set gene by
#' gene: within each bin of `ctrlMatchKey` (by default the mean-variance
#' bins made by [scdrsPreprocess()]) that contains \eqn{d} disease genes,
#' each control set gets \eqn{d} genes drawn at random, without
#' replacement, from all genes of that bin (disease genes included). A
#' control gene takes the weight of a disease gene of the same bin. This is
#' the R port of `scdrs.method._select_ctrl_geneset()`.
#'
#' Random numbers come from R's generator, so call [set.seed()] first for
#' reproducible results. The control sets are not the same as the ones
#' Python draws with the same seed (see the vignette).
#'
#' @param geneStats A [SingleCellExperiment::SingleCellExperiment]
#'     processed by [scdrsPreprocess()] (the statistics are taken from
#'     `rowData`, without the `scdrs_` prefix), or a `data.frame` /
#'     `DataFrame` of gene statistics with gene names as row names.
#' @param genes Disease genes: a character vector, or a named numeric
#'     vector of weights (names = genes), as returned by
#'     [scdrsReadGeneSets()]. Genes not in `geneStats` are dropped; a gene
#'     given twice keeps its last weight.
#' @param weights Optional gene weights when `genes` is a character vector
#'     (default 1).
#' @param ctrlMatchKey Gene statistic used for matching, e.g. `"mean_var"`
#'     (default) or `"mean"`. A statistic with fewer distinct values than
#'     10% of the genes is used as categories; otherwise genes are split
#'     into `nGeneBin` quantile bins.
#' @param nCtrl Number of control gene sets.
#' @param nGeneBin Number of quantile bins for a continuous `ctrlMatchKey`.
#'
#' @return A list of `nCtrl` named numeric vectors (names = control genes,
#'     values = weights).
#'
#' @seealso [scdrsScoreCells()]
#' @references Zhang MJ, Hou K, et al. (2022). Polygenic enrichment
#'     distinguishes disease associations of individual cells in single-cell
#'     RNA-seq data. *Nature Genetics* 54, 1572-1580.
#'
#' @examples
#' sce <- scdrsPreprocess(scdrsLoadToyData())
#' gs <- scdrsReadGeneSets(system.file("extdata", "toydata_mouse.gs",
#'     package = "RscDRS"))
#' set.seed(1)
#' ctrl <- scdrsSelectControlGenes(sce, gs[[1]], nCtrl = 5)
#' lengths(ctrl)
#' head(ctrl[[1]])
#'
#' @export
scdrsSelectControlGenes <- function(geneStats, genes, weights = NULL,
    ctrlMatchKey = "mean_var", nCtrl = 1000L, nGeneBin = 200L) {
    stats <- .geneStatsTable(geneStats)
    if (!is.character(ctrlMatchKey) || length(ctrlMatchKey) != 1L ||
        !ctrlMatchKey %in% names(stats)) {
        stop("'ctrlMatchKey' = ", paste(ctrlMatchKey, collapse = ", "),
            " not found in the gene statistics")
    }
    nCtrl <- as.integer(nCtrl)
    if (length(nCtrl) != 1L || is.na(nCtrl) || nCtrl < 1L) {
        stop("'nCtrl' must be a positive integer")
    }
    universe <- rownames(stats)
    w <- .prepGenes(genes, weights, universe)
    if (length(w) == 0L) {
        stop("none of the disease genes are in the gene statistics")
    }
    bins <- .controlBins(stats[[ctrlMatchKey]], nGeneBin)
    draws <- list()
    ctrlWeights <- list()
    for (b in bins$order) {
        binGenes <- sort(universe[!is.na(bins$bin) & bins$bin == b],
            method = "radix")
        binDisease <- binGenes[binGenes %in% names(w)]
        d <- length(binDisease)
        if (d == 0L) {
            next
        }
        nb <- length(binGenes)
        m <- vapply(seq_len(nCtrl), function(k) sample.int(nb, d),
            integer(d))
        draws[[length(draws) + 1L]] <- matrix(binGenes[m], nrow = d)
        ctrlWeights[[length(ctrlWeights) + 1L]] <- unname(w[binDisease])
    }
    draws <- do.call(rbind, draws)
    ctrlWeights <- unlist(ctrlWeights)
    lapply(seq_len(nCtrl), function(k) {
        stats::setNames(ctrlWeights, draws[, k])
    })
}

#' Compute scDRS disease scores for each cell
#'
#' For one disease gene set, computes each cell's raw disease score (a
#' weighted mean of the covariate-corrected expression of the disease
#' genes), the same score for `nCtrl` matched control gene sets, the
#' normalized score, and a cell-level p-value from the empirical
#' distribution of the control scores pooled over all cells. This is the R
#' port of `scdrs.score_cell()`.
#'
#' Run [scdrsPreprocess()] first. Control gene sets are drawn at random
#' (see [scdrsSelectControlGenes()]); call [set.seed()] first for
#' reproducible results.
#'
#' @param sce A [SingleCellExperiment::SingleCellExperiment] processed by
#'     [scdrsPreprocess()].
#' @param genes Disease genes: a character vector, or a named numeric
#'     vector of weights (names = genes), e.g. one element of
#'     [scdrsReadGeneSets()]. Genes not in `sce` are dropped.
#' @param weights Optional gene weights when `genes` is a character vector
#'     (default 1).
#' @param ctrlMatchKey Gene statistic used to match control genes
#'     (`"mean_var"`, the default, or e.g. `"mean"`); see
#'     [scdrsSelectControlGenes()].
#' @param nCtrl Number of control gene sets.
#' @param nGeneBin Number of quantile bins for a continuous `ctrlMatchKey`.
#' @param weightOpt How the raw score weights genes: `"vs"` (default,
#'     \eqn{1/\sqrt{\mathrm{var\_tech} + 0.01}}{1/sqrt(var_tech + 0.01)}),
#'     `"uniform"`, `"inv_std"`
#'     (\eqn{1/\sqrt{\mathrm{var} + 0.01}}{1/sqrt(var + 0.01)}), or `"od"`
#'     (overdispersion score, a benchmarking option); always multiplied by
#'     the gene weights.
#' @param returnCtrlRawScore,returnCtrlNormScore Also return the raw /
#'     normalized control scores? The normalized control scores are needed
#'     by the downstream analyses ([scdrsTestGroups()] etc.).
#' @param controlSets Optional precomputed control gene sets, either a
#'     list of named numeric vectors as returned by
#'     [scdrsSelectControlGenes()], or a `data.frame` with columns
#'     `ctrl_id`, `gene` and `weight`. If given, `nCtrl`, `ctrlMatchKey`
#'     and `nGeneBin` only matter for the variance ratio (computed when
#'     `ctrlMatchKey = "mean_var"`), and no random numbers are used.
#'
#' @return A [S4Vectors::DataFrame] with one row per cell (row names = cell
#'     names) and columns
#' \describe{
#'   \item{`raw_score`}{Raw disease score.}
#'   \item{`norm_score`}{Normalized disease score.}
#'   \item{`mc_pval`}{Monte Carlo p-value from the cell's own control
#'     scores.}
#'   \item{`pval`}{scDRS cell-level p-value (control scores pooled over all
#'     cells).}
#'   \item{`nlog10_pval`}{`-log10(pval)`.}
#'   \item{`zscore`}{One-sided z-score of `pval`, clipped to
#'     \[-10, 10\].}
#'   \item{`ctrl_raw_score`, `ctrl_norm_score`}{Only if requested: cells x
#'     `nCtrl` matrices of control scores (matrix columns, columns named
#'     `0`, `1`, ... as in Python's `ctrl_norm_score_0`, ...).}
#' }
#' `metadata()` of the result holds the disease genes and weights used
#' (`genes`), the control sets (`controlSets`) and the settings.
#'
#' @details Let \eqn{s_c} be the raw disease score of cell \eqn{c} and
#' \eqn{S_{ck}} that of control set \eqn{k}. The scores are centered per
#' gene set, the control scores are scaled to the same variance as the
#' disease score (assuming independent genes), each cell is standardized
#' by the mean and SD of its control scores, and the scores are centered
#' again. Cells with a raw score of exactly 0 get the smallest score. The
#' p-value is \eqn{(1 + \#\{(c',k): Z_{c'k} \ge z_c\}) / (1 + nK)}.
#'
#' @seealso [scdrsPreprocess()], [scdrsComputeScore()] to score several
#'     gene sets, [scdrsTestGroups()].
#' @references Zhang MJ, Hou K, et al. (2022). Polygenic enrichment
#'     distinguishes disease associations of individual cells in single-cell
#'     RNA-seq data. *Nature Genetics* 54, 1572-1580.
#'
#' @examples
#' sce <- scdrsPreprocess(scdrsLoadToyData(), cov = "covariate")
#' gs <- scdrsReadGeneSets(system.file("extdata", "toydata_mouse.gs",
#'     package = "RscDRS"))
#' set.seed(1)
#' res <- scdrsScoreCells(sce, gs[["toydata_gs_mouse"]], nCtrl = 100)
#' head(res)
#' tapply(res$norm_score, sce$cell_type, mean)
#'
#' @export
scdrsScoreCells <- function(sce, genes, weights = NULL,
    ctrlMatchKey = "mean_var", nCtrl = 1000L, nGeneBin = 200L,
    weightOpt = c("vs", "uniform", "inv_std", "od"),
    returnCtrlRawScore = FALSE, returnCtrlNormScore = FALSE,
    controlSets = NULL) {
    weightOpt <- match.arg(weightOpt)
    dat <- .scoringSetup(sce)
    stats <- .geneStatsTable(sce)
    universe <- rownames(sce)
    w <- .prepGenes(genes, weights, universe)
    if (length(w) == 0L) {
        stop("none of the disease genes are in 'sce'")
    }
    if (is.null(controlSets)) {
        controlSets <- scdrsSelectControlGenes(stats, w,
            ctrlMatchKey = ctrlMatchKey, nCtrl = nCtrl, nGeneBin = nGeneBin)
    } else {
        controlSets <- .asControlSets(controlSets, universe)
    }
    nCtrl <- length(controlSets)
    W <- .scoreWeightMatrix(c(list(w), controlSets), stats, universe,
        weightOpt)
    raw <- if (weightOpt == "od") {
        .odScores(dat, W, stats)
    } else {
        .rawScores(dat, W)
    }
    varRatio <- .varRatio(W, stats, ctrlMatchKey, weightOpt)
    bg <- .correctBackground(raw[, 1L], raw[, -1L, drop = FALSE], varRatio)
    res <- .scoreTable(raw, bg, colnames(sce), returnCtrlRawScore,
        returnCtrlNormScore)
    S4Vectors::metadata(res)$scdrs <- list(genes = w,
        controlSets = controlSets, nCtrl = nCtrl,
        ctrlMatchKey = ctrlMatchKey, weightOpt = weightOpt)
    res
}

## Variance of each control score relative to the disease score
## (mean_var matching, not od); 1 otherwise.
.varRatio <- function(W, stats, ctrlMatchKey, weightOpt) {
    if (ctrlMatchKey != "mean_var" || weightOpt == "od") {
        return(rep(1, ncol(W) - 1L))
    }
    v <- as.vector(Matrix::crossprod(W^2, stats$var))
    v[-1L] / v[1L]
}

## The result table of scdrsScoreCells() from the raw scores (disease
## first) and the background-corrected scores.
.scoreTable <- function(raw, bg, cells, returnRaw, returnNorm) {
    z <- bg$score
    Z <- bg$ctrl
    nCtrl <- ncol(Z)
    mcp <- (1 + rowSums(Z >= z, na.rm = TRUE)) / (1 + nCtrl)
    p <- .empiricalPvalue(z, as.vector(Z))
    res <- S4Vectors::DataFrame(raw_score = raw[, 1L], norm_score = z,
        mc_pval = mcp, pval = p, nlog10_pval = -log10(p),
        zscore = -pmin(pmax(stats::qnorm(p), -10), 10),
        row.names = cells)
    ctrlNames <- as.character(seq_len(nCtrl) - 1L)
    if (returnRaw) {
        res$ctrl_raw_score <- matrix(raw[, -1L], ncol = nCtrl,
            dimnames = list(NULL, ctrlNames))
    }
    if (returnNorm) {
        res$ctrl_norm_score <- matrix(Z, ncol = nCtrl,
            dimnames = list(NULL, ctrlNames))
    }
    res
}
