## Internal helpers for scoring (port of the subroutines of
## scdrs/method.py: _select_ctrl_geneset, _compute_raw_score,
## _compute_overdispersion_score, _correct_background,
## _get_p_from_empi_null).
##
## Orientation: the expression matrix is genes x cells; score matrices are
## cells x gene sets.

## Preprocessing results needed for scoring, with a helpful error if
## scdrsPreprocess() has not been run.
.scoringSetup <- function(sce) {
    if (!is(sce, "SingleCellExperiment")) {
        stop("'sce' must be a SingleCellExperiment")
    }
    md <- S4Vectors::metadata(sce)$scdrs
    rd <- SummarizedExperiment::rowData(sce)
    need <- paste0("scdrs_", c("mean", "var", "var_tech"))
    if (is.null(md) || !all(need %in% names(rd))) {
        stop("preprocessing results not found in 'sce'; ",
            "run scdrsPreprocess() first")
    }
    x <- .getAssay(sce, md$assay)
    implicit <- isTRUE(md$implicit)
    out <- list(x = x, implicit = implicit)
    if (implicit) {
        covMat <- SingleCellExperiment::reducedDim(sce, "scdrs_cov")
        out$covMat <- covMat
        out$covBeta <- as.matrix(rd$scdrs_cov_beta)[, colnames(covMat),
            drop = FALSE]
        out$geneMean <- rd$scdrs_cov_gene_mean
    }
    out
}

## Gene statistics from rowData(sce) (columns without the "scdrs_" prefix)
## or from a data.frame / DataFrame with genes as row names.
.geneStatsTable <- function(geneStats) {
    if (is(geneStats, "SummarizedExperiment")) {
        rd <- SummarizedExperiment::rowData(geneStats)
        cols <- grep("^scdrs_", names(rd), value = TRUE)
        cols <- setdiff(cols, c("scdrs_cov_beta", "scdrs_cov_gene_mean"))
        df <- as.data.frame(rd[, cols, drop = FALSE], optional = TRUE)
        names(df) <- sub("^scdrs_", "", cols)
        rownames(df) <- rownames(geneStats)
        return(df)
    }
    if (is.data.frame(geneStats) || is(geneStats, "DataFrame")) {
        if (is.null(rownames(geneStats))) {
            stop("'geneStats' must have gene names as row names")
        }
        return(as.data.frame(geneStats, optional = TRUE))
    }
    stop("'geneStats' must be a SingleCellExperiment or a data.frame")
}

## Disease genes and weights as score_cell() prepares them: weights default
## to 1; a gene given twice keeps its last weight (Python dict); genes not
## in `universe` are dropped; genes are sorted by code point (Python
## sorted()). `genes` may be a named numeric vector (names = genes, values
## = weights). Returns a named numeric vector.
.prepGenes <- function(genes, weights = NULL, universe) {
    if (is.numeric(genes) && !is.null(names(genes))) {
        if (!is.null(weights)) {
            stop("give weights either as the values of 'genes' or in ",
                "'weights', not both")
        }
        weights <- unname(genes)
        genes <- names(genes)
    }
    if (!is.character(genes)) {
        stop("'genes' must be a character vector or a named numeric vector")
    }
    if (is.null(weights)) {
        weights <- rep(1, length(genes))
    }
    if (length(weights) != length(genes)) {
        stop("'weights' must have one value per gene")
    }
    last <- !duplicated(genes, fromLast = TRUE)
    w <- stats::setNames(as.numeric(weights[last]), genes[last])
    keep <- sort(intersect(names(w), universe), method = "radix")
    w[keep]
}

## Control gene sets given by the user: a list of named numeric vectors
## (names = genes, values = weights), or a data.frame in the long format of
## the Python reference files (columns ctrl_id, gene, weight).
.asControlSets <- function(controlSets, universe) {
    if (is.data.frame(controlSets) || is(controlSets, "DataFrame")) {
        df <- as.data.frame(controlSets)
        if (!all(c("ctrl_id", "gene", "weight") %in% names(df))) {
            stop("a control set table needs columns ctrl_id, gene, weight")
        }
        ids <- sort(unique(df$ctrl_id))
        controlSets <- lapply(ids, function(i) {
            s <- df[df$ctrl_id == i, , drop = FALSE]
            stats::setNames(as.numeric(s$weight), as.character(s$gene))
        })
    }
    if (!is.list(controlSets) || length(controlSets) == 0L) {
        stop("'controlSets' must be a non-empty list or a data.frame")
    }
    ok <- vapply(controlSets, function(s) {
        is.numeric(s) && !is.null(names(s)) && all(names(s) %in% universe)
    }, logical(1))
    if (!all(ok)) {
        stop("every control set must be a named numeric vector of genes ",
            "found in 'sce'")
    }
    unname(controlSets)
}

## Gene bins for matching control genes (part of _select_ctrl_geneset):
## a key with fewer distinct values than 10% of the genes is used as it is
## (categorical, e.g. "mean_var"); otherwise genes are put in nGeneBin
## quantile bins (pd.qcut). Returns one bin label per gene and the bins in
## the order Python iterates over them (sorted).
.controlBins <- function(key, nGeneBin) {
    if (length(unique(key)) < length(key) / 10) {
        bins <- key
    } else {
        bins <- .qcut(as.numeric(key), nGeneBin)
    }
    order <- if (is.character(bins)) {
        sort(unique(bins), method = "radix")
    } else {
        sort(unique(bins))
    }
    list(bin = bins, order = order)
}

## Sparse genes x sets matrix of normalized score weights (part of
## _compute_raw_score / _compute_overdispersion_score).
## sets: list of named numeric vectors (gene weights); stats: gene table
## with var, var_tech; genes: rownames of the expression matrix.
.scoreWeightMatrix <- function(sets, stats, genes, weightOpt) {
    idx <- lapply(sets, function(s) match(names(s), genes))
    vals <- lapply(seq_along(sets), function(k) {
        i <- idx[[k]]
        u <- switch(weightOpt,
            uniform = rep(1, length(i)),
            vs = 1 / sqrt(stats$var_tech[i] + 1e-2),
            inv_std = 1 / sqrt(stats$var[i] + 1e-2),
            od = 1 / (stats$var_tech[i] + 1e-2))
        w <- u * unname(sets[[k]])
        w / sum(w)
    })
    Matrix::sparseMatrix(i = unlist(idx), j = rep(seq_along(sets),
        lengths(idx)), x = unlist(vals), dims = c(length(genes),
        length(sets)))
}

## Raw scores X~' W for all sets at once (cells x sets). In implicit mode
## X~ = X + covBeta covMat' + geneMean, so
## X~' W = X' W + covMat (covBeta' W) + 1 (geneMean' W).
.rawScores <- function(dat, W) {
    s <- as.matrix(Matrix::crossprod(dat$x, W))
    if (dat$implicit) {
        s <- s + dat$covMat %*% as.matrix(Matrix::crossprod(dat$covBeta, W))
        s <- s + rep(as.vector(Matrix::crossprod(W, dat$geneMean)),
            each = nrow(s))
    }
    unname(s)
}

## Overdispersion scores (weightOpt = "od", a benchmarking option of
## scDRS): s_c = sum_g w_g [(x~_cg - mean_g)^2 - var_tech_g], computed as
## quadratic + linear + constant terms as in Python. In implicit mode the
## corrected data is built in chunks of genes.
.odScores <- function(dat, W, stats) {
    mu <- stats$mean
    vt <- stats$var_tech
    lin <- W * (2 * mu)
    const <- as.vector(Matrix::crossprod(W, mu^2 - vt))
    if (!dat$implicit) {
        s <- as.matrix(Matrix::crossprod(dat$x^2, W)) -
            as.matrix(Matrix::crossprod(dat$x, lin))
    } else {
        used <- which(Matrix::rowSums(W != 0) > 0)
        s <- matrix(0, ncol(dat$x), ncol(W))
        tCov <- t(dat$covMat)
        for (start in seq(1L, length(used), by = 500L)) {
            g <- used[seq.int(start, min(start + 499L, length(used)))]
            xt <- as.matrix(dat$x[g, , drop = FALSE]) +
                dat$covBeta[g, , drop = FALSE] %*% tCov + dat$geneMean[g]
            s <- s + as.matrix(Matrix::crossprod(xt^2, W[g, , drop = FALSE])) -
                as.matrix(Matrix::crossprod(xt, lin[g, , drop = FALSE]))
        }
    }
    unname(s + rep(const, each = nrow(s)))
}

## Port of _correct_background(): first gene-set alignment, cell-wise
## standardization (population SD), second gene-set alignment, and cells
## with a raw score of exactly 0 set to the minimum. rawScore: length n;
## ctrlRawScore: n x K; varRatio: length K. With returnIntermediate = TRUE
## the four stages are returned too (for the tests).
.correctBackground <- function(rawScore, ctrlRawScore, varRatio,
    returnIntermediate = FALSE) {
    zero <- rawScore == 0
    zeroCtrl <- ctrlRawScore == 0
    s <- rawScore - mean(rawScore)
    S <- sweep(ctrlRawScore, 2L, colMeans(ctrlRawScore))
    S <- sweep(S, 2L, sqrt(varRatio), "/")
    stage1 <- list(score = s, ctrl = S)
    m <- rowMeans(S)
    sd0 <- sqrt(rowMeans((S - m)^2))
    z <- (s - m) / sd0
    Z <- (S - m) / sd0
    stage2 <- list(score = z, ctrl = Z)
    z <- z - mean(z)
    Z <- sweep(Z, 2L, colMeans(Z))
    stage3 <- list(score = z, ctrl = Z)
    minScore <- min(min(z), min(Z))
    z[zero] <- minScore - 1e-3
    Z[zeroCtrl] <- minScore
    out <- list(score = z, ctrl = Z)
    if (returnIntermediate) {
        out$intermediate <- list(firstAlignment = stage1,
            cellwiseStandardization = stage2, secondAlignment = stage3)
    }
    out
}

## Port of _get_p_from_empi_null(): p_i = (1 + #{null >= t_i}) / (1 + N),
## in O(N log N).
.empiricalPvalue <- function(t, tNull) {
    n <- length(tNull)
    pos <- findInterval(t, sort(tNull), left.open = TRUE)
    (n - pos + 1) / (n + 1)
}
