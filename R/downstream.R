## Downstream analyses (port of scdrs/method.py: downstream_group_analysis,
## downstream_corr_analysis, downstream_gene_analysis, test_gearysc,
## gearys_c, _pearson_corr, _pearson_corr_sparse; and of the
## `scdrs perform-downstream` command).

## Population standard deviation (numpy's default ddof = 0).
.sd0 <- function(x) {
    sqrt(mean((x - mean(x))^2))
}

## Normalized control scores (cells x K) from a score table: the matrix
## column `ctrl_norm_score` of scdrsScoreCells(), or the columns
## ctrl_norm_score_0, ctrl_norm_score_1, ... of a Python .full_score file.
.ctrlNormScores <- function(scores) {
    nm <- names(scores)
    if ("ctrl_norm_score" %in% nm) {
        m <- as.matrix(scores[["ctrl_norm_score"]])
    } else {
        cols <- grep("^ctrl_norm_score_", nm, value = TRUE)
        if (length(cols) == 0L) {
            stop("'scores' has no normalized control scores; use ",
                "scdrsScoreCells(..., returnCtrlNormScore = TRUE)")
        }
        m <- as.matrix(as.data.frame(scores[, cols, drop = FALSE]))
    }
    storage.mode(m) <- "double"
    unname(m)
}

## Check a score table and return the sorted cells shared with `sce`, as
## Python's sorted(set(adata.obs_names) & set(df_full_score.index)).
.scoredCells <- function(sce, scores) {
    if (!is(sce, "SingleCellExperiment")) {
        stop("'sce' must be a SingleCellExperiment")
    }
    if (!(is.data.frame(scores) || is(scores, "DataFrame")) ||
        is.null(rownames(scores)) || !"norm_score" %in% names(scores)) {
        stop("'scores' must be a score table with cell names as row ",
            "names, e.g. from scdrsScoreCells()")
    }
    cells <- sort(intersect(colnames(sce), rownames(scores)),
        method = "radix")
    if (length(cells) == 0L) {
        stop("no cells in common between 'sce' and 'scores'")
    }
    cells
}

## The cell-cell graph (n x n sparse matrix, dimnames = cell names): a
## colPair() of `sce` by name, or a matrix given directly.
.getGraph <- function(sce, graph) {
    if (is.character(graph) && length(graph) == 1L) {
        if (!graph %in% SingleCellExperiment::colPairNames(sce)) {
            stop("graph '", graph, "' not found in colPair(sce); ",
                "run scdrsBuildGraph() first or give a graph")
        }
        g <- SingleCellExperiment::colPair(sce, graph, asSparse = TRUE)
        dimnames(g) <- list(colnames(sce), colnames(sce))
    } else {
        g <- as(as(as(graph, "CsparseMatrix"), "generalMatrix"), "dMatrix")
        if (!identical(dim(g), rep(ncol(sce), 2L))) {
            stop("'graph' must be a cells x cells matrix")
        }
        if (!is.null(rownames(g))) {
            g <- g[colnames(sce), colnames(sce)]
        }
        dimnames(g) <- list(colnames(sce), colnames(sce))
    }
    g
}

## Sorted distinct values of a grouping variable, as Python's
## sorted(set(x)); missing values are dropped.
.groupLevels <- function(x) {
    x <- x[!is.na(x)]
    if (is.numeric(x)) {
        as.character(sort(unique(x)))
    } else {
        sort(unique(as.character(x)), method = "radix")
    }
}

## Quantile of each column with numpy's default (linear) method, i.e.
## R's type 7.
.colQuantile <- function(m, prob) {
    s <- apply(m, 2L, sort)
    s <- matrix(s, nrow = nrow(m))
    h <- (nrow(m) - 1) * prob
    lo <- floor(h)
    hi <- min(lo + 1, nrow(m) - 1)
    s[lo + 1, ] + (h - lo) * (s[hi + 1, ] - s[lo + 1, ])
}

#' Geary's C statistic
#'
#' Computes Geary's C, a measure of spatial autocorrelation, of one or
#' several cell-level variables on a cell-cell graph. Small values mean
#' that neighbouring cells have similar values. This is the R port of
#' `scdrs.method.gearys_c()`.
#'
#' @param graph A cells x cells matrix of edge weights (sparse or dense),
#'     e.g. from [scdrsBuildGraph()].
#' @param x A numeric vector with one value per cell, or a matrix with one
#'     row per cell (one statistic per column).
#'
#' @return A numeric vector with one value per column of `x`.
#'
#' @details
#' \deqn{C = \frac{(N - 1)\sum_{i,j} w_{ij} (x_i - x_j)^2}{2W \sum_i (x_i -
#' \bar{x})^2},}{C = (N - 1) sum_ij w_ij (x_i - x_j)^2 / (2 W sum_i (x_i -
#' mean(x))^2),}
#' with \eqn{W = \sum_{ij} w_{ij}}. It is computed for all columns at once
#' from \eqn{\sum_{ij} w_{ij}(x_i - x_j)^2 = \sum_i r_i x_i^2 + \sum_j k_j
#' x_j^2 - 2x^\top W x} (row sums \eqn{r}, column sums \eqn{k}), on the
#' centered values.
#'
#' @examples
#' g <- Matrix::sparseMatrix(i = c(1, 2, 2, 3), j = c(2, 1, 3, 2),
#'     x = 1, dims = c(3, 3))
#' scdrsComputeGearysC(g, cbind(smooth = c(1, 2, 3), rough = c(1, 3, 1)))
#'
#' @export
scdrsComputeGearysC <- function(graph, x) {
    w <- as(as(as(graph, "CsparseMatrix"), "generalMatrix"), "dMatrix")
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    if (nrow(w) != ncol(w) || nrow(w) != nrow(x)) {
        stop("'graph' must be square with one row per row of 'x'")
    }
    n <- nrow(w)
    xc <- sweep(x, 2L, colMeans(x))
    xc2 <- xc^2
    num <- colSums(Matrix::rowSums(w) * xc2) +
        colSums(Matrix::colSums(w) * xc2) -
        2 * colSums(xc * as.matrix(w %*% xc))
    unname((n - 1) * num / (2 * sum(w) * colSums(xc2)))
}

#' Test the heterogeneity of disease scores within cell groups
#'
#' For each group of cells, computes Geary's C of the normalized disease
#' score on the group's part of the cell-cell graph and compares it with
#' null statistics, to test whether disease scores are more similar among
#' neighbouring cells than expected (heterogeneity within the group). This
#' is the R port of `scdrs.method.test_gearysc()`.
#'
#' @param graph A cells x cells matrix of edge weights, with the cells in
#'     the order of the rows of `scores` (or with cell names as dimnames).
#' @param scores A score table with normalized control scores, from
#'     [scdrsScoreCells()] with `returnCtrlNormScore = TRUE` (or read from
#'     a Python `.full_score.gz` file with [scdrsReadScores()]).
#' @param groups A vector with the group of each cell (rows of `scores`).
#' @param method How the null statistics are made:
#'     `"control_distribution_match"` (default: the disease scores of the
#'     group, rearranged in the order of each control score),
#'     `"permutation"` (random permutations of the disease scores; uses
#'     random numbers) or `"control"` (the control scores themselves).
#'
#' @return A [S4Vectors::DataFrame] with one row per group (sorted):
#'     `pval` (\eqn{(1 + \#\{k: C > C_k\}) / (1 + K)}; small values mean
#'     heterogeneity), `trait` (Geary's C of the disease score),
#'     `ctrl_mean`, `ctrl_std` (mean and sample SD of the null statistics)
#'     and `zsc` (\eqn{-(C - \mathrm{ctrl\_mean}) /
#'     \mathrm{ctrl\_std}}{-(C - ctrl_mean) / ctrl_std}).
#'
#' @seealso [scdrsTestGroups()], which calls this function.
#'
#' @examples
#' sce <- scdrsPreprocess(scdrsBuildGraph(scdrsLoadToyData()))
#' gs <- scdrsReadGeneSets(system.file("extdata", "toydata_mouse.gs",
#'     package = "scdRs"))
#' set.seed(1)
#' res <- scdrsScoreCells(sce, gs[[1]], nCtrl = 100,
#'     returnCtrlNormScore = TRUE)
#' g <- SingleCellExperiment::colPair(sce, "connectivities", asSparse = TRUE)
#' scdrsTestGearysC(g, res, sce$cell_type)
#'
#' @export
scdrsTestGearysC <- function(graph, scores, groups,
    method = c("control_distribution_match", "permutation", "control")) {
    method <- match.arg(method)
    z <- as.numeric(scores$norm_score)
    Z <- .ctrlNormScores(scores)
    g <- as(as(as(graph, "CsparseMatrix"), "generalMatrix"), "dMatrix")
    if (!identical(dim(g), rep(length(z), 2L))) {
        stop("'graph' must be a cells x cells matrix matching 'scores'")
    }
    if (!is.null(rownames(g)) && !is.null(rownames(scores))) {
        g <- g[rownames(scores), rownames(scores)]
    }
    if (length(groups) != length(z)) {
        stop("'groups' must have one value per row of 'scores'")
    }
    nNull <- ncol(Z)
    lev <- .groupLevels(groups)
    groups <- as.character(groups)
    res <- vapply(lev, function(grp) {
        idx <- which(groups == grp)
        zg <- z[idx]
        Zg <- Z[idx, , drop = FALSE]
        nulls <- switch(method,
            control_distribution_match = {
                r <- matrix(apply(Zg, 2L, rank, ties.method = "first"),
                    nrow = length(idx))
                matrix(sort(zg)[r], nrow = length(idx))
            },
            permutation = vapply(seq_len(nNull),
                function(k) zg[sample.int(length(zg))], numeric(length(zg))),
            control = Zg)
        cc <- scdrsComputeGearysC(g[idx, idx, drop = FALSE],
            cbind(zg, matrix(nulls, nrow = length(idx))))
        trait <- cc[1L]
        null <- cc[-1L]
        pval <- if (is.na(trait)) {
            NA_real_
        } else {
            (sum(trait > null, na.rm = TRUE) + 1) / (nNull + 1)
        }
        m <- mean(null, na.rm = TRUE)
        s <- stats::sd(null, na.rm = TRUE)
        c(pval = pval, trait = trait, ctrl_mean = m, ctrl_std = s,
            zsc = -(trait - m) / s)
    }, numeric(5))
    res <- matrix(res, nrow = 5L, dimnames = list(c("pval", "trait",
        "ctrl_mean", "ctrl_std", "zsc"), lev))
    S4Vectors::DataFrame(t(res), check.names = FALSE)
}

#' Group-level scDRS analysis
#'
#' For each annotation in `groupCols` (e.g. cell types) and each group of
#' cells, counts the significant cells and tests whether the group as a
#' whole is associated with the disease (association) and whether
#' disease scores vary among the group's cells beyond chance
#' (heterogeneity). This is the R port of
#' `scdrs.method.downstream_group_analysis()`.
#'
#' @param sce A [SingleCellExperiment::SingleCellExperiment] with the
#'     annotations in `colData` and a cell-cell graph (see `graph`).
#' @param scores A score table with normalized control scores: the output
#'     of [scdrsScoreCells()] with `returnCtrlNormScore = TRUE`, one
#'     element of [scdrsComputeScore()], or a table read with
#'     [scdrsReadScores()].
#' @param groupCols Names of `colData(sce)` columns defining the groups.
#' @param fdrThresholds FDR levels at which significant cells are counted.
#' @param graph Name of a `colPair(sce)` graph (default
#'     `"connectivities"`, made by [scdrsBuildGraph()]), or a cells x cells
#'     matrix.
#'
#' @return A named [S4Vectors::SimpleList] with one [S4Vectors::DataFrame]
#'     per element of `groupCols`, with one row per group (sorted) and
#'     columns:
#' \describe{
#'   \item{`n_cell`, `n_ctrl`}{Numbers of cells and of control sets.}
#'   \item{`assoc_mcp`, `assoc_mcz`}{Association: Monte Carlo p-value and
#'     z-score of the 95% quantile of the group's disease scores against
#'     the same quantile of each control score.}
#'   \item{`hetero_mcp`, `hetero_mcz`}{Heterogeneity: p-value and z-score
#'     of Geary's C (see [scdrsTestGearysC()]).}
#'   \item{`n_fdr_<t>`}{Number of cells with FDR (Benjamini-Hochberg, over
#'     all cells) below `t`, for each of `fdrThresholds`.}
#' }
#'
#' @seealso [scdrsTestCorrelation()], [scdrsRankGenes()],
#'     [scdrsPerformDownstream()].
#' @references Zhang MJ, Hou K, et al. (2022). Polygenic enrichment
#'     distinguishes disease associations of individual cells in single-cell
#'     RNA-seq data. *Nature Genetics* 54, 1572-1580.
#'
#' @examples
#' sce <- scdrsBuildGraph(scdrsLoadToyData())
#' sce <- scdrsPreprocess(sce, cov = "covariate")
#' gs <- scdrsReadGeneSets(system.file("extdata", "toydata_mouse.gs",
#'     package = "scdRs"))
#' set.seed(1)
#' res <- scdrsScoreCells(sce, gs[[1]], nCtrl = 100,
#'     returnCtrlNormScore = TRUE)
#' scdrsTestGroups(sce, res, "cell_type")[["cell_type"]]
#'
#' @export
scdrsTestGroups <- function(sce, scores, groupCols,
    fdrThresholds = c(0.05, 0.1, 0.2), graph = "connectivities") {
    cells <- .scoredCells(sce, scores)
    cd <- SummarizedExperiment::colData(sce)
    miss <- setdiff(groupCols, names(cd))
    if (length(miss) > 0L) {
        stop("'groupCols' not in colData(sce): ", paste(miss, collapse = ", "))
    }
    rows <- match(cells, rownames(scores))
    sub <- scores[rows, , drop = FALSE]
    z <- as.numeric(sub$norm_score)
    Z <- .ctrlNormScores(sub)
    nCtrl <- ncol(Z)
    fdr <- stats::p.adjust(as.numeric(sub$pval), "BH")
    g <- .getGraph(sce, graph)[cells, cells]
    out <- lapply(groupCols, function(col) {
        lev <- .groupLevels(cd[[col]])
        grp <- as.character(cd[cells, col])
        hetero <- scdrsTestGearysC(g, sub, grp)
        res <- vapply(lev, function(l) {
            idx <- which(grp == l)
            nFdr <- vapply(fdrThresholds, function(t) sum(fdr[idx] < t),
                numeric(1))
            if (length(idx) == 0L) {
                warning("group '", l, "' of '", col, "' has no scored cells")
                return(c(0, nCtrl, NA, NA, NA, NA, nFdr))
            }
            q <- .colQuantile(matrix(z[idx]), 0.95)
            qCtrl <- .colQuantile(Z[idx, , drop = FALSE], 0.95)
            mcp <- (sum(qCtrl >= q) + 1) / (nCtrl + 1)
            mcz <- (q - mean(qCtrl)) / .sd0(qCtrl)
            h <- hetero[l, ]
            c(length(idx), nCtrl, mcp, mcz, h$pval, h$zsc, nFdr)
        }, numeric(6 + length(fdrThresholds)))
        res <- matrix(res, ncol = length(lev), dimnames = list(c("n_cell",
            "n_ctrl", "assoc_mcp", "assoc_mcz", "hetero_mcp", "hetero_mcz",
            paste0("n_fdr_", fdrThresholds)), lev))
        S4Vectors::DataFrame(t(res), check.names = FALSE)
    })
    S4Vectors::SimpleList(stats::setNames(out, groupCols))
}

#' Correlation of disease scores with cell-level variables
#'
#' For each continuous cell-level variable (e.g. a T cell effector
#' gradient), tests the Pearson correlation between the variable and the
#' normalized disease score against the correlations with the control
#' scores. This is the R port of
#' `scdrs.method.downstream_corr_analysis()`.
#'
#' @inheritParams scdrsTestGroups
#' @param varCols Names of numeric `colData(sce)` columns.
#'
#' @return A [S4Vectors::DataFrame] with one row per variable and columns
#'     `n_ctrl`, `corr_mcp` (Monte Carlo p-value) and `corr_mcz` (z-score
#'     against the control correlations).
#'
#' @examples
#' sce <- scdrsPreprocess(scdrsLoadToyData(), cov = "covariate")
#' gs <- scdrsReadGeneSets(system.file("extdata", "toydata_mouse.gs",
#'     package = "scdRs"))
#' set.seed(1)
#' res <- scdrsScoreCells(sce, gs[[1]], nCtrl = 100,
#'     returnCtrlNormScore = TRUE)
#' scdrsTestCorrelation(sce, res, c("causal_variable",
#'     "non_causal_variable"))
#'
#' @export
scdrsTestCorrelation <- function(sce, scores, varCols) {
    cells <- .scoredCells(sce, scores)
    cd <- SummarizedExperiment::colData(sce)
    miss <- setdiff(varCols, names(cd))
    if (length(miss) > 0L) {
        stop("'varCols' not in colData(sce): ", paste(miss, collapse = ", "))
    }
    sub <- scores[match(cells, rownames(scores)), , drop = FALSE]
    z <- as.numeric(sub$norm_score)
    Z <- .ctrlNormScores(sub)
    res <- vapply(varCols, function(col) {
        v <- as.numeric(cd[cells, col])
        r <- stats::cor(v, z)
        rCtrl <- as.vector(stats::cor(v, Z))
        c(n_ctrl = ncol(Z), corr_mcp = (sum(rCtrl >= r) + 1) / (ncol(Z) + 1),
            corr_mcz = (r - mean(rCtrl)) / .sd0(rCtrl))
    }, numeric(3))
    S4Vectors::DataFrame(t(matrix(res, nrow = 3L, dimnames = list(
        c("n_ctrl", "corr_mcp", "corr_mcz"), varCols))), check.names = FALSE)
}

## Port of _pearson_corr() / _pearson_corr_sparse(): Pearson correlation
## of every column of x (n x p) with every column of y (n x q); SDs are
## population SDs clipped below at 1e-8. If either input is sparse, means
## and variances use E[x^2] - E[x]^2 as Python's sparse version does.
## Returns a p x q matrix, or a vector if p or q is 1.
.pearsonCor <- function(x, y) {
    if (is.null(dim(x))) x <- matrix(x)
    if (is.null(dim(y))) y <- matrix(y)
    n <- nrow(x)
    if (is(x, "sparseMatrix") || is(y, "sparseMatrix")) {
        mx <- Matrix::colMeans(x)
        my <- Matrix::colMeans(y)
        sx <- pmax(sqrt(Matrix::colMeans(x^2) - mx^2), 1e-8)
        sy <- pmax(sqrt(Matrix::colMeans(y^2) - my^2), 1e-8)
        r <- (as.matrix(Matrix::crossprod(x, y)) / n - outer(mx, my)) /
            outer(sx, sy)
    } else {
        sc <- function(m) {
            mu <- colMeans(m)
            m <- sweep(m, 2L, mu)
            sweep(m, 2L, pmax(sqrt(colMeans(m^2)), 1e-8), "/")
        }
        r <- crossprod(sc(x), sc(y)) / n
    }
    r <- unname(as.matrix(r))
    if (ncol(x) == 1L || ncol(y) == 1L) as.vector(r) else r
}

#' Rank genes by their correlation with the disease score
#'
#' Computes the Pearson correlation across cells between each gene's
#' expression and the normalized disease score, to prioritize genes. This
#' is the R port of `scdrs.method.downstream_gene_analysis()`.
#'
#' @inheritParams scdrsTestGroups
#' @param assay Expression assay to use. As in Python, this is the
#'     log-normalized data, not the covariate-corrected data.
#'
#' @return A [S4Vectors::DataFrame] with one row per gene (row names),
#'     sorted by decreasing correlation, with columns `CORR` and `RANK`
#'     (1 = highest correlation; Python's `RANK` starts at 0).
#'
#' @examples
#' sce <- scdrsPreprocess(scdrsLoadToyData(), cov = "covariate")
#' gs <- scdrsReadGeneSets(system.file("extdata", "toydata_mouse.gs",
#'     package = "scdRs"))
#' set.seed(1)
#' res <- scdrsScoreCells(sce, gs[[1]], nCtrl = 100)
#' head(scdrsRankGenes(sce, res))
#'
#' @export
scdrsRankGenes <- function(sce, scores, assay = "logcounts") {
    cells <- .scoredCells(sce, scores)
    x <- .getAssay(sce, assay)[, cells, drop = FALSE]
    z <- as.numeric(scores[cells, "norm_score"])
    r <- .pearsonCor(Matrix::t(x), z)
    o <- order(-r, method = "radix")
    S4Vectors::DataFrame(CORR = r[o], RANK = seq_along(o),
        row.names = rownames(x)[o])
}

#' Run the scDRS downstream analyses
#'
#' Runs the group-level analysis ([scdrsTestGroups()]), the correlation
#' analysis ([scdrsTestCorrelation()]) and the gene analysis
#' ([scdrsRankGenes()]) for one or several traits. This is the R port of
#' the `scdrs perform-downstream` command.
#'
#' @inheritParams scdrsTestGroups
#' @param scores One score table, or a named list of them (e.g. from
#'     [scdrsComputeScore()] or [scdrsReadScores()]), with normalized
#'     control scores.
#' @param groupCols,varCols `NULL`, or `colData(sce)` columns for the group
#'     analysis and the correlation analysis.
#' @param geneAnalysis Run the gene analysis?
#' @param assay Assay used by the gene analysis.
#'
#' @return For one score table, a list with elements `group`, `corr` and
#'     `gene` (`NULL` for analyses not run); for a list of score tables, a
#'     named list of such lists.
#'
#' @examples
#' sce <- scdrsBuildGraph(scdrsLoadToyData())
#' gs <- scdrsReadGeneSets(system.file("extdata", "toydata_mouse.gs",
#'     package = "scdRs"))
#' set.seed(1)
#' res <- scdrsComputeScore(sce, gs, cov = "covariate", nCtrl = 100)
#' out <- scdrsPerformDownstream(sce, res, groupCols = "cell_type",
#'     varCols = "causal_variable")
#' out$toydata_gs_mouse$group$cell_type
#'
#' @export
scdrsPerformDownstream <- function(sce, scores, groupCols = NULL,
    varCols = NULL, geneAnalysis = FALSE, graph = "connectivities",
    assay = "logcounts") {
    if (is.null(groupCols) && is.null(varCols) && !geneAnalysis) {
        stop("expect at least one of 'groupCols', 'varCols' and ",
            "'geneAnalysis = TRUE'")
    }
    one <- function(s) {
        list(group = if (!is.null(groupCols)) {
                scdrsTestGroups(sce, s, groupCols, graph = graph)
            },
            corr = if (!is.null(varCols)) {
                scdrsTestCorrelation(sce, s, varCols)
            },
            gene = if (geneAnalysis) scdrsRankGenes(sce, s, assay))
    }
    if (is.data.frame(scores) || is(scores, "DataFrame")) {
        return(one(scores))
    }
    if (is.null(names(scores))) {
        stop("'scores' must be a score table or a named list of them")
    }
    lapply(scores, one)
}
