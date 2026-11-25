## Internal helpers for preprocessing (port of scdrs/pp.py).
##
## Orientation: matrices are genes x cells (SingleCellExperiment), the
## transpose of AnnData. "gene" statistics are computed across cells (Python
## axis = 0) and "cell" statistics across genes (Python axis = 1).

## Fetch an assay as a dgCMatrix (sparse input) or a double matrix.
.getAssay <- function(sce, assay) {
    if (!is(sce, "SingleCellExperiment")) {
        stop("'sce' must be a SingleCellExperiment")
    }
    if (!is.character(assay) || length(assay) != 1L ||
        !assay %in% SummarizedExperiment::assayNames(sce)) {
        stop("assay '", paste(assay, collapse = ", "), "' not found in 'sce'")
    }
    if (is.null(rownames(sce)) || anyDuplicated(rownames(sce))) {
        stop("'sce' must have unique gene names (rownames)")
    }
    if (is.null(colnames(sce)) || anyDuplicated(colnames(sce))) {
        stop("'sce' must have unique cell names (colnames)")
    }
    x <- SummarizedExperiment::assay(sce, assay, withDimnames = TRUE)
    if (is(x, "sparseMatrix")) {
        x <- as(as(as(x, "CsparseMatrix"), "generalMatrix"), "dMatrix")
    } else {
        x <- as.matrix(x)
        storage.mode(x) <- "double"
    }
    x
}

## Port of pp.reg_out(): residuals of the least-squares regression of y on
## x, solved with the normal equations as in Python (not with qr()).
.regOut <- function(y, x) {
    y <- as.matrix(y)
    x <- as.matrix(x)
    n <- nrow(y)
    coef <- solve(crossprod(x) / n, crossprod(x, y) / n)
    y - x %*% coef
}

## Port of pp.category2dummy(): non-numeric columns become 0/1 dummy
## columns named "<col>_<level>", dropping the first level, as
## pd.get_dummies(drop_first = TRUE). Levels of character columns are
## sorted by code point (radix sort), as Python's sorted(); factors keep
## their level order (like a pandas Categorical). Missing values stay NA in
## every dummy column. Logical columns count as numeric (0/1), as in pandas.
.categoryToDummy <- function(df, cols = NULL, verbose = FALSE) {
    df <- as.data.frame(df, optional = TRUE)
    isNum <- vapply(df, function(v) is.numeric(v) || is.logical(v),
        logical(1))
    if (is.null(cols)) {
        cols <- names(df)[!isNum]
    }
    if (!all(cols %in% names(df))) {
        stop("'cols' must be a subset of the column names")
    }
    out <- df[, setdiff(names(df), cols), drop = FALSE]
    for (j in seq_along(out)) {
        if (is.logical(out[[j]])) {
            out[[j]] <- as.numeric(out[[j]])
        }
    }
    added <- character(0)
    for (col in cols) {
        v <- df[[col]]
        lev <- if (is.factor(v)) {
            levels(v)
        } else {
            sort(unique(as.character(v[!is.na(v)])), method = "radix")
        }
        v <- as.character(v)
        for (l in lev[-1]) {
            nm <- paste0(col, "_", l)
            out[[nm]] <- ifelse(is.na(v), NA_real_, as.numeric(v == l))
            added <- c(added, nm)
        }
    }
    if (verbose && length(added) > 0L) {
        message("Detected categorical columns: ", paste(cols, collapse = ","),
            ". Added dummy columns: ", paste(added, collapse = ","), ".")
    }
    out
}

## Build the covariate matrix (cells x covariates) as pp.preprocess() does:
## match to cells, dummy-code, fill NA with column means, and add a constant
## column "SCDRS_CONST" unless the covariates already span the constant.
.covariateMatrix <- function(sce, cov) {
    cells <- colnames(sce)
    if (is.character(cov)) {
        miss <- setdiff(cov, names(SummarizedExperiment::colData(sce)))
        if (length(miss) > 0L) {
            stop("covariate column(s) not in colData(sce): ",
                paste(miss, collapse = ", "))
        }
        df <- as.data.frame(SummarizedExperiment::colData(sce)[, cov,
            drop = FALSE], optional = TRUE)
    } else if (is.data.frame(cov) || is(cov, "DataFrame")) {
        df <- as.data.frame(cov, optional = TRUE)
        if (is.null(rownames(cov))) {
            stop("'cov' must have cell names as row names")
        }
        if (!(sum(unique(rownames(df)) %in% cells) > 0.75 * length(cells))) {
            stop("'cov' does not match the cells in 'sce' (fewer than 75% ",
                "of the cells found)")
        }
        df <- df[match(cells, rownames(df)), , drop = FALSE]
    } else {
        stop("'cov' must be a character vector of colData columns or a ",
            "data.frame with cell names as row names")
    }
    m <- as.matrix(.categoryToDummy(df, verbose = TRUE))
    storage.mode(m) <- "double"
    for (j in seq_len(ncol(m))) {
        na <- is.na(m[, j])
        if (any(na)) {
            m[na, j] <- mean(m[!na, j])
        }
    }
    if (anyNA(m)) {
        stop("covariate(s) with no non-missing values: ",
            paste(colnames(m)[colSums(is.na(m)) > 0], collapse = ", "))
    }
    if (ncol(m) == 0L ||
        mean(.regOut(rep(1, nrow(m)), m)^2) > 0.01) {
        m <- cbind(m, SCDRS_CONST = 1)
    }
    rownames(m) <- cells
    m
}

## Cell weights for adjusting for cell-group proportions (adj_prop): each
## cell is weighted by n / (size of its group); group sizes are clipped
## below at floor(1% of the largest group); weights have mean 1.
.cellWeights <- function(sce, adjProp) {
    cd <- SummarizedExperiment::colData(sce)
    if (!is.character(adjProp) || length(adjProp) != 1L ||
        !adjProp %in% names(cd)) {
        stop("'adjProp' = ", paste(adjProp, collapse = ", "),
            " not in colData(sce)")
    }
    grp <- as.character(cd[[adjProp]])
    if (anyNA(grp)) {
        stop("'adjProp' column '", adjProp, "' has missing values")
    }
    n <- length(grp)
    if (!(length(unique(grp)) < 0.1 * n)) {
        stop("On average <10 cells per group, maybe '", adjProp,
            "' is not categorical?")
    }
    tab <- table(grp)
    size <- pmax(as.numeric(tab), floor(0.01 * max(tab)))
    names(size) <- names(tab)
    w <- n / size[grp]
    unname(w / mean(w))
}

## Port of pp._get_mean_var(): mean and variance per gene (across cells) or
## per cell (across genes). Unweighted dense input uses a two-pass variance
## (np.var); sparse or weighted input uses E[x^2] - E[x]^2, as in Python.
.meanVar <- function(x, margin = c("gene", "cell"), weights = NULL) {
    margin <- match.arg(margin)
    byGene <- margin == "gene"
    if (is.null(weights)) {
        m <- if (byGene) Matrix::rowMeans(x) else Matrix::colMeans(x)
        if (is(x, "sparseMatrix")) {
            m2 <- if (byGene) Matrix::rowMeans(x^2) else Matrix::colMeans(x^2)
            v <- m2 - m^2
        } else if (byGene) {
            v <- rowMeans((x - m)^2)
        } else {
            v <- colMeans((x - rep(m, each = nrow(x)))^2)
        }
    } else {
        w <- as.numeric(weights)
        if (length(w) != (if (byGene) ncol(x) else nrow(x))) {
            stop("'weights' has the wrong length")
        }
        if (byGene) {
            m <- as.vector(x %*% w) / sum(w)
            v <- as.vector(x^2 %*% w) / sum(w) - m^2
        } else {
            m <- as.vector(Matrix::crossprod(x, w)) / sum(w)
            v <- as.vector(Matrix::crossprod(x^2, w)) / sum(w) - m^2
        }
    }
    list(mean = unname(as.vector(m)), var = unname(as.vector(v)))
}

## Port of pp._get_mean_var_implicit_cov_corr(): statistics of
## fun(x + covBeta %*% t(covMat) + geneMean) without building the whole
## corrected matrix; works on chunks of genes ("gene") or cells ("cell").
## covMat: cells x k; covBeta: genes x k (already negated, = COV_BETA).
## Python's chunk size n // n_chunk can be 0 (infinite loop); R uses at
## least 1. Results do not depend on the chunking.
.meanVarImplicit <- function(x, covMat, covBeta, geneMean,
    margin = c("gene", "cell"), weights = NULL, fun = identity,
    nChunk = 20L) {
    margin <- match.arg(margin)
    byGene <- margin == "gene"
    n <- if (byGene) nrow(x) else ncol(x)
    size <- max(1L, n %/% as.integer(nChunk))
    tCov <- t(covMat)
    m <- v <- numeric(n)
    for (start in seq(1L, n, by = size)) {
        idx <- seq.int(start, min(start + size - 1L, n))
        if (byGene) {
            chunk <- as.matrix(x[idx, , drop = FALSE]) +
                covBeta[idx, , drop = FALSE] %*% tCov + geneMean[idx]
        } else {
            chunk <- as.matrix(x[, idx, drop = FALSE]) +
                covBeta %*% tCov[, idx, drop = FALSE] + geneMean
        }
        r <- .meanVar(fun(chunk), margin, weights)
        m[idx] <- r$mean
        v[idx] <- r$var
    }
    list(mean = m, var = v)
}

## expm1() that keeps a sparse matrix sparse (expm1(0) = 0).
.expm1 <- function(x) {
    if (is(x, "sparseMatrix")) {
        x@x <- expm1(x@x)
        x
    } else {
        expm1(x)
    }
}

## Quantiles exactly as numpy's default ("linear") method, which pandas
## uses; probabilities as np.linspace(0, 1, q + 1) and the * 100 / 100 round
## trip of pandas -> np.percentile, so cut points match to the last bit.
.npQuantile <- function(x, q) {
    x <- sort(x, method = "radix")
    n <- length(x)
    p <- c(seq(0, q - 1) * (1 / q), 1)
    p <- (p * 100) / 100
    vi <- (n - 1) * p
    lo <- floor(vi)
    g <- vi - lo
    above <- vi >= n - 1
    lo[above] <- n - 1
    hi <- pmin(lo + 1, n - 1)
    a <- x[lo + 1]
    b <- x[hi + 1]
    d <- b - a
    out <- ifelse(g >= 0.5, b - d * (1 - g), a + d * g)
    out[above] <- x[n]
    out
}

## Port of pd.qcut(x, q, labels = False, duplicates = "drop"): 0-based bin
## labels; bins are right-closed, the lowest edge is included.
.qcut <- function(x, q) {
    edges <- unique(.npQuantile(x, q))
    if (length(edges) < 2L) {
        return(integer(length(x)))
    }
    findInterval(x, edges, left.open = TRUE, all.inside = TRUE) - 1L
}

## Technical variance, the Seurat v3 recipe used by scDRS: LOESS of
## log10(ct_var) on log10(ct_mean) (span 0.3, degree 2); genes left out of
## the fit get ct_var_tech = 1 as in Python.
.techVar <- function(ctMean, ctVar) {
    use <- ctVar > 0
    if (any(ctMean <= 0)) {
        warning(sum(ctMean < 0), " genes with ct_mean<0")
        use <- use & ctMean > 0
    }
    lx <- log10(ctMean[use])
    ly <- log10(ctVar[use])
    fit <- stats::loess(ly ~ lx, span = 0.3, degree = 2,
        family = "gaussian",
        control = stats::loess.control(surface = "interpolate",
            statistics = "approximate", cell = 0.2))
    est <- numeric(length(ctMean))
    est[use] <- stats::fitted(fit)
    10^est
}

## Mean-variance bins: qcut of the mean, then qcut of the variance within
## each mean bin; labels "<meanBin>.<varBin>". If either bin count is above
## floor(sqrt(nGene / 10)), both are set to that value, as in Python.
.meanVarBins <- function(mean, var, nMeanBin, nVarBin) {
    nBinMax <- floor(sqrt(length(mean) / 10))
    if (nMeanBin > nBinMax || nVarBin > nBinMax) {
        message("Too few genes for ", nMeanBin, "*", nVarBin,
            " bins, setting nMeanBin = nVarBin = ", nBinMax)
        nMeanBin <- nVarBin <- nBinMax
    }
    meanBin <- .qcut(mean, nMeanBin)
    out <- character(length(mean))
    for (b in sort(unique(meanBin))) {
        sel <- meanBin == b
        out[sel] <- sprintf("%d.%d", b, .qcut(var[sel], nVarBin))
    }
    out
}
