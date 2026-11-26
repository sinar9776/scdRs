## Provenance: R port of the Python scDRS 1.0.4
## (https://github.com/martinjzhang/scDRS; MIT License, Copyright (c)
## 2021 Martin Jinye Zhang and Kangcheng Hou; see inst/COPYRIGHTS).
## Assisted-by: Claude (Anthropic), reviewed by the maintainer.

#' Convert between one-sided p-values and z-scores
#'
#' `scdrsZToP()` returns \eqn{p = \Phi(-z)} and `scdrsPToZ()` returns
#' \eqn{z = -\Phi^{-1}(p)}, accurate far in the tails. These are the R
#' ports of `scdrs.util.zsc2pval()` and `scdrs.util.pval2zsc()`.
#'
#' @param z Numeric vector of z-scores.
#' @param p Numeric vector of p-values.
#' @return A numeric vector of the same length.
#'
#' @examples
#' scdrsZToP(c(0, 1.96, 10))
#' scdrsPToZ(c(0.5, 0.025, 1e-30))
#'
#' @export
scdrsZToP <- function(z) {
    stats::pnorm(-z)
}

#' @rdname scdrsZToP
#' @export
scdrsPToZ <- function(p) {
    -stats::qnorm(p)
}

#' Meta-analysis of effect sizes
#'
#' Fixed-effect or DerSimonian-Laird random-effects meta-analysis, e.g. to
#' combine scDRS results across data sets. This is the R port of
#' `scdrs.util.meta_analysis()`.
#'
#' @param effects Numeric vector of effect sizes.
#' @param se Numeric vector of their standard errors.
#' @param method `"random"` (default) or `"fixed"`.
#' @param weights Optional weights to use instead of the inverse
#'     variances.
#'
#' @return A named numeric vector: `summary` (the combined effect) and
#'     `se` (its standard error).
#'
#' @details With \eqn{v_i = 1/\mathrm{se}_i^2}, the between-study variance
#' is \eqn{\tau^2 = \max(0, (Q - (m-1)) / (\sum v - \sum v^2 / \sum v))},
#' where \eqn{Q} is Cochran's Q. Random-effects weights are
#' \eqn{1 / (\mathrm{se}^2 + \tau^2)}.
#'
#' @examples
#' scdrsMetaAnalyze(c(0.2, 0.5, -0.1), c(0.1, 0.2, 0.15))
#' scdrsMetaAnalyze(c(0.2, 0.5, -0.1), c(0.1, 0.2, 0.15), method = "fixed")
#'
#' @export
scdrsMetaAnalyze <- function(effects, se, method = c("random", "fixed"),
    weights = NULL) {
    method <- match.arg(method)
    if (length(effects) != length(se)) {
        stop("'effects' and 'se' must have the same length")
    }
    v <- se^2
    vw <- 1 / v
    fixed <- sum(vw * effects) / sum(vw)
    q <- sum((effects - fixed)^2 / v)
    tau2 <- max(0, (q - (length(effects) - 1)) /
        (sum(vw) - sum(vw^2) / sum(vw)))
    total <- if (method == "fixed") v else v + tau2
    wt <- if (is.null(weights)) 1 / total else weights
    summ <- sum(wt * effects) / sum(wt)
    c(summary = summ, se = sqrt(sum(wt^2 * total) / sum(wt)^2))
}

#' Test the overlap of two gene sets
#'
#' Fisher's exact test for the overlap of two gene sets within a
#' background, with the sample odds ratio and its 95% confidence interval.
#' This is the R port of `scdrs.util.test_overlap()`.
#'
#' @param set1,set2 Character vectors of genes.
#' @param background Character vector of all genes considered.
#'
#' @return A named numeric vector: `pval` (two-sided Fisher's exact test),
#'     `odds_ratio` (sample odds ratio), `or_ub` and `or_lb` (upper and
#'     lower 95% Wald limits on the log scale; both 0 if a cell of the 2 x
#'     2 table is 0, as in Python).
#'
#' @details The table counts background genes in neither set, in only one
#' of them, and in both. Unlike [stats::fisher.test()], the odds ratio is
#' the sample odds ratio (as `scipy.stats.fisher_exact()`), not the
#' conditional maximum likelihood estimate.
#'
#' @examples
#' bg <- paste0("g", seq_len(100))
#' scdrsTestOverlap(bg[1:30], bg[21:45], bg)
#'
#' @export
scdrsTestOverlap <- function(set1, set2, background) {
    set1 <- unique(set1)
    set2 <- unique(set2)
    n1 <- length(set1)
    n2 <- length(set2)
    nOverlap <- length(intersect(set1, set2))
    nOther <- length(setdiff(unique(background), union(set1, set2)))
    tab <- matrix(c(nOther, n2 - nOverlap, n1 - nOverlap, nOverlap), 2L)
    pval <- stats::fisher.test(tab)$p.value
    odds <- if (n1 > nOverlap && n2 > nOverlap) {
        (nOther * nOverlap) / ((n1 - nOverlap) * (n2 - nOverlap))
    } else {
        Inf
    }
    if (nOverlap == 0 || nOther == 0 || n2 == nOverlap || n1 == nOverlap) {
        return(c(pval = pval, odds_ratio = odds, or_ub = 0, or_lb = 0))
    }
    se <- sqrt(1 / (n1 - nOverlap) + 1 / (n2 - nOverlap) + 1 / nOverlap +
        1 / nOther)
    c(pval = pval, odds_ratio = odds, or_ub = exp(log(odds) + 1.96 * se),
        or_lb = exp(log(odds) - 1.96 * se))
}
