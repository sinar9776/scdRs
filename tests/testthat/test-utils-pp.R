## Internal helpers of the preprocessing step, each compared with the
## Python function it ports, on small fixed inputs exported by
## make_reference_outputs.py (export_units). No random numbers involved.

test_that(".regOut matches pp.reg_out", {
    y <- as.matrix(readReference("units.reg_out.Y", rowNames = FALSE))
    x <- as.matrix(readReference("units.reg_out.X", rowNames = FALSE))
    ref <- as.matrix(readReference("units.reg_out.resid", rowNames = FALSE))
    res <- scdRs:::.regOut(y, x)
    expect_equal(unname(res), unname(ref), tolerance = 1e-10)
})

test_that(".categoryToDummy matches pp.category2dummy", {
    df <- readReference("units.category2dummy.input")
    ref <- readReference("units.category2dummy.output")
    res <- scdRs:::.categoryToDummy(df)
    ## Python sorts levels by code point: "C" < "a" < "b", so "C" is the
    ## dropped first level.
    expect_identical(names(res), names(ref))
    expect_equal(as.matrix(res), as.matrix(ref), ignore_attr = TRUE)
    ## a factor keeps its own level order
    res2 <- scdRs:::.categoryToDummy(data.frame(
        g = factor(c("x", "y", "z"), levels = c("z", "y", "x"))))
    expect_identical(names(res2), c("g_y", "g_x"))
})

test_that(".meanVar matches pp._get_mean_var (dense/sparse, weights)", {
    m <- as.matrix(readReference("units.mean_var.X", rowNames = FALSE))
    w <- readReference("units.mean_var.w", rowNames = FALSE)$w
    ref <- readReference("units.mean_var.output", rowNames = FALSE)
    x <- t(unname(m))  # Python is cells x genes; R is genes x cells
    inputs <- list(dense = x, sparse = Matrix::Matrix(x, sparse = TRUE))
    for (i in seq_len(nrow(unique(ref[, c("kind", "axis", "weights")])))) {
        key <- unique(ref[, c("kind", "axis", "weights")])[i, ]
        sub <- ref[ref$kind == key$kind & ref$axis == key$axis &
            ref$weights == key$weights, ]
        margin <- if (key$axis == 0) "gene" else "cell"
        wt <- if (key$weights == "w") w else NULL
        res <- scdRs:::.meanVar(inputs[[key$kind]], margin, wt)
        expect_equal(res$mean, sub$mean, tolerance = 1e-10)
        expect_equal(res$var, sub$var, tolerance = 1e-10)
    }
})

test_that(".qcut matches pd.qcut(labels=False, duplicates='drop')", {
    ref <- readReference("units.qcut", rowNames = FALSE)
    expect_identical(scdRs:::.qcut(ref$x, 7), as.integer(ref$bin))
})

test_that(".meanVarImplicit equals statistics of the corrected matrix", {
    ## pp._get_mean_var_implicit_cov_corr is checked against Python in
    ## test-scdrsPreprocess.R (sparse + cov scenarios); here: chunking and
    ## both margins agree with the explicitly corrected matrix.
    set.seed(1)
    x <- Matrix::rsparsematrix(12, 9, density = 0.4)
    covMat <- cbind(a = rnorm(9), SCDRS_CONST = 1)
    covBeta <- matrix(rnorm(24), 12, 2)
    mu <- rnorm(12)
    full <- as.matrix(x) + covBeta %*% t(covMat) + mu
    for (nChunk in c(1L, 5L, 100L)) {
        for (margin in c("gene", "cell")) {
            res <- scdRs:::.meanVarImplicit(x, covMat, covBeta, mu, margin,
                fun = expm1, nChunk = nChunk)
            expect_equal(res, scdRs:::.meanVar(expm1(full), margin),
                tolerance = 1e-12)
        }
    }
})
