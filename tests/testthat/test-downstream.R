## Downstream analyses compared with Python scdrs.method on the toy data
## (scdrs 1.0.4, seed 0). Inputs are Python's: the full score file of the
## covariate-corrected mouse gene set with 100 control sets
## (fixtures/downstream.toydata_gs_mouse.full_score.gz) and scanpy's kNN
## graph (fixtures/downstream.connectivities), so the comparisons are exact
## (no random numbers). Python returns float32: tolerance 1e-5.

pyScores <- function() {
    f <- referencePath("downstream.toydata_gs_mouse.full_score.gz")
    scdrsReadScores(f)[["downstream.toydata_gs_mouse"]]
}

pyGraph <- function(cells) {
    e <- readReference("downstream.connectivities", rowNames = FALSE)
    Matrix::sparseMatrix(i = match(e$cell_i, cells),
        j = match(e$cell_j, cells), x = e$weight,
        dims = c(length(cells), length(cells)),
        dimnames = list(cells, cells))
}

## Toy data with Python's graph and the 3-group annotation the reference
## script adds (cells taken in turn).
toyWithGraph <- function() {
    sce <- scdrsLoadToyData()
    sce$test_group <- paste0("g", (seq_len(ncol(sce)) - 1L) %% 3L)
    colPair(sce, "connectivities") <- unname(pyGraph(colnames(sce)))
    sce
}

test_that("scdrsReadScores reads Python's .full_score.gz", {
    s <- pyScores()
    expect_s4_class(s, "DataFrame")
    expect_identical(dim(s$ctrl_norm_score), c(30L, 100L))
    expect_identical(dim(s$ctrl_raw_score), c(30L, 100L))
    ref <- readReference("cov.toydata_gs_mouse.score")
    expect_equal(s[rownames(ref), "norm_score"], ref$norm_score)
    nrm <- readReference("cov.toydata_gs_mouse.ctrl_norm_score")
    expect_equal(unname(s$ctrl_norm_score[match(rownames(nrm),
        rownames(s)), ]), unname(as.matrix(nrm)))
})

for (col in c("cell_type", "test_group")) {
    test_that(paste("scdrsTestGroups matches downstream_group_analysis:",
        col), {
        sce <- toyWithGraph()
        res <- scdrsTestGroups(sce, pyScores(), col)
        expect_named(res, col)
        res <- res[[col]]
        ref <- readReference(paste0("downstream.toydata_gs_mouse.group.",
            col))
        expect_identical(rownames(res), rownames(ref))
        expect_identical(names(res), names(ref))
        for (nm in names(ref)) {
            expect_equal(res[[nm]], ref[[nm]], tolerance = 1e-5,
                label = paste(col, nm))
        }
    })
}

test_that("scdrsTestGroups takes a graph matrix and handles empty groups", {
    sce <- scdrsLoadToyData()
    g <- pyGraph(colnames(sce))
    res <- scdrsTestGroups(sce, pyScores(), "cell_type", graph = g)
    ref <- readReference("downstream.toydata_gs_mouse.group.cell_type")
    expect_equal(res$cell_type$assoc_mcz, ref$assoc_mcz, tolerance = 1e-5)
    expect_error(scdrsTestGroups(sce, pyScores(), "cell_type"),
        "scdrsBuildGraph")
    ## a group whose cells were not scored gets NA (Python fails)
    s <- pyScores()
    causal <- colnames(sce)[sce$cell_type == "causal_cell"]
    expect_warning(res2 <- scdrsTestGroups(sce, s[!rownames(s) %in% causal, ],
        "cell_type", graph = g), "no scored cells")
    expect_true(is.na(res2$cell_type["causal_cell", "assoc_mcp"]))
    expect_identical(res2$cell_type["causal_cell", "n_cell"], 0)
})

test_that("scdrsTestCorrelation matches downstream_corr_analysis", {
    vars <- c("causal_variable", "non_causal_variable", "covariate")
    res <- scdrsTestCorrelation(scdrsLoadToyData(), pyScores(), vars)
    ref <- readReference("downstream.toydata_gs_mouse.corr")
    expect_identical(rownames(res), vars)
    for (nm in names(ref)) {
        expect_equal(res[[nm]], ref[[nm]], tolerance = 1e-5, label = nm)
    }
})

test_that("scdrsRankGenes matches downstream_gene_analysis", {
    res <- scdrsRankGenes(scdrsLoadToyData(), pyScores())
    ref <- readReference("downstream.toydata_gs_mouse.gene")
    expect_setequal(rownames(res), rownames(ref))
    expect_equal(res[rownames(ref), "CORR"], ref$CORR, tolerance = 1e-5)
    ## ranks are 1-based in R; compare where float32 CORR has no ties
    single <- !ref$CORR %in% ref$CORR[duplicated(ref$CORR)]
    expect_identical(res[rownames(ref)[single], "RANK"] - 1L,
        as.integer(ref$RANK[single]))
    expect_false(is.unsorted(rev(res$CORR)))
})

test_that("scdrsComputeGearysC matches gearys_c", {
    sce <- scdrsLoadToyData()
    s <- pyScores()[colnames(sce), ]
    ref <- readReference("downstream.gearys_c", rowNames = FALSE)
    vals <- cbind(s$norm_score, s$ctrl_norm_score[, 1:5])
    expect_equal(scdrsComputeGearysC(pyGraph(colnames(sce)), vals), ref$C,
        tolerance = 1e-7)
    expect_equal(scdrsComputeGearysC(as.matrix(pyGraph(colnames(sce))),
        s$norm_score), ref$C[1], tolerance = 1e-7)
    ## the direct definition, on a small asymmetric graph
    set.seed(4)
    w <- matrix(runif(36) * (runif(36) > 0.5), 6)
    x <- rnorm(6)
    direct <- 5 * sum(w * outer(x, x, "-")^2) /
        (2 * sum(w) * sum((x - mean(x))^2))
    expect_equal(scdrsComputeGearysC(w, x), direct)
})

for (opt in c("control_distribution_match", "control")) {
    test_that(paste("scdrsTestGearysC matches test_gearysc:", opt), {
        sce <- scdrsLoadToyData()
        s <- pyScores()[colnames(sce), ]
        res <- scdrsTestGearysC(pyGraph(colnames(sce)), s, sce$cell_type,
            method = opt)
        ref <- readReference(paste0("downstream.gearysc.", opt))
        expect_identical(rownames(res), rownames(ref))
        for (nm in names(ref)) {
            expect_equal(res[[nm]], ref[[nm]], tolerance = 1e-6,
                label = paste(opt, nm))
        }
    })
}

test_that("scdrsTestGearysC permutation option behaves as a null", {
    ## Random permutations (no Python parity possible): the trait statistic
    ## is the same as with the other options and the p-value is in range.
    sce <- scdrsLoadToyData()
    s <- pyScores()[colnames(sce), ]
    g <- pyGraph(colnames(sce))
    set.seed(5)
    res <- scdrsTestGearysC(g, s, sce$cell_type, method = "permutation")
    ref <- readReference("downstream.gearysc.control")
    expect_equal(res$trait, ref$trait, tolerance = 1e-6)
    expect_true(all(res$pval > 0 & res$pval <= 1))
    ## a permutation null has mean close to 1 (Geary's C of random values)
    expect_true(all(abs(res$ctrl_mean - 1) < 0.15))
})

test_that(".pearsonCor matches _pearson_corr (dense and sparse)", {
    x <- as.matrix(readReference("units.pearson.X", rowNames = FALSE))
    y <- as.matrix(readReference("units.pearson.Y", rowNames = FALSE))
    ref <- readReference("units.pearson.output", rowNames = FALSE)
    for (kind in c("dense", "sparse")) {
        r <- ref[ref$kind == kind, ]
        xx <- if (kind == "sparse") Matrix::Matrix(x, sparse = TRUE) else x
        res <- RscDRS:::.pearsonCor(xx, y)
        expect_equal(res[cbind(r$i + 1, r$j + 1)], r$corr, tolerance = 1e-6,
            label = kind)
    }
})

test_that("scdrsPerformDownstream runs all analyses like the CLI", {
    sce <- toyWithGraph()
    s <- pyScores()
    out <- scdrsPerformDownstream(sce, list(BMI = s), groupCols =
        c("cell_type", "test_group"), varCols = "causal_variable",
        geneAnalysis = TRUE)
    expect_named(out, "BMI")
    expect_named(out$BMI, c("group", "corr", "gene"))
    ref <- readReference("downstream.toydata_gs_mouse.group.test_group")
    expect_equal(out$BMI$group$test_group$hetero_mcz, ref$hetero_mcz,
        tolerance = 1e-5)
    one <- scdrsPerformDownstream(sce, s, varCols = "causal_variable")
    expect_null(one$group)
    expect_equal(one$corr, out$BMI$corr)
    expect_error(scdrsPerformDownstream(sce, s), "at least one")
})
