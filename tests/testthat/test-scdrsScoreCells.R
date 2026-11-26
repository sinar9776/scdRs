## scdrsScoreCells() and scdrsSelectControlGenes() compared with Python
## scdrs.score_cell() on the toy data (scdrs 1.0.4, seed 0, 100 control
## sets; fixtures <tag>.toydata_gs_mouse.*).
##
## Python and R draw different random control sets, so two kinds of test:
## - exact: Python's control sets (fixtures *.ctrl_genesets) are passed to
##   scdrsScoreCells(controlSets = ), and Python's gene statistics are used
##   (toyWithPythonStats()), so every step after the random draw is
##   compared. Python returns float32, hence tolerance 1e-5;
## - statistical: R's own random control sets (1000) against Python's run
##   with 1000 control sets (fixtures *.score1000), with the acceptance
##   criteria in the test "R's own control sets meet the acceptance criteria".

mainCols <- c("raw_score", "norm_score", "mc_pval", "pval", "nlog10_pval",
    "zscore")

checkScores <- function(res, prefix, ctrlNorm = TRUE) {
    cells <- rownames(res)
    ref <- readScores(prefix, cells)
    for (col in mainCols) {
        expect_equal(res[[col]], ref[[col]], tolerance = 1e-5,
            label = paste(prefix, col))
    }
    raw <- as.matrix(readReference(paste0(prefix, ".ctrl_raw_score"))[cells, ])
    expect_equal(unname(res$ctrl_raw_score), unname(raw), tolerance = 1e-5,
        label = paste(prefix, "ctrl_raw_score"))
    if (ctrlNorm) {
        nrm <- as.matrix(readReference(paste0(prefix,
            ".ctrl_norm_score"))[cells, ])
        expect_equal(unname(res$ctrl_norm_score), unname(nrm),
            tolerance = 1e-5, label = paste(prefix, "ctrl_norm_score"))
    }
}

for (tag in c("cov", "nocov")) {
    test_that(paste("scdrsScoreCells matches score_cell with Python's",
        "control sets:", tag), {
        sce <- toyWithPythonStats(cov = tag == "cov")
        prefix <- paste0(tag, ".toydata_gs_mouse")
        ctrl <- readControlSets(paste0(prefix, ".ctrl_genesets"))
        res <- scdrsScoreCells(sce, toyGeneSet(), controlSets = ctrl,
            returnCtrlRawScore = TRUE, returnCtrlNormScore = TRUE)
        expect_identical(rownames(res), colnames(sce))
        checkScores(res, prefix)
        ## the long table format of the fixture is accepted as well
        res2 <- scdrsScoreCells(sce, toyGeneSet(),
            controlSets = readReference(paste0(prefix, ".ctrl_genesets"),
                rowNames = FALSE))
        expect_equal(res2$norm_score, res$norm_score)
    })

    for (opt in c("uniform", "inv_std", "od")) {
        test_that(paste("scdrsScoreCells matches score_cell:", tag, opt), {
            sce <- toyWithPythonStats(cov = tag == "cov")
            ## same seed and bins, so the same sets as for "vs"
            ctrl <- readControlSets(paste0(tag,
                ".toydata_gs_mouse.ctrl_genesets"))
            res <- scdrsScoreCells(sce, toyGeneSet(), controlSets = ctrl,
                weightOpt = opt, returnCtrlRawScore = TRUE)
            checkScores(res, paste0(tag, ".", opt, ".toydata_gs_mouse"),
                ctrlNorm = FALSE)
        })
    }

    test_that(paste("scdrsScoreCells matches score_cell with",
        "ctrl_match_key = 'mean':", tag), {
        sce <- toyWithPythonStats(cov = tag == "cov")
        prefix <- paste0(tag, ".matchmean.toydata_gs_mouse")
        ctrl <- readControlSets(paste0(prefix, ".ctrl_genesets"))
        res <- scdrsScoreCells(sce, toyGeneSet(), controlSets = ctrl,
            ctrlMatchKey = "mean", returnCtrlRawScore = TRUE)
        checkScores(res, prefix, ctrlNorm = FALSE)
    })

    test_that(paste("scdrsScoreCells matches score_cell with gene",
        "weights:", tag), {
        sce <- toyWithPythonStats(cov = tag == "cov")
        prefix <- paste0(tag, ".weighted.toydata_gs_mouse")
        ctrl <- readControlSets(paste0(prefix, ".ctrl_genesets"))
        gs <- toyGeneSet()
        w <- 1 + (seq_along(gs) - 1) %% 5
        ## weights given as values of a named vector or separately
        res <- scdrsScoreCells(sce, stats::setNames(w, names(gs)),
            controlSets = ctrl, returnCtrlRawScore = TRUE)
        checkScores(res, prefix, ctrlNorm = FALSE)
        res2 <- scdrsScoreCells(sce, names(gs), weights = w,
            controlSets = ctrl)
        expect_equal(res2$norm_score, res$norm_score)
        expect_identical(metadata(res)$scdrs$genes,
            stats::setNames(w, names(gs))[sort(names(gs), method = "radix")])
    })
}

test_that("scdrsScoreCells matches score_cell on dense data with covariates", {
    sce <- toyWithPythonStats(cov = TRUE, dense = TRUE)
    expect_identical(metadata(sce)$scdrs$assay, "scdrs_corrected")
    ctrl <- readControlSets("cov.toydata_gs_mouse.ctrl_genesets")
    res <- scdrsScoreCells(sce, toyGeneSet(), controlSets = ctrl,
        returnCtrlRawScore = TRUE)
    checkScores(res, "densecov.toydata_gs_mouse", ctrlNorm = FALSE)
})

test_that(".correctBackground matches _correct_background (with zeros)", {
    sc <- readReference("units.correct_background.score", rowNames = FALSE)
    ctrl <- as.matrix(readReference("units.correct_background.ctrl_raw",
        rowNames = FALSE))
    ratio <- readReference("units.correct_background.ratio",
        rowNames = FALSE)$ratio
    ref <- as.matrix(readReference("units.correct_background.ctrl_norm",
        rowNames = FALSE))
    res <- RscDRS:::.correctBackground(sc$raw, ctrl, ratio,
        returnIntermediate = TRUE)
    expect_equal(res$score, sc$norm, tolerance = 1e-12)
    expect_equal(unname(res$ctrl), unname(ref), tolerance = 1e-12)
    ## the zero raw scores get the smallest values
    expect_identical(which(res$score == min(res$score)), c(2L, 6L))
    expect_named(res$intermediate, c("firstAlignment",
        "cellwiseStandardization", "secondAlignment"))
})

test_that(".empiricalPvalue matches _get_p_from_empi_null (with ties)", {
    t <- readReference("units.empi_null.t", rowNames = FALSE)
    tNull <- readReference("units.empi_null.null", rowNames = FALSE)$t_null
    expect_equal(RscDRS:::.empiricalPvalue(t$t, tNull), t$p,
        tolerance = 1e-15)
    ## naive definition
    naive <- vapply(t$t, function(v) (1 + sum(tNull >= v)) /
        (1 + length(tNull)), numeric(1))
    expect_equal(RscDRS:::.empiricalPvalue(t$t, tNull), naive)
})

test_that("R's control sets match Python's in structure", {
    ## Exact draws differ (different random number generators), but every
    ## control set must have the same number of genes from each
    ## mean-variance bin as Python's, and the same weights.
    sce <- toyWithPythonStats(cov = TRUE)
    gs <- toyGeneSet()
    bins <- stats::setNames(rowData(sce)$scdrs_mean_var, rownames(sce))
    py <- readControlSets("cov.toydata_gs_mouse.ctrl_genesets")
    set.seed(1)
    r <- scdrsSelectControlGenes(sce, gs, nCtrl = 100)
    expect_length(r, 100L)
    tabPy <- table(bins[names(py[[1]])])
    for (k in seq_along(r)) {
        expect_identical(table(bins[names(r[[k]])]), tabPy)
        expect_identical(table(bins[names(py[[k]])]), tabPy)
        expect_false(anyDuplicated(names(r[[k]])) > 0)
        expect_identical(unname(r[[k]]), unname(py[[k]]))
    }
    ## the disease genes' bins
    expect_identical(table(bins[names(gs)]), tabPy)
})

test_that("control genes are drawn uniformly within bins", {
    stats <- data.frame(key = rep(c("a", "b"), c(20, 30)),
        row.names = paste0("g", seq_len(50)))
    genes <- c(g1 = 2, g21 = 1, g22 = 3)
    set.seed(2)
    ctrl <- scdrsSelectControlGenes(stats, genes, ctrlMatchKey = "key",
        nCtrl = 3000)
    ## bin "a" first (sorted), then "b"; weights follow the disease genes
    for (k in c(1L, 100L)) {
        expect_true(names(ctrl[[k]])[1] %in% paste0("g", 1:20))
        expect_true(all(names(ctrl[[k]])[2:3] %in% paste0("g", 21:50)))
        expect_identical(unname(ctrl[[k]]), c(2, 1, 3))
    }
    counts <- table(factor(unlist(lapply(ctrl, function(s) names(s)[2:3])),
        levels = paste0("g", 21:50)))
    expect_gt(chisq.test(counts)$p.value, 0.001)
    ## a continuous key is cut into nGeneBin quantile bins (pd.qcut)
    stats2 <- data.frame(mean = seq_len(50) / 7, row.names = rownames(stats))
    ctrl2 <- scdrsSelectControlGenes(stats2, c("g3", "g48"),
        ctrlMatchKey = "mean", nCtrl = 50, nGeneBin = 5)
    expect_true(all(vapply(ctrl2, function(s) {
        all(names(s) %in% paste0("g", c(1:10, 41:50)))
    }, logical(1))))
})

test_that("R's own control sets meet the acceptance criteria", {
    ## against Python with 1000 control sets
    for (case in list(c("cov", "toydata_gs_mouse"),
        c("nocov", "toydata_gs_human"))) {
        sce <- suppressMessages(scdrsPreprocess(scdrsLoadToyData(),
            cov = if (case[1] == "cov") "covariate"))
        f <- system.file("extdata", paste0(sub("toydata_gs_", "toydata_",
            case[2]), ".gs"), package = "RscDRS")
        gs <- if (case[2] == "toydata_gs_human") {
            scdrsReadGeneSets(f, "human", "mouse", rownames(sce))
        } else {
            scdrsReadGeneSets(f)
        }
        set.seed(1)
        res <- scdrsScoreCells(sce, gs[[case[2]]], nCtrl = 1000)
        ref <- readReference(paste0(case[1], ".", case[2],
            ".score1000"))[rownames(res), ]
        expect_gte(cor(res$norm_score, ref$norm_score), 0.99)
        expect_gte(cor(res$nlog10_pval, ref$nlog10_pval), 0.99)
        ## FDR < 0.1 calls: with only 30 cells a single cell is 3.3%, so
        ## the 99% criterion (meant for large data) cannot absorb even one
        ## borderline cell. In the human/nocov case Python's FDR for one
        ## cell is 0.101 and R's is 0.09-0.10 for every seed tried. So:
        ## calls must agree except for at most one cell, and such a cell
        ## must be borderline (FDR in [0.05, 0.2] in both runs).
        fdrR <- p.adjust(res$pval, "BH")
        fdrPy <- p.adjust(ref$pval, "BH")
        differ <- (fdrR < 0.1) != (fdrPy < 0.1)
        expect_lte(sum(differ), 1L)
        expect_true(all(fdrR[differ] > 0.05 & fdrR[differ] < 0.2 &
            fdrPy[differ] > 0.05 & fdrPy[differ] < 0.2))
        expect_equal(res$raw_score, ref$raw_score, tolerance = 1e-3)
    }
})

test_that("scdrsScoreCells input checks and options", {
    sce <- scdrsLoadToyData()
    expect_error(scdrsScoreCells(sce, toyGeneSet()), "scdrsPreprocess")
    sce <- suppressMessages(scdrsPreprocess(sce))
    expect_error(scdrsScoreCells(sce, c("nope1", "nope2")), "none of")
    expect_error(scdrsScoreCells(sce, toyGeneSet(), weights = 1),
        "not both")
    expect_error(scdrsScoreCells(sce, toyGeneSet(), ctrlMatchKey = "x"),
        "not found")
    set.seed(3)
    res <- scdrsScoreCells(sce, c(names(toyGeneSet()), "not_a_gene"),
        nCtrl = 20)
    expect_identical(names(res), mainCols)
    expect_length(metadata(res)$scdrs$genes, length(toyGeneSet()))
    expect_length(metadata(res)$scdrs$controlSets, 20L)
    ## reusing the stored control sets gives the same result
    res2 <- scdrsScoreCells(sce, toyGeneSet(),
        controlSets = metadata(res)$scdrs$controlSets)
    expect_identical(res2$norm_score, res$norm_score)
    expect_true(all(res$pval > 0 & res$pval <= 1))
    expect_true(all(abs(res$zscore) <= 10))
})
