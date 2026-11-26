## Gene set functions compared with Python scdrs.util (load_gs, save_gs,
## load_homolog_mapping) and the `scdrs munge-gs` command (scdrs 1.0.4).
## Fixtures util.load_gs.*, util.save_gs.gs, util.munge.*. No random
## numbers in the functions; the inputs were simulated in Python (seed 3).

## Python's gene set table (trait, pos, gene, weight) as a named list.
readGsTable <- function(name) {
    df <- readReference(name, rowNames = FALSE)
    traits <- unique(df$trait)
    stats::setNames(lapply(traits, function(t) {
        s <- df[df$trait == t, ]
        s <- s[order(s$pos), ]
        stats::setNames(as.numeric(s$weight), s$gene)
    }), traits)
}

ext <- function(f) system.file("extdata", f, package = "RscDRS")

test_that("scdrsReadGeneSets matches load_gs", {
    expect_identical(scdrsReadGeneSets(ext("toydata_mouse.gs")),
        readGsTable("util.load_gs.mouse"))
    expect_identical(scdrsReadGeneSets(ext("toydata_mouse.gs"), "mouse",
        "human"), readGsTable("util.load_gs.mouse2human"))
    expect_identical(scdrsReadGeneSets(ext("toydata_human.gs"), "human",
        "mouse"), readGsTable("util.load_gs.human2mouse"))
    sce <- scdrsLoadToyData()
    expect_identical(scdrsReadGeneSets(ext("toydata_human.gs"), "Homo_sapiens",
        "mmusculus", intersect = rownames(sce)[1:1500]),
        readGsTable("util.load_gs.human2mouse.intersect"))
    ## same species: no mapping
    expect_identical(scdrsReadGeneSets(ext("toydata_mouse.gs"), "mouse",
        "mmusculus"), readGsTable("util.load_gs.mouse"))
})

test_that("scdrsReadGeneSets handles weights and duplicated genes as Python", {
    ## Gene_b given twice: first position, last weight
    res <- scdrsReadGeneSets(referencePath("util.weighted.gs.gz"))
    expect_identical(res, readGsTable("util.load_gs.weighted"))
    expect_identical(res$T1, c(Gene_b = 0.25, Gene_a = 2, Gene_c = -0.001))
    bad <- tempfile(fileext = ".gs")
    writeLines(c("TRAIT\tGENESET", "X\ta:1,b"), bad)
    expect_error(scdrsReadGeneSets(bad), "invalid format")
    expect_error(scdrsReadGeneSets(ext("toydata_mouse.gs"), "mouse"),
        "both")
    expect_error(scdrsReadGeneSets(ext("toydata_mouse.gs"), "rat", "mouse"),
        "not supported")
})

test_that("scdrsWriteGeneSets writes the same file as save_gs", {
    gs <- readGsTable("util.load_gs.weighted")
    f <- tempfile(fileext = ".gs")
    scdrsWriteGeneSets(gs, f)
    ref <- readLines(gzfile(referencePath("util.save_gs.gs.gz")))
    expect_identical(readLines(f), ref)
    ## round trip, and gene-only sets
    expect_identical(scdrsReadGeneSets(f), gs)
    scdrsWriteGeneSets(list(A = c("x", "y")), f)
    expect_identical(readLines(f), c("TRAIT\tGENESET", "A\tx,y"))
})

test_that("scdrsGetHomologs matches load_homolog_mapping", {
    h2m <- scdrsGetHomologs("human", "mouse")
    m2h <- scdrsGetHomologs("mouse", "human")
    expect_length(h2m, 16466L)
    expect_identical(unname(h2m[c("ACADM", "ACADVL")]), c("Acadm", "Acadvl"))
    ## the Python mapping applied to the toy gene set
    ref <- readGsTable("util.load_gs.mouse2human")[[1]]
    gs <- scdrsReadGeneSets(ext("toydata_mouse.gs"))[[1]]
    expect_identical(unique(unname(m2h[names(gs)[names(gs) %in%
        names(m2h)]])), names(ref))
    expect_error(scdrsGetHomologs("human", "Homo_sapiens"), "same")
})

## Python's munge-gs output file as gene sets (weights as written, %.5g).
readMunged <- function(name) {
    scdrsReadGeneSets(referencePath(paste0("util.munge.", name, ".gs.gz")))
}

test_that("scdrsMungeGeneSets matches scdrs munge-gs", {
    pv <- referencePath("util.munge.pval.tsv.gz")
    zs <- referencePath("util.munge.zscore.tsv.gz")
    runs <- list(
        fdr = list(pvals = pv, fdr = 0.05, nMin = 10, nMax = 100),
        fwer = list(pvals = pv, fwer = 0.05, nMin = 5, nMax = 100),
        top = list(pvals = pv, nMin = 10, nMax = 40),
        zscore_uniform = list(zscores = zs, weight = "uniform", fdr = 0.1,
            nMin = 10, nMax = 1000))
    for (nm in names(runs)) {
        f <- tempfile(fileext = ".gs")
        res <- do.call(scdrsMungeGeneSets, c(runs[[nm]], file = f))
        ref <- readMunged(nm)
        expect_identical(names(res), c("TRAIT_A", "TRAIT_B"))
        expect_equal(res, ref, label = nm)
        ## the file written is the same as Python's
        expect_identical(readLines(f), readLines(gzfile(referencePath(
            paste0("util.munge.", nm, ".gs.gz")))), label = nm)
    }
    ## a data.frame input gives the same result as the file
    df <- read.delim(pv, row.names = 1)
    expect_equal(scdrsMungeGeneSets(pvals = df, fdr = 0.05, nMin = 10,
        nMax = 100), readMunged("fdr"))
    expect_error(scdrsMungeGeneSets(), "exactly one")
    expect_error(scdrsMungeGeneSets(pvals = df, fdr = 0.1, fwer = 0.1),
        "at most one")
})
