# The example data in inst/extdata is the Python scDRS toy data. Values
# are compared with toydata_mouse.h5ad as exported by
# inst/scripts/python/make_reference_outputs.py (scdrs 1.0.4).
test_that("toy data in extdata matches the Python toy data", {
    ext <- system.file("extdata", package = "scdRs")
    expr <- as.matrix(read.csv(file.path(ext, "toydata_logexpr.csv.gz"),
        row.names = 1, check.names = FALSE))
    cells <- read.csv(file.path(ext, "toydata_cell_metadata.csv"),
        row.names = 1)
    cov <- read.csv(file.path(ext, "toydata_covariates.csv"),
        row.names = 1)
    expect_identical(dim(expr), c(2500L, 30L))
    expect_setequal(colnames(expr), rownames(cells))
    expect_setequal(rownames(cov), rownames(cells))
    expect_true(all(expr >= 0))
    expect_identical(sort(unique(cells$cell_type)),
        c("causal_cell", "non_causal_cell"))
    gs <- read.delim(file.path(ext, "toydata_mouse.gs"))
    expect_identical(gs$TRAIT, "toydata_gs_mouse")
})
