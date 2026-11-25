## Small statistics utilities compared with Python scdrs.util
## (test_overlap, meta_analysis, zsc2pval, pval2zsc; scdrs 1.0.4).
## Fixtures util.test_overlap, util.meta_analysis, util.zsc_pval. No random
## numbers.

test_that("scdrsTestOverlap matches test_overlap", {
    bg <- paste0("g", 0:99)
    cases <- list(basic = list(bg[1:30], bg[21:45]),
        no_overlap = list(bg[1:10], bg[51:60]),
        subset = list(bg[1:10], bg[1:40]))
    ref <- readReference("util.test_overlap")
    for (nm in names(cases)) {
        res <- scdrsTestOverlap(cases[[nm]][[1]], cases[[nm]][[2]], bg)
        expect_named(res, c("pval", "odds_ratio", "or_ub", "or_lb"))
        r <- unlist(ref[nm, ])
        expect_equal(unname(res), unname(r), tolerance = 1e-8, label = nm)
    }
})

test_that("scdrsMetaAnalyze matches meta_analysis", {
    inp <- readReference("util.meta_analysis.input", rowNames = FALSE)
    ref <- readReference("util.meta_analysis", rowNames = FALSE)
    for (i in seq_len(nrow(ref))) {
        w <- if (ref$weights[i] == "w") c(1, 2, 1, 3, 1)
        res <- scdrsMetaAnalyze(inp$effect, inp$se, method = ref$method[i],
            weights = w)
        expect_equal(unname(res), c(ref$summary[i], ref$se[i]),
            tolerance = 1e-12)
    }
    expect_named(scdrsMetaAnalyze(inp$effect, inp$se), c("summary", "se"))
})

test_that("scdrsZToP and scdrsPToZ match zsc2pval and pval2zsc", {
    ref <- readReference("util.zsc_pval", rowNames = FALSE)
    expect_equal(scdrsZToP(ref$z), ref$p, tolerance = 1e-12)
    expect_equal(scdrsPToZ(ref$p), ref$z_back, tolerance = 1e-10)
})
