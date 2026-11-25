# scdRs

A native R implementation of [scDRS](https://github.com/martinjzhang/scDRS)
(single-cell disease-relevance score; Zhang, Hou et al., *Nature Genetics*
2022) for `SingleCellExperiment` objects, intended for Bioconductor.

**Status:** under development. The full scDRS workflow (gene sets,
preprocessing, cell scoring and the downstream analyses) is implemented
and tested against the Python package; see the vignette
(`vignette("scdRs")`). Not yet submitted to Bioconductor.

```r
sce <- scdrsPreprocess(sce, cov = "covariate")
set.seed(1)
res <- scdrsScoreCells(sce, geneSets[["BMI"]], returnCtrlNormScore = TRUE)
sce <- scdrsBuildGraph(sce)
scdrsTestGroups(sce, res, groupCols = "cell_type")
```

scdRs is a port of the original MIT-licensed Python package by Martin
Jinye Zhang and Kangcheng Hou; see `inst/COPYRIGHTS`. If you use it,
please cite the original paper (`citation("scdRs")`).
