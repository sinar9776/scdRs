# scdRs

A native R implementation of [scDRS](https://github.com/martinjzhang/scDRS)
(single-cell disease-relevance score; Zhang, Hou et al., *Nature Genetics*
2022) for `SingleCellExperiment` objects.

scDRS links single-cell RNA-seq data to polygenic disease risk. It scores
every cell for excess expression of a disease gene set, such as the top
genes of a MAGMA analysis of GWAS results, compared with matched control
gene sets. It then tests cell groups for disease association and for
heterogeneity. scdRs follows the Python package (scdrs 1.0.4) and is tested
against it, without calling Python.

## Installation

scdRs has been submitted to Bioconductor. Once it is accepted:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("scdRs")
```

Until then, install the development version from GitHub. It requires
R >= 4.6 and Bioconductor devel:

```r
BiocManager::install("sinar9776/scdRs")
```

## Example

The package includes the scDRS toy data: 2,500 mouse genes and 30 cells,
10 of which overexpress the disease genes.

```r
library(scdRs)

sce <- scdrsLoadToyData()
geneSets <- scdrsReadGeneSets(system.file("extdata", "toydata_mouse.gs",
    package = "scdRs"))

## covariate correction and gene statistics
sce <- scdrsPreprocess(sce, cov = "covariate")

## score each cell against 1,000 matched control gene sets
set.seed(1)
res <- scdrsScoreCells(sce, geneSets[["toydata_gs_mouse"]],
    returnCtrlNormScore = TRUE)
table(sce$cell_type, FDR_below_0.1 = p.adjust(res$pval, "BH") < 0.1)

## group-level association and heterogeneity
sce <- scdrsBuildGraph(sce)
scdrsTestGroups(sce, res, groupCols = "cell_type")
```

The vignette has the full workflow, including gene sets from GWAS
results, several traits at once, and the correlation and gene analyses:

```r
vignette("scdRs")
```

## Development

scdRs was developed with the assistance of an AI coding assistant
(Claude, Anthropic). All code has been reviewed by the maintainer, and
every exported function is tested against the output of the original
Python scDRS (see `tests/testthat/` and `inst/scripts/python/`).

## Citation

scdRs is a port of the MIT-licensed Python package by Martin Jinye Zhang
and Kangcheng Hou (see `inst/COPYRIGHTS`). If you use it, please cite the
original paper (`citation("scdRs")`):

> Zhang MJ, Hou K, et al. (2022). Polygenic enrichment distinguishes
> disease associations of individual cells in single-cell RNA-seq data.
> *Nature Genetics* 54, 1572-1580. doi:10.1038/s41588-022-01167-z
