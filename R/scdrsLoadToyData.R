## Provenance: R port of the Python scDRS 1.0.4
## (https://github.com/martinjzhang/scDRS; MIT License, Copyright (c)
## 2021 Martin Jinye Zhang and Kangcheng Hou; see inst/COPYRIGHTS).
## Assisted-by: Claude (Anthropic), reviewed by the maintainer.

#' Load the scDRS toy data
#'
#' Loads the small example data set shipped with the original scDRS (a
#' subset of Tabula Muris Senis: 2,500 mouse genes and 30 cells, 10
#' "causal" and 20 "non-causal") as a `SingleCellExperiment`.
#'
#' @return A [SingleCellExperiment::SingleCellExperiment] with a sparse
#'     `logcounts` assay (size-factor normalized, `log1p`), cell metadata
#'     in `colData` (`label`, `cell_type`, `causal_variable`,
#'     `non_causal_variable`, `covariate`) and gene metadata in `rowData`
#'     (`label`). The files are in `system.file("extdata", package =
#'     "RscDRS")`, including the disease gene sets `toydata_mouse.gs` and
#'     `toydata_human.gs`; see `inst/scripts/python/README.md` for their
#'     origin.
#'
#' @examples
#' sce <- scdrsLoadToyData()
#' sce
#' table(sce$cell_type)
#'
#' @export
scdrsLoadToyData <- function() {
    ext <- system.file("extdata", package = "RscDRS")
    read <- function(f) {
        utils::read.csv(file.path(ext, f), row.names = 1,
            check.names = FALSE)
    }
    expr <- as.matrix(read("toydata_logexpr.csv.gz"))
    cells <- read("toydata_cell_metadata.csv")
    genes <- read("toydata_gene_metadata.csv")
    SingleCellExperiment::SingleCellExperiment(
        assays = list(logcounts = Matrix::Matrix(expr, sparse = TRUE)),
        colData = S4Vectors::DataFrame(cells[colnames(expr), , drop = FALSE],
            check.names = FALSE),
        rowData = S4Vectors::DataFrame(genes[rownames(expr), , drop = FALSE],
            check.names = FALSE))
}
