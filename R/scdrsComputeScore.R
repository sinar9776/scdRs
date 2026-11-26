## Provenance: R port of the Python scDRS 1.0.4
## (https://github.com/martinjzhang/scDRS; MIT License, Copyright (c)
## 2021 Martin Jinye Zhang and Kangcheng Hou; see inst/COPYRIGHTS).
## Assisted-by: Claude (Anthropic), reviewed by the maintainer.

#' Filter and normalize raw counts as scDRS does
#'
#' Optionally removes cells with fewer than `minGenes` expressed genes and
#' then genes expressed in fewer than `minCells` cells, and normalizes each
#' cell to 10,000 total counts followed by a natural-log `log1p`. This
#' copies what `scdrs compute-score` does with `--flag-filter-data` and
#' `--flag-raw-count` (scanpy's `filter_cells`, `filter_genes`,
#' `normalize_per_cell(counts_per_cell_after = 1e4)` and `log1p`; R port
#' of `scdrs.util.load_h5ad()`).
#'
#' Note that `scuttle::logNormCounts()` gives different values (log2 and
#' centred size factors); scDRS expects natural-log data.
#'
#' @param sce A [SingleCellExperiment::SingleCellExperiment] with raw
#'     counts.
#' @param assay Name of the counts assay.
#' @param filter Filter cells and genes first?
#' @param minGenes,minCells Filtering thresholds.
#' @param normalize Normalize and log-transform? If `FALSE`, only the
#'     filtering is done.
#' @param name Name of the new assay.
#'
#' @return `sce`, possibly with fewer cells and genes, with the normalized
#'     data in assay `name` (sparse if the counts are sparse). Cells
#'     without any counts are removed when normalizing, as in scanpy.
#'
#' @examples
#' set.seed(1)
#' counts <- matrix(rpois(2000, 2), 50, 40,
#'     dimnames = list(paste0("g", 1:50), paste0("c", 1:40)))
#' sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
#' sce <- scdrsNormalize(sce, minGenes = 10, minCells = 5)
#' sce
#'
#' @export
scdrsNormalize <- function(sce, assay = "counts", filter = TRUE,
    minGenes = 250L, minCells = 50L, normalize = TRUE, name = "logcounts") {
    x <- .getAssay(sce, assay)
    if (anyNA(if (is(x, "sparseMatrix")) x@x else x)) {
        stop("the counts must not contain missing values")
    }
    if (any((if (is(x, "sparseMatrix")) x@x else x) < 0)) {
        stop("the counts must not contain negative values")
    }
    if (filter) {
        keep <- Matrix::colSums(x > 0) >= minGenes
        sce <- sce[, keep]
        x <- x[, keep, drop = FALSE]
        keep <- Matrix::rowSums(x > 0) >= minCells
        sce <- sce[keep, ]
        x <- x[keep, , drop = FALSE]
    }
    if (normalize) {
        if (name == assay) {
            stop("'name' must differ from 'assay': the input assay is ",
                "never overwritten")
        }
        total <- Matrix::colSums(x)
        keep <- total >= 1
        sce <- sce[, keep]
        x <- x[, keep, drop = FALSE]
        f <- 1e4 / total[keep]
        if (is(x, "sparseMatrix")) {
            x@x <- log1p(x@x * rep(f, diff(x@p)))
        } else {
            x <- log1p(sweep(x, 2L, f, "*"))
        }
        SummarizedExperiment::assay(sce, name, withDimnames = FALSE) <- x
    }
    sce
}

#' Compute scDRS scores for several gene sets
#'
#' Preprocesses the data if needed ([scdrsPreprocess()]) and scores every
#' gene set with at least `minGenes` genes in the data
#' ([scdrsScoreCells()]). This is the R port of the `scdrs compute-score`
#' command, with the same defaults.
#'
#' @param sce A [SingleCellExperiment::SingleCellExperiment] with
#'     log-normalized data (see [scdrsNormalize()] for raw counts).
#' @param geneSets A named list of gene sets as returned by
#'     [scdrsReadGeneSets()] (named numeric vectors of weights, or
#'     character vectors), or the path of a `.gs` file with genes named as
#'     in `sce`.
#' @param assay Assay to use.
#' @param cov,adjProp Covariates and cell-group column for
#'     [scdrsPreprocess()]. The data is preprocessed when these are given
#'     or when `sce` has not been preprocessed yet; otherwise the existing
#'     preprocessing is used.
#' @param ctrlMatchKey,weightOpt,nCtrl See [scdrsScoreCells()].
#' @param returnCtrlRawScore,returnCtrlNormScore Return the control
#'     scores? The normalized ones are returned by default, because the
#'     downstream analyses need them.
#' @param minGenes Gene sets with fewer genes in `sce` are skipped.
#'
#' @return A named [S4Vectors::SimpleList] of score tables (see
#'     [scdrsScoreCells()]), one per scored gene set.
#'
#' @seealso [scdrsWriteScores()], [scdrsPerformDownstream()].
#'
#' @examples
#' sce <- scdrsLoadToyData()
#' gs <- scdrsReadGeneSets(system.file("extdata", "toydata_human.gs",
#'     package = "RscDRS"), srcSpecies = "human", dstSpecies = "mouse")
#' set.seed(1)
#' res <- scdrsComputeScore(sce, gs, cov = "covariate", nCtrl = 100)
#' names(res)
#' head(res[[1]][, 1:6])
#'
#' @export
scdrsComputeScore <- function(sce, geneSets, assay = "logcounts",
    cov = NULL, adjProp = NULL, ctrlMatchKey = "mean_var",
    weightOpt = "vs", nCtrl = 1000L, returnCtrlRawScore = FALSE,
    returnCtrlNormScore = TRUE, minGenes = 10L) {
    if (is.character(geneSets) && length(geneSets) == 1L) {
        geneSets <- scdrsReadGeneSets(geneSets)
    }
    if (!is.list(geneSets) || is.null(names(geneSets))) {
        stop("'geneSets' must be a named list of gene sets or a .gs file")
    }
    if (is.null(S4Vectors::metadata(sce)$scdrs) || !is.null(cov) ||
        !is.null(adjProp)) {
        sce <- scdrsPreprocess(sce, assay = assay, cov = cov,
            adjProp = adjProp)
    }
    out <- list()
    for (trait in names(geneSets)) {
        gs <- geneSets[[trait]]
        genes <- if (is.numeric(gs)) names(gs) else gs
        n <- length(unique(intersect(genes, rownames(sce))))
        if (n < minGenes) {
            message("trait=", trait, ": skipped due to small size (n_gene=",
                n, ")")
            next
        }
        res <- scdrsScoreCells(sce, gs, ctrlMatchKey = ctrlMatchKey,
            nCtrl = nCtrl, weightOpt = weightOpt,
            returnCtrlRawScore = returnCtrlRawScore,
            returnCtrlNormScore = returnCtrlNormScore)
        fdr <- stats::p.adjust(res$pval, "BH")
        message("Trait=", trait, ", n_gene=", n, ": ", sum(fdr < 0.1), "/",
            nrow(res), " FDR<0.1 cells, ", sum(fdr < 0.2), "/", nrow(res),
            " FDR<0.2 cells")
        out[[trait]] <- res
    }
    S4Vectors::SimpleList(out)
}

## Score table -> plain data.frame with Python's column names
## (ctrl_raw_score_0, ..., ctrl_norm_score_0, ...).
.flattenScores <- function(scores) {
    main <- c("raw_score", "norm_score", "mc_pval", "pval", "nlog10_pval",
        "zscore")
    df <- as.data.frame(scores[, intersect(main, names(scores)),
        drop = FALSE], optional = TRUE)
    for (nm in c("ctrl_raw_score", "ctrl_norm_score")) {
        if (nm %in% names(scores)) {
            m <- as.matrix(scores[[nm]])
            colnames(m) <- paste0(nm, "_", seq_len(ncol(m)) - 1L)
            df <- cbind(df, m)
        } else {
            cols <- grep(paste0("^", nm, "_[0-9]+$"), names(scores),
                value = TRUE)
            if (length(cols) > 0L) {
                df <- cbind(df, as.data.frame(scores[, cols, drop = FALSE],
                    optional = TRUE))
            }
        }
    }
    rownames(df) <- rownames(scores)
    df
}

#' Read and write scDRS score files
#'
#' `scdrsWriteScores()` writes score tables in the format of `scdrs
#' compute-score`: `<trait>.score.gz` (the six main columns) and, if the
#' tables contain control scores, `<trait>.full_score.gz` (all columns,
#' with `ctrl_raw_score_0`, ..., `ctrl_norm_score_0`, ...).
#' `scdrsReadScores()` reads such files, written by R or by Python, back
#' into score tables that the downstream functions accept. This is the R
#' port of `scdrs.util.load_scdrs_score()`.
#'
#' @param scores A named list of score tables (e.g. from
#'     [scdrsComputeScore()]).
#' @param dir Output directory (created if needed).
#'
#' @return `scdrsWriteScores()`: the paths of the files written,
#'     invisibly.
#'
#' @examples
#' sce <- scdrsLoadToyData()
#' gs <- scdrsReadGeneSets(system.file("extdata", "toydata_mouse.gs",
#'     package = "RscDRS"))
#' set.seed(1)
#' res <- scdrsComputeScore(sce, gs, nCtrl = 20)
#' dir <- file.path(tempdir(), "scdrs_out")
#' scdrsWriteScores(res, dir)
#' list.files(dir)
#' back <- scdrsReadScores(file.path(dir, "@.full_score.gz"))
#' all.equal(back$toydata_gs_mouse$norm_score, res[[1]]$norm_score)
#'
#' @export
scdrsWriteScores <- function(scores, dir) {
    if (is.data.frame(scores) || is(scores, "DataFrame") ||
        is.null(names(scores))) {
        stop("'scores' must be a named list of score tables")
    }
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    files <- character(0)
    write <- function(df, f) {
        con <- gzfile(f, "w")
        on.exit(close(con))
        utils::write.table(df, con, sep = "\t", quote = FALSE,
            col.names = NA)
    }
    for (trait in names(scores)) {
        df <- .flattenScores(scores[[trait]])
        f <- file.path(dir, paste0(trait, ".score.gz"))
        write(df[, seq_len(min(6L, ncol(df))), drop = FALSE], f)
        files <- c(files, f)
        if (ncol(df) > 6L) {
            f <- file.path(dir, paste0(trait, ".full_score.gz"))
            write(df, f)
            files <- c(files, f)
        }
    }
    invisible(files)
}

#' @rdname scdrsWriteScores
#' @param file Path of a `.full_score.gz` or `.score.gz` file. Use `@` in
#'     the file name to read several files at once, e.g.
#'     `"out/@.full_score.gz"` (the directory itself must not contain
#'     `@`).
#' @param cells Optional cell names: files containing fewer than 10% of
#'     them are skipped, with a message.
#' @return `scdrsReadScores()`: a named list (traits from the file names)
#'     of [S4Vectors::DataFrame]s with cell names as row names; control
#'     scores are returned as the matrix columns `ctrl_raw_score` and
#'     `ctrl_norm_score`.
#' @export
scdrsReadScores <- function(file, cells = NULL) {
    dir <- dirname(file)
    pattern <- basename(file)
    if (!grepl("(full_score|score)\\.gz$", pattern)) {
        stop("expect scDRS .full_score.gz or .score.gz files")
    }
    rx <- paste0("^", gsub("@", ".*", gsub(".", "\\.", pattern,
        fixed = TRUE), fixed = TRUE), "$")
    files <- sort(list.files(dir, pattern = rx), method = "radix")
    if (length(files) == 0L) {
        stop("no score files found for '", file, "'")
    }
    out <- list()
    for (f in files) {
        df <- utils::read.delim(file.path(dir, f), check.names = FALSE,
            row.names = 1L, colClasses = c("character"))
        df[] <- lapply(df, as.numeric)
        if (!is.null(cells)) {
            n <- length(intersect(cells, rownames(df)))
            if (n < 0.1 * length(cells)) {
                message(f, " skipped, containing ", n, "/", length(cells),
                    " target cells")
                next
            }
        }
        res <- S4Vectors::DataFrame(df[, !grepl("^ctrl_(raw|norm)_score_",
            names(df)), drop = FALSE], check.names = FALSE)
        for (nm in c("ctrl_raw_score", "ctrl_norm_score")) {
            cols <- grep(paste0("^", nm, "_[0-9]+$"), names(df), value = TRUE)
            if (length(cols) > 0L) {
                id <- as.integer(sub(paste0(nm, "_"), "", cols))
                cols <- cols[order(id)]
                res[[nm]] <- matrix(as.matrix(df[, cols, drop = FALSE]),
                    ncol = length(cols),
                    dimnames = list(NULL, as.character(sort(id))))
            }
        }
        out[[sub("\\.(full_score|score)\\.gz$", "", f)]] <- res
    }
    out
}
