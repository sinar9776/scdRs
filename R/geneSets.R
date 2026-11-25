## Gene set files (.gs), mouse-human homologs and munge-gs
## (port of scdrs/util.py load_gs, save_gs, load_homolog_mapping,
## convert_species_name, and of the `scdrs munge-gs` command).

## Port of util.convert_species_name().
.speciesName <- function(species) {
    if (species %in% c("Mouse", "mouse", "Mus_musculus", "mus_musculus",
        "mmusculus")) {
        return("mmusculus")
    }
    if (species %in% c("Human", "human", "Homo_sapiens", "homo_sapiens",
        "hsapiens")) {
        return("hsapiens")
    }
    stop("species name '", species, "' is not supported")
}

## Keep the first position of each name and its last value (what a Python
## dict does when a key is assigned twice).
.dictDedup <- function(x) {
    nm <- names(x)
    if (!anyDuplicated(nm)) {
        return(x)
    }
    u <- unique(nm)
    last <- length(nm) + 1L - match(u, rev(nm))
    stats::setNames(x[last], u)
}

#' Mouse-human gene homologs
#'
#' Returns the mouse-human gene symbol map shipped with scDRS (16,466
#' pairs, `system.file("extdata", "mouse_human_homologs.txt", package =
#' "scdRs")`). This is the R port of `scdrs.util.load_homolog_mapping()`.
#'
#' @param src,dst Source and destination species: `"mouse"` or `"human"`
#'     (also `"mmusculus"`, `"Mus_musculus"`, `"hsapiens"`,
#'     `"Homo_sapiens"`, ...). They must differ.
#'
#' @return A named character vector: names are `src` gene symbols, values
#'     the `dst` homologs.
#'
#' @examples
#' map <- scdrsGetHomologs("human", "mouse")
#' map[c("ACADM", "APOE")]
#'
#' @export
scdrsGetHomologs <- function(src, dst) {
    src <- .speciesName(src)
    dst <- .speciesName(dst)
    if (src == dst) {
        stop("'src' and 'dst' cannot be the same species")
    }
    tab <- utils::read.delim(system.file("extdata",
        "mouse_human_homologs.txt", package = "scdRs"),
        colClasses = "character", quote = "", na.strings = character(0))
    from <- if (src == "hsapiens") "HUMAN_GENE_SYM" else "MOUSE_GENE_SYM"
    to <- if (src == "hsapiens") "MOUSE_GENE_SYM" else "HUMAN_GENE_SYM"
    .dictDedup(stats::setNames(tab[[to]], tab[[from]]))
}

#' Read and write scDRS gene set files
#'
#' `scdrsReadGeneSets()` reads a `.gs` file: a tab-separated file with
#' columns `TRAIT` and `GENESET`, where `GENESET` is a comma-separated
#' list of genes (`"g1,g2,g3"`) or of genes and weights
#' (`"g1:w1,g2:w2"`). Genes can be mapped between mouse and human and
#' intersected with the genes of a data set. `scdrsWriteGeneSets()` writes
#' such a file. These are the R ports of `scdrs.util.load_gs()` and
#' `scdrs.util.save_gs()`.
#'
#' @param file Path of the `.gs` file.
#' @param srcSpecies,dstSpecies `NULL` (default, no mapping), or the
#'     species of the gene set file and of the data (`"mouse"` or
#'     `"human"`). Genes without a homolog are dropped.
#' @param intersect `NULL`, or a character vector of genes (e.g.
#'     `rownames(sce)`) to keep.
#'
#' @return `scdrsReadGeneSets()`: a named list (one element per trait) of
#'     named numeric vectors (names = genes, values = weights; 1 for an
#'     unweighted gene set). A gene listed twice keeps its first position
#'     and its last weight, as in Python.
#'
#' @seealso [scdrsMungeGeneSets()] to make gene sets from GWAS gene
#'     statistics, [scdrsComputeScore()].
#'
#' @examples
#' f <- system.file("extdata", "toydata_human.gs", package = "scdRs")
#' sce <- scdrsLoadToyData()
#' gs <- scdrsReadGeneSets(f, srcSpecies = "human", dstSpecies = "mouse",
#'     intersect = rownames(sce))
#' lengths(gs)
#' out <- tempfile(fileext = ".gs")
#' scdrsWriteGeneSets(gs, out)
#' identical(scdrsReadGeneSets(out), gs)
#'
#' @export
scdrsReadGeneSets <- function(file, srcSpecies = NULL, dstSpecies = NULL,
    intersect = NULL) {
    if (is.null(srcSpecies) != is.null(dstSpecies)) {
        stop("'srcSpecies' and 'dstSpecies' must be both NULL or both given")
    }
    map <- NULL
    if (!is.null(srcSpecies) &&
        .speciesName(srcSpecies) != .speciesName(dstSpecies)) {
        map <- scdrsGetHomologs(srcSpecies, dstSpecies)
    }
    df <- utils::read.delim(file, colClasses = "character",
        check.names = FALSE, na.strings = character(0))
    if (!all(c("TRAIT", "GENESET") %in% names(df))) {
        stop("'file' must have columns TRAIT and GENESET")
    }
    out <- lapply(seq_len(nrow(df)), function(i) {
        items <- strsplit(strsplit(df$GENESET[i], ",", fixed = TRUE)[[1]],
            ":", fixed = TRUE)
        len <- lengths(items)
        if (all(len == 1L)) {
            w <- rep(1, length(items))
        } else if (all(len == 2L)) {
            w <- as.numeric(vapply(items, `[`, character(1), 2L))
        } else {
            stop("gene set ", df$TRAIT[i],
                " contains genes with invalid format")
        }
        w <- .dictDedup(stats::setNames(w,
            vapply(items, `[`, character(1), 1L)))
        if (!is.null(map)) {
            keep <- names(w) %in% names(map)
            w <- .dictDedup(stats::setNames(w[keep], map[names(w)[keep]]))
        }
        if (!is.null(intersect)) {
            w <- w[names(w) %in% intersect]
        }
        w
    })
    .dictDedup(stats::setNames(out, df$TRAIT))
}

## Python's str() of a float, for the weights written by save_gs: R's
## 15-significant-digit form, with ".0" added to whole numbers.
.pyFloat <- function(x) {
    s <- as.character(x)
    plain <- grepl("^-?[0-9]+$", s)
    s[plain] <- paste0(s[plain], ".0")
    s
}

#' @rdname scdrsReadGeneSets
#' @param geneSets A named list of gene sets: named numeric vectors (genes
#'     and weights, written as `"g:w"`) or character vectors (genes only).
#' @return `scdrsWriteGeneSets()`: `file`, invisibly.
#' @export
scdrsWriteGeneSets <- function(geneSets, file) {
    if (!is.list(geneSets) || is.null(names(geneSets))) {
        stop("'geneSets' must be a named list")
    }
    sets <- vapply(geneSets, function(s) {
        if (is.numeric(s)) {
            paste0(names(s), ":", .pyFloat(unname(s)), collapse = ",")
        } else {
            paste(s, collapse = ",")
        }
    }, character(1))
    .writeGsTable(names(geneSets), sets, file)
}

.writeGsTable <- function(traits, sets, file) {
    utils::write.table(data.frame(TRAIT = traits, GENESET = unname(sets)),
        file, sep = "\t", quote = FALSE, row.names = FALSE)
    invisible(file)
}

#' Make scDRS gene sets from gene-level GWAS statistics
#'
#' For each trait, selects the disease genes with the smallest gene-level
#' p-values (e.g. from MAGMA) and weights them by their z-scores. This is
#' the R port of the `scdrs munge-gs` command.
#'
#' The number of genes per trait is the number of discoveries at
#' `fdr` (Benjamini-Hochberg) or `fwer` (Bonferroni), or `nMax` if neither
#' is given, and is then restricted to `[nMin, nMax]`.
#'
#' @param pvals,zscores Exactly one of: gene-level p-values or one-sided
#'     z-scores, as a `data.frame` (genes as row names, one column per
#'     trait) or the path of a tab-separated file with the genes in the
#'     first column (header `GENE`). Z-scores are converted to p-values
#'     with `pnorm(-z)`. Missing values are dropped per trait.
#' @param weight `"zscore"` (default): weight
#'     \eqn{\min(-\Phi^{-1}(\max(p, 10^{-100})), 10)}{min(-qnorm(max(p,
#'     1e-100)), 10)}; or `"uniform"`: weight 1.
#' @param fdr,fwer At most one of them: the FDR or FWER level used to
#'     choose the number of genes.
#' @param nMin,nMax Minimum and maximum number of genes per trait.
#' @param file Optional path: also write the gene sets to this `.gs` file,
#'     in the format of `scdrs munge-gs` (weights with 5 significant
#'     digits).
#'
#' @return A named list (traits in sorted order) of named numeric vectors
#'     (names = genes, ordered by p-value; values = weights rounded to 5
#'     significant digits, as in the file), like [scdrsReadGeneSets()].
#'
#' @examples
#' set.seed(1)
#' p <- data.frame(BMI = runif(500)^4, HEIGHT = runif(500)^2,
#'     row.names = paste0("GENE", seq_len(500)))
#' gs <- scdrsMungeGeneSets(pvals = p, fdr = 0.05, nMin = 10, nMax = 100)
#' lengths(gs)
#' head(gs$BMI)
#'
#' @export
scdrsMungeGeneSets <- function(pvals = NULL, zscores = NULL,
    weight = c("zscore", "uniform"), fdr = NULL, fwer = NULL,
    nMin = 100L, nMax = 1000L, file = NULL) {
    weight <- match.arg(weight)
    if (is.null(pvals) == is.null(zscores)) {
        stop("exactly one of 'pvals' and 'zscores' must be given")
    }
    if (!is.null(fdr) && !is.null(fwer)) {
        stop("at most one of 'fdr' and 'fwer' is allowed")
    }
    tab <- .readGeneTable(if (is.null(pvals)) zscores else pvals)
    if (!is.null(zscores)) {
        tab[] <- lapply(tab, function(z) scdrsZToP(z))
    }
    nMin <- min(nMin, nMax)
    traits <- sort(names(tab), method = "radix")
    out <- lapply(traits, function(trait) {
        p <- stats::setNames(tab[[trait]], rownames(tab))
        p <- p[!is.na(p)]
        nGene <- if (!is.null(fdr)) {
            sum(stats::p.adjust(p, "BH") <= fdr)
        } else if (!is.null(fwer)) {
            sum(p <= fwer / length(p))
        } else {
            nMax
        }
        nGene <- max(min(nGene, nMax), nMin)
        p <- p[order(p, method = "radix")]
        p <- p[seq_len(min(nGene, length(p)))]
        w <- if (weight == "zscore") {
            pmin(scdrsPToZ(pmax(p, 1e-100)), 10)
        } else {
            rep(1, length(p))
        }
        stats::setNames(as.numeric(sprintf("%.5g", w)), names(p))
    })
    names(out) <- traits
    if (!is.null(file)) {
        sets <- vapply(out, function(s) {
            paste0(names(s), ":", sprintf("%.5g", unname(s)), collapse = ",")
        }, character(1))
        .writeGsTable(traits, sets, file)
    }
    out
}

## A gene x trait table from a data.frame (genes as row names) or a
## tab-separated file (genes in the first column). Checks that the gene
## names are unique and drops rows without a gene name.
.readGeneTable <- function(x) {
    if (is.character(x) && length(x) == 1L) {
        df <- utils::read.delim(x, check.names = FALSE,
            colClasses = c("character"), na.strings = c("", "NA", "NaN"))
        genes <- df[[1L]]
        df <- df[, -1L, drop = FALSE]
        df[] <- lapply(df, as.numeric)
    } else if (is.data.frame(x) || is(x, "DataFrame")) {
        df <- as.data.frame(x, optional = TRUE)
        genes <- rownames(df)
    } else {
        stop("gene statistics must be a data.frame or a file path")
    }
    ok <- !is.na(genes) & genes != ""
    if (!all(ok)) {
        message("Removing ", sum(!ok), " rows without a gene name")
    }
    df <- df[ok, , drop = FALSE]
    genes <- genes[ok]
    if (anyDuplicated(genes)) {
        stop("gene names are not unique: ",
            paste(unique(genes[duplicated(genes)]), collapse = ","))
    }
    if (!all(vapply(df, is.numeric, logical(1)))) {
        stop("all trait columns must be numeric")
    }
    rownames(df) <- genes
    df
}
