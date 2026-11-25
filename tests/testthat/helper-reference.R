## Helpers for comparing scdRs output with reference output from the
## Python scDRS (v1.0.4, commit 644a5c7). Reference files are made by
## inst/scripts/python/make_reference_outputs.py (seed 0) and stored,
## gzipped, in tests/testthat/fixtures/.

referencePath <- function(...) {
    testthat::test_path("fixtures", ...)
}

## Read reference table `name` (file fixtures/<name>.csv.gz). With
## rowNames = TRUE the first column holds row names. `mean_var` bin labels
## are read as character ("3.10" must not become 3.1), Python's
## "True"/"False" become logical, and empty fields are NA. Skips the test
## if the file is missing.
readReference <- function(name, rowNames = TRUE) {
    path <- referencePath(paste0(name, ".csv.gz"))
    if (!file.exists(path)) {
        testthat::skip(paste("reference file missing:", name))
    }
    header <- names(utils::read.csv(path, nrows = 1, check.names = FALSE))
    cc <- if ("mean_var" %in% header) c(mean_var = "character") else NA
    df <- utils::read.csv(path, row.names = if (rowNames) 1 else NULL,
        check.names = FALSE, colClasses = cc, na.strings = c("NA", ""))
    for (j in seq_along(df)) {
        v <- df[[j]]
        if (is.character(v) && all(v %in% c("True", "False", NA))) {
            df[[j]] <- v == "True"
        }
    }
    df
}

## Fraction of equal elements.
agreement <- function(x, y) {
    mean(x == y)
}
