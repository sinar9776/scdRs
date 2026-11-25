## Helpers for comparing scdRs output with reference output from the
## Python scDRS. Reference files are produced by the scripts in
## inst/scripts/python/ and stored in tests/testthat/fixtures/.

referencePath <- function(...) {
    testthat::test_path("fixtures", ...)
}

## Read a reference table written by the Python scripts (TSV, first
## column = row names). Skips the test if the file is not there yet.
readReference <- function(file) {
    path <- referencePath(file)
    if (!file.exists(path)) {
        testthat::skip(paste("reference file missing:", file))
    }
    utils::read.delim(path, row.names = 1, check.names = FALSE)
}
