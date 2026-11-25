# Python reference outputs

These scripts run the original Python scDRS (version 1.0.4) to produce
the reference files in `tests/testthat/fixtures/` that the R tests
compare against. They are not used by the R package itself.

Run on the HPC / Jupyter environment:

```bash
pip install "scdrs==1.0.4"
python inst/scripts/python/<script>.py --out tests/testthat/fixtures
```

Every script fixes its random seeds and writes plain TSV (or Matrix
Market `.mtx`) files so R can read them without Python.
