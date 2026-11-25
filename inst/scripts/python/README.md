# Python reference outputs and example data

`make_reference_outputs.py` runs the original Python scDRS (v1.0.4,
GitHub commit 644a5c7) on its bundled toy data and writes plain CSV
files. They are used for two things:

- the reference ("answer key") files in `tests/testthat/fixtures/` that
  the R tests compare against, and
- the example data in `inst/extdata/`.

The script is not used by the R package itself.

```bash
python -m venv venv
venv/bin/pip install "git+https://github.com/martinjzhang/scDRS@644a5c7" \
    "numpy<2"
venv/bin/python -I inst/scripts/python/make_reference_outputs.py \
    --out outputs --n-ctrl 1000 --seed 0
```

## Provenance of `inst/extdata/`

All files come from the scDRS repository (`scdrs/data/`, commit
644a5c7; MIT licence, see `inst/COPYRIGHTS`).

| File | Source |
|---|---|
| `toydata_logexpr.csv.gz` | `toydata_mouse.h5ad`, `X` (log-normalised), written genes x cells by the script, then `gzip -9` |
| `toydata_cell_metadata.csv` | `toydata_mouse.h5ad`, `obs` |
| `toydata_gene_metadata.csv` | `toydata_mouse.h5ad`, `var` |
| `toydata_covariates.csv` | `toydata_mouse.cov` |
| `toydata_mouse.gs`, `toydata_human.gs` | copied unchanged |
| `mouse_human_homologs.txt` | copied unchanged |

The toy data has 2,500 genes and 30 cells (a subset of Tabula Muris
Senis).
