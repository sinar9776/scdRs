# Python reference outputs and example data

`make_reference_outputs.py` runs the original Python scDRS (v1.0.4,
GitHub commit 644a5c7) on its bundled toy data and writes plain CSV
files. They are used for two things:

- the reference ("answer key") files in `tests/testthat/fixtures/` that
  the R tests compare against, and
- the example data in `inst/extdata/`.

The script is not used by the R package itself. Main outputs:

- `pp.<sparse|dense>.<cov|nocov>.<none|celltype>.*`: `preprocess()` in
  all 8 modes (gene/cell statistics, covariate matrix and coefficients);
  `pp.*.cov.*.nonconst.*`: the same without the 10 toy genes that are
  constant after covariate correction (`pp.constant_genes.csv`);
- `units.*`: small inputs and outputs of the internal helpers
  (`reg_out`, `category2dummy`, `_get_mean_var`, `pd.qcut`);
- `<cov|nocov|densecov>[.<option>].<trait>.*`: `score_cell()` results
  (main columns, control raw/normalized scores) and the control gene sets
  it used, with `--n-ctrl` (default 100) control sets, for the exact
  scoring tests; options `uniform`, `inv_std`, `od` (weighting),
  `matchmean` (`ctrl_match_key = "mean"`) and `weighted` (gene weights);
- `<cov|nocov>.<trait>.score1000.csv`: main columns with 1000 control
  sets, for the statistical comparison with R's own control sets;
- `units.correct_background.*`, `units.empi_null.*`, `units.pearson.*`:
  internal scoring/downstream helpers on small inputs;
- `downstream.*`: scanpy PCA and kNN graph, the `.full_score.gz` file as
  `scdrs compute-score` writes it, group (2 annotations), correlation and
  gene analyses, `gearys_c` and `test_gearysc`;
- `util.*`: `load_gs` (species mapping, duplicates), `save_gs`,
  `scdrs munge-gs` (fdr / fwer / top-n / z-score input), `test_overlap`,
  `meta_analysis`, `zsc2pval`/`pval2zsc`, and `load_h5ad` filtering and
  normalization on simulated counts.

The test fixtures are gzipped copies (`gzip -9`) of these files, except
copies that would be byte-identical to another fixture (the control sets
of the weighting options and of the human gene set equal those of the
mouse gene set; the tests reuse those).

```bash
python -m venv venv
venv/bin/pip install "git+https://github.com/martinjzhang/scDRS@644a5c7" \
    "numpy<2"
venv/bin/python -I inst/scripts/python/make_reference_outputs.py \
    --out outputs --n-ctrl 100 --seed 0
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
