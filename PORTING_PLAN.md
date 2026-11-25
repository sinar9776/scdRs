# scdRs porting plan

Port of scDRS v1.0.4 (commit `644a5c7`) to native R.

**Status (2026-11-25):** steps 1-6 of section 6 are implemented and tested
against Python; step 7 has the vignette (toy data end to end). Still open
in step 7: `R CMD check` / `BiocCheck` on the HPC, and the
statistical-equivalence run on a larger dataset. Deviations from this plan,
decided during implementation:

- No `BiocParallel` and no `BPPARAM`: all control scores of a trait are one
  sparse matrix product, so a trait takes seconds; parallel workers would
  also make `set.seed()` reproducibility harder.
- `scdrsSelectControlGenes(geneStats, genes, weights, ...)` accepts a
  preprocessed `SingleCellExperiment` or a gene-statistics table, and
  sorts / intersects / de-duplicates the genes as `score_cell()` does.
- `scdrsBuildGraph()` reimplements scanpy's exact path (dense SVD PCA,
  exact Euclidean kNN with the cell itself first, umap-learn's
  `smooth_knn_dist` + fuzzy union) instead of `uwot`; on the toy data it
  equals scanpy's graph to 1e-4. `irlba` (large PCA) and `BiocNeighbors`
  (large kNN) are Suggests.
- `scdrsTestGearysC(graph, scores, groups, method)` returns groups in
  sorted order (Python: order of appearance).
- Acceptance criterion "FDR < 0.1 calls agree for >= 99% of cells": on the
  30 toy cells one cell is 3.3%, and one cell has Python FDR 0.101; the toy
  test allows at most one disagreeing cell, which must be borderline (FDR
  in [0.05, 0.2] in both). The 99% criterion stays for the larger dataset.
- `scdrsPerformDownstream()` errors if the graph is missing (pointing to
  `scdrsBuildGraph()`) rather than building it silently.

Sources read: `scdrs/pp.py`, `scdrs/method.py`, `scdrs/util.py`,
`bin/scdrs` (the CLI), plus `scdrs/data_loader.py`, the top-level
`compute_score.py` / `compute_downstream.py`, and `tests/`.

---

## 0. Conventions used below

- **Orientation.** AnnData is cells × genes; a `SingleCellExperiment` is
  **genes × cells**. All formulas below use $X_{cg}$ (cell $c$, gene $g$),
  $n$ cells, $G$ genes; the R code will work on the transposed matrix.
- $\tilde X$ = covariate-corrected log expression (see §2.1).
- "sd₀" = population SD (divide by $n$, numpy default `ddof=0`);
  "sd₁" = sample SD (divide by $n-1$, pandas default `ddof=1`, R `sd()`).
  scDRS uses both, in different places; the R port must copy each one.
- $K$ = number of control gene sets (`n_ctrl`).

---

## 1. Normal analysis: call order

The CLI (`bin/scdrs`) defines the canonical pipeline.

```
                  ┌───────────── scdrs munge-gs ─────────────┐
 MAGMA gene       │ util.zsc2pval (if z-scores)               │
 p-values/z  ───► │ multipletests (fdr/fwer) → top genes      │ ──► .gs file
                  │ util.pval2zsc → weights                    │
                  └───────────────────────────────────────────┘

                  ┌───────────── scdrs compute-score ────────┐
 .h5ad ─────────► │ util.load_h5ad (filter, normalize, log1p) │
 .cov  ─────────► │ read covariates                           │
 .gs   ─────────► │ util.load_gs (+ load_homolog_mapping)     │
                  │ pp.preprocess                              │
                  │   ├─ category2dummy, reg_out (covariates)  │
                  │   ├─ cell weights (adj_prop)               │
                  │   └─ compute_stats (+ _get_mean_var*)      │
                  │ for each trait (≥ 10 genes):               │
                  │   method.score_cell                        │
                  │     ├─ _select_ctrl_geneset                │
                  │     ├─ _compute_raw_score (disease + K ctrl)│
                  │     │    └─ _compute_overdispersion_score   │
                  │     ├─ _correct_background                 │
                  │     └─ _get_p_from_empi_null               │
                  └───────────────────────────────────────────┘ ──► .score.gz,
                                                                     .full_score.gz
                  ┌──────────── scdrs perform-downstream ────┐
 .h5ad + full ──► │ util.load_h5ad, util.load_scdrs_score     │
 scores           │ (sc.pp.pca + sc.pp.neighbors if no graph) │
                  │ method.downstream_group_analysis          │
                  │   └─ test_gearysc → gearys_c              │
                  │ method.downstream_corr_analysis           │
                  │ method.downstream_gene_analysis           │
                  │   └─ _pearson_corr(_sparse)               │
                  └───────────────────────────────────────────┘ ──► .scdrs_group.*,
                                                                     .scdrs_cell_corr,
                                                                     .scdrs_gene
```

The equivalent R pipeline:

```r
gs  <- scdrsReadGeneSets("traits.gs", srcSpecies = "human",
                         dstSpecies = "mouse",
                         intersect = rownames(sce)) # or scdrsMungeGeneSets()
sce <- scdrsNormalize(sce, filter = TRUE)           # only if raw counts
sce <- scdrsPreprocess(sce, cov = c("n_genes", "age"))
res <- scdrsComputeScore(sce, gs, nCtrl = 1000)     # list of DataFrames
grp <- scdrsTestGroups(sce, res[["BMI"]], groupCols = "cell_type")
cor <- scdrsTestCorrelation(sce, res[["BMI"]], varCols = "effector_score")
gen <- scdrsRankGenes(sce, res[["BMI"]])
```

---

## 2. Function inventory

Each entry: what it does → inputs/outputs → math → proposed R function.
*Exported* = user-facing in R. *Internal* = not exported, but still tested
against Python (CLAUDE.md rule).

### 2.1 `pp.py`

#### `preprocess(data, cov, adj_prop, n_mean_bin=20, n_var_bin=20, n_chunk, copy)`

Prepares the data: covariate correction, then gene- and cell-level stats.

- **In:** AnnData (log-normalized), covariate DataFrame (cells × covariates),
  optional cell-group column for proportion adjustment.
- **Out (in `uns["SCDRS_PARAM"]`):** `FLAG_SPARSE`, `FLAG_COV`,
  `FLAG_ADJ_PROP`, `COV_MAT` (n × k), `COV_BETA` (G × k, **negated**),
  `COV_GENE_MEAN` (G), `GENE_STATS` (G × 7), `CELL_STATS` (n × 2).
  In dense mode `X` is overwritten by $\tilde X$; in sparse+cov mode
  ("implicit covariate correction") `X` is left sparse and $\tilde X$ is
  only ever formed in chunks.

**Math.**

1. Covariates: join to cells (cells absent from `cov` become NA; at least
   75 % must match), dummy-code non-numeric columns (`category2dummy`),
   replace NA by column mean → $C$ ($n \times k$).
2. Constant term: regress $\mathbf 1$ on $C$; if
   $\frac1n\lVert \mathbf 1 - C(C^\top C)^{-1}C^\top \mathbf 1\rVert^2 > 0.01$,
   append a column `SCDRS_CONST = 1`.
3. Gene means $\mu_g = \frac1n\sum_c X_{cg}$.
4. $\beta = (C^\top C / n)^{-1}(C^\top X / n)$ ($k \times G$).
5. Corrected data: $\tilde X = X - C\beta + \mathbf 1\mu^\top$
   (stored as `COV_BETA` $= -\beta^\top$, so
   $\tilde X = X + C\,\mathrm{COV\_BETA}^\top + \mu$).
6. Cell weights (only if `adj_prop`): group sizes $m_j$ clipped below at
   $\lfloor 0.01 \max_j m_j\rfloor$; $w_c = n / m_{j(c)}$; then
   $w \leftarrow w/\bar w$. Requires fewer than $0.1n$ groups.
7. `n_chunk` default: $5nG/\mathrm{nnz} + 1$ in implicit mode, else 20.
8. Calls `compute_stats`.

**R:** *exported*

```r
scdrsPreprocess(sce, assay = "logcounts", cov = NULL, adjProp = NULL,
                nMeanBin = 20L, nVarBin = 20L, nChunk = NULL,
                implicit = NULL)
```

- `cov`: character vector of `colData` column names, **or** a
  `DataFrame`/`data.frame` with rownames = cell names (mirrors the `.cov`
  file).
- `adjProp`: one `colData` column name.
- `implicit`: `NULL` = follow Python's rule (sparse + cov → implicit);
  `TRUE`/`FALSE` forces a mode.
- Returns `sce` with results stored **subset-safely**:

  | Python | R location |
  |---|---|
  | `GENE_STATS` | `rowData(sce)` columns `scdrs_mean`, `scdrs_var`, `scdrs_var_tech`, `scdrs_ct_mean`, `scdrs_ct_var`, `scdrs_ct_var_tech`, `scdrs_mean_var` |
  | `COV_BETA`, `COV_GENE_MEAN` | `rowData(sce)$scdrs_cov_beta` (nested G × k matrix), `rowData(sce)$scdrs_cov_gene_mean` |
  | `CELL_STATS` | `colData(sce)` columns `scdrs_mean`, `scdrs_var` |
  | `COV_MAT` | `reducedDim(sce, "scdrs_cov")` (n × k) |
  | flags, parameters | `metadata(sce)$scdrs` |
  | dense $\tilde X$ | **new** assay `"scdrs_corrected"` (the input assay is never overwritten) |

  Storing per-gene and per-cell values in `rowData`/`colData`/`reducedDim`
  means they stay correct when the user subsets `sce`, which Python
  achieves by indexing by name.

#### `category2dummy(df, cols, verbose)`

Turns non-numeric columns into 0/1 dummies, dropping the first level
(`pd.get_dummies(drop_first=True)`); NA rows stay NA in all dummy columns.
Column names are `<col>_<level>`.

**R:** *internal* `.categoryToDummy(df, cols = NULL)`. Levels are sorted
with `sort(method = "radix")` (Python order, see §4.6). Uses
`model.matrix`-free code, so NA handling matches.

#### `compute_stats(adata, implicit_cov_corr, cell_weight, n_mean_bin, n_var_bin, n_chunk)`

Gene and cell statistics on $\tilde X$.

**Math.** With weights $w$ (all 1 if none):

- $\mathrm{mean}_g = \sum_c w_c \tilde X_{cg} / \sum_c w_c$
- $\mathrm{var}_g = \sum_c w_c \tilde X_{cg}^2 / \sum_c w_c - \mathrm{mean}_g^2$
  (computed as $E[x^2]-E[x]^2$, **not** a two-pass variance)
- `ct_mean`, `ct_var`: same two formulas on $e^{\tilde X}-1$.
- Technical variance (Seurat v3 recipe): on genes with `ct_var > 0` (and
  `ct_mean > 0` if any `ct_mean ≤ 0`, with a warning), fit LOESS
  $\log_{10}\mathrm{ct\_var} \sim \log_{10}\mathrm{ct\_mean}$, `span = 0.3`,
  `degree = 2`. $\mathrm{ct\_var\_tech}_g = 10^{\hat y_g}$; excluded genes
  get $10^0 = 1$.
- $\mathrm{var\_tech}_g = \mathrm{var}_g \cdot \mathrm{ct\_var\_tech}_g /
  \mathrm{ct\_var}_g$; NaN (0/0) → 0.
- Mean–variance bins: $b_{\max} = \lfloor\sqrt{G/10}\rfloor$; if either bin
  count exceeds it, **both** become $b_{\max}$. `qcut(mean, n_mean_bin)`
  (quantile bins, duplicate edges dropped, 0-based labels), then within each
  mean bin `qcut(var, n_var_bin)`; label `"<i>.<j>"`.
- Cell stats: **unweighted** mean and var of $\tilde X_{c\cdot}$ over genes.

**R:** *exported* (advanced) `scdrsComputeGeneStats(sce, assay, cellWeight = NULL,
nMeanBin = 20L, nVarBin = 20L, nChunk = 20L)`, called by
`scdrsPreprocess()`. Returns a list of two `DataFrame`s (`gene`, `cell`).

#### `reg_out(mat_Y, mat_X)`

Residuals of a least-squares regression:
$Y - X(X^\top X/n)^{-1}(X^\top Y/n)$.

**R:** *internal* `.regOut(y, x)`, using `solve(crossprod(x)/n,
crossprod(x, y)/n)` (same normal-equations route as numpy, not `qr()`, so
the results match).

#### `_get_mean_var(X, axis, weights)`, `_weighted_sparse_average(X, weights, axis)`

Row or column mean and $E[x^2]-E[x]^2$ variance, weighted or unweighted,
dense or sparse. The sparse weighted mean simplifies to
$\sum w x / \sum w$.

**R:** *internal* `.meanVar(x, margin = c("gene", "cell"), weights = NULL)`
using `MatrixGenerics::rowMeans2/colMeans2` (works for base matrices and
`dgCMatrix` via `sparseMatrixStats`); weighted version via
`x %*% w / sum(w)`.

#### `_get_mean_var_implicit_cov_corr(adata, axis, weights, transform_func, n_chunk)`

Same statistics for $f(\tilde X)$ without materializing $\tilde X$: loops over
chunks of genes (axis 0) or cells (axis 1), builds the dense chunk
$X_{\cdot S} + C\,\beta^\top_{S} + \mu_S$, applies $f$ (identity or
`expm1`), calls `_get_mean_var`.

**R:** *internal* `.meanVarImplicit(x, covMat, covBeta, geneMean, margin,
weights = NULL, fun = identity, nChunk = 20L)`.
Chunk size `max(1, n %/% nChunk)`: Python's `n // n_chunk` can be 0 and
loop forever; R guards this (§4.11).

### 2.2 `method.py`: scoring

#### `score_cell(data, gene_list, gene_weight, ctrl_match_key="mean_var", n_ctrl=1000, n_genebin=200, weight_opt="vs", return_ctrl_raw_score, return_ctrl_norm_score, random_seed=0, ...)`

The core scDRS score.

- **In:** preprocessed data, disease genes, optional gene weights.
- **Out:** per-cell table (float32) with `raw_score`, `norm_score`,
  `mc_pval`, `pval`, `nlog10_pval`, `zscore`, optional
  `ctrl_raw_score_0..K-1`, `ctrl_norm_score_0..K-1`.

**Steps / math.**

1. Keep disease genes present in the data; **sort** them; weights $a_g$
   (default 1). Duplicate gene names are dropped.
2. Control sets: `_select_ctrl_geneset`.
3. Raw disease score and $K$ raw control scores: `_compute_raw_score`.
4. Variance ratio (only if `ctrl_match_key == "mean_var"` and weighting is
   not `od`; else all 1):
   $r_k = \dfrac{\sum_{g \in \mathcal C_k} \mathrm{var}_g\, w_{g,k}^2}
   {\sum_{g \in \mathcal D} \mathrm{var}_g\, w_g^2}$,
   where $w$ are the normalized score weights returned by
   `_compute_raw_score`.
5. `_correct_background` → normalized scores $z_c$, $Z_{ck}$.
6. Per-cell Monte Carlo p-value:
   $p^{\mathrm{mc}}_c = \dfrac{1 + \#\{k: Z_{ck} \ge z_c\}}{1+K}$.
7. Pooled p-value (main scDRS p-value), using the null pooled over all
   cells and controls:
   $p_c = \dfrac{1 + \#\{(c',k): Z_{c'k} \ge z_c\}}{1 + nK}$.
8. `nlog10_pval` $= -\log_{10} p_c$;
   `zscore` $= \mathrm{clip}(-\Phi^{-1}(p_c), -10, 10)$.

**R:** *exported*

```r
scdrsScoreCells(sce, genes, weights = NULL, ctrlMatchKey = "mean_var",
           nCtrl = 1000L, nGeneBin = 200L,
           weightOpt = c("vs", "uniform", "inv_std", "od"),
           returnCtrlRawScore = FALSE, returnCtrlNormScore = FALSE,
           controlSets = NULL, BPPARAM = BiocParallel::SerialParam())
```

- `genes`: character vector, **or** a named numeric vector (names = genes,
  values = weights; the format `scdrsReadGeneSets()` returns).
- `controlSets`: optional precomputed control sets (list of $K$ character
  vectors, or the output of `scdrsSelectControlGenes()`). If given, no random
  draws happen. This is how the R tests reproduce Python exactly (§4.1),
  and it lets users reuse control sets.
- Randomness follows R conventions: no `seed` argument. Users call
  `set.seed()` first; the vignette and examples do.
- Returns a `DataFrame` (one row per cell, rownames = cell names) with
  double-precision columns `raw_score`, `norm_score`, `mc_pval`, `pval`,
  `nlog10_pval`, `zscore`. Control scores, if requested, are stored as
  **nested matrix columns** `ctrl_raw_score` and `ctrl_norm_score`
  ($n \times K$), not as $K$ separate columns.
- Speed: all $K$ control raw scores are one sparse matrix product
  $\tilde X W$ with $W$ a $G \times K$ sparse weight matrix (a gene appears at
  most once per control set). This replaces the Python loop of $K$
  mat-vec products. For implicit mode:
  $\tilde X W = XW + C(\beta W) + \mathbf 1(\mu^\top W)$.

#### `_select_ctrl_geneset(df_gene, gene_list, gene_weight, ctrl_match_key, n_ctrl, n_genebin, random_seed)`

- If `ctrl_match_key` has fewer than $G/10$ distinct values it is treated
  as categorical (bins = its values, e.g. `mean_var`); otherwise
  `qcut(·, n_genebin)`.
- For each bin $b$ (in sorted bin order) that contains $d_b > 0$ disease
  genes, and for each of the $K$ sets: draw $d_b$ genes **uniformly
  without replacement from all genes in bin $b$, disease genes included**.
  The control genes inherit the weights of that bin's disease genes (in
  sorted gene order).

**R:** *exported* (advanced)
`scdrsSelectControlGenes(geneStats, genes, weights, ctrlMatchKey = "mean_var",
nCtrl = 1000L, nGeneBin = 200L)`. Returns a list of $K$ named numeric
vectors (names = control genes, values = weights).

#### `_compute_raw_score(adata, gene_list, gene_weight, weight_opt)`

Weighted mean expression over a gene set.

- Per-gene factor $u_g$: `uniform` 1; `vs`
  $1/\sqrt{\mathrm{var\_tech}_g + 0.01}$; `inv_std`
  $1/\sqrt{\mathrm{var}_g + 0.01}$.
- $w_g = u_g a_g / \sum_{g'} u_{g'} a_{g'}$.
- $s_c = \sum_{g} w_g \tilde X_{cg}$ (implicit mode: computed as
  $X w + C(\beta w) + \mu^\top w$).
- `od`: delegates to `_compute_overdispersion_score`.
- Returns $(s, w)$.

**R:** *internal* `.rawScore(sce, genes, weights, weightOpt)`, plus a
matrix version `.rawScoreMatrix()` used for all controls at once.

#### `_compute_overdispersion_score(adata, gene_list, gene_weight)`

Benchmarking option `od`.
$w_g \propto a_g / (\mathrm{var\_tech}_g + 0.01)$, normalized to sum 1;
$s_c = \sum_g w_g\left[(\tilde X_{cg} - \mathrm{mean}_g)^2 -
\mathrm{var\_tech}_g\right]$, computed as a quadratic + linear + constant
term. Returns weights = all ones.

**R:** *internal* `.overdispersionScore(...)`.

#### `_correct_background(v_raw_score, mat_ctrl_raw_score, v_var_ratio_c2t)`

1. Remember which raw scores are exactly 0.
2. First gene-set alignment: $s \leftarrow s - \bar s$;
   $S_{\cdot k} \leftarrow (S_{\cdot k} - \bar S_{\cdot k}) / \sqrt{r_k}$.
3. Cell-wise standardization: $m_c = \mathrm{mean}_k S_{ck}$,
   $\sigma_c = \mathrm{sd}_{0,k}(S_{ck})$;
   $z_c = (s_c - m_c)/\sigma_c$, $Z_{ck} = (S_{ck} - m_c)/\sigma_c$.
4. Second gene-set alignment: $z \leftarrow z - \bar z$;
   $Z_{\cdot k} \leftarrow Z_{\cdot k} - \bar Z_{\cdot k}$.
5. Zero raw scores: let $\min = \min(\min z, \min Z)$ (taken **before**
   the replacement); $z_c = \min - 10^{-3}$ where $s_c$ was 0;
   $Z_{ck} = \min$ where $S_{ck}$ was 0.

**R:** *internal* `.correctBackground(rawScore, ctrlRawScore, varRatio)`;
the optional `save_intermediate` stages become a `returnIntermediate =
FALSE` argument, used by the tests to compare each stage.

#### `_get_p_from_empi_null(v_t, v_t_null)`

$p_i = \dfrac{1 + \#\{j: T^{\mathrm{null}}_j \ge t_i\}}{1+N}$, in
$O(N\log N)$ via sort + `searchsorted(side="left")`.

**R:** *internal* `.empiricalPvalue(t, tNull)`:
`(N - findInterval(t, sort(tNull), left.open = TRUE) + 1) / (N + 1)`.

### 2.3 `method.py`: downstream

#### `downstream_group_analysis(adata, df_full_score, group_cols, fdr_thresholds=[0.05, 0.1, 0.2])`

Per annotation column and per group: number of significant cells, group
association, and within-group heterogeneity.

- Cells = sorted intersection of data and score table.
- $q_c$ = BH-adjusted $p_c$ over **all** cells (once per column).
- For each group $G_j$:
  - `n_cell`, `n_ctrl` $=K$, `n_fdr_<t>` $= \#\{c \in G_j: q_c < t\}$.
  - Association: $Q = \mathrm{quantile}_{0.95}(z_{G_j})$,
    $Q_k = \mathrm{quantile}_{0.95}(Z_{G_j,k})$ (linear interpolation);
    `assoc_mcp` $= (1 + \#\{k: Q_k \ge Q\})/(K+1)$,
    `assoc_mcz` $= (Q - \bar Q_k)/\mathrm{sd}_0(Q_k)$.
  - Heterogeneity: `test_gearysc` →
    `hetero_mcp`, `hetero_mcz`.
- Needs `adata.obsp["connectivities"]` (the scanpy kNN graph).

**R:** *exported*

```r
scdrsTestGroups(sce, scores, groupCols,
                fdrThresholds = c(0.05, 0.1, 0.2),
                graph = "connectivities")
```

- `scores`: the `DataFrame` from `scdrsScoreCells()` (must contain
  `ctrl_norm_score`), or a data.frame read from a Python `.full_score.gz`
  (columns `ctrl_norm_score_*` are recognized).
- `graph`: name of a `colPair(sce)` entry, or an $n \times n$ sparse matrix.
- Returns a named `List` of `DataFrame`s, one per `groupCols` entry, with
  Python's column names (`n_cell`, `n_ctrl`, `assoc_mcp`, `assoc_mcz`,
  `hetero_mcp`, `hetero_mcz`, `n_fdr_0.05`, …).

#### `test_gearysc(adata, df_full_score, groupby, subset, opt="control_distribution_match")`

For each group: Geary's C of the disease score on the group's induced
subgraph, versus $K$ null statistics.

- Default option `control_distribution_match`: null $k$ uses the values of
  $z$ rearranged to follow the order of $Z_{\cdot k}$:
  $\mathrm{null}_k = \mathrm{sort}(z_{G})[\mathrm{rank}_{\text{ordinal}}(Z_{G,k})]$.
  Other options: `permutation` (random permutations of $z$) and `control`
  (raw $Z_{\cdot k}$).
- `pval` $= (1 + \#\{k: C_{\mathrm{trait}} > C_k\})/(K+1)$ (strict `>`);
  NaN if $C_{\mathrm{trait}}$ is NaN.
- `zsc` $= -(C_{\mathrm{trait}} - \bar C_k)/\mathrm{sd}_1(C_k)$
  (**sample** SD here, unlike `assoc_mcz`).
- Out: per group `pval`, `trait`, `ctrl_mean`, `ctrl_std`, `zsc`.

**R:** *exported* (advanced)
`scdrsTestGearysC(graph, scores, groups, method = c("control_distribution_match",
"permutation", "control"))`.

#### `gearys_c(adata, vals)`

$C = \dfrac{(N-1)\sum_{i,j} w_{ij}(x_i - x_j)^2}{2W\sum_i (x_i - \bar x)^2}$,
$W = \sum_{ij} w_{ij}$, vectorized over columns of `vals`.

**R:** *exported* `scdrsComputeGearysC(graph, x)`. Instead of Python's per-row loop,
vectorized with Matrix:
$\sum_{ij} w_{ij}(x_i-x_j)^2 = \sum_i r_i x_i^2 + \sum_j k_j x_j^2 - 2x^\top W x$
($r$ = row sums, $k$ = column sums), valid for asymmetric $W$ too, for all
columns at once.

#### `downstream_corr_analysis(adata, df_full_score, var_cols)`

For each continuous cell variable $v$:
$\rho = \mathrm{cor}(v, z)$, $\rho_k = \mathrm{cor}(v, Z_{\cdot k})$;
`corr_mcp` $= (1 + \#\{k: \rho_k \ge \rho\})/(K+1)$,
`corr_mcz` $= (\rho - \bar\rho_k)/\mathrm{sd}_0(\rho_k)$; plus `n_ctrl`.

**R:** *exported* `scdrsTestCorrelation(sce, scores, varCols)` → `DataFrame`
(one row per variable).

#### `downstream_gene_analysis(adata, df_full_score)`

Pearson correlation of each gene's expression with $z$ across cells; table
sorted by `CORR` descending, `RANK` = 0-based position. It uses `adata.X`
as loaded (log-normalized, **not** covariate-corrected; the CLI reloads the
data without running `preprocess`).

**R:** *exported* `scdrsRankGenes(sce, scores, assay = "logcounts")`.
`RANK` is 1-based in R (R convention; the tests compare `RANK - 1`).

#### `_pearson_corr(X, Y)`, `_pearson_corr_sparse(X, Y)`

Column-wise Pearson correlation with population SD clipped at $10^{-8}$;
the sparse version uses $\frac1n X^\top Y - \bar x\bar y^\top$ to keep $X$
sparse. Output cast to float32.

**R:** *internal* `.pearsonCor(x, y)`, sparse-aware
(`crossprod(x, y)/n - outer(mx, my)`).

#### Not ported (stated in the vignette)

- `score_cell_vision`, `score_cell_scanpy`: comparison methods from the
  paper's benchmarks. `scdrsScoreVision()` could be added later (about 10
  lines); `score_cell_scanpy` wraps scanpy's `score_genes`.
- `correlate_gene`, `_spearman_corr`, `_get_rank`: marked "Archived" in the
  source.

### 2.4 `util.py`

| Python | What / math | R proposal |
|---|---|---|
| `load_h5ad(file, flag_filter_data, flag_raw_count, min_genes=250, min_cells=50)` | Read .h5ad; error on NaN or negative values; optional `filter_cells(min_genes)` **then** `filter_genes(min_cells)`; optional `normalize_per_cell(counts_per_cell_after=1e4)` + natural `log1p`: $X_{cg} \leftarrow \log(1 + 10^4 X_{cg}/\sum_{g'}X_{cg'})$ | *exported* `scdrsNormalize(sce, assay = "counts", filter = TRUE, minGenes = 250L, minCells = 50L, normalize = TRUE)` → adds `"logcounts"`. Reading .h5ad is left to `zellkonverter::readH5AD(reader = "R")` (pure-R reader; Suggests only). **Not** `scuttle::logNormCounts` (log2 + centred size factors give different numbers). |
| `load_scdrs_score(score_file, obs_names)` | Read one or many (`@` wildcard) `.full_score.gz`; skip files covering < 10 % of cells | *exported* `scdrsReadScores(file, cells = NULL)` → named list of `DataFrame`s with nested `ctrl_*` matrices; reads Python output too |
| `load_gs(gs_path, src_species, dst_species, to_intersect)` | Parse `TRAIT\tGENESET` with `g1,g2` or `g1:w1,g2:w2`; optional mouse↔human mapping (unmapped genes dropped); optional intersect | *exported* `scdrsReadGeneSets(file, srcSpecies = NULL, dstSpecies = NULL, intersect = NULL)` → named list of named numeric vectors |
| `save_gs(gs_path, dict_gs)` | Write `.gs` | *exported* `scdrsWriteGeneSets(geneSets, file)` |
| `load_homolog_mapping(src, dst)` | Mouse↔human symbol map from the bundled `mouse_human_homologs.txt` (16,466 pairs) | *exported* `scdrsGetHomologs(src, dst)`; table shipped in `inst/extdata/` (MIT, attributed) |
| `convert_species_name(species)` | Normalize "mouse"/"Mus_musculus"/… → `mmusculus`/`hsapiens` | *internal* `.speciesName()` |
| `str_or_list_like(x)` | CLI argument helper | not needed |
| `test_overlap(l1, l2, bg)` | Fisher's exact test on 2×2 table; sample odds ratio with 95 % Wald CI on log OR (0 bounds if any cell is 0) | *exported* `scdrsTestOverlap(set1, set2, background)` (see §4.12) |
| `meta_analysis(effects, se, method, weights)` | Fixed / DerSimonian–Laird random-effects: $\tau^2 = \max\left(0, \frac{Q - (m-1)}{\sum v - \sum v^2/\sum v}\right)$, $v = 1/\mathrm{se}^2$ | *exported* `scdrsMetaAnalyze(effects, se, method = c("random", "fixed"), weights = NULL)` |
| `zsc2pval`, `pval2zsc` | $p = \Phi(-z)$; $z = -\Phi^{-1}(p)$ | *exported* `scdrsZToP()`, `scdrsPToZ()` (`pnorm(-z)`, `-qnorm(p)`, accurate in the tails) |
| `plot_group_stats`, `plot_heatmap`, `discrete_cmap`, `small_squares` | Heatmap: proportion of FDR<0.1 cells per trait × group; squares = significant association (BH across all trait × group), × = significant heterogeneity | *exported* `scdrsPlotGroupStats(groupResults, assocFdr = 0.05, heteroFdr = 0.05)` with ggplot2 (Suggests) — **phase 2** |
| `plot_qq`, `qqplot`, `plot_score_umap`, `plot_assoc_matrix` | Plotting conveniences | phase 2, optional (`scdrsPlotQQ()`, `scdrsPlotUmap()`) |
| `p_2_str`, `p_2_str_num`, `reorder_col` | Formatting / clustering helpers for figures | not ported |

### 2.5 CLI (`bin/scdrs`)

Bioconductor packages don't need a command-line tool; each CLI subcommand
becomes a high-level R wrapper with the same defaults.

| CLI | R wrapper | Notes |
|---|---|---|
| `munge-gs` | `scdrsMungeGeneSets(pvals = NULL, zscores = NULL, weight = c("zscore", "uniform"), fdr = NULL, fwer = NULL, nMin = 100L, nMax = 1000L)` | Input: data.frame with gene rownames, one column per trait. Per trait: drop NA; $n$ = #BH rejections at `fdr` / #Bonferroni rejections at `fwer` / `nMax`; clamp to [`nMin`, `nMax`]; take the $n$ smallest p; weights $\min(-\Phi^{-1}(\max(p, 10^{-100})), 10)$ or 1; traits sorted. `.gs` weights are written with `%.5g`. |
| `compute-score` | `scdrsComputeScore(sce, geneSets, cov = NULL, adjProp = NULL, ctrlMatchKey = "mean_var", weightOpt = "vs", nCtrl = 1000L, returnCtrlRawScore = FALSE, returnCtrlNormScore = TRUE, minGenes = 10L, BPPARAM = SerialParam())` | Runs `scdrsPreprocess()` if not done; skips traits with < 10 genes (message); returns named `List` of score `DataFrame`s; `scdrsWriteScores()` writes Python-compatible `.score.gz` / `.full_score.gz`. |
| `perform-downstream` | `scdrsPerformDownstream(sce, scores, groupCols = NULL, varCols = NULL, geneAnalysis = FALSE, graph = "connectivities")` | Calls the three downstream functions for each trait. If the graph is missing: error with a message pointing to `scdrsBuildGraph()` (below). |

The legacy argparse scripts `compute_score.py` / `compute_downstream.py`
are superseded by the CLI and are not ported. `data_loader.py` holds
loaders for specific datasets from the paper and is not ported.

**kNN graph.** Python builds it with `sc.pp.pca(n_comps=20)` +
`sc.pp.neighbors(n_neighbors=15)` (UMAP fuzzy-simplicial-set
connectivities). R: *exported* `scdrsBuildGraph(sce, assay = "logcounts",
nNeighbors = 15L, nPCs = 20L)` using `BiocSingular` PCA +
`uwot::similarity_graph()` (same UMAP membership-strength formula), stored
in `colPair(sce, "connectivities")`. Users can also supply any graph
(e.g. from `bluster`/`scran`). The result is close to scanpy's but not
identical (§4.9).

---

## 3. Dependencies: Python → R / Bioconductor

| Python | Used for | R / Bioconductor | Field |
|---|---|---|---|
| anndata | data container, .h5ad | `SingleCellExperiment`, `SummarizedExperiment`, `S4Vectors` | Imports |
| anndata `read_h5ad` | reading .h5ad | `zellkonverter::readH5AD(reader = "R")` (pure R via rhdf5) | Suggests |
| numpy | arrays, linear algebra, sort/searchsorted | base R (`solve`, `crossprod`, `findInterval`) | base |
| scipy.sparse | sparse CSR/CSC matrices | `Matrix` (`dgCMatrix`) | Imports |
| numpy row/col reductions | means, variances | `MatrixGenerics` + `sparseMatrixStats` (`rowMeans2`, `colMeans2`, …); `DelayedMatrixStats` for `DelayedArray` assays | Imports / Suggests |
| — | large on-disk data | `DelayedArray` / `HDF5Array` via block processing (`DelayedArray::blockApply`) | Suggests, phase 2 |
| pandas | data frames, `qcut`, `groupby`, `get_dummies` | `S4Vectors::DataFrame`, base `cut`/`quantile(type = 7)`, `split`, own dummy coder | Imports |
| scipy.stats `norm` | `cdf`, `ppf` | `stats::pnorm`, `stats::qnorm` | base |
| scipy.stats `rankdata(method="ordinal")` | distribution matching | `rank(ties.method = "first")` | base |
| scipy.stats `fisher_exact` | `test_overlap` | `stats::fisher.test` (p-value) + own sample odds ratio | base |
| skmisc.loess | technical variance trend | `stats::loess` (same netlib `dloess` Fortran code underneath) | base |
| statsmodels `multipletests` | BH / Bonferroni | `stats::p.adjust("BH" / "bonferroni")` | base |
| scanpy `normalize_per_cell`, `log1p` | normalization | own 3-line implementation (`Matrix` ops) | — |
| scanpy `filter_cells`, `filter_genes` | QC filter | own implementation (`colSums(x > 0)`, `rowSums(x > 0)`) | — |
| scanpy `pca`, `neighbors` | kNN graph for Geary's C | `BiocSingular::runPCA` (or `scater::runPCA`) + `uwot::similarity_graph`; alternatives `BiocNeighbors` + `bluster` | Suggests |
| tqdm | progress bar | `BiocParallel` progress bar / `message()` | Imports |
| (none; loops) | parallelism over traits | `BiocParallel` (`bplapply`, `BPPARAM`) | Imports |
| fire | CLI | not needed | — |
| matplotlib, seaborn | plots | `ggplot2` | Suggests, phase 2 |
| pytest | tests | `testthat` (3rd edition) | Suggests |

`Imports` stays small: `SingleCellExperiment`, `SummarizedExperiment`,
`S4Vectors`, `Matrix`, `MatrixGenerics`, `sparseMatrixStats`, `BiocParallel`,
`methods`, `stats`, `utils`. Rcpp only if profiling shows a need
(candidates: the implicit-mode chunk loop, Geary's C on very large groups).

---

## 4. What will be hard to reproduce exactly, and the plan for each

### 4.1 Random control gene sets (the big one)

Python draws control genes with `np.random.seed(0)` +
`np.random.choice(bin_genes, size, replace=False)` (numpy legacy MT19937 and
a full Fisher–Yates permutation per call). R's `set.seed()` + `sample()` uses
a different seeding scramble and a different sampling algorithm, so **R will
never draw the same control genes as Python**, and the final scores will
differ by Monte Carlo noise.

**Plan:**

1. **Exact tests through injection.** The Python reference scripts export
   the control sets they used (long TSV: `ctrl_id, gene, weight`).
   `scdrsScoreCells(..., controlSets = )` accepts them, so every deterministic
   step after sampling is compared at tight tolerance (`raw_score`, control
   raw scores, the four `_correct_background` stages, `norm_score`,
   `mc_pval`, `pval`, `zscore`).
2. **Algorithm tests for `scdrsSelectControlGenes()`** (no Python parity
   possible): every control gene lies in the same bin as the disease gene
   it replaces; each bin contributes the same number of genes; weights are
   copied correctly; genes are distinct within a set; frequencies are
   uniform within bins (χ² test over many draws).
3. **Statistical equivalence tests** end to end with R's own random draws
   (toy data and one larger dataset on the HPC), using the acceptance
   criteria in §7.4, plus a KS test that p-value distributions agree
   under the null.
4. Bit-identical RNG (reimplementing numpy's legacy MT19937 + `choice` in
   Rcpp) is **skipped** as long as the acceptance criteria in §7.4 are
   met.

### 4.2 Quantile bin boundaries (`mean_var`)

`pd.qcut` puts genes into bins at quantile cut points. Tiny floating-point
differences in `mean`/`var` (summation order, BLAS) can move a gene that
sits exactly on a cut point into the neighbouring bin, which then changes
which control genes can be drawn. Python's own tests have the same problem:
they overwrite `mean_var` with a stored reference "to avoid numerical
instability" and compare final scores at `rtol = 1e-2`.

**Plan:** reproduce `qcut` exactly (`quantile(type = 7)` edges, drop
duplicate edges, right-closed bins with the lowest edge included, 0-based
consecutive labels). Tests: (a) numeric gene stats at tight tolerance;
(b) bin labels agree for ≥ 99.5 % of genes, and every mismatch is a gene
within 1e-10 of a cut point; (c) downstream tests inject Python's
`mean_var` labels, as upstream does.

### 4.3 LOESS (`skmisc.loess` vs `stats::loess`)

Both wrap Cleveland's netlib `dloess`, so they should agree to about 1e-6,
**if** the options match: `span = 0.3`, `degree = 2`, `family = "gaussian"`,
`surface = "interpolate"`, `statistics = "approximate"`, `cell = 0.2`,
`normalize` (irrelevant with one predictor).

**Plan:** test `ct_var_tech` at `tolerance = 1e-6`. If it fails, compare
with `surface = "direct"` on both sides to check whether the k-d tree
interpolation is the cause, and record the result.

**Result (preprocessing port):** given the same inputs, `stats::loess`
and `skmisc.loess` agree to ~1e-13. One data quirk: genes that are an
exact linear function of the covariates are constant after correction,
and their `ct_var` is floating-point noise (~1e-28, different in R and
Python) that still enters the LOESS fit in Python. R copies this rule
(no special-casing), so on such data `var_tech` can differ by up to ~20%
for genes with similar means. The toy data has 10 such genes; the tests
compare those scenarios loosely and, strictly (1e-6), the same scenarios
with the 10 genes removed (`pp.*.nonconst`).

### 4.4 Variance formula and SD conventions

- Variances are $E[x^2] - E[x]^2$, which can be slightly negative or lose
  precision for nearly constant genes. R must use the **same** formula, not
  `var()`, so the numbers match (including the edge cases).
- SDs: `_correct_background`, `assoc_mcz` and `corr_mcz` use sd₀; the
  Geary `zsc` uses sd₁ (pandas). Each gets a dedicated test.

### 4.5 Float32 output

Python casts all `score_cell` and downstream results to float32 (about 7
significant digits). R keeps doubles.

**Plan:** compare with `tolerance = 1e-5` (relative) on reference values;
p-values are compared on the `nlog10_pval` scale as well.

### 4.6 String ordering

Python sorts strings by Unicode code point; R's `sort()` follows the
locale (e.g. `"a" < "B"` in `en_US`, but `"B" < "a"` in Python). It affects:
the sorted disease-gene list (weight order), bin iteration order, dummy-coded
covariate levels, `sorted(set(group))` in the downstream output, and sorted
cell lists.

**Plan:** use `sort(method = "radix")` (locale-independent, C ordering) and
`order(method = "radix")` everywhere Python calls `sorted`, plus a test with
mixed-case names.

### 4.7 Covariate handling

- Dummy coding must match `pd.get_dummies(drop_first = True)` (first level
  in sorted order dropped, NA propagated, then NA → column mean).
- The "is there already a constant?" check (mean residual² > 0.01) must be
  copied, not replaced by "always add an intercept", because covariates
  such as a complete set of dummies already span the constant.
- `np.linalg.solve` and R `solve()` both use LAPACK `gesv`; differences of
  about 1e-12 are expected. Near-collinear covariates may differ more;
  tested with a collinear example.

### 4.8 Normalization

`normalize_per_cell(counts_per_cell_after=1e4)` + `log1p` uses the natural
log and a fixed 10,000 target. Bioconductor's usual `logNormCounts` uses
log2 and centred size factors, which gives different values.
`scdrsNormalize()` copies scanpy exactly, and the order of filtering (cells,
then genes) is preserved. `normalize_per_cell` is deprecated in current
scanpy; the reference scripts pin the scanpy version used with
scdrs 1.0.4.

### 4.9 kNN graph for Geary's C

scanpy's graph depends on PCA (ARPACK), the approximate neighbour search
(pynndescent for larger data), and UMAP's smoothing. An R graph will be
similar, not identical, so `hetero_*` results will differ slightly.

**Plan:** the reference scripts export `obsp["connectivities"]` as Matrix
Market; tests put it into `colPair(sce)` so Geary's C and `hetero_*` are
compared exactly. `scdrsBuildGraph()` is tested separately for neighbour
overlap with scanpy (e.g. ≥ 90 % shared neighbours), not equality. Upstream
itself only reaches `rtol = 0.05` for `hetero_mcz`.

### 4.10 Ties

- `mc_pval`/`pval` use `≥`; Geary `pval` uses strict `>`; copied
  exactly.
- `rankdata(method = "ordinal")` breaks ties by position →
  `rank(ties.method = "first")`.
- Exact zero raw scores (`== 0`) are rare after covariate correction but
  common without it; tested on cells with no expression of the gene set.

### 4.11 Python edge cases and quirks

Copied unless they are plain bugs; each one is documented in the code and
tested:

- `chunk_size = n // n_chunk` can be 0 → infinite loop in Python. R uses
  `max(1, ·)` (bug fix, no change to results otherwise).
- Control genes may include disease genes (sampling from the whole bin).
  Copied: this is the published method.
- Genes excluded from the LOESS fit get `ct_var_tech = 1`, not NA. Copied.
- Duplicate genes in a `.gs` line: the last weight wins, the first position
  is kept. Copied.
- `downstream_group_analysis` lists groups from all cells in `adata`, not
  just scored cells; a group with no scored cells fails in Python. R gives
  NA for that group with a warning.
- Duplicate gene names: Python silently drops duplicates. R requires
  unique `rownames(sce)` and errors with a helpful message.
- `munge-gs --fwer` prints the FDR value in its header (cosmetic, not
  copied).

### 4.12 Fisher's exact test odds ratio

`scipy.stats.fisher_exact` returns the **sample** odds ratio $ad/bc$; R's
`fisher.test()` returns the conditional MLE. `scdrsTestOverlap()` takes the
p-value from `fisher.test()` and computes the sample OR and its Wald CI
itself, as Python does.

### 4.13 Memory and scale (not a parity issue, but a design constraint)

- Control score matrices are $n \times K$ doubles (100k cells × 1000
  controls ≈ 800 MB). Defaults follow Python (`returnCtrlNormScore = FALSE`
  in `scdrsScoreCells()`, `TRUE` in `scdrsComputeScore()`, which downstream needs).
- Implicit covariate correction keeps the assay sparse; dense mode needs a
  dense $n \times G$ copy, as in Python.
- `DelayedArray`/HDF5 assays: phase 2, via block processing.

---

## 5. Testing strategy (all functions)

- **Reference data.** Use the toy data bundled with scDRS
  (`toydata_mouse.h5ad`: 30 cells, 10 "causal" and 20 "non_causal";
  already log-normalized; plus `.cov` and `.gs` files and reference
  outputs). Export it to plain files (`.mtx`, TSV) and ship it in
  `inst/extdata/` (< 1 MB) for tests and the vignette.
- **Reference scripts** (`inst/scripts/python/`, run on the HPC):
  1. `export_toydata.py`: matrix, obs, var, cov, gs.
  2. `ref_preprocess.py`: `GENE_STATS`, `CELL_STATS`, `COV_BETA`,
     `COV_GENE_MEAN` for all 8 combinations of
     dense/sparse × cov/no cov × adj_prop/none; also `reg_out`,
     `_get_mean_var`, `_get_mean_var_implicit_cov_corr`, `category2dummy`
     on small inputs.
  3. `ref_score_cell.py`: control sets; disease and control raw scores for
     each `weight_opt`; all intermediate `_correct_background` stages
     (`save_intermediate`); final scores; `_get_p_from_empi_null`.
  4. `ref_downstream.py`: connectivities (`.mtx`); group / corr / gene
     outputs; `gearys_c`, `test_gearysc` (all three options);
     `_pearson_corr`.
  5. `ref_util.py`: `load_gs` with/without species mapping, `munge-gs`
     (p-value and z-score inputs, fdr / fwer / top-n), `test_overlap`,
     `meta_analysis`, `zsc2pval`/`pval2zsc`, normalization
     (`load_h5ad` with filtering + normalization on a raw-count matrix).
  6. One larger public dataset (on the HPC only, not shipped) for the
     statistical-equivalence tests in §4.1.
- **Tolerances.** Deterministic math: `1e-6` relative (double vs double) or
  `1e-5` (vs float32 output); LOESS `1e-6`; bins ≥ 99.5 % agreement; end to
  end with R's own random draws: statistical criteria (§4.1).
- **Test files** map one to one onto R source files, e.g.
  `test-scdrsPreprocess.R`, `test-scdrsScoreCells.R`, `test-scdrsTestGroups.R`.

---

## 6. Proposed implementation order

1. Reference scripts 1–2 + toy data export (you run them on the HPC).
2. `pp.py` port: `.regOut`, `.meanVar`, `.meanVarImplicit`,
   `.categoryToDummy`, `scdrsComputeGeneStats`, `scdrsPreprocess`.
3. Script 3 + `method.py` scoring: `scdrsSelectControlGenes`, `.rawScore`,
   `.overdispersionScore`, `.correctBackground`, `.empiricalPvalue`,
   `scdrsScoreCells`.
4. `util.py` I/O: `scdrsReadGeneSets`, `scdrsWriteGeneSets`, `scdrsGetHomologs`,
   `scdrsMungeGeneSets`, `scdrsNormalize`, `scdrsReadScores`,
   `scdrsWriteScores`, `scdrsComputeScore`.
5. Script 4 + downstream: `scdrsComputeGearysC`, `scdrsTestGearysC`, `.pearsonCor`,
   `scdrsTestGroups`, `scdrsTestCorrelation`, `scdrsRankGenes`,
   `scdrsPerformDownstream`, `scdrsBuildGraph`.
6. Small utilities: `scdrsZToP`, `scdrsPToZ`, `scdrsMetaAnalyze`, `scdrsTestOverlap`.
7. Vignette (toy data end to end), `R CMD check` + `BiocCheck` clean,
   statistical-equivalence run on a larger dataset.
8. Phase 2: plotting (`scdrsPlotGroupStats`, …), `DelayedArray` support, Rcpp
   where profiling shows it is needed, then Bioconductor submission.

---

## 7. Decisions (reviewed)

1. **Function names:** every exported function is `scdrs` + verb
   (+ object), camelCase: `scdrsPreprocess`, `scdrsScoreCells`,
   `scdrsComputeScore`, `scdrsTestGroups`, … (table in §8). Internal
   helpers keep a leading dot (`.rawScore`, …) and are not exported.
2. **Where results live:** per-gene results in `rowData(sce)` and per-cell
   results in `colData(sce)`, all with the `scdrs_` prefix; the covariate
   matrix in `reducedDim(sce, "scdrs_cov")`; only flags and parameters in
   `metadata(sce)$scdrs` (§2.1 table).
3. **RNG:** no `seed` argument. Users call `set.seed()` themselves; the
   examples, vignette and tests do.
4. **Exact-RNG helper:** skipped, on the condition that R's own random
   control sets reproduce Python's final results (acceptance criteria
   below). If they don't, this decision is revisited.

   **Acceptance criteria** (R with its own random control sets vs. Python,
   same data and settings, `nCtrl = 1000`):
   - Pearson correlation of `norm_score` ≥ 0.99;
   - Pearson correlation of `-log10(pval)` ≥ 0.99;
   - cells called at FDR < 0.1 agree for ≥ 99 % of cells;
   - group-level `assoc_mcz`, `hetero_mcz` and `corr_mcz` correlate
     ≥ 0.99 across groups/variables.

   Checked on the toy data and on one larger dataset run on the HPC.
5. **Scope of phase 1:** everything in §8 marked phase 1. Plotting and
   `DelayedArray`/HDF5 support are phase 2.

---

## 8. Exported function names (decided: `scdrs` + verb)

| Python | R | Phase |
|---|---|---|
| `pp.preprocess` | `scdrsPreprocess()` | 1 |
| `pp.compute_stats` | `scdrsComputeGeneStats()` | 1 |
| `method.score_cell` | `scdrsScoreCells()` | 1 |
| `method._select_ctrl_geneset` | `scdrsSelectControlGenes()` | 1 |
| CLI `compute-score` | `scdrsComputeScore()` | 1 |
| `method.downstream_group_analysis` | `scdrsTestGroups()` | 1 |
| `method.downstream_corr_analysis` | `scdrsTestCorrelation()` | 1 |
| `method.downstream_gene_analysis` | `scdrsRankGenes()` | 1 |
| CLI `perform-downstream` | `scdrsPerformDownstream()` | 1 |
| `method.test_gearysc` | `scdrsTestGearysC()` | 1 |
| `method.gearys_c` | `scdrsComputeGearysC()` | 1 |
| `sc.pp.pca` + `sc.pp.neighbors` | `scdrsBuildGraph()` | 1 |
| `util.load_h5ad` (filter + normalize) | `scdrsNormalize()` | 1 |
| `util.load_gs` | `scdrsReadGeneSets()` | 1 |
| `util.save_gs` | `scdrsWriteGeneSets()` | 1 |
| CLI `munge-gs` | `scdrsMungeGeneSets()` | 1 |
| `util.load_homolog_mapping` | `scdrsGetHomologs()` | 1 |
| `util.load_scdrs_score` | `scdrsReadScores()` | 1 |
| (CLI file output) | `scdrsWriteScores()` | 1 |
| `util.zsc2pval` / `util.pval2zsc` | `scdrsZToP()` / `scdrsPToZ()` | 1 |
| `util.meta_analysis` | `scdrsMetaAnalyze()` | 1 |
| `util.test_overlap` | `scdrsTestOverlap()` | 1 |
| `method.score_cell_vision` | `scdrsScoreVision()` | optional |
| `util.plot_group_stats` | `scdrsPlotGroupStats()` | 2 |
| `util.plot_qq` | `scdrsPlotQQ()` | 2 |
| `util.plot_score_umap` | `scdrsPlotUmap()` | 2 |
