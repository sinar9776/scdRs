#!/usr/bin/env python
"""
Make Python scDRS reference outputs for the R port (RscDRS).

Runs scDRS (v1.0.4, commit 644a5c7) on its own bundled toy data and writes
plain CSV files that R can read without Python.

Usage:
    reference/venv/bin/python reference/make_reference_outputs.py \
        [--out reference/outputs] [--n-ctrl 100] [--seed 0]

Orientation: every expression-like matrix is written genes x cells (rows =
genes), i.e. the SingleCellExperiment orientation, not AnnData's
cells x genes. Per-cell tables have one row per cell.

Note: `mean_var` bin labels look like numbers ("3.12"); read them in R as
character, e.g. read.csv(..., colClasses = c(mean_var = "character")).
"""

import argparse
import os
import shutil
import sys
from importlib.metadata import version

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from anndata import read_h5ad

import scdrs

SCENARIOS = {
    # name: use covariates? (toy data is sparse, so cov -> implicit
    # covariate-correction mode, the same as the CLI test)
    "cov": True,
    "nocov": False,
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    here = os.path.dirname(os.path.abspath(__file__))
    p.add_argument("--out", default=os.path.join(here, "outputs"))
    p.add_argument("--n-ctrl", type=int, default=100)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def write_csv(df, out, name, index=True):
    path = os.path.join(out, name)
    df.to_csv(path, index=index)
    return path


def export_toydata(adata, df_cov, data_dir, out):
    """Toy data in R-readable form."""
    X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
    expr = pd.DataFrame(X.T, index=adata.var_names, columns=adata.obs_names)
    expr.index.name = "gene"
    write_csv(expr, out, "toydata_logexpr.csv")
    obs = adata.obs.copy()
    obs.index.name = "cell"
    write_csv(obs, out, "toydata_cell_metadata.csv")
    var = adata.var.copy()
    var.index.name = "gene"
    write_csv(var, out, "toydata_gene_metadata.csv")
    cov = df_cov.copy()
    cov.index.name = "cell"
    write_csv(cov, out, "toydata_covariates.csv")
    # Original .gs files, unchanged (tab-separated TRAIT / GENESET)
    for f in ["toydata_mouse.gs", "toydata_human.gs"]:
        shutil.copy(os.path.join(data_dir, f), os.path.join(out, f))


def export_gene_sets(dict_gs, out):
    """Gene sets as loaded by scdrs.util.load_gs (after species mapping)."""
    rows = []
    for trait, (genes, weights) in dict_gs.items():
        rows += [(trait, g, w) for g, w in zip(genes, weights)]
    df = pd.DataFrame(rows, columns=["trait", "gene", "weight"])
    write_csv(df, out, "genesets_loaded.csv", index=False)


def export_preprocess(adata, tag, out):
    param = adata.uns["SCDRS_PARAM"]
    gs = param["GENE_STATS"].copy()
    gs.index.name = "gene"
    write_csv(gs, out, f"{tag}.gene_stats.csv")
    cs = param["CELL_STATS"].copy()
    cs.index.name = "cell"
    write_csv(cs, out, f"{tag}.cell_stats.csv")
    if "COV_BETA" in param:  # only stored in sparse + cov mode
        beta = param["COV_BETA"].copy()
        beta.index.name = "gene"
        write_csv(beta, out, f"{tag}.cov_beta.csv")
        mu = param["COV_GENE_MEAN"].rename("cov_gene_mean").to_frame()
        mu.index.name = "gene"
        write_csv(mu, out, f"{tag}.cov_gene_mean.csv")
        cmat = param["COV_MAT"].copy()
        cmat.index.name = "cell"
        write_csv(cmat, out, f"{tag}.cov_mat.csv")
    flags = pd.DataFrame(
        {k: [param[k]] for k in ["FLAG_SPARSE", "FLAG_COV", "FLAG_ADJ_PROP"]}
    )
    write_csv(flags, out, f"{tag}.flags.csv", index=False)


def control_sets(adata, gene_list, gene_weight, n_ctrl, seed,
                 ctrl_match_key="mean_var"):
    """Re-create the control gene sets exactly as score_cell() does.

    score_cell() does not return its control sets, so this repeats its
    gene-table preparation and calls the same subroutine with the same seed.
    """
    df_gene = adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[adata.var_names].copy()
    df_gene["gene"] = df_gene.index
    df_gene.drop_duplicates(subset="gene", inplace=True)
    dic_w = dict(zip(gene_list, gene_weight))
    genes = sorted(set(gene_list) & set(df_gene["gene"]))
    weights = [dic_w[g] for g in genes]
    dic_list, dic_weight = scdrs.method._select_ctrl_geneset(
        df_gene, genes, weights, ctrl_match_key, n_ctrl, 200, seed
    )
    return genes, weights, dic_list, dic_weight


MAIN_COLS = ["raw_score", "norm_score", "mc_pval", "pval", "nlog10_pval",
             "zscore"]


def export_scores(adata, tag, trait, gene_list, gene_weight, n_ctrl, seed,
                  out, weight_opt="vs", ctrl_match_key="mean_var",
                  norm=True):
    """score_cell() output and the control gene sets it used.

    Files <tag>.<trait>.{score,ctrl_raw_score,ctrl_norm_score,
    ctrl_genesets}.csv (ctrl_norm_score only with norm=True).
    """
    df = scdrs.score_cell(
        adata,
        gene_list,
        gene_weight=gene_weight,
        ctrl_match_key=ctrl_match_key,
        n_ctrl=n_ctrl,
        weight_opt=weight_opt,
        return_ctrl_raw_score=True,
        return_ctrl_norm_score=True,
        random_seed=seed,
    )
    df.index.name = "cell"
    prefix = f"{tag}.{trait}"
    raw_cols = [f"ctrl_raw_score_{i}" for i in range(n_ctrl)]
    norm_cols = [f"ctrl_norm_score_{i}" for i in range(n_ctrl)]
    write_csv(df[MAIN_COLS], out, f"{prefix}.score.csv")
    write_csv(df[raw_cols], out, f"{prefix}.ctrl_raw_score.csv")
    if norm:
        write_csv(df[norm_cols], out, f"{prefix}.ctrl_norm_score.csv")

    # Control gene sets, long format: ctrl_id (0-based), gene, weight
    genes, weights, dic_list, dic_weight = control_sets(
        adata, gene_list, gene_weight, n_ctrl, seed, ctrl_match_key
    )
    rows = []
    for i in range(n_ctrl):
        rows += [(i, g, w) for g, w in zip(dic_list[i], dic_weight[i])]
    ctrl = pd.DataFrame(rows, columns=["ctrl_id", "gene", "weight"])
    write_csv(ctrl, out, f"{prefix}.ctrl_genesets.csv", index=False)

    # Check: the exported control sets reproduce score_cell's control raw
    # scores, so they are exactly the sets score_cell used.
    for i in [0, n_ctrl // 2, n_ctrl - 1]:
        v, _ = scdrs.method._compute_raw_score(
            adata, dic_list[i], dic_weight[i], weight_opt
        )
        if not np.allclose(v, df[f"ctrl_raw_score_{i}"].values, rtol=1e-5):
            sys.exit(f"control set {i} does not reproduce score_cell output")
    return df


def export_downstream(adata, df_full, trait, out):
    """Downstream analyses, the same steps as `scdrs perform-downstream`."""
    n_pcs = min(20, min(adata.shape) - 1)
    sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)

    pca = pd.DataFrame(adata.obsm["X_pca"], index=adata.obs_names,
                       columns=[f"PC{i + 1}" for i in range(n_pcs)])
    pca.index.name = "cell"
    write_csv(pca, out, "downstream.pca.csv")

    g = adata.obsp["connectivities"].tocoo()
    conn = pd.DataFrame({
        "cell_i": adata.obs_names[g.row],
        "cell_j": adata.obs_names[g.col],
        "weight": g.data,
    })
    write_csv(conn, out, "downstream.connectivities.csv", index=False)

    # The full score file exactly as `scdrs compute-score` writes it
    df_full.to_csv(os.path.join(out, f"downstream.{trait}.full_score.gz"),
                   sep="\t", index=True, compression="gzip")

    # A second annotation with 3 groups (cells taken in turn), to test more
    # than two groups; R adds the same column.
    adata.obs["test_group"] = [f"g{i % 3}" for i in range(adata.n_obs)]
    for col in ["cell_type", "test_group"]:
        group = scdrs.method.downstream_group_analysis(
            adata=adata, df_full_score=df_full, group_cols=[col]
        )[col]
        write_csv(group, out, f"downstream.{trait}.group.{col}.csv")

    # test_gearysc, deterministic options, on all cells
    cells = list(adata.obs_names)
    for opt in ["control_distribution_match", "control"]:
        rls = scdrs.method.test_gearysc(
            adata, df_full.loc[cells], groupby="cell_type", opt=opt
        )
        rls.index.name = "group"
        write_csv(rls, out, f"downstream.gearysc.{opt}.csv")

    # gearys_c on the whole graph: disease score and the first 5 controls
    vals = df_full.loc[cells, ["norm_score"] + [
        f"ctrl_norm_score_{i}" for i in range(5)]]
    c = scdrs.method.gearys_c(adata, vals.values)
    write_csv(pd.DataFrame({"column": vals.columns, "C": c}), out,
              "downstream.gearys_c.csv", index=False)

    corr = scdrs.method.downstream_corr_analysis(
        adata=adata, df_full_score=df_full,
        var_cols=["causal_variable", "non_causal_variable", "covariate"],
    )
    corr.index.name = "variable"
    write_csv(corr, out, f"downstream.{trait}.corr.csv")

    gene = scdrs.method.downstream_gene_analysis(
        adata=adata, df_full_score=df_full
    )
    gene.index.name = "gene"
    write_csv(gene, out, f"downstream.{trait}.gene.csv")


def export_pp_scenarios(adata0, df_cov, out):
    """preprocess() in all 8 modes: sparse/dense x cov/nocov x adj_prop.

    Files: pp.<x>.<cov>.<adj>.{gene_stats,cell_stats,cov_beta,
    cov_gene_mean,cov_mat,flags}.csv. Dense + cov also writes the
    covariate-corrected matrix (genes x cells) as .corrected.csv.
    """
    for x_mode in ["sparse", "dense"]:
        for cov_mode in ["cov", "nocov"]:
            for adj in ["none", "celltype"]:
                adata = adata0.copy()
                if x_mode == "dense":
                    adata.X = adata.X.toarray()
                scdrs.preprocess(
                    adata,
                    cov=df_cov if cov_mode == "cov" else None,
                    adj_prop="cell_type" if adj == "celltype" else None,
                    n_mean_bin=20, n_var_bin=20, copy=False,
                )
                tag = f"pp.{x_mode}.{cov_mode}.{adj}"
                export_preprocess(adata, tag, out)
                if x_mode == "dense" and cov_mode == "cov":
                    corr = pd.DataFrame(
                        np.asarray(adata.X).T,
                        index=adata.var_names, columns=adata.obs_names,
                    )
                    corr.index.name = "gene"
                    write_csv(corr, out, f"{tag}.corrected.csv")


def export_pp_nonconstant(adata0, df_cov, out):
    """preprocess() with covariates on the toy data without the genes that
    are constant after covariate correction.

    In the toy data 10 genes are an exact linear function of the covariate.
    After correction their variance is floating-point noise (~1e-28), which
    differs between implementations but still enters the LOESS fit. Without
    them, the results are free of that noise.
    Files: pp.<x>.cov.<adj>.nonconst.*.csv and pp.constant_genes.csv.
    """
    adata = adata0.copy()
    scdrs.preprocess(adata, cov=df_cov, n_mean_bin=20, n_var_bin=20)
    gs = adata.uns["SCDRS_PARAM"]["GENE_STATS"]
    const = gs.index[gs["ct_var"].values < 1e-12 * gs["ct_mean"].values ** 2]
    write_csv(pd.DataFrame({"gene": const}), out, "pp.constant_genes.csv",
              index=False)
    keep = [g for g in adata0.var_names if g not in set(const)]
    for x_mode in ["sparse", "dense"]:
        for adj in ["none", "celltype"]:
            adata = adata0[:, keep].copy()
            if x_mode == "dense":
                adata.X = adata.X.toarray()
            scdrs.preprocess(
                adata, cov=df_cov,
                adj_prop="cell_type" if adj == "celltype" else None,
                n_mean_bin=20, n_var_bin=20, copy=False,
            )
            export_preprocess(adata, f"pp.{x_mode}.cov.{adj}.nonconst", out)


def export_units(out):
    """Small inputs/outputs for the internal pp.py helpers."""
    rng = np.random.RandomState(1)

    # reg_out
    Y = rng.normal(size=(20, 3))
    X = np.column_stack([rng.normal(size=20), np.ones(20)])
    R = scdrs.pp.reg_out(Y, X)
    write_csv(pd.DataFrame(Y, columns=["y1", "y2", "y3"]), out,
              "units.reg_out.Y.csv", index=False)
    write_csv(pd.DataFrame(X, columns=["x1", "const"]), out,
              "units.reg_out.X.csv", index=False)
    write_csv(pd.DataFrame(R, columns=["y1", "y2", "y3"]), out,
              "units.reg_out.resid.csv", index=False)

    # category2dummy: strings (mixed case, NA), numbers with NA, bool
    df = pd.DataFrame({
        "num": [1.0, 2.0, np.nan, 4.0, 5.0, 6.0],
        "grp": ["b", "a", "C", None, "a", "b"],
        "flag": [True, False, True, True, False, False],
    }, index=[f"c{i}" for i in range(6)])
    df.index.name = "cell"
    write_csv(df, out, "units.category2dummy.input.csv")
    dummy = scdrs.pp.category2dummy(df).astype(float)
    dummy = dummy[["num", "flag"] + sorted(
        c for c in dummy.columns if c.startswith("grp_"))]
    write_csv(dummy, out, "units.category2dummy.output.csv")

    # _get_mean_var: dense / sparse, weighted / unweighted, both axes.
    # Matrix is cells x genes (Python orientation).
    M = rng.poisson(1.0, size=(8, 5)).astype(float)
    w = rng.uniform(0.5, 2.0, size=8)
    write_csv(pd.DataFrame(M), out, "units.mean_var.X.csv", index=False)
    write_csv(pd.DataFrame({"w": w}), out, "units.mean_var.w.csv",
              index=False)
    rows = []
    for kind, mat in [("dense", M), ("sparse", sp.csr_matrix(M))]:
        for axis in [0, 1]:
            for wt in ["none", "w"]:
                if wt == "w" and axis == 1:
                    continue
                m, v = scdrs.pp._get_mean_var(
                    mat, axis=axis, weights=None if wt == "none" else w)
                for i, (mi, vi) in enumerate(zip(m, v)):
                    rows.append((kind, axis, wt, i, mi, vi))
    write_csv(pd.DataFrame(rows, columns=[
        "kind", "axis", "weights", "index", "mean", "var"]), out,
        "units.mean_var.output.csv", index=False)

    # pd.qcut(labels=False, duplicates="drop"), incl. ties and values
    # lying exactly on a cut point
    x = np.concatenate([np.repeat(0.0, 7), np.arange(1, 25) / 4.0,
                        rng.normal(size=20)])
    q = pd.qcut(x, 7, labels=False, duplicates="drop")
    write_csv(pd.DataFrame({"x": x, "bin": q}), out, "units.qcut.csv",
              index=False)


def export_score_units(out):
    """Small inputs/outputs for the internal scoring and downstream helpers:
    _correct_background (with exact zeros), _get_p_from_empi_null (with
    ties), _pearson_corr (dense and sparse)."""
    rng = np.random.RandomState(2)
    n, k = 12, 6
    raw = rng.normal(size=n)
    ctrl = rng.normal(size=(n, k))
    raw[[1, 5]] = 0.0
    ctrl[[0, 3, 7], [2, 4, 0]] = 0.0
    ratio = rng.uniform(0.5, 2.0, size=k)
    norm, ctrl_norm = scdrs.method._correct_background(raw, ctrl.copy(), ratio)
    cols = [f"c{i}" for i in range(k)]
    write_csv(pd.DataFrame({"raw": raw, "norm": norm}), out,
              "units.correct_background.score.csv", index=False)
    write_csv(pd.DataFrame(ctrl, columns=cols), out,
              "units.correct_background.ctrl_raw.csv", index=False)
    write_csv(pd.DataFrame(ctrl_norm, columns=cols), out,
              "units.correct_background.ctrl_norm.csv", index=False)
    write_csv(pd.DataFrame({"ratio": ratio}), out,
              "units.correct_background.ratio.csv", index=False)

    t = np.array([-1.0, 0.0, 0.5, 0.5, 2.0, 10.0])
    t_null = np.array([0.5, 0.5, -3.0, 1.0, 0.0, 2.0, 2.0, 0.25])
    p = scdrs.method._get_p_from_empi_null(t, t_null)
    write_csv(pd.DataFrame({"t": t, "p": p}), out, "units.empi_null.t.csv",
              index=False)
    write_csv(pd.DataFrame({"t_null": t_null}), out,
              "units.empi_null.null.csv", index=False)

    # _pearson_corr: X has a constant column (sd clipped at 1e-8)
    X = rng.poisson(0.7, size=(15, 4)).astype(float)
    X[:, 3] = 0.0
    Y = rng.normal(size=(15, 2))
    write_csv(pd.DataFrame(X), out, "units.pearson.X.csv", index=False)
    write_csv(pd.DataFrame(Y), out, "units.pearson.Y.csv", index=False)
    rows = []
    for kind, mx in [("dense", X), ("sparse", sp.csr_matrix(X))]:
        r = scdrs.method._pearson_corr(mx, Y)
        for i in range(r.shape[0]):
            for j in range(r.shape[1]):
                rows.append((kind, i, j, float(r[i, j])))
    write_csv(pd.DataFrame(rows, columns=["kind", "i", "j", "corr"]), out,
              "units.pearson.output.csv", index=False)


def load_cli():
    """The scdrs command-line script (bin/scdrs) as a module."""
    import importlib.machinery
    import importlib.util
    import shutil as sh
    path = os.path.join(os.path.dirname(sys.executable), "scdrs")
    if not os.path.exists(path):
        path = sh.which("scdrs")
    loader = importlib.machinery.SourceFileLoader("scdrs_cli", path)
    spec = importlib.util.spec_from_loader("scdrs_cli", loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


def write_gs_table(dict_gs, out, name):
    rows = []
    for trait, (genes, weights) in dict_gs.items():
        rows += [(trait, i, g, w) for i, (g, w) in
                 enumerate(zip(genes, weights))]
    write_csv(pd.DataFrame(rows, columns=["trait", "pos", "gene", "weight"]),
              out, name, index=False)


def export_util(adata0, data_dir, out):
    """util.py functions and `scdrs munge-gs`."""
    mouse_gs = os.path.join(data_dir, "toydata_mouse.gs")
    human_gs = os.path.join(data_dir, "toydata_human.gs")

    # load_gs: as is, mouse -> human, human -> mouse, with/without intersect
    write_gs_table(scdrs.util.load_gs(mouse_gs), out, "util.load_gs.mouse.csv")
    write_gs_table(scdrs.util.load_gs(mouse_gs, "mouse", "human"), out,
                   "util.load_gs.mouse2human.csv")
    write_gs_table(scdrs.util.load_gs(human_gs, "human", "mouse"), out,
                   "util.load_gs.human2mouse.csv")
    write_gs_table(scdrs.util.load_gs(human_gs, "human", "mouse",
                                      to_intersect=adata0.var_names[:1500]),
                   out, "util.load_gs.human2mouse.intersect.csv")

    # Weighted .gs with a duplicated gene (last weight wins, first position
    # kept) and a second, unweighted trait
    gs_file = os.path.join(out, "util.weighted.gs")
    with open(gs_file, "w") as f:
        f.write("TRAIT\tGENESET\n")
        f.write("T1\tGene_b:1.5,Gene_a:2,Gene_b:0.25,Gene_c:-1e-3\n")
        f.write("T2\tGene_c,Gene_a\n")
    dict_w = scdrs.util.load_gs(gs_file)
    write_gs_table(dict_w, out, "util.load_gs.weighted.csv")
    scdrs.util.save_gs(os.path.join(out, "util.save_gs.gs"), dict_w)

    # munge-gs on simulated p-values (genes x 2 traits, with NAs and ties)
    rng = np.random.RandomState(3)
    genes = [f"G{i:03d}" for i in range(300)]
    pv = pd.DataFrame({
        "TRAIT_B": rng.uniform(size=300) ** 3,
        "TRAIT_A": rng.uniform(size=300) ** 6,
    }, index=pd.Index(genes, name="GENE"))
    pv.iloc[[3, 50, 77], 0] = np.nan
    pv.iloc[:5, 1] = 1e-120  # clipped at 1e-100 for the weights
    pv_file = os.path.join(out, "util.munge.pval.tsv")
    pv.to_csv(pv_file, sep="\t")
    zs = pd.DataFrame(scdrs.util.pval2zsc(pv.values), index=pv.index,
                      columns=pv.columns)
    zs_file = os.path.join(out, "util.munge.zscore.tsv")
    zs.to_csv(zs_file, sep="\t")
    cli = load_cli()
    runs = {
        "fdr": dict(pval_file=pv_file, fdr=0.05, n_min=10, n_max=100),
        "fwer": dict(pval_file=pv_file, fwer=0.05, n_min=5, n_max=100),
        "top": dict(pval_file=pv_file, n_min=10, n_max=40),
        "zscore_uniform": dict(zscore_file=zs_file, weight="uniform",
                               fdr=0.1, n_min=10, n_max=1000),
    }
    for name, kw in runs.items():
        cli.munge_gs(out_file=os.path.join(out, f"util.munge.{name}.gs"),
                     **kw)

    # test_overlap
    bg = [f"g{i}" for i in range(100)]
    cases = {
        "basic": (bg[:30], bg[20:45]),
        "no_overlap": (bg[:10], bg[50:60]),
        "subset": (bg[:10], bg[:40]),
    }
    rows = []
    for name, (a, b) in cases.items():
        pval, odds, ub, lb = scdrs.util.test_overlap(a, b, bg)
        rows.append((name, pval, odds, ub, lb))
    write_csv(pd.DataFrame(rows, columns=["case", "pval", "odds_ratio",
                                          "or_ub", "or_lb"]),
              out, "util.test_overlap.csv", index=False)

    # meta_analysis
    eff = np.array([0.2, 0.5, -0.1, 0.35, 0.8])
    se = np.array([0.1, 0.2, 0.15, 0.05, 0.3])
    rows = []
    for method in ["random", "fixed"]:
        s, s_se = scdrs.util.meta_analysis(eff, se, method=method)
        rows.append((method, "none", s, s_se))
    s, s_se = scdrs.util.meta_analysis(eff, se, method="random",
                                       weights=np.array([1, 2, 1, 3, 1.0]))
    rows.append(("random", "w", s, s_se))
    write_csv(pd.DataFrame(rows, columns=["method", "weights", "summary",
                                          "se"]),
              out, "util.meta_analysis.csv", index=False)
    write_csv(pd.DataFrame({"effect": eff, "se": se}), out,
              "util.meta_analysis.input.csv", index=False)

    # zsc2pval / pval2zsc, including the far tail
    z = np.array([-5.0, -1.0, 0.0, 1.0, 2.5, 8.0, 30.0])
    p = scdrs.util.zsc2pval(z)
    write_csv(pd.DataFrame({"z": z, "p": p, "z_back":
                            scdrs.util.pval2zsc(p)}), out,
              "util.zsc_pval.csv", index=False)

    # load_h5ad: filtering + normalize_per_cell + log1p on raw counts
    # (cells x genes Poisson counts; a few cells / genes fall below the
    # thresholds; one cell has no counts at all)
    counts = rng.poisson(rng.uniform(0.05, 3.0, size=(1, 60)),
                         size=(40, 60)).astype(np.float32)
    counts[5, :] = 0
    counts[7, 30:] = 0
    counts[:, 11] = 0
    counts[:35, 12] = 0
    cells = [f"cell{i}" for i in range(40)]
    gnames = [f"gene{j}" for j in range(60)]
    raw = sc.AnnData(sp.csr_matrix(counts),
                     obs=pd.DataFrame(index=cells),
                     var=pd.DataFrame(index=gnames))
    df = pd.DataFrame(counts.T, index=gnames, columns=cells)
    df.index.name = "gene"
    write_csv(df, out, "util.normalize.counts.csv")
    h5 = os.path.join(out, "util.normalize.h5ad")
    raw.write_h5ad(h5)
    for name, flt in [("filter", True), ("nofilter", False)]:
        ad = scdrs.util.load_h5ad(h5, flag_filter_data=flt,
                                  flag_raw_count=True, min_genes=20,
                                  min_cells=6)
        X = ad.X.toarray() if sp.issparse(ad.X) else ad.X
        res = pd.DataFrame(X.T, index=ad.var_names, columns=ad.obs_names)
        res.index.name = "gene"
        write_csv(res, out, f"util.normalize.{name}.csv")
    os.remove(h5)


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    data_dir = os.path.join(scdrs.__path__[0], "data")

    adata0 = read_h5ad(os.path.join(data_dir, "toydata_mouse.h5ad"))
    df_cov = pd.read_csv(
        os.path.join(data_dir, "toydata_mouse.cov"), sep="\t", index_col=0
    )
    df_cov.index = [str(x) for x in df_cov.index]

    export_toydata(adata0, df_cov, data_dir, args.out)

    # Gene sets: mouse set as is; human set mapped to mouse (as the CLI does)
    dict_gs = {}
    dict_gs.update(scdrs.util.load_gs(
        os.path.join(data_dir, "toydata_mouse.gs"),
        src_species="mouse", dst_species="mouse",
        to_intersect=adata0.var_names,
    ))
    dict_gs.update(scdrs.util.load_gs(
        os.path.join(data_dir, "toydata_human.gs"),
        src_species="human", dst_species="mouse",
        to_intersect=adata0.var_names,
    ))
    export_gene_sets(dict_gs, args.out)

    full_scores = {}
    mouse = "toydata_gs_mouse"
    for tag, use_cov in SCENARIOS.items():
        adata = adata0.copy()
        scdrs.preprocess(
            adata, cov=df_cov if use_cov else None,
            n_mean_bin=20, n_var_bin=20, copy=False,
        )
        export_preprocess(adata, tag, args.out)
        for trait, (genes, weights) in dict_gs.items():
            # Full output (control scores, control sets) with n_ctrl
            # controls, for the exact tests with injected control sets
            full_scores[(tag, trait)] = export_scores(
                adata, tag, trait, genes, weights,
                args.n_ctrl, args.seed, args.out,
            )
            # Main columns only, with 1000 controls, for the statistical
            # comparison with R's own random control sets
            df = scdrs.score_cell(
                adata, genes, gene_weight=weights, ctrl_match_key="mean_var",
                n_ctrl=1000, weight_opt="vs", random_seed=args.seed,
            )
            df.index.name = "cell"
            write_csv(df[MAIN_COLS], args.out,
                      f"{tag}.{trait}.score1000.csv")
        # Other weighting options and ctrl_match_key = "mean"
        genes, weights = dict_gs[mouse]
        for opt in ["uniform", "inv_std", "od"]:
            export_scores(adata, f"{tag}.{opt}", mouse, genes, weights,
                          args.n_ctrl, args.seed, args.out,
                          weight_opt=opt, norm=False)
        export_scores(adata, f"{tag}.matchmean", mouse, genes, weights,
                      args.n_ctrl, args.seed, args.out,
                      ctrl_match_key="mean", norm=False)
        # Weighted gene set: weights 1, 2, ..., 5 in turn, in the order of
        # the .gs file
        w5 = [1.0 + (i % 5) for i in range(len(genes))]
        export_scores(adata, f"{tag}.weighted", mouse, genes, w5,
                      args.n_ctrl, args.seed, args.out, norm=False)

    # Dense data with covariates (corrected data stored in X)
    adata = adata0.copy()
    adata.X = adata.X.toarray()
    scdrs.preprocess(adata, cov=df_cov, n_mean_bin=20, n_var_bin=20)
    genes, weights = dict_gs[mouse]
    export_scores(adata, "densecov", mouse, genes, weights, args.n_ctrl,
                  args.seed, args.out, norm=False)

    # Preprocessing in every mode, and the internal helpers
    export_pp_scenarios(adata0, df_cov, args.out)
    export_pp_nonconstant(adata0, df_cov, args.out)
    export_units(args.out)
    export_score_units(args.out)
    export_util(adata0, data_dir, args.out)

    # Downstream on a fresh copy of the data (the CLI reloads the .h5ad and
    # does not run preprocess), using the covariate-corrected mouse scores.
    trait = "toydata_gs_mouse"
    export_downstream(adata0.copy(), full_scores[("cov", trait)], trait,
                      args.out)

    pkgs = ["scdrs", "scanpy", "anndata", "numpy", "scipy", "pandas",
            "scikit-misc", "statsmodels", "scikit-learn", "umap-learn",
            "pynndescent"]
    info = pd.DataFrame({
        "item": ["python", "seed", "n_ctrl", "scdrs_commit"] + pkgs,
        "value": [sys.version.split()[0], args.seed, args.n_ctrl,
                  "644a5c7"] + [version(p) for p in pkgs],
    })
    write_csv(info, args.out, "session_info.csv", index=False)
    print(f"Done: {len(os.listdir(args.out))} files in {args.out}")


if __name__ == "__main__":
    main()
