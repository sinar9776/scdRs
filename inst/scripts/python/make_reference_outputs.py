#!/usr/bin/env python
"""
Make Python scDRS reference outputs for the R port (scdRs).

Runs scDRS (v1.0.4, commit 644a5c7) on its own bundled toy data and writes
plain CSV files that R can read without Python.

Usage:
    reference/venv/bin/python reference/make_reference_outputs.py \
        [--out reference/outputs] [--n-ctrl 1000] [--seed 0]

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
    p.add_argument("--n-ctrl", type=int, default=1000)
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


def control_sets(adata, gene_list, gene_weight, n_ctrl, seed):
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
        df_gene, genes, weights, "mean_var", n_ctrl, 200, seed
    )
    return genes, weights, dic_list, dic_weight


def export_scores(adata, tag, trait, gene_list, gene_weight, n_ctrl, seed, out):
    df = scdrs.score_cell(
        adata,
        gene_list,
        gene_weight=gene_weight,
        ctrl_match_key="mean_var",
        n_ctrl=n_ctrl,
        weight_opt="vs",
        return_ctrl_raw_score=True,
        return_ctrl_norm_score=True,
        random_seed=seed,
    )
    df.index.name = "cell"
    prefix = f"{tag}.{trait}"
    main_cols = ["raw_score", "norm_score", "mc_pval", "pval",
                 "nlog10_pval", "zscore"]
    raw_cols = [f"ctrl_raw_score_{i}" for i in range(n_ctrl)]
    norm_cols = [f"ctrl_norm_score_{i}" for i in range(n_ctrl)]
    write_csv(df[main_cols], out, f"{prefix}.score.csv")
    write_csv(df[raw_cols], out, f"{prefix}.ctrl_raw_score.csv")
    write_csv(df[norm_cols], out, f"{prefix}.ctrl_norm_score.csv")

    # Control gene sets, long format: ctrl_id (0-based), gene, weight
    genes, weights, dic_list, dic_weight = control_sets(
        adata, gene_list, gene_weight, n_ctrl, seed
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
            adata, dic_list[i], dic_weight[i], "vs"
        )
        if not np.allclose(v, df[f"ctrl_raw_score_{i}"].values, rtol=1e-5):
            sys.exit(f"control set {i} does not reproduce score_cell output")
    return df


def export_downstream(adata, df_full, trait, out):
    """Downstream analyses, the same steps as `scdrs perform-downstream`."""
    n_pcs = min(20, min(adata.shape) - 1)
    sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)

    g = adata.obsp["connectivities"].tocoo()
    conn = pd.DataFrame({
        "cell_i": adata.obs_names[g.row],
        "cell_j": adata.obs_names[g.col],
        "weight": g.data,
    })
    write_csv(conn, out, "downstream.connectivities.csv", index=False)

    group = scdrs.method.downstream_group_analysis(
        adata=adata, df_full_score=df_full, group_cols=["cell_type"]
    )["cell_type"]
    write_csv(group, out, f"downstream.{trait}.group.cell_type.csv")

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
    for tag, use_cov in SCENARIOS.items():
        adata = adata0.copy()
        scdrs.preprocess(
            adata, cov=df_cov if use_cov else None,
            n_mean_bin=20, n_var_bin=20, copy=False,
        )
        export_preprocess(adata, tag, args.out)
        for trait, (genes, weights) in dict_gs.items():
            full_scores[(tag, trait)] = export_scores(
                adata, tag, trait, genes, weights,
                args.n_ctrl, args.seed, args.out,
            )

    # Preprocessing in every mode, and the internal helpers
    export_pp_scenarios(adata0, df_cov, args.out)
    export_pp_nonconstant(adata0, df_cov, args.out)
    export_units(args.out)

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
