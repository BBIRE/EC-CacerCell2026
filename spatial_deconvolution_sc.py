#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import re
import sys
import json
import argparse
import warnings
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
import h5py

from scipy.optimize import nnls
from joblib import Parallel, delayed

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings(
    "ignore",
    message=".*DataFrame is highly fragmented.*",
    category=pd.errors.PerformanceWarning,
)

def normalize_key(s: str) -> str:
    return "".join(str(s).strip().lower().split())

def resolve_celltype(requested: str, available: list[str]) -> str:
    if requested is None:
        return None
    if requested in available:
        return requested
    req = normalize_key(requested)
    avail_norm = {normalize_key(a): a for a in available}
    if req in avail_norm:
        return avail_norm[req]
    for a in available:
        an = normalize_key(a)
        if req in an or an in req:
            return a
    syn = {
        "macrophage": ["em1", "em2", "monocyte", "myeloid", "immune_myeloid", "immune myeloid"],
        "nk": ["unk", "unk1", "unk2", "unk3", "ilc", "ilc3"],
        "tcell": ["t_cell", "t cell", "t_reg", "treg", "cd4", "cd8"],
        "bcell": ["b_cell", "plasma_b_cell", "plasma cell"],
        "dc": ["cdc", "cdc1", "cdc2", "pdc", "dendritic"],
    }
    if req in syn:
        for token in syn[req]:
            tok = normalize_key(token)
            for a in available:
                if tok in normalize_key(a):
                    return a
    raise ValueError(
        f"Celltype '{requested}' not found. Available examples: {available[:30]} (total={len(available)})."
    )

def safe_mkdir(p: str | Path):
    Path(p).mkdir(parents=True, exist_ok=True)

def sanitize_name(s: str) -> str:
    s = str(s).strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_.-]+", "", s)
    return s or "sample"

def assert_h5ad_valid(path: str, label: str):
    if not os.path.exists(path):
        raise FileNotFoundError(f"[{label}] not found: {path}")
    if os.path.isdir(path):
        raise IsADirectoryError(f"[{label}] is a directory, expected .h5ad file: {path}")
    size = os.path.getsize(path)
    if size < 1024:
        raise OSError(f"[{label}] file too small (<1KB), likely corrupt: {path}")
    try:
        with h5py.File(path, "r") as f:
            _ = list(f.keys())
    except Exception as e:
        raise OSError(
            f"[{label}] not a valid HDF5/.h5ad (or corrupted): {path}\nOriginal error: {e}"
        )

def ensure_csr_float32(adata: ad.AnnData):
    X = adata.X
    if sp.issparse(X):
        if not sp.isspmatrix_csr(X):
            X = X.tocsr()
        if X.dtype != np.float32:
            X = X.astype(np.float32)
        adata.X = X
    else:
        adata.X = sp.csr_matrix(np.asarray(X, dtype=np.float32))

def ensure_csc_float32(adata: ad.AnnData):
    X = adata.X
    if sp.issparse(X):
        if not sp.isspmatrix_csc(X):
            X = X.tocsc()
        if X.dtype != np.float32:
            X = X.astype(np.float32)
        adata.X = X
    else:
        adata.X = sp.csc_matrix(np.asarray(X, dtype=np.float32))

def safe_hvg(adata: ad.AnnData, n_top: int = 3000):
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat_v3")
    except Exception:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="cell_ranger")

def _norm(s: str) -> str:
    return re.sub(r"[\s_]+", "", str(s).strip().lower())

def find_obs_column(adata: ad.AnnData, candidates):
    norm_map = {_norm(c): c for c in candidates}
    for col in adata.obs.columns:
        if _norm(col) in norm_map:
            return col
    return None

def get_coords(st_adata: ad.AnnData) -> pd.DataFrame:
    if "spatial" in st_adata.obsm_keys():
        arr = st_adata.obsm["spatial"]
        if arr is not None and arr.shape[0] == st_adata.n_obs and arr.shape[1] >= 2:
            return pd.DataFrame(arr[:, :2], index=st_adata.obs_names, columns=["x", "y"]).astype("float32")
    if {"x", "y"}.issubset(st_adata.obs.columns):
        df = st_adata.obs[["x", "y"]].copy()
        df.index = st_adata.obs_names
        return df.astype("float32")
    if {"array_row", "array_col"}.issubset(st_adata.obs.columns):
        df = st_adata.obs.rename(columns={"array_row": "x", "array_col": "y"})[["x", "y"]].copy()
        df.index = st_adata.obs_names
        return df.astype("float32")
    raise KeyError("No spatial coordinates found (need obsm['spatial'] or obs x/y or array_row/array_col).")

def pick_symbol_column(var: pd.DataFrame):
    for c in ["real_gene_name", "gene_symbols", "gene_symbol", "symbol", "feature_name", "gene_name", "name"]:
        if c in var.columns:
            return c
    return None

def gene_id_diagnostics(st_adata: ad.AnnData, sc_adata: ad.AnnData, outdir: str):
    def ens_like(vn):
        vn = list(map(str, vn[:200]))
        return sum(x.startswith("ENSG") for x in vn), len(vn)
    st_ens, st_n = ens_like(st_adata.var_names)
    sc_ens, sc_n = ens_like(sc_adata.var_names)
    st_sym_col = pick_symbol_column(st_adata.var)
    sc_sym_col = pick_symbol_column(sc_adata.var)
    lines = [
        f"ST shape={st_adata.shape}",
        f"SC shape={sc_adata.shape}",
        f"ST var_names example={list(map(str, st_adata.var_names[:10]))}",
        f"SC var_names example={list(map(str, sc_adata.var_names[:10]))}",
        f"ST ENSG-like(first{st_n})={st_ens}/{st_n}",
        f"SC ENSG-like(first{sc_n})={sc_ens}/{sc_n}",
        f"ST symbol column candidate={st_sym_col}",
        f"SC symbol column candidate={sc_sym_col}",
        f"ST var columns(first30)={list(st_adata.var.columns)[:30]}",
        f"SC var columns(first30)={list(sc_adata.var.columns)[:30]}",
        f"SC obs columns(first60)={list(sc_adata.obs.columns)[:60]}",
    ]
    with open(os.path.join(outdir, "gene_id_diagnostics.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")

def harmonize_to_symbols(st_adata: ad.AnnData, sc_adata: ad.AnnData, st_symbol_col: str = "real_gene_name", sc_symbol_col: str | None = None):
    info = {}
    st = st_adata.copy()
    scd = sc_adata.copy()
    st.var["original_var_names"] = st.var_names.astype(str).values
    scd.var["original_var_names"] = scd.var_names.astype(str).values
    st_col = st_symbol_col if st_symbol_col in st.var.columns else pick_symbol_column(st.var)
    sc_col = sc_symbol_col if (sc_symbol_col and sc_symbol_col in scd.var.columns) else pick_symbol_column(scd.var)
    info["st_symbol_col_used"] = st_col
    info["sc_symbol_col_used"] = sc_col
    if st_col is not None:
        sym = st.var[st_col].astype(object)
        sym = sym.where(~pd.isna(sym), other=st.var["original_var_names"].astype(object))
        sym = sym.astype(str)
        st.var["gene_symbol"] = sym.astype(object)
        st.var_names = pd.Index(st.var["gene_symbol"].astype(str))
        st.var_names_make_unique()
    else:
        st.var["gene_symbol"] = st.var_names.astype(str).astype(object)
    if sc_col is not None:
        sym = scd.var[sc_col].astype(object)
        sym = sym.where(~pd.isna(sym), other=scd.var["original_var_names"].astype(object))
        sym = sym.astype(str)
        scd.var["gene_symbol"] = sym.astype(object)
        scd.var_names = pd.Index(scd.var["gene_symbol"].astype(str))
        scd.var_names_make_unique()
    else:
        scd.var["gene_symbol"] = scd.var_names.astype(str).astype(object)
    st.var.index.name = "gene"
    scd.var.index.name = "gene"
    return st, scd, info

def filter_sc_celltypes_by_min_cells(sc_adata: ad.AnnData, celltype_col: str, min_cells: int):
    ct = sc_adata.obs[celltype_col].astype(str)
    counts = ct.value_counts(dropna=False)
    keep_types = counts[counts >= min_cells].index.tolist()
    drop_types = counts[counts <  min_cells].index.tolist()
    keep_mask = ct.isin(keep_types).to_numpy()
    sc_filt = sc_adata[keep_mask].copy()
    return sc_filt, counts, keep_types, drop_types

def build_marker_gene_list(sc_adata: ad.AnnData, celltype_col: str, n_top_markers: int = 50, max_cells_per_type: int = 5000, seed: int = 0) -> list[str]:
    rng = np.random.default_rng(seed)
    obs = sc_adata.obs[celltype_col].astype(str)
    keep_idx = []
    for ct, idx in obs.groupby(obs).groups.items():
        idx = np.asarray(list(idx))
        if idx.size > max_cells_per_type:
            idx = rng.choice(idx, size=max_cells_per_type, replace=False)
        keep_idx.append(idx)
    keep_idx = np.concatenate(keep_idx)
    sc_sub = sc_adata[keep_idx].copy()
    ensure_csr_float32(sc_sub)
    try:
        xmax = float(sc_sub.X.max())
    except Exception:
        xmax = 0.0
    if xmax > 50:
        sc.pp.normalize_total(sc_sub, target_sum=1e4)
        sc.pp.log1p(sc_sub)
    safe_hvg(sc_sub, n_top=min(3000, sc_sub.n_vars))
    sc_sub = sc_sub[:, sc_sub.var["highly_variable"]].copy()
    sc.tl.rank_genes_groups(sc_sub, groupby=celltype_col, method="wilcoxon", n_genes=min(200, sc_sub.n_vars), use_raw=False)
    mdf = sc.get.rank_genes_groups_df(sc_sub, group=None).copy()
    mdf = (mdf.sort_values(["group", "pvals_adj", "logfoldchanges"], ascending=[True, True, False]).groupby("group").head(n_top_markers))
    genes = mdf["names"].astype(str).dropna().unique().tolist()
    return genes

def _nnls_one(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    x, _ = nnls(A, b)
    s = x.sum()
    if s > 0:
        x = x / s
    return x.astype(np.float32, copy=False)

def build_reference_matrix(sc_adata: ad.AnnData, celltype_col: str, genes: list[str]):
    sc_ref = sc_adata[:, genes].copy()
    ensure_csr_float32(sc_ref)
    try:
        xmax = float(sc_ref.X.max())
    except Exception:
        xmax = 0.0
    if xmax > 50:
        sc.pp.normalize_total(sc_ref, target_sum=1e4)
        sc.pp.log1p(sc_ref)
    celltypes = sc_ref.obs[celltype_col].astype(str)
    ct_names = sorted(celltypes.unique().tolist())
    A = np.zeros((sc_ref.n_vars, len(ct_names)), dtype=np.float32)
    X = sc_ref.X
    for j, ct in enumerate(ct_names):
        idx = np.where(celltypes.to_numpy() == ct)[0]
        A[:, j] = np.asarray(X[idx, :].mean(axis=0)).ravel().astype(np.float32)
    keep = (A.sum(axis=0) > 1e-8)
    A = A[:, keep]
    ct_names = [c for c, ok in zip(ct_names, keep) if ok]
    g = np.sqrt((A * A).mean(axis=1))
    g[g == 0] = 1.0
    A = (A.T / g).T
    return A, ct_names, g

def run_nnls_deconvolution(st_adata: ad.AnnData, sc_adata: ad.AnnData, celltype_col: str, outdir: str, n_top_markers: int, max_cells_per_type: int, n_jobs: int, seed: int, min_overlap_genes: int) -> pd.DataFrame:
    print("[NNLS] building marker genes from scRNA...")
    markers = build_marker_gene_list(sc_adata=sc_adata, celltype_col=celltype_col, n_top_markers=n_top_markers, max_cells_per_type=max_cells_per_type, seed=seed)
    pd.DataFrame({"marker_gene": markers}).to_csv(os.path.join(outdir, "marker_genes.csv"), index=False)
    print(f"[NNLS] markers: {len(markers)}")
    marker_set = set(map(str, markers))
    genes_use = sorted(list(marker_set & set(map(str, st_adata.var_names)) & set(map(str, sc_adata.var_names))))
    print(f"[NNLS] overlap genes (markers ∩ ST ∩ SC): {len(genes_use)}")
    if len(genes_use) < min_overlap_genes:
        with open(os.path.join(outdir, "nnls_overlap_debug.txt"), "w") as f:
            f.write(f"markers={len(markers)}\n")
            f.write(f"overlap={len(genes_use)}\n\n")
            f.write("Example ST var_names:\n" + "\n".join(list(map(str, st_adata.var_names[:40]))) + "\n\n")
            f.write("Example SC var_names:\n" + "\n".join(list(map(str, sc_adata.var_names[:40]))) + "\n\n")
            f.write("Example markers:\n" + "\n".join(list(map(str, markers[:40]))) + "\n")
        raise ValueError(f"[NNLS] too few marker genes overlap between sc and st: {len(genes_use)}")
    print("[NNLS] building reference profiles...")
    A, ct_names, g = build_reference_matrix(sc_adata, celltype_col, genes_use)
    print(f"[NNLS] A matrix: genes={A.shape[0]} celltypes={A.shape[1]}")
    st_ref = st_adata[:, genes_use].copy()
    ensure_csr_float32(st_ref)
    try:
        xmax = float(st_ref.X.max())
    except Exception:
        xmax = 0.0
    if xmax > 50:
        sc.pp.normalize_total(st_ref, target_sum=1e4)
        sc.pp.log1p(st_ref)
    Xst = st_ref.X
    spot_names = st_ref.obs_names.to_numpy()
    print(f"[NNLS] running NNLS on {len(spot_names)} spots (n_jobs={n_jobs})...")
    def solve_i(i: int):
        b = Xst[i, :].toarray().ravel().astype(np.float32, copy=False)
        b = b / g
        return _nnls_one(A, b)
    frac = Parallel(n_jobs=n_jobs, prefer="processes", verbose=0)(delayed(solve_i)(i) for i in range(Xst.shape[0]))
    frac = np.vstack(frac).astype(np.float32, copy=False)
    frac_df = pd.DataFrame(frac, index=spot_names, columns=ct_names)
    frac_df.to_csv(os.path.join(outdir, "celltype_fractions_nnls.csv"))
    st_adata.obsm["celltype_fractions_nnls"] = frac_df.to_numpy(dtype=np.float32)
    st_adata.uns["celltype_fraction_cols_nnls"] = ct_names
    st_adata.obs["dominant_celltype"] = frac_df.idxmax(axis=1).astype("category")
    print("[NNLS] wrote celltype_fractions_nnls.csv")
    return frac_df

def explicit_aggregation_map_cellxgene_endometrium():
    m = {}
    m.update({
        "epithelial cell": "Endometrial",
        "luminal endometrial epithelial cell": "Endometrial",
        "glandular secretory epithelial cell": "Endometrial",
        "ciliated epithelial cell": "Endometrial",
        "myoepithelial cell": "Endometrial",
    })
    m.update({
        "stromal cell": "Stromal",
        "stromal cell of endometrium": "Stromal",
        "fibroblast": "Stromal",
        "fibrocyte": "Stromal",
        "fibro/adipogenic progenitor cell": "Stromal",
    })
    m.update({
        "blood vessel endothelial cell": "Vascular",
        "capillary endothelial cell": "Vascular",
        "endothelial cell of artery": "Vascular",
        "endothelial cell of venule": "Vascular",
        "endothelial cell of lymphatic vessel": "Vascular",
        "endothelial cell of uterus": "Vascular",
    })
    m.update({
        "perivascular cell": "Perivascular",
        "uterine smooth muscle cell": "Perivascular",
    })
    m.update({
        "B cell": "B_cell",
        "plasma cell": "B_cell",
        "natural killer cell": "NK",
        "innate lymphoid cell": "NK",
        "macrophage": "Macrophage",
        "monocyte": "Myeloid",
        "leukocyte": "Myeloid",
        "hematopoietic precursor cell": "Myeloid",
        "memory T cell": "T_cell",
        "regulatory T cell": "T_cell",
        "CD4-positive, alpha-beta memory T cell": "T_cell",
        "CD8-positive, alpha-beta memory T cell": "T_cell",
    })
    m.update({
        "dendritic cell": "DC",
        "plasmacytoid dendritic cell": "DC",
    })
    m.update({"mast cell": "Other"})
    return m

def aggregate_fraction_df_explicit(frac_df: pd.DataFrame, agg_map: dict, fallback="Other") -> pd.DataFrame:
    cols = [str(c) for c in frac_df.columns]
    mapped = [agg_map.get(c, fallback) for c in cols]
    out = frac_df.copy()
    out.columns = mapped
    out = out.groupby(level=0, axis=1).sum()
    s = out.sum(axis=1).replace(0, np.nan)
    out = out.div(s, axis=0).fillna(0.0)
    return out

def apply_aggregation(st_adata, frac_df: pd.DataFrame, outdir: str):
    agg_map = explicit_aggregation_map_cellxgene_endometrium()
    present_cols = [str(c) for c in frac_df.columns]
    map_df = pd.DataFrame({"original": present_cols, "aggregated": [agg_map.get(c, "Other") for c in present_cols]})
    map_df.to_csv(os.path.join(outdir, "celltype_aggregation_map.csv"), index=False)
    frac_agg = aggregate_fraction_df_explicit(frac_df, agg_map, fallback="Other")
    frac_agg.to_csv(os.path.join(outdir, "celltype_fractions_nnls_aggregated.csv"))
    st_adata.obsm["celltype_fractions_nnls_aggregated"] = frac_agg.to_numpy(dtype=np.float32)
    st_adata.uns["celltype_fraction_cols_nnls_aggregated"] = frac_agg.columns.astype(str).tolist()
    st_adata.obsm["celltype_fractions_nnls"] = st_adata.obsm["celltype_fractions_nnls_aggregated"]
    st_adata.uns["celltype_fraction_cols_nnls"] = st_adata.uns["celltype_fraction_cols_nnls_aggregated"]
    st_adata.obs["dominant_celltype"] = frac_agg.idxmax(axis=1).astype("category")
    print(f"[Aggregate] Applied aggregation: {len(frac_df.columns)} -> {len(frac_agg.columns)} columns")
    return st_adata, frac_agg

def compute_regions_leiden(st_adata: ad.AnnData, coords_df: pd.DataFrame, resolution: float, hvg_n_top: int, n_pcs: int, spatial_n_neigh: int, fallback_n_neighbors: int, use_spatial_neighbors: bool):
    try:
        import squidpy as sq
        HAVE_SQUIDPY = True
    except Exception:
        HAVE_SQUIDPY = False
    st_for_regions = st_adata.copy()
    ensure_csr_float32(st_for_regions)
    try:
        xmax = float(st_for_regions.X.max())
    except Exception:
        xmax = 0.0
    if xmax > 50:
        sc.pp.normalize_total(st_for_regions, target_sum=1e4)
        sc.pp.log1p(st_for_regions)
    safe_hvg(st_for_regions, n_top=min(hvg_n_top, st_for_regions.n_vars))
    st_for_regions = st_for_regions[:, st_for_regions.var["highly_variable"]].copy()
    sc.pp.scale(st_for_regions, max_value=10)
    sc.tl.pca(st_for_regions, n_comps=n_pcs)
    st_for_regions.obsm["spatial"] = coords_df.loc[st_for_regions.obs_names, ["x", "y"]].to_numpy()
    if HAVE_SQUIDPY and use_spatial_neighbors:
        sq.gr.spatial_neighbors(st_for_regions, coord_type="generic", n_neigh=spatial_n_neigh)
        sc.tl.leiden(st_for_regions, adjacency=st_for_regions.obsp["spatial_connectivities"], key_added="region", resolution=resolution)
        method = f"squidpy_spatial_neighbors(k={spatial_n_neigh}) + leiden(res={resolution})"
    else:
        sc.pp.neighbors(st_for_regions, n_neighbors=fallback_n_neighbors, n_pcs=min(30, n_pcs))
        sc.tl.leiden(st_for_regions, key_added="region", resolution=resolution)
        method = f"scanpy_neighbors(n_neighbors={fallback_n_neighbors}) + leiden(res={resolution})"
    return st_for_regions, method

def _to_builtin(obj):
    """Recursively convert numpy scalars/arrays in .uns to HDF5-friendly Python types."""
    if isinstance(obj, (np.generic,)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, dict):
        return {str(k): _to_builtin(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [ _to_builtin(v) for v in obj ]
    return obj

def make_h5ad_write_safe(adata: ad.AnnData) -> ad.AnnData:
    if adata.var.index.name is None:
        adata.var.index.name = "gene"
    if adata.var.index.name in adata.var.columns:
        adata.var.index.name = "gene"
    if adata.obs.index.name is None:
        adata.obs.index.name = "obs"
    adata.var_names = pd.Index(adata.var_names.astype(str))
    adata.obs_names = pd.Index(adata.obs_names.astype(str))
    def to_object_str(s: pd.Series) -> pd.Series:
        return s.where(s.isna(), other=s.astype(str)).astype(object)
    for df in (adata.var, adata.obs):
        for c in list(df.columns):
            if df[c].dtype == "object" or pd.api.types.is_string_dtype(df[c]):
                df[c] = to_object_str(df[c])
            if isinstance(df[c].dtype, pd.CategoricalDtype):
                df[c] = df[c].astype("category")
    # sanitize .uns
    try:
        adata.uns = _to_builtin(dict(adata.uns))
    except Exception:
        pass
    return adata

def write_h5ad_atomic(adata: ad.AnnData, out_path: str):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # temp file in same dir (important on NFS for atomic replace)
    fd, tmp = tempfile.mkstemp(prefix=out_path.name + ".", suffix=".tmp", dir=str(out_path.parent))
    os.close(fd)
    try:
        adata.write(tmp)
        os.replace(tmp, out_path)
    finally:
        try:
            if os.path.exists(tmp):
                os.remove(tmp)
        except Exception:
            pass

def main():
    ap = argparse.ArgumentParser(description="Spatial NNLS deconvolution + regions (cellxgene scRNA reference).")
    ap.add_argument("--st", required=True, help="Path to spatial .h5ad")
    ap.add_argument("--sc", required=True, help="Path to scRNA reference .h5ad")
    ap.add_argument("--out_base", default="/data/scripts/results", help="Base output folder")
    ap.add_argument("--sample", required=True, help="Sample name -> creates <out_base>/<sample>/")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing sample folder")
    ap.add_argument("--celltype_col", default=None, help="sc_adata.obs column for cell type labels (default auto-detect)")
    ap.add_argument("--min_cells_per_celltype", type=int, default=20, help="Exclude scRNA cell types with fewer than this many cells")
    ap.add_argument("--harmonize_genes", choices=["none", "symbols"], default="symbols", help="Gene harmonization across ST/SC")
    ap.add_argument("--st_symbol_col", default="real_gene_name", help="ST var column for symbols")
    ap.add_argument("--sc_symbol_col", default=None, help="SC var column for symbols (auto-detect if None)")
    ap.add_argument("--n_top_markers", type=int, default=50, help="Markers per celltype used for NNLS")
    ap.add_argument("--max_cells_per_type", type=int, default=5000, help="SC downsample per celltype for marker finding")
    ap.add_argument("--min_overlap_genes", type=int, default=50, help="Minimum overlap genes required")
    ap.add_argument("--n_jobs", type=int, default=8, help="Parallel jobs for NNLS")
    ap.add_argument("--seed", type=int, default=0, help="Random seed")
    ap.add_argument("--aggregated", action="store_true", help="Aggregate into macro-groups")
    ap.add_argument("--resolution", type=float, default=0.6, help="Leiden resolution")
    ap.add_argument("--hvg_n_top", type=int, default=3000, help="HVGs for region calling")
    ap.add_argument("--n_pcs", type=int, default=50, help="PCA components for region calling")
    ap.add_argument("--spatial_n_neigh", type=int, default=8, help="Squidpy spatial neighbors")
    ap.add_argument("--fallback_n_neighbors", type=int, default=15, help="Scanpy neighbors if squidpy not used")
    ap.add_argument("--no_spatial_neighbors", action="store_true", help="Disable squidpy spatial neighbors")
    args = ap.parse_args()

    sample = sanitize_name(args.sample)
    outdir = os.path.join(args.out_base, sample)

    if os.path.exists(outdir) and os.listdir(outdir) and not args.overwrite:
        print(f"[ERROR] Output folder exists and is not empty: {outdir}")
        print("Use --overwrite or choose a different --sample")
        sys.exit(2)

    safe_mkdir(outdir)
    print("[IO] outdir:", outdir)

    assert_h5ad_valid(args.st, "ST")
    assert_h5ad_valid(args.sc, "SC")

    print("[Load] ST:", args.st)
    st_adata = sc.read_h5ad(args.st); st_adata.var_names_make_unique()
    print("[Load] SC:", args.sc)
    sc_adata = sc.read_h5ad(args.sc); sc_adata.var_names_make_unique()

    ensure_csc_float32(st_adata)
    ensure_csc_float32(sc_adata)

    coords_df = get_coords(st_adata)
    coords_df.to_csv(os.path.join(outdir, "coords_df.csv"))

    gene_id_diagnostics(st_adata, sc_adata, outdir)
    print("[Genes] wrote gene_id_diagnostics.txt")

    harmonize_info = {"st_symbol_col_used": None, "sc_symbol_col_used": None}
    if args.harmonize_genes == "symbols":
        st_adata, sc_adata, harmonize_info = harmonize_to_symbols(
            st_adata=st_adata, sc_adata=sc_adata, st_symbol_col=args.st_symbol_col, sc_symbol_col=args.sc_symbol_col
        )
        print("[Genes] harmonized to symbols:", "ST col=", harmonize_info["st_symbol_col_used"], "SC col=", harmonize_info["sc_symbol_col_used"])
    else:
        st_adata.var["original_var_names"] = st_adata.var_names.astype(str).astype(object)
        sc_adata.var["original_var_names"] = sc_adata.var_names.astype(str).astype(object)
        st_adata.var.index.name = "gene"
        sc_adata.var.index.name = "gene"

    if args.celltype_col is None:
        celltype_aliases = [
            "cell_type", "celltype", "cell_type_label", "celltype_major",
            "cell_annotation", "annotation", "labels", "label",
            "cell_class", "cell_subtype"
        ]
        ctcol = find_obs_column(sc_adata, celltype_aliases)
        if ctcol is None:
            raise KeyError(
                "Cannot auto-detect celltype column in sc_adata.obs. "
                f"Available columns: {list(sc_adata.obs.columns)[:120]}"
            )
        args.celltype_col = ctcol
    print("[SC] using celltype_col:", args.celltype_col)

    sc_adata, ct_counts, keep_types, drop_types = filter_sc_celltypes_by_min_cells(
        sc_adata, args.celltype_col, args.min_cells_per_celltype
    )
    ct_counts.to_csv(os.path.join(outdir, "sc_celltype_counts_before_filter.csv"))
    pd.DataFrame({"kept_celltypes": keep_types}).to_csv(os.path.join(outdir, "sc_celltypes_kept.csv"), index=False)
    pd.DataFrame({"dropped_celltypes": drop_types}).to_csv(os.path.join(outdir, "sc_celltypes_dropped.csv"), index=False)
    print(f"[SC] kept cell types: {len(keep_types)}; dropped (<{args.min_cells_per_celltype}): {len(drop_types)}; cells remaining: {sc_adata.n_obs}")
    if sc_adata.n_obs == 0 or sc_adata.obs[args.celltype_col].nunique() < 2:
        raise ValueError("[SC] After filtering, too few cells/cell types remain for deconvolution.")

    frac_csv = os.path.join(outdir, "celltype_fractions_nnls.csv")
    if os.path.exists(frac_csv) and (not args.overwrite):
        print("[NNLS] found existing fractions:", frac_csv)
        frac_df = pd.read_csv(frac_csv, index_col=0).reindex(st_adata.obs_names).fillna(0.0)
        st_adata.obsm["celltype_fractions_nnls"] = frac_df.to_numpy(dtype=np.float32)
        st_adata.uns["celltype_fraction_cols_nnls"] = frac_df.columns.astype(str).tolist()
        st_adata.obs["dominant_celltype"] = frac_df.idxmax(axis=1).astype("category")
    else:
        frac_df = run_nnls_deconvolution(
            st_adata=st_adata, sc_adata=sc_adata, celltype_col=args.celltype_col, outdir=outdir,
            n_top_markers=args.n_top_markers, max_cells_per_type=args.max_cells_per_type,
            n_jobs=args.n_jobs, seed=args.seed, min_overlap_genes=args.min_overlap_genes
        )

    if args.aggregated:
        st_adata, frac_df = apply_aggregation(st_adata, frac_df, outdir)

    # --- export #1: NNLS outputs only (dominant cell type + fractions), no regions ---
    # Useful for downstream plotting where you do not want region annotations.
    out_h5ad_nnls_only = os.path.join(outdir, "mirxes_celltypes_nnls_only.h5ad")
    st_nnls_only = make_h5ad_write_safe(st_adata.copy())
    write_h5ad_atomic(st_nnls_only, out_h5ad_nnls_only)
    print("[Export] wrote NNLS-only h5ad:", out_h5ad_nnls_only)

    print("[Regions] computing regions...")
    st_for_regions, region_method = compute_regions_leiden(
        st_adata=st_adata, coords_df=coords_df, resolution=args.resolution,
        hvg_n_top=args.hvg_n_top, n_pcs=args.n_pcs, spatial_n_neigh=args.spatial_n_neigh,
        fallback_n_neighbors=args.fallback_n_neighbors, use_spatial_neighbors=(not args.no_spatial_neighbors),
    )
    st_adata.obs["region"] = st_for_regions.obs["region"].reindex(st_adata.obs_names).astype("category")

    region_means = frac_df.join(st_adata.obs["region"]).groupby("region").mean()
    region_label = region_means.idxmax(axis=1).to_dict()
    st_adata.obs["region_celltype_label"] = st_adata.obs["region"].map(region_label).astype("category")
    region_means.to_csv(os.path.join(outdir, "region_mean_composition.csv"))
    if args.aggregated:
        region_means.to_csv(os.path.join(outdir, "region_mean_composition_aggregated.csv"))
    st_adata.obs[["region", "region_celltype_label", "dominant_celltype"]].to_csv(os.path.join(outdir, "spot_region_labels.csv"))

    with open(os.path.join(outdir, "run_params.txt"), "w") as f:
        f.write(f"st={args.st}\n")
        f.write(f"sc={args.sc}\n")
        f.write(f"sample={sample}\n")
        f.write(f"outdir={outdir}\n")
        f.write(f"celltype_col={args.celltype_col}\n")
        f.write(f"min_cells_per_celltype={args.min_cells_per_celltype}\n")
        f.write(f"harmonize_genes={args.harmonize_genes}\n")
        f.write(f"st_symbol_col={args.st_symbol_col}\n")
        f.write(f"sc_symbol_col={args.sc_symbol_col}\n")
        f.write(f"st_symbol_col_used={harmonize_info.get('st_symbol_col_used')}\n")
        f.write(f"sc_symbol_col_used={harmonize_info.get('sc_symbol_col_used')}\n")
        f.write(f"n_top_markers={args.n_top_markers}\n")
        f.write(f"max_cells_per_type={args.max_cells_per_type}\n")
        f.write(f"min_overlap_genes={args.min_overlap_genes}\n")
        f.write(f"n_jobs={args.n_jobs}\n")
        f.write(f"seed={args.seed}\n")
        f.write(f"aggregated={args.aggregated}\n")
        f.write(f"region_method={region_method}\n")
        f.write(f"resolution={args.resolution}\n")
        f.write(f"hvg_n_top={args.hvg_n_top}\n")
        f.write(f"n_pcs={args.n_pcs}\n")
        f.write(f"spatial_n_neigh={args.spatial_n_neigh}\n")
        f.write(f"fallback_n_neighbors={args.fallback_n_neighbors}\n")

    # --- robust H5AD export ---
    out_h5ad = os.path.join(outdir, "mirxes_regions_with_celltypes_nnls.h5ad")
    st_adata = make_h5ad_write_safe(st_adata)
    write_h5ad_atomic(st_adata, out_h5ad)

    print("[Done] wrote:", out_h5ad)
    print("[Done] outputs in:", outdir)
    print("[Regions] method:", region_method)

if __name__ == "__main__":
    main()
