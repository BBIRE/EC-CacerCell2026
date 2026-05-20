#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl

from scipy import sparse
from scipy.stats import mannwhitneyu
from sklearn.neighbors import kneighbors_graph
from matplotlib.colors import to_rgba

mpl.rcParams["savefig.bbox"] = "standard"

SIGNATURES = {
    "macrophage_phagocytosis_negative_regulation": [
        "PTPN6", "PTPN11", "FCGR2B", "LILRB1", "LILRB2", "CD300A",
        "TGFB1", "IL10", "SOCS1", "SOCS3",
    ]
}

CELLTYPE_COLORS = {
    "B_cell": "#69B3A2", "B cell": "#69B3A2", "Bcells": "#69B3A2", "B cells": "#69B3A2",
    "Myeloid": "#E7B07A", "Dendritic": "#E7B07A", "Dendritic cells": "#E7B07A",
    "Endometrial": "#9d8189",
    "Macrophage": "#E6A43A", "Macrophages": "#E6A43A",
    "NK": "#7FB77E", "NK cell": "#7FB77E", "NK cells": "#7FB77E",
    "T_cell": "#D95F5F", "Tcell": "#D95F5F", "T cell": "#D95F5F", "T cells": "#D95F5F",
    "T cells CD4": "#B08CCB", "T cells CD8": "#D95F5F", "CD4 T cell": "#B08CCB", "CD8 T cell": "#D95F5F",
    "Stromal": "#fdf0d5", "Vascular": "#952406", "Perivascular": "#f9b5ac",
    "Others": "#D98C8C", "Other": "#d9d9d9",
}

ORIGIN_COLORS = {
    "endometrial": CELLTYPE_COLORS["Endometrial"],
    "macrophage": CELLTYPE_COLORS["Macrophage"],
    "NK": CELLTYPE_COLORS["NK"],
    "Tcell": CELLTYPE_COLORS["T_cell"],
}


def origin_rgba(origin, positive=True):
    return to_rgba(ORIGIN_COLORS[origin], alpha=1.0 if positive else 0.35)


def normalize_celltype_label(label):
    label = str(label)
    if label in CELLTYPE_COLORS:
        return label
    ll = label.lower()
    if "endo" in ll: return "Endometrial"
    if "macro" in ll: return "Macrophage"
    if "nk" in ll or "natural killer" in ll: return "NK"
    if "cd8" in ll: return "T cells CD8"
    if "cd4" in ll: return "T cells CD4"
    if "t_cell" in ll or "tcell" in ll or "t cell" in ll: return "T_cell"
    if "strom" in ll: return "Stromal"
    if "vascular" in ll and "peri" in ll: return "Perivascular"
    if "vascular" in ll: return "Vascular"
    if "myeloid" in ll: return "Myeloid"
    if "dend" in ll: return "Dendritic"
    if "bcell" in ll or "b cell" in ll or "b_cell" in ll: return "B_cell"
    return label


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def to_dense(x):
    return x.toarray() if sparse.issparse(x) else np.asarray(x)


def get_xy(adata, spatial_key="spatial"):
    if spatial_key in adata.obsm:
        xy = np.asarray(adata.obsm[spatial_key])
        if xy.ndim != 2 or xy.shape[1] != 2:
            raise ValueError(f"adata.obsm['{spatial_key}'] must be n_obs x 2.")
        return xy.astype(float)
    for xcol, ycol in [("x", "y"), ("pxl_col_in_fullres", "pxl_row_in_fullres"), ("array_col", "array_row")]:
        if xcol in adata.obs.columns and ycol in adata.obs.columns:
            return np.asarray(adata.obs[[xcol, ycol]], dtype=float)
    raise KeyError("No spatial coordinates found.")


def get_gene_vector(adata, gene, layer=None, use_raw=False):
    if use_raw and adata.raw is not None:
        if gene not in adata.raw.var_names:
            return np.full(adata.n_obs, np.nan)
        x = adata.raw[:, gene].X
    else:
        if gene not in adata.var_names:
            return np.full(adata.n_obs, np.nan)
        x = adata[:, gene].layers[layer] if layer is not None and layer in adata.layers else adata[:, gene].X
    return to_dense(x).reshape(-1).astype(float)


def get_total_transcripts(adata, layer=None, use_raw=False):
    if "total_counts" in adata.obs.columns:
        return pd.to_numeric(adata.obs["total_counts"], errors="coerce").to_numpy(float)
    if use_raw and adata.raw is not None:
        X = adata.raw.X
    elif layer is not None and layer in adata.layers:
        X = adata.layers[layer]
    else:
        X = adata.X
    return np.asarray(X.sum(axis=1)).reshape(-1).astype(float) if sparse.issparse(X) else np.asarray(X.sum(axis=1)).reshape(-1).astype(float)


def safe_pct(n, d):
    return 100.0 * n / d if d else np.nan


def significance_stars(p):
    if p is None or not np.isfinite(p): return "n.s."
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "n.s."


def build_spatial_graph(xy, n_neighbors=8):
    n_neighbors = min(int(n_neighbors), xy.shape[0] - 1)
    graph = kneighbors_graph(xy, n_neighbors=n_neighbors, mode="connectivity", include_self=False)
    graph = graph.maximum(graph.T)
    graph.setdiag(0)
    graph.eliminate_zeros()
    return graph.tocsr()


def reachable_hit_within_hop(graph, start_idx, allowed_mask, max_hop):
    graph = graph.tocsr()
    allowed_mask = np.asarray(allowed_mask, dtype=bool)
    hit = np.zeros(graph.shape[0], dtype=bool)
    for start in np.asarray(start_idx, dtype=int):
        visited = {int(start)}
        frontier = np.array([start], dtype=int)
        found = False
        for _ in range(1, int(max_hop) + 1):
            if frontier.size == 0: break
            neigh = np.unique(graph[frontier].indices)
            neigh = np.array([x for x in neigh if int(x) not in visited], dtype=int)
            if neigh.size == 0: break
            visited.update([int(x) for x in neigh])
            if np.any(allowed_mask[neigh]):
                found = True
                break
            frontier = neigh
        hit[start] = found
    return hit


def score_signature(adata, genes, score_name, layer=None, use_raw=False):
    var_names = adata.raw.var_names if (use_raw and adata.raw is not None) else adata.var_names
    present = [g for g in genes if g in var_names]
    if len(present) == 0:
        adata.obs[score_name] = np.nan
        return []
    if use_raw and adata.raw is not None:
        X = adata.raw[:, present].X
    elif layer is not None and layer in adata.layers:
        X = adata[:, present].layers[layer]
    else:
        X = adata[:, present].X
    X = to_dense(X).astype(float)
    gene_mean = np.nanmean(X, axis=0)
    gene_sd = np.nanstd(X, axis=0)
    gene_sd[gene_sd == 0] = 1.0
    adata.obs[score_name] = np.nanmean((X - gene_mean) / gene_sd, axis=1)
    return present


def add_expression_and_signature_calls(adata, layer=None, use_raw=False, transcript_threshold=0.0):
    for gene, name in [("EPCAM", "EPCAM"), ("CD47", "CD47"), ("SIRPA", "SIRPA"), ("NECTIN2", "NECTIN2"), ("TIGIT", "TIGIT")]:
        vec = get_gene_vector(adata, gene, layer=layer, use_raw=use_raw)
        adata.obs[f"{name}_expr"] = vec
        adata.obs[f"{name}_pos"] = vec > transcript_threshold
    rows = []
    for sig, genes in SIGNATURES.items():
        present = score_signature(adata, genes, sig, layer=layer, use_raw=use_raw)
        rows.extend({"signature": sig, "gene": g} for g in present)
    adata.obs["total_transcripts_pipeline"] = get_total_transcripts(adata, layer=layer, use_raw=use_raw)
    return pd.DataFrame(rows)


def build_cell_masks(labels, endometrial_label, macrophage_label):
    labels_s = labels.astype(str)
    labels_l = labels_s.str.lower()
    labels_norm = labels_l.str.replace("_", " ", regex=False).str.replace("-", " ", regex=False).str.replace("+", "", regex=False).str.strip()
    endometrial_mask = labels_s.eq(endometrial_label).to_numpy()
    macrophage_mask = labels_s.eq(macrophage_label).to_numpy()
    nk_mask = (labels_s.eq("NK") | labels_s.eq("NK_cell") | labels_s.eq("NK_cells") | labels_norm.eq("nk") | labels_norm.eq("nk cell") | labels_norm.eq("nk cells") | labels_norm.eq("natural killer cell") | labels_norm.eq("natural killer cells")).to_numpy()
    tcell_mask = (labels_s.eq("T_cell") | labels_s.eq("T_cells") | labels_s.eq("Tcell") | labels_s.eq("Tcells") | labels_l.eq("t_cell") | labels_l.eq("t_cells") | labels_l.eq("tcell") | labels_l.eq("tcells") | labels_norm.eq("t cell") | labels_norm.eq("t cells") | labels_norm.eq("cd4 t cell") | labels_norm.eq("cd8 t cell") | labels_norm.eq("cd4 t cells") | labels_norm.eq("cd8 t cells") | labels_norm.eq("cd8 positive t cell")).to_numpy()
    return endometrial_mask, macrophage_mask, nk_mask, tcell_mask


def composition_table(adata, label_col):
    out = adata.obs[label_col].astype(str).value_counts().rename_axis("cell_type").reset_index(name="n")
    out["cell_type_norm"] = out["cell_type"].map(normalize_celltype_label)
    out["pct"] = 100.0 * out["n"] / out["n"].sum()
    return out


def two_by_two_table(adata, source_mask, gene_a_col, gene_b_col, a_name, b_name):
    obs = adata.obs
    a_pos = obs[f"{gene_a_col}_pos"].to_numpy(bool)
    b_pos = obs[f"{gene_b_col}_pos"].to_numpy(bool)
    rows = []
    for a_class, a_mask in {f"{a_name}-": ~a_pos, f"{a_name}+": a_pos}.items():
        base = source_mask & a_mask
        total = int(base.sum())
        for b_class, b_mask in {f"{b_name}+": b_pos, f"{b_name}-": ~b_pos}.items():
            n = int((base & b_mask).sum())
            rows.append({"x_class": a_class, "stack_class": b_class, "n": n, "pct": safe_pct(n, total), "total": total})
    return pd.DataFrame(rows)


def macrophage_sirpa_table(adata, macrophage_mask):
    obs = adata.obs
    sirpa_pos = obs["SIRPA_pos"].to_numpy(bool)
    total = int(macrophage_mask.sum())
    rows = []
    for cls, mask in {"SIRPA+": macrophage_mask & sirpa_pos, "SIRPA-": macrophage_mask & ~sirpa_pos}.items():
        n = int(mask.sum())
        rows.append({"bar": "Macrophages", "sirpa_class": cls, "n": n, "pct": safe_pct(n, total), "total": total})
    return pd.DataFrame(rows)


def immune_tigit_status_table(adata, target_mask, target_label):
    obs = adata.obs
    tigit_pos = obs["TIGIT_pos"].to_numpy(bool)
    total = int(target_mask.sum())
    rows = []
    for cls, mask in {"TIGIT+": target_mask & tigit_pos, "TIGIT-": target_mask & ~tigit_pos}.items():
        n = int(mask.sum())
        rows.append({"target_cell_type": target_label, "tigit_class": cls, "n": n, "pct": safe_pct(n, total), "total": total})
    return pd.DataFrame(rows)


def classify_macrophage_cd47_exposure(adata, graph, endometrial_mask, macrophage_mask, local_hop):
    cd47_pos = adata.obs["CD47_pos"].to_numpy(bool)
    mac_idx = np.where(macrophage_mask)[0]
    has_cd47_pos = reachable_hit_within_hop(graph, mac_idx, endometrial_mask & cd47_pos, local_hop)
    exposure = np.full(adata.n_obs, "not_macrophage", dtype=object)
    exposure[macrophage_mask & has_cd47_pos] = "CD47+"
    exposure[macrophage_mask & ~has_cd47_pos] = "CD47-"
    return exposure


def macrophage_signature_table(adata, graph, endometrial_mask, macrophage_mask, local_hop):
    obs = adata.obs
    sirpa_pos = obs["SIRPA_pos"].to_numpy(bool)
    cd47_exposure = classify_macrophage_cd47_exposure(adata, graph, endometrial_mask, macrophage_mask, local_hop)
    rows = []
    for sig in SIGNATURES:
        tmp = pd.DataFrame({"obs_name": obs.index.astype(str), "signature": sig, "signature_score": pd.to_numeric(obs[sig], errors="coerce").to_numpy(float), "sirpa_class": np.where(sirpa_pos, "SIRPA+", "SIRPA-"), "cd47_exposure": cd47_exposure, "is_macrophage": macrophage_mask})
        tmp = tmp[tmp["is_macrophage"]].dropna(subset=["signature_score"]).copy()
        tmp["group"] = tmp["cd47_exposure"] + " / " + tmp["sirpa_class"]
        rows.append(tmp)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def endometrial_epcam_cd47_depth_table(adata, endometrial_mask):
    obs = adata.obs
    total = pd.to_numeric(obs["total_transcripts_pipeline"], errors="coerce").to_numpy(float)
    epcam_pos = obs["EPCAM_pos"].to_numpy(bool)
    cd47_pos = obs["CD47_pos"].to_numpy(bool)
    groups = {"CD47+ / EPCAM+": endometrial_mask & cd47_pos & epcam_pos, "CD47+ / EPCAM-": endometrial_mask & cd47_pos & ~epcam_pos, "CD47- / EPCAM+": endometrial_mask & ~cd47_pos & epcam_pos, "CD47- / EPCAM-": endometrial_mask & ~cd47_pos & ~epcam_pos}
    rows = []
    for group, mask in groups.items():
        rows.append(pd.DataFrame({"obs_name": obs.index.astype(str), "group": group, "total_transcripts": total, "keep": mask}).query("keep"))
    return pd.concat(rows, ignore_index=True).replace([np.inf, -np.inf], np.nan).dropna(subset=["total_transcripts"])


def tigit_nectin2_depth_table(adata, graph, endometrial_mask, target_mask, target_label, local_hop):
    obs = adata.obs
    total = pd.to_numeric(obs["total_transcripts_pipeline"], errors="coerce").to_numpy(float)
    tigit_pos = obs["TIGIT_pos"].to_numpy(bool)
    nectin2_pos = obs["NECTIN2_pos"].to_numpy(bool)
    target_idx = np.where(target_mask)[0]
    has_nectin2_pos = reachable_hit_within_hop(graph, target_idx, endometrial_mask & nectin2_pos, local_hop)
    exposure = np.zeros(adata.n_obs, dtype=bool)
    exposure[target_mask] = has_nectin2_pos[target_mask]
    groups = {"TIGIT+ / NECTIN2+": target_mask & tigit_pos & exposure, "TIGIT+ / NECTIN2-": target_mask & tigit_pos & ~exposure, "TIGIT- / NECTIN2+": target_mask & ~tigit_pos & exposure, "TIGIT- / NECTIN2-": target_mask & ~tigit_pos & ~exposure}
    rows = []
    for group, mask in groups.items():
        rows.append(pd.DataFrame({"target_cell_type": target_label, "group": group, "total_transcripts": total, "keep": mask}).query("keep"))
    return pd.concat(rows, ignore_index=True).dropna(subset=["total_transcripts"])


def plot_celltype_composition(ax, df):
    bottom = 0.0
    for _, r in df.iterrows():
        color = CELLTYPE_COLORS.get(str(r["cell_type_norm"]), "#CCCCCC")
        ax.bar([0], [r["pct"]], bottom=bottom, color=color, edgecolor="black", linewidth=0.5, label=f"{r['cell_type']} (n={int(r['n'])})")
        bottom += r["pct"]
    ax.set_xticks([0]); ax.set_xticklabels(["All cells"]); ax.set_ylim(0, 100); ax.set_ylabel("% cells"); ax.set_title("Cell-type composition")
    ax.legend(frameon=False, fontsize=6, loc="upper right")


def plot_stacked_two_bar(ax, df, x_order, stack_order, origin, ylabel, title):
    pivot = df.pivot_table(index="x_class", columns="stack_class", values="n", fill_value=0).reindex(index=x_order, columns=stack_order, fill_value=0)
    totals = pivot.sum(axis=1).replace(0, np.nan)
    x = np.arange(len(x_order)); bottom = np.zeros(len(x_order))
    for cls in stack_order:
        positive = cls.endswith("+")
        vals = 100.0 * pivot[cls].to_numpy(float) / totals.to_numpy(float)
        vals = np.nan_to_num(vals)
        ax.bar(x, vals, bottom=bottom, color=origin_rgba(origin, positive), edgecolor="black", linewidth=0.6, label=cls)
        bottom += vals
    ax.set_xticks(x); ax.set_xticklabels([f"{g}\nn={int(pivot.loc[g].sum())}" for g in x_order])
    ax.set_ylim(0, 100); ax.set_ylabel(ylabel); ax.set_title(title); ax.legend(frameon=False, fontsize=8)


def plot_single_stacked_bar(ax, df, class_col, order, origin, title):
    total = int(df["n"].sum()); bottom = 0.0
    for cls in order:
        n = int(df.loc[df[class_col] == cls, "n"].sum())
        pct = safe_pct(n, total)
        ax.bar([0], [pct], bottom=bottom, color=origin_rgba(origin, cls.endswith("+")), edgecolor="black", linewidth=0.6, label=f"{cls} (n={n})")
        bottom += pct
    ax.set_xticks([0]); ax.set_xticklabels([f"Total\nn={total}"]); ax.set_ylim(0, 100); ax.set_ylabel("% cells"); ax.set_title(title); ax.legend(frameon=False, fontsize=8)


def plot_signature_boxplot(ax, df):
    group_order = ["CD47+ / SIRPA+", "CD47+ / SIRPA-", "CD47- / SIRPA+", "CD47- / SIRPA-"]
    sig = "macrophage_phagocytosis_negative_regulation"
    sub = df[df["signature"] == sig].copy()
    data = [sub[sub["group"] == g]["signature_score"].dropna().to_numpy(float) for g in group_order]
    plot_data = [x if len(x) else np.array([np.nan]) for x in data]
    positions = np.arange(1, len(group_order) + 1)
    ax.boxplot(plot_data, positions=positions, widths=0.45, showfliers=False, patch_artist=True, boxprops=dict(facecolor=origin_rgba("macrophage", True), alpha=0.35, linewidth=0.8), medianprops=dict(color="black", linewidth=1.2), whiskerprops=dict(color="black", linewidth=0.8), capprops=dict(color="black", linewidth=0.8))
    ref_vals = data[0]
    finite_vals = np.concatenate([x[np.isfinite(x)] for x in data if len(x)])
    ymax = np.nanmax(finite_vals) if finite_vals.size else 1.0
    ymin = np.nanmin(finite_vals) if finite_vals.size else 0.0
    ypad = 0.08 * (ymax - ymin + 1e-6)
    for j in range(1, len(group_order)):
        vals = data[j]
        p = mannwhitneyu(ref_vals, vals, alternative="two-sided").pvalue if len(ref_vals) >= 3 and len(vals) >= 3 else np.nan
        ax.text(positions[j], ymax + ypad, significance_stars(p), ha="center", va="bottom", fontsize=9)
    ax.axhline(0, linestyle="--", linewidth=0.7, alpha=0.45)
    ax.set_xticks(positions); ax.set_xticklabels([f"{g}\nn={len(v)}" for g, v in zip(group_order, data)], rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("Signature score"); ax.set_title("Negative phagocytosis regulation\nReference: CD47+ / SIRPA+")


def plot_depth_trajectory(ax, df, group_order, origin, title):
    medians, errors, ns = [], [], []
    for group in group_order:
        vals = df[df["group"] == group]["total_transcripts"].dropna().to_numpy(float)
        ns.append(len(vals))
        if len(vals) == 0:
            medians.append(np.nan); errors.append(np.nan)
        else:
            medians.append(np.nanmedian(vals)); q1, q3 = np.nanpercentile(vals, [25, 75]); errors.append((q3 - q1) / 2)
    x = np.arange(len(group_order))
    ax.errorbar(x, medians, yerr=errors, fmt="o-", color=ORIGIN_COLORS[origin], linewidth=2, capsize=4, markersize=6)
    ax.set_xticks(x); ax.set_xticklabels([f"{g}\nn={n}" for g, n in zip(group_order, ns)], rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Median total transcripts"); ax.set_title(title)


def plot_signature_violin_separate(df, out_png):
    group_order = ["CD47+ / SIRPA+", "CD47+ / SIRPA-", "CD47- / SIRPA+", "CD47- / SIRPA-"]
    for sig in df["signature"].unique():
        sub = df[df["signature"] == sig].copy()
        data = [sub[sub["group"] == g]["signature_score"].dropna().to_numpy(float) for g in group_order]
        plot_data = [x if len(x) else np.array([np.nan]) for x in data]
        positions = np.arange(1, len(group_order) + 1)
        fig, ax = plt.subplots(figsize=(7, 5))
        parts = ax.violinplot(plot_data, positions=positions, showmeans=False, showmedians=False, showextrema=False)
        for body in parts["bodies"]:
            body.set_alpha(0.35); body.set_edgecolor("black"); body.set_linewidth(0.6)
        ax.boxplot(plot_data, positions=positions, widths=0.22, showfliers=False, patch_artist=True, boxprops=dict(facecolor="white", alpha=0.75, linewidth=0.8), medianprops=dict(color="black", linewidth=1.2))
        ref_vals = data[0]
        finite_vals = np.concatenate([x[np.isfinite(x)] for x in data if len(x)])
        ymax = np.nanmax(finite_vals) if finite_vals.size else 1.0
        for j in range(1, len(group_order)):
            vals = data[j]
            p = mannwhitneyu(ref_vals, vals, alternative="two-sided").pvalue if len(ref_vals) >= 3 and len(vals) >= 3 else np.nan
            ax.text(j + 1, ymax * 1.05, significance_stars(p), ha="center")
        ax.set_xticks(positions); ax.set_xticklabels([f"{g}\nn={len(v)}" for g, v in zip(group_order, data)], rotation=35, ha="right")
        ax.set_ylabel("Signature score"); ax.set_title(sig)
        fig.subplots_adjust(left=0.15, right=0.95, bottom=0.28, top=0.88)
        fig.savefig(out_png.replace(".png", f".{sig}.png"), dpi=300, bbox_inches=None, pad_inches=0.1)
        plt.close(fig)


def plot_cd47_sirpa_high_signature_zoom(adata, graph, endometrial_mask, macrophage_mask, mac_sig_df, out_png, n_regions=5, crop_radius=120):
    xy = get_xy(adata)
    obs = adata.obs
    cd47_pos = obs["CD47_pos"].to_numpy(bool); sirpa_pos = obs["SIRPA_pos"].to_numpy(bool)
    lookup = pd.Series(np.arange(adata.n_obs), index=adata.obs_names.astype(str))
    sig = "macrophage_phagocytosis_negative_regulation"
    df = mac_sig_df[(mac_sig_df["signature"] == sig) & (mac_sig_df["group"] == "CD47+ / SIRPA+")].dropna(subset=["signature_score"]).copy()
    df = df[df["obs_name"].isin(lookup.index)]
    if df.empty:
        print("No CD47+/SIRPA+ high-signature macrophages for zoom plot.", flush=True); return
    df = df.sort_values("signature_score", ascending=False).head(n_regions)
    center_idx = lookup.loc[df["obs_name"]].to_numpy(dtype=int); centers = xy[center_idx]
    fig, axes = plt.subplots(1, len(centers), figsize=(4.2 * len(centers), 4.2), squeeze=False); axes = axes.ravel()
    for ax, center, idx in zip(axes, centers, center_idx):
        in_crop = (xy[:, 0] >= center[0] - crop_radius) & (xy[:, 0] <= center[0] + crop_radius) & (xy[:, 1] >= center[1] - crop_radius) & (xy[:, 1] <= center[1] + crop_radius)
        ax.scatter(xy[in_crop, 0], xy[in_crop, 1], s=4, color="#E6E6E6", alpha=0.35, linewidths=0)
        ax.scatter(xy[in_crop & endometrial_mask & ~cd47_pos, 0], xy[in_crop & endometrial_mask & ~cd47_pos, 1], s=12, color=origin_rgba("endometrial", False), linewidths=0, label="Endometrial CD47-")
        ax.scatter(xy[in_crop & endometrial_mask & cd47_pos, 0], xy[in_crop & endometrial_mask & cd47_pos, 1], s=16, color=origin_rgba("endometrial", True), linewidths=0, label="Endometrial CD47+")
        ax.scatter(xy[in_crop & macrophage_mask & ~sirpa_pos, 0], xy[in_crop & macrophage_mask & ~sirpa_pos, 1], s=18, color=origin_rgba("macrophage", False), linewidths=0, label="Macrophage SIRPA-")
        ax.scatter(xy[in_crop & macrophage_mask & sirpa_pos, 0], xy[in_crop & macrophage_mask & sirpa_pos, 1], s=22, color=origin_rgba("macrophage", True), edgecolor="black", linewidths=0.3, label="Macrophage SIRPA+")
        ax.scatter(xy[idx, 0], xy[idx, 1], s=90, facecolors="none", edgecolors="black", linewidths=1.5)
        score_val = df.loc[df["obs_name"] == adata.obs_names[idx], "signature_score"].iloc[0]
        ax.set_xlim(center[0] - crop_radius, center[0] + crop_radius); ax.set_ylim(center[1] + crop_radius, center[1] - crop_radius)
        ax.set_title(f"CD47+/SIRPA+\nscore={score_val:.2f}", fontsize=9); ax.set_xticks([]); ax.set_yticks([])
    handles, labels = axes[0].get_legend_handles_labels(); fig.legend(handles, labels, frameon=False, fontsize=8, loc="lower center", ncol=4)
    fig.suptitle("High negative-phagocytosis-regulation CD47+/SIRPA+ macrophage niches", fontsize=13)
    fig.subplots_adjust(left=0.03, right=0.98, bottom=0.18, top=0.82, wspace=0.15)
    fig.savefig(out_png, dpi=300, bbox_inches=None, pad_inches=0.1); plt.close(fig)


def run_publication_pipeline(h5ad_paths, outdir, label_col="dominant_celltype", endometrial_label="Endometrial", macrophage_label="Macrophage", layer=None, use_raw=False, transcript_threshold=0.0, spatial_k_neighbors=8, local_hop=3, min_transcripts=None):
    ensure_dir(outdir)
    for h5ad_path in h5ad_paths:
        sample_id = os.path.splitext(os.path.basename(h5ad_path))[0]
        print(f"Processing {sample_id}", flush=True)
        sample_dir = ensure_dir(os.path.join(outdir, sample_id)); files_dir = ensure_dir(os.path.join(sample_dir, "files")); plots_dir = ensure_dir(os.path.join(sample_dir, "plots"))
        adata = sc.read_h5ad(h5ad_path)
        if min_transcripts is not None:
            total_pre = get_total_transcripts(adata, layer=layer, use_raw=use_raw)
            keep = np.isfinite(total_pre) & (total_pre >= float(min_transcripts))
            print(f"{sample_id}: keeping {keep.sum():,}/{adata.n_obs:,} cells with >= {min_transcripts} transcripts", flush=True)
            adata = adata[keep].copy()
        if label_col not in adata.obs.columns:
            raise KeyError(f"{sample_id}: missing adata.obs['{label_col}'].")
        present_df = add_expression_and_signature_calls(adata=adata, layer=layer, use_raw=use_raw, transcript_threshold=transcript_threshold)
        xy = get_xy(adata); graph = build_spatial_graph(xy, n_neighbors=spatial_k_neighbors)
        labels = adata.obs[label_col].astype(str)
        endometrial_mask, macrophage_mask, nk_mask, tcell_mask = build_cell_masks(labels=labels, endometrial_label=endometrial_label, macrophage_label=macrophage_label)
        print(f"{sample_id}: n_endometrial={endometrial_mask.sum()}, n_macrophage={macrophage_mask.sum()}, n_NK={nk_mask.sum()}, n_Tcell={tcell_mask.sum()}", flush=True)
        celltype_df = composition_table(adata, label_col)
        endo_cd47_df = two_by_two_table(adata, endometrial_mask, "EPCAM", "CD47", "EPCAM", "CD47")
        endo_nectin2_df = two_by_two_table(adata, endometrial_mask, "EPCAM", "NECTIN2", "EPCAM", "NECTIN2")
        mac_sirpa_df = macrophage_sirpa_table(adata, macrophage_mask)
        mac_sig_df = macrophage_signature_table(adata, graph, endometrial_mask, macrophage_mask, local_hop)
        endo_depth_df = endometrial_epcam_cd47_depth_table(adata, endometrial_mask)
        nk_tigit_df = immune_tigit_status_table(adata, nk_mask, "NK cells")
        tcell_tigit_df = immune_tigit_status_table(adata, tcell_mask, "T cells")
        nk_depth_df = tigit_nectin2_depth_table(adata, graph, endometrial_mask, nk_mask, "NK cells", local_hop)
        tcell_depth_df = tigit_nectin2_depth_table(adata, graph, endometrial_mask, tcell_mask, "T cells", local_hop)
        tables = {"signature_genes_present": present_df, "celltype_composition": celltype_df, "endometrial_EPCAM_CD47": endo_cd47_df, "endometrial_EPCAM_NECTIN2": endo_nectin2_df, "macrophage_SIRPA": mac_sirpa_df, "macrophage_signature": mac_sig_df, "endometrial_depth_EPCAM_CD47": endo_depth_df, "NK_TIGIT_status": nk_tigit_df, "Tcell_TIGIT_status": tcell_tigit_df, "NK_depth_TIGIT_NECTIN2": nk_depth_df, "Tcell_depth_TIGIT_NECTIN2": tcell_depth_df}
        for name, df in tables.items():
            df.to_csv(os.path.join(files_dir, f"{sample_id}.{name}.csv"), index=False)
        plot_cd47_sirpa_high_signature_zoom(adata, graph, endometrial_mask, macrophage_mask, mac_sig_df, os.path.join(plots_dir, f"{sample_id}.CD47_SIRPA_high_signature_spatial_zooms.png"), n_regions=5, crop_radius=120)
        fig = plt.figure(figsize=(20, 10)); gs = fig.add_gridspec(nrows=2, ncols=4, width_ratios=[0.85, 0.85, 0.85, 2.3], height_ratios=[1, 1])
        ax1 = fig.add_subplot(gs[0, 0]); ax2 = fig.add_subplot(gs[0, 1]); ax3 = fig.add_subplot(gs[0, 2]); ax5 = fig.add_subplot(gs[1, 0:3]); ax4 = fig.add_subplot(gs[:, 3])
        plot_celltype_composition(ax1, celltype_df)
        plot_stacked_two_bar(ax2, endo_cd47_df, ["EPCAM-", "EPCAM+"], ["CD47+", "CD47-"], "endometrial", "% endometrial cells", "Endometrial EPCAM/CD47")
        plot_single_stacked_bar(ax3, mac_sirpa_df, "sirpa_class", ["SIRPA+", "SIRPA-"], "macrophage", "Macrophage SIRPA status")
        plot_depth_trajectory(ax5, endo_depth_df, ["CD47+ / EPCAM+", "CD47+ / EPCAM-", "CD47- / EPCAM+", "CD47- / EPCAM-"], "endometrial", "Endometrial CD47/EPCAM transcript-depth trajectory")
        plot_signature_boxplot(ax4, mac_sig_df)
        fig.suptitle(f"{sample_id}: CD47/SIRPA macrophage axis", fontsize=15); fig.subplots_adjust(left=0.06, right=0.96, bottom=0.12, top=0.88, wspace=0.60, hspace=0.45)
        fig.savefig(os.path.join(plots_dir, f"{sample_id}.CD47_SIRPA_publication_panel.png"), dpi=300, bbox_inches=None, pad_inches=0.1); plt.close(fig)
        fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
        plot_nectin2_tigit_publication_bar(axes[0], endo_nectin2_df, nk_tigit_df, tcell_tigit_df)
        plot_depth_trajectory(axes[1], nk_depth_df, ["TIGIT+ / NECTIN2+", "TIGIT+ / NECTIN2-", "TIGIT- / NECTIN2+", "TIGIT- / NECTIN2-"], "NK", "NK TIGIT/NECTIN2 transcript-depth trajectory")
        plot_depth_trajectory(axes[2], tcell_depth_df, ["TIGIT+ / NECTIN2+", "TIGIT+ / NECTIN2-", "TIGIT- / NECTIN2+", "TIGIT- / NECTIN2-"], "Tcell", "T-cell TIGIT/NECTIN2 transcript-depth trajectory")
        fig.suptitle(f"{sample_id}: NECTIN2/TIGIT lymphoid axis", fontsize=15); fig.subplots_adjust(left=0.06, right=0.96, bottom=0.22, top=0.84, wspace=0.45)
        fig.savefig(os.path.join(plots_dir, f"{sample_id}.NECTIN2_TIGIT_publication_panel.png"), dpi=300, bbox_inches=None, pad_inches=0.1); plt.close(fig)
        plot_signature_violin_separate(mac_sig_df, os.path.join(plots_dir, f"{sample_id}.CD47_SIRPA_signature_violin.png"))
        print(f"Saved outputs in: {sample_dir}", flush=True)
    return outdir

def plot_nectin2_tigit_publication_bar(ax, endo_nectin2_df, nk_tigit_df, tcell_tigit_df):
    """
    Combined NECTIN2/TIGIT stacked barplot.

    Left:
      Endometrial cells split into EPCAM- and EPCAM+ bars,
      stacked by NECTIN2+ / NECTIN2-.

    Right:
      NK cells and T cells as separate bars,
      stacked by TIGIT+ / TIGIT-.
    """

    # -----------------------------
    # Endometrial EPCAM / NECTIN2
    # -----------------------------
    epcam_order = ["EPCAM-", "EPCAM+"]
    nectin_order = ["NECTIN2+", "NECTIN2-"]

    endo_pivot = (
        endo_nectin2_df
        .pivot_table(
            index="epcam_class",
            columns="nectin2_class",
            values="n",
            aggfunc="sum",
            fill_value=0,
        )
        .reindex(
            index=epcam_order,
            columns=nectin_order,
            fill_value=0,
        )
    )

    # -----------------------------
    # Immune TIGIT
    # -----------------------------
    immune_df = pd.concat([nk_tigit_df, tcell_tigit_df], ignore_index=True)

    immune_order = ["NK cells", "T cells"]
    tigit_order = ["TIGIT+", "TIGIT-"]

    immune_pivot = (
        immune_df
        .pivot_table(
            index="target_cell_type",
            columns="tigit_class",
            values="n",
            aggfunc="sum",
            fill_value=0,
        )
        .reindex(
            index=immune_order,
            columns=tigit_order,
            fill_value=0,
        )
    )

    # Layout: two endometrial bars, gap, two immune bars
    x_endo = np.array([0, 1])
    x_immune = np.array([3, 4])

    # -----------------------------
    # Plot endometrial bars
    # -----------------------------
    endo_totals = endo_pivot.sum(axis=1).replace(0, np.nan)
    bottom = np.zeros(len(epcam_order), dtype=float)

    endo_colors = {
        "NECTIN2+": origin_rgba("endometrial", True),
        "NECTIN2-": origin_rgba("endometrial", False),
    }

    for cls in nectin_order:
        vals = 100.0 * endo_pivot[cls].to_numpy(float) / endo_totals.to_numpy(float)
        vals = np.nan_to_num(vals, nan=0.0)

        ax.bar(
            x_endo,
            vals,
            bottom=bottom,
            color=endo_colors[cls],
            edgecolor="black",
            linewidth=0.6,
            label=cls,
        )

        bottom += vals

    # -----------------------------
    # Plot NK/T-cell TIGIT bars
    # -----------------------------
    immune_totals = immune_pivot.sum(axis=1).replace(0, np.nan)
    bottom = np.zeros(len(immune_order), dtype=float)

    for cls in tigit_order:
        vals = 100.0 * immune_pivot[cls].to_numpy(float) / immune_totals.to_numpy(float)
        vals = np.nan_to_num(vals, nan=0.0)

        colors = [
            origin_rgba("NK", cls == "TIGIT+"),
            origin_rgba("Tcell", cls == "TIGIT+"),
        ]

        ax.bar(
            x_immune,
            vals,
            bottom=bottom,
            color=colors,
            edgecolor="black",
            linewidth=0.6,
            label=cls,
        )

        bottom += vals

    # -----------------------------
    # Formatting
    # -----------------------------
    xticks = np.concatenate([x_endo, x_immune])
    xticklabels = (
        [f"{g}\nn={int(endo_pivot.loc[g].sum())}" for g in epcam_order] +
        [f"{g}\nn={int(immune_pivot.loc[g].sum())}" for g in immune_order]
    )

    ax.set_xticks(xticks)
    ax.set_xticklabels(
        xticklabels,
        rotation=20,
        ha="right",
        fontsize=8,
    )

    ax.set_ylim(0, 100)
    ax.set_ylabel("% cells")
    ax.set_title("Endometrial NECTIN2 and lymphoid TIGIT status")

    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))

    ax.legend(
        unique.values(),
        unique.keys(),
        frameon=False,
        fontsize=8,
        loc="upper right",
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", nargs="+", required=True)
    parser.add_argument("--outdir", default="spatial_publication_plots_out")
    parser.add_argument("--label_col", default="dominant_celltype")
    parser.add_argument("--endometrial_label", default="Endometrial")
    parser.add_argument("--macrophage_label", default="Macrophage")
    parser.add_argument("--layer", default=None)
    parser.add_argument("--use_raw", action="store_true")
    parser.add_argument("--transcript_threshold", type=float, default=0.0)
    parser.add_argument("--spatial_k_neighbors", type=int, default=8)
    parser.add_argument("--local_hop", type=int, default=3)
    parser.add_argument("--min_transcripts", type=float, default=None)
    args = parser.parse_args()
    run_publication_pipeline(h5ad_paths=args.h5ad, outdir=args.outdir, label_col=args.label_col, endometrial_label=args.endometrial_label, macrophage_label=args.macrophage_label, layer=args.layer, use_raw=args.use_raw, transcript_threshold=args.transcript_threshold, spatial_k_neighbors=args.spatial_k_neighbors, local_hop=args.local_hop, min_transcripts=args.min_transcripts)


if __name__ == "__main__":
    main()
