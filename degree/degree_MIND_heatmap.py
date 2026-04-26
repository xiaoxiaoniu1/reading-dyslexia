#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
按 group_age 和 group_d_or_c 将被试分成四组（Adult/Child × DD/TD），
分别计算 MIND 平均矩阵，并单独输出四张热力图。
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

# =========================
# 1) 基本配置
# =========================
demo_file = "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/all_data_cqt_mean_1.5.xlsx"

# 选择矩阵类型： "DK68" 或 "DK318"
MIND_TYPE = "DK318"

# DK68 配置
DK68_MIND_DIR = "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK68_combat"

# DK318 配置（按已有脚本习惯）
DK318_MIND_DIR = "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_combat/"

out_dir = "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_group_mean_heatmaps"
os.makedirs(out_dir, exist_ok=True)

# 是否显示全部 ROI 标签
SHOW_ALL_LABELS = True

# 单张图尺寸（英寸）
FIGSIZE = (18, 18)

# DK318 仅显示 part1 标签
DK318_LABEL_SUFFIX = "_part1"

# 如果想固定色条范围，可手动设置；否则自动从四组平均矩阵中取 min/max
FIXED_VMIN = None
FIXED_VMAX = None

# =========================
# 2) 读取人口学信息并构建分组
# =========================
df = pd.read_excel(demo_file, sheet_name="Sheet1")

df["Diagnosis"] = df["group_d_or_c"].map({0: "TD", 1: "DD"})
df["AgeGroup"] = df["group_age"].map({1: "Adult", 2: "Child"})
df["Group"] = df["AgeGroup"] + "-" + df["Diagnosis"]

# 构建文件路径
if MIND_TYPE == "DK68":
    df["id"] = df["id"].astype(str)
    df["mind_file"] = df["id"].apply(lambda x: os.path.join(DK68_MIND_DIR, f"{x}_combat.csv"))
elif MIND_TYPE == "DK318":
    df["original_project"] = df["original-project"].astype(str)
    df["id_old"] = df["id_old"].astype(str)
    df["mind_file"] = df.apply(lambda row: os.path.join(DK318_MIND_DIR, f"{row['original_project']}_{row['id_old']}_MIND_DK318_combat.csv"), axis=1)
else:
    raise ValueError("MIND_TYPE 只能是 'DK68' 或 'DK318'")

df["has_file"] = df["mind_file"].apply(os.path.exists)
df = df.dropna(subset=["Group"])
df = df[df["has_file"]]

print(f"Total subjects with matrix: {len(df)}")
print(df["Group"].value_counts())

# =========================
# 3) 读取矩阵并累计求均值
# =========================
def read_mind_matrix(fp):
    mat = pd.read_csv(fp, index_col=0)
    mat = mat.apply(pd.to_numeric, errors="coerce")
    return mat

# 用第一份矩阵确定 ROI 顺序
first_file = df["mind_file"].iloc[0]
first_mat = read_mind_matrix(first_file)
roi_names = first_mat.columns.tolist()
n_roi = len(roi_names)
upper_tri_idx = np.triu_indices(n_roi, k=1)

if MIND_TYPE == "DK318":
    label_names = [name if name.endswith(DK318_LABEL_SUFFIX) else "" for name in roi_names]
else:
    label_names = roi_names

group_order = ["Adult-DD", "Adult-TD", "Child-DD", "Child-TD"]
sum_mats = {g: np.zeros((n_roi, n_roi), dtype=float) for g in group_order}
counts = {g: 0 for g in group_order}
group_edges = {g: [] for g in group_order}

for _, row in df.iterrows():
    group = row["Group"]
    if group not in sum_mats:
        continue
    try:
        mat = read_mind_matrix(row["mind_file"]).values
        if mat.shape != (n_roi, n_roi):
            print(f"Skip (shape mismatch): {row['mind_file']}")
            continue
        sum_mats[group] += mat
        counts[group] += 1
        group_edges[group].append(mat[upper_tri_idx])
    except Exception as e:
        print(f"Skip (read error): {row['mind_file']} | {e}")

mean_mats = {}
for g in group_order:
    if counts[g] > 0:
        mean_mats[g] = sum_mats[g] / counts[g]
    else:
        mean_mats[g] = None
        print(f"Warning: {g} has 0 subjects")

# 自动设置色条范围
if FIXED_VMIN is None or FIXED_VMAX is None:
    all_vals = []
    for g in group_order:
        if mean_mats[g] is not None:
            all_vals.append(mean_mats[g].ravel())
    all_vals = np.concatenate(all_vals)
    vmin = np.nanmin(all_vals) if FIXED_VMIN is None else FIXED_VMIN
    vmax = np.nanmax(all_vals) if FIXED_VMAX is None else FIXED_VMAX
else:
    vmin, vmax = FIXED_VMIN, FIXED_VMAX

# =========================
# 4) 分别绘制四张热力图
# =========================
sns.set(style="white")

for group in group_order:
    mat = mean_mats[group]
    if mat is None:
        print(f"Skip plot (no data): {group}")
        continue

    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=300)
    mask = np.eye(n_roi, dtype=bool)

    hm = sns.heatmap(
        mat,
        ax=ax,
        cmap="Reds",
        vmin=vmin,
        vmax=vmax,
        square=True,
        cbar=False,
        mask=mask,
        xticklabels=label_names if SHOW_ALL_LABELS else False,
        yticklabels=label_names if SHOW_ALL_LABELS else False
    )

    fig.subplots_adjust(right=0.80)
    pos = ax.get_position()
    cax = fig.add_axes([0.83, pos.y0, 0.02, pos.height])
    fig.colorbar(hm.collections[0], cax=cax)

    if SHOW_ALL_LABELS:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=5)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=5)

    ax.set_title(group, fontsize=18, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")

    plt.tight_layout(rect=[0, 0, 0.80, 1])
    out_file = os.path.join(out_dir, f"MIND_mean_heatmap_{MIND_TYPE}_{group}.png")
    plt.savefig(out_file, bbox_inches="tight")
    plt.close()

    print(f"Saved: {out_file}")

# =========================
# 5) 组间差异矩阵（DD - TD）与 FDR 显著差异矩阵
# =========================

def fdr_bh(pvals, alpha=0.05):
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    reject = np.zeros_like(p, dtype=bool)

    valid = np.isfinite(p)
    if valid.sum() == 0:
        return reject, q

    pv = p[valid]
    m = pv.size
    order = np.argsort(pv)
    pv_sorted = pv[order]

    q_sorted = pv_sorted * m / (np.arange(1, m + 1))
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0, 1)

    rej_sorted = q_sorted <= alpha

    q_valid = np.empty_like(pv)
    rej_valid = np.empty_like(pv, dtype=bool)
    q_valid[order] = q_sorted
    rej_valid[order] = rej_sorted

    q[valid] = q_valid
    reject[valid] = rej_valid

    return reject, q


def save_matrix_csv(mat, names, out_csv):
    pd.DataFrame(mat, index=names, columns=names).to_csv(out_csv)


def plot_diff_heatmap(diff_mat, title, out_png, labels, mask=None):
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=300)
    diag_mask = np.eye(diff_mat.shape[0], dtype=bool)
    use_mask = diag_mask if mask is None else (diag_mask | mask)

    vals = diff_mat[~use_mask & np.isfinite(diff_mat)]
    max_abs = np.nanmax(np.abs(vals)) if vals.size > 0 else 1.0
    if not np.isfinite(max_abs) or max_abs < 1e-12:
        max_abs = 1.0

    hm = sns.heatmap(
        diff_mat,
        ax=ax,
        cmap="coolwarm",
        center=0,
        vmin=-max_abs,
        vmax=max_abs,
        square=True,
        cbar=False,
        mask=use_mask,
        xticklabels=labels if SHOW_ALL_LABELS else False,
        yticklabels=labels if SHOW_ALL_LABELS else False,
    )

    fig.subplots_adjust(right=0.80)
    pos = ax.get_position()
    cax = fig.add_axes([0.83, pos.y0, 0.02, pos.height])
    fig.colorbar(hm.collections[0], cax=cax)

    if SHOW_ALL_LABELS:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=5)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=5)

    ax.set_title(title, fontsize=18, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")

    plt.tight_layout(rect=[0, 0, 0.80, 1])
    plt.savefig(out_png, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_png}")


comparisons = [
    ("Adult", "Adult-DD", "Adult-TD"),
    ("Child", "Child-DD", "Child-TD"),
]

for tag, g_dd, g_td in comparisons:
    if mean_mats.get(g_dd) is None or mean_mats.get(g_td) is None:
        print(f"Skip diff ({tag}): missing group mean matrix")
        continue

    diff_mat = mean_mats[g_dd] - mean_mats[g_td]

    out_csv = os.path.join(out_dir, f"MIND_diff_{MIND_TYPE}_{tag}_DD_minus_TD.csv")
    save_matrix_csv(diff_mat, roi_names, out_csv)

    out_png = os.path.join(out_dir, f"MIND_diff_heatmap_{MIND_TYPE}_{tag}_DD_minus_TD.png")
    plot_diff_heatmap(diff_mat, f"{tag} DD - TD", out_png, label_names)

    dd_edges = np.asarray(group_edges.get(g_dd, []), dtype=float)
    td_edges = np.asarray(group_edges.get(g_td, []), dtype=float)

    if dd_edges.size == 0 or td_edges.size == 0:
        print(f"Skip FDR ({tag}): missing subject edges")
        continue

    _, pvals = stats.ttest_ind(dd_edges, td_edges, axis=0, equal_var=False, nan_policy='omit')
    reject, qvals = fdr_bh(pvals, alpha=0.05)

    p_mat = np.full((n_roi, n_roi), np.nan, dtype=float)
    q_mat = np.full((n_roi, n_roi), np.nan, dtype=float)
    sig_mat = np.zeros((n_roi, n_roi), dtype=int)

    p_mat[upper_tri_idx] = pvals
    q_mat[upper_tri_idx] = qvals
    sig_mat[upper_tri_idx] = reject.astype(int)

    p_mat = p_mat + p_mat.T
    q_mat = q_mat + q_mat.T
    sig_mat = sig_mat + sig_mat.T

    out_p_csv = os.path.join(out_dir, f"MIND_diff_pvalues_{MIND_TYPE}_{tag}_DD_vs_TD.csv")
    out_q_csv = os.path.join(out_dir, f"MIND_diff_qvalues_FDR_{MIND_TYPE}_{tag}_DD_vs_TD.csv")
    out_sig_csv = os.path.join(out_dir, f"MIND_diff_significant_FDR_{MIND_TYPE}_{tag}_DD_vs_TD.csv")

    save_matrix_csv(p_mat, roi_names, out_p_csv)
    save_matrix_csv(q_mat, roi_names, out_q_csv)
    save_matrix_csv(sig_mat, roi_names, out_sig_csv)

    nonsig_mask = (sig_mat == 0)
    out_png_sig = os.path.join(out_dir, f"MIND_diff_heatmap_{MIND_TYPE}_{tag}_DD_minus_TD_FDRsig.png")
    plot_diff_heatmap(diff_mat, f"{tag} DD - TD (FDR<0.05)", out_png_sig, label_names, mask=nonsig_mask)