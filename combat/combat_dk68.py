#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
使用【Jfortin1 官方 NeuroComBat】对 MIND edge（DK-68 上三角 2278 维）做批次校正
- 基于: https://github.com/Jfortin1/neuroCombat
- batch: site（作为 covars 中的一列）
- covariates: site, group_d_or_c, age_month, sex
- 自动过滤仅含单一 group 的站点
- 启用 mean_only=True 确保稳定性
- 输入为完整的 68×68 CSV 矩阵，提取上三角进行校正
"""

import os
import numpy as np
import pandas as pd
from neuroCombat import neuroCombat

# ====================== 路径配置 ======================
EXCEL_PATH = "/data/home/tqi/data1/share/after_freesurfer/FILE/all_data_cqt.xlsx"
MIND_DIR = "/data/home/tqi/data1/share/after_freesurfer/FILE/MIND_DK68"
OUT_DIR = "/data/home/tqi/data1/share/after_freesurfer/FILE/MIND_DK68_combat"
os.makedirs(OUT_DIR, exist_ok=True)

N_ROI = 68
TRIU_IDX = np.triu_indices(N_ROI, k=1)  # 上三角索引，不含对角线


def calc_degree(mat: np.ndarray) -> np.ndarray:
    mat = mat.copy()
    np.fill_diagonal(mat, np.nan)
    return np.nanmean(mat, axis=1)

# ====================== 1. 读协变量表 ======================
cov_df = pd.read_excel(EXCEL_PATH)
print("协变量表总行数：", len(cov_df))

mind_rows = []
ages = []
sexs = []
groups = []
sites = []
used_ids = []

# ====================== 2. 遍历被试，加载 MIND 矩阵 ======================
for _, row in cov_df.iterrows():
    sid = str(int(row["id"]))  # 使用 id 列对齐
    
    csv_path = os.path.join(MIND_DIR, f"{sid}.csv")
    if not os.path.exists(csv_path):
        print(f"[跳过] 找不到 MIND：{csv_path}")
        continue
    
    # 读取 CSV 矩阵（第一行和第一列是标签）
    mat_df = pd.read_csv(csv_path, index_col=0)
    mat = mat_df.values
    
    if mat.shape != (N_ROI, N_ROI):
        print(f"[跳过] ID={sid} 的矩阵形状错误：{mat.shape}")
        continue
    
    # 提取上三角（不含对角线）
    edges = mat[TRIU_IDX]  # 68×67/2 = 2278 维
    mind_rows.append(edges)
    ages.append(row["age_month"])
    sexs.append(row["sex"])
    groups.append(row["group_d_or_c"])
    sites.append(row["site"])
    used_ids.append(sid)

mind_rows = np.array(mind_rows)  # (n_subj, 2278)
n_subj, n_edge = mind_rows.shape
print(f"\n有效被试数: {n_subj}, 每个被试 edge 数: {n_edge}")

if n_subj == 0:
    raise ValueError("没有有效被试！")

# ====================== 3. 构建 covars（官方要求：包含 site） ======================
dat = mind_rows.T  # (2278, n_subj)

covars_df = pd.DataFrame({
    "site": sites,
    "age": ages,
    "sex": sexs,
    "group": groups
})

# 强制分类变量为字符串（避免数值型被误判）
covars_df["site"] = covars_df["site"].astype(str)
covars_df["sex"] = covars_df["sex"].astype(str)
covars_df["group"] = covars_df["group"].astype(str)
covars_df["age"] = covars_df["age"].astype(float)


# ====================== 4. 过滤仅含单一 group 的站点 ======================
print("\n【诊断】site × group 交叉表:")
cross_tab = pd.crosstab(covars_df['site'], covars_df['group'])
print(cross_tab)

# 仅保留同时包含 group=0 和 group=1 的站点
valid_sites = cross_tab[(cross_tab.iloc[:, 0] > 0) & (cross_tab.iloc[:, 1] > 0)].index
print(f"\n✅ 有效站点: {sorted(valid_sites.tolist())}")

mask = covars_df['site'].isin(valid_sites)
dat = dat[:, mask]
covars_df = covars_df[mask].reset_index(drop=True)
used_ids = [sid for i, sid in enumerate(used_ids) if mask.iloc[i]]

print(f"\n过滤后被试数: {len(used_ids)}")

# ====================== 5. 调用【Jfortin1 官方 NeuroComBat】 ======================
print("\n[ComBat] 使用官方 NeuroComBat (mean_only=True) 开始校正...")
res = neuroCombat(
    dat=dat,
    covars=covars_df,           # 包含 "site" 列
    batch_col="site",           # 指定哪一列是 batch
    categorical_cols=["sex", "group"],
    continuous_cols=["age"]
)

dat_combat = res["data"]  # (2278, n_subj)
print("[ComBat] 校正完成！")

# ====================== 6. 重建并保存 68×68 矩阵 + degree ======================
for i, sid in enumerate(used_ids):
    edges_h = dat_combat[:, i]
    
    mat_h = np.zeros((N_ROI, N_ROI), dtype=np.float32)
    mat_h[TRIU_IDX] = edges_h
    mat_h.T[TRIU_IDX] = edges_h  # 对称填充
    
    # 保留原始对角线
    orig_df = pd.read_csv(os.path.join(MIND_DIR, f"{sid}.csv"), index_col=0)
    orig_mat = orig_df.values
    np.fill_diagonal(mat_h, np.diagonal(orig_mat))
    
    # 保存 combat 后矩阵
    mat_h_df = pd.DataFrame(mat_h, index=orig_df.index, columns=orig_df.columns)
    matrix_out_path = os.path.join(OUT_DIR, f"{sid}_combat.csv")
    mat_h_df.to_csv(matrix_out_path)
    
    # 保存该被试的 degree
    degree = calc_degree(mat_h)
    degree_df = pd.DataFrame({
        "ROI": orig_df.index,
        "degree": degree
    })
    degree_out_path = os.path.join(OUT_DIR, f"{sid}_combat_degree.csv")
    degree_df.to_csv(degree_out_path, index=False)
    
    print(f"[保存] ID={sid} | matrix={os.path.basename(matrix_out_path)} | degree={os.path.basename(degree_out_path)}")

print(f"\n✅ 全部完成！输出目录：{OUT_DIR}")
