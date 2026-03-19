#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
批量计算 MIND 矩阵和 degree（DK-318）
- 使用原始 MIND 仓库中的 compute_MIND 函数
- 对每个被试输出：
    * <sub>_MIND_DK318.npy：318x318 MIND 矩阵（numpy 数组）
    * <sub>_MIND_DK318.csv：318x318 MIND 矩阵（带行列标签的 csv）
    * <sub>_MIND_DK318_edges.npy：上三角（去对角线）101124 维 edge 向量
    * <sub>_MIND_DK318_degree.npy：318 维 degree
    * <sub>_MIND_DK318_degree.csv：带 ROI 名称的 degree 表
"""

import os
import sys
import numpy as np
import pandas as pd

# =================== 路径设置 ===================

# 1）MIND 代码所在目录（你 git clone 的位置）
MIND_CODE_DIR = "/data1/tqi/share/after_freesurfer/CODE/MIND"  # 如有变化改这里

# 2）FreeSurfer 输出的根目录（包含所有被试）
BASE_FS_DIR = "/data1/tqi/share/after_freesurfer/fs_subjects_all"


# 3）MIND 结果输出目录
MIND_OUT_DIR = os.path.join(BASE_FS_DIR, "MIND_out")
os.makedirs(MIND_OUT_DIR, exist_ok=True)

# =================== 加载 MIND 函数 ===================

sys.path.insert(1, MIND_CODE_DIR)
from MIND import compute_MIND  # 来自原作者仓库

# 使用的特征（MIND 原文默认组合）
FEATURES = ["CT", "MC", "Vol", "SD", "SA"]

# 使用的脑区分区：必须和 label/ 里的 .annot 文件后缀一致
# 若你生成的是 lh.DK318.annot / rh.DK318.annot，就用 "DK318"
# 若是 lh.DK-318.annot / rh.DK-318.annot，就改成 "DK-318"
PARCELLATION = "DK318"

# =================== 工具函数 ===================

def find_subjects(base_dir: str):
    subs = []
    for name in os.listdir(base_dir):
        full_path = os.path.join(base_dir, name)
        surf_dir = os.path.join(full_path, "surf")

        # 必须是目录 + 有 surf 子目录
        if os.path.isdir(full_path) and os.path.isdir(surf_dir):
            subs.append(name)

    subs.sort()
    print(f"找到的被试数量: {len(subs)}")
    print(f"被试列表: {subs}")
    return subs


def compute_degree_from_matrix(mat: np.ndarray) -> np.ndarray:
    """
    按你文档要求：
    对每一行，把 318*318 矩阵中该行除对角线以外的元素求和，再除以 (N-1)
    得到每个脑区与其他脑区的平均相似度（degree）
    """
    n = mat.shape[0]
    # 行和减去对角线
    row_sum = mat.sum(axis=1) - np.diag(mat)
    degree = row_sum / (n - 1)
    return degree


def main():
    subjects = find_subjects(BASE_FS_DIR)
    print(f"在 {BASE_FS_DIR} 下找到 {len(subjects)} 个被试：{subjects}")

    # 预先算好上三角索引（k=1 去掉对角线）
    n_roi = 318  # DK-318
    triu_idx = np.triu_indices(n_roi, k=1)

    for sub in subjects:
        surf_dir = os.path.join(BASE_FS_DIR, sub)
        print(f"\n====== 处理被试 {sub} ======")
        print(f"FreeSurfer 目录：{surf_dir}")
        
        #====断点续算：如果文件已经存在，就跳过这个被试====
        npy_path=os.path.join(MIND_OUT_DIR,f"{sub}_MIND_DK318.npy")
        if os.path.exists(npy_path):
          print(f"发现该被试已经有结果，跳过：{npy_path}")
          continue
          
        try:
        # ---- 1) 调用原始 compute_MIND 计算 318x318 网络 ----
          mind_df = compute_MIND(
            surf_dir,
            FEATURES,
            PARCELLATION
          )  # 返回 DataFrame: ROI x ROI
        except Exception as e:
          print(f"❌ 跳过被试 {sub}，原因：{e}")
          # 记录错误到日志文件，便于排查
          with open("bad_subjects.log", "a") as f:
            f.write(f"{sub}\t{repr(e)}\n")
          continue  # 跳过当前被试，继续下一个

        # DataFrame -> numpy 矩阵
        mat = mind_df.to_numpy(dtype=np.float32)  # 318 x 318

        # ---- 2) 保存完整的 MIND 矩阵 ----
        npy_path = os.path.join(MIND_OUT_DIR, f"{sub}_MIND_DK318.npy")
        csv_path = os.path.join(MIND_OUT_DIR, f"{sub}_MIND_DK318.csv")
        np.save(npy_path, mat)
        mind_df.to_csv(csv_path)
        print(f"已保存 MIND 矩阵到：\n  {npy_path}\n  {csv_path}")

        # ---- 3) 展开上三角得到 edge 向量（50403 维）----
        edges = mat[triu_idx]
        edges_path = os.path.join(MIND_OUT_DIR, f"{sub}_MIND_DK318_edges.npy")
        np.save(edges_path, edges.astype(np.float32))
        print(f"已保存 edge 向量到：\n  {edges_path}")

        # ---- 4) 计算 degree（按你说的“每一行除对角线加和再除以 317”）----
        degree = compute_degree_from_matrix(mat)

        degree_npy_path = os.path.join(MIND_OUT_DIR, f"{sub}_MIND_DK318_degree.npy")
        degree_csv_path = os.path.join(MIND_OUT_DIR, f"{sub}_MIND_DK318_degree.csv")

        np.save(degree_npy_path, degree.astype(np.float32))

        degree_df = pd.DataFrame({
            "region": mind_df.index,   # ROI 名称
            "degree": degree
        })
        degree_df.to_csv(degree_csv_path, index=False)

        print(f"已保存 degree 到：\n  {degree_npy_path}\n  {degree_csv_path}")

    print("\n全部被试处理完成！")


if __name__ == "__main__":
    main()
