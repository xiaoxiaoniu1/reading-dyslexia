#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
从现成（ComBat后）的 318x318 MIND 矩阵（.npy）
1) 计算 degree
2) 导出带 ROI 名称的 degree.csv
3) 同时导出带 ROI 名称的 318x318 matrix.csv
"""

import os
import numpy as np
import pandas as pd


# =================== 路径设置 ===================

IN_DIR = "/data1/tqi/share/after_freesurfer/fs_subjects_all/MIND_out_combat_group"
OUT_DIR = "/data1/tqi/share/after_freesurfer/fs_subjects_all/MIND_out_combat_degree"
os.makedirs(OUT_DIR, exist_ok=True)

ROI_NAME_FILE = "/data1/tqi/share/after_freesurfer/FILE/DK318_roi_names.csv"


NPY_SUFFIX = "_MIND_DK318_combat.npy"

# ★ 一定要 True
SAVE_LABELED_MATRIX_CSV = True


# =================== 工具函数 ===================

def load_roi_names(path: str, n_roi: int = 318):
    df = pd.read_csv(path)
    names = df["region"].astype(str).tolist()
    if len(names) != n_roi:
        raise ValueError(f"ROI 名称数量不等于 {n_roi}")
    return names


def read_matrix(fp: str) -> np.ndarray:
    mat = np.load(fp)
    mat = np.asarray(mat, dtype=np.float32)
    if mat.shape != (318, 318):
        raise ValueError(f"{fp} 不是 318x318")
    return mat


def compute_degree_from_matrix(mat: np.ndarray) -> np.ndarray:
    n = mat.shape[0]
    row_sum = mat.sum(axis=1) - np.diag(mat)
    return (row_sum / (n - 1)).astype(np.float32)


def find_subject_files(in_dir: str):
    files = [f for f in os.listdir(in_dir) if f.endswith(NPY_SUFFIX)]
    files.sort()
    return files


# =================== 主流程 ===================

def main():
    roi_names = load_roi_names(ROI_NAME_FILE)

    files = find_subject_files(IN_DIR)
    print(f"在 {IN_DIR} 找到 {len(files)} 个 npy 文件")

    for fname in files:
        sub = fname.replace(NPY_SUFFIX, "")
        fp = os.path.join(IN_DIR, fname)

        try:
            mat = read_matrix(fp)

            # ---------- degree ----------
            degree = compute_degree_from_matrix(mat)
            degree_df = pd.DataFrame({
                "region": roi_names,
                "degree": degree
            })
            degree_df.to_csv(
                os.path.join(OUT_DIR, f"{sub}_degree.csv"),
                index=False
            )

            # ---------- ★ 318×318 matrix csv ----------
            if SAVE_LABELED_MATRIX_CSV:
                mat_df = pd.DataFrame(
                    mat,
                    index=roi_names,
                    columns=roi_names
                )
                mat_df.to_csv(
                    os.path.join(OUT_DIR, f"{sub}_MIND_DK318_combat_labeled.csv")
                )

            print(f"✅ {sub} 完成")

        except Exception as e:
            print(f"❌ {sub} 失败：{e}")

    print("\n全部被试处理完成")


if __name__ == "__main__":
    main()
