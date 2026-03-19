import os
import numpy as np
import pandas as pd

import mne
from nibabel.freesurfer.io import read_annot, read_geometry
from nilearn import plotting
import matplotlib.pyplot as plt


# ==========================
# 基本路径配置
# ==========================
BASE_DIR = r"/data/home/tqi/data1/share/after_freesurfer"

LH_ANNOT = os.path.join(BASE_DIR, r"FILE/DK-318", "lh.DK318.annot")
RH_ANNOT = os.path.join(BASE_DIR, r"FILE/DK-318", "rh.DK318.annot")

DIAG_SIG_FDR_CSV = os.path.join(
    BASE_DIR,
    r"FILE/MIND_ANOVA/dk318_degree_brainmaps(ANOVA_significant)/p_Diagnosis_FDR",
    "significant_regions_p_Diagnosis_FDR.csv"
)


# ==========================
# 工具函数：加载 fsaverage + annot
# ==========================
def load_fsaverage_and_annot():
    # 将 fsaverage 下载到 FILE 目录，避免权限问题
    subjects_dir = os.path.join(BASE_DIR, "FILE")
    fs_dir = mne.datasets.fetch_fsaverage(subjects_dir=subjects_dir, verbose=True)
    print("fsaverage dir:", fs_dir)

    lh_surf_path = os.path.join(fs_dir, "surf", "lh.inflated")
    rh_surf_path = os.path.join(fs_dir, "surf", "rh.inflated")

    lh_coords, lh_faces = read_geometry(lh_surf_path)
    rh_coords, rh_faces = read_geometry(rh_surf_path)

    print("LH surface vertices:", lh_coords.shape[0])
    print("RH surface vertices:", rh_coords.shape[0])

    lh_labels, lh_ctab, lh_names = read_annot(LH_ANNOT)
    rh_labels, rh_ctab, rh_names = read_annot(RH_ANNOT)

    print("LH annot vertices:", len(lh_labels))
    print("RH annot vertices:", len(rh_labels))

    if lh_coords.shape[0] != len(lh_labels) or rh_coords.shape[0] != len(rh_labels):
        raise RuntimeError(
            "表面顶点数与annot不匹配！\n"
            "你的annot是 163842/半球（fsaverage），但你当前使用的surface不是同一版本。\n"
            "请确认你用的是 fsaverage（不是 fsaverage5/6）。"
        )

    lh_mesh = (lh_coords, lh_faces)
    rh_mesh = (rh_coords, rh_faces)

    return (lh_mesh, rh_mesh,
            lh_labels, lh_ctab, lh_names,
            rh_labels, rh_ctab, rh_names)


# ==========================
# 工具函数：名字处理 + 顶点赋值
# ==========================
def strip_hemi_prefix(s: str) -> str:
    # 把 lh_ / rh_ 前缀去掉，便于匹配 annot 名称
    return s.replace("lh_", "").replace("rh_", "")


def make_vertex_data(labels, names, roi_map):
    """
    labels: annot 返回的每个顶点的标签（整数）
    names : annot 里的 ROI 名称（bytes 列表）
    roi_map: dict, key=roi_norm(去掉lh_/rh_)，value=数值（p或FDR）
    """
    name_list = [n.decode("utf-8") for n in names]
    v = np.full(labels.shape[0], np.nan, dtype=float)

    hit_count = 0
    for lab_id, roi in enumerate(name_list):
        if roi in roi_map:
            v[labels == lab_id] = float(roi_map[roi])
            hit_count += 1
    return v, name_list, hit_count


# ==========================
# 核心函数：画一类 degree 脑图
# ==========================
def plot_degree_brainmaps(
    csv_path,
    out_dir,
    cols,
    lh_mesh,
    rh_mesh,
    lh_labels,
    lh_names,
    rh_labels,
    rh_names,
    title_prefix="DK318",
    fixed_vmin=1.3,
    fixed_vmax=3.0
):
    """
    csv_path : 含 318 个 ROI 结果的表（必须有 'feature' + 若干 p列）
    out_dir  : 输出目录
    cols     : 想要画的列名列表
    """

    os.makedirs(out_dir, exist_ok=True)

    # 3) 读取 degree 结果表（318行），建立 ROI->值 的映射
    df = pd.read_csv(csv_path)
    need_cols = {"feature"} | set(cols)
    missing = need_cols - set(df.columns)
    if missing:
        raise ValueError(f"CSV缺少列：{missing}。当前列名：{list(df.columns)}")

    df["roi_norm"] = df["feature"].astype(str).map(strip_hemi_prefix)

    plots = [
        ("lh", "lateral", lh_mesh, "left"),
        ("lh", "medial",  lh_mesh, "left"),
        ("rh", "lateral", rh_mesh, "right"),
        ("rh", "medial",  rh_mesh, "right"),
    ]

    # ===== 循环每个列名，逐列绘图 =====
    for VALUE_COL in cols:
        print("\n============================")
        print("绘制列:", VALUE_COL)

        # 默认配置
        use_symmetric_cbar = False
        cmap_name = 'viridis'
        is_signed = False
        
        # 辅助文件目录（假设与主csv在同一目录）
        sig_file_dir = os.path.dirname(csv_path)

        # 只显示显著脑区：p<0.05 或 FDR<0.05，其余ROI不上色（NaN）
        # 并转换为 -log10(p)
        if VALUE_COL.lower().startswith("p"):
            df_sig = df[df[VALUE_COL] < 0.05].copy()
            base_val = -np.log10(
                np.clip(df_sig[VALUE_COL].to_numpy(), a_min=np.finfo(float).tiny, a_max=None)
            )
            
            # --- 尝试加载方向信息 (Mean Values) ---
            # 1. Diagnosis: TD vs DD (TD > DD 为正/红)
            if "Diagnosis" in VALUE_COL:
                diag_sig_csv = os.path.join(sig_file_dir, "Significant_Diagnosis_DK318_degree_results.csv")
                try:
                    if os.path.exists(diag_sig_csv):
                        df_dir = pd.read_csv(diag_sig_csv)
                        if "feature" in df_dir.columns:
                            sig_features = df_dir["feature"].astype(str).dropna().unique()
                            df_sig = df[df["feature"].astype(str).isin(sig_features)].copy()
                            base_val = -np.log10(
                                np.clip(df_sig[VALUE_COL].to_numpy(), a_min=np.finfo(float).tiny, a_max=None)
                            )
                        if {"mean_TD", "mean_DD"}.issubset(df_dir.columns):
                            df_dir_sub = df_dir[["feature", "mean_TD", "mean_DD"]]
                            df_sig = pd.merge(df_sig, df_dir_sub, on="feature", how="left")
                except Exception as e:
                    print(f"Warning: Failed to load significant Diagnosis file: {e}")

                if "mean_TD" in df_sig.columns and "mean_DD" in df_sig.columns:
                    diff = df_sig["mean_TD"] - df_sig["mean_DD"]
                    signs = np.sign(diff)
                    signs[diff == 0] = 1
                    signs = signs.fillna(1)

                    df_sig[VALUE_COL] = base_val * signs
                    plot_title_prefix = f"{title_prefix} Signed -log10({VALUE_COL})\n(Red: TD>DD, Blue: DD>TD)"
                    use_symmetric_cbar = True
                    cmap_name = 'coolwarm'
                    is_signed = True
                else:
                    df_sig[VALUE_COL] = base_val
                    plot_title_prefix = f"{title_prefix} -log10({VALUE_COL})"

            # 2. AgeGroup: Adult vs Child (Adult > Child 为正/红)
            elif "AgeGroup" in VALUE_COL:
                age_sig_csv = os.path.join(sig_file_dir, "Significant_AgeGroup_DK318_degree_results.csv")
                
                if ("mean_Adult" not in df_sig.columns) and os.path.exists(age_sig_csv):
                    print(f"Loading directional info from: {age_sig_csv}")
                    try:
                        df_dir = pd.read_csv(age_sig_csv)
                        if "mean_Adult" in df_dir.columns and "mean_Child" in df_dir.columns:
                            df_dir_sub = df_dir[["feature", "mean_Adult", "mean_Child"]]
                            df_sig = pd.merge(df_sig, df_dir_sub, on="feature", how="left")
                    except Exception as e:
                        print(f"Warning: Failed to load direction file: {e}")

                if "mean_Adult" in df_sig.columns and "mean_Child" in df_sig.columns:
                    diff = df_sig["mean_Adult"] - df_sig["mean_Child"]
                    signs = np.sign(diff)
                    signs[diff == 0] = 1
                    signs = signs.fillna(1)
                    
                    df_sig[VALUE_COL] = base_val * signs
                    plot_title_prefix = f"{title_prefix} Signed -log10({VALUE_COL})\n(Red: Adult>Child, Blue: Child>Adult)"
                    use_symmetric_cbar = True
                    cmap_name = 'coolwarm'
                    is_signed = True
                else:
                    df_sig[VALUE_COL] = base_val
                    plot_title_prefix = f"{title_prefix} -log10({VALUE_COL})"

            else:
                # 其他情况 (Interaction等)
                df_sig[VALUE_COL] = base_val
                plot_title_prefix = f"{title_prefix} -log10({VALUE_COL})"
                
        else:
            df_sig = df.copy()
            plot_title_prefix = f"{title_prefix} {VALUE_COL}"

        if df_sig.empty:
            print(f"{VALUE_COL}: 没有显著脑区 (阈值<0.05)，跳过绘图。")
            continue

        # 左/右半球 ROI -> 值
        lh_map = dict(
            zip(
                df_sig.loc[df_sig["feature"].str.startswith("lh_"), "roi_norm"],
                df_sig.loc[df_sig["feature"].str.startswith("lh_"), VALUE_COL],
            )
        )
        rh_map = dict(
            zip(
                df_sig.loc[df_sig["feature"].str.startswith("rh_"), "roi_norm"],
                df_sig.loc[df_sig["feature"].str.startswith("rh_"), VALUE_COL],
            )
        )

        lh_vtx, lh_roi_names, lh_hit = make_vertex_data(
            lh_labels, lh_names, lh_map
        )
        rh_vtx, rh_roi_names, rh_hit = make_vertex_data(
            rh_labels, rh_names, rh_map
        )

        lh_match_rate = np.isfinite(lh_vtx).mean() * 100
        rh_match_rate = np.isfinite(rh_vtx).mean() * 100
        print(f"LH vertex filled (%): {lh_match_rate:.2f}")
        print(f"RH vertex filled (%): {rh_match_rate:.2f}")

        # 5) 设定颜色范围
        all_vals = np.concatenate(
            [lh_vtx[np.isfinite(lh_vtx)], rh_vtx[np.isfinite(rh_vtx)]]
        )
        if all_vals.size == 0:
            print(f"{VALUE_COL}: 顶点全是 NaN，跳过绘图。")
            continue

        if is_signed:
            max_abs = np.max(np.abs(all_vals))
            if max_abs < 1e-6: max_abs = 1.0
            vmin = -max_abs
            vmax = max_abs
        else:
            vmin = fixed_vmin
            vmax = fixed_vmax
            
        print("Color range vmin/vmax:", vmin, vmax)

        # 6) 画四张图：LH/RH × lateral/medial
        col_dir = os.path.join(out_dir, VALUE_COL)
        os.makedirs(col_dir, exist_ok=True)

        # --- 新增：保存显著脑区列表到 CSV ---
        sig_table_path = os.path.join(col_dir, f"significant_regions_{VALUE_COL}.csv")
        # 挑选出显著的行，并按显著性排序
        df_export = df_sig[["feature", VALUE_COL]].copy()
        # 如果是 P 值列，增加一列原始 P 值方便查看
        if "plot_title_prefix" in locals() and "-log10" in plot_title_prefix:
            df_export["p_raw"] = 10**(-df_export[VALUE_COL])
            df_export = df_export.rename(columns={VALUE_COL: "minus_log10_p"})
        
        df_export = df_export.sort_values(by=df_export.columns[1], ascending=False)
        df_export.to_csv(sig_table_path, index=False)
        print(f"显著脑区列表已保存至: {sig_table_path}")
        # -----------------------------------

        for hemi_tag, view, mesh, hemi_name in plots:
            out_file = os.path.join(
                col_dir, f"{VALUE_COL}_{hemi_tag}_{view}.png"
            )
            vtx = lh_vtx if hemi_tag == "lh" else rh_vtx

            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(111, projection='3d')
            
            display = plotting.plot_surf_stat_map(
                surf_mesh=mesh,
                stat_map=vtx,
                hemi=hemi_name,
                view=view,
                colorbar=True,
                symmetric_cbar=use_symmetric_cbar,
                vmin=vmin,
                vmax=vmax,
                cmap=cmap_name,
                title=f"{plot_title_prefix} | {hemi_tag.upper()} {view}",
                axes=ax,
                darkness=None
            )

            # 终极兼容方案：使用 add_contours，并处理可能的 AttributeError
            curr_labels = lh_labels if hemi_tag == "lh" else rh_labels
            try:
                # 尝试使用最高效的方法
                display.add_contours(curr_labels, colors='k', linewidths=0.5)
            except AttributeError:
                # 如果上述方法失败，使用手动匹配长度的 plot_surf_contours
                unique_labels = np.unique(curr_labels)
                plotting.plot_surf_contours(
                    surf_mesh=mesh,
                    roi_map=curr_labels,
                    hemi=hemi_name,
                    view=view,
                    axes=ax,
                    levels=unique_labels,
                    colors=['black'] * len(unique_labels),
                    linewidths=0.2
                )

            # 保存图片
            fig.savefig(out_file, bbox_inches="tight", pad_inches=0.2)
            plt.close(fig)
            print("Saved:", out_file)

        print(f"Done. Four brain maps for {VALUE_COL} saved in: {col_dir}")


# ==========================
# 主程序：分别画 ANOVA / Adult / Child
# ==========================
if __name__ == "__main__":
    # 确保当前工作目录是 BASE_DIR，方便相对路径
    os.chdir(BASE_DIR)
    print("当前工作目录：", os.getcwd())

    (lh_mesh, rh_mesh,
     lh_labels, lh_ctab, lh_names,
     rh_labels, rh_ctab, rh_names) = load_fsaverage_and_annot()

    # ---- ANOVA 结果 ----
    anova_csv = os.path.join(
        BASE_DIR, r"FILE/MIND_ANOVA", "ANOVA_DK318_degree_results.csv"
    )
    anova_out = os.path.join(BASE_DIR, "FILE/MIND_ANOVA/dk318_degree_brainmaps(ANOVA_significant)")
    anova_cols = [
        "p_Interaction", "p_Interaction_FDR",
        "p_Diagnosis_FDR","p_AgeGroup_FDR"
    ]
    plot_degree_brainmaps(
        csv_path=anova_csv,
        out_dir=anova_out,
        cols=anova_cols,
        lh_mesh=lh_mesh,
        rh_mesh=rh_mesh,
        lh_labels=lh_labels,
        lh_names=lh_names,
        rh_labels=rh_labels,
        rh_names=rh_names,
        title_prefix="DK318 ANOVA",
    )

    