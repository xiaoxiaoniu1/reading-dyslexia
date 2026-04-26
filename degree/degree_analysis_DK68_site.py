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

# ==========================
# 工具函数：加载 fsaverage + aparc (DK-68)
# ==========================
def load_fsaverage_and_aparc():
    # 获取 fsaverage 数据集路径
    subjects_dir = os.path.join(BASE_DIR, "FILE")
    fs_dir = mne.datasets.fetch_fsaverage(subjects_dir=subjects_dir, verbose=True)
    print("fsaverage dir:", fs_dir)

    # 表面几何文件
    lh_surf_path = os.path.join(fs_dir, "surf", "lh.inflated")
    rh_surf_path = os.path.join(fs_dir, "surf", "rh.inflated")

    lh_coords, lh_faces = read_geometry(lh_surf_path)
    rh_coords, rh_faces = read_geometry(rh_surf_path)

    # 默认 DK-68 (aparc) 标注文件路径
    LH_ANNOT = os.path.join(fs_dir, "label", "lh.aparc.annot")
    RH_ANNOT = os.path.join(fs_dir, "label", "rh.aparc.annot")

    lh_labels, lh_ctab, lh_names = read_annot(LH_ANNOT)
    rh_labels, rh_ctab, rh_names = read_annot(RH_ANNOT)

    print(f"LH: {lh_coords.shape[0]} vertices, {len(lh_names)} ROIs")
    print(f"RH: {rh_coords.shape[0]} vertices, {len(rh_names)} ROIs")

    lh_mesh = (lh_coords, lh_faces)
    rh_mesh = (rh_coords, rh_faces)

    return (lh_mesh, rh_mesh,
            lh_labels, lh_ctab, lh_names,
            rh_labels, rh_ctab, rh_names)

# ==========================
# 工具函数：名字处理 + 顶点赋值
# ==========================
def strip_hemi_prefix(s: str) -> str:
    return str(s).replace("lh_", "").replace("rh_", "")

def make_vertex_data(labels, names, roi_map):
    name_list = [n.decode("utf-8") for n in names]
    v = np.full(labels.shape[0], np.nan, dtype=float)
    hit_count = 0
    for lab_id, roi in enumerate(name_list):
        if roi in roi_map:
            v[labels == lab_id] = float(roi_map[roi])
            hit_count += 1
    return v, name_list, hit_count

# ==========================
# 核心函数：画脑图
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
    title_prefix="DK68",
    fixed_vmin=1.3,
    fixed_vmax=3.0
):
    os.makedirs(out_dir, exist_ok=True)
    df = pd.read_csv(csv_path)
    need_cols = {"feature"} | set(cols)
    missing = need_cols - set(df.columns)
    if missing:
        raise ValueError(f"CSV缺少列：{missing}。当前列名：{list(df.columns)}")
    df["roi_norm"] = df["feature"].map(strip_hemi_prefix)

    plots = [
        ("lh", "lateral", lh_mesh, "left"),
        ("lh", "medial",  lh_mesh, "left"),
        ("rh", "lateral", rh_mesh, "right"),
        ("rh", "medial",  rh_mesh, "right"),
    ]

    for VALUE_COL in cols:
        print(f"\n绘制列: {VALUE_COL}")
        use_symmetric_cbar = False
        cmap_name = 'viridis'
        is_signed = False
        sig_file_dir = os.path.dirname(csv_path)

        if VALUE_COL.lower().startswith("p"):
            df_sig = df[df[VALUE_COL] < 0.05].copy()

            if "Diagnosis" in VALUE_COL:
                diag_sig_csv = os.path.join(sig_file_dir, "Significant_Diagnosis_DK68_degree_results_with_site.csv")
                if ("mean_TD" not in df_sig.columns) and os.path.exists(diag_sig_csv):
                    try:
                        df_dir = pd.read_csv(diag_sig_csv)
                        if "mean_TD" in df_dir.columns and "mean_DD" in df_dir.columns:
                            df_dir_sub = df_dir[["feature", "mean_TD", "mean_DD"]]
                            df_sig = pd.merge(df_sig, df_dir_sub, on="feature", how="left")
                    except Exception:
                        pass

                if "mean_TD" in df_sig.columns and "mean_DD" in df_sig.columns:
                    diff = df_sig["mean_TD"] - df_sig["mean_DD"]
                    signs = np.sign(diff)
                    signs[diff == 0] = 1
                    signs = signs.fillna(1)
                    base_val = -np.log10(df_sig[VALUE_COL].to_numpy())
                    df_sig[VALUE_COL] = base_val * signs.to_numpy()
                    plot_title_prefix = f"{title_prefix} Signed -log10({VALUE_COL})\n(Red: TD>DD, Blue: DD>TD)"
                    use_symmetric_cbar = True
                    cmap_name = 'coolwarm'
                    is_signed = True
                else:
                    df_sig[VALUE_COL] = -np.log10(df_sig[VALUE_COL].to_numpy())
                    plot_title_prefix = f"{title_prefix} -log10({VALUE_COL})"

            elif "AgeGroup" in VALUE_COL:
                age_sig_csv = os.path.join(sig_file_dir, "Significant_AgeGroup_DK68_degree_results_with_site.csv")
                if ("mean_Adult" not in df_sig.columns) and os.path.exists(age_sig_csv):
                    try:
                        df_dir = pd.read_csv(age_sig_csv)
                        if "mean_Adult" in df_dir.columns and "mean_Child" in df_dir.columns:
                            df_dir_sub = df_dir[["feature", "mean_Adult", "mean_Child"]]
                            df_sig = pd.merge(df_sig, df_dir_sub, on="feature", how="left")
                    except Exception:
                        pass

                if "mean_Adult" in df_sig.columns and "mean_Child" in df_sig.columns:
                    diff = df_sig["mean_Adult"] - df_sig["mean_Child"]
                    signs = np.sign(diff)
                    signs[diff == 0] = 1
                    signs = signs.fillna(1)
                    base_val = -np.log10(df_sig[VALUE_COL].to_numpy())
                    df_sig[VALUE_COL] = base_val * signs.to_numpy()
                    plot_title_prefix = f"{title_prefix} Signed -log10({VALUE_COL})\n(Red: Adult>Child, Blue: Child>Adult)"
                    use_symmetric_cbar = True
                    cmap_name = 'coolwarm'
                    is_signed = True
                else:
                    df_sig[VALUE_COL] = -np.log10(df_sig[VALUE_COL].to_numpy())
                    plot_title_prefix = f"{title_prefix} -log10({VALUE_COL})"

            else:
                df_sig[VALUE_COL] = -np.log10(df_sig[VALUE_COL].to_numpy())
                plot_title_prefix = f"{title_prefix} -log10({VALUE_COL})"
        else:
            df_sig = df.copy()
            plot_title_prefix = f"{title_prefix} {VALUE_COL}"

        if df_sig.empty:
            print(f"{VALUE_COL}: 无显著脑区 (p<0.05)，跳过。")
            continue

        lh_map = dict(zip(df_sig.loc[df_sig["feature"].str.startswith("lh_"), "roi_norm"],
                          df_sig.loc[df_sig["feature"].str.startswith("lh_"), VALUE_COL]))
        rh_map = dict(zip(df_sig.loc[df_sig["feature"].str.startswith("rh_"), "roi_norm"],
                          df_sig.loc[df_sig["feature"].str.startswith("rh_"), VALUE_COL]))

        lh_vtx, _, _ = make_vertex_data(lh_labels, lh_names, lh_map)
        rh_vtx, _, _ = make_vertex_data(rh_labels, rh_names, rh_map)

        all_vals = np.concatenate([
            lh_vtx[np.isfinite(lh_vtx)],
            rh_vtx[np.isfinite(rh_vtx)]
        ])
        if all_vals.size == 0:
            print(f"{VALUE_COL}: 顶点全是 NaN，跳过绘图。")
            continue

        if is_signed:
            max_abs = np.max(np.abs(all_vals))
            if max_abs < 1e-6:
                max_abs = 1.0
            vmin = -max_abs
            vmax = max_abs
        else:
            vmin = fixed_vmin
            vmax = fixed_vmax

        col_dir = os.path.join(out_dir, VALUE_COL)
        os.makedirs(col_dir, exist_ok=True)

        for hemi_tag, view, mesh, hemi_name in plots:
            out_file = os.path.join(col_dir, f"{VALUE_COL}_{hemi_tag}_{view}.png")
            vtx = lh_vtx if hemi_tag == "lh" else rh_vtx
            curr_labels = lh_labels if hemi_tag == "lh" else rh_labels

            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(111, projection='3d')
            
            display = plotting.plot_surf_stat_map(
                surf_mesh=mesh, stat_map=vtx, hemi=hemi_name, view=view,
                colorbar=True, symmetric_cbar=use_symmetric_cbar, vmin=vmin, vmax=vmax,
                cmap=cmap_name, title=f"{plot_title_prefix} | {hemi_tag.upper()} {view}",
                axes=ax
            )

            # 添加脑区轮廓线
            try:
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

            fig.savefig(out_file, bbox_inches="tight", pad_inches=0.2)
            plt.close(fig)
            print("已保存:", out_file)

# ==========================
# 主程序
# ==========================
if __name__ == "__main__":
    os.chdir(BASE_DIR)
    
    # 1. 加载 DK-68 标注和 fsaverage 表面
    (lh_mesh, rh_mesh, 
     lh_labels, _, lh_names, 
     rh_labels, _, rh_names) = load_fsaverage_and_aparc()

    # 2. 配置 DK-68 ANOVA 结果路径
    dk68_csv = os.path.join(BASE_DIR, "FILE/test_mean_1.5/MIND_DK68_ANOVA", "ANOVA_DK68_degree_results.csv")
    dk68_out = os.path.join(BASE_DIR, "FILE/test_mean_1.5/MIND_DK68_ANOVA", "dk68_degree_brainmaps(ANOVA_site)")
    
    # 想要绘制的 P 值列名
    dk68_cols = [
        "p_Interaction","p_Interaction_FDR",
        "p_Diagnosis_FDR", "p_AgeGroup_FDR"
    ]

    # 3. 绘图
    plot_degree_brainmaps(
        csv_path=dk68_csv,
        out_dir=dk68_out,
        cols=dk68_cols,
        lh_mesh=lh_mesh,
        rh_mesh=rh_mesh,
        lh_labels=lh_labels,
        lh_names=lh_names,
        rh_labels=rh_labels,
        rh_names=rh_names,
        title_prefix="DK68 ANOVA"
    )
