import os
import numpy as np
import pandas as pd
import mne
from nibabel.freesurfer.io import read_annot, read_geometry
from nilearn import plotting
import matplotlib.pyplot as plt

BASE_DIR = r"/data/home/tqi/data1/share/after_freesurfer"

LH_ANNOT = os.path.join(BASE_DIR, r"FILE/DK-318", "lh.DK318.annot")
RH_ANNOT = os.path.join(BASE_DIR, r"FILE/DK-318", "rh.DK318.annot")


def load_fsaverage_and_annot():
    subjects_dir = os.path.join(BASE_DIR, "FILE")
    fs_dir = mne.datasets.fetch_fsaverage(subjects_dir=subjects_dir, verbose=True)

    lh_surf_path = os.path.join(fs_dir, "surf", "lh.inflated")
    rh_surf_path = os.path.join(fs_dir, "surf", "rh.inflated")

    lh_coords, lh_faces = read_geometry(lh_surf_path)
    rh_coords, rh_faces = read_geometry(rh_surf_path)

    lh_labels, lh_ctab, lh_names = read_annot(LH_ANNOT)
    rh_labels, rh_ctab, rh_names = read_annot(RH_ANNOT)

    lh_mesh = (lh_coords, lh_faces)
    rh_mesh = (rh_coords, rh_faces)

    return (lh_mesh, rh_mesh,
            lh_labels, lh_ctab, lh_names,
            rh_labels, rh_ctab, rh_names)


def strip_hemi_prefix(s: str) -> str:
    return str(s).replace("lh_", "").replace("rh_", "")


def make_vertex_data(labels, names, roi_map):
    name_list = [n.decode("utf-8") for n in names]
    v = np.full(labels.shape[0], np.nan, dtype=float)
    for lab_id, roi in enumerate(name_list):
        if roi in roi_map:
            v[labels == lab_id] = float(roi_map[roi])
    return v


def apply_direction(df_sig, value_col, direction_csv):
    base_val = -np.log10(df_sig[value_col].to_numpy())
    if not direction_csv or not os.path.exists(direction_csv):
        return df_sig, base_val

    df_dir = pd.read_csv(direction_csv)
    df_dir.columns = df_dir.columns.str.strip()
    df_sig.columns = df_sig.columns.str.strip()

    if "feature" not in df_dir.columns:
        return df_sig, base_val

    df_dir["feature"] = df_dir["feature"].astype(str).str.strip()
    df_sig["feature"] = df_sig["feature"].astype(str).str.strip()

    if "direction" in df_dir.columns:
        dir_map = df_dir.set_index("feature")["direction"].astype(str)
        directions = df_sig["feature"].map(dir_map).fillna("")
        signs = np.where(directions.str.contains("DD<TD", na=False), 1,
                np.where(directions.str.contains("TD<DD", na=False) | directions.str.contains("DD<TD", na=False), -1, 1))
        return df_sig, base_val * signs

    if "mean_DD" in df_dir.columns and "mean_TD" in df_dir.columns:
        df_dir = df_dir.set_index("feature")
        diff = df_sig["feature"].map(df_dir["mean_DD"]) - df_sig["feature"].map(df_dir["mean_TD"])
        signs = np.sign(diff).fillna(1).to_numpy()
        return df_sig, base_val * signs

    return df_sig, base_val


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
    fixed_vmax=3.0,
    direction_csv=None,
    sig_region_csv=None
):
    os.makedirs(out_dir, exist_ok=True)

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

    for VALUE_COL in cols:
        use_symmetric_cbar = False
        cmap_name = "viridis"
        is_signed = False

        if VALUE_COL.lower().startswith("p"):
            if VALUE_COL == "p_fdr" and sig_region_csv and os.path.exists(sig_region_csv):
                df_sig_list = pd.read_csv(sig_region_csv)
                if "feature" in df_sig_list.columns:
                    sig_features = df_sig_list["feature"].astype(str).str.strip().tolist()
                    df_sig = df[df["feature"].astype(str).str.strip().isin(sig_features)].copy()
                else:
                    df_sig = df[df[VALUE_COL] < 0.05].copy()
            else:
                df_sig = df[df[VALUE_COL] < 0.05].copy()
            df_sig, signed_val = apply_direction(df_sig, VALUE_COL, direction_csv)
            df_sig[VALUE_COL] = signed_val
            plot_title_prefix = f"{title_prefix} Signed -log10({VALUE_COL})\n(Red: DD<TD, Blue: TD<DD)"
            use_symmetric_cbar = True
            cmap_name = "coolwarm"
            is_signed = True
        else:
            df_sig = df.copy()
            plot_title_prefix = f"{title_prefix} {VALUE_COL}"

        if df_sig.empty:
            print(f"{VALUE_COL}: 没有显著脑区 (阈值<0.05)，跳过绘图。")
            continue

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

        lh_vtx = make_vertex_data(lh_labels, lh_names, lh_map)
        rh_vtx = make_vertex_data(rh_labels, rh_names, rh_map)

        all_vals = np.concatenate(
            [lh_vtx[np.isfinite(lh_vtx)], rh_vtx[np.isfinite(rh_vtx)]]
        )
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

            try:
                display.add_contours(curr_labels, colors='k', linewidths=0.5)
            except AttributeError:
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
            print("Saved:", out_file)


if __name__ == "__main__":
    os.chdir(BASE_DIR)

    (lh_mesh, rh_mesh,
     lh_labels, lh_ctab, lh_names,
     rh_labels, rh_ctab, rh_names) = load_fsaverage_and_annot()

    adult_csv = os.path.join(
        BASE_DIR, r"FILE/MIND_t_test", "Stratified_Adult_degree_results.csv"
    )
    adult_out = os.path.join(BASE_DIR, "FILE/MIND_t_test/dk318_degree_brainmaps(t-test_Adult_significant)")
    adult_cols = ["p_fdr","p_value"]
    adult_dir_csv = os.path.join(
        BASE_DIR, r"FILE/MIND_t_test", "Stratified_Adult_degree_results.csv"
    )
    adult_sig_p_fdr_csv = os.path.join(
        BASE_DIR, r"FILE/MIND_t_test/dk318_degree_brainmaps(t-test_Adult_significant)",
        "p_fdr", "significant_regions_p_fdr.csv"
    )
    plot_degree_brainmaps(
        csv_path=adult_csv,
        out_dir=adult_out,
        cols=adult_cols,
        lh_mesh=lh_mesh,
        rh_mesh=rh_mesh,
        lh_labels=lh_labels,
        lh_names=lh_names,
        rh_labels=rh_labels,
        rh_names=rh_names,
        title_prefix="DK318 t-test Adult",
        direction_csv=adult_dir_csv,
        sig_region_csv=adult_sig_p_fdr_csv
    )

    child_csv = os.path.join(
        BASE_DIR, r"FILE/MIND_t_test", "Stratified_Child_degree_results.csv"
    )
    child_out = os.path.join(BASE_DIR, "FILE/MIND_t_test/dk318_degree_brainmaps(t-test_Child_significant)")
    child_cols = ["p_fdr"]
    child_dir_csv = os.path.join(
        BASE_DIR, r"FILE/MIND_t_test", "Significant_Child_degree_results.csv"
    )
    plot_degree_brainmaps(
        csv_path=child_csv,
        out_dir=child_out,
        cols=child_cols,
        lh_mesh=lh_mesh,
        rh_mesh=rh_mesh,
        lh_labels=lh_labels,
        lh_names=lh_names,
        rh_labels=rh_labels,
        rh_names=rh_names,
        title_prefix="DK318 t-test Child",
        direction_csv=child_dir_csv
    )