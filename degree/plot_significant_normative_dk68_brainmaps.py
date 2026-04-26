import os
import numpy as np
import pandas as pd
import mne
from nibabel.freesurfer.io import read_annot, read_geometry
from nilearn import plotting
import matplotlib.pyplot as plt

BASE_DIR = "/data/home/tqi/data1/share/after_freesurfer"
SUMMARY_CSV = os.path.join(
    BASE_DIR,
    "FILE/test_mean_1.5/Normative_model/Results/Individual_area_volume",
    "DGLM_Zscore_AgeGroup_Diagnosis_significant_summary_fdr.csv",
)
RESULT_CSV = os.path.join(
    BASE_DIR,
    "FILE/test_mean_1.5/Normative_model/Results/Individual_area_volume",
    "DGLM_Zscore_AgeGroup_Diagnosis.csv",
)
OUTPUT_DIR = os.path.join(
    BASE_DIR,
    "FILE/test_mean_1.5/Normative_model/Results/Individual_area_volume",
    "DGLM_Zscore_AgeGroup_Diagnosis_significant_brainmaps_fdr",
)


PLOT_CONFIGS = [
    {
        "test": "mu_agegroup",
        "anatomy_group": "Cortical area",
        "title": "mu_agegroup Cortical area (FDR significant)",
        "out_name": "mu_agegroup_cortical_area",
        "signed": False,
    },
    {
        "test": "mu_agegroup",
        "anatomy_group": "Cortical thickness",
        "title": "mu_agegroup Cortical thickness (FDR significant)",
        "out_name": "mu_agegroup_cortical_thickness",
        "signed": False,
    },
    {
        "test": "mu_agegroup",
        "anatomy_group": "Cortical volume",
        "title": "mu_agegroup Cortical volume (FDR significant)",
        "out_name": "mu_agegroup_cortical_volume",
        "signed": False,
    },
    {
        "test": "sigma_agegroup",
        "anatomy_group": "Cortical thickness",
        "title": "sigma_agegroup Cortical thickness (FDR significant)",
        "out_name": "sigma_agegroup_cortical_thickness",
        "signed": False,
    },
    {
        "test": "sigma_agegroup",
        "anatomy_group": "Cortical volume",
        "title": "sigma_agegroup Cortical volume (FDR significant)",
        "out_name": "sigma_agegroup_cortical_volume",
        "signed": False,
    },
    {
        "test": "sigma_diagnosis",
        "anatomy_group": "Cortical thickness",
        "title": "sigma_diagnosis Cortical thickness (FDR significant)",
        "out_name": "sigma_diagnosis_cortical_thickness",
        "signed": False,
    },
    {
        "test": "sigma_diagnosis",
        "anatomy_group": "Cortical volume",
        "title": "sigma_diagnosis Cortical volume (FDR significant)",
        "out_name": "sigma_diagnosis_cortical_volume",
        "signed": False,
    },
    {
        "test": "sigma_interaction",
        "anatomy_group": "Cortical thickness",
        "title": "sigma_interaction Cortical thickness (FDR significant)",
        "out_name": "sigma_interaction_cortical_thickness",
        "signed": False,
    },
]


def load_fsaverage_and_aparc():
    subjects_dir = os.path.join(BASE_DIR, "FILE")
    fs_dir = mne.datasets.fetch_fsaverage(subjects_dir=subjects_dir, verbose=True)
    print("fsaverage dir:", fs_dir)

    lh_surf_path = os.path.join(fs_dir, "surf", "lh.inflated")
    rh_surf_path = os.path.join(fs_dir, "surf", "rh.inflated")
    lh_coords, lh_faces = read_geometry(lh_surf_path)
    rh_coords, rh_faces = read_geometry(rh_surf_path)

    lh_annot = os.path.join(fs_dir, "label", "lh.aparc.annot")
    rh_annot = os.path.join(fs_dir, "label", "rh.aparc.annot")
    lh_labels, _, lh_names = read_annot(lh_annot)
    rh_labels, _, rh_names = read_annot(rh_annot)

    return (
        (lh_coords, lh_faces),
        (rh_coords, rh_faces),
        lh_labels,
        lh_names,
        rh_labels,
        rh_names,
    )


def normalize_roi_name(value: str) -> str:
    value = str(value)
    value = value.replace("lh_", "").replace("rh_", "")
    for suffix in ("_area", "_thickness", "_volume"):
        if value.endswith(suffix):
            value = value[: -len(suffix)]
    return value


def make_vertex_data(labels, names, roi_map):
    name_list = [name.decode("utf-8") for name in names]
    vertex_data = np.full(labels.shape[0], np.nan, dtype=float)
    hit_count = 0
    for label_id, roi_name in enumerate(name_list):
        if roi_name in roi_map:
            vertex_data[labels == label_id] = float(roi_map[roi_name])
            hit_count += 1
    return vertex_data, hit_count


def parse_features(features_text: str) -> list[str]:
    if pd.isna(features_text):
        return []
    return [item.strip() for item in str(features_text).split(";") if item.strip()]

def add_roi_contours(display, mesh, labels, hemi_name, view, axis):
    try:
        display.add_contours(labels, colors="k", linewidths=0.5)
    except AttributeError:
        unique_labels = np.unique(labels)
        plotting.plot_surf_contours(
            surf_mesh=mesh,
            roi_map=labels,
            hemi=hemi_name,
            view=view,
            axes=axis,
            levels=unique_labels,
            colors=["black"] * len(unique_labels),
            linewidths=0.2,
        )


def get_q_column(test: str) -> str:
    q_map = {
        "mu_agegroup": "q_mu_agegroup",
        "mu_diagnosis": "q_mu_diagnosis",
        "mu_interaction": "q_mu_interaction",
        "sigma_agegroup": "q_sigma_agegroup",
        "sigma_diagnosis": "q_sigma_diagnosis",
        "sigma_interaction": "q_sigma_interaction",
    }
    if test not in q_map:
        raise ValueError(f"未定义 {test} 对应的 FDR 列。")
    return q_map[test]


def build_plot_dataframe(
    summary_df: pd.DataFrame,
    result_df: pd.DataFrame,
    test: str,
    anatomy_group: str,
) -> pd.DataFrame:
    subset = summary_df.loc[
        (summary_df["test"].astype(str) == test)
        & (summary_df["anatomy_group"].astype(str) == anatomy_group)
    ].copy()
    if subset.empty:
        raise ValueError(f"未找到 {test} / {anatomy_group} 对应记录。")

    features = parse_features(subset.iloc[0]["features"])
    plot_df = pd.DataFrame({"feature": features})
    plot_df = plot_df.loc[plot_df["feature"].str.startswith(("lh_", "rh_"))].copy()
    if plot_df.empty:
        raise ValueError(f"{test} / {anatomy_group} 没有可绘制的皮层脑区。")

    q_col = get_q_column(test)
    value_df = result_df[["Feature", q_col]].copy()
    value_df = value_df.rename(columns={"Feature": "feature", q_col: "q_value"})
    value_df["feature"] = value_df["feature"].astype(str)
    value_df["q_value"] = pd.to_numeric(value_df["q_value"], errors="coerce")

    plot_df = plot_df.merge(value_df, on="feature", how="left")
    plot_df = plot_df.loc[np.isfinite(plot_df["q_value"]) & (plot_df["q_value"] > 0)].copy()
    if plot_df.empty:
        raise ValueError(f"{test} / {anatomy_group} 没有可用的 FDR p 值。")

    plot_df["roi_norm"] = plot_df["feature"].map(normalize_roi_name)
    plot_df["plot_value"] = -np.log10(
        np.clip(plot_df["q_value"].to_numpy(dtype=float), a_min=np.finfo(float).tiny, a_max=None)
    )
    return plot_df


def plot_single_map(
    plot_df,
    out_dir,
    title,
    signed,
    lh_mesh,
    rh_mesh,
    lh_labels,
    lh_names,
    rh_labels,
    rh_names,
):
    os.makedirs(out_dir, exist_ok=True)

    lh_df = plot_df.loc[plot_df["feature"].str.startswith("lh_")].copy()
    rh_df = plot_df.loc[plot_df["feature"].str.startswith("rh_")].copy()

    lh_map = dict(zip(lh_df["roi_norm"], lh_df["plot_value"]))
    rh_map = dict(zip(rh_df["roi_norm"], rh_df["plot_value"]))

    lh_vtx, lh_hit = make_vertex_data(lh_labels, lh_names, lh_map)
    rh_vtx, rh_hit = make_vertex_data(rh_labels, rh_names, rh_map)

    print(title)
    print(f"LH ROI matched: {lh_hit}")
    print(f"RH ROI matched: {rh_hit}")

    all_values = np.concatenate([
        lh_vtx[np.isfinite(lh_vtx)],
        rh_vtx[np.isfinite(rh_vtx)],
    ])
    if all_values.size == 0:
        raise ValueError(f"{title} 没有成功映射到脑表面的 ROI。")

    vmin = 1.3
    vmax = float(np.nanmax(all_values))
    if not np.isfinite(vmax) or vmax < vmin:
        vmax = vmin
    cmap = "YlOrRd"

    export_path = os.path.join(out_dir, "significant_regions_for_plot.csv")
    plot_df[["feature", "roi_norm", "q_value", "plot_value"]].to_csv(export_path, index=False)
    print("Saved:", export_path)

    plots = [
        ("lh", "lateral", lh_mesh, "left"),
        ("lh", "medial", lh_mesh, "left"),
        ("rh", "lateral", rh_mesh, "right"),
        ("rh", "medial", rh_mesh, "right"),
    ]

    for hemi_tag, view, mesh, hemi_name in plots:
        out_file = os.path.join(out_dir, f"{hemi_tag}_{view}.png")
        vertex_data = lh_vtx if hemi_tag == "lh" else rh_vtx

        fig = plt.figure(figsize=(6, 6))
        axis = fig.add_subplot(111, projection="3d")
        display = plotting.plot_surf_stat_map(
            surf_mesh=mesh,
            stat_map=vertex_data,
            hemi=hemi_name,
            view=view,
            colorbar=True,
            symmetric_cbar=False,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            title=f"{title} | -log10(FDR p) | {hemi_tag.upper()} {view}",
            axes=axis,
            darkness=None,
        )
        add_roi_contours(
            display,
            mesh,
            lh_labels if hemi_tag == "lh" else rh_labels,
            hemi_name,
            view,
            axis,
        )
        fig.savefig(out_file, bbox_inches="tight", pad_inches=0.2)
        plt.close(fig)
        print("Saved:", out_file)


if __name__ == "__main__":
    os.chdir(BASE_DIR)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    summary_df = pd.read_csv(SUMMARY_CSV)
    result_df = pd.read_csv(RESULT_CSV)
    (
        lh_mesh,
        rh_mesh,
        lh_labels,
        lh_names,
        rh_labels,
        rh_names,
    ) = load_fsaverage_and_aparc()

    for config in PLOT_CONFIGS:
        print("\n========================================")
        print("Processing:", config["out_name"])
        print("========================================")
        plot_df = build_plot_dataframe(
            summary_df=summary_df,
            result_df=result_df,
            test=config["test"],
            anatomy_group=config["anatomy_group"],
        )
        plot_single_map(
            plot_df=plot_df,
            out_dir=os.path.join(OUTPUT_DIR, config["out_name"]),
            title=config["title"],
            signed=config["signed"],
            lh_mesh=lh_mesh,
            rh_mesh=rh_mesh,
            lh_labels=lh_labels,
            lh_names=lh_names,
            rh_labels=rh_labels,
            rh_names=rh_names,
        )

    print("\n全部 FDR 显著脑区脑图绘制完成。")
