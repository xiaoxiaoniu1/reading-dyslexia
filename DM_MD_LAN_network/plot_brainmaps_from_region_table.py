import os
import argparse
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import mne
from nibabel.freesurfer.io import read_annot, read_geometry, read_morph_data
from nilearn import plotting
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


BASE_DIR = "/data/home/tqi/data1/share/after_freesurfer"
DEFAULT_CSV_DIR = os.path.join(
    BASE_DIR,
    "FILE/test_mean_1.5/DM_MD_LAN_network",
)
DEFAULT_OUTPUT_DIR = os.path.join(
    DEFAULT_CSV_DIR,
    "brain_region_table_brainmaps",
)
LH_ANNOT = os.path.join(BASE_DIR, "FILE/DK-318", "lh.DK318.annot")
RH_ANNOT = os.path.join(BASE_DIR, "FILE/DK-318", "rh.DK318.annot")
DEFAULT_CSVS = [
    os.path.join(DEFAULT_CSV_DIR, "brain_regions_at_least_one_task_mean_gt_2.csv"),
    os.path.join(DEFAULT_CSV_DIR, "brain_regions_at_least_two_tasks_mean_gt_2.csv"),
    os.path.join(DEFAULT_CSV_DIR, "alphabetic_structural.csv"),
    os.path.join(DEFAULT_CSV_DIR, "morphosyllabic_functional.csv"),
    os.path.join(DEFAULT_CSV_DIR, "morphosyllabic_structural.csv"),
    os.path.join(DEFAULT_CSV_DIR, "alphabetic_functional.csv"),
]
PLOT_VIEWS = [
    ("lh", "lateral", "left"),
    ("lh", "medial", "left"),
    ("rh", "lateral", "right"),
    ("rh", "medial", "right"),
]


def load_fsaverage_and_annot() -> Dict[str, object]:
    subjects_dir = os.path.join(BASE_DIR, "FILE")
    fs_dir = mne.datasets.fetch_fsaverage(subjects_dir=subjects_dir, verbose=True)
    print("fsaverage dir:", fs_dir)

    lh_mesh = read_geometry(os.path.join(fs_dir, "surf", "lh.inflated"))
    rh_mesh = read_geometry(os.path.join(fs_dir, "surf", "rh.inflated"))
    lh_sulc = read_morph_data(os.path.join(fs_dir, "surf", "lh.sulc"))
    rh_sulc = read_morph_data(os.path.join(fs_dir, "surf", "rh.sulc"))

    lh_labels, _, lh_names = read_annot(LH_ANNOT)
    rh_labels, _, rh_names = read_annot(RH_ANNOT)

    if lh_mesh[0].shape[0] != len(lh_labels) or rh_mesh[0].shape[0] != len(rh_labels):
        raise RuntimeError(
            "Surface 顶点数与 annot 不匹配，请确认使用的是 fsaverage 而不是 fsaverage5/fsaverage6。"
        )

    return {
        "lh_mesh": lh_mesh,
        "rh_mesh": rh_mesh,
        "lh_sulc": lh_sulc,
        "rh_sulc": rh_sulc,
        "lh_labels": lh_labels,
        "rh_labels": rh_labels,
        "lh_names": lh_names,
        "rh_names": rh_names,
    }


def strip_hemi_prefix(value: str) -> str:
    value = str(value)
    return value.replace("lh.", "").replace("rh.", "").replace("lh_", "").replace("rh_", "")


def build_roi_lookup_keys(feature: str) -> List[str]:
    feature = str(feature).strip()
    normalized = feature.replace(".", "_")
    return [
        feature,
        normalized,
        strip_hemi_prefix(feature),
        strip_hemi_prefix(normalized),
    ]


def make_vertex_data(labels, names, roi_map: Dict[str, float]) -> Tuple[np.ndarray, List[str]]:
    name_list = [name.decode("utf-8") for name in names]
    vertex_data = np.full(labels.shape[0], np.nan, dtype=float)

    matched_names = []
    for label_id, roi_name in enumerate(name_list):
        if roi_name in roi_map:
            vertex_data[labels == label_id] = float(roi_map[roi_name])
            matched_names.append(roi_name)

    return vertex_data, matched_names


def build_single_color_cmap(color: str):
    return LinearSegmentedColormap.from_list(
        "table_region_single_color",
        ["#d8d8d8", color],
    )


def has_any_finite_value(vertex_data: np.ndarray) -> bool:
    return bool(np.isfinite(vertex_data).any())


def sanitize_name(path: str) -> str:
    return os.path.splitext(os.path.basename(path))[0]


def load_table_regions(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    if df.shape[1] == 0:
        raise ValueError(f"CSV 没有任何列: {csv_path}")

    first_col = df.columns[0]
    plot_df = pd.DataFrame()
    plot_df["feature"] = df.iloc[:, 0].astype(str).str.strip()
    plot_df = plot_df[plot_df["feature"].ne("")].copy()
    plot_df = plot_df.drop_duplicates(subset=["feature"]).reset_index(drop=True)
    plot_df["value"] = 1.0
    plot_df["source_first_column"] = first_col
    return plot_df


def build_hemi_maps(plot_df: pd.DataFrame):
    lh_df = plot_df.loc[
        plot_df["feature"].str.startswith("lh_") | plot_df["feature"].str.startswith("lh.")
    ].copy()
    rh_df = plot_df.loc[
        plot_df["feature"].str.startswith("rh_") | plot_df["feature"].str.startswith("rh.")
    ].copy()

    lh_map: Dict[str, float] = {}
    for _, row in lh_df.iterrows():
        for key in build_roi_lookup_keys(row["feature"]):
            lh_map[key] = row["value"]

    rh_map: Dict[str, float] = {}
    for _, row in rh_df.iterrows():
        for key in build_roi_lookup_keys(row["feature"]):
            rh_map[key] = row["value"]

    return lh_df, rh_df, lh_map, rh_map


def render_region_table_brainmaps(csv_path: str, fs_data: Dict[str, object], output_dir: str, color: str):
    table_name = sanitize_name(csv_path)
    table_out_dir = os.path.join(output_dir, table_name)
    os.makedirs(table_out_dir, exist_ok=True)

    plot_df = load_table_regions(csv_path)
    lh_df, rh_df, lh_map, rh_map = build_hemi_maps(plot_df)

    lh_vtx, lh_matched_names = make_vertex_data(fs_data["lh_labels"], fs_data["lh_names"], lh_map)
    rh_vtx, rh_matched_names = make_vertex_data(fs_data["rh_labels"], fs_data["rh_names"], rh_map)

    lh_unmatched = sorted(set(lh_map) - set(lh_matched_names))
    rh_unmatched = sorted(set(rh_map) - set(rh_matched_names))

    print(f"\n========== {table_name} ==========")
    print(f"CSV: {csv_path}")
    print(f"First column used: {plot_df['source_first_column'].iloc[0] if not plot_df.empty else 'N/A'}")
    print(f"Total requested regions: {len(plot_df)}")
    print(f"LH requested ROIs: {int(lh_df.shape[0])}, matched annot names: {len(lh_matched_names)}")
    print(f"RH requested ROIs: {int(rh_df.shape[0])}, matched annot names: {len(rh_matched_names)}")
    print(f"LH unmatched examples: {lh_unmatched[:10]}")
    print(f"RH unmatched examples: {rh_unmatched[:10]}")

    export_path = os.path.join(table_out_dir, f"{table_name}_regions_used.csv")
    plot_df[["feature"]].to_csv(export_path, index=False)
    print("Saved:", export_path)

    cmap = build_single_color_cmap(color)
    vertex_data_by_hemi = {"lh": lh_vtx, "rh": rh_vtx}
    mesh_by_hemi = {"lh": fs_data["lh_mesh"], "rh": fs_data["rh_mesh"]}
    sulc_by_hemi = {"lh": fs_data["lh_sulc"], "rh": fs_data["rh_sulc"]}

    for hemi_tag, view, hemi_name in PLOT_VIEWS:
        out_file = os.path.join(table_out_dir, f"{table_name}_{hemi_tag}_{view}.png")
        fig = plt.figure(figsize=(4.8, 4.2), facecolor="white")
        axis = fig.add_subplot(111, projection="3d")

        vertex_data = vertex_data_by_hemi[hemi_tag]
        if has_any_finite_value(vertex_data):
            plotting.plot_surf_stat_map(
                surf_mesh=mesh_by_hemi[hemi_tag],
                stat_map=vertex_data,
                hemi=hemi_name,
                view=view,
                bg_map=sulc_by_hemi[hemi_tag],
                bg_on_data=True,
                colorbar=False,
                symmetric_cbar=False,
                threshold=0.5,
                vmin=0.0,
                vmax=1.0,
                cmap=cmap,
                axes=axis,
            )
        else:
            print(f"No matched regions for {table_name} {hemi_tag}; saving background-only surface.")
            plotting.plot_surf(
                surf_mesh=mesh_by_hemi[hemi_tag],
                hemi=hemi_name,
                view=view,
                bg_map=sulc_by_hemi[hemi_tag],
                colorbar=False,
                axes=axis,
            )

        fig.savefig(out_file, bbox_inches="tight", pad_inches=0.0, dpi=300, facecolor="white")
        plt.close(fig)
        print("Saved:", out_file)

    summary_row = {
        "table_name": table_name,
        "csv_path": csv_path,
        "first_column_name": plot_df["source_first_column"].iloc[0] if not plot_df.empty else "",
        "n_regions": int(plot_df.shape[0]),
        "n_lh": int(lh_df.shape[0]),
        "n_rh": int(rh_df.shape[0]),
        "n_lh_matched": len(lh_matched_names),
        "n_rh_matched": len(rh_matched_names),
    }
    return summary_row


def parse_args():
    parser = argparse.ArgumentParser(
        description="根据 CSV 第一列中的 brain region 列表绘制 DK318 脑图。"
    )
    parser.add_argument(
        "--csvs",
        nargs="+",
        default=DEFAULT_CSVS,
        help="一个或多个 CSV 路径。脚本会直接读取第一列作为脑区名。",
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help="输出目录。每个 CSV 会生成一个同名子目录。",
    )
    parser.add_argument(
        "--color",
        default="#ff5a7a",
        help="脑区高亮颜色，默认是粉红色。",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print("当前工作目录:", BASE_DIR)
    fs_data = load_fsaverage_and_annot()

    summary_rows = []
    for csv_path in args.csvs:
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"CSV 不存在: {csv_path}")

        summary = render_region_table_brainmaps(
            csv_path=csv_path,
            fs_data=fs_data,
            output_dir=args.output_dir,
            color=args.color,
        )
        summary_rows.append(summary)

    summary_df = pd.DataFrame(summary_rows)
    summary_csv = os.path.join(args.output_dir, "brain_region_table_brainmaps_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    print("Saved:", summary_csv)
    print("\n脑图绘制完成。")


if __name__ == "__main__":
    main()
