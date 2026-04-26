import os
import numpy as np
import pandas as pd

import mne
from nibabel.freesurfer.io import read_annot, read_geometry, read_morph_data
from nilearn import plotting
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


BASE_DIR = "/data/home/tqi/data1/share/after_freesurfer"
NETWORK_DIR = os.path.join(
    BASE_DIR,
    "FILE/test_mean_1.5/DM_MD_LAN_network",
)
OUTPUT_DIR = os.path.join(NETWORK_DIR, "DK318_network_brainmaps")
LH_ANNOT = os.path.join(BASE_DIR, "FILE/DK-318", "lh.DK318.annot")
RH_ANNOT = os.path.join(BASE_DIR, "FILE/DK-318", "rh.DK318.annot")

NETWORK_CONFIGS = [
    {
        "csv_path": os.path.join(NETWORK_DIR, "DK318_Kong17_Language_hard_max4.csv"),
        "display_name": "Language network",
        "short_name": "LAN",
        "subtitle": "Language knowledge and processing",
        "color": "#c78cff",
    },
    {
        "csv_path": os.path.join(NETWORK_DIR, "DK318_Kong17_MD_hard_max4.csv"),
        "display_name": "Multiple demand network",
        "short_name": "MD",
        "subtitle": "Response inhibition, working memory and goal-directed behavior",
        "color": "#22d9ff",
    },
    {
        "csv_path": os.path.join(NETWORK_DIR, "DK318_Kong17_DM_hard_max4.csv"),
        "display_name": "Default mode network",
        "short_name": "DM",
        "subtitle": "Mind-wandering, self-referencing, autobiographical memories, and abstract thoughts",
        "color": "#2fe6b0",
    },
]

SUMMARY_PLOTS = [
    ("lh", "lateral", "left"),
    ("rh", "lateral", "right"),
]


def load_fsaverage_and_annot():
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
            "表面顶点数与 annot 不匹配，请确认使用的是 fsaverage 而不是 fsaverage5/fsaverage6。"
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


def build_roi_lookup_keys(feature: str):
    feature = str(feature)
    normalized = feature.replace(".", "_")
    return [
        feature,
        normalized,
        strip_hemi_prefix(feature),
        strip_hemi_prefix(normalized),
    ]


def make_vertex_data(labels, names, roi_map):
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
        "network_single_color",
        ["#d8d8d8", color],
    )


def prepare_network_dataframe(csv_path: str, short_name: str):
    df = pd.read_csv(csv_path)

    required_columns = {"parcel_id"}
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        raise ValueError(f"CSV 缺少必要列: {sorted(missing_columns)}")

    plot_df = df.copy()
    plot_df["feature"] = plot_df["parcel_id"].astype(str)
    plot_df["network_display_name"] = short_name
    plot_df["value"] = 1.0
    plot_df = plot_df.drop_duplicates(subset=["feature"])

    return plot_df[["network_display_name", "feature", "value"]].copy()


def build_hemi_maps(plot_df):
    lh_df = plot_df.loc[
        plot_df["feature"].str.startswith("lh_") | plot_df["feature"].str.startswith("lh.")
    ].copy()
    rh_df = plot_df.loc[
        plot_df["feature"].str.startswith("rh_") | plot_df["feature"].str.startswith("rh.")
    ].copy()

    lh_map = {}
    for _, row in lh_df.iterrows():
        for key in build_roi_lookup_keys(row["feature"]):
            lh_map[key] = row["value"]

    rh_map = {}
    for _, row in rh_df.iterrows():
        for key in build_roi_lookup_keys(row["feature"]):
            rh_map[key] = row["value"]

    return lh_df, rh_df, lh_map, rh_map


def render_network_views(config, plot_df, fs_data, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    lh_df, rh_df, lh_map, rh_map = build_hemi_maps(plot_df)

    lh_vtx, lh_matched_names = make_vertex_data(fs_data["lh_labels"], fs_data["lh_names"], lh_map)
    rh_vtx, rh_matched_names = make_vertex_data(fs_data["rh_labels"], fs_data["rh_names"], rh_map)

    lh_unmatched = sorted(set(lh_map) - set(lh_matched_names))
    rh_unmatched = sorted(set(rh_map) - set(rh_matched_names))

    print(f"\n========== {config['short_name']} ==========")
    print(f"CSV: {config['csv_path']}")
    print(f"LH requested ROIs: {int(lh_df.shape[0])}, matched annot names: {len(lh_matched_names)}")
    print(f"RH requested ROIs: {int(rh_df.shape[0])}, matched annot names: {len(rh_matched_names)}")
    print(f"LH unmatched examples: {lh_unmatched[:10]}")
    print(f"RH unmatched examples: {rh_unmatched[:10]}")

    export_path = os.path.join(output_dir, f"{config['short_name']}_regions.csv")
    plot_df[["network_display_name", "feature"]].to_csv(export_path, index=False)
    print("Saved:", export_path)

    cmap = build_single_color_cmap(config["color"])

    vertex_data_by_hemi = {
        "lh": lh_vtx,
        "rh": rh_vtx,
    }
    mesh_by_hemi = {
        "lh": fs_data["lh_mesh"],
        "rh": fs_data["rh_mesh"],
    }
    sulc_by_hemi = {
        "lh": fs_data["lh_sulc"],
        "rh": fs_data["rh_sulc"],
    }

    for hemi_tag, view, hemi_name in SUMMARY_PLOTS:
        out_file = os.path.join(output_dir, f"{config['short_name']}_{hemi_tag}_{view}.png")
        fig = plt.figure(figsize=(4.6, 4.0), facecolor="white")
        axis = fig.add_subplot(111, projection="3d")

        plotting.plot_surf_stat_map(
            surf_mesh=mesh_by_hemi[hemi_tag],
            stat_map=vertex_data_by_hemi[hemi_tag],
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
            darkness=0.55,
        )

        fig.savefig(out_file, bbox_inches="tight", pad_inches=0.0, dpi=300, facecolor="white")
        plt.close(fig)
        print("Saved:", out_file)

    return {
        "config": config,
        "plot_df": plot_df,
        "lh_vtx": lh_vtx,
        "rh_vtx": rh_vtx,
    }


def create_summary_figure(network_results, fs_data):
    summary_path = os.path.join(OUTPUT_DIR, "DK318_DM_MD_LAN_networks_summary.png")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    fig = plt.figure(figsize=(9.6, 8.4), facecolor="white")
    grid = fig.add_gridspec(
        nrows=3,
        ncols=2,
        left=0.28,
        right=0.98,
        top=0.97,
        bottom=0.06,
        wspace=0.02,
        hspace=0.16,
    )

    for row_idx, result in enumerate(network_results):
        config = result["config"]
        row_top = 0.97 - row_idx * 0.305

        fig.text(
            0.03,
            row_top - 0.018,
            config["display_name"],
            fontsize=17,
            fontweight="bold",
            color="black",
            ha="left",
            va="top",
        )
        fig.text(
            0.03,
            row_top - 0.050,
            config["subtitle"],
            fontsize=10.5,
            color="black",
            ha="left",
            va="top",
        )

        for col_idx, (hemi_tag, view, hemi_name) in enumerate(SUMMARY_PLOTS):
            axis = fig.add_subplot(grid[row_idx, col_idx], projection="3d")
            mesh = fs_data["lh_mesh"] if hemi_tag == "lh" else fs_data["rh_mesh"]
            sulc = fs_data["lh_sulc"] if hemi_tag == "lh" else fs_data["rh_sulc"]
            vertex_data = result["lh_vtx"] if hemi_tag == "lh" else result["rh_vtx"]

            plotting.plot_surf_stat_map(
                surf_mesh=mesh,
                stat_map=vertex_data,
                hemi=hemi_name,
                view=view,
                bg_map=sulc,
                bg_on_data=True,
                colorbar=False,
                symmetric_cbar=False,
                threshold=0.5,
                vmin=0.0,
                vmax=1.0,
                cmap=build_single_color_cmap(config["color"]),
                axes=axis,
                darkness=0.55,
            )

    fig.savefig(summary_path, dpi=300, bbox_inches="tight", pad_inches=0.03, facecolor="white")
    plt.close(fig)
    print("Saved:", summary_path)


def main():
    print("当前工作目录：", BASE_DIR)
    fs_data = load_fsaverage_and_annot()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    summary_rows = []
    network_results = []

    for config in NETWORK_CONFIGS:
        plot_df = prepare_network_dataframe(
            csv_path=config["csv_path"],
            short_name=config["short_name"],
        )

        network_out_dir = os.path.join(OUTPUT_DIR, config["short_name"])
        result = render_network_views(
            config=config,
            plot_df=plot_df,
            fs_data=fs_data,
            output_dir=network_out_dir,
        )
        network_results.append(result)

        summary_rows.append(
            {
                "network": config["short_name"],
                "display_name": config["display_name"],
                "csv_path": config["csv_path"],
                "n_regions": int(plot_df.shape[0]),
                "n_lh": int(
                    plot_df["feature"].str.startswith("lh_").sum()
                    + plot_df["feature"].str.startswith("lh.").sum()
                ),
                "n_rh": int(
                    plot_df["feature"].str.startswith("rh_").sum()
                    + plot_df["feature"].str.startswith("rh.").sum()
                ),
            }
        )

    create_summary_figure(network_results, fs_data)

    summary_df = pd.DataFrame(summary_rows)
    summary_csv = os.path.join(OUTPUT_DIR, "network_region_counts.csv")
    summary_df.to_csv(summary_csv, index=False)
    print("Saved:", summary_csv)
    print("\n三个网络脑图绘制完成。")


if __name__ == "__main__":
    main()
