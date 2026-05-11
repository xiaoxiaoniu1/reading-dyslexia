#!/usr/bin/env python3
"""
DD-only NMF subtype brain maps.

This script visualizes DD-only NMF component ROI weights on the DK318 surface.
It reads DD-only NMF outputs and writes per-component brain maps, summary maps,
and top-ROI barplots.
"""

import os
import re
import numpy as np
import pandas as pd
import mne
from nibabel.freesurfer.io import read_annot, read_geometry, read_morph_data
from nilearn import plotting
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

BASE_DIR = "/data/home/tqi/data1/share/after_freesurfer"
NMF_DIR = os.environ.get(
    "DD_NMF_OUT_DIR",
    os.path.join(BASE_DIR, "FILE/test_mean_1.5/MIND_DK318_NMF_DD_only"),
)
H_MATRIX_FILE = os.path.join(NMF_DIR, "NMF_DD_only_H_component_roi_weights_long.csv")
TOP_ROIS_FILE = os.path.join(NMF_DIR, "NMF_DD_only_component_top_ROIs.csv")
OUTPUT_DIR = os.path.join(NMF_DIR, "brain_maps")
LH_ANNOT = os.path.join(BASE_DIR, "FILE/DK-318", "lh.DK318.annot")
RH_ANNOT = os.path.join(BASE_DIR, "FILE/DK-318", "rh.DK318.annot")

VIEWS = [
    ("lh", "lateral", "left"),
    ("rh", "lateral", "right"),
    ("lh", "medial", "left"),
    ("rh", "medial", "right"),
]
COLORS = ["#e74c3c", "#3498db", "#2ecc71", "#f39c12", "#9b59b6", "#1abc9c", "#34495e", "#e67e22"]


def component_sort_key(component_name: str):
    match = re.search(r"(\d+)$", str(component_name))
    return int(match.group(1)) if match else str(component_name)


def sanitize_name(value: str) -> str:
    value = str(value)
    value = re.sub(r"[^A-Za-z0-9_\-]+", "_", value)
    return value.strip("_")


def build_component_configs(components):
    configs = []
    for idx, component in enumerate(sorted(components, key=component_sort_key)):
        comp_num = component_sort_key(component)
        if isinstance(comp_num, int):
            short_name = f"DD_Subtype_{comp_num}"
            display_name = f"DD-only Subtype Component {comp_num}"
        else:
            short_name = sanitize_name(component)
            display_name = f"DD-only {component}"
        configs.append(
            {
                "component_name": component,
                "display_name": display_name,
                "short_name": short_name,
                "color": COLORS[idx % len(COLORS)],
            }
        )
    return configs


def load_fsaverage_and_annot():
    print("Loading fsaverage and DK318 annotations...")
    subjects_dir = os.path.join(BASE_DIR, "FILE")
    fs_dir = mne.datasets.fetch_fsaverage(subjects_dir=subjects_dir, verbose=False)

    lh_mesh = read_geometry(os.path.join(fs_dir, "surf", "lh.inflated"))
    rh_mesh = read_geometry(os.path.join(fs_dir, "surf", "rh.inflated"))
    lh_sulc = read_morph_data(os.path.join(fs_dir, "surf", "lh.sulc"))
    rh_sulc = read_morph_data(os.path.join(fs_dir, "surf", "rh.sulc"))
    lh_labels, _, lh_names = read_annot(LH_ANNOT)
    rh_labels, _, rh_names = read_annot(RH_ANNOT)

    print(f"  LH vertices: {lh_mesh[0].shape[0]}, RH vertices: {rh_mesh[0].shape[0]}")
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
    return str(value).replace("lh.", "").replace("rh.", "").replace("lh_", "").replace("rh_", "")


def build_roi_lookup_keys(roi_name: str):
    roi_name = str(roi_name)
    normalized = roi_name.replace(".", "_")
    return [roi_name, normalized, strip_hemi_prefix(roi_name), strip_hemi_prefix(normalized)]


def make_vertex_data(labels, names, roi_map):
    name_list = [name.decode("utf-8") for name in names]
    vertex_data = np.zeros(labels.shape[0], dtype=float)
    matched_names = []
    for label_id, roi_name in enumerate(name_list):
        if roi_name in roi_map:
            vertex_data[labels == label_id] = float(roi_map[roi_name])
            matched_names.append(roi_name)
    return vertex_data, matched_names


def build_gradient_cmap(color: str):
    return LinearSegmentedColormap.from_list("component_gradient", ["#ffffff", color])


def load_nmf_data():
    print(f"\nLoading DD-only NMF data from: {NMF_DIR}")
    if not os.path.exists(H_MATRIX_FILE):
        raise FileNotFoundError(f"Missing H matrix file: {H_MATRIX_FILE}")
    if not os.path.exists(TOP_ROIS_FILE):
        raise FileNotFoundError(f"Missing top ROIs file: {TOP_ROIS_FILE}")

    h_df = pd.read_csv(H_MATRIX_FILE)
    top_df = pd.read_csv(TOP_ROIS_FILE)
    components = sorted(h_df["component"].unique(), key=component_sort_key)
    print(f"  Components: {components}")
    for comp in components:
        count = len(h_df[h_df["component"] == comp])
        print(f"    {comp}: {count} ROIs")
    return h_df, top_df, components


def prepare_component_dataframe(h_df, component_name, weight_threshold=0.05):
    comp_data = h_df[h_df["component"] == component_name].copy()
    if weight_threshold > 0:
        comp_data = comp_data[comp_data["weight"] >= weight_threshold]
    print(f"  {component_name}: {len(comp_data)} ROIs (threshold={weight_threshold})")
    return comp_data


def build_hemi_maps(comp_data):
    roi_series = comp_data["ROI"].astype(str)
    lh_rois = comp_data[roi_series.str.startswith("lh_") | roi_series.str.startswith("lh.")].copy()
    rh_rois = comp_data[roi_series.str.startswith("rh_") | roi_series.str.startswith("rh.")].copy()

    lh_map = {}
    for _, row in lh_rois.iterrows():
        for key in build_roi_lookup_keys(row["ROI"]):
            lh_map[key] = row["weight"]

    rh_map = {}
    for _, row in rh_rois.iterrows():
        for key in build_roi_lookup_keys(row["ROI"]):
            rh_map[key] = row["weight"]

    return lh_rois, rh_rois, lh_map, rh_map


def render_component_views(config, comp_data, fs_data, output_dir, vmax=None):
    print(f"\n{'=' * 60}\nRendering {config['display_name']}\n{'=' * 60}")
    os.makedirs(output_dir, exist_ok=True)

    lh_rois, rh_rois, lh_map, rh_map = build_hemi_maps(comp_data)
    lh_vtx, lh_matched = make_vertex_data(fs_data["lh_labels"], fs_data["lh_names"], lh_map)
    rh_vtx, rh_matched = make_vertex_data(fs_data["rh_labels"], fs_data["rh_names"], rh_map)

    print(f"  LH ROIs: {len(lh_rois)}, matched: {len(lh_matched)}")
    print(f"  RH ROIs: {len(rh_rois)}, matched: {len(rh_matched)}")

    if vmax is None:
        vmax = max(lh_vtx.max(), rh_vtx.max())

    comp_data_sorted = comp_data.sort_values("weight", ascending=False)
    roi_list_path = os.path.join(output_dir, f"{config['short_name']}_ROIs.csv")
    comp_data_sorted[["ROI", "weight", "rank_within_component"]].to_csv(roi_list_path, index=False)
    print(f"  Saved ROI list: {roi_list_path}")

    cmap = build_gradient_cmap(config["color"])
    vertex_data_by_hemi = {"lh": lh_vtx, "rh": rh_vtx}
    mesh_by_hemi = {"lh": fs_data["lh_mesh"], "rh": fs_data["rh_mesh"]}
    sulc_by_hemi = {"lh": fs_data["lh_sulc"], "rh": fs_data["rh_sulc"]}

    for hemi_tag, view, hemi_name in VIEWS:
        out_file = os.path.join(output_dir, f"{config['short_name']}_{hemi_tag}_{view}.png")
        fig = plt.figure(figsize=(6, 5), facecolor="white")
        axis = fig.add_subplot(111, projection="3d")
        plotting.plot_surf_stat_map(
            surf_mesh=mesh_by_hemi[hemi_tag],
            stat_map=vertex_data_by_hemi[hemi_tag],
            hemi=hemi_name,
            view=view,
            bg_map=sulc_by_hemi[hemi_tag],
            bg_on_data=True,
            colorbar=True,
            symmetric_cbar=False,
            threshold=0.01,
            vmin=0.0,
            vmax=vmax,
            cmap=cmap,
            axes=axis,
            darkness=0.55,
        )
        fig.suptitle(f"{config['display_name']}\n{hemi_tag.upper()} {view.capitalize()} View", fontsize=11, y=0.98)
        fig.savefig(out_file, bbox_inches="tight", pad_inches=0.1, dpi=300, facecolor="white")
        plt.close(fig)
        print(f"  Saved: {os.path.basename(out_file)}")

    return {"config": config, "lh_vtx": lh_vtx, "rh_vtx": rh_vtx, "vmax": vmax}


def create_summary_figure(component_results, fs_data):
    print(f"\n{'=' * 60}\nCreating summary figure...\n{'=' * 60}")
    summary_path = os.path.join(OUTPUT_DIR, "DD_only_All_Subtype_Components_Summary.png")
    n_components = len(component_results)
    fig = plt.figure(figsize=(12, 3 * n_components), facecolor="white")

    for row_idx, result in enumerate(component_results):
        config = result["config"]
        vmax = result["vmax"]
        ax_lh = fig.add_subplot(n_components, 2, row_idx * 2 + 1, projection="3d")
        plotting.plot_surf_stat_map(
            surf_mesh=fs_data["lh_mesh"],
            stat_map=result["lh_vtx"],
            hemi="left",
            view="lateral",
            bg_map=fs_data["lh_sulc"],
            bg_on_data=True,
            colorbar=True,
            threshold=0.01,
            vmin=0.0,
            vmax=vmax,
            cmap=build_gradient_cmap(config["color"]),
            axes=ax_lh,
            darkness=0.55,
        )
        ax_lh.set_title(f"{config['display_name']} - LH Lateral", fontsize=10)

        ax_rh = fig.add_subplot(n_components, 2, row_idx * 2 + 2, projection="3d")
        plotting.plot_surf_stat_map(
            surf_mesh=fs_data["rh_mesh"],
            stat_map=result["rh_vtx"],
            hemi="right",
            view="lateral",
            bg_map=fs_data["rh_sulc"],
            bg_on_data=True,
            colorbar=True,
            threshold=0.01,
            vmin=0.0,
            vmax=vmax,
            cmap=build_gradient_cmap(config["color"]),
            axes=ax_rh,
            darkness=0.55,
        )
        ax_rh.set_title(f"{config['display_name']} - RH Lateral", fontsize=10)

    plt.tight_layout()
    fig.savefig(summary_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {summary_path}")


def create_top_rois_barplots(top_df, configs, output_dir):
    print(f"\n{'=' * 60}\nCreating DD-only top ROI bar plots...\n{'=' * 60}")
    top20_df = top_df[top_df["top_type"] == "top20"].copy()

    for config in configs:
        comp_name = config["component_name"]
        comp_top = top20_df[top20_df["component"] == comp_name].copy()
        if len(comp_top) == 0:
            print(f"  No top20 data for {comp_name}")
            continue

        comp_top = comp_top.sort_values("weight", ascending=True)
        fig, ax = plt.subplots(figsize=(8, 10), facecolor="white")
        ax.barh(range(len(comp_top)), comp_top["weight"], color=config["color"], alpha=0.7)
        ax.set_yticks(range(len(comp_top)))
        ax.set_yticklabels(comp_top["ROI"], fontsize=9)
        ax.set_xlabel("Weight", fontsize=11)
        ax.set_title(f"{config['display_name']}\nTop 20 ROIs by Weight", fontsize=12, fontweight="bold")
        ax.grid(axis="x", alpha=0.3)

        for i, (_, row) in enumerate(comp_top.iterrows()):
            ax.text(row["weight"], i, f" {row['weight']:.3f}", va="center", fontsize=8)

        plt.tight_layout()
        out_file = os.path.join(output_dir, f"{config['short_name']}_top20_barplot.png")
        fig.savefig(out_file, dpi=300, facecolor="white", bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {os.path.basename(out_file)}")


def main():
    print("=" * 60 + "\nDD-only NMF Subtype Brain Visualization\n" + "=" * 60)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    fs_data = load_fsaverage_and_annot()
    h_df, top_df, components = load_nmf_data()
    configs = build_component_configs(components)

    global_vmax = h_df["weight"].max()
    print(f"\nGlobal weight max: {global_vmax:.3f}")

    component_results = []
    for config in configs:
        comp_data = prepare_component_dataframe(h_df, config["component_name"], weight_threshold=0.05)
        comp_output_dir = os.path.join(OUTPUT_DIR, config["short_name"])
        result = render_component_views(config, comp_data, fs_data, comp_output_dir, vmax=global_vmax)
        component_results.append(result)

    if component_results:
        create_summary_figure(component_results, fs_data)

    create_top_rois_barplots(top_df, configs, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("DD-only brain visualizations completed!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
