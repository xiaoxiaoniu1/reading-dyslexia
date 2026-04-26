import argparse
import os
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import mne
import numpy as np
from nibabel.freesurfer.io import read_annot, read_geometry
from nilearn import plotting


BASE_DIR = "/data/home/tqi/data1/share/after_freesurfer"
DEFAULT_OUT_DIR = os.path.join(BASE_DIR, "FILE", "test_mean_1.5", "region_locator_plots")
DK318_LH_ANNOT = os.path.join(BASE_DIR, "FILE", "DK-318", "lh.DK318.annot")
DK318_RH_ANNOT = os.path.join(BASE_DIR, "FILE", "DK-318", "rh.DK318.annot")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="输入 DK68 或 DK318 的脑区名，直接高亮画出该脑区位置。"
    )
    parser.add_argument("--atlas", choices=["DK68", "DK318"], required=True)
    parser.add_argument(
        "--region",
        required=True,
        help="脑区名。可输入完整名称如 lh_precuneus_part7，也可只输入 precuneus_part7。",
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUT_DIR,
        help="输出目录，默认放在 FILE/test_mean_1.5/region_locator_plots",
    )
    parser.add_argument(
        "--list-only",
        action="store_true",
        help="只列出该 atlas 的全部脑区名，不画图。",
    )
    return parser.parse_args()


def ensure_exists(path: str) -> None:
    if not os.path.exists(path):
        raise FileNotFoundError(f"文件不存在: {path}")


def sanitize_name(text: str) -> str:
    keep = []
    for ch in str(text):
        if ch.isalnum() or ch in {"-", "_", "."}:
            keep.append(ch)
        else:
            keep.append("_")
    return "".join(keep).strip("_")


def normalize_region_query(region: str) -> Tuple[str, str]:
    text = str(region).strip()
    if text.startswith("lh_"):
        return "lh", text[3:]
    if text.startswith("rh_"):
        return "rh", text[3:]
    return "both", text


def fetch_fsaverage_surfaces() -> Tuple[Tuple[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray], str]:
    subjects_dir = os.path.join(BASE_DIR, "FILE")
    fs_dir = mne.datasets.fetch_fsaverage(subjects_dir=subjects_dir, verbose=True)
    lh_surf_path = os.path.join(fs_dir, "surf", "lh.inflated")
    rh_surf_path = os.path.join(fs_dir, "surf", "rh.inflated")
    lh_coords, lh_faces = read_geometry(lh_surf_path)
    rh_coords, rh_faces = read_geometry(rh_surf_path)
    return (lh_coords, lh_faces), (rh_coords, rh_faces), fs_dir


def load_atlas(atlas: str) -> Dict[str, object]:
    lh_mesh, rh_mesh, fs_dir = fetch_fsaverage_surfaces()

    if atlas == "DK68":
        lh_annot = os.path.join(fs_dir, "label", "lh.aparc.annot")
        rh_annot = os.path.join(fs_dir, "label", "rh.aparc.annot")
    else:
        lh_annot = DK318_LH_ANNOT
        rh_annot = DK318_RH_ANNOT

    ensure_exists(lh_annot)
    ensure_exists(rh_annot)

    lh_labels, _, lh_names = read_annot(lh_annot)
    rh_labels, _, rh_names = read_annot(rh_annot)

    if len(lh_labels) != lh_mesh[0].shape[0] or len(rh_labels) != rh_mesh[0].shape[0]:
        raise RuntimeError("annot 顶点数与 fsaverage 表面顶点数不匹配。")

    return {
        "lh_mesh": lh_mesh,
        "rh_mesh": rh_mesh,
        "lh_labels": lh_labels,
        "rh_labels": rh_labels,
        "lh_names": [n.decode("utf-8") for n in lh_names],
        "rh_names": [n.decode("utf-8") for n in rh_names],
    }


def list_all_regions(atlas_data: Dict[str, object]) -> List[str]:
    regions = [f"lh_{name}" for name in atlas_data["lh_names"]]
    regions.extend(f"rh_{name}" for name in atlas_data["rh_names"])
    return sorted(regions)


def find_matching_regions(atlas_data: Dict[str, object], region_query: str) -> List[Tuple[str, str]]:
    hemi_query, region_name = normalize_region_query(region_query)
    matches: List[Tuple[str, str]] = []

    if hemi_query in {"lh", "both"}:
        for name in atlas_data["lh_names"]:
            if name == region_name:
                matches.append(("lh", name))

    if hemi_query in {"rh", "both"}:
        for name in atlas_data["rh_names"]:
            if name == region_name:
                matches.append(("rh", name))

    return matches


def build_vertex_map(labels: np.ndarray, names: List[str], target_region: str) -> np.ndarray:
    vertex_map = np.full(labels.shape[0], np.nan, dtype=float)
    for lab_id, roi_name in enumerate(names):
        if roi_name == target_region:
            vertex_map[labels == lab_id] = 1.0
            break
    return vertex_map


def render_single_region(
    atlas_data: Dict[str, object],
    hemi: str,
    region_name: str,
    output_prefix: str,
) -> List[str]:
    if hemi == "lh":
        mesh = atlas_data["lh_mesh"]
        labels = atlas_data["lh_labels"]
        names = atlas_data["lh_names"]
        hemi_name = "left"
    else:
        mesh = atlas_data["rh_mesh"]
        labels = atlas_data["rh_labels"]
        names = atlas_data["rh_names"]
        hemi_name = "right"

    vertex_map = build_vertex_map(labels, names, region_name)
    if not np.isfinite(vertex_map).any():
        raise ValueError(f"annot 中没有找到脑区: {hemi}_{region_name}")

    saved_files: List[str] = []
    for view in ["lateral", "medial"]:
        out_file = f"{output_prefix}_{hemi}_{view}.png"
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection="3d")
        display = plotting.plot_surf_stat_map(
            surf_mesh=mesh,
            stat_map=vertex_map,
            hemi=hemi_name,
            view=view,
            colorbar=False,
            vmin=0.0,
            vmax=1.0,
            cmap="autumn",
            title=f"{hemi}_{region_name} | {view}",
            axes=ax,
            darkness=None,
        )
        try:
            display.add_contours(labels, colors="k", linewidths=0.4)
        except AttributeError:
            unique_labels = np.unique(labels)
            plotting.plot_surf_contours(
                surf_mesh=mesh,
                roi_map=labels,
                hemi=hemi_name,
                view=view,
                axes=ax,
                levels=unique_labels,
                colors=["black"] * len(unique_labels),
                linewidths=0.15,
            )
        fig.savefig(out_file, bbox_inches="tight", pad_inches=0.2)
        plt.close(fig)
        saved_files.append(out_file)

    return saved_files


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    atlas_data = load_atlas(args.atlas)

    if args.list_only:
        atlas_out_dir = os.path.join(args.output_dir, args.atlas)
        os.makedirs(atlas_out_dir, exist_ok=True)
        region_list = list_all_regions(atlas_data)
        list_path = os.path.join(atlas_out_dir, f"{args.atlas}_all_regions.txt")
        with open(list_path, "w", encoding="utf-8") as f:
            for item in region_list:
                f.write(f"{item}\n")
        print(f"已保存全部脑区名: {list_path}")
        for item in region_list:
            print(item)
        return

    matches = find_matching_regions(atlas_data, args.region)
    if not matches:
        print(f"没有找到脑区: {args.region}")
        print("你可以先用 --list-only 看全部可用脑区名。")
        return

    atlas_out_dir = os.path.join(args.output_dir, args.atlas, sanitize_name(args.region))
    os.makedirs(atlas_out_dir, exist_ok=True)

    print("找到以下脑区：")
    for hemi, region_name in matches:
        print(f"- {hemi}_{region_name}")

    for hemi, region_name in matches:
        output_prefix = os.path.join(atlas_out_dir, sanitize_name(f"{hemi}_{region_name}"))
        saved_files = render_single_region(
            atlas_data=atlas_data,
            hemi=hemi,
            region_name=region_name,
            output_prefix=output_prefix,
        )
        for file_path in saved_files:
            print(f"已保存: {file_path}")

    print(f"输出目录: {atlas_out_dir}")


if __name__ == "__main__":
    main()
