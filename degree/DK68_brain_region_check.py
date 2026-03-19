import os
import argparse
import numpy as np
import pandas as pd
import mne
from nibabel.freesurfer.io import read_annot, read_geometry
from nilearn import plotting
import matplotlib.pyplot as plt

BASE_DIR = r"/data/home/tqi/data1/share/after_freesurfer"
DEFAULT_OUT_DIR = os.path.join(BASE_DIR, "FILE/DK68-degree")
SELECTED_ROIS = ["lh_bankssts", "rh_bankssts"]

def load_fsaverage_and_aparc():
    subjects_dir = os.path.join(BASE_DIR, "FILE")
    fs_dir = mne.datasets.fetch_fsaverage(subjects_dir=subjects_dir, verbose=True)
    lh_surf_path = os.path.join(fs_dir, "surf", "lh.inflated")
    rh_surf_path = os.path.join(fs_dir, "surf", "rh.inflated")
    lh_coords, lh_faces = read_geometry(lh_surf_path)
    rh_coords, rh_faces = read_geometry(rh_surf_path)
    lh_annot = os.path.join(fs_dir, "label", "lh.aparc.annot")
    rh_annot = os.path.join(fs_dir, "label", "rh.aparc.annot")
    lh_labels, lh_ctab, lh_names = read_annot(lh_annot)
    rh_labels, rh_ctab, rh_names = read_annot(rh_annot)
    if lh_coords.shape[0] != len(lh_labels) or rh_coords.shape[0] != len(rh_labels):
        raise RuntimeError("Surface vertices do not match annot labels. Ensure fsaverage is used.")
    lh_mesh = (lh_coords, lh_faces)
    rh_mesh = (rh_coords, rh_faces)
    return (lh_mesh, rh_mesh, lh_labels, lh_names, rh_labels, rh_names)

def strip_hemi_prefix(s: str) -> str:
    return str(s).replace("lh_", "").replace("rh_", "")

def make_vertex_data(labels, names, roi_map):
    name_list = [n.decode("utf-8") for n in names]
    v = np.full(labels.shape[0], np.nan, dtype=float)
    hit = 0
    for lab_id, roi in enumerate(name_list):
        if roi in roi_map:
            v[labels == lab_id] = float(roi_map[roi])
            hit += 1
    return v, hit

def load_rois_from_csv(path):
    df = pd.read_csv(path)
    if "feature" in df.columns:
        return df["feature"].astype(str).tolist()
    return df.iloc[:, 0].astype(str).tolist()

def plot_selected_rois(selected_rois, out_dir=DEFAULT_OUT_DIR, cmap_name="Reds"):
    os.makedirs(out_dir, exist_ok=True)
    (lh_mesh, rh_mesh, lh_labels, lh_names, rh_labels, rh_names) = load_fsaverage_and_aparc()
    selected_norm = [strip_hemi_prefix(s) for s in selected_rois]
    roi_label = ", ".join([str(s) for s in selected_rois])
    safe_tag = "".join([c if (c.isalnum() or c in ("_", "-")) else "_" for c in roi_label])
    if len(safe_tag) == 0:
        safe_tag = "selected"
    roi_map = {n: 1.0 for n in selected_norm}
    lh_vtx, lh_hit = make_vertex_data(lh_labels, lh_names, roi_map)
    rh_vtx, rh_hit = make_vertex_data(rh_labels, rh_names, roi_map)
    if lh_hit + rh_hit == 0:
        raise ValueError("No ROI matched annot names. Check ROI labels.")
    plots = [
        ("lh", "lateral", lh_mesh, "left"),
        ("lh", "medial",  lh_mesh, "left"),
        ("rh", "lateral", rh_mesh, "right"),
        ("rh", "medial",  rh_mesh, "right"),
    ]
    vmin, vmax = 0.5, 1.0
    for hemi_tag, view, mesh, hemi_name in plots:
        out_file = os.path.join(out_dir, f"{safe_tag}_{hemi_tag}_{view}.png")
        vtx = lh_vtx if hemi_tag == "lh" else rh_vtx
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection="3d")
        display = plotting.plot_surf_stat_map(
            surf_mesh=mesh,
            stat_map=vtx,
            hemi=hemi_name,
            view=view,
            colorbar=True,
            symmetric_cbar=False,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap_name,
            title=f"DK68 {roi_label} | {hemi_tag.upper()} {view}",
            axes=ax,
            darkness=None
        )
        curr_labels = lh_labels if hemi_tag == "lh" else rh_labels
        try:
            display.add_contours(curr_labels, colors="k", linewidths=0.5)
        except AttributeError:
            unique_labels = np.unique(curr_labels)
            plotting.plot_surf_contours(
                surf_mesh=mesh,
                roi_map=curr_labels,
                hemi=hemi_name,
                view=view,
                axes=ax,
                levels=unique_labels,
                colors=["black"] * len(unique_labels),
                linewidths=0.2
            )
        fig.savefig(out_file, bbox_inches="tight", pad_inches=0.2)
        plt.close(fig)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rois", type=str, default="", help="Comma-separated ROI names")
    parser.add_argument("--roi-file", type=str, default="", help="CSV file with ROI names (feature column or first column)")
    parser.add_argument("--out-dir", type=str, default=DEFAULT_OUT_DIR)
    parser.add_argument("--cmap", type=str, default="Reds")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    roi_list = []
    if args.rois:
        roi_list.extend([s.strip() for s in args.rois.split(",") if s.strip()])
    if args.roi_file:
        roi_list.extend(load_rois_from_csv(args.roi_file))
    if len(roi_list) == 0 and len(SELECTED_ROIS) > 0:
        roi_list = list(SELECTED_ROIS)
    roi_list = list(dict.fromkeys(roi_list))
    if len(roi_list) == 0:
        raise ValueError("No ROI provided. Use --rois or --roi-file, or set SELECTED_ROIS in script.")
    plot_selected_rois(roi_list, out_dir=args.out_dir, cmap_name=args.cmap)