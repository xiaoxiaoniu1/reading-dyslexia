import os
import sys
import subprocess
import importlib
import argparse
from typing import Optional, List, Dict, Tuple


def _ensure_module(module_name: str, pip_name: Optional[str] = None):
    try:
        return importlib.import_module(module_name)
    except ModuleNotFoundError:
        pkg = pip_name or module_name
        print(f"Missing Python package '{pkg}', installing it now...")
        cmd = [sys.executable, "-m", "pip", "install", pkg]
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError:
            subprocess.check_call(cmd + ["--user"])
        return importlib.import_module(module_name)


_ensure_module("numpy")
_ensure_module("pandas")
_ensure_module("mne")
_ensure_module("nibabel")
_ensure_module("nilearn")
_ensure_module("matplotlib")

import numpy as np
import pandas as pd
import mne
from nibabel.freesurfer.io import read_annot, read_geometry
from nilearn import plotting
import matplotlib.pyplot as plt


P_THRESHOLD = 0.05
DEFAULT_PLOT_VMIN = 1.3
DEFAULT_PLOT_VMAX = 3.0
DEFAULT_SIGNED_PLOT_VMAX = DEFAULT_PLOT_VMAX
DEFAULT_SIGNED_PLOT_VMIN = -DEFAULT_SIGNED_PLOT_VMAX

BASE_DIR = r"/data/home/tqi/data1/share/after_freesurfer"
RESULT_BASE_DIR = os.path.join(BASE_DIR, "FILE", "test_mean_1.5")

DEFAULT_ROOTS = {
    "default": os.path.join(RESULT_BASE_DIR, "MIND_DK318_DGLM"),
    "main": os.path.join(RESULT_BASE_DIR, "MIND_DK318_DGLM_main"),
    "sensitivity": os.path.join(RESULT_BASE_DIR, "MIND_DK318_DGLM_sensitivity"),
    "main_sex": os.path.join(RESULT_BASE_DIR, "MIND_DK318_DGLM_main_sex"),
    "sensitivity_sex": os.path.join(RESULT_BASE_DIR, "MIND_DK318_DGLM_sensitivity_sex"),
}

LH_ANNOT = os.path.join(BASE_DIR, r"FILE/DK-318", "lh.DK318.annot")
RH_ANNOT = os.path.join(BASE_DIR, r"FILE/DK-318", "rh.DK318.annot")

DEGREE_CSV = "DGLM_DK318_degree_results.csv"
NETWORK_DEGREE_CSV = "DGLM_DK318_degree_network.csv"

EFFECTS = [
    "Diagnosis",
    "AgeGroup",
    "Sex",
    "Diagnosis_AgeGroup",
    "Diagnosis_Sex",
    "AgeGroup_Sex",
    "Diagnosis_AgeGroup_Sex",
]
MAIN_EFFECTS = {"Diagnosis", "AgeGroup", "Sex"}
DOMAINS = ["mean", "disp"]


def p_adjust_bh(pvals):
    p = pd.to_numeric(pd.Series(pvals), errors="coerce").to_numpy(dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    vals = p[ok]
    n = vals.size
    order = np.argsort(vals)
    ranked = vals[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    tmp = np.empty_like(adj)
    tmp[order] = adj
    out[ok] = tmp
    return out


def safe_name(x: str) -> str:
    import re
    x = re.sub(r"[^A-Za-z0-9_]+", "_", str(x))
    x = re.sub(r"_+", "_", x).strip("_")
    return x or "NA"


def strip_hemi_prefix(s: str) -> str:
    return str(s).replace("lh_", "").replace("rh_", "")


def parse_p_col(col: str) -> Optional[Tuple[str, str, str]]:
    if not col.startswith("p_"):
        return None
    rest = col[2:]
    for domain in DOMAINS:
        prefix = f"{domain}_"
        if not rest.startswith(prefix):
            continue
        tail = rest[len(prefix):]
        for effect in sorted(EFFECTS, key=len, reverse=True):
            if tail == effect:
                return domain, effect, "uncorrected"
            if tail == f"{effect}_FDR":
                return domain, effect, "FDR"
            net_prefix = f"{effect}_networkFDR_"
            if tail.startswith(net_prefix):
                return domain, effect, "networkFDR"
    return None


def augment_missing_plot_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for domain in DOMAINS:
        for effect in EFFECTS:
            raw_p = f"p_{domain}_{effect}"
            fdr_p = f"{raw_p}_FDR"
            if raw_p in df.columns and fdr_p not in df.columns:
                df[fdr_p] = p_adjust_bh(df[raw_p])

            est_col = f"estimate_{domain}_{effect}"
            sign_col = f"sign_{domain}_{effect}"
            if est_col in df.columns and sign_col not in df.columns:
                est = pd.to_numeric(df[est_col], errors="coerce").to_numpy(dtype=float)
                s = np.sign(est)
                s[~np.isfinite(s)] = 1.0
                s[s == 0] = 1.0
                df[sign_col] = s
    return df


def is_main_effect(effect: str) -> bool:
    return effect in MAIN_EFFECTS


def direction_annotation(effect: str) -> str:
    if effect == "Diagnosis":
        return "Red: TD>DD; Blue: DD>TD"
    if effect == "AgeGroup":
        return "Positive: Adult>Child; Negative: Child>Adult"
    if effect == "Sex":
        return "Positive: Male>Female; Negative: Female>Male"
    return ""


def load_fsaverage_and_annot():
    subjects_dir = os.path.join(BASE_DIR, "FILE")
    fs_dir = mne.datasets.fetch_fsaverage(subjects_dir=subjects_dir, verbose=True)
    print("fsaverage dir:", fs_dir)

    lh_surf_path = os.path.join(fs_dir, "surf", "lh.inflated")
    rh_surf_path = os.path.join(fs_dir, "surf", "rh.inflated")

    lh_coords, lh_faces = read_geometry(lh_surf_path)
    rh_coords, rh_faces = read_geometry(rh_surf_path)

    lh_labels, _, lh_names = read_annot(LH_ANNOT)
    rh_labels, _, rh_names = read_annot(RH_ANNOT)

    if lh_coords.shape[0] != len(lh_labels) or rh_coords.shape[0] != len(rh_labels):
        raise RuntimeError("Surface vertex count and DK318 annot vertex count do not match.")

    return (
        (lh_coords, lh_faces),
        (rh_coords, rh_faces),
        lh_labels,
        lh_names,
        rh_labels,
        rh_names,
    )


def make_vertex_data(labels, names, roi_map: Dict[str, float]):
    name_list = [n.decode("utf-8") if isinstance(n, bytes) else str(n) for n in names]
    v = np.full(labels.shape[0], np.nan, dtype=float)
    for lab_id, roi in enumerate(name_list):
        if roi in roi_map:
            v[labels == lab_id] = float(roi_map[roi])
    return v


def resolve_degree_dir(path: str) -> Optional[str]:
    if path is None:
        return None
    path = os.path.normpath(path)
    if os.path.exists(os.path.join(path, DEGREE_CSV)):
        return path
    nested = os.path.join(path, "degree")
    if os.path.exists(os.path.join(nested, DEGREE_CSV)):
        return nested
    return None


def find_network_dirs(degree_dir: str) -> List[str]:
    net_root = os.path.join(degree_dir, "network_FDR")
    if not os.path.isdir(net_root):
        return []
    out = []
    for name in sorted(os.listdir(net_root)):
        p = os.path.join(net_root, name)
        if os.path.isdir(p) and os.path.exists(os.path.join(p, NETWORK_DEGREE_CSV)):
            out.append(p)
    return out


def sex_dirs_for_root(root: str, sex: str) -> List[str]:
    labels = ["male", "female"] if sex == "all" else [sex.lower()]
    out = []
    for s in labels:
        d = os.path.join(root, f"{s}_sex", "degree")
        if os.path.exists(os.path.join(d, DEGREE_CSV)):
            out.append(d)
        d2 = os.path.join(root, f"{s}_sex")
        if os.path.exists(os.path.join(d2, DEGREE_CSV)):
            out.append(d2)
    return out


def discover_degree_dirs(model: str, sex: str = "all", result_dir: Optional[str] = None) -> List[str]:
    if result_dir:
        d = resolve_degree_dir(result_dir)
        if d:
            return [d]
        # result_dir may be sex root containing male_sex/female_sex
        out = sex_dirs_for_root(result_dir, sex)
        if out:
            return out
        raise FileNotFoundError(f"Cannot find {DEGREE_CSV} under: {result_dir}")

    models = ["default", "main", "sensitivity", "main_sex", "sensitivity_sex"] if model == "all" else [model]
    out = []

    for m in models:
        root = DEFAULT_ROOTS[m]
        if m.endswith("_sex"):
            out.extend(sex_dirs_for_root(root, sex))
        else:
            d = resolve_degree_dir(root)
            if d:
                out.append(d)

    # Preserve order, remove duplicates
    seen = set()
    unique = []
    for d in out:
        nd = os.path.normpath(d)
        if nd not in seen:
            unique.append(nd)
            seen.add(nd)
    return unique


def label_for_dir(degree_dir: str, is_network: bool = False) -> str:
    low = degree_dir.lower()
    if "main_sex" in low:
        model = "main sex"
    elif "sensitivity_sex" in low:
        model = "sensitivity sex"
    elif "dglm_main" in low or "mind_dk318_dglm_main" in low:
        model = "main"
    elif "dglm_sensitivity" in low or "mind_dk318_dglm_sensitivity" in low:
        model = "sensitivity"
    else:
        model = "default"

    sex = ""
    if "male_sex" in low:
        sex = " Male"
    elif "female_sex" in low:
        sex = " Female"

    if is_network:
        return f"DK318 {model}{sex} | {os.path.basename(degree_dir)}"
    return f"DK318 {model}{sex}"


def build_plot_columns(df: pd.DataFrame, include_uncorrected: bool, include_fdr: bool,
                       is_network: bool = False) -> List[str]:
    cols = []
    for c in df.columns:
        parsed = parse_p_col(c)
        if parsed is None:
            continue
        _, effect, corr = parsed

        if is_network:
            if corr == "networkFDR" and include_fdr and not c.endswith("_networkFDR_min"):
                cols.append(c)
            continue

        # For main effects, only keep whole-brain FDR results.
        # Do not plot uncorrected or network-specific results.
        if effect in MAIN_EFFECTS:
            if corr == "FDR":
                cols.append(c)
            continue

        if corr == "uncorrected" and include_uncorrected:
            cols.append(c)
        elif corr == "FDR" and include_fdr:
            cols.append(c)
    return cols


def plot_one_p_column(df, value_col, out_root, surfaces, title_prefix,
                      plot_index: Optional[int] = None,
                      total_plots: Optional[int] = None):
    parsed = parse_p_col(value_col)
    if parsed is None:
        return
    domain, effect, correction = parsed

    if value_col not in df.columns:
        return
    if "feature" not in df.columns:
        raise ValueError("Result table lacks 'feature' column.")

    plot_prefix = f"[{plot_index}/{total_plots}] " if plot_index is not None and total_plots is not None else ""
    print(f"{plot_prefix}Preparing brainmap column: {value_col}")

    p = pd.to_numeric(df[value_col], errors="coerce")
    sig = df[np.isfinite(p) & (p < P_THRESHOLD)].copy()
    if sig.empty:
        print(f"{plot_prefix}{value_col}: no significant ROI at p<{P_THRESHOLD}; skip.")
        return

    print(f"{plot_prefix}{value_col}: {len(sig)} significant ROI rows before hemisphere mapping.")

    sig["roi_norm"] = sig["feature"].astype(str).map(strip_hemi_prefix)
    p_sig = pd.to_numeric(sig[value_col], errors="coerce").to_numpy(dtype=float)
    values = -np.log10(np.clip(p_sig, np.finfo(float).tiny, None))

    signed = False
    cmap_name = "viridis"
    symmetric_cbar = False
    vmin = DEFAULT_PLOT_VMIN
    vmax = DEFAULT_PLOT_VMAX

    if is_main_effect(effect):
        sign_col = f"sign_{domain}_{effect}"
        est_col = f"estimate_{domain}_{effect}"
        if sign_col in sig.columns:
            signs = pd.to_numeric(sig[sign_col], errors="coerce").fillna(1.0).to_numpy(dtype=float)
        elif est_col in sig.columns:
            est = pd.to_numeric(sig[est_col], errors="coerce").to_numpy(dtype=float)
            signs = np.sign(est)
            signs[~np.isfinite(signs)] = 1.0
            signs[signs == 0] = 1.0
        else:
            signs = np.ones(sig.shape[0], dtype=float)

        signs[~np.isfinite(signs)] = 1.0
        signs[signs == 0] = 1.0
        if effect == "Diagnosis":
            signs = -signs
        values = values * signs
        signed = True
        cmap_name = "RdBu_r"
        symmetric_cbar = True
        vmin = DEFAULT_SIGNED_PLOT_VMIN
        vmax = DEFAULT_SIGNED_PLOT_VMAX

    sig["_plot_value"] = values

    lh_map = dict(zip(
        sig.loc[sig["feature"].astype(str).str.startswith("lh_"), "roi_norm"],
        sig.loc[sig["feature"].astype(str).str.startswith("lh_"), "_plot_value"],
    ))
    rh_map = dict(zip(
        sig.loc[sig["feature"].astype(str).str.startswith("rh_"), "roi_norm"],
        sig.loc[sig["feature"].astype(str).str.startswith("rh_"), "_plot_value"],
    ))

    print(
        f"{plot_prefix}{value_col}: ROI mapped to hemispheres | "
        f"lh={len(lh_map)} | rh={len(rh_map)}"
    )

    (lh_mesh, rh_mesh, lh_labels, lh_names, rh_labels, rh_names) = surfaces
    lh_vtx = make_vertex_data(lh_labels, lh_names, lh_map)
    rh_vtx = make_vertex_data(rh_labels, rh_names, rh_map)

    lh_n = int(np.isfinite(lh_vtx).sum())
    rh_n = int(np.isfinite(rh_vtx).sum())
    total_n = lh_n + rh_n
    print(
        f"{plot_prefix}{value_col}: finite surface vertices | "
        f"lh={lh_n} | rh={rh_n} | total={total_n}"
    )

    if total_n == 0:
        print(f"{plot_prefix}{value_col}: all mapped surface values are NaN; skip entire brainmap.")
        return

    col_dir = os.path.join(out_root, safe_name(value_col))
    os.makedirs(col_dir, exist_ok=True)

    export = sig[["feature", value_col]].copy()
    export["minus_log10_p"] = -np.log10(np.clip(pd.to_numeric(export[value_col], errors="coerce"), np.finfo(float).tiny, None))
    if signed:
        export["signed_minus_log10_p"] = sig["_plot_value"]
    export_csv = os.path.join(col_dir, f"significant_regions_{safe_name(value_col)}.csv")
    export.to_csv(export_csv, index=False)
    print(f"{plot_prefix}Saved significant ROI table: {export_csv}")

    title = f"{title_prefix} | {domain} {effect} | {correction}"
    if signed:
        note = direction_annotation(effect)
        if note:
            title = f"{title}\n{note}"

    plots = [
        ("lh", "lateral", lh_mesh, lh_vtx, lh_labels, "left"),
        ("lh", "medial",  lh_mesh, lh_vtx, lh_labels, "left"),
        ("rh", "lateral", rh_mesh, rh_vtx, rh_labels, "right"),
        ("rh", "medial",  rh_mesh, rh_vtx, rh_labels, "right"),
    ]

    for hemi_tag, view, mesh, vtx, labels, hemi_name in plots:
        finite_n = int(np.isfinite(vtx).sum())
        if finite_n == 0:
            print(
                f"{plot_prefix}{value_col}: skip {hemi_tag.upper()} {view} "
                f"because all vertex values are NaN."
            )
            continue

        out_file = os.path.join(col_dir, f"{safe_name(value_col)}_{hemi_tag}_{view}.png")
        print(
            f"{plot_prefix}{value_col}: plotting {hemi_tag.upper()} {view} "
            f"with {finite_n} finite vertices -> {out_file}"
        )

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection="3d")

        kwargs = dict(
            surf_mesh=mesh,
            stat_map=vtx,
            hemi=hemi_name,
            view=view,
            colorbar=True,
            symmetric_cbar=symmetric_cbar,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap_name,
            title=f"{title} | {hemi_tag.upper()} {view}",
            axes=ax,
        )

        try:
            display = plotting.plot_surf_stat_map(**kwargs, darkness=None)
        except TypeError as e:
            if "darkness" not in str(e):
                raise
            display = plotting.plot_surf_stat_map(**kwargs)

        try:
            display.add_contours(labels, colors="k", linewidths=0.5)
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
                linewidths=0.2,
            )

        fig.savefig(out_file, bbox_inches="tight", pad_inches=0.2)
        plt.close(fig)
        print("Saved:", out_file)


def plot_degree_dir(degree_dir: str, surfaces: tuple, is_network: bool,
                    include_uncorrected: bool, include_fdr: bool):
    csv_name = NETWORK_DEGREE_CSV if is_network else DEGREE_CSV
    csv_path = os.path.join(degree_dir, csv_name)
    if not os.path.exists(csv_path):
        print("Skip missing:", csv_path)
        return

    print("\n========================================")
    print("Plotting:", csv_path)
    print("Network:", is_network)
    print("========================================")

    df = pd.read_csv(csv_path)
    print(f"Loaded result table: {csv_path} | rows={len(df)} | cols={len(df.columns)}")
    df = augment_missing_plot_columns(df)

    plot_cols = build_plot_columns(
        df,
        include_uncorrected=include_uncorrected,
        include_fdr=include_fdr,
        is_network=is_network,
    )
    if not plot_cols:
        print("No plottable p columns:", csv_path)
        return

    out_root = os.path.join(degree_dir, "brainmaps")
    os.makedirs(out_root, exist_ok=True)

    title_prefix = label_for_dir(degree_dir, is_network=is_network)
    print(f"Plottable columns in {degree_dir}: {len(plot_cols)}")
    for i, col in enumerate(plot_cols, start=1):
        print(f"Starting plot column {i}/{len(plot_cols)}: {col}")
        plot_one_p_column(
            df,
            col,
            out_root,
            surfaces,
            title_prefix,
            plot_index=i,
            total_plots=len(plot_cols),
        )
        print(f"Finished plot column {i}/{len(plot_cols)}: {col}")


def run_all(model: str, sex: str, result_dir: Optional[str],
            include_network: bool, include_uncorrected: bool, include_fdr: bool):
    degree_dirs = discover_degree_dirs(model=model, sex=sex, result_dir=result_dir)
    if not degree_dirs:
        raise FileNotFoundError("No degree result directories found for selected model.")

    print("Selected degree result directories:")
    for d in degree_dirs:
        print(" -", d)

    os.chdir(BASE_DIR)
    surfaces = load_fsaverage_and_annot()

    for di, d in enumerate(degree_dirs, start=1):
        print(f"\n===== Degree result directory {di}/{len(degree_dirs)} =====")
        print(d)
        plot_degree_dir(d, surfaces, is_network=False,
                        include_uncorrected=include_uncorrected,
                        include_fdr=include_fdr)
        if include_network:
            network_dirs = find_network_dirs(d)
            print(f"Network result directories under {d}: {len(network_dirs)}")
            for ni, nd in enumerate(network_dirs, start=1):
                print(f"\n----- Network directory {ni}/{len(network_dirs)} for degree dir {di}/{len(degree_dirs)} -----")
                print(nd)
                plot_degree_dir(nd, surfaces, is_network=True,
                                include_uncorrected=include_uncorrected,
                                include_fdr=include_fdr)

    print("\n========== DK318 DGLM degree brainmap plotting complete ==========")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model",
        default="default",
        choices=["default", "main", "sensitivity", "main_sex", "sensitivity_sex", "all"],
        help="Which model result to plot.",
    )
    parser.add_argument(
        "--sex",
        default="all",
        choices=["all", "male", "female"],
        help="For sex-stratified model results only.",
    )
    parser.add_argument(
        "--result-dir",
        default=None,
        help="Directly specify a DGLM root, degree directory, or sex-stratified root.",
    )
    parser.add_argument("--include-network", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--include-uncorrected", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--include-fdr", action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args()

    run_all(
        model=args.model,
        sex=args.sex,
        result_dir=args.result_dir,
        include_network=args.include_network,
        include_uncorrected=args.include_uncorrected,
        include_fdr=args.include_fdr,
    )
