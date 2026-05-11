#!/usr/bin/env python3
"""绘制 ROI-level t 值分布，并绘制 subject-level cluster mean degree 分布。"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

BASE_DIR = "/data/home/tqi/data1/share/after_freesurfer"
FILE_DIR = os.path.join(BASE_DIR, "FILE/test_mean_1.5")
DEMO_FILE = os.path.join(FILE_DIR, "all_data_cqt_mean_1.5.xlsx")
MIND_COMBAT_DIR = os.path.join(FILE_DIR, "MIND_DK318_combat")
OUTPUT_DIR = os.path.join(BASE_DIR, "CODE/Clustering/tvalue_distribution_plots")
os.makedirs(OUTPUT_DIR, exist_ok=True)

CONFIGS = {
    3: {
        "dir": "Clustering_3",
        "clusters": [1, 2, 3],
        "names": {
            1: "C1: Adult Emergence",
            2: "C2: Child Predominant",
            3: "C3: Consistently Weak",
        },
    },
    4: {
        "dir": "Clustering_4",
        "clusters": [1, 2, 3, 4],
        "names": {
            1: "C1: Adult Emergence",
            2: "C2: Persistent Strong",
            3: "C3: Reversal",
            4: "C4: Child Predominant",
        },
    },
    5: {
        "dir": "Clustering_5",
        "clusters": [1, 2, 3, 4, 5],
        "names": {1: "C1", 2: "C2", 3: "C3", 4: "C4", 5: "C5"},
    },
}

GROUP_ORDER = ["Child TD", "Child DD", "Adult TD", "Adult DD"]
GROUP_COLORS = {
    "Child TD": "#4C9BE8",
    "Child DD": "#F28E2B",
    "Adult TD": "#59A14F",
    "Adult DD": "#E15759",
}

SIGNIFICANCE_COMPARISONS = [
    ("Child TD", "Child DD"),
    ("Adult TD", "Adult DD"),
]


def p_to_stars(p_value):
    if pd.isna(p_value):
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""


def add_significance_annotations(ax, arrays_by_group, positions_by_group, y_padding_fraction=0.05):
    finite_values = np.concatenate(
        [arr[np.isfinite(arr)] for arr in arrays_by_group.values() if np.isfinite(arr).any()]
    )
    if len(finite_values) == 0:
        return

    y_min = np.nanmin(finite_values)
    y_max = np.nanmax(finite_values)
    y_range = y_max - y_min
    if y_range == 0:
        y_range = max(abs(y_max), 1.0) * 0.1

    annotation_idx = 0
    for group_a, group_b in SIGNIFICANCE_COMPARISONS:
        arr_a = arrays_by_group[group_a]
        arr_b = arrays_by_group[group_b]
        arr_a = arr_a[np.isfinite(arr_a)]
        arr_b = arr_b[np.isfinite(arr_b)]
        if len(arr_a) < 2 or len(arr_b) < 2:
            continue

        _, p_value = stats.ttest_ind(arr_a, arr_b, equal_var=False, nan_policy="omit")
        stars = p_to_stars(p_value)
        if not stars:
            continue

        x1 = positions_by_group[group_a]
        x2 = positions_by_group[group_b]
        y = y_max + y_range * (y_padding_fraction + annotation_idx * 0.14)
        h = y_range * 0.04
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", lw=1.2, clip_on=False)
        ax.text(
            (x1 + x2) / 2,
            y + h,
            stars,
            ha="center",
            va="bottom",
            color="black",
            fontsize=13,
            fontweight="bold",
            clip_on=False,
        )
        annotation_idx += 1

    if annotation_idx:
        ax.set_ylim(top=y_max + y_range * (y_padding_fraction + (annotation_idx - 1) * 0.14 + 0.16))


def format_id(value):
    if pd.isna(value):
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def cluster_path(k):
    cfg = CONFIGS[k]
    return os.path.join(FILE_DIR, cfg["dir"], "clustering", f"degree_Tmap_cluster_K{k}.csv")


def load_cluster_data(k):
    return pd.read_csv(cluster_path(k))


def plot_bars(k):
    cfg = CONFIGS[k]
    df = load_cluster_data(k)
    print(f"\nK={k}: {len(df)} ROIs")

    n = len(cfg["clusters"])
    fig, axes = plt.subplots(n, 2, figsize=(14, 4 * n))
    if n == 1:
        axes = axes.reshape(1, -1)
    fig.suptitle(f"T-value Distribution K={k}", fontsize=16, fontweight="bold")

    for i, cid in enumerate(cfg["clusters"]):
        cd = df[df["cluster"] == cid]
        mc = cd["t_child_TD_minus_DD"].mean()
        ma = cd["t_adult_TD_minus_DD"].mean()

        ax = axes[i, 0]
        cd_sort = cd.sort_values("t_child_TD_minus_DD", ascending=False)
        ax.bar(
            range(len(cd)),
            cd_sort["t_child_TD_minus_DD"],
            alpha=0.7,
            color="skyblue",
            edgecolor="black",
            lw=0.5,
        )
        ax.axhline(0, color="black", lw=1.5, alpha=0.5)
        ax.axhline(mc, color="red", ls="--", lw=2, label=f"Mean={mc:.2f}")
        ax.set_xlabel("ROI Index (sorted)")
        ax.set_ylabel("t-value (Child: TD-DD)")
        ax.set_title(f'{cfg["names"][cid]} - Child (n={len(cd)})', fontweight="bold")
        ax.legend()
        ax.grid(axis="y", alpha=0.3)

        ax = axes[i, 1]
        cd_sort = cd.sort_values("t_adult_TD_minus_DD", ascending=False)
        ax.bar(
            range(len(cd)),
            cd_sort["t_adult_TD_minus_DD"],
            alpha=0.7,
            color="salmon",
            edgecolor="black",
            lw=0.5,
        )
        ax.axhline(0, color="black", lw=1.5, alpha=0.5)
        ax.axhline(ma, color="red", ls="--", lw=2, label=f"Mean={ma:.2f}")
        ax.set_xlabel("ROI Index (sorted)")
        ax.set_ylabel("t-value (Adult: TD-DD)")
        ax.set_title(f'{cfg["names"][cid]} - Adult (n={len(cd)})', fontweight="bold")
        ax.legend()
        ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, f"tvalue_bars_K{k}.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved: {out}")


def plot_hist(k):
    cfg = CONFIGS[k]
    df = load_cluster_data(k)

    n = len(cfg["clusters"])
    fig, axes = plt.subplots(n, 2, figsize=(14, 4 * n))
    if n == 1:
        axes = axes.reshape(1, -1)
    fig.suptitle(f"T-value Histogram K={k}", fontsize=16, fontweight="bold")

    for i, cid in enumerate(cfg["clusters"]):
        cd = df[df["cluster"] == cid]
        mc = cd["t_child_TD_minus_DD"].mean()
        ma = cd["t_adult_TD_minus_DD"].mean()
        sc = cd["t_child_TD_minus_DD"].std()
        sa = cd["t_adult_TD_minus_DD"].std()

        ax = axes[i, 0]
        ax.hist(cd["t_child_TD_minus_DD"], bins=20, alpha=0.7, color="skyblue", edgecolor="black")
        ax.axvline(0, color="black", lw=2, alpha=0.5)
        ax.axvline(mc, color="red", ls="--", lw=2, label=f"Mean={mc:.2f}")
        ax.set_xlabel("t-value (Child: TD-DD)")
        ax.set_ylabel("Frequency")
        ax.set_title(f'{cfg["names"][cid]} - Child (n={len(cd)}, μ={mc:.2f}, σ={sc:.2f})', fontweight="bold")
        ax.legend()
        ax.grid(axis="y", alpha=0.3)

        ax = axes[i, 1]
        ax.hist(cd["t_adult_TD_minus_DD"], bins=20, alpha=0.7, color="salmon", edgecolor="black")
        ax.axvline(0, color="black", lw=2, alpha=0.5)
        ax.axvline(ma, color="red", ls="--", lw=2, label=f"Mean={ma:.2f}")
        ax.set_xlabel("t-value (Adult: TD-DD)")
        ax.set_ylabel("Frequency")
        ax.set_title(f'{cfg["names"][cid]} - Adult (n={len(cd)}, μ={ma:.2f}, σ={sa:.2f})', fontweight="bold")
        ax.legend()
        ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, f"tvalue_hist_K{k}.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved: {out}")


def load_subject_metadata():
    df = pd.read_excel(DEMO_FILE, sheet_name="Sheet1")
    df = df.copy()
    df["original_project"] = df["original-project"].apply(format_id)
    df["id_old_str"] = df["id_old"].apply(format_id)
    df["subj_prefix"] = df["original_project"] + "_" + df["id_old_str"]
    df["file_base"] = df["subj_prefix"] + "_MIND_DK318_combat"
    df["degree_file"] = df["file_base"].apply(lambda x: os.path.join(MIND_COMBAT_DIR, f"{x}_degree.csv"))
    df["Diagnosis"] = np.where(df["group_d_or_c"] == 0, "TD", "DD")
    df["AgeGroup"] = np.where(df["group_age"] == 1, "Adult", "Child")
    df["Sex"] = np.where(df["sex"] == 1, "Male", "Female")
    df["group_label"] = df["AgeGroup"] + " " + df["Diagnosis"]
    df["has_file"] = df["degree_file"].apply(os.path.exists)

    keep = (
        df["original_project"].ne("")
        & df["id_old_str"].ne("")
        & df["has_file"]
        & df["group_label"].isin(GROUP_ORDER)
    )
    df = df.loc[keep].reset_index(drop=True)
    print(f"Loaded metadata: {len(df)} subjects with degree files")
    print(df["group_label"].value_counts().reindex(GROUP_ORDER, fill_value=0).to_string())
    return df


def read_degree_series(path):
    d = pd.read_csv(path)
    if "ROI" not in d.columns or "degree" not in d.columns:
        raise ValueError(f"Degree file must contain ROI and degree columns: {path}")
    return pd.Series(pd.to_numeric(d["degree"], errors="coerce").values, index=d["ROI"].astype(str).values)


def compute_subject_cluster_means(k, metadata):
    cfg = CONFIGS[k]
    cluster_df = load_cluster_data(k)
    cluster_rois = {
        cid: cluster_df.loc[cluster_df["cluster"] == cid, "feature"].astype(str).tolist()
        for cid in cfg["clusters"]
    }

    rows = []
    for idx, subj in metadata.iterrows():
        if (idx + 1) % 50 == 0 or idx == 0 or idx + 1 == len(metadata):
            print(f"  Subject {idx + 1}/{len(metadata)}")
        degrees = read_degree_series(subj["degree_file"])
        for cid, rois in cluster_rois.items():
            vals = degrees.reindex(rois).dropna()
            rows.append(
                {
                    "subject": subj["subj_prefix"],
                    "file_base": subj["file_base"],
                    "Diagnosis": subj["Diagnosis"],
                    "AgeGroup": subj["AgeGroup"],
                    "Sex": subj["Sex"],
                    "group_label": subj["group_label"],
                    "cluster": cid,
                    "cluster_label": cfg["names"][cid],
                    "mean_degree": vals.mean() if len(vals) else np.nan,
                    "n_rois_available": int(len(vals)),
                    "n_rois_cluster": int(len(rois)),
                }
            )

    out_df = pd.DataFrame(rows)
    out_csv = os.path.join(OUTPUT_DIR, f"subject_cluster_mean_degree_K{k}.csv")
    out_df.to_csv(out_csv, index=False)
    print(f"  Saved: {out_csv}")
    return out_df


def plot_subject_cluster_violin_box(k, subject_df):
    cfg = CONFIGS[k]
    n = len(cfg["clusters"])
    fig, axes = plt.subplots(n, 1, figsize=(12, 4.2 * n), sharex=True)
    if n == 1:
        axes = np.array([axes])
    fig.suptitle(f"Subject-level Cluster Mean Degree K={k}", fontsize=16, fontweight="bold")

    rng = np.random.default_rng(12345)
    positions = np.arange(1, len(GROUP_ORDER) + 1)

    for i, cid in enumerate(cfg["clusters"]):
        ax = axes[i]
        cd = subject_df[subject_df["cluster"] == cid]
        arrays = [cd.loc[cd["group_label"] == group, "mean_degree"].dropna().values for group in GROUP_ORDER]
        arrays_for_violin = [arr if len(arr) else np.array([np.nan]) for arr in arrays]

        valid_positions = [pos for pos, arr in zip(positions, arrays_for_violin) if np.isfinite(arr).any()]
        valid_arrays = [arr[np.isfinite(arr)] for arr in arrays_for_violin if np.isfinite(arr).any()]
        if valid_arrays:
            violins = ax.violinplot(valid_arrays, positions=valid_positions, widths=0.72, showextrema=False)
            for body, pos in zip(violins["bodies"], valid_positions):
                group = GROUP_ORDER[pos - 1]
                body.set_facecolor(GROUP_COLORS[group])
                body.set_edgecolor("black")
                body.set_alpha(0.35)

        box_arrays = [arr if len(arr) else np.array([np.nan]) for arr in arrays]
        ax.boxplot(
            box_arrays,
            positions=positions,
            widths=0.28,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 2},
            boxprops={"facecolor": "white", "edgecolor": "black", "linewidth": 1.2},
            whiskerprops={"color": "black", "linewidth": 1.0},
            capprops={"color": "black", "linewidth": 1.0},
        )

        for pos, group, arr in zip(positions, GROUP_ORDER, arrays):
            if len(arr) == 0:
                continue
            jitter = rng.normal(0, 0.055, size=len(arr))
            ax.scatter(
                np.full(len(arr), pos) + jitter,
                arr,
                s=18,
                alpha=0.65,
                color=GROUP_COLORS[group],
                edgecolors="white",
                linewidths=0.35,
                zorder=3,
            )

        group_ns = [len(arr) for arr in arrays]
        label_with_n = [f"{group}\n(n={n_subj})" for group, n_subj in zip(GROUP_ORDER, group_ns)]
        arrays_by_group = {group: arr for group, arr in zip(GROUP_ORDER, arrays)}
        positions_by_group = {group: pos for group, pos in zip(GROUP_ORDER, positions)}
        add_significance_annotations(ax, arrays_by_group, positions_by_group)
        ax.set_xticks(positions)
        ax.set_xticklabels(label_with_n)
        ax.set_ylabel("Mean MIND degree")
        ax.set_title(f'{cfg["names"][cid]} - cluster mean degree across ROIs', fontweight="bold")
        ax.grid(axis="y", alpha=0.3)

    axes[-1].set_xlabel("AgeGroup + Diagnosis")
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, f"subject_cluster_mean_degree_violin_box_K{k}.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved: {out}")


def plot_subject_level(k, metadata):
    print(f"\nComputing subject-level cluster mean degree for K={k}")
    subject_df = compute_subject_cluster_means(k, metadata)
    plot_subject_cluster_violin_box(k, subject_df)


def main():
    print("=" * 60)
    print("T-value and Subject-level Degree Distribution Visualization")
    print("=" * 60)
    print(f"Output: {OUTPUT_DIR}\n")

    metadata = load_subject_metadata()

    for k in [3, 4, 5]:
        print(f"\n{'=' * 60}")
        print(f"Processing K={k}")
        print(f"{'=' * 60}")
        plot_bars(k)
        plot_hist(k)
        plot_subject_level(k, metadata)

    print(f"\n{'=' * 60}")
    print("All visualizations completed!")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
