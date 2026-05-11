#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Visualize DK318 DGLM results computed within clustering-defined ROI sets."""

import argparse
import importlib.util
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE_DIR = Path("/data/home/tqi/data1/share/after_freesurfer")
DEFAULT_RESULT_ROOT = BASE_DIR / "FILE/test_mean_1.5/MIND_DK318_DGLM_by_cluster/degree_Tmap_cluster_K3"
DEFAULT_DEMO_FILE = BASE_DIR / "FILE/test_mean_1.5/all_data_cqt_mean_1.5.xlsx"
DEFAULT_MIND_DIR = BASE_DIR / "FILE/test_mean_1.5/MIND_DK318_combat"
DEGREE_VIS = BASE_DIR / "CODE/degree/degree_analysis_DK318_DGLM.py"
RESULT_CSV = "DGLM_DK318_degree_within_cluster.csv"
P_THRESHOLD = 0.05

EFFECTS = [
    "Diagnosis", "AgeGroup", "Sex", "Diagnosis_AgeGroup", "Diagnosis_Sex",
    "AgeGroup_Sex", "Diagnosis_AgeGroup_Sex",
]
INTERACTION_EFFECTS = {
    "Diagnosis_AgeGroup", "Diagnosis_Sex", "AgeGroup_Sex", "Diagnosis_AgeGroup_Sex",
}
DOMAINS = ["mean", "disp"]
GROUP_ORDER = [
    ("TD", "Child", "Female"), ("DD", "Child", "Female"),
    ("TD", "Adult", "Female"), ("DD", "Adult", "Female"),
    ("TD", "Child", "Male"), ("DD", "Child", "Male"),
    ("TD", "Adult", "Male"), ("DD", "Adult", "Male"),
]
GROUP_COLORS = {
    "Child-TD-Female": "#4C78A8", "Child-DD-Female": "#F58518",
    "Adult-TD-Female": "#54A24B", "Adult-DD-Female": "#E45756",
    "Child-TD-Male": "#72B7B2", "Child-DD-Male": "#FF9DA6",
    "Adult-TD-Male": "#B279A2", "Adult-DD-Male": "#9D755D",
}
FACTOR_LEVELS = {
    "Diagnosis": ["TD", "DD"],
    "AgeGroup": ["Child", "Adult"],
    "Sex": ["Female", "Male"],
}
INTERACTION_FACTOR_ORDER = {
    "Diagnosis_AgeGroup": ["AgeGroup", "Diagnosis"],
    "Diagnosis_Sex": ["Sex", "Diagnosis"],
    "AgeGroup_Sex": ["Sex", "AgeGroup"],
    "Diagnosis_AgeGroup_Sex": ["AgeGroup", "Diagnosis", "Sex"],
}
INTERACTION_GROUP_COLORS = {
    "Child-TD": "#4C78A8", "Child-DD": "#F58518",
    "Adult-TD": "#54A24B", "Adult-DD": "#E45756",
    "Female-TD": "#4C78A8", "Female-DD": "#F58518",
    "Male-TD": "#72B7B2", "Male-DD": "#FF9DA6",
    "Female-Child": "#4C78A8", "Female-Adult": "#54A24B",
    "Male-Child": "#72B7B2", "Male-Adult": "#B279A2",
}


def load_degree_brainmap_module():
    spec = importlib.util.spec_from_file_location("degree_dglm_vis", DEGREE_VIS)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["degree_dglm_vis"] = mod
    spec.loader.exec_module(mod)
    return mod


def safe_name(x):
    x = re.sub(r"[^A-Za-z0-9_]+", "_", str(x))
    return re.sub(r"_+", "_", x).strip("_") or "NA"


def p_adjust_bh(pvals):
    p = pd.to_numeric(pd.Series(pvals), errors="coerce").to_numpy(float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order]
    adj = ranked * len(vals) / np.arange(1, len(vals) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    tmp = np.empty_like(adj)
    tmp[order] = np.clip(adj, 0, 1)
    out[ok] = tmp
    return out


def normalize_for_brainmap(df):
    df = df.copy()
    for domain in DOMAINS:
        for effect in EFFECTS:
            raw = f"p_{domain}_{effect}"
            fdr_cluster = f"{raw}_FDR_cluster"
            fdr = f"{raw}_FDR"
            if fdr_cluster in df.columns:
                df[fdr] = pd.to_numeric(df[fdr_cluster], errors="coerce")
            elif raw in df.columns and fdr not in df.columns:
                df[fdr] = p_adjust_bh(df[raw])

            est = f"estimate_{domain}_{effect}"
            sign = f"sign_{domain}_{effect}"
            if est in df.columns and sign not in df.columns:
                s = np.sign(pd.to_numeric(df[est], errors="coerce").to_numpy(float))
                s[~np.isfinite(s)] = 1
                s[s == 0] = 1
                df[sign] = s
    return df


def cluster_dirs(result_root):
    return sorted([p for p in Path(result_root).glob("cluster_*") if (p / RESULT_CSV).exists()])


def plot_cols(df, correction="FDR_cluster", domains=None, effects=None):
    domains = domains or DOMAINS
    effects = effects or EFFECTS
    cols = []
    for d in domains:
        for e in effects:
            if correction == "FDR_cluster":
                c = f"p_{d}_{e}_FDR_cluster"
            elif correction == "raw":
                c = f"p_{d}_{e}"
            else:
                c = f"p_{d}_{e}_FDR"
            if c in df.columns:
                cols.append(c)
    return cols


def sig_rows(df, p_col, threshold):
    p = pd.to_numeric(df[p_col], errors="coerce")
    return df[np.isfinite(p) & (p < threshold)].copy()


def raw_group_table(row):
    rows = []
    for diagnosis, age, sex in GROUP_ORDER:
        raw_key = f"{diagnosis}_{age}_{sex}"
        label = f"{age}-{diagnosis}-{sex}"
        n = row.get(f"raw_n_{raw_key}", np.nan)
        mean = row.get(f"raw_mean_{raw_key}", np.nan)
        sd = row.get(f"raw_sd_{raw_key}", np.nan)
        n = pd.to_numeric(n, errors="coerce")
        mean = pd.to_numeric(mean, errors="coerce")
        sd = pd.to_numeric(sd, errors="coerce")
        se = sd / np.sqrt(n) if np.isfinite(sd) and np.isfinite(n) and n > 0 else np.nan
        rows.append({"Group": label, "Diagnosis": diagnosis, "AgeGroup": age, "Sex": sex, "n": n, "mean": mean, "se": se})
    return pd.DataFrame(rows)


def plot_roi_raw_bars(row, p_col, out_file, title):
    tab = raw_group_table(row)
    x = np.arange(len(tab))
    colors = [GROUP_COLORS.get(g, "#777777") for g in tab["Group"]]
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(x, tab["mean"], yerr=tab["se"], color=colors, edgecolor="black", linewidth=0.6, capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels(tab["Group"], rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("Raw degree mean ± SE")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.grid(axis="y", linestyle="--", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_top_roi_bars(df, p_col, out_dir, threshold, top_n):
    sig = sig_rows(df, p_col, threshold)
    if sig.empty:
        return 0
    sig["_p"] = pd.to_numeric(sig[p_col], errors="coerce")
    sig = sig.sort_values("_p").head(top_n)
    bar_dir = Path(out_dir) / safe_name(p_col) / "raw_group_bars"
    bar_dir.mkdir(parents=True, exist_ok=True)
    for _, row in sig.iterrows():
        roi = str(row["feature"])
        pval = row["_p"]
        est_col = p_col.replace("p_", "estimate_").replace("_FDR_cluster", "")
        est = row.get(est_col, np.nan)
        title = f"{roi}\n{p_col}={pval:.3g}; estimate={pd.to_numeric(est, errors='coerce'):.3g}"
        plot_roi_raw_bars(row, p_col, bar_dir / f"{safe_name(roi)}.png", title)
    sig.drop(columns=["_p"]).to_csv(bar_dir / f"top{top_n}_significant_roi_table.csv", index=False)
    return len(sig)


def parse_p_col_effect(p_col):
    for domain in DOMAINS:
        prefix = f"p_{domain}_"
        if not p_col.startswith(prefix):
            continue
        effect = p_col[len(prefix):]
        for suffix in ("_FDR_cluster", "_FDR"):
            if effect.endswith(suffix):
                effect = effect[:-len(suffix)]
        return domain, effect
    return None, None


def is_interaction_p_col(p_col):
    _, effect = parse_p_col_effect(p_col)
    return effect in INTERACTION_EFFECTS


def load_subject_degree_long(demo_file=DEFAULT_DEMO_FILE, mind_dir=DEFAULT_MIND_DIR):
    demo_file = Path(demo_file)
    mind_dir = Path(mind_dir)
    if not demo_file.exists():
        raise FileNotFoundError(f"Cannot find subject metadata file: {demo_file}")
    if not mind_dir.exists():
        raise FileNotFoundError(f"Cannot find MIND/degree directory: {mind_dir}")

    subj = pd.read_excel(demo_file, sheet_name="Sheet1")
    required = ["original-project", "id_old", "group_d_or_c", "group_age", "sex"]
    missing = [c for c in required if c not in subj.columns]
    if missing:
        raise ValueError(f"Subject metadata is missing columns: {missing}")

    subj = subj.copy()
    subj["original_project"] = subj["original-project"].astype(str)
    subj["id_old"] = subj["id_old"].astype(str)
    subj["subj_prefix"] = subj["original_project"] + "_" + subj["id_old"]
    subj["file_base"] = subj["subj_prefix"] + "_MIND_DK318_combat"
    subj["degree_file"] = subj["file_base"].map(lambda x: mind_dir / f"{x}_degree.csv")
    subj["Diagnosis"] = np.where(pd.to_numeric(subj["group_d_or_c"], errors="coerce") == 0, "TD", "DD")
    subj["AgeGroup"] = np.where(pd.to_numeric(subj["group_age"], errors="coerce") == 1, "Adult", "Child")
    subj["Sex"] = np.where(pd.to_numeric(subj["sex"], errors="coerce") == 1, "Male", "Female")
    subj = subj[
        subj["degree_file"].map(lambda p: Path(p).exists())
        & subj["original_project"].notna()
        & subj["id_old"].notna()
        & (subj["original_project"] != "")
        & (subj["id_old"] != "")
    ].copy()
    if subj.empty:
        raise FileNotFoundError(f"No subject degree files found under: {mind_dir}")
    return subj


def read_subject_roi_degree(subj, roi):
    rows = []
    for _, s in subj.iterrows():
        deg = pd.read_csv(s["degree_file"], usecols=["ROI", "degree"])
        hit = deg.loc[deg["ROI"].astype(str) == roi, "degree"]
        if hit.empty:
            continue
        value = pd.to_numeric(hit.iloc[0], errors="coerce")
        if not np.isfinite(value):
            continue
        label = f"{s['AgeGroup']}-{s['Diagnosis']}-{s['Sex']}"
        rows.append({
            "Subject": s["subj_prefix"],
            "Diagnosis": s["Diagnosis"],
            "AgeGroup": s["AgeGroup"],
            "Sex": s["Sex"],
            "Group": label,
            "degree": float(value),
        })
    return pd.DataFrame(rows)


def interaction_group_spec(effect):
    factors = INTERACTION_FACTOR_ORDER.get(effect)
    if not factors:
        return None, None
    labels = [""]
    for factor in factors:
        labels = [f"{prefix}-{level}" if prefix else level for prefix in labels for level in FACTOR_LEVELS[factor]]
    return factors, labels


def add_interaction_group(values, effect):
    factors, labels = interaction_group_spec(effect)
    if factors is None:
        return values, []
    values = values.copy()
    values["InteractionGroup"] = values[factors].astype(str).agg("-".join, axis=1)
    return values, labels


def plot_roi_violin(values, row, p_col, out_file, title):
    if values.empty:
        return False
    _, effect = parse_p_col_effect(p_col)
    values, labels = add_interaction_group(values, effect)
    if not labels:
        labels = [f"{age}-{diagnosis}-{sex}" for diagnosis, age, sex in GROUP_ORDER]
        values = values.copy()
        values["InteractionGroup"] = values["Group"]
    data = [values.loc[values["InteractionGroup"] == lab, "degree"].to_numpy(float) for lab in labels]
    if not any(len(x) > 0 for x in data):
        return False

    x = np.arange(1, len(labels) + 1)
    colors = [INTERACTION_GROUP_COLORS.get(g, GROUP_COLORS.get(g, "#777777")) for g in labels]
    fig_w = 8 if len(labels) <= 4 else 12.5
    fig, ax = plt.subplots(figsize=(fig_w, 5.5))
    nonempty = [i for i, arr in enumerate(data) if len(arr) > 0]
    parts = ax.violinplot(
        [data[i] for i in nonempty],
        positions=[x[i] for i in nonempty],
        widths=0.75,
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )
    for body, i in zip(parts["bodies"], nonempty):
        body.set_facecolor(colors[i])
        body.set_edgecolor("black")
        body.set_alpha(0.45)
        body.set_linewidth(0.8)
    if "cmedians" in parts:
        parts["cmedians"].set_color("black")
        parts["cmedians"].set_linewidth(1.2)

    rng = np.random.default_rng(20260511)
    for i, arr in enumerate(data):
        if len(arr) == 0:
            continue
        jitter = rng.normal(0, 0.055, size=len(arr))
        ax.scatter(
            np.full(len(arr), x[i]) + jitter,
            arr,
            s=18,
            color=colors[i],
            edgecolor="white",
            linewidth=0.35,
            alpha=0.75,
            zorder=3,
        )
        ax.scatter(
            x[i],
            np.nanmean(arr),
            s=80,
            marker="D",
            color="black",
            edgecolor="white",
            linewidth=0.8,
            zorder=4,
        )

    pval = pd.to_numeric(row.get(p_col), errors="coerce")
    est_col = p_col.replace("p_", "estimate_").replace("_FDR_cluster", "").replace("_FDR", "")
    est = pd.to_numeric(row.get(est_col, np.nan), errors="coerce")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("Subject-level degree")
    ax.set_title(f"{title}\n{p_col}={pval:.3g}; estimate={est:.3g}", fontsize=11, fontweight="bold")
    ax.grid(axis="y", linestyle="--", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_interaction_roi_violins(df, p_col, out_dir, threshold, subj):
    if not is_interaction_p_col(p_col):
        return 0
    sig = sig_rows(df, p_col, threshold)
    if sig.empty:
        return 0
    sig = sig.copy()
    sig["_p"] = pd.to_numeric(sig[p_col], errors="coerce")
    sig = sig.sort_values("_p")
    violin_dir = Path(out_dir) / safe_name(p_col) / "subject_violin_plots_by_interaction"
    violin_dir.mkdir(parents=True, exist_ok=True)
    made = 0
    all_values = []
    for _, row in sig.iterrows():
        roi = str(row["feature"])
        values = read_subject_roi_degree(subj, roi)
        if values.empty:
            continue
        _, effect = parse_p_col_effect(p_col)
        export_values, _ = add_interaction_group(values, effect)
        export_values.insert(0, "feature", roi)
        export_values.insert(0, "effect_col", p_col)
        all_values.append(export_values)
        title = f"{roi} | {effect} interaction"
        if plot_roi_violin(values, row, p_col, violin_dir / f"{safe_name(roi)}.png", title):
            made += 1
    if all_values:
        pd.concat(all_values, ignore_index=True).to_csv(violin_dir / "subject_level_degree_values.csv", index=False)
        sig.drop(columns=["_p"]).to_csv(violin_dir / "interaction_significant_roi_table.csv", index=False)
    return made


def plot_cluster_summary(summary, out_file):
    if summary.empty:
        return
    pivot = summary.pivot_table(index="cluster", columns="effect_col", values="n_sig", fill_value=0, aggfunc="sum")
    fig_w = max(10, 0.35 * len(pivot.columns))
    fig_h = max(4, 0.55 * len(pivot.index))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(pivot.values, aspect="auto", cmap="Reds")
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels(pivot.index)
    ax.set_title("Number of significant ROIs within each cluster")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = int(pivot.values[i, j])
            if val > 0:
                ax.text(j, i, str(val), ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, label="n significant ROI")
    fig.tight_layout()
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_cluster_brainmaps(brain, surfaces, df, cluster_name, cluster_dir, out_dir, threshold, cols):
    brain.P_THRESHOLD = threshold
    brain_df = normalize_for_brainmap(df)
    brain_out = Path(out_dir) / "brainmaps"
    brain_out.mkdir(parents=True, exist_ok=True)
    for c in cols:
        brain_col = c.replace("_FDR_cluster", "_FDR")
        title = f"DK318 cluster DGLM | {cluster_name}"
        brain.plot_one_p_column(brain_df, brain_col, str(brain_out), surfaces, title)


def run(result_root, threshold, top_n, correction, domains, effects, skip_brainmaps, skip_violins, demo_file, mind_dir):
    result_root = Path(result_root)
    cdirs = cluster_dirs(result_root)
    if not cdirs:
        raise FileNotFoundError(f"No cluster_* directories with {RESULT_CSV} under {result_root}")

    all_summary = []
    subj = None
    if not skip_violins:
        subj = load_subject_degree_long(demo_file=demo_file, mind_dir=mind_dir)
        print(f"Loaded subject-level degree metadata: n={len(subj)}")

    brain = None
    surfaces = None
    if not skip_brainmaps:
        brain = load_degree_brainmap_module()
        os.chdir(BASE_DIR)
        surfaces = brain.load_fsaverage_and_annot()

    for cdir in cdirs:
        cluster_name = cdir.name
        csv_path = cdir / RESULT_CSV
        out_dir = cdir / "visualization"
        out_dir.mkdir(exist_ok=True)
        df = pd.read_csv(csv_path)
        cols = plot_cols(df, correction=correction, domains=domains, effects=effects)
        print(f"\n{cluster_name}: rows={len(df)}, plot columns={len(cols)}")
        if not cols:
            continue

        for col in cols:
            sig = sig_rows(df, col, threshold)
            all_summary.append({
                "cluster": cluster_name, "effect_col": col, "n_roi_total": len(df),
                "n_sig": len(sig), "min_p": pd.to_numeric(df[col], errors="coerce").min(),
            })
            sig.to_csv(out_dir / f"significant_regions_{safe_name(col)}.csv", index=False)
            n_bar = plot_top_roi_bars(df, col, out_dir, threshold, top_n)
            n_violin = 0 if skip_violins else plot_interaction_roi_violins(df, col, out_dir, threshold, subj)
            suffix = f", interaction violins={n_violin}" if is_interaction_p_col(col) else ""
            print(f"  {col}: n_sig={len(sig)}, raw barplots={n_bar}{suffix}")

        if not skip_brainmaps:
            plot_cluster_brainmaps(brain, surfaces, df, cluster_name, cdir, out_dir, threshold, cols)

    summary = pd.DataFrame(all_summary)
    summary.to_csv(result_root / "cluster_DGLM_visualization_summary.csv", index=False)
    plot_cluster_summary(summary, result_root / "cluster_DGLM_significant_count_heatmap.png")
    print(f"\nDONE. Summary saved to: {result_root / 'cluster_DGLM_visualization_summary.csv'}")


def parse_list_arg(x, default):
    if x in (None, "all", ""):
        return default
    return [v.strip() for v in x.split(",") if v.strip()]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize within-cluster DK318 DGLM degree results.")
    parser.add_argument("--result-root", default=str(DEFAULT_RESULT_ROOT))
    parser.add_argument("--threshold", type=float, default=P_THRESHOLD)
    parser.add_argument("--top-n", type=int, default=10, help="Top significant ROIs per effect for raw group barplots.")
    parser.add_argument("--correction", choices=["FDR_cluster", "raw", "FDR"], default="FDR_cluster")
    parser.add_argument("--domains", default="mean,disp", help="Comma list: mean,disp")
    parser.add_argument("--effects", default="all", help="Comma list or all")
    parser.add_argument("--skip-brainmaps", action="store_true", help="Only export tables/barplots/summary heatmap.")
    parser.add_argument("--skip-violins", action="store_true", help="Do not make subject-level violin plots for significant interaction ROIs.")
    parser.add_argument("--demo-file", default=str(DEFAULT_DEMO_FILE), help="Subject metadata Excel file.")
    parser.add_argument("--mind-dir", default=str(DEFAULT_MIND_DIR), help="Directory containing subject DK318 degree CSV files.")
    args = parser.parse_args()

    run(
        result_root=args.result_root,
        threshold=args.threshold,
        top_n=args.top_n,
        correction=args.correction,
        domains=parse_list_arg(args.domains, DOMAINS),
        effects=parse_list_arg(args.effects, EFFECTS),
        skip_brainmaps=args.skip_brainmaps,
        skip_violins=args.skip_violins,
        demo_file=args.demo_file,
        mind_dir=args.mind_dir,
    )
