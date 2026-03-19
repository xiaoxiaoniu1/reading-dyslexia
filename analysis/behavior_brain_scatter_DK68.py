import os
import re
import glob
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

behavior_file = "/data/home/tqi/data1/share/after_freesurfer/FILE/caofan_behavior_data.xlsx"
sheet_name = 0

parcellation_tag = "DK68"
expected_n_roi = 68

use_labeled_csv = True
mind_dir_csv = "/data/home/tqi/data1/share/after_freesurfer/fs_subjects_all/MIND_DK68_combat"
mind_dir_npy = None

output_dir = "/data/home/tqi/data1/share/after_freesurfer/FILE/behavior_brain_scatter_dk68"

diagnostic_behaviors = [
    "CharacterNaming(150)",
    "ReadingFluency(100)",
    "One-minuteA(150)",
    "One-minuteB(150)"
]

nondiagnostic_behavior = "Initial sound deletion(30)"

def norm_name(s):
    return re.sub(r"[^a-z0-9]+", "", str(s).lower())

def find_col(df, candidates):
    norm_map = {norm_name(c): c for c in df.columns}
    for c in candidates:
        key = norm_name(c)
        if key in norm_map:
            return norm_map[key]
    return None

def slugify(s):
    s = re.sub(r"[^\w\-]+", "_", str(s))
    s = re.sub(r"_+", "_", s).strip("_")
    return s

def map_diagnosis(series):
    vals = series.dropna().unique()
    if len(vals) == 0:
        return series
    if all(str(v).isdigit() for v in vals):
        u = sorted(set(int(v) for v in vals))
        if set(u).issubset({0, 1}):
            return series.map({0: "TD", 1: "DD"})
        if set(u).issubset({1, 2}):
            return series.map({1: "TD", 2: "DD"})
    def _map_one(v):
        if pd.isna(v):
            return np.nan
        v = str(v).strip().lower()
        if "dd" in v or "dys" in v:
            return "DD"
        if "td" in v or "control" in v:
            return "TD"
        return np.nan
    return series.map(_map_one)

def map_agegroup(series):
    vals = series.dropna().unique()
    if len(vals) == 0:
        return series
    if all(str(v).isdigit() for v in vals):
        u = sorted(set(int(v) for v in vals))
        if set(u).issubset({1, 2}):
            return series.map({1: "Adult", 2: "Child"})
    def _map_one(v):
        if pd.isna(v):
            return np.nan
        v = str(v).strip().lower()
        if "adult" in v:
            return "Adult"
        if "child" in v or "kid" in v:
            return "Child"
        return np.nan
    return series.map(_map_one)

def map_sex(series):
    vals = series.dropna().unique()
    if len(vals) == 0:
        return series
    if all(str(v).isdigit() for v in vals):
        u = sorted(set(int(v) for v in vals))
        if set(u).issubset({1, 2}):
            return series.map({1: "Male", 2: "Female"})
        if set(u).issubset({0, 1}):
            return series.map({0: "Female", 1: "Male"})
    def _map_one(v):
        if pd.isna(v):
            return np.nan
        v = str(v).strip().lower()
        if v.startswith("m") or "male" in v:
            return "Male"
        if v.startswith("f") or "female" in v:
            return "Female"
        return np.nan
    return series.map(_map_one)

def resolve_degree_file(row, id_col):
    if id_col:
        sid = str(row[id_col]).strip()
        if sid != "nan" and len(sid) > 0:
            f = os.path.join(mind_dir_csv, f"{sid}_combat.csv")
            if os.path.exists(f):
                return f
            hits = glob.glob(os.path.join(mind_dir_csv, f"*{sid}*_combat.csv"))
            if len(hits) == 1:
                return hits[0]
    return np.nan

def load_matrix(path):
    if use_labeled_csv:
        mat = pd.read_csv(path, index_col=0)
        return mat.values, mat.columns.tolist()
    mat = np.load(path)
    return mat, None

def get_roi_names_from_any_labeled():
    files = glob.glob(os.path.join(mind_dir_csv, f"*_MIND_{parcellation_tag}_combat_labeled.csv"))
    if not files:
        return None
    mat = pd.read_csv(files[0], index_col=0)
    return mat.columns.tolist()

def calc_degree(mat):
    mat_copy = mat.copy()
    np.fill_diagonal(mat_copy, np.nan)
    return np.nanmean(mat_copy, axis=1)

def format_p(p):
    if not np.isfinite(p):
        return "NA"
    return f"{p:.3g}"

def extract_p(pvalues, pattern):
    hits = [p for name, p in pvalues.items() if pattern in name]
    if not hits:
        return np.nan
    return float(np.min(hits))

def zscore_array(arr):
    mean = np.nanmean(arr)
    std = np.nanstd(arr)
    if not np.isfinite(std) or std <= 0:
        return np.full_like(arr, np.nan, dtype=float)
    return (arr - mean) / std

def add_fdr(res_df, alpha=0.05):
    res = res_df.copy()
    res["p_target"] = res["p_behavior"]
    valid = res["p_target"].notna()
    res["p_fdr"] = np.nan
    if valid.sum() > 0:
        res.loc[valid, "p_fdr"] = multipletests(res.loc[valid, "p_target"], method="fdr_bh")[1]
    res["significant_fdr"] = res["p_fdr"].notna() & (res["p_fdr"] < alpha)
    return res

def add_fdr_interaction(res_df, p_col, suffix, alpha=0.05):
    res = res_df.copy()
    fdr_col = f"p_fdr_{suffix}"
    sig_col = f"significant_fdr_{suffix}"
    valid = res[p_col].notna()
    res[fdr_col] = np.nan
    if valid.sum() > 0:
        res.loc[valid, fdr_col] = multipletests(res.loc[valid, p_col], method="fdr_bh")[1]
    res[sig_col] = res[fdr_col].notna() & (res[fdr_col] < alpha)
    return res

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def plot_simple(x, y, title, r, p_r, out_path, extra_text=None):
    plt.figure(figsize=(5, 4))
    plt.scatter(y, x, s=18, alpha=0.7)
    if len(x) > 1:
        b1, b0 = np.polyfit(y, x, 1)
        xs = np.linspace(np.min(y), np.max(y), 100)
        ys = b1 * xs + b0
        plt.plot(xs, ys, color="red", linewidth=1.5)
    plt.title(title)
    plt.xlabel("Behavior (z)")
    plt.ylabel("Brain")
    text = extra_text if extra_text is not None else f"r={r:.3f}, p_r={format_p(p_r)}"
    plt.text(0.02, 0.98, text, transform=plt.gca().transAxes, va="top")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

def plot_interaction(x, y, group, title, r, p_r, out_path, extra_text=None, legend_outside=False):
    plt.figure(figsize=(5.5, 4.5))
    for g in sorted(group.unique()):
        mask = group == g
        plt.scatter(y[mask], x[mask], s=18, alpha=0.7, label=str(g))
        if mask.sum() > 1:
            b1, b0 = np.polyfit(y[mask], x[mask], 1)
            xs = np.linspace(np.min(y[mask]), np.max(y[mask]), 100)
            ys = b1 * xs + b0
            plt.plot(xs, ys, linewidth=1.3)
    plt.title(title)
    plt.xlabel("Behavior (z)")
    plt.ylabel("Brain")
    if legend_outside:
        plt.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    else:
        plt.legend(frameon=False)
    text = extra_text if extra_text is not None else f"r={r:.3f}, p_r={format_p(p_r)}"
    plt.text(0.02, 0.98, text, transform=plt.gca().transAxes, va="top")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

def run_behavior_models(df, brain_mat, brain_names, behavior_col, formula, required_cols, group_col=None, out_base=None, tag=None, p_check_terms=None, alpha=0.05, legend_outside=False):
    results = []
    y_all = pd.to_numeric(df[behavior_col], errors="coerce").values
    y_all = zscore_array(y_all)
    for k, brain_name in enumerate(brain_names):
        x_all = brain_mat[:, k]
        valid = np.isfinite(x_all) & np.isfinite(y_all)
        for col in required_cols:
            valid = valid & df[col].notna().values
        if group_col is not None:
            valid = valid & df[group_col].notna().values
        if valid.sum() < 3:
            results.append({
                "roi": brain_name,
                "n": int(valid.sum()),
                "r": np.nan,
                "p_r": np.nan,
                "p_behavior": np.nan,
                "p_interaction_agegroup": np.nan,
                "p_interaction_diagnosis": np.nan,
                "p_interaction_dx_age": np.nan
            })
            continue
        x = x_all[valid]
        y = y_all[valid]
        r, p_r = stats.pearsonr(x, y)
        data = df.loc[valid, required_cols].copy()
        data["brain"] = x
        data["behavior"] = y
        if group_col is not None:
            data[group_col] = df.loc[valid, group_col].values
        for col in ["Diagnosis", "AgeGroup", "Sex", "DxAgeGroup"]:
            if col in data.columns:
                data[col] = data[col].astype("category")
        fit = smf.ols(formula, data=data).fit()
        p_behavior = float(fit.pvalues.get("behavior", np.nan))
        p_int_age = extract_p(fit.pvalues, "behavior:AgeGroup")
        p_int_dx = extract_p(fit.pvalues, "behavior:Diagnosis")
        p_int_dx_age = extract_p(fit.pvalues, "behavior:Diagnosis:AgeGroup")
        results.append({
            "roi": brain_name,
            "n": int(valid.sum()),
            "r": float(r),
            "p_r": float(p_r),
            "p_behavior": p_behavior,
            "p_interaction_agegroup": p_int_age,
            "p_interaction_diagnosis": p_int_dx,
            "p_interaction_dx_age": p_int_dx_age
        })
        plot_sig = False
        if p_check_terms:
            p_vals = []
            for term in p_check_terms:
                if term == "behavior":
                    p_vals.append(p_behavior)
                elif term == "behavior:AgeGroup":
                    p_vals.append(p_int_age)
                elif term == "behavior:Diagnosis":
                    p_vals.append(p_int_dx)
                elif term == "behavior:Diagnosis:AgeGroup":
                    p_vals.append(p_int_dx_age)
            plot_sig = any(np.isfinite(p) and p < alpha for p in p_vals)
        if out_base and plot_sig:
            out_dir = out_base
            ensure_dir(out_dir)
            title = f"{behavior_col} | {brain_name}"
            text_lines = [f"r={r:.3f}, p_r={format_p(p_r)}", f"p_behavior={format_p(p_behavior)}"]
            if np.isfinite(p_int_age):
                text_lines.append(f"p_behavior:AgeGroup={format_p(p_int_age)}")
            if np.isfinite(p_int_dx):
                text_lines.append(f"p_behavior:Diagnosis={format_p(p_int_dx)}")
            if np.isfinite(p_int_dx_age):
                text_lines.append(f"p_behavior:Diagnosis:AgeGroup={format_p(p_int_dx_age)}")
            extra_text = "\n".join(text_lines)
            if group_col is None:
                out_path = os.path.join(out_dir, f"{slugify(brain_name)}_{tag}.png")
                plot_simple(x, y, title, r, p_r, out_path, extra_text=extra_text)
            else:
                out_path = os.path.join(out_dir, f"{slugify(brain_name)}_{tag}.png")
                plot_interaction(x, y, data[group_col], title, r, p_r, out_path, extra_text=extra_text, legend_outside=legend_outside)
    res_df = pd.DataFrame(results)
    return res_df

def main():
    ensure_dir(output_dir)
    df = pd.read_excel(behavior_file, sheet_name=sheet_name)

    diagnosis_col = find_col(df, ["Diagnosis", "diagnosis", "group_d_or_c", "group", "DD_TD", "dx"])
    agegroup_col = find_col(df, ["AgeGroup", "agegroup", "group_age", "age_group", "GroupAge"])
    sex_col = find_col(df, ["Sex", "sex", "gender", "Gender", "SEX"])
    project_col = None
    id_col = find_col(df, ["id", "ID", "subject_id", "participant_id"]) 

    if diagnosis_col:
        df["Diagnosis"] = map_diagnosis(df[diagnosis_col])
    else:
        df["Diagnosis"] = np.nan
    if agegroup_col:
        df["AgeGroup"] = map_agegroup(df[agegroup_col])
    else:
        df["AgeGroup"] = np.nan
    if sex_col:
        df["Sex"] = map_sex(df[sex_col])
    else:
        df["Sex"] = np.nan

    df["DxAgeGroup"] = np.where(
        df["Diagnosis"].notna() & df["AgeGroup"].notna(),
        df["Diagnosis"] + "_" + df["AgeGroup"],
        np.nan
    )

    df["mind_file"] = df.apply(lambda r: resolve_degree_file(r, id_col), axis=1)
    df = df[df["mind_file"].notna()].copy()
    df["mind_file"] = df["mind_file"].astype(str)
    df = df[df["mind_file"].apply(os.path.exists)].copy()

    if df.empty:
        raise RuntimeError("No valid subjects matched to DK68 degree CSVs by id. Please check the 'id' column and folder path.")

    first_path = df["mind_file"].iloc[0]
    mat0, roi_names = load_matrix(first_path)
    if roi_names is None:
        roi_names = get_roi_names_from_any_labeled()
    if roi_names is None:
        roi_names = [f"ROI{i+1}" for i in range(mat0.shape[0])]

    if mat0.shape[0] != expected_n_roi:
        raise RuntimeError(f"Detected n_roi={mat0.shape[0]}, expected {expected_n_roi}. Please check DK68 matrix path.")

    degree_mat = []
    for fp in df["mind_file"]:
        mat, _ = load_matrix(fp)
        degree_mat.append(calc_degree(mat))
    degree_mat = np.asarray(degree_mat)

    for beh in diagnostic_behaviors:
        beh_col = find_col(df, [beh])
        if beh_col is None:
            continue
        for dx in ["DD", "TD"]:
            df_dx = df[df["Diagnosis"] == dx].copy()
            if df_dx.empty:
                continue
            degree_dx = degree_mat[df["Diagnosis"] == dx]
            out_base = os.path.join(output_dir, "diagnostic", slugify(beh), dx)
            res_dx = run_behavior_models(
                df_dx,
                degree_dx,
                roi_names,
                beh_col,
                formula="brain ~ behavior * AgeGroup + Sex",
                required_cols=["AgeGroup", "Sex"],
                group_col="AgeGroup",
                out_base=os.path.join(out_base, "degree_agegroup"),
                tag="behavior_agegroup",
                p_check_terms=["behavior", "behavior:AgeGroup"]
            )
            res_dx = add_fdr(res_dx)
            res_dx = add_fdr_interaction(res_dx, "p_interaction_agegroup", "agegroup")
            res_dx.to_csv(os.path.join(out_base, "behavior_degree_agegroup_results.csv"), index=False)
            res_dx[res_dx["significant_fdr"]].to_csv(
                os.path.join(out_base, "behavior_degree_agegroup_results_fdr_sig.csv"),
                index=False
            )
            res_dx[res_dx["significant_fdr_agegroup"]].to_csv(
                os.path.join(out_base, "behavior_degree_agegroup_results_fdr_sig_agegroup.csv"),
                index=False
            )

    beh_col = find_col(df, [nondiagnostic_behavior])
    if beh_col is not None:
        out_base = os.path.join(output_dir, "nondiagnostic", slugify(nondiagnostic_behavior))
        res_overall = run_behavior_models(
            df,
            degree_mat,
            roi_names,
            beh_col,
            formula="brain ~ behavior * Diagnosis * AgeGroup + Sex",
            required_cols=["Diagnosis", "AgeGroup", "Sex"],
            group_col="DxAgeGroup",
            out_base=os.path.join(out_base, "degree_diagnosis_agegroup"),
            tag="behavior_diagnosis_agegroup",
            p_check_terms=["behavior", "behavior:Diagnosis", "behavior:AgeGroup", "behavior:Diagnosis:AgeGroup"],
            legend_outside=True
        )
        res_overall = add_fdr(res_overall)
        res_overall = add_fdr_interaction(res_overall, "p_interaction_agegroup", "agegroup")
        res_overall = add_fdr_interaction(res_overall, "p_interaction_diagnosis", "diagnosis")
        res_overall = add_fdr_interaction(res_overall, "p_interaction_dx_age", "dx_age")
        res_overall.to_csv(os.path.join(out_base, "behavior_degree_diagnosis_agegroup_results.csv"), index=False)
        res_overall[res_overall["significant_fdr"]].to_csv(
            os.path.join(out_base, "behavior_degree_diagnosis_agegroup_results_fdr_sig.csv"),
            index=False
        )
        res_overall[res_overall["significant_fdr_agegroup"]].to_csv(
            os.path.join(out_base, "behavior_degree_diagnosis_agegroup_results_fdr_sig_agegroup.csv"),
            index=False
        )
        res_overall[res_overall["significant_fdr_diagnosis"]].to_csv(
            os.path.join(out_base, "behavior_degree_diagnosis_agegroup_results_fdr_sig_diagnosis.csv"),
            index=False
        )
        res_overall[res_overall["significant_fdr_dx_age"]].to_csv(
            os.path.join(out_base, "behavior_degree_diagnosis_agegroup_results_fdr_sig_dx_age.csv"),
            index=False
        )

        for dx in ["DD", "TD"]:
            df_dx = df[df["Diagnosis"] == dx].copy()
            if df_dx.empty:
                continue
            degree_dx = degree_mat[df["Diagnosis"] == dx]
            out_dx = os.path.join(out_base, dx)
            res_dx = run_behavior_models(
                df_dx,
                degree_dx,
                roi_names,
                beh_col,
                formula="brain ~ behavior * AgeGroup + Sex",
                required_cols=["AgeGroup", "Sex"],
                group_col="AgeGroup",
                out_base=os.path.join(out_dx, "degree_agegroup"),
                tag="behavior_agegroup",
                p_check_terms=["behavior", "behavior:AgeGroup"]
            )
            res_dx = add_fdr(res_dx)
            res_dx = add_fdr_interaction(res_dx, "p_interaction_agegroup", "agegroup")
            res_dx.to_csv(os.path.join(out_dx, "behavior_degree_agegroup_results.csv"), index=False)
            res_dx[res_dx["significant_fdr"]].to_csv(
                os.path.join(out_dx, "behavior_degree_agegroup_results_fdr_sig.csv"),
                index=False
            )
            res_dx[res_dx["significant_fdr_agegroup"]].to_csv(
                os.path.join(out_dx, "behavior_degree_agegroup_results_fdr_sig_agegroup.csv"),
                index=False
            )

if __name__ == "__main__":
    main()