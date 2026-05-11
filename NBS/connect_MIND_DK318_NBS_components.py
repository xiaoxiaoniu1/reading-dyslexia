import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from nilearn import plotting

DEFAULT_NBS_DIR = "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NBS_100"
DEFAULT_TOP_N = 50
EDGE_SUFFIXES = (
    "_edges_in_FWE_significant_components.csv",
    "_suprathreshold_component_edges.csv",
    "_FWE_significant_edges.csv",
)
PREFERRED_ALL = [
    "NBR_ALL_edges_in_FWE_significant_components.csv",
    "NBR_ALL_suprathreshold_component_edges.csv",
]

# Network color palette (matching Yeo7 networks)
NETWORK_COLORS = {
    "VIS": "#781286",      # Purple
    "SMN": "#4682B4",      # Steel Blue
    "DAN": "#00760E",      # Green
    "VAN": "#C43AFA",      # Light Purple
    "Limbic": "#DCDC00",   # Yellow
    "DMN": "#CD3E4E",      # Red
    "FPN": "#E69422",      # Orange
    "Other": "#808080",    # Gray
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Plot per-component top-N connectome figures from MIND_DK318_NBS outputs."
    )
    p.add_argument("--nbs-dir", default=DEFAULT_NBS_DIR)
    p.add_argument("--top-n", type=int, default=DEFAULT_TOP_N)
    p.add_argument("--display-mode", default="lyrz")
    p.add_argument("--node-size", type=float, default=2.0)
    p.add_argument("--edge-linewidth", type=float, default=1.2)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--analysis", nargs="*", default=None)
    p.add_argument("--include-overall", action="store_true")
    return p.parse_args()


def require_path(path: Path, desc: str) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"{desc} not found: {path}")
    return path


def discover_roi_metadata(nbs_dir: Path) -> Path:
    for path in [
        nbs_dir / "roi_metadata_with_mni_coordinates.csv",
        nbs_dir / "ROI_metadata_with_mni_coordinates.csv",
    ]:
        if path.exists():
            return path
    matches = sorted(nbs_dir.glob("*roi*metadata*.csv"))
    if matches:
        return matches[0]
    raise FileNotFoundError("Could not find ROI metadata with MNI coordinates.")


def discover_edge_tables(tables_dir: Path, include_overall: bool):
    out, seen = [], set()
    if include_overall:
        for name in PREFERRED_ALL:
            p = tables_dir / name
            if p.exists():
                out.append(p)
                seen.add(p)
    for p in sorted(tables_dir.glob("*.csv")):
        if not p.name.endswith(EDGE_SUFFIXES):
            continue
        if not include_overall and p.name.startswith("NBR_ALL_"):
            continue
        if p not in seen:
            out.append(p)
            seen.add(p)
    if not out:
        raise FileNotFoundError("Could not find any supported edge table.")
    return out


def sanitize_name(text: str) -> str:
    safe = "".join(ch if ch.isalnum() or ch in "_-" else "_" for ch in str(text))
    while "__" in safe:
        safe = safe.replace("__", "_")
    return safe.strip("_")


def infer_analysis_name(path: Path) -> str:
    name = path.name
    for s in EDGE_SUFFIXES:
        if name.endswith(s):
            return name[: -len(s)]
    return path.stem


def roi_lookup(roi_df: pd.DataFrame) -> pd.DataFrame:
    need = {"ROI", "x", "y", "z"}
    miss = need.difference(roi_df.columns)
    if miss:
        raise ValueError(f"ROI metadata missing columns: {sorted(miss)}")
    # Include network information if available
    cols = ["ROI", "x", "y", "z"]
    if "network" in roi_df.columns:
        cols.append("network")
    return roi_df[cols].drop_duplicates("ROI").set_index("ROI")


def normalize_table(df: pd.DataFrame, fallback: str, coords: pd.DataFrame) -> pd.DataFrame:
    if "analysis" not in df.columns:
        df["analysis"] = df["contrast"] if "contrast" in df.columns else fallback
    df["analysis"] = df["analysis"].fillna(fallback).astype(str)
    if "ROI1" not in df.columns and "ROI1.x" in df.columns:
        df = df.rename(columns={"ROI1.x": "ROI1"})
    if "ROI2" not in df.columns and "ROI2.x" in df.columns:
        df = df.rename(columns={"ROI2.x": "ROI2"})
    for c in ["analysis", "Component", "ROI1", "ROI2"]:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")
    
    # Add network information for ROI1 and ROI2
    if "network" in coords.columns:
        net1 = coords["network"].reindex(df["ROI1"]).to_numpy()
        net2 = coords["network"].reindex(df["ROI2"]).to_numpy()
        if "network1" not in df.columns:
            df["network1"] = net1
        if "network2" not in df.columns:
            df["network2"] = net2
    
    for roi_col, k in [("ROI1", "1"), ("ROI2", "2")]:
        x, y, z = f"x{k}", f"y{k}", f"z{k}"
        if not all(c in df.columns for c in [x, y, z]):
            xyz = coords.reindex(df[roi_col])
            df[x], df[y], df[z] = xyz["x"].to_numpy(), xyz["y"].to_numpy(), xyz["z"].to_numpy()
    if "t" in df.columns:
        df["t"] = pd.to_numeric(df["t"], errors="coerce")
        df["abs_t"] = df["t"].abs()
    elif "abs_t" in df.columns:
        df["abs_t"] = pd.to_numeric(df["abs_t"], errors="coerce")
    else:
        raise ValueError("Need 't' or 'abs_t' column")
    if "plot_weight" not in df.columns:
        df["plot_weight"] = df["abs_t"]
    df["plot_weight"] = pd.to_numeric(df["plot_weight"], errors="coerce")
    if "direction_for_plot" not in df.columns:
        if "direction" in df.columns:
            df["direction_for_plot"] = df["direction"]
        elif "t" in df.columns:
            df["direction_for_plot"] = np.where(df["t"] >= 0, "TD > DD", "DD > TD")
        else:
            df["direction_for_plot"] = np.where(
                df["analysis"].str.contains("DD_gt_TD|DD > TD", case=False, regex=True),
                "DD > TD",
                "TD > DD",
            )
    if "component_sig" in df.columns:
        keep = df["component_sig"]
        if keep.dtype != bool:
            keep = keep.astype(str).str.lower().isin({"true", "t", "1", "yes"})
        df = df.loc[keep].copy()
    elif "ncompFWE" in df.columns:
        fwe = pd.to_numeric(df["ncompFWE"], errors="coerce")
        df = df.loc[fwe < 0.05].copy()
    else:
        return pd.DataFrame(columns=df.columns)
    for c in ["x1", "y1", "z1", "x2", "y2", "z2", "plot_weight"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["Component"] = pd.to_numeric(df["Component"], errors="coerce")
    df = df.dropna(subset=["Component", "ROI1", "ROI2", "x1", "y1", "z1", "x2", "y2", "z2", "plot_weight"])
    df["Component"] = df["Component"].astype(int)
    return df.sort_values(["analysis", "Component", "plot_weight"], ascending=[True, True, False]).reset_index(drop=True)


def edge_color_by_network(net1: str, net2: str) -> str:
    """
    Determine edge color based on the networks of the two connected nodes.
    If both nodes are in the same network, use that network's color.
    If they are in different networks, use a blended or neutral color.
    """
    net1 = str(net1).strip() if pd.notna(net1) else "Other"
    net2 = str(net2).strip() if pd.notna(net2) else "Other"
    
    # If same network, use that network's color
    if net1 == net2:
        return NETWORK_COLORS.get(net1, NETWORK_COLORS["Other"])
    
    # If different networks, use a gray color to indicate inter-network connection
    return "#999999"


def node_color_by_network(network: str) -> str:
    """Get node color based on network."""
    network = str(network).strip() if pd.notna(network) else "Other"
    return NETWORK_COLORS.get(network, NETWORK_COLORS["Other"])


def build_graph(df: pd.DataFrame):
    rois = pd.Index(pd.unique(pd.concat([df["ROI1"], df["ROI2"]], ignore_index=True)))
    idx = {r: i for i, r in enumerate(rois)}
    adj = np.zeros((len(rois), len(rois)), dtype=float)
    edge_colors, node_colors, xyz = {}, {}, {}
    node_networks = {}
    
    for r in df.itertuples(index=False):
        i, j = idx[r.ROI1], idx[r.ROI2]
        w = float(r.plot_weight)
        adj[i, j] = adj[j, i] = w
        
        # Get network information
        net1 = getattr(r, "network1", "Other") if hasattr(r, "network1") else "Other"
        net2 = getattr(r, "network2", "Other") if hasattr(r, "network2") else "Other"
        
        # Determine edge color based on networks
        edge_colors[(i, j)] = edge_colors[(j, i)] = edge_color_by_network(net1, net2)
        
        # Store node networks and colors
        node_networks[r.ROI1] = net1
        node_networks[r.ROI2] = net2
        node_colors[i] = node_color_by_network(net1)
        node_colors[j] = node_color_by_network(net2)
        
        xyz[r.ROI1] = (float(r.x1), float(r.y1), float(r.z1))
        xyz[r.ROI2] = (float(r.x2), float(r.y2), float(r.z2))
    
    coords = np.array([xyz[r] for r in rois], dtype=float)
    node_color_array = [node_colors.get(i, NETWORK_COLORS["Other"]) for i in range(len(rois))]
    
    return rois.tolist(), adj, coords, edge_colors, node_color_array, node_networks


def plot_component(df: pd.DataFrame, args, out_png: Path):
    analysis = str(df["analysis"].iloc[0])
    comp = int(df["Component"].iloc[0])
    rois, adj, coords, edge_colors, node_colors, node_networks = build_graph(df)
    
    fig = plt.figure(figsize=(16, 9))
    
    # Plot nodes with network colors
    disp = plotting.plot_connectome(
        np.zeros_like(adj),
        coords,
        node_size=args.node_size,
        node_color=node_colors,
        edge_threshold=0,
        display_mode=args.display_mode,
        figure=fig,
        colorbar=False,
        black_bg=False,
        annotate=False,
    )
    
    vals = adj[np.triu_indices_from(adj, 1)]
    vmax = float(vals[vals > 0].max()) if np.any(vals > 0) else 1.0
    n = len(rois)
    
    # Plot edges with network-based colors
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i, j] <= 0:
                continue
            one = np.zeros_like(adj)
            one[i, j] = one[j, i] = adj[i, j]
            lw = args.edge_linewidth * (0.6 + 1.8 * adj[i, j] / vmax)
            disp.add_graph(
                one,
                coords,
                node_size=0,
                edge_threshold=0,
                edge_kwargs={"color": edge_colors[(i, j)], "linewidth": lw, "alpha": 0.85},
                colorbar=False,
            )
    
    # Add legend for networks
    from matplotlib.patches import Patch
    unique_networks = sorted(set(node_networks.values()))
    legend_elements = [
        Patch(facecolor=NETWORK_COLORS.get(net, NETWORK_COLORS["Other"]), 
              edgecolor='white', label=net)
        for net in unique_networks if net in NETWORK_COLORS
    ]
    if legend_elements:
        fig.legend(handles=legend_elements, loc='lower center', ncol=len(legend_elements), 
                  frameon=True, fontsize=10, bbox_to_anchor=(0.5, -0.02))
    
    fig.suptitle(f"{analysis} | Component {comp} | Top {len(df)} edges", fontsize=16, y=0.98)
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    fig.savefig(out_png, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)


def save_edges(df: pd.DataFrame, out_csv: Path):
    cols = [
        c
        for c in [
            "analysis",
            "Component",
            "ROI1",
            "ROI2",
            "i",
            "j",
            "t",
            "abs_t",
            "direction_for_plot",
            "plot_weight",
            "p",
            "p_fdr",
            "ncompFWE",
            "strnFWE",
        ]
        if c in df.columns
    ]
    df.loc[:, cols].to_csv(out_csv, index=False)


def main():
    args = parse_args()
    nbs_dir = require_path(Path(args.nbs_dir), "NBS output directory")
    tables_dir = require_path(nbs_dir / "tables", "NBS tables directory")
    fig_dir = nbs_dir / "figures" / "component_top_connectome"
    fig_dir.mkdir(parents=True, exist_ok=True)
    coords = roi_lookup(pd.read_csv(discover_roi_metadata(nbs_dir)))
    tables = discover_edge_tables(tables_dir, args.include_overall)
    keep_analysis = set(args.analysis) if args.analysis else None
    written = 0
    for table in tables:
        df = normalize_table(pd.read_csv(table), infer_analysis_name(table), coords)
        if keep_analysis is not None:
            df = df.loc[df["analysis"].isin(keep_analysis)].copy()
        if df.empty:
            continue
        for (analysis, comp), sub in df.groupby(["analysis", "Component"], sort=True):
            top = sub.sort_values("plot_weight", ascending=False).head(args.top_n).copy()
            if top.empty:
                continue
            stem = f"{sanitize_name(analysis)}_Component_{int(comp)}_top{len(top)}"
            save_edges(top, fig_dir / f"{stem}_edges.csv")
            plot_component(top, args, fig_dir / f"{stem}_connectome.png")
            written += 1
            print(f"Saved: {stem}")
    if written == 0:
        raise RuntimeError("No plottable components were found.")
    print(f"Finished. Wrote {written} figure(s) to {fig_dir}")


if __name__ == "__main__":
    main()
