#!/usr/bin/env python3
"""
绘制 T-map 聚类结果的脑区图
基于 K=3 聚类结果，为每个聚类生成脑区可视化
"""

import os
import numpy as np
import pandas as pd
import mne
from nibabel.freesurfer.io import read_annot, read_geometry, read_morph_data
from nilearn import plotting
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ============================================================
# 路径配置
# ============================================================
BASE_DIR = "/data/home/tqi/data1/share/after_freesurfer"
TMAP_DIR = os.path.join(BASE_DIR, "FILE/test_mean_1.5/Clustering_3")
CLUSTER_FILE = os.path.join(TMAP_DIR, "clustering/degree_Tmap_cluster_K3.csv")
OUTPUT_DIR = os.path.join(TMAP_DIR, "clustering/brain_maps")

LH_ANNOT = os.path.join(BASE_DIR, "FILE/DK-318", "lh.DK318.annot")
RH_ANNOT = os.path.join(BASE_DIR, "FILE/DK-318", "rh.DK318.annot")

# ============================================================
# 聚类配置
# ============================================================
CLUSTER_CONFIGS = [
    {
        "cluster_id": 1,
        "display_name": "Cluster 1: Adult Emergence",
        "short_name": "C1_Adult_Emergence",
        "subtitle": "Weak → Strong (Child: t=1.06, Adult: t=2.32)",
        "description": "115 ROIs showing adult enhancement pattern",
        "color": "#e74c3c",  # 红色
    },
    {
        "cluster_id": 2,
        "display_name": "Cluster 2: Child Predominant",
        "short_name": "C2_Child_Predominant",
        "subtitle": "Strong → Moderate (Child: t=2.44, Adult: t=1.23)",
        "description": "111 ROIs showing child predominant pattern",
        "color": "#3498db",  # 蓝色
    },
    {
        "cluster_id": 3,
        "display_name": "Cluster 3: Consistently Weak",
        "short_name": "C3_Consistently_Weak",
        "subtitle": "Weak → Weak (Child: t=0.48, Adult: t=0.31)",
        "description": "92 ROIs showing minimal TD-DD difference",
        "color": "#9b59b6",  # 紫色
    },
]

VIEWS = [
    ("lh", "lateral", "left"),
    ("rh", "lateral", "right"),
    ("lh", "medial", "left"),
    ("rh", "medial", "right"),
]


# ============================================================
# 辅助函数
# ============================================================
def load_fsaverage_and_annot():
    """加载 fsaverage 表面和 DK318 注释"""
    print("Loading fsaverage and DK318 annotations...")
    
    subjects_dir = os.path.join(BASE_DIR, "FILE")
    fs_dir = mne.datasets.fetch_fsaverage(subjects_dir=subjects_dir, verbose=False)
    
    lh_mesh = read_geometry(os.path.join(fs_dir, "surf", "lh.inflated"))
    rh_mesh = read_geometry(os.path.join(fs_dir, "surf", "rh.inflated"))
    lh_sulc = read_morph_data(os.path.join(fs_dir, "surf", "lh.sulc"))
    rh_sulc = read_morph_data(os.path.join(fs_dir, "surf", "rh.sulc"))
    
    lh_labels, _, lh_names = read_annot(LH_ANNOT)
    rh_labels, _, rh_names = read_annot(RH_ANNOT)
    
    print(f"  LH vertices: {lh_mesh[0].shape[0]}, labels: {len(lh_labels)}")
    print(f"  RH vertices: {rh_mesh[0].shape[0]}, labels: {len(rh_labels)}")
    
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
    """移除半球前缀"""
    value = str(value)
    return value.replace("lh.", "").replace("rh.", "").replace("lh_", "").replace("rh_", "")


def build_roi_lookup_keys(feature: str):
    """构建 ROI 查找键的多种变体"""
    feature = str(feature)
    normalized = feature.replace(".", "_")
    return [
        feature,
        normalized,
        strip_hemi_prefix(feature),
        strip_hemi_prefix(normalized),
    ]


def make_vertex_data(labels, names, roi_map):
    """将 ROI 映射到顶点数据"""
    name_list = [name.decode("utf-8") for name in names]
    vertex_data = np.full(labels.shape[0], np.nan, dtype=float)
    
    matched_names = []
    for label_id, roi_name in enumerate(name_list):
        if roi_name in roi_map:
            vertex_data[labels == label_id] = float(roi_map[roi_name])
            matched_names.append(roi_name)
    
    return vertex_data, matched_names


def build_single_color_cmap(color: str):
    """创建单色渐变色图"""
    return LinearSegmentedColormap.from_list(
        "cluster_color",
        ["#d8d8d8", color],
    )


# ============================================================
# 主要函数
# ============================================================
def load_cluster_data():
    """加载聚类结果"""
    print(f"\nLoading cluster data from: {CLUSTER_FILE}")
    df = pd.read_csv(CLUSTER_FILE)
    print(f"  Total ROIs: {len(df)}")
    print(f"  Clusters: {sorted(df['cluster'].unique())}")
    
    # 统计每个聚类的 ROI 数量
    for cluster_id in sorted(df['cluster'].unique()):
        count = len(df[df['cluster'] == cluster_id])
        print(f"    Cluster {cluster_id}: {count} ROIs")
    
    return df


def prepare_cluster_dataframe(cluster_df, cluster_id):
    """准备特定聚类的数据"""
    cluster_rois = cluster_df[cluster_df['cluster'] == cluster_id].copy()
    cluster_rois['value'] = 1.0
    return cluster_rois


def build_hemi_maps(cluster_rois):
    """构建左右半球的 ROI 映射"""
    lh_rois = cluster_rois[
        cluster_rois['feature'].str.startswith('lh_') | 
        cluster_rois['feature'].str.startswith('lh.')
    ].copy()
    
    rh_rois = cluster_rois[
        cluster_rois['feature'].str.startswith('rh_') | 
        cluster_rois['feature'].str.startswith('rh.')
    ].copy()
    
    lh_map = {}
    for _, row in lh_rois.iterrows():
        for key in build_roi_lookup_keys(row['feature']):
            lh_map[key] = row['value']
    
    rh_map = {}
    for _, row in rh_rois.iterrows():
        for key in build_roi_lookup_keys(row['feature']):
            rh_map[key] = row['value']
    
    return lh_rois, rh_rois, lh_map, rh_map


def render_cluster_views(config, cluster_rois, fs_data, output_dir):
    """渲染单个聚类的脑区图"""
    print(f"\n{'='*60}")
    print(f"Rendering {config['display_name']}")
    print(f"{'='*60}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    lh_rois, rh_rois, lh_map, rh_map = build_hemi_maps(cluster_rois)
    
    lh_vtx, lh_matched = make_vertex_data(
        fs_data['lh_labels'], fs_data['lh_names'], lh_map
    )
    rh_vtx, rh_matched = make_vertex_data(
        fs_data['rh_labels'], fs_data['rh_names'], rh_map
    )
    
    print(f"  LH ROIs: {len(lh_rois)}, matched: {len(lh_matched)}")
    print(f"  RH ROIs: {len(rh_rois)}, matched: {len(rh_matched)}")
    
    # 保存 ROI 列表
    roi_list_path = os.path.join(output_dir, f"{config['short_name']}_ROIs.csv")
    cluster_rois[['feature', 't_child_TD_minus_DD', 't_adult_TD_minus_DD', 'cluster']].to_csv(
        roi_list_path, index=False
    )
    print(f"  Saved ROI list: {roi_list_path}")
    
    # 创建色图
    cmap = build_single_color_cmap(config['color'])
    
    vertex_data_by_hemi = {'lh': lh_vtx, 'rh': rh_vtx}
    mesh_by_hemi = {'lh': fs_data['lh_mesh'], 'rh': fs_data['rh_mesh']}
    sulc_by_hemi = {'lh': fs_data['lh_sulc'], 'rh': fs_data['rh_sulc']}
    
    # 渲染各个视图
    for hemi_tag, view, hemi_name in VIEWS:
        out_file = os.path.join(
            output_dir, 
            f"{config['short_name']}_{hemi_tag}_{view}.png"
        )
        
        fig = plt.figure(figsize=(5, 4.5), facecolor='white')
        axis = fig.add_subplot(111, projection='3d')
        
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
        
        # 添加标题
        n_rois = len(cluster_rois)
        fig.suptitle(
            f"{config['display_name']}\n{config['subtitle']}\n{config['description']}", 
            fontsize=9, 
            y=0.98
        )
        
        fig.savefig(
            out_file, 
            bbox_inches='tight', 
            pad_inches=0.1, 
            dpi=300, 
            facecolor='white'
        )
        plt.close(fig)
        print(f"  Saved: {os.path.basename(out_file)}")
    
    return {
        'config': config,
        'cluster_rois': cluster_rois,
        'lh_vtx': lh_vtx,
        'rh_vtx': rh_vtx,
    }


def create_summary_figure(cluster_results, fs_data):
    """创建 3 个聚类的汇总图"""
    print(f"\n{'='*60}")
    print("Creating summary figure...")
    print(f"{'='*60}")
    
    summary_path = os.path.join(OUTPUT_DIR, "All_Clusters_Summary.png")
    
    fig = plt.figure(figsize=(14, 10), facecolor='white')
    
    # 添加总标题
    fig.suptitle(
        'K-means Clustering (K=3) of TD-DD Degree T-maps\nDevelopmental Patterns Across Child and Adult Groups',
        fontsize=14, fontweight='bold', y=0.98
    )
    
    # 创建 3x2 的网格（3个聚类 × 2个半球）
    for row_idx, result in enumerate(cluster_results):
        config = result['config']
        n_rois = len(result['cluster_rois'])
        
        # 左半球 lateral view
        ax_lh = fig.add_subplot(3, 2, row_idx*2 + 1, projection='3d')
        plotting.plot_surf_stat_map(
            surf_mesh=fs_data['lh_mesh'],
            stat_map=result['lh_vtx'],
            hemi='left',
            view='lateral',
            bg_map=fs_data['lh_sulc'],
            bg_on_data=True,
            colorbar=False,
            threshold=0.5,
            vmin=0.0,
            vmax=1.0,
            cmap=build_single_color_cmap(config['color']),
            axes=ax_lh,
            darkness=0.55,
        )
        ax_lh.set_title(
            f"{config['display_name']} - LH\n{config['subtitle']}\n(n={n_rois} ROIs)", 
            fontsize=9, pad=10
        )
        
        # 右半球 lateral view
        ax_rh = fig.add_subplot(3, 2, row_idx*2 + 2, projection='3d')
        plotting.plot_surf_stat_map(
            surf_mesh=fs_data['rh_mesh'],
            stat_map=result['rh_vtx'],
            hemi='right',
            view='lateral',
            bg_map=fs_data['rh_sulc'],
            bg_on_data=True,
            colorbar=False,
            threshold=0.5,
            vmin=0.0,
            vmax=1.0,
            cmap=build_single_color_cmap(config['color']),
            axes=ax_rh,
            darkness=0.55,
        )
        ax_rh.set_title(
            f"{config['display_name']} - RH\n{config['subtitle']}\n(n={n_rois} ROIs)", 
            fontsize=9, pad=10
        )
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(summary_path, dpi=300, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    
    print(f"  Saved: {summary_path}")


# ============================================================
# 主程序
# ============================================================
def main():
    print("="*60)
    print("T-map Cluster Brain Visualization (K=3)")
    print("="*60)
    print("\nCluster Characteristics:")
    print("  Cluster 1: Adult Emergence (n=115)")
    print("    - Child: t=1.06, Adult: t=2.32")
    print("    - Pattern: Weak → Strong difference")
    print("  Cluster 2: Child Predominant (n=111)")
    print("    - Child: t=2.44, Adult: t=1.23")
    print("    - Pattern: Strong → Moderate difference")
    print("  Cluster 3: Consistently Weak (n=92)")
    print("    - Child: t=0.48, Adult: t=0.31")
    print("    - Pattern: Minimal TD-DD difference")
    print("="*60)
    
    # 加载数据
    fs_data = load_fsaverage_and_annot()
    cluster_df = load_cluster_data()
    
    # 为每个聚类生成脑区图
    cluster_results = []
    for config in CLUSTER_CONFIGS:
        cluster_rois = prepare_cluster_dataframe(cluster_df, config['cluster_id'])
        
        cluster_output_dir = os.path.join(OUTPUT_DIR, config['short_name'])
        
        result = render_cluster_views(
            config, 
            cluster_rois, 
            fs_data, 
            cluster_output_dir
        )
        cluster_results.append(result)
    
    # 创建汇总图
    create_summary_figure(cluster_results, fs_data)
    
    print("\n" + "="*60)
    print("All visualizations completed!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("="*60)


if __name__ == "__main__":
    main()
