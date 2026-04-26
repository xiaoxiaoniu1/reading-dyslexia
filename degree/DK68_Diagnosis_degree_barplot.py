#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
绘制显著 Diagnosis degree 脑区的柱状图
每个脑区分为四组：Adult-TD, Adult-DD, Child-TD, Child-DD
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# =============================================================================
# 1) 路径配置
# =============================================================================
demo_file = "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/all_data_cqt_mean_1.5.xlsx"
mind_csv_dir = "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK68_combat/"
sig_degree_file = "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK68_ANOVA/Significant_Diagnosis_DK68_degree_results.csv"
out_dir = "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK68_ANOVA/degree_diagnosis_barplot"

os.makedirs(out_dir, exist_ok=True)

# =============================================================================
# 2) 读取人口学数据
# =============================================================================
print("Loading demographic data...")
df = pd.read_excel(demo_file, sheet_name="Sheet1")

# 构建文件路径
df['id'] = df['id'].astype(str)
df['file_base'] = df['id'] + '_combat'
df['mind_file'] = df['file_base'].apply(lambda x: os.path.join(mind_csv_dir, f"{x}.csv"))

# 编码分组变量
df['Diagnosis'] = df['group_d_or_c'].map({0: 'TD', 1: 'DD'})
df['AgeGroup'] = df['group_age'].map({1: 'Adult', 2: 'Child'})
df['Sex'] = df['sex'].map({1: 'Male', 2: 'Female'})
df['Site'] = df['site'].astype('category')

# 创建四组标签
df['Group'] = df['AgeGroup'] + '-' + df['Diagnosis']

# 检查文件是否存在
df['has_file'] = df['mind_file'].apply(os.path.exists)

# 清理数据
df = df.dropna(subset=['id', 'Diagnosis', 'AgeGroup', 'Sex'])
df = df[df['has_file']]

print(f"Total subjects: {len(df)}")
print(f"  - Adult-TD: {len(df[(df['AgeGroup']=='Adult') & (df['Diagnosis']=='TD')])}")
print(f"  - Adult-DD: {len(df[(df['AgeGroup']=='Adult') & (df['Diagnosis']=='DD')])}")
print(f"  - Child-TD: {len(df[(df['AgeGroup']=='Child') & (df['Diagnosis']=='TD')])}")
print(f"  - Child-DD: {len(df[(df['AgeGroup']=='Child') & (df['Diagnosis']=='DD')])}")

# =============================================================================
# 3) 读取显著 degree 结果
# =============================================================================
print("\nLoading significant degree results...")
sig_df = pd.read_csv(sig_degree_file)
print(f"Number of significant ROIs: {len(sig_df)}")

# 获取显著 ROI 列表
sig_rois = sig_df['feature'].tolist()

# =============================================================================
# 4) 读取 degree 数据
# =============================================================================
def read_mind_matrix(filepath):
    """读取 MIND 矩阵"""
    mat = pd.read_csv(filepath, index_col=0)
    return mat

def calc_degree(mat):
    """计算度中心性（排除对角线）"""
    mat_copy = mat.values.copy()
    np.fill_diagonal(mat_copy, np.nan)
    return np.nanmean(mat_copy, axis=1)

# 读取第一个矩阵获取 ROI 信息
first_file = df['mind_file'].iloc[0]
first_mat = read_mind_matrix(first_file)
roi_names = first_mat.columns.tolist()
n_roi = len(roi_names)

print(f"\nTotal ROIs in matrix: {n_roi}")

# 构建 degree 矩阵
print("Loading degree data...")
degree_data = []
for idx, row in df.iterrows():
    mat = read_mind_matrix(row['mind_file'])
    deg = calc_degree(mat)
    degree_data.append(deg)

degree_mat = np.array(degree_data)  # shape: (n_subj, n_roi)
print(f"Degree matrix shape: {degree_mat.shape}")

# =============================================================================
# 5) 为每个显著 ROI 准备长表数据
# =============================================================================
print("\nPreparing long-form data for significant ROIs...")

groups = ['Adult-TD', 'Adult-DD', 'Child-TD', 'Child-DD']
plot_data = []

for roi in sig_rois:
    if roi not in roi_names:
        print(f"Warning: ROI {roi} not found in matrix")
        continue
    
    roi_idx = roi_names.index(roi)
    roi_degree = degree_mat[:, roi_idx]
    
    for group in groups:
        group_mask = df['Group'] == group
        group_degree = roi_degree[group_mask]
        if len(group_degree) == 0:
            continue
        for val in group_degree:
            plot_data.append({
                'ROI': roi,
                'Group': group,
                'Degree': val
            })

plot_df = pd.DataFrame(plot_data)
print(f"\nPlot data prepared: {len(plot_df)} rows")

# =============================================================================
# 6) 为前 10 个显著 ROI 分别绘制小提琴图
# =============================================================================
print("\nGenerating individual violin plots for top 10 ROIs...")

# 选择前 10 个最显著的 ROI
n_rois_to_plot = min(10, len(sig_rois))
rois_to_plot = sig_rois[:n_rois_to_plot]

# 定义颜色（加重饱和度）
colors = {
    'Adult-TD': '#2E5C99',   # 深蓝色
    'Adult-DD': '#D96519',   # 深橙色
    'Child-TD': '#5A8F35',   # 深绿色
    'Child-DD': '#E6A800'    # 深黄色
}

# 为每个 ROI 生成单独的图
for roi in rois_to_plot:
    # 获取该 ROI 的数据
    roi_data = plot_df[plot_df['ROI'] == roi]
    
    if len(roi_data) == 0:
        print(f"Warning: No data for ROI {roi}")
        continue
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(8, 6))
    
    sns.violinplot(
        data=roi_data,
        x='Group',
        y='Degree',
        order=groups,
        palette=colors,
        cut=0,
        inner=None,
        linewidth=1,
        ax=ax
    )
    sns.stripplot(
        data=roi_data,
        x='Group',
        y='Degree',
        order=groups,
        color='black',
        size=3,
        alpha=0.5,
        jitter=0.2,
        ax=ax
    )
    
    # 设置 x 轴标签
    ax.set_xticklabels(groups, fontsize=11, fontweight='bold')
    
    # 设置标签和标题
    roi_label = roi.replace('lh_', 'L-').replace('rh_', 'R-').replace('_part', '-')
    ax.set_ylabel('Degree Centrality', fontsize=12, fontweight='bold')
    ax.set_title(f'{roi_label}\nDegree Centrality by Group', 
                 fontsize=13, fontweight='bold', pad=15)
    
    # 添加网格
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)
    ax.set_axisbelow(True)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    roi_safe_name = roi.replace('/', '_').replace('\\', '_')
    output_file = os.path.join(out_dir, f"Significant_Diagnosis_DK68_degree_{roi_safe_name}.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {roi_safe_name}.png")

# 保存绘图数据
plot_data_file = os.path.join(out_dir, "Significant_Diagnosis_DK68_degree_barplot_data.csv")
plot_df[plot_df['ROI'].isin(rois_to_plot)].to_csv(plot_data_file, index=False)
print(f"\nPlot data saved: {plot_data_file}")

print(f"\nGenerated {n_rois_to_plot} individual barplots!")
print("DONE!")
