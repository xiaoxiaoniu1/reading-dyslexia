#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分层分析：在每个年龄组内单独比较 DD vs TD
检查诊断效应在 Child 和 Adult 组中的模式
"""

import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
import os
from pathlib import Path

# =============================================================================
# 1) 配置路径
# =============================================================================
demo_file = "/data/home/tqi/data1/share/after_freesurfer/FILE/test_any_2/all_data_cqt_any_2.xlsx"
mind_csv_dir = "/data/home/tqi/data1/share/after_freesurfer/FILE/test_any_2/MIND_DK68_combat/"
out_dir = "/data/home/tqi/data1/share/after_freesurfer/FILE/test_any_2/MIND_t_test_DK68"

os.makedirs(out_dir, exist_ok=True)

out_degree_child = os.path.join(out_dir, "Stratified_Child_degree_results.csv")
out_degree_adult = os.path.join(out_dir, "Stratified_Adult_degree_results.csv")
out_edge_child = os.path.join(out_dir, "Stratified_Child_edge_results.csv")
out_edge_adult = os.path.join(out_dir, "Stratified_Adult_edge_results.csv")
out_summary = os.path.join(out_dir, "Stratified_Analysis_Summary.txt")

# =============================================================================
# 2) 读取人口学数据
# =============================================================================
print("Loading demographic data...")
df = pd.read_excel(demo_file, sheet_name="Sheet1")

# 构建文件路径（DK68：id_combat.csv）
df['subj_id'] = df['id'].astype(str).str.strip()
df['file_base'] = df['subj_id'] + '_combat'
df['mind_file'] = df['file_base'].apply(lambda x: os.path.join(mind_csv_dir, f"{x}.csv"))

# 编码分组变量
df['Diagnosis'] = df['group_d_or_c'].map({0: 'TD', 1: 'DD'})
df['AgeGroup'] = df['group_age'].map({1: 'Adult', 2: 'Child'})
df['Sex'] = df['sex'].map({1: 'Male', 2: 'Female'})
df['Site'] = df['site'].astype('category')

# 检查文件是否存在
df['has_file'] = df['mind_file'].apply(os.path.exists)

# 清理数据
df = df.dropna(subset=['id', 'Diagnosis', 'AgeGroup', 'Sex'])
df = df[df['subj_id'] != ""]
df = df[df['has_file']]

print(f"Total subjects with valid data: {len(df)}")
print(f"  - Child TD: {len(df[(df['AgeGroup']=='Child') & (df['Diagnosis']=='TD')])}")
print(f"  - Child DD: {len(df[(df['AgeGroup']=='Child') & (df['Diagnosis']=='DD')])}")
print(f"  - Adult TD: {len(df[(df['AgeGroup']=='Adult') & (df['Diagnosis']=='TD')])}")
print(f"  - Adult DD: {len(df[(df['AgeGroup']=='Adult') & (df['Diagnosis']=='DD')])}")

# =============================================================================
# 3) 读取矩阵数据
# =============================================================================
def read_mind_matrix(filepath):
    """读取 MIND 矩阵"""
    mat = pd.read_csv(filepath, index_col=0)
    return mat.values

def calc_degree(mat):
    """计算度中心性（排除对角线）"""
    mat_copy = mat.copy()
    np.fill_diagonal(mat_copy, np.nan)
    return np.nanmean(mat_copy, axis=1)

# 读取第一个矩阵获取 ROI 信息
first_file = df['mind_file'].iloc[0]
first_mat = pd.read_csv(first_file, index_col=0)
roi_names = first_mat.columns.tolist()
n_roi = len(roi_names)

print(f"\nDetected {n_roi} ROIs")

# 获取上三角索引
upper_tri_idx = np.triu_indices(n_roi, k=1)
n_edge = len(upper_tri_idx[0])
print(f"Number of edges (upper triangle): {n_edge}")

# 构建 degree 和 edge 矩阵
print("\nLoading matrices...")
degree_data = []
edge_data = []

for idx, row in df.iterrows():
    mat = read_mind_matrix(row['mind_file'])
    
    # Degree
    deg = calc_degree(mat)
    degree_data.append(deg)
    
    # Edges
    edges = mat[upper_tri_idx]
    edge_data.append(edges)

degree_mat = np.array(degree_data)  # shape: (n_subj, n_roi)
edge_mat = np.array(edge_data)      # shape: (n_subj, n_edge)

print(f"Degree matrix shape: {degree_mat.shape}")
print(f"Edge matrix shape: {edge_mat.shape}")

# =============================================================================
# 4) 分层分析函数：在特定年龄组内比较 DD vs TD
# =============================================================================
def run_stratified_analysis(feature_mat, df_demo, feature_names, group_name, 
                           min_n=10, var_eps=1e-12):
    """
    在指定年龄组内进行 DD vs TD 比较
    
    Parameters:
    -----------
    feature_mat : ndarray, shape (n_subj, n_features)
    df_demo : DataFrame with demographic info
    feature_names : list of feature names
    group_name : 'Child' or 'Adult'
    
    Returns:
    --------
    DataFrame with results for each feature
    """
    
    # 筛选该年龄组
    group_mask = df_demo['AgeGroup'] == group_name
    df_group = df_demo[group_mask].reset_index(drop=True)
    feat_group = feature_mat[group_mask]
    
    print(f"\n=== Analyzing {group_name} group (n={len(df_group)}) ===")
    
    results = []
    
    for k, feat_name in enumerate(feature_names):
        y = feat_group[:, k]
        
        # 排除 NaN/Inf
        valid_mask = np.isfinite(y)
        y_valid = y[valid_mask]
        df_valid = df_group[valid_mask].copy()
        
        # 检查样本量
        if len(y_valid) < min_n:
            results.append({
                'feature': feat_name,
                'n_total': len(y_valid),
                'n_TD': np.nan,
                'n_DD': np.nan,
                'mean_TD': np.nan,
                'mean_DD': np.nan,
                'std_TD': np.nan,
                'std_DD': np.nan,
                'cohens_d': np.nan,
                'p_value': np.nan,
                'direction': 'NA',
                'note': 'insufficient_sample'
            })
            continue
        
        # 检查方差
        if np.var(y_valid) < var_eps:
            results.append({
                'feature': feat_name,
                'n_total': len(y_valid),
                'n_TD': np.nan,
                'n_DD': np.nan,
                'mean_TD': np.nan,
                'mean_DD': np.nan,
                'std_TD': np.nan,
                'std_DD': np.nan,
                'cohens_d': np.nan,
                'p_value': np.nan,
                'direction': 'NA',
                'note': 'zero_variance'
            })
            continue
        
        # 准备数据
        df_valid['y'] = y_valid
        
        # 检查两组是否都有数据
        n_td = (df_valid['Diagnosis'] == 'TD').sum()
        n_dd = (df_valid['Diagnosis'] == 'DD').sum()
        
        if n_td < 3 or n_dd < 3:
            results.append({
                'feature': feat_name,
                'n_total': len(y_valid),
                'n_TD': n_td,
                'n_DD': n_dd,
                'mean_TD': np.nan,
                'mean_DD': np.nan,
                'std_TD': np.nan,
                'std_DD': np.nan,
                'cohens_d': np.nan,
                'p_value': np.nan,
                'direction': 'NA',
                'note': 'insufficient_group_size'
            })
            continue
        
        # 计算描述统计
        mean_td = df_valid[df_valid['Diagnosis']=='TD']['y'].mean()
        mean_dd = df_valid[df_valid['Diagnosis']=='DD']['y'].mean()
        std_td = df_valid[df_valid['Diagnosis']=='TD']['y'].std()
        std_dd = df_valid[df_valid['Diagnosis']=='DD']['y'].std()
        
        # Cohen's d 效应量
        pooled_std = np.sqrt(((n_td-1)*std_td**2 + (n_dd-1)*std_dd**2) / (n_td + n_dd - 2))
        cohens_d = (mean_dd - mean_td) / pooled_std if pooled_std > 0 else np.nan
        
        # 效应方向
        if mean_dd > mean_td:
            direction = 'DD>TD'
        elif mean_dd < mean_td:
            direction = 'DD<TD'
        else:
            direction = 'equal'
        
        # 线性模型：y ~ Diagnosis + Sex (控制性别协变量)
        # 注意：如果要控制 age_month 和 Site，可以添加
        try:
            # 编码 Diagnosis 为数值 (TD=0, DD=1)
            df_valid['Diagnosis_num'] = (df_valid['Diagnosis'] == 'DD').astype(int)
            df_valid['Sex_num'] = (df_valid['Sex'] == 'Male').astype(int)
            
            # 构建设计矩阵
            X = df_valid[['Diagnosis_num', 'Sex_num']].values
            X = sm.add_constant(X)  # 添加截距
            
            model = sm.OLS(y_valid, X)
            fit_result = model.fit()
            
            # Diagnosis 的 p 值（第1个系数，因为第0个是截距）
            p_diagnosis = fit_result.pvalues[1]
            
        except Exception as e:
            p_diagnosis = np.nan
        
        results.append({
            'feature': feat_name,
            'n_total': len(y_valid),
            'n_TD': n_td,
            'n_DD': n_dd,
            'mean_TD': mean_td,
            'mean_DD': mean_dd,
            'std_TD': std_td,
            'std_DD': std_dd,
            'cohens_d': cohens_d,
            'p_value': p_diagnosis,
            'direction': direction,
            'note': 'success'
        })
    
    # 转为 DataFrame
    res_df = pd.DataFrame(results)
    
    # FDR 校正
    valid_p = res_df['p_value'].notna()
    if valid_p.sum() > 0:
        from statsmodels.stats.multitest import multipletests
        res_df.loc[valid_p, 'p_fdr'] = multipletests(
            res_df.loc[valid_p, 'p_value'], 
            method='fdr_bh'
        )[1]
    else:
        res_df['p_fdr'] = np.nan
    
    # 显著性标记
    res_df['significant'] = (res_df['p_fdr'] < 0.05) & res_df['p_fdr'].notna()
    
    return res_df

# =============================================================================
# 5) 运行分层分析 - Degree
# =============================================================================
print("\n" + "="*70)
print("DEGREE ANALYSIS")
print("="*70)

# Child 组
degree_child_res = run_stratified_analysis(
    degree_mat, df, roi_names, 'Child'
)
degree_child_res.to_csv(out_degree_child, index=False)
print(f"\nChild group - Significant ROIs (FDR<0.05): {degree_child_res['significant'].sum()}")
print(f"Saved: {out_degree_child}")

# Adult 组
degree_adult_res = run_stratified_analysis(
    degree_mat, df, roi_names, 'Adult'
)
degree_adult_res.to_csv(out_degree_adult, index=False)
print(f"\nAdult group - Significant ROIs (FDR<0.05): {degree_adult_res['significant'].sum()}")
print(f"Saved: {out_degree_adult}")

# =============================================================================
# 6) 运行分层分析 - Edge (分块处理)
# =============================================================================
print("\n" + "="*70)
print("EDGE ANALYSIS")
print("="*70)

# 创建 edge 名称
edge_names = [f"Edge_{upper_tri_idx[0][i]}_{upper_tri_idx[1][i]}" 
              for i in range(n_edge)]

def run_edge_analysis_chunked(edge_mat, df_demo, edge_names, upper_tri_idx, 
                               group_name, chunk_size=2000):
    """分块处理 edge 分析"""
    n_edge = len(edge_names)
    n_chunks = int(np.ceil(n_edge / chunk_size))
    
    all_results = []
    
    for chunk_idx in range(n_chunks):
        start_idx = chunk_idx * chunk_size
        end_idx = min((chunk_idx + 1) * chunk_size, n_edge)
        
        print(f"Processing edge chunk {chunk_idx+1}/{n_chunks} "
              f"(edges {start_idx}-{end_idx})...")
        
        chunk_mat = edge_mat[:, start_idx:end_idx]
        chunk_names = edge_names[start_idx:end_idx]
        
        chunk_res = run_stratified_analysis(
            chunk_mat, df_demo, chunk_names, group_name
        )
        
        # 添加 edge 索引信息
        for i, idx in enumerate(range(start_idx, end_idx)):
            chunk_res.loc[i, 'edge_id'] = idx
            chunk_res.loc[i, 'roi_i'] = upper_tri_idx[0][idx]
            chunk_res.loc[i, 'roi_j'] = upper_tri_idx[1][idx]
        
        all_results.append(chunk_res)
    
    return pd.concat(all_results, ignore_index=True)

# Child 组 - Edge
edge_child_res = run_edge_analysis_chunked(
    edge_mat, df, edge_names, upper_tri_idx, 'Child', chunk_size=2000
)
edge_child_res.to_csv(out_edge_child, index=False)
print(f"\nChild group - Significant Edges (FDR<0.05): {edge_child_res['significant'].sum()}")
print(f"Saved: {out_edge_child}")

# Adult 组 - Edge
edge_adult_res = run_edge_analysis_chunked(
    edge_mat, df, edge_names, upper_tri_idx, 'Adult', chunk_size=2000
)
edge_adult_res.to_csv(out_edge_adult, index=False)
print(f"\nAdult group - Significant Edges (FDR<0.05): {edge_adult_res['significant'].sum()}")
print(f"Saved: {out_edge_adult}")

# =============================================================================
# 7) 生成汇总报告
# =============================================================================
with open(out_summary, 'w', encoding='utf-8') as f:
    f.write("="*70 + "\n")
    f.write("分层分析汇总报告：DD vs TD 在不同年龄组中的比较\n")
    f.write("="*70 + "\n\n")
    
    # Degree 结果
    f.write("【一】Degree (节点度中心性) 分析\n")
    f.write("-"*70 + "\n\n")
    
    f.write("1. Child 组 (儿童)\n")
    f.write(f"   总 ROI 数: {len(degree_child_res)}\n")
    f.write(f"   显著 ROI (FDR<0.05): {degree_child_res['significant'].sum()}\n")
    if degree_child_res['significant'].sum() > 0:
        sig_child = degree_child_res[degree_child_res['significant']].sort_values('p_fdr')
        f.write(f"   Top 10 显著 ROI:\n")
        for _, row in sig_child.head(10).iterrows():
            f.write(f"     - {row['feature']}: p_fdr={row['p_fdr']:.2e}, "
                   f"Cohen's d={row['cohens_d']:.3f}, {row['direction']}\n")
    f.write("\n")
    
    f.write("2. Adult 组 (成人)\n")
    f.write(f"   总 ROI 数: {len(degree_adult_res)}\n")
    f.write(f"   显著 ROI (FDR<0.05): {degree_adult_res['significant'].sum()}\n")
    if degree_adult_res['significant'].sum() > 0:
        sig_adult = degree_adult_res[degree_adult_res['significant']].sort_values('p_fdr')
        f.write(f"   Top 10 显著 ROI:\n")
        for _, row in sig_adult.head(10).iterrows():
            f.write(f"     - {row['feature']}: p_fdr={row['p_fdr']:.2e}, "
                   f"Cohen's d={row['cohens_d']:.3f}, {row['direction']}\n")
    f.write("\n")
    
    # 比较两组的效应方向
    f.write("3. 效应方向比较\n")
    if degree_child_res['significant'].sum() > 0 and degree_adult_res['significant'].sum() > 0:
        sig_child_rois = set(degree_child_res[degree_child_res['significant']]['feature'])
        sig_adult_rois = set(degree_adult_res[degree_adult_res['significant']]['feature'])
        common_rois = sig_child_rois & sig_adult_rois
        
        f.write(f"   两组共同显著的 ROI: {len(common_rois)}\n")
        if len(common_rois) > 0:
            for roi in common_rois:
                child_dir = degree_child_res[degree_child_res['feature']==roi]['direction'].iloc[0]
                adult_dir = degree_adult_res[degree_adult_res['feature']==roi]['direction'].iloc[0]
                consistency = "一致" if child_dir == adult_dir else "不一致"
                f.write(f"     - {roi}: Child={child_dir}, Adult={adult_dir} ({consistency})\n")
    f.write("\n\n")
    
    # Edge 结果
    f.write("【二】Edge (连接边) 分析\n")
    f.write("-"*70 + "\n\n")
    
    f.write("1. Child 组 (儿童)\n")
    f.write(f"   总 Edge 数: {len(edge_child_res)}\n")
    f.write(f"   显著 Edge (FDR<0.05): {edge_child_res['significant'].sum()}\n")
    if edge_child_res['significant'].sum() > 0:
        sig_edge_child = edge_child_res[edge_child_res['significant']].sort_values('p_fdr')
        f.write(f"   Top 10 显著 Edge:\n")
        for _, row in sig_edge_child.head(10).iterrows():
            f.write(f"     - {row['feature']}: p_fdr={row['p_fdr']:.2e}, "
                   f"Cohen's d={row['cohens_d']:.3f}, {row['direction']}\n")
    f.write("\n")
    
    f.write("2. Adult 组 (成人)\n")
    f.write(f"   总 Edge 数: {len(edge_adult_res)}\n")
    f.write(f"   显著 Edge (FDR<0.05): {edge_adult_res['significant'].sum()}\n")
    if edge_adult_res['significant'].sum() > 0:
        sig_edge_adult = edge_adult_res[edge_adult_res['significant']].sort_values('p_fdr')
        f.write(f"   Top 10 显著 Edge:\n")
        for _, row in sig_edge_adult.head(10).iterrows():
            f.write(f"     - {row['feature']}: p_fdr={row['p_fdr']:.2e}, "
                   f"Cohen's d={row['cohens_d']:.3f}, {row['direction']}\n")
    f.write("\n")
    
    # Edge 效应方向比较
    f.write("3. 效应方向比较\n")
    if edge_child_res['significant'].sum() > 0 and edge_adult_res['significant'].sum() > 0:
        sig_child_edges = set(edge_child_res[edge_child_res['significant']]['feature'])
        sig_adult_edges = set(edge_adult_res[edge_adult_res['significant']]['feature'])
        common_edges = sig_child_edges & sig_adult_edges
        
        f.write(f"   两组共同显著的 Edge: {len(common_edges)}\n")
        if len(common_edges) > 0:
            for edge in list(common_edges)[:20]:  # 只显示前20个
                child_dir = edge_child_res[edge_child_res['feature']==edge]['direction'].iloc[0]
                adult_dir = edge_adult_res[edge_adult_res['feature']==edge]['direction'].iloc[0]
                consistency = "一致" if child_dir == adult_dir else "不一致"
                f.write(f"     - {edge}: Child={child_dir}, Adult={adult_dir} ({consistency})\n")
    f.write("\n")

print(f"\n{'='*70}")
print(f"Summary report saved: {out_summary}")
print(f"{'='*70}")

# =============================================================================
# 8) 输出显著结果表（degree 和 edge）
# =============================================================================
print("\n" + "="*70)
print("GENERATING SIGNIFICANT RESULTS TABLES")
print("="*70)

# Degree 显著结果 - Child
sig_degree_child = degree_child_res[degree_child_res['significant']].copy()
if len(sig_degree_child) > 0:
    sig_degree_child_output = sig_degree_child[['feature', 'n_total', 'n_TD', 'n_DD', 
                                                  'mean_TD', 'mean_DD', 'std_TD', 'std_DD',
                                                  'cohens_d', 'p_value', 'p_fdr', 'direction']].copy()
    sig_degree_child_output = sig_degree_child_output.sort_values('p_fdr')
    out_sig_degree_child = os.path.join(out_dir, "Significant_Child_degree_results.csv")
    sig_degree_child_output.to_csv(out_sig_degree_child, index=False)
    print(f"Significant Child degree results (n={len(sig_degree_child)}): {out_sig_degree_child}")

# Degree 显著结果 - Adult
sig_degree_adult = degree_adult_res[degree_adult_res['significant']].copy()
if len(sig_degree_adult) > 0:
    sig_degree_adult_output = sig_degree_adult[['feature', 'n_total', 'n_TD', 'n_DD',
                                                  'mean_TD', 'mean_DD', 'std_TD', 'std_DD',
                                                  'cohens_d', 'p_value', 'p_fdr', 'direction']].copy()
    sig_degree_adult_output = sig_degree_adult_output.sort_values('p_fdr')
    out_sig_degree_adult = os.path.join(out_dir, "Significant_Adult_degree_results.csv")
    sig_degree_adult_output.to_csv(out_sig_degree_adult, index=False)
    print(f"Significant Adult degree results (n={len(sig_degree_adult)}): {out_sig_degree_adult}")

# Edge 显著结果 - Child
sig_edge_child = edge_child_res[edge_child_res['significant']].copy()
if len(sig_edge_child) > 0:
    sig_edge_child_output = sig_edge_child[['feature', 'roi_i', 'roi_j', 'n_total', 'n_TD', 'n_DD',
                                              'mean_TD', 'mean_DD', 'std_TD', 'std_DD',
                                              'cohens_d', 'p_value', 'p_fdr', 'direction']].copy()
    sig_edge_child_output = sig_edge_child_output.sort_values('p_fdr')
    out_sig_edge_child = os.path.join(out_dir, "Significant_Child_edge_results.csv")
    sig_edge_child_output.to_csv(out_sig_edge_child, index=False)
    print(f"Significant Child edge results (n={len(sig_edge_child)}): {out_sig_edge_child}")

# Edge 显著结果 - Adult
sig_edge_adult = edge_adult_res[edge_adult_res['significant']].copy()
if len(sig_edge_adult) > 0:
    sig_edge_adult_output = sig_edge_adult[['feature', 'roi_i', 'roi_j', 'n_total', 'n_TD', 'n_DD',
                                              'mean_TD', 'mean_DD', 'std_TD', 'std_DD',
                                              'cohens_d', 'p_value', 'p_fdr', 'direction']].copy()
    sig_edge_adult_output = sig_edge_adult_output.sort_values('p_fdr')
    out_sig_edge_adult = os.path.join(out_dir, "Significant_Adult_edge_results.csv")
    sig_edge_adult_output.to_csv(out_sig_edge_adult, index=False)
    print(f"Significant Adult edge results (n={len(sig_edge_adult)}): {out_sig_edge_adult}")

sig_degree_child_uncorr = degree_child_res[(degree_child_res['p_value'] < 0.05) & degree_child_res['p_value'].notna()].copy()
if len(sig_degree_child_uncorr) > 0:
    sig_degree_child_uncorr_output = sig_degree_child_uncorr[['feature', 'n_total', 'n_TD', 'n_DD', 
                                                                'mean_TD', 'mean_DD', 'std_TD', 'std_DD',
                                                                'cohens_d', 'p_value', 'p_fdr', 'direction']].copy()
    sig_degree_child_uncorr_output = sig_degree_child_uncorr_output.sort_values('p_value')
    out_sig_degree_child_uncorr = os.path.join(out_dir, "Significant_Child_degree_results_uncorrected.csv")
    sig_degree_child_uncorr_output.to_csv(out_sig_degree_child_uncorr, index=False)
    print(f"Significant Child degree results (uncorrected, n={len(sig_degree_child_uncorr)}): {out_sig_degree_child_uncorr}")

sig_degree_adult_uncorr = degree_adult_res[(degree_adult_res['p_value'] < 0.05) & degree_adult_res['p_value'].notna()].copy()
if len(sig_degree_adult_uncorr) > 0:
    sig_degree_adult_uncorr_output = sig_degree_adult_uncorr[['feature', 'n_total', 'n_TD', 'n_DD',
                                                                'mean_TD', 'mean_DD', 'std_TD', 'std_DD',
                                                                'cohens_d', 'p_value', 'p_fdr', 'direction']].copy()
    sig_degree_adult_uncorr_output = sig_degree_adult_uncorr_output.sort_values('p_value')
    out_sig_degree_adult_uncorr = os.path.join(out_dir, "Significant_Adult_degree_results_uncorrected.csv")
    sig_degree_adult_uncorr_output.to_csv(out_sig_degree_adult_uncorr, index=False)
    print(f"Significant Adult degree results (uncorrected, n={len(sig_degree_adult_uncorr)}): {out_sig_degree_adult_uncorr}")

sig_edge_child_uncorr = edge_child_res[(edge_child_res['p_value'] < 0.05) & edge_child_res['p_value'].notna()].copy()
if len(sig_edge_child_uncorr) > 0:
    sig_edge_child_uncorr_output = sig_edge_child_uncorr[['feature', 'roi_i', 'roi_j', 'n_total', 'n_TD', 'n_DD',
                                                            'mean_TD', 'mean_DD', 'std_TD', 'std_DD',
                                                            'cohens_d', 'p_value', 'p_fdr', 'direction']].copy()
    sig_edge_child_uncorr_output = sig_edge_child_uncorr_output.sort_values('p_value')
    out_sig_edge_child_uncorr = os.path.join(out_dir, "Significant_Child_edge_results_uncorrected.csv")
    sig_edge_child_uncorr_output.to_csv(out_sig_edge_child_uncorr, index=False)
    print(f"Significant Child edge results (uncorrected, n={len(sig_edge_child_uncorr)}): {out_sig_edge_child_uncorr}")

sig_edge_adult_uncorr = edge_adult_res[(edge_adult_res['p_value'] < 0.05) & edge_adult_res['p_value'].notna()].copy()
if len(sig_edge_adult_uncorr) > 0:
    sig_edge_adult_uncorr_output = sig_edge_adult_uncorr[['feature', 'roi_i', 'roi_j', 'n_total', 'n_TD', 'n_DD',
                                                            'mean_TD', 'mean_DD', 'std_TD', 'std_DD',
                                                            'cohens_d', 'p_value', 'p_fdr', 'direction']].copy()
    sig_edge_adult_uncorr_output = sig_edge_adult_uncorr_output.sort_values('p_value')
    out_sig_edge_adult_uncorr = os.path.join(out_dir, "Significant_Adult_edge_results_uncorrected.csv")
    sig_edge_adult_uncorr_output.to_csv(out_sig_edge_adult_uncorr, index=False)
    print(f"Significant Adult edge results (uncorrected, n={len(sig_edge_adult_uncorr)}): {out_sig_edge_adult_uncorr}")

print("\nANALYSIS COMPLETE!")
