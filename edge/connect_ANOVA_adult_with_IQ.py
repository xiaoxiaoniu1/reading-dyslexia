import pandas as pd
import numpy as np
from nilearn import plotting
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ==========================================
# 1. 读取坐标文件
# ==========================================
coords_file = '/data/home/tqi/data1/share/after_freesurfer/FILE/DK-318/DK318_MNI_Coordinates.csv'
coords_df = pd.read_csv(coords_file)
node_coords = coords_df[['x', 'y', 'z']].values
roi_names = coords_df['label'].values  # 保存脑区名称

n_nodes = len(node_coords)
print(f"加载了 {n_nodes} 个脑区坐标")

# ==========================================
# 2. 读取 ANOVA 结果并转换为矩阵
# ==========================================
anova_file = '/data/home/tqi/data1/share/after_freesurfer/FILE/MIND_ANOVA_adult_with_IQ/ANOVA_DK318_edge_results.csv'
edge_df = pd.read_csv(anova_file)

# 读取 combat 后的矩阵数据用于计算效应方向
demo_file = '/data/home/tqi/data1/share/after_freesurfer/FILE/all_data_cqt_adult_with_IQ.xlsx'
demo_df = pd.read_excel(demo_file)
iq_candidates = [c for c in demo_df.columns if 'iq' in c.lower()]
if len(iq_candidates) == 0:
    raise ValueError('未找到 IQ 列，请检查列名。')
if len(iq_candidates) > 1:
    print('Multiple IQ columns detected, using:', iq_candidates[0])
iq_col = iq_candidates[0]
mind_dir = '/data1/tqi/share/after_freesurfer/fs_subjects_all/MIND_out_combat_degree/'

# 初始化一个全零的矩阵 (318 x 318)
# 我们用 -log10(P值) 作为连接强度：P值越小，值越大，连线越亮
matrix = np.zeros((n_nodes, n_nodes))

# 设置显著性阈值 (例如只画 FDR p < 0.05 的线)
# 如果你想画所有线，可以设为 1.0
p_threshold = 0.05 

# 设置边数限制 (二选一，注释掉不用的)
# 方式1：限制边的绝对数量
top_n_edges = 20  # 只画显著性最高的前20条边
#top_n_edges = None  # 取消注释这行来禁用边数限制

# 方式2：限制边的百分比
# top_percent_edges = 20  # 只画显著性最高的前20%的边
top_percent_edges = None  # 取消注释来禁用百分比限制

# 定义要分析的6列
target_p_cols = [
    'p_Diagnosis_FDR',
    "p_Diagnosis",
    'p_IQ_FDR',
    "p_IQ"
]

# 输出目录
output_dir = '/data/home/tqi/data1/share/after_freesurfer/FILE/MIND_ANOVA_adult_with_IQ/'

# ==========================================
# 3. 计算每条边的效应方向（TD vs DD）
# ==========================================
def load_matrices_and_info(demo_df, mind_dir, n_nodes, iq_col):
    demo_df = demo_df.copy()
    demo_df['subj_prefix'] = demo_df['original-project'].astype(str) + '_' + demo_df['id_old'].astype(str) + '_MIND_DK318'
    demo_df['file_base'] = demo_df['subj_prefix'] + '_combat_labeled'
    demo_df['mind_file'] = mind_dir + demo_df['file_base'] + '.csv'
    demo_df['exists'] = demo_df['mind_file'].apply(lambda x: pd.io.common.file_exists(x))
    if 'group_age' in demo_df.columns:
        demo_df = demo_df[demo_df['group_age'] == 1]
    demo_df = demo_df[demo_df['exists']].copy()
    demo_df = demo_df[demo_df[iq_col].notna()].copy()

    all_matrices = []
    diagnoses = []
    iq_values = []

    for _, row in demo_df.iterrows():
        try:
            mat = pd.read_csv(row['mind_file'], index_col=0).values
            if mat.shape == (n_nodes, n_nodes):
                all_matrices.append(mat)
                diagnoses.append(row['group_d_or_c'])
                iq_values.append(row[iq_col])
        except:
            continue

    return np.array(all_matrices), np.array(diagnoses), np.array(iq_values, dtype=float)


def compute_td_dd_diff(all_matrices, diagnoses):
    if all_matrices.size == 0:
        return np.zeros((n_nodes, n_nodes))
    td_mask = (diagnoses == 0)
    dd_mask = (diagnoses == 1)
    if td_mask.sum() == 0 or dd_mask.sum() == 0:
        return np.zeros((n_nodes, n_nodes))
    td_mean = all_matrices[td_mask].mean(axis=0)
    dd_mean = all_matrices[dd_mask].mean(axis=0)
    return td_mean - dd_mean


def compute_iq_corr_matrix(all_matrices, iq_values):
    if all_matrices.size == 0:
        return np.zeros((n_nodes, n_nodes))
    n_subj, n_nodes_local, _ = all_matrices.shape
    upper = np.triu_indices(n_nodes_local, 1)
    X = all_matrices[:, upper[0], upper[1]].astype(float)
    iq = iq_values.astype(float)

    iq_mean = np.nanmean(iq)
    iq_std = np.nanstd(iq, ddof=1)
    if not np.isfinite(iq_std) or iq_std == 0:
        corr = np.zeros(X.shape[1])
    else:
        iq_z = (iq - iq_mean) / iq_std
        X_mean = np.nanmean(X, axis=0)
        X_std = np.nanstd(X, axis=0, ddof=1)
        X_z = (X - X_mean) / X_std
        corr = np.nanmean(X_z * iq_z[:, None], axis=0)
        corr = np.where(np.isfinite(corr), corr, 0)
        corr = np.where(X_std == 0, 0, corr)

    corr_matrix = np.zeros((n_nodes_local, n_nodes_local))
    corr_matrix[upper] = corr
    corr_matrix = corr_matrix + corr_matrix.T
    return corr_matrix

print("正在加载矩阵并计算效应方向...")
all_matrices, diagnoses, iq_values = load_matrices_and_info(demo_df, mind_dir, n_nodes, iq_col)
effect_direction_matrix = compute_td_dd_diff(all_matrices, diagnoses)
iq_corr_matrix = compute_iq_corr_matrix(all_matrices, iq_values)
print("效应方向计算完成！")

# ==========================================
# 4. 对每一列分别画图
# ==========================================
for target_p_col in target_p_cols:
    print(f"\n{'='*50}")
    print(f"正在处理列: {target_p_col} (阈值 p < {p_threshold})")
    print(f"{'='*50}")
    
    # 筛选显著的边
    if target_p_col == "p_Diagnosis":
        diag_sig_csv = f"{output_dir}Significant_Diagnosis_DK318_edge_results.csv"
        if pd.io.common.file_exists(diag_sig_csv):
            print(f"Using significant diagnosis table: {diag_sig_csv}")
            significant_edges = pd.read_csv(diag_sig_csv)
            if "sig_raw" in significant_edges.columns:
                significant_edges = significant_edges[significant_edges["sig_raw"] == True]
        else:
            significant_edges = edge_df[edge_df[target_p_col] < p_threshold].copy()
    else:
        significant_edges = edge_df[edge_df[target_p_col] < p_threshold].copy()
    
    if len(significant_edges) == 0:
        print(f"⚠️ 提示：{target_p_col} 在阈值 p < {p_threshold} 下没有发现显著连接")
        continue
    
    # 按p值排序（升序，越小越显著）
    significant_edges = significant_edges.sort_values(by=target_p_col)
    
    # 应用边数限制
    if top_n_edges is not None:
        significant_edges = significant_edges.head(top_n_edges)
        print(f"应用边数限制：保留显著性最高的前 {top_n_edges} 条边")
    elif top_percent_edges is not None:
        n_keep = int(len(significant_edges) * top_percent_edges / 100)
        significant_edges = significant_edges.head(n_keep)
        print(f"应用百分比限制：保留显著性最高的前 {top_percent_edges}% ({n_keep} 条边)")
    
    # 初始化矩阵
    matrix = np.zeros((n_nodes, n_nodes))
    
    count = 0
    for index, row in significant_edges.iterrows():
        # 提取节点索引 (注意：CSV里通常是1-based，Python是0-based，所以要减1)
        i = int(row['i']) - 1
        j = int(row['j']) - 1
        
        p_val = row[target_p_col]
        if target_p_col == "p_Diagnosis" and "direction" in row.index:
            dir_map = {"TD>DD": 1, "TD<DD": -1, "DD>TD": -1, "DD<TD": 1, "equal": 1}
            effect_sign = dir_map.get(str(row["direction"]), 1)
        elif "IQ" in target_p_col:
            effect_sign = np.sign(iq_corr_matrix[i, j])
        else:
            effect_sign = np.sign(effect_direction_matrix[i, j])
        strength = -np.log10(p_val) * effect_sign
        
        # 矩阵是对称的，所以 i,j 和 j,i 都要填
        matrix[i, j] = strength
        matrix[j, i] = strength
        count += 1
    
    print(f"矩阵构建完成！共有 {count} 条显著连接。")
    print(f"  TD更强（正值）: {np.sum(matrix > 0) // 2} 条")
    print(f"  DD更强（负值）: {np.sum(matrix < 0) // 2} 条")
    
    # 保存显著连接的边表格（包含脑区名称和效应方向）
    edge_table = significant_edges.copy()
    edge_table['ROI_i_name'] = (edge_table['i'] - 1).apply(lambda x: roi_names[int(x)])
    edge_table['ROI_j_name'] = (edge_table['j'] - 1).apply(lambda x: roi_names[int(x)])
    
    # 添加效应方向列
    if "IQ" in target_p_col:
        edge_table['effect_direction'] = edge_table.apply(
            lambda row: iq_corr_matrix[int(row['i'])-1, int(row['j'])-1], axis=1
        )
        edge_table['stronger_in'] = edge_table['effect_direction'].apply(
            lambda x: 'IQ+' if x > 0 else ('IQ-' if x < 0 else 'Equal')
        )
    elif target_p_col == "p_Diagnosis" and "direction" in edge_table.columns:
        edge_table['effect_direction'] = edge_table['direction']
        edge_table['stronger_in'] = edge_table['direction'].apply(
            lambda x: 'TD' if str(x) in ('TD>DD', 'DD<TD') else ('DD' if str(x) in ('DD>TD', 'TD<DD') else 'Equal')
        )
    else:
        edge_table['effect_direction'] = edge_table.apply(
            lambda row: effect_direction_matrix[int(row['i'])-1, int(row['j'])-1], axis=1
        )
        edge_table['stronger_in'] = edge_table['effect_direction'].apply(
            lambda x: 'TD' if x > 0 else ('DD' if x < 0 else 'Equal')
        )
    
    # 获取对应的效应量列
    if 'Diagnosis' in target_p_col:
        eta2_col = 'eta2_Diagnosis'
    elif 'IQ' in target_p_col:
        eta2_col = 'eta2_IQ'
    else:
        eta2_col = 'eta2_Diagnosis'
    
    # 选择输出的列
    output_columns = ['ROI_i_name', 'ROI_j_name', 'i', 'j', target_p_col, eta2_col, 'effect_direction', 'stronger_in']
    edge_table_output = edge_table[output_columns]
    
    # 保存为CSV文件
    edge_table_path = f"{output_dir}ANOVA_{target_p_col}_top{top_n_edges}_edges.csv"
    edge_table_output.to_csv(edge_table_path, index=False)
    print(f"✓ 显著连接表格已保存至: {edge_table_path}")
    
    # 计算颜色范围（对称）
    max_abs_val = max(abs(matrix.min()), abs(matrix.max()))
    
    # 画图
    fig = plotting.plot_connectome(
        adjacency_matrix=matrix,    # 带符号的矩阵：正=TD强，负=DD强
        node_coords=node_coords,    # 坐标
        
        # 视图设置
        display_mode='lyrz',        # 左、后、右、顶 四视图
        node_size=1,                # 节点尺寸
        node_color='black',         # 节点颜色
        
        # 连线设置 - 使用红蓝双色渐变
        edge_threshold='10%',       # 显示所有非零值
        edge_kwargs={'linewidth': 1.5},
        edge_cmap='RdBu_r',         # 红蓝配色：红=正值(TD强)，蓝=负值(DD强)
        edge_vmin=-max_abs_val,     # 对称的负值范围
        edge_vmax=max_abs_val,      # 对称的正值范围
        colorbar=True,              # 显示颜色条
        
        title=(
            f"Significant Edges: {target_p_col} -log10(p) "
            f"(Red=IQ+, Blue=IQ-)" if "IQ" in target_p_col else
            f"Significant Edges: {target_p_col} -log10(p) (Red=TD>DD, Blue=DD>TD)"
        )
    )
    
    # 保存图片
    output_path = f"{output_dir}ANOVA_{target_p_col}_plot.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ 图片已保存至: {output_path}")
    
    plt.close()

print(f"\n{'='*50}")
print("所有图片生成完成！")
print(f"{'='*50}")