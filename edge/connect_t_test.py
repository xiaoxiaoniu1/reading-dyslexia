import pandas as pd
import numpy as np
from nilearn import plotting
import matplotlib.pyplot as plt

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
# 2. 读取 t-test 结果并转换为矩阵
# ==========================================
ttest_file = '/data/home/tqi/data1/share/after_freesurfer/FILE/MIND_t_test/Stratified_Child_edge_results.csv'
edge_df = pd.read_csv(ttest_file)

# 从文件名提取类别（Adult 或 Child）
import os
file_basename = os.path.basename(ttest_file)  # Stratified_Adult_edge_results.csv
group_name = file_basename.split('_')[1]  # Adult 或 Child

# 设置边数限制 (二选一，注释掉不用的)
# 方式1：限制边的绝对数量
top_n_edges = 20  # 只画显著性最高的前200条边
#top_n_edges = None  # 取消注释这行来禁用边数限制

# 方式2：限制边的百分比
# top_percent_edges = 20  # 只画显著性最高的前20%的边
top_percent_edges = None  # 取消注释来禁用百分比限制

# 定义要分析的2列
target_p_cols = [
    #'p_value',
    'p_fdr'
]

# 输出目录
output_dir = '/data/home/tqi/data1/share/after_freesurfer/FILE/MIND_t_test/'

# ==========================================
# 3. 对每一列分别画图
# ==========================================
for target_p_col in target_p_cols:
    print(f"\n{'='*50}")
    print(f"正在处理列: {target_p_col}")
    print(f"{'='*50}")
    
    # 筛选显著的边（只保留significant=True的行）
    significant_edges = edge_df[edge_df['significant'] == True].copy()
    
    if len(significant_edges) == 0:
        print(f"⚠️ 提示：{target_p_col} 没有发现显著连接")
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
        # 提取节点索引 (t-test文件中是0-based，直接使用)
        i = int(row['roi_i'])
        j = int(row['roi_j'])
        
        p_val = row[target_p_col]
        direction = row['direction']
        if direction == 'DD>TD':
            strength = -abs(-np.log10(p_val))
        else:
            strength = abs(-np.log10(p_val))
        
        # 矩阵是对称的，所以 i,j 和 j,i 都要填
        matrix[i, j] = strength
        matrix[j, i] = strength
        count += 1
    
    print(f"矩阵构建完成！共有 {count} 条显著连接。")
    print(f"  TD更强（正值）: {np.sum(matrix > 0) // 2} 条")
    print(f"  DD更强（负值）: {np.sum(matrix < 0) // 2} 条")
    
    # 保存显著连接的边表格（包含脑区名称）
    edge_table = significant_edges.copy()
    edge_table['ROI_i_name'] = edge_table['roi_i'].apply(lambda x: roi_names[int(x)])
    edge_table['ROI_j_name'] = edge_table['roi_j'].apply(lambda x: roi_names[int(x)])
    
    # 选择输出的列
    output_columns = ['ROI_i_name', 'ROI_j_name', 'roi_i', 'roi_j', 'cohens_d', target_p_col, 'direction']
    edge_table_output = edge_table[output_columns]
    
    # 保存为CSV文件
    edge_table_path = f"{output_dir}ttest_{group_name}_{target_p_col}_significant_edges.csv"
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
        
        title=f"Significant Edges: {target_p_col} -log10(p) (Red=TD>DD, Blue=DD>TD)"
    )
    
    # 保存图片
    output_path = f"{output_dir}ttest_{group_name}_{target_p_col}_plot.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ 图片已保存至: {output_path}")
    
    plt.close()

print(f"\n{'='*50}")
print("所有图片生成完成！")
print(f"{'='*50}")