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
anova_file = '/data/home/tqi/data1/share/after_freesurfer/FILE/MIND_ANOVA/ANOVA_DK318_edge_results.csv'
edge_df = pd.read_csv(anova_file)

# 读取 combat 后的矩阵数据用于计算效应方向
demo_file = '/data/home/tqi/data1/share/after_freesurfer/FILE/all_data_cqt.xlsx'
demo_df = pd.read_excel(demo_file)
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
    #'p_Diagnosis',
    'p_Diagnosis_FDR',
    #'p_AgeGroup',
    'p_AgeGroup_FDR',
    'p_Interaction',
    #'p_Interaction_FDR'
]

# 输出目录
output_dir = '/data/home/tqi/data1/share/after_freesurfer/FILE/MIND_ANOVA/'

# ==========================================
# 3. 计算每条边的效应方向
# ==========================================
def compute_edge_effect_direction(demo_df, mind_dir, n_nodes, group_col, group_a, group_b):
    demo_df = demo_df.copy()
    demo_df['subj_prefix'] = demo_df['original-project'].astype(str) + '_' + demo_df['id_old'].astype(str) + '_MIND_DK318'
    demo_df['file_base'] = demo_df['subj_prefix'] + '_combat_labeled'
    demo_df['mind_file'] = mind_dir + demo_df['file_base'] + '.csv'
    demo_df['exists'] = demo_df['mind_file'].apply(lambda x: pd.io.common.file_exists(x))
    demo_df = demo_df[demo_df['exists']].copy()

    all_matrices = []
    groups = []

    for _, row in demo_df.iterrows():
        try:
            mat = pd.read_csv(row['mind_file'], index_col=0).values
            if mat.shape == (n_nodes, n_nodes):
                all_matrices.append(mat)
                groups.append(row[group_col])
        except:
            continue

    all_matrices = np.array(all_matrices)
    groups = np.array(groups)

    a_mean = all_matrices[groups == group_a].mean(axis=0)
    b_mean = all_matrices[groups == group_b].mean(axis=0)

    diff_matrix = a_mean - b_mean

    return diff_matrix

print("正在计算效应方向（TD vs DD / Adult vs Child）...")
effect_direction_diagnosis = compute_edge_effect_direction(demo_df, mind_dir, n_nodes, "group_d_or_c", 0, 1)
effect_direction_age = compute_edge_effect_direction(demo_df, mind_dir, n_nodes, "group_age", 1, 2)
print("效应方向计算完成！")

# ==========================================
# 4. 对每一列分别画图
# ==========================================
for target_p_col in target_p_cols:
    print(f"\n{'='*50}")
    print(f"正在处理列: {target_p_col} (阈值 p < {p_threshold})")
    print(f"{'='*50}")
    
    # 筛选显著的边
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
    
    if 'AgeGroup' in target_p_col:
        effect_direction_matrix = effect_direction_age
        pos_label = 'Adult'
        neg_label = 'Child'
    else:
        effect_direction_matrix = effect_direction_diagnosis
        pos_label = 'TD'
        neg_label = 'DD'

    count = 0
    for index, row in significant_edges.iterrows():
        i = int(row['i']) - 1
        j = int(row['j']) - 1
        
        p_val = row[target_p_col]
        effect_sign = np.sign(effect_direction_matrix[i, j])
        strength = -np.log10(p_val) * effect_sign
        
        matrix[i, j] = strength
        matrix[j, i] = strength
        count += 1
    
    print(f"矩阵构建完成！共有 {count} 条显著连接。")
    print(f"  {pos_label}更强（正值）: {np.sum(matrix > 0) // 2} 条")
    print(f"  {neg_label}更强（负值）: {np.sum(matrix < 0) // 2} 条")
    
    # 保存显著连接的边表格（包含脑区名称和效应方向）
    edge_table = significant_edges.copy()
    edge_table['ROI_i_name'] = (edge_table['i'] - 1).apply(lambda x: roi_names[int(x)])
    edge_table['ROI_j_name'] = (edge_table['j'] - 1).apply(lambda x: roi_names[int(x)])
    
    # 添加效应方向列
    edge_table['effect_direction'] = edge_table.apply(
        lambda row: effect_direction_matrix[int(row['i'])-1, int(row['j'])-1], axis=1
    )
    edge_table['stronger_in'] = edge_table['effect_direction'].apply(
        lambda x: pos_label if x > 0 else neg_label
    )
    
    # 获取对应的效应量列
    if 'Diagnosis' in target_p_col:
        eta2_col = 'eta2_Diagnosis'
    elif 'AgeGroup' in target_p_col:
        eta2_col = 'eta2_AgeGroup'
    elif 'Interaction' in target_p_col:
        eta2_col = 'eta2_Interaction'
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
        edge_cmap='RdBu_r',         # 红蓝配色：红=正值，蓝=负值
        edge_vmin=-max_abs_val,     # 对称的负值范围
        edge_vmax=max_abs_val,      # 对称的正值范围
        colorbar=True,              # 显示颜色条
        
        title=f"Significant Edges: {target_p_col} -log10(p) (Red={pos_label}>{neg_label}, Blue={neg_label}>{pos_label})"
    )
    
    # 保存图片
    output_path = f"{output_dir}ANOVA_{target_p_col}_plot.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ 图片已保存至: {output_path}")
    
    plt.close()

print(f"\n{'='*50}")
print("所有图片生成完成！")
print(f"{'='*50}")