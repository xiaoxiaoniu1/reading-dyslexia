import pandas as pd
import numpy as np
from nilearn import plotting
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os

# ==========================================
# 1. 读取坐标文件
# ==========================================
coords_file = '/data/home/tqi/data1/share/after_freesurfer/FILE/DK-68/DK68_MNI_coordinates.csv'
coords_df = pd.read_csv(coords_file)
node_coords = coords_df[['x', 'y', 'z']].values
roi_names = coords_df['label'].values  # 保存脑区名称

n_nodes = len(node_coords)
print(f"加载了 {n_nodes} 个脑区坐标")

# ==========================================
# 2. 读取 ANOVA 结果并转换为矩阵
# ==========================================
anova_file = '/data/home/tqi/data1/share/after_freesurfer/FILE/MIND_DK68_ANOVA/ANOVA_DK68_edge_results.csv'
edge_df = pd.read_csv(anova_file)

# 读取 combat 后的矩阵数据用于计算效应方向 (TD vs DD)
demo_file = '/data/home/tqi/data1/share/after_freesurfer/FILE/all_data_cqt.xlsx'
demo_df = pd.read_excel(demo_file)
mind_dir = '/data/home/tqi/data1/share/after_freesurfer/fs_subjects_all/MIND_DK68_combat/'

# ==========================================
# 3. 计算每条边的效应方向
# ==========================================
def compute_edge_effect_direction(demo_df, mind_dir, n_nodes, group_col, group_a, group_b):
    demo_df = demo_df.copy()
    demo_df['subj_id'] = demo_df['id'].astype(int).astype(str)
    demo_df['mind_file'] = demo_df['subj_id'].apply(lambda x: os.path.join(mind_dir, f"{x}_combat.csv"))

    demo_df['exists'] = demo_df['mind_file'].apply(lambda x: os.path.exists(x))
    valid_df = demo_df[demo_df['exists']].copy()
    print(f"找到 {len(valid_df)} 个有效的矩阵文件用于计算效应方向")

    all_matrices = []
    groups = []

    for _, row in valid_df.iterrows():
        try:
            mat = pd.read_csv(row['mind_file'], index_col=0).values
            if mat.shape == (n_nodes, n_nodes):
                all_matrices.append(mat)
                groups.append(row[group_col])
        except Exception as e:
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
# 4. 配置绘图参数
# ==========================================
p_threshold = 0.05 
top_n_edges = 20  # 只画显著性最高的前20条边
target_p_cols = ['p_Diagnosis_FDR', 'p_AgeGroup_FDR', 'p_Interaction']

# 输出目录
output_dir = '/data/home/tqi/data1/share/after_freesurfer/FILE/MIND_DK68_ANOVA/'
os.makedirs(output_dir, exist_ok=True)

# ==========================================
# 5. 对每一列分别画图
# ==========================================
for target_p_col in target_p_cols:
    print(f"\n处理: {target_p_col}")
    
    # 筛选显著边
    significant_edges = edge_df[edge_df[target_p_col] < p_threshold].copy()
    if len(significant_edges) == 0:
        print(f"⚠️ 提示：{target_p_col} 无显著结果")
        continue
    
    significant_edges = significant_edges.sort_values(by=target_p_col).head(top_n_edges)
    
    if 'AgeGroup' in target_p_col:
        effect_direction_matrix = effect_direction_age
        pos_label = 'Adult'
        neg_label = 'Child'
    else:
        effect_direction_matrix = effect_direction_diagnosis
        pos_label = 'TD'
        neg_label = 'DD'

    matrix = np.zeros((n_nodes, n_nodes))
    for _, row in significant_edges.iterrows():
        i, j = int(row['i']) - 1, int(row['j']) - 1
        
        p_val = row[target_p_col]
        effect_sign = np.sign(effect_direction_matrix[i, j])
        strength = -np.log10(p_val) * effect_sign
        
        matrix[i, j] = matrix[j, i] = strength
    
    # 保存显著边表格
    edge_table = significant_edges.copy()
    edge_table['ROI_i_name'] = edge_table['i'].apply(lambda x: roi_names[int(x)-1])
    edge_table['ROI_j_name'] = edge_table['j'].apply(lambda x: roi_names[int(x)-1])
    edge_table['stronger_in'] = edge_table.apply(
        lambda r: pos_label if effect_direction_matrix[int(r['i'])-1, int(r['j'])-1] > 0 else neg_label,
        axis=1
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
    
    edge_table_path = f"{output_dir}ANOVA_DK68_{target_p_col}_top{top_n_edges}_edges.csv"
    edge_table[['ROI_i_name', 'ROI_j_name', 'i', 'j', target_p_col, eta2_col, 'stronger_in']].to_csv(edge_table_path, index=False)
    
    # 绘图
    max_abs = max(abs(matrix.min()), abs(matrix.max()), 1e-6)
    
    fig = plotting.plot_connectome(
        adjacency_matrix=matrix,
        node_coords=node_coords,
        display_mode='lyrz',
        node_size=5,
        node_color='black',
        edge_cmap='RdBu_r',
        edge_vmin=-max_abs,
        edge_vmax=max_abs,
        edge_threshold=0,  # 显示所有非零边
        edge_kwargs={'linewidth': 1.5},
        colorbar=True,
        title=f"DK68 {target_p_col} -log10(p) (Red:{pos_label}>{neg_label}, Blue:{neg_label}>{pos_label})"
    )
    
    output_path = f"{output_dir}ANOVA_DK68_{target_p_col}_plot.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ 已保存图片: {output_path}")

print("\n全部完成！")
