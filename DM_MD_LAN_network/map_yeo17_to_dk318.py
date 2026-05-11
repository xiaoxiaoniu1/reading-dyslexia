#!/usr/bin/env python3
"""
Map Yeo 17 Networks to DK318 parcellation
将Yeo17网络标签映射到DK318脑区上

使用方法 (Usage):
    python3 map_yeo17_to_dk318.py
    
    或指定参数:
    python3 map_yeo17_to_dk318.py --freesurfer_dir /path/to/freesurfer --dk318_dir /path/to/dk318

输出文件 (Output files):
    - lh_DK318_to_Yeo17_mapping.csv: 左半球简化映射
    - lh_DK318_to_Yeo17_mapping_detailed.csv: 左半球详细映射
    - rh_DK318_to_Yeo17_mapping.csv: 右半球简化映射
    - rh_DK318_to_Yeo17_mapping_detailed.csv: 右半球详细映射

依赖 (Dependencies):
    - nibabel
    - numpy
    - pandas
"""

import os
import sys
import numpy as np
import pandas as pd
import nibabel as nib
from collections import Counter
import argparse


def get_yeo17_network_names():
    """
    获取Yeo17网络的标准名称映射
    
    Returns:
    --------
    dict: 从annot文件中的标签名到标准网络名的映射
    
    Notes:
    ------
    Yeo17网络是Yeo7网络的细分版本，提供更精细的功能网络划分
    参考: Yeo et al. (2011) J Neurophysiol
    根据官方文档的编号对应关系
    """
    # Yeo17网络的标准名称（根据官方对应表）
    yeo17_name_mapping = {
        '17Networks_1': 'VisCent',           # Visual Central / Visual A
        '17Networks_2': 'VisPeri',           # Visual Peripheral / Visual B
        '17Networks_3': 'SomMotA',           # Somatomotor A
        '17Networks_4': 'SomMotB',           # Somatomotor B
        '17Networks_5': 'DorsAttnA',         # Dorsal Attention A
        '17Networks_6': 'DorsAttnB',         # Dorsal Attention B
        '17Networks_7': 'SalVentAttnA',      # Salience/Ventral Attention A
        '17Networks_8': 'SalVentAttnB',      # Salience/Ventral Attention B
        '17Networks_9': 'LimbicA',           # Limbic A
        '17Networks_10': 'LimbicB',          # Limbic B
        '17Networks_11': 'ContC',            # Control C
        '17Networks_12': 'ContA',            # Control A
        '17Networks_13': 'ContB',            # Control B
        '17Networks_14': 'TempPar',          # Temporal Parietal
        '17Networks_15': 'DefaultC',         # Default C
        '17Networks_16': 'DefaultA',         # Default A
        '17Networks_17': 'DefaultB',         # Default B
        'FreeSurfer_Defined_Medial_Wall': 'Unknown'
    }
    return yeo17_name_mapping


def get_yeo17_to_yeo7_mapping():
    """
    获取Yeo17到Yeo7的对应关系
    
    Returns:
    --------
    dict: Yeo17网络名到Yeo7网络名的映射
    
    Notes:
    ------
    根据官方对应表：
    Yeo7-1 (Visual) → Yeo17-1,2
    Yeo7-2 (Somatomotor) → Yeo17-3,4
    Yeo7-3 (Dorsal Attention) → Yeo17-5,6
    Yeo7-4 (Salience/Ventral Attention) → Yeo17-7,8
    Yeo7-5 (Limbic) → Yeo17-9,10
    Yeo7-6 (Control) → Yeo17-11,12,13
    Yeo7-7 (Default) → Yeo17-15,16,17
    Yeo17-14 (Temporal Parietal) 是独立的
    """
    yeo17_to_yeo7 = {
        'VisCent': 'Visual',
        'VisPeri': 'Visual',
        'SomMotA': 'Somatomotor',
        'SomMotB': 'Somatomotor',
        'DorsAttnA': 'Dorsal Attention',
        'DorsAttnB': 'Dorsal Attention',
        'SalVentAttnA': 'Salience/Ventral Attention',
        'SalVentAttnB': 'Salience/Ventral Attention',
        'LimbicA': 'Limbic',
        'LimbicB': 'Limbic',
        'ContC': 'Control',
        'ContA': 'Control',
        'ContB': 'Control',
        'TempPar': 'Temporal Parietal',
        'DefaultC': 'Default',
        'DefaultA': 'Default',
        'DefaultB': 'Default',
        'Unknown': 'Unknown'
    }
    return yeo17_to_yeo7


def read_annot(annot_file):
    """读取FreeSurfer annot文件"""
    try:
        labels, ctab, names = nib.freesurfer.read_annot(annot_file)
        return labels, ctab, names
    except Exception as e:
        print(f"Error reading {annot_file}: {e}")
        sys.exit(1)


def map_yeo17_to_dk318(yeo17_annot, dk318_annot, surf_file=None):
    """
    将Yeo17网络映射到DK318脑区
    
    Parameters:
    -----------
    yeo17_annot : str
        Yeo17 annot文件路径
    dk318_annot : str
        DK318 annot文件路径
    surf_file : str, optional
        表面文件路径 (用于获取顶点数)
    
    Returns:
    --------
    mapping_df : pd.DataFrame
        映射结果的DataFrame
    """
    print(f"Reading Yeo17 annot: {yeo17_annot}")
    yeo17_labels, yeo17_ctab, yeo17_names = read_annot(yeo17_annot)
    yeo17_names = [name.decode('utf-8') if isinstance(name, bytes) else name 
                   for name in yeo17_names]
    
    # 获取Yeo17网络名称映射
    yeo17_name_mapping = get_yeo17_network_names()
    yeo17_to_yeo7 = get_yeo17_to_yeo7_mapping()
    
    # 将annot中的标签名映射到标准网络名
    yeo17_standard_names = []
    yeo7_parent_names = []
    for name in yeo17_names:
        standard_name = yeo17_name_mapping.get(name, name)
        yeo17_standard_names.append(standard_name)
        parent_name = yeo17_to_yeo7.get(standard_name, 'Unknown')
        yeo7_parent_names.append(parent_name)
    
    print(f"Reading DK318 annot: {dk318_annot}")
    dk318_labels, dk318_ctab, dk318_names = read_annot(dk318_annot)
    dk318_names = [name.decode('utf-8') if isinstance(name, bytes) else name 
                   for name in dk318_names]
    
    print(f"\nYeo17 networks found: {len(yeo17_names)}")
    print(f"Original Yeo17 labels: {yeo17_names}")
    print(f"Standard Yeo17 names: {yeo17_standard_names}")
    print(f"\nDK318 regions found: {len(dk318_names)}")
    print(f"Total vertices: {len(yeo17_labels)}")
    
    # 检查顶点数是否匹配
    if len(yeo17_labels) != len(dk318_labels):
        print(f"\nWarning: Vertex count mismatch!")
        print(f"  Yeo17: {len(yeo17_labels)} vertices")
        print(f"  DK318: {len(dk318_labels)} vertices")
        print(f"  Using minimum: {min(len(yeo17_labels), len(dk318_labels))}")
        min_len = min(len(yeo17_labels), len(dk318_labels))
        yeo17_labels = yeo17_labels[:min_len]
        dk318_labels = dk318_labels[:min_len]
    
    # 创建映射字典
    mapping_results = []
    
    # 遍历每个DK318脑区
    print(f"\nMapping DK318 regions to Yeo17 networks...")
    for dk_idx, dk_name in enumerate(dk318_names):
        # 找到属于当前DK318脑区的所有顶点
        dk_vertices = np.where(dk318_labels == dk_idx)[0]
        
        if len(dk_vertices) == 0:
            continue
        
        # 获取这些顶点对应的Yeo17标签
        yeo17_labels_in_dk = yeo17_labels[dk_vertices]
        
        # 统计每个Yeo17网络的顶点数
        yeo17_counter = Counter(yeo17_labels_in_dk)
        
        # 找到最多的Yeo17网络
        if len(yeo17_counter) > 0:
            most_common_yeo17_idx, most_common_count = yeo17_counter.most_common(1)[0]
            
            # 获取Yeo17网络的原始名称和标准名称
            if most_common_yeo17_idx < len(yeo17_names):
                most_common_yeo17_original = yeo17_names[most_common_yeo17_idx]
                most_common_yeo17_name = yeo17_standard_names[most_common_yeo17_idx]
                most_common_yeo7_parent = yeo7_parent_names[most_common_yeo17_idx]
            else:
                most_common_yeo17_original = f"Unknown_{most_common_yeo17_idx}"
                most_common_yeo17_name = f"Unknown_{most_common_yeo17_idx}"
                most_common_yeo7_parent = "Unknown"
            
            # 计算占比
            percentage = (most_common_count / len(dk_vertices)) * 100
            
            # 获取所有Yeo17网络的分布
            yeo17_distribution = {}
            for yeo_idx, count in yeo17_counter.items():
                if yeo_idx < len(yeo17_names):
                    yeo_original = yeo17_names[yeo_idx]
                    yeo_standard = yeo17_standard_names[yeo_idx]
                    yeo7_parent = yeo7_parent_names[yeo_idx]
                else:
                    yeo_original = f"Unknown_{yeo_idx}"
                    yeo_standard = f"Unknown_{yeo_idx}"
                    yeo7_parent = "Unknown"
                yeo17_distribution[yeo_standard] = {
                    'count': count,
                    'percentage': (count / len(dk_vertices)) * 100,
                    'original_label': yeo_original,
                    'yeo7_parent': yeo7_parent
                }
            
            mapping_results.append({
                'DK318_index': dk_idx,
                'DK318_region': dk_name,
                'num_vertices': len(dk_vertices),
                'primary_Yeo17_network': most_common_yeo17_name,
                'primary_Yeo17_original_label': most_common_yeo17_original,
                'primary_Yeo17_index': most_common_yeo17_idx,
                'primary_Yeo7_parent': most_common_yeo7_parent,
                'primary_count': most_common_count,
                'primary_percentage': percentage,
                'yeo17_distribution': yeo17_distribution
            })
    
    # 创建DataFrame
    mapping_df = pd.DataFrame(mapping_results)
    print(f"Successfully mapped {len(mapping_df)} DK318 regions")
    
    return mapping_df


def save_detailed_mapping(mapping_df, output_file):
    """保存详细的映射结果"""
    detailed_rows = []
    
    for _, row in mapping_df.iterrows():
        dk_region = row['DK318_region']
        num_vertices = row['num_vertices']
        
        # 添加所有网络信息（按百分比排序）
        yeo17_dist_sorted = sorted(row['yeo17_distribution'].items(), 
                                   key=lambda x: x[1]['percentage'], 
                                   reverse=True)
        
        for yeo_name, info in yeo17_dist_sorted:
            is_primary = (yeo_name == row['primary_Yeo17_network'])
            detailed_rows.append({
                'DK318_region': dk_region,
                'DK318_index': row['DK318_index'],
                'total_vertices': num_vertices,
                'Yeo17_network': yeo_name,
                'Yeo7_parent': info['yeo7_parent'],
                'vertex_count': info['count'],
                'percentage': round(info['percentage'], 2),
                'is_primary': is_primary
            })
    
    detailed_df = pd.DataFrame(detailed_rows)
    detailed_df.to_csv(output_file, index=False)
    print(f"  Detailed mapping saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Map Yeo 17 Networks to DK318 parcellation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 (Examples):
  python3 map_yeo17_to_dk318.py
  python3 map_yeo17_to_dk318.py --freesurfer_dir /usr/local/freesurfer
        """
    )
    parser.add_argument('--freesurfer_dir', type=str,
                       default='/data/software/freesurfer',
                       help='FreeSurfer installation directory')
    parser.add_argument('--dk318_dir', type=str,
                       default='/data/home/tqi/data1/share/after_freesurfer/FILE/DK-318',
                       help='DK318 annot files directory')
    parser.add_argument('--output_dir', type=str,
                       default='/data/home/tqi/data1/share/after_freesurfer/FILE/DK-318',
                       help='Output directory')
    
    args = parser.parse_args()
    
    print("="*70)
    print("Yeo17 to DK318 Mapping Tool")
    print("="*70)
    
    # 设置路径
    fsaverage_label_dir = os.path.join(args.freesurfer_dir, 
                                       'subjects/fsaverage/label')
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 处理左右半球
    for hemi in ['lh', 'rh']:
        print(f"\n{'='*70}")
        print(f"Processing {hemi.upper()} hemisphere")
        print(f"{'='*70}")
        
        # 文件路径
        yeo17_annot = os.path.join(fsaverage_label_dir, 
                                   f'{hemi}.Yeo2011_17Networks_N1000.annot')
        dk318_annot = os.path.join(args.dk318_dir, f'{hemi}.DK318.annot')
        
        # 检查文件是否存在
        if not os.path.exists(yeo17_annot):
            print(f"ERROR: Yeo17 annot file not found: {yeo17_annot}")
            continue
        if not os.path.exists(dk318_annot):
            print(f"ERROR: DK318 annot file not found: {dk318_annot}")
            continue
        
        # 执行映射
        mapping_df = map_yeo17_to_dk318(yeo17_annot, dk318_annot)
        
        # 保存简化版映射结果
        simple_output = os.path.join(args.output_dir, 
                                     f'{hemi}_DK318_to_Yeo17_mapping.csv')
        mapping_df[['DK318_region', 'DK318_index', 'num_vertices', 
                    'primary_Yeo17_network', 'primary_Yeo17_original_label', 
                    'primary_Yeo17_index', 'primary_Yeo7_parent',
                    'primary_count', 'primary_percentage']].to_csv(
            simple_output, index=False
        )
        print(f"\n  Simple mapping saved to: {simple_output}")
        
        # 保存详细版映射结果
        detailed_output = os.path.join(args.output_dir, 
                                       f'{hemi}_DK318_to_Yeo17_mapping_detailed.csv')
        save_detailed_mapping(mapping_df, detailed_output)
        
        # 打印统计信息
        print(f"\n  Summary for {hemi}:")
        print(f"  - Total DK318 regions: {len(mapping_df)}")
        print(f"\n  Yeo17 network distribution (number of DK318 regions):")
        yeo17_counts = mapping_df['primary_Yeo17_network'].value_counts()
        for network, count in yeo17_counts.items():
            print(f"    {network}: {count} regions")
        
        print(f"\n  Yeo7 parent network distribution (number of DK318 regions):")
        yeo7_counts = mapping_df['primary_Yeo7_parent'].value_counts()
        for network, count in yeo7_counts.items():
            print(f"    {network}: {count} regions")
    
    print(f"\n{'='*70}")
    print("Mapping completed successfully!")
    print(f"Output directory: {args.output_dir}")
    print(f"{'='*70}\n")


if __name__ == '__main__':
    main()
