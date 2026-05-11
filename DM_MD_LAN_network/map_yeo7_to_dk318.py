#!/usr/bin/env python3
"""
Map Yeo 7 Networks to DK318 parcellation
将Yeo7网络标签映射到DK318脑区上

使用方法 (Usage):
    python3 map_yeo7_to_dk318.py
    
    或指定参数:
    python3 map_yeo7_to_dk318.py --freesurfer_dir /path/to/freesurfer --dk318_dir /path/to/dk318

输出文件 (Output files):
    - lh_DK318_to_Yeo7_mapping.csv: 左半球简化映射
    - lh_DK318_to_Yeo7_mapping_detailed.csv: 左半球详细映射
    - rh_DK318_to_Yeo7_mapping.csv: 右半球简化映射
    - rh_DK318_to_Yeo7_mapping_detailed.csv: 右半球详细映射

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


def get_yeo7_network_names():
    """
    获取Yeo7网络的标准名称映射
    
    Returns:
    --------
    dict: 从annot文件中的标签名到标准网络名的映射
    
    Notes:
    ------
    根据Yeo et al. 2011的定义和官方对应表
    """
    # Yeo7网络的标准名称
    yeo7_name_mapping = {
        '7Networks_1': 'Visual',
        '7Networks_2': 'Somatomotor', 
        '7Networks_3': 'Dorsal Attention',
        '7Networks_4': 'Salience/Ventral Attention',
        '7Networks_5': 'Limbic',
        '7Networks_6': 'Control',
        '7Networks_7': 'Default',
        'FreeSurfer_Defined_Medial_Wall': 'Unknown'
    }
    return yeo7_name_mapping


def read_annot(annot_file):
    """读取FreeSurfer annot文件"""
    try:
        labels, ctab, names = nib.freesurfer.read_annot(annot_file)
        return labels, ctab, names
    except Exception as e:
        print(f"Error reading {annot_file}: {e}")
        sys.exit(1)


def map_yeo7_to_dk318(yeo7_annot, dk318_annot, surf_file=None):
    """
    将Yeo7网络映射到DK318脑区
    
    Parameters:
    -----------
    yeo7_annot : str
        Yeo7 annot文件路径
    dk318_annot : str
        DK318 annot文件路径
    surf_file : str, optional
        表面文件路径 (用于获取顶点数)
    
    Returns:
    --------
    mapping_df : pd.DataFrame
        映射结果的DataFrame
    """
    print(f"Reading Yeo7 annot: {yeo7_annot}")
    yeo7_labels, yeo7_ctab, yeo7_names = read_annot(yeo7_annot)
    yeo7_names = [name.decode('utf-8') if isinstance(name, bytes) else name 
                  for name in yeo7_names]
    
    # 获取Yeo7网络名称映射
    yeo7_name_mapping = get_yeo7_network_names()
    
    # 将annot中的标签名映射到标准网络名
    yeo7_standard_names = []
    for name in yeo7_names:
        standard_name = yeo7_name_mapping.get(name, name)
        yeo7_standard_names.append(standard_name)
    
    print(f"Reading DK318 annot: {dk318_annot}")
    dk318_labels, dk318_ctab, dk318_names = read_annot(dk318_annot)
    dk318_names = [name.decode('utf-8') if isinstance(name, bytes) else name 
                   for name in dk318_names]
    
    print(f"\nYeo7 networks found: {len(yeo7_names)}")
    print(f"Original Yeo7 labels: {yeo7_names}")
    print(f"Standard Yeo7 names: {yeo7_standard_names}")
    print(f"\nDK318 regions found: {len(dk318_names)}")
    print(f"Total vertices: {len(yeo7_labels)}")
    
    # 检查顶点数是否匹配
    if len(yeo7_labels) != len(dk318_labels):
        print(f"\nWarning: Vertex count mismatch!")
        print(f"  Yeo7: {len(yeo7_labels)} vertices")
        print(f"  DK318: {len(dk318_labels)} vertices")
        print(f"  Using minimum: {min(len(yeo7_labels), len(dk318_labels))}")
        min_len = min(len(yeo7_labels), len(dk318_labels))
        yeo7_labels = yeo7_labels[:min_len]
        dk318_labels = dk318_labels[:min_len]
    
    # 创建映射字典
    mapping_results = []
    
    # 遍历每个DK318脑区
    print(f"\nMapping DK318 regions to Yeo7 networks...")
    for dk_idx, dk_name in enumerate(dk318_names):
        # 找到属于当前DK318脑区的所有顶点
        dk_vertices = np.where(dk318_labels == dk_idx)[0]
        
        if len(dk_vertices) == 0:
            continue
        
        # 获取这些顶点对应的Yeo7标签
        yeo7_labels_in_dk = yeo7_labels[dk_vertices]
        
        # 统计每个Yeo7网络的顶点数
        yeo7_counter = Counter(yeo7_labels_in_dk)
        
        # 找到最多的Yeo7网络
        if len(yeo7_counter) > 0:
            most_common_yeo7_idx, most_common_count = yeo7_counter.most_common(1)[0]
            
            # 获取Yeo7网络的原始名称和标准名称
            if most_common_yeo7_idx < len(yeo7_names):
                most_common_yeo7_original = yeo7_names[most_common_yeo7_idx]
                most_common_yeo7_name = yeo7_standard_names[most_common_yeo7_idx]
            else:
                most_common_yeo7_original = f"Unknown_{most_common_yeo7_idx}"
                most_common_yeo7_name = f"Unknown_{most_common_yeo7_idx}"
            
            # 计算占比
            percentage = (most_common_count / len(dk_vertices)) * 100
            
            # 获取所有Yeo7网络的分布
            yeo7_distribution = {}
            for yeo_idx, count in yeo7_counter.items():
                if yeo_idx < len(yeo7_names):
                    yeo_original = yeo7_names[yeo_idx]
                    yeo_standard = yeo7_standard_names[yeo_idx]
                else:
                    yeo_original = f"Unknown_{yeo_idx}"
                    yeo_standard = f"Unknown_{yeo_idx}"
                yeo7_distribution[yeo_standard] = {
                    'count': count,
                    'percentage': (count / len(dk_vertices)) * 100,
                    'original_label': yeo_original
                }
            
            mapping_results.append({
                'DK318_index': dk_idx,
                'DK318_region': dk_name,
                'num_vertices': len(dk_vertices),
                'primary_Yeo7_network': most_common_yeo7_name,
                'primary_Yeo7_original_label': most_common_yeo7_original,
                'primary_Yeo7_index': most_common_yeo7_idx,
                'primary_count': most_common_count,
                'primary_percentage': percentage,
                'yeo7_distribution': yeo7_distribution
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
        yeo7_dist_sorted = sorted(row['yeo7_distribution'].items(), 
                                 key=lambda x: x[1]['percentage'], 
                                 reverse=True)
        
        for yeo_name, info in yeo7_dist_sorted:
            is_primary = (yeo_name == row['primary_Yeo7_network'])
            detailed_rows.append({
                'DK318_region': dk_region,
                'DK318_index': row['DK318_index'],
                'total_vertices': num_vertices,
                'Yeo7_network': yeo_name,
                'vertex_count': info['count'],
                'percentage': round(info['percentage'], 2),
                'is_primary': is_primary
            })
    
    detailed_df = pd.DataFrame(detailed_rows)
    detailed_df.to_csv(output_file, index=False)
    print(f"  Detailed mapping saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Map Yeo 7 Networks to DK318 parcellation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 (Examples):
  python3 map_yeo7_to_dk318.py
  python3 map_yeo7_to_dk318.py --freesurfer_dir /usr/local/freesurfer
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
    print("Yeo7 to DK318 Mapping Tool")
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
        yeo7_annot = os.path.join(fsaverage_label_dir, 
                                  f'{hemi}.Yeo2011_7Networks_N1000.annot')
        dk318_annot = os.path.join(args.dk318_dir, f'{hemi}.DK318.annot')
        
        # 检查文件是否存在
        if not os.path.exists(yeo7_annot):
            print(f"ERROR: Yeo7 annot file not found: {yeo7_annot}")
            continue
        if not os.path.exists(dk318_annot):
            print(f"ERROR: DK318 annot file not found: {dk318_annot}")
            continue
        
        # 执行映射
        mapping_df = map_yeo7_to_dk318(yeo7_annot, dk318_annot)
        
        # 保存简化版映射结果
        simple_output = os.path.join(args.output_dir, 
                                     f'{hemi}_DK318_to_Yeo7_mapping.csv')
        mapping_df[['DK318_region', 'DK318_index', 'num_vertices', 
                    'primary_Yeo7_network', 'primary_Yeo7_original_label', 'primary_Yeo7_index',
                    'primary_count', 'primary_percentage']].to_csv(
            simple_output, index=False
        )
        print(f"\n  Simple mapping saved to: {simple_output}")
        
        # 保存详细版映射结果
        detailed_output = os.path.join(args.output_dir, 
                                       f'{hemi}_DK318_to_Yeo7_mapping_detailed.csv')
        save_detailed_mapping(mapping_df, detailed_output)
        
        # 打印统计信息
        print(f"\n  Summary for {hemi}:")
        print(f"  - Total DK318 regions: {len(mapping_df)}")
        print(f"\n  Yeo7 network distribution (number of DK318 regions):")
        yeo7_counts = mapping_df['primary_Yeo7_network'].value_counts()
        for network, count in yeo7_counts.items():
            print(f"    {network}: {count} regions")
    
    print(f"\n{'='*70}")
    print("Mapping completed successfully!")
    print(f"Output directory: {args.output_dir}")
    print(f"{'='*70}\n")


if __name__ == '__main__':
    main()
