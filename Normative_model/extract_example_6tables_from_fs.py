#!/usr/bin/env python3
"""
从 FreeSurfer recon-all 结果提取 example 所需的 6 个形态学表：
- cortical-thickness-female.xlsx
- cortical-thickness-male.xlsx
- surface-area-female.xlsx
- surface-area-male.xlsx
- subcortical-volume-female.xlsx
- subcortical-volume-male.xlsx
"""

from pathlib import Path
import re

import pandas as pd

DK68_REGIONS = [
    'bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal', 'cuneus',
    'entorhinal', 'fusiform', 'inferiorparietal', 'inferiortemporal',
    'isthmuscingulate', 'lateraloccipital', 'lateralorbitofrontal',
    'lingual', 'medialorbitofrontal', 'middletemporal', 'parahippocampal',
    'paracentral', 'parsopercularis', 'parsorbitalis', 'parstriangularis',
    'pericalcarine', 'postcentral', 'posteriorcingulate', 'precentral',
    'precuneus', 'rostralanteriorcingulate', 'rostralmiddlefrontal',
    'superiorfrontal', 'superiorparietal', 'superiortemporal',
    'supramarginal', 'frontalpole', 'temporalpole', 'transversetemporal', 'insula'
]

APARC_METRIC_INDEX = {'area': 2, 'volume': 3, 'thickness': 4}
PROJECT_ALIAS = {
    'adulec': 'adolec',
    'aduld': 'adoled',
    'controldataadultscontrol': 'controldataadultscontrol',
    'controldataagecontrol': 'controldataagecontrol'
}

SUBCORTICAL_FEATURES = [
    ('ICV', ['eTIV', 'EstimatedTotalIntraCranialVol']),
    ('Lthal', ['Left-Thalamus', 'Left-Thalamus-Proper']),
    ('Rthal', ['Right-Thalamus', 'Right-Thalamus-Proper']),
    ('Lcaud', ['Left-Caudate']),
    ('Rcaud', ['Right-Caudate']),
    ('Lput', ['Left-Putamen']),
    ('Rput', ['Right-Putamen']),
    ('Lpal', ['Left-Pallidum']),
    ('Rpal', ['Right-Pallidum']),
    ('Lhippo', ['Left-Hippocampus']),
    ('Rhippo', ['Right-Hippocampus']),
    ('Lamyg', ['Left-Amygdala']),
    ('Ramyg', ['Right-Amygdala']),
    ('Laccumb', ['Left-Accumbens-area']),
    ('Raccumb', ['Right-Accumbens-area']),
]

META_COLUMNS = ['SITE', 'SubjectID', 'Vendor', 'FreeSurfer_Version', 'age', 'sex']

REGION_NAME_MAP = {
    'entorhinal': 'entorhil',
    'supramarginal': 'supramargil',
}

BASE_DIR = Path('/data/home/tqi/data1/share/after_freesurfer')
FS_SUBJECTS_DIR = BASE_DIR / 'fs_subjects_all'
CLINICAL_CSV = BASE_DIR / 'FILE/test_mean_1.5/Normative_model/Datasets/Datasets-new/Clinical_vars_merged.csv'
OUTPUT_DIR = BASE_DIR / 'FILE/test_mean_1.5/Normative_model/example_extracted'


def normalize_id_value(value):
    if pd.isna(value):
        return ''
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if value.is_integer():
            return str(int(value))
        return str(value).strip()
    text = str(value).strip()
    if text.endswith('.0') and text[:-2].isdigit():
        return text[:-2]
    return text


def normalize_project(value):
    text = re.sub(r'[^a-z0-9]+', '', str(value).strip().lower())
    return PROJECT_ALIAS.get(text, text)


def normalize_numeric_tail(value):
    text = normalize_id_value(value)
    match = re.search(r'(\d+)$', text)
    if match:
        return str(int(match.group(1)))
    return text.lower()


def build_subject_key(project, sid):
    return f'{normalize_project(project)}_{normalize_numeric_tail(sid)}'


def split_subject_name(subject_name):
    left, sep, right = subject_name.rpartition('_')
    if sep:
        return left, right
    match = re.match(r'^(.*?)(\d+)$', subject_name)
    if match:
        return match.group(1), match.group(2)
    return subject_name, subject_name


def to_float(text):
    return float(str(text).strip())


def parse_measure_line(line):
    content = line[len('# Measure '):].strip()
    parts = [part.strip() for part in content.split(',')]
    if len(parts) < 4:
        return None, None
    field = parts[1]
    value_text = parts[3].split()[0]
    try:
        return field, float(value_text)
    except ValueError:
        return field, None


def parse_aparc_stats(stats_file):
    summary = {}
    regions = {}
    data_start = None
    with open(stats_file, 'r', encoding='utf-8', errors='ignore') as handle:
        lines = handle.readlines()
    for index, line in enumerate(lines):
        if line.startswith('# Measure '):
            field, value = parse_measure_line(line)
            if field and value is not None:
                summary[field] = value
        if line.startswith('# ColHeaders'):
            data_start = index + 1
            break
    if data_start is None:
        raise ValueError(f'无法在 {stats_file} 中找到 ColHeaders')
    for line in lines[data_start:]:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) < 5:
            continue
        region = parts[0]
        if region not in DK68_REGIONS:
            continue
        regions[region] = {
            'area': to_float(parts[APARC_METRIC_INDEX['area']]),
            'volume': to_float(parts[APARC_METRIC_INDEX['volume']]),
            'thickness': to_float(parts[APARC_METRIC_INDEX['thickness']]),
        }
    return summary, regions


def parse_aseg_stats(stats_file):
    summary = {}
    struct_vol = {}
    data_start = None
    with open(stats_file, 'r', encoding='utf-8', errors='ignore') as handle:
        lines = handle.readlines()
    for index, line in enumerate(lines):
        if line.startswith('# Measure '):
            field, value = parse_measure_line(line)
            if field and value is not None:
                summary[field] = value
        if line.startswith('# ColHeaders'):
            data_start = index + 1
            break
    if data_start is None:
        raise ValueError(f'无法在 {stats_file} 中找到 ColHeaders')
    for line in lines[data_start:]:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) < 5:
            continue
        struct_name = parts[4]
        struct_vol[struct_name] = to_float(parts[3])
    return summary, struct_vol


def build_mri_database(subjects_dir):
    records = {}
    failed = []
    for subject_dir in sorted(path for path in Path(subjects_dir).iterdir() if path.is_dir()):
        lh_file = subject_dir / 'stats' / 'lh.aparc.stats'
        rh_file = subject_dir / 'stats' / 'rh.aparc.stats'
        aseg_file = subject_dir / 'stats' / 'aseg.stats'
        if not (lh_file.exists() and rh_file.exists() and aseg_file.exists()):
            continue
        project, sid = split_subject_name(subject_dir.name)
        subject_key = build_subject_key(project, sid)
        try:
            lh_summary, lh_regions = parse_aparc_stats(lh_file)
            rh_summary, rh_regions = parse_aparc_stats(rh_file)
            aseg_summary, aseg_struct = parse_aseg_stats(aseg_file)
            missing_lh = [region for region in DK68_REGIONS if region not in lh_regions]
            missing_rh = [region for region in DK68_REGIONS if region not in rh_regions]
            if missing_lh or missing_rh:
                raise ValueError(
                    f'DK 区域不完整，lh缺失{len(missing_lh)}个，rh缺失{len(missing_rh)}个'
                )
            records[subject_key] = {
                'ID': subject_key,
                'SubjectName': subject_dir.name,
                'lh_summary': lh_summary,
                'rh_summary': rh_summary,
                'aseg_summary': aseg_summary,
                'aseg_struct': aseg_struct,
                'lh_regions': lh_regions,
                'rh_regions': rh_regions,
            }
        except Exception as exc:
            failed.append((subject_dir.name, str(exc)))
    print(f'可解析 MRI 被试: {len(records)}')
    if failed:
        print(f'解析失败被试: {len(failed)}')
    return records, failed


def load_clinical_df(clinical_csv):
    clinical_df = pd.read_csv(clinical_csv)
    required = ['ID', 'Age', 'Sex', 'Center', 'Scanner']
    missing = [col for col in required if col not in clinical_df.columns]
    if missing:
        raise ValueError(f'临床表缺少必要列: {missing}')
    clinical_df = clinical_df.copy()
    clinical_df['ID'] = clinical_df['ID'].astype(str)
    clinical_df['Sex'] = clinical_df['Sex'].astype(str).str.strip()
    clinical_df = clinical_df.drop_duplicates('ID', keep='first')
    return clinical_df


def get_meta_row(clinical_row):
    return {
        'SITE': clinical_row['Center'],
        'SubjectID': clinical_row['ID'],
        'Vendor': 'Siemens',
        'FreeSurfer_Version': '7.4.1',
        'age': clinical_row['Age'],
        'sex': clinical_row['Sex'],
    }


def get_first_available(mapping, keys):
    for key in keys:
        if key in mapping and pd.notna(mapping[key]):
            return mapping[key]
    return pd.NA


def build_cortical_table(aligned_df, mri_db, metric):
    metric_suffix = 'thickavg' if metric == 'thickness' else 'surfavg'
    total_key = 'MeanThickness' if metric == 'thickness' else 'WhiteSurfArea'
    total_left_col = 'LThickness' if metric == 'thickness' else 'LSurfArea'
    total_right_col = 'RThickness' if metric == 'thickness' else 'RSurfArea'

    rows = []
    for _, clinical_row in aligned_df.iterrows():
        sid = clinical_row['ID']
        record = mri_db[sid]
        row = get_meta_row(clinical_row)
        row[total_left_col] = record['lh_summary'].get(total_key, pd.NA)
        row[total_right_col] = record['rh_summary'].get(total_key, pd.NA)
        for region in DK68_REGIONS:
            region_name = REGION_NAME_MAP.get(region, region)
            row[f'L_{region_name}_{metric_suffix}'] = record['lh_regions'][region][metric]
        for region in DK68_REGIONS:
            region_name = REGION_NAME_MAP.get(region, region)
            row[f'R_{region_name}_{metric_suffix}'] = record['rh_regions'][region][metric]
        rows.append(row)
    return pd.DataFrame(rows)


def build_subcortical_table(aligned_df, mri_db):
    rows = []
    for _, clinical_row in aligned_df.iterrows():
        sid = clinical_row['ID']
        record = mri_db[sid]
        row = get_meta_row(clinical_row)
        for output_col, source_keys in SUBCORTICAL_FEATURES:
            if output_col == 'ICV':
                row[output_col] = get_first_available(record['aseg_summary'], source_keys)
            else:
                row[output_col] = get_first_available(record['aseg_struct'], source_keys)
        rows.append(row)
    return pd.DataFrame(rows)


def sort_output(df):
    feature_cols = [col for col in df.columns if col not in META_COLUMNS]
    id_col = 'SubjectID'
    return df[META_COLUMNS + feature_cols].sort_values(id_col).reset_index(drop=True)


def export_tables(thickness_df, area_df, subcortical_df, suffix):
    thickness_path = OUTPUT_DIR / f'cortical-thickness-{suffix}.xlsx'
    area_path = OUTPUT_DIR / f'surface-area-{suffix}.xlsx'
    subcortical_path = OUTPUT_DIR / f'subcortical-volume-{suffix}.xlsx'

    sort_output(thickness_df).to_excel(thickness_path, index=False)
    sort_output(area_df).to_excel(area_path, index=False)
    sort_output(subcortical_df).to_excel(subcortical_path, index=False)

    print(f'已导出: {thickness_path}')
    print(f'已导出: {area_path}')
    print(f'已导出: {subcortical_path}')


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print('=' * 60)
    print('第一步：解析 FreeSurfer stats')
    print('=' * 60)
    mri_db, failed = build_mri_database(FS_SUBJECTS_DIR)

    print('\n' + '=' * 60)
    print('第二步：读取临床表并匹配被试')
    print('=' * 60)
    clinical_df = load_clinical_df(CLINICAL_CSV)
    aligned_df = clinical_df[clinical_df['ID'].isin(set(mri_db.keys()))].copy()
    aligned_df = aligned_df[aligned_df['Sex'].isin(['Male', 'Female'])]
    aligned_df = aligned_df.sort_values('ID').reset_index(drop=True)

    if aligned_df.empty:
        raise ValueError('Clinical_vars_merged.csv 与 FreeSurfer 结果没有可用交集')

    male_df = aligned_df[aligned_df['Sex'] == 'Male'].copy()
    female_df = aligned_df[aligned_df['Sex'] == 'Female'].copy()

    print(f'交集被试数: {len(aligned_df)}')
    print(f'Male: {len(male_df)}')
    print(f'Female: {len(female_df)}')

    print('\n' + '=' * 60)
    print('第三步：生成 6 个 example 风格表')
    print('=' * 60)

    export_tables(
        build_cortical_table(female_df, mri_db, 'thickness'),
        build_cortical_table(female_df, mri_db, 'area'),
        build_subcortical_table(female_df, mri_db),
        'female',
    )

    export_tables(
        build_cortical_table(male_df, mri_db, 'thickness'),
        build_cortical_table(male_df, mri_db, 'area'),
        build_subcortical_table(male_df, mri_db),
        'male',
    )

    if failed:
        failed_df = pd.DataFrame(failed, columns=['SubjectName', 'Error'])
        failed_path = OUTPUT_DIR / 'example_extraction_failed_subjects.csv'
        failed_df.to_csv(failed_path, index=False)
        print(f'解析失败记录: {failed_path}')

    print('\n完成。输出目录: ' + str(OUTPUT_DIR))


if __name__ == '__main__':
    main()
