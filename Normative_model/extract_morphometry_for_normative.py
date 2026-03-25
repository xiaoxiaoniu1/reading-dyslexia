#!/usr/bin/env python3
"""
按官方 Dataset-norms 风格导出 FreeSurfer 形态学表与临床表
"""

import json
import os
import re
from pathlib import Path

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

HEMIS = ['lh', 'rh']
APARC_METRIC_INDEX = {'area': 2, 'volume': 3, 'thickness': 4}
SHEET_ORDER = [
    'Global.table',
    'aseg.vol.table',
    'lh.aparc.volume.table',
    'rh.aparc.volume.table',
    'lh.aparc.thickness.table',
    'rh.aparc.thickness.table',
    'lh.aparc.area.table',
    'rh.aparc.area.table'
]
PROJECT_ALIAS = {
    'adulec': 'adolec',
    'aduld': 'adoled',
    'controldataadultscontrol': 'controldataadultscontrol',
    'controldataagecontrol': 'controldataagecontrol'
}

ASEG_STRUCT_OUTPUT_MAP = {
    'Left-Accumbens-area': 'Left.Accumbens.area',
    'Right-Accumbens-area': 'Right.Accumbens.area',
    'Left-Amygdala': 'Left.Amygdala',
    'Right-Amygdala': 'Right.Amygdala',
    'Left-Caudate': 'Left.Caudate',
    'Right-Caudate': 'Right.Caudate',
    'Left-Hippocampus': 'Left.Hippocampus',
    'Right-Hippocampus': 'Right.Hippocampus',
    'Left-Pallidum': 'Left.Pallidum',
    'Right-Pallidum': 'Right.Pallidum',
    'Left-Putamen': 'Left.Putamen',
    'Right-Putamen': 'Right.Putamen',
    'Left-Thalamus': 'Left.Thalamus',
    'Right-Thalamus': 'Right.Thalamus',
    'Left-VentralDC': 'Left.VentralDC',
    'Right-VentralDC': 'Right.VentralDC',
    'Brain-Stem': 'Brain.Stem'
}

CLINICAL_COLUMNS = [
    'ID', 'Freesufer_Path1', 'Freesufer_Path2', 'Freesufer_Path3',
    'Clinical_Path1', 'Clinical_Path2', 'Clinical_Path3', 'Age', 'Sex',
    'Diagnosis', 'Scanner', 'Age_y', 'Age_m', 'Age_d', 'Filed_strength',
    'Manufacturer', 'euler_number_l', 'euler_number_r', 'City', 'Province',
    'Country', 'Center'
]


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
    m = re.match(r'^(.*?)(\d+)$', subject_name)
    if m:
        return m.group(1), m.group(2)
    return subject_name, subject_name


def to_float(text):
    return float(str(text).strip())


def parse_measure_line(line):
    content = line[len('# Measure '):].strip()
    parts = [p.strip() for p in content.split(',')]
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
    with open(stats_file, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith('# Measure '):
            field, value = parse_measure_line(line)
            if field and value is not None:
                summary[field] = value
        if line.startswith('# ColHeaders'):
            data_start = i + 1
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
            'thickness': to_float(parts[APARC_METRIC_INDEX['thickness']])
        }
    return summary, regions


def parse_aseg_stats(stats_file):
    summary = {}
    struct_vol = {}
    data_start = None
    with open(stats_file, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith('# Measure '):
            field, value = parse_measure_line(line)
            if field and value is not None:
                summary[field] = value
        if line.startswith('# ColHeaders'):
            data_start = i + 1
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


def map_sex(value):
    if pd.isna(value):
        return pd.NA
    text = str(value).strip().lower()
    if text.isdigit():
        v = int(text)
        if v == 1:
            return 'Male'
        if v == 2:
            return 'Female'
        if v == 0:
            return 'Female'
    if text.startswith('m') or 'male' in text:
        return 'Male'
    if text.startswith('f') or 'female' in text:
        return 'Female'
    return pd.NA


def map_diagnosis(value):
    if pd.isna(value):
        return pd.NA
    text = str(value).strip().lower()
    if text.isdigit():
        v = int(text)
        if v in (0, 1):
            return 'TD' if v == 0 else 'DD'
        if v in (1, 2):
            return 'TD' if v == 1 else 'DD'
    if 'dd' in text or 'dys' in text:
        return 'DD'
    if 'td' in text or 'control' in text:
        return 'TD'
    return pd.NA


def map_center(value):
    if pd.isna(value):
        return pd.NA
    text = str(value).strip()
    if text == '':
        return pd.NA
    if text.isdigit():
        return f'Site_{text}'
    return text


def map_scanner(value):
    if pd.isna(value):
        return pd.NA
    text = str(value).strip()
    if text == '':
        return pd.NA
    return text


def infer_field_strength(scanner_value):
    if pd.isna(scanner_value):
        return pd.NA
    text = str(scanner_value).lower()
    if '3t' in text or '3.0' in text:
        return '3T'
    if '1.5t' in text or '1.5' in text:
        return '1.5T'
    return pd.NA


def infer_manufacturer(scanner_value):
    if pd.isna(scanner_value):
        return pd.NA
    text = str(scanner_value).lower()
    if 'siemens' in text:
        return 'Siemens'
    if 'ge' in text:
        return 'GE'
    if 'philips' in text:
        return 'Philips'
    return pd.NA


def build_mri_database(subjects_dir):
    records = {}
    failed = []
    subjects_path = Path(subjects_dir)
    for subject_dir in sorted([x for x in subjects_path.iterdir() if x.is_dir()]):
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
            missing_lh = [r for r in DK68_REGIONS if r not in lh_regions]
            missing_rh = [r for r in DK68_REGIONS if r not in rh_regions]
            if missing_lh or missing_rh:
                raise ValueError(
                    f'DK 区域不完整，lh缺失{len(missing_lh)}个，rh缺失{len(missing_rh)}个'
                )
            records[subject_key] = {
                'ID': subject_key,
                'SubjectName': subject_dir.name,
                'Freesufer_Path1': str(subject_dir),
                'Freesufer_Path2': str(lh_file),
                'Freesufer_Path3': str(rh_file),
                'aseg_path': str(aseg_file),
                'lh_summary': lh_summary,
                'rh_summary': rh_summary,
                'aseg_summary': aseg_summary,
                'aseg_struct': aseg_struct,
                'lh_regions': lh_regions,
                'rh_regions': rh_regions
            }
        except Exception as exc:
            failed.append((subject_dir.name, str(exc)))
    print(f'可解析 MRI 被试: {len(records)}')
    if failed:
        print(f'解析失败被试: {len(failed)}')
    return records, failed


def build_clinical_df(clinical_file):
    clinical_df = pd.read_excel(clinical_file)
    required_cols = ['original-project', 'id_old', 'age_month', 'sex', 'site', 'group_d_or_c']
    missing_required = [c for c in required_cols if c not in clinical_df.columns]
    if missing_required:
        raise ValueError(f'临床表缺少必要列: {missing_required}')
    clinical_df = clinical_df.copy()
    clinical_df['subject_key'] = clinical_df.apply(
        lambda row: build_subject_key(row['original-project'], row['id_old']),
        axis=1
    )
    age_years = pd.to_numeric(clinical_df['age_month'], errors='coerce')
    clinical_df['Age_y'] = age_years
    clinical_df['Age_m'] = age_years * 12.0
    clinical_df['Age_d'] = age_years * 365.25
    clinical_df['Sex'] = clinical_df['sex'].map(map_sex)
    clinical_df['Diagnosis'] = clinical_df['group_d_or_c'].map(map_diagnosis)
    clinical_df['Center'] = clinical_df['site'].map(map_center)
    clinical_df['Scanner'] = clinical_df['scanner'].map(map_scanner) if 'scanner' in clinical_df.columns else pd.NA
    clinical_df['Filed_strength'] = clinical_df['Scanner'].map(infer_field_strength)
    clinical_df['Manufacturer'] = clinical_df['Scanner'].map(infer_manufacturer)
    clinical_df['City'] = clinical_df['location'] if 'location' in clinical_df.columns else pd.NA
    clinical_df['Province'] = pd.NA
    clinical_df['Country'] = pd.NA
    clinical_df['Clinical_Path1'] = clinical_file
    clinical_df['Clinical_Path2'] = clinical_df['original-project'].astype(str)
    clinical_df['Clinical_Path3'] = clinical_df['id_old'].map(normalize_id_value)
    clinical_df['Age'] = age_years
    clinical_df = clinical_df.drop_duplicates('subject_key', keep='first')
    return clinical_df


def build_global_sheet(aligned_ids, mri_db):
    rows = []
    for sid in aligned_ids:
        m = mri_db[sid]
        aseg_summary = m['aseg_summary']
        lh_summary = m['lh_summary']
        rh_summary = m['rh_summary']
        lh_mean = lh_summary.get('MeanThickness', pd.NA)
        rh_mean = rh_summary.get('MeanThickness', pd.NA)
        row = {
            'ID': sid,
            'GMV': aseg_summary.get('TotalGrayVol', pd.NA),
            'sGMV': aseg_summary.get('SubCortGrayVol', pd.NA),
            'WMV': aseg_summary.get('CerebralWhiteMatterVol', pd.NA),
            'Ventricles': aseg_summary.get('VentricleChoroidVol', pd.NA),
            'lhMeanThickness': lh_mean,
            'rhMeanThickness': rh_mean,
            'meanCT2': (lh_mean + rh_mean) / 2.0 if pd.notna(lh_mean) and pd.notna(rh_mean) else pd.NA,
            'lh_totaISA2': lh_summary.get('WhiteSurfArea', pd.NA),
            'rh_totaISA2': rh_summary.get('WhiteSurfArea', pd.NA),
            'TCV': aseg_summary.get('eTIV', pd.NA),
            'euler_number_l': aseg_summary.get('lhSurfaceHoles', pd.NA),
            'euler_number_r': aseg_summary.get('rhSurfaceHoles', pd.NA),
            'Freesufer_Path1': m['Freesufer_Path1'],
            'Freesufer_Path2': m['Freesufer_Path2'],
            'Freesufer_Path3': m['Freesufer_Path3']
        }
        rows.append(row)
    return pd.DataFrame(rows)


def build_aseg_sheet(aligned_ids, mri_db):
    rows = []
    for sid in aligned_ids:
        m = mri_db[sid]
        aseg_summary = m['aseg_summary']
        struct = m['aseg_struct']
        row = {'ID': sid}
        for fs_name, out_name in ASEG_STRUCT_OUTPUT_MAP.items():
            row[out_name] = struct.get(fs_name, pd.NA)
        cerebellum_parts = [
            struct.get('Left-Cerebellum-Cortex', 0.0),
            struct.get('Right-Cerebellum-Cortex', 0.0),
            struct.get('Left-Cerebellum-White-Matter', 0.0),
            struct.get('Right-Cerebellum-White-Matter', 0.0)
        ]
        row['cerebellum_total'] = sum(x for x in cerebellum_parts if pd.notna(x))
        row['BrainSegVol'] = aseg_summary.get('BrainSegVol', pd.NA)
        row['BrainSegVolNotVent'] = aseg_summary.get('BrainSegVolNotVent', pd.NA)
        row['lhCortexVol'] = aseg_summary.get('lhCortexVol', pd.NA)
        row['rhCortexVol'] = aseg_summary.get('rhCortexVol', pd.NA)
        row['CortexVol'] = aseg_summary.get('CortexVol', pd.NA)
        row['CerebralWhiteMatterVol'] = aseg_summary.get('CerebralWhiteMatterVol', pd.NA)
        row['SubCortGrayVol'] = aseg_summary.get('SubCortGrayVol', pd.NA)
        row['TotalGrayVol'] = aseg_summary.get('TotalGrayVol', pd.NA)
        row['EstimatedTotalIntraCranialVol'] = aseg_summary.get('eTIV', pd.NA)
        row['lhSurfaceHoles'] = aseg_summary.get('lhSurfaceHoles', pd.NA)
        row['rhSurfaceHoles'] = aseg_summary.get('rhSurfaceHoles', pd.NA)
        row['Freesufer_Path1'] = m['Freesufer_Path1']
        row['Freesufer_Path2'] = m['aseg_path']
        row['Freesufer_Path3'] = m['Freesufer_Path3']
        rows.append(row)
    return pd.DataFrame(rows)


def build_aparc_sheet(aligned_ids, mri_db, hemi, metric):
    rows = []
    for sid in aligned_ids:
        m = mri_db[sid]
        hemi_regions = m[f'{hemi}_regions']
        aseg_summary = m['aseg_summary']
        hemi_summary = m[f'{hemi}_summary']
        row = {'ID': sid}
        for region in DK68_REGIONS:
            row[f'{hemi}_{region}_{metric}'] = hemi_regions[region][metric]
        row['MeanThickness'] = hemi_summary.get('MeanThickness', pd.NA)
        row['BrainSegVolNotVent'] = aseg_summary.get('BrainSegVolNotVent', pd.NA)
        row['eTIV'] = aseg_summary.get('eTIV', pd.NA)
        row['Freesufer_Path1'] = m['Freesufer_Path1']
        row['Freesufer_Path2'] = m['Freesufer_Path2']
        row['Freesufer_Path3'] = m['Freesufer_Path3']
        rows.append(row)
    return pd.DataFrame(rows)


def build_clinical_output(aligned_ids, clinical_df, mri_db):
    clinical_idx = clinical_df.set_index('subject_key')
    rows = []
    for sid in aligned_ids:
        c = clinical_idx.loc[sid]
        m = mri_db[sid]
        row = {
            'ID': sid,
            'Freesufer_Path1': m['Freesufer_Path1'],
            'Freesufer_Path2': m['Freesufer_Path2'],
            'Freesufer_Path3': m['Freesufer_Path3'],
            'Clinical_Path1': c['Clinical_Path1'],
            'Clinical_Path2': c['Clinical_Path2'],
            'Clinical_Path3': c['Clinical_Path3'],
            'Age': c['Age'],
            'Sex': c['Sex'],
            'Diagnosis': c['Diagnosis'],
            'Scanner': c['Scanner'],
            'Age_y': c['Age_y'],
            'Age_m': c['Age_m'],
            'Age_d': c['Age_d'],
            'Filed_strength': c['Filed_strength'],
            'Manufacturer': c['Manufacturer'],
            'euler_number_l': m['aseg_summary'].get('lhSurfaceHoles', pd.NA),
            'euler_number_r': m['aseg_summary'].get('rhSurfaceHoles', pd.NA),
            'City': c['City'],
            'Province': c['Province'],
            'Country': c['Country'],
            'Center': c['Center']
        }
        rows.append(row)
    out = pd.DataFrame(rows)
    for col in CLINICAL_COLUMNS:
        if col not in out.columns:
            out[col] = pd.NA
    return out[CLINICAL_COLUMNS]


def validate_outputs(clinical_out, sheets):
    clinical_ids = clinical_out['ID'].tolist()
    clinical_cols = set(clinical_out.columns) - {'ID'}
    allowed_overlap = {
        'Freesufer_Path1', 'Freesufer_Path2', 'Freesufer_Path3',
        'euler_number_l', 'euler_number_r'
    }
    for sheet_name, df in sheets.items():
        if len(df) != len(clinical_out):
            raise ValueError(f'{sheet_name} 行数 {len(df)} 与 Clinical_vars.csv {len(clinical_out)} 不一致')
        if df['ID'].tolist() != clinical_ids:
            raise ValueError(f'{sheet_name} 的被试顺序与 Clinical_vars.csv 不一致')
        overlap_cols = (clinical_cols & set(df.columns)) - allowed_overlap
        if overlap_cols:
            raise ValueError(f'{sheet_name} 出现临床列混入: {sorted(overlap_cols)}')
    for hemi in HEMIS:
        for metric in ['thickness', 'area', 'volume']:
            sheet_name = f'{hemi}.aparc.{metric}.table'
            expected = [f'{hemi}_{region}_{metric}' for region in DK68_REGIONS]
            miss = [c for c in expected if c not in sheets[sheet_name].columns]
            if miss:
                raise ValueError(f'{sheet_name} 缺失 DK 字段: {miss[:5]} ...')
            for region in ['frontalpole', 'temporalpole']:
                col = f'{hemi}_{region}_{metric}'
                if col not in sheets[sheet_name].columns:
                    raise ValueError(f'{sheet_name} 缺失关键区域字段 {col}')
    if len(sheets) != 8:
        raise ValueError(f'MR_measures.xlsx sheet 数量错误，当前 {len(sheets)}')


def build_field_source_map(sheets):
    source_map = {}
    for sheet_name, df in sheets.items():
        col_source = {}
        for col in df.columns:
            if col == 'ID':
                col_source[col] = 'derived:subject_key'
            elif col.startswith('lh_') or col.startswith('rh_'):
                col_source[col] = 'direct:aparc_table'
            elif col in {'GMV', 'sGMV', 'WMV', 'Ventricles', 'TCV'}:
                col_source[col] = 'direct:aseg_summary'
            elif col in {'lhMeanThickness', 'rhMeanThickness', 'lh_totaISA2', 'rh_totaISA2'}:
                col_source[col] = 'direct:aparc_summary'
            elif col in {'meanCT2', 'cerebellum_total'}:
                col_source[col] = 'derived:aggregated'
            elif col in {'BrainSegVol', 'BrainSegVolNotVent', 'lhCortexVol', 'rhCortexVol', 'CortexVol',
                         'CerebralWhiteMatterVol', 'SubCortGrayVol', 'TotalGrayVol',
                         'EstimatedTotalIntraCranialVol', 'lhSurfaceHoles', 'rhSurfaceHoles',
                         'eTIV', 'MeanThickness'}:
                col_source[col] = 'direct:aseg_or_aparc_summary'
            elif col.startswith('Left.') or col.startswith('Right.') or col == 'Brain.Stem':
                col_source[col] = 'direct:aseg_table'
            elif col.startswith('Freesufer_Path'):
                col_source[col] = 'direct:path'
            elif col.startswith('euler_number'):
                col_source[col] = 'direct:aseg_summary'
            else:
                col_source[col] = 'unknown'
        source_map[sheet_name] = col_source
    return source_map


def build_log(clinical_out, sheets):
    log = {}
    for name, df in sheets.items():
        log[name] = {
            'rows': int(len(df)),
            'columns': list(df.columns),
            'missing_count': {k: int(v) for k, v in df.isna().sum().to_dict().items()}
        }
    log['Clinical_vars.csv'] = {
        'rows': int(len(clinical_out)),
        'columns': list(clinical_out.columns),
        'missing_count': {k: int(v) for k, v in clinical_out.isna().sum().to_dict().items()},
        'field_source': {c: 'clinical_or_derived' for c in clinical_out.columns}
    }
    for sheet_name, mapping in build_field_source_map(sheets).items():
        log[sheet_name]['field_source'] = mapping
    return log


def main():
    base_dir = '/data/home/tqi/data1/share/after_freesurfer'
    fs_subjects_dir = f'{base_dir}/fs_subjects_all'
    clinical_file = f'{base_dir}/FILE/all_data_cqt_without_age_control.xlsx'
    output_dir = f'{base_dir}/FILE/Normative_data'

    os.makedirs(output_dir, exist_ok=True)

    print('=' * 60)
    print('第一步：解析 FreeSurfer stats（lh/rh.aparc.stats + aseg.stats）')
    print('=' * 60)
    mri_db, failed = build_mri_database(fs_subjects_dir)
    mri_ids = set(mri_db.keys())

    print('\n' + '=' * 60)
    print('第二步：构建临床表并生成 subject key')
    print('=' * 60)
    clinical_df = build_clinical_df(clinical_file)
    clinical_ids = set(clinical_df['subject_key'].tolist())

    overlap_ids = sorted(mri_ids & clinical_ids)
    print(f'MRI 被试数: {len(mri_ids)}')
    print(f'临床被试数: {len(clinical_ids)}')
    print(f'交集被试数: {len(overlap_ids)}')
    print(f'仅MRI: {len(mri_ids - clinical_ids)}')
    print(f'仅临床: {len(clinical_ids - mri_ids)}')
    if not overlap_ids:
        raise ValueError('MRI 与临床无交集，无法导出')

    clinical_out = build_clinical_output(overlap_ids, clinical_df, mri_db)
    clinical_out_path = f'{output_dir}/Clinical_vars.csv'
    clinical_out.to_csv(clinical_out_path, index=False)

    sheets = {
        'Global.table': build_global_sheet(overlap_ids, mri_db),
        'aseg.vol.table': build_aseg_sheet(overlap_ids, mri_db),
        'lh.aparc.volume.table': build_aparc_sheet(overlap_ids, mri_db, 'lh', 'volume'),
        'rh.aparc.volume.table': build_aparc_sheet(overlap_ids, mri_db, 'rh', 'volume'),
        'lh.aparc.thickness.table': build_aparc_sheet(overlap_ids, mri_db, 'lh', 'thickness'),
        'rh.aparc.thickness.table': build_aparc_sheet(overlap_ids, mri_db, 'rh', 'thickness'),
        'lh.aparc.area.table': build_aparc_sheet(overlap_ids, mri_db, 'lh', 'area'),
        'rh.aparc.area.table': build_aparc_sheet(overlap_ids, mri_db, 'rh', 'area')
    }

    validate_outputs(clinical_out, sheets)

    mr_xlsx = f'{output_dir}/MR_measures.xlsx'
    with pd.ExcelWriter(mr_xlsx, engine='openpyxl') as writer:
        for sheet_name in SHEET_ORDER:
            sheets[sheet_name].to_excel(writer, sheet_name=sheet_name, index=False)

    log_data = build_log(clinical_out, sheets)
    log_data['summary'] = {
        'subjects_in_mri': len(mri_ids),
        'subjects_in_clinical': len(clinical_ids),
        'subjects_intersection': len(overlap_ids),
        'failed_subjects_count': len(failed),
        'failed_subjects': failed[:100]
    }
    log_file = f'{output_dir}/extraction_log.json'
    with open(log_file, 'w', encoding='utf-8') as f:
        json.dump(log_data, f, ensure_ascii=False, indent=2)

    print('\n' + '=' * 60)
    print('导出完成')
    print('=' * 60)
    print(f'Clinical_vars.csv: {clinical_out_path}')
    print(f'MR_measures.xlsx: {mr_xlsx}')
    print(f'日志文件: {log_file}')


if __name__ == '__main__':
    main()
