#!/bin/bash
#SBATCH --job-name=extract_example_6tables_from_fs
#SBATCH --output=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs/extract_example_6tables_from_fs_%j.out
#SBATCH --error=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs/extract_example_6tables_from_fs_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""

set -euo pipefail

LOG_DIR=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs
SCRIPT_DIR=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model
PY_SCRIPT=${SCRIPT_DIR}/extract_example_6tables_from_fs.py
OUTPUT_DIR=/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/example_extracted

mkdir -p "${LOG_DIR}"
echo "Host: $(hostname)"
echo "Start: $(date)"
echo "Python script: ${PY_SCRIPT}"
echo "Output dir: ${OUTPUT_DIR}"

set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

cd "${SCRIPT_DIR}"

python3 "${PY_SCRIPT}"

python3 - <<'PY'
from pathlib import Path
import pandas as pd

output_dir = Path('/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/example_extracted')
files = {
    'female': output_dir / 'cortical-thickness-female.xlsx',
    'male': output_dir / 'cortical-thickness-male.xlsx',
}

print('--- Simple extraction summary ---')
for sex, path in files.items():
    if path.exists():
        df = pd.read_excel(path)
        print(f'{sex}_subjects: {len(df)}')
    else:
        print(f'{sex}_subjects: file not found -> {path}')
PY

echo "End: $(date)"
