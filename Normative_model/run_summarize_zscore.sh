#!/bin/bash
#SBATCH --job-name=summarize_zscore
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs/summarize_zscore_%j.log
#SBATCH --error=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs/summarize_zscore_%j.err

set -euo pipefail

mkdir -p /data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs

echo "========================================"
echo "Starting Z-score Summary"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "R Environment: rd_env_r"
echo "R Path: $(which R || true)"
echo "Start time: $(date)"
echo "Script: Summarize_Zscore.R"
echo "========================================"

set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

SCRIPT_DIR="/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model"
R_SCRIPT="Summarize_Zscore.R"

cd "${SCRIPT_DIR}"
Rscript "${R_SCRIPT}"

EXIT_CODE=$?

echo "========================================"
echo "End time: $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "========================================"

exit ${EXIT_CODE}
