#!/bin/bash
#SBATCH --job-name=evaluate_calibration_holdout_all_features
#SBATCH --output=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs/evaluate_calibration_holdout_all_features_%j.out
#SBATCH --error=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs/evaluate_calibration_holdout_all_features_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --account=""

set -euo pipefail

mkdir -p /data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs
echo "Host: $(hostname)"
echo "Start: $(date)"

set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

cd /data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model

Rscript /data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/Evaluate_Calibration_Holdout_All_Features.R

echo "End: $(date)"
