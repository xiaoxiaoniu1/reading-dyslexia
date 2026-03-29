#!/bin/bash
#SBATCH --job-name=DK_normative_calibration_newsite
#SBATCH --output=Normative_model/logs/DK_normative_calibration_newsite_%j.out
#SBATCH --error=Normative_model/logs/DK_normative_calibration_newsite_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""

set -euo pipefail

mkdir -p logs
echo "Host: $(hostname)"
echo "Start: $(date)"

set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

cd /data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model

Rscript /data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/Calibration-normative-model-using-new-dataset.R

echo "End: $(date)"
