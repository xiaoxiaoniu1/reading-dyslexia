#!/bin/bash
#SBATCH --job-name=check_uncalibrated_TD_bias
#SBATCH --output=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs/check_uncalibrated_TD_bias_%j.out
#SBATCH --error=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs/check_uncalibrated_TD_bias_%j.err
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

Rscript /data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/check_uncalibrated_TD_bias.R

echo "End: $(date)"
