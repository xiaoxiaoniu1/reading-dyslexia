#!/bin/bash
#SBATCH --job-name=DK318_sig_brainmaps
#SBATCH --output=logs/DK318_sig_brainmaps_%j.out
#SBATCH --error=logs/DK318_sig_brainmaps_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --account=""
set -euo pipefail

cd /data/home/tqi/data1/share/after_freesurfer/CODE/degree

echo "Start: $(date)"
echo "Script: plot_significant_DK318_stratified_brainmaps.py"
echo "Job ID: ${SLURM_JOB_ID:-NA}"

set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env
set -u

python plot_significant_DK318_stratified_brainmaps.py

echo "End: $(date)"
