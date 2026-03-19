#!/bin/bash
#SBATCH --job-name=DK318_brain_region_check
#SBATCH --output=logs/DK318_brain_region_check_%j.out
#SBATCH --error=logs/DK318_brain_region_check_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --account=""

set -euo pipefail

cd /data/home/tqi/data1/share/after_freesurfer/CODE/degree

echo "Start: $(date)"

set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env
set -u

python DK318_brain_region_check.py

echo "End: $(date)"