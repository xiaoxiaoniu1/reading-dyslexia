#!/bin/bash
#SBATCH --job-name=DK318_edge_ANOVA_adult_IQ
#SBATCH --output=logs/DK318_edge_ANOVA_adult_IQ_%j.out
#SBATCH --error=logs/DK318_edge_ANOVA_adult_IQ_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""
set -euo pipefail

cd /data/home/tqi/data1/share/after_freesurfer/CODE/edge

echo "Start: $(date)"

# 激活 Python 环境
set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env
set -u

# 运行 Python 脚本
python connect_DK318_ANOVA_adult_with_IQ.py

echo "End: $(date)"
