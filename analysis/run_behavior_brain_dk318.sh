#!/bin/bash
#SBATCH --job-name=behavior_brain_scatter
#SBATCH --output=analysis/logs/behavior_brain_scatter_%j.out
#SBATCH --error=analysis/logs/behavior_brain_scatter_%j.err
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

# 如果用 conda 管理 Python，按你自己的环境名修改
# 例如：conda run -n rd_env_py python stratified_analysis_dd_td.py
source ~/.bashrc

# 进入代码目录（按实际路径调整）
cd /data/home/tqi/data1/share/after_freesurfer/CODE/analysis

conda run -n rd_env python /data/home/tqi/data1/share/after_freesurfer/CODE/analysis/behavior_brain_scatter.py

echo "End: $(date)"