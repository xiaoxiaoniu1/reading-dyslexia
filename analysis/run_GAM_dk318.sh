#!/bin/bash
#SBATCH --job-name=DK318_GAM
#SBATCH --output=logs/DK318_GAM_%j.out
#SBATCH --error=logs/DK318_GAM_%j.err
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

# 方式A：如果你用 conda 管 R（常见）
source ~/.bashrc   # 可留可不留

cd /data1/tqi/share/after_freesurfer/CODE/analysis

# 跑【连续年龄】版本的 R 脚本（文件名带括号，一定要加引号）
conda run -n rd_env_r Rscript "GAM_dk318_degree_edge.R"

echo "End: $(date)"