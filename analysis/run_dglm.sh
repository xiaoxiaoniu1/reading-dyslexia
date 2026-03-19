#!/bin/bash
#SBATCH --job-name=DK318_DGLM
#SBATCH --output=logs/DK318_DGLM_%j.out
#SBATCH --error=logs/DK318_DGLM_%j.err
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

# 激活 conda 环境（临时关闭 -u 以避免 conda 激活脚本的未定义变量错误）
set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

# 进入代码目录
cd /data/home/tqi/data1/share/after_freesurfer/CODE/analysis

# 跑 DGLM R 脚本
Rscript dglm_dk318_degree_edge.R

echo "End: $(date)"
