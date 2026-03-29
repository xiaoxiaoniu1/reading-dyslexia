#!/bin/bash
#SBATCH -J combat_dk318               # 作业名
#SBATCH -p partition_1              # 队列名：按你们集群实际改
#SBATCH -N 1
#SBATCH --cpus-per-task=4           # 用几个核
#SBATCH --mem=16G                   # 内存
#SBATCH -t 72:00:00                 # 最长运行时间
#SBATCH -o combat/logs/combat_dk318_%j.out
#SBATCH -e combat/logs/combat_dk318_%j.err

set -euo pipefail

echo "=========================================="
echo "开始执行 DK-318 Combat 批次校正"
echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================="

# 激活 conda 环境
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env

# 运行 combat_dk318.py
python3 /data/home/tqi/data1/share/after_freesurfer/CODE/combat/combat_dk318.py

echo "=========================================="
echo "DK-318 Combat 批次校正完成"
echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================="
