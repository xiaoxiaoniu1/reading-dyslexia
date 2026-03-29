#!/bin/bash
#SBATCH -J DK68_site            # 作业名
#SBATCH -p partition_1              # 队列名：按你们集群实际改
#SBATCH -N 1
#SBATCH --cpus-per-task=4           # 用几个核
#SBATCH --mem=16G                   # 内存
#SBATCH -t 12:00:00                 # 最长运行时间
#SBATCH -o analysis_logs/DK68_site_%j.out
#SBATCH -e analysis_logs/DK68_site_%j.err

set -euo pipefail

echo "=========================================="
echo "开始执行 DK-68 ANOVA 分析 (原始数据 + Site作为协变量)"
echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================="

# 激活 conda 环境（临时关闭 -u 以避免 conda 激活脚本的未定义变量错误）
set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

# 运行 anova_dk68_degree_edge_with_site.R
Rscript /data/home/tqi/data1/share/after_freesurfer/CODE/analysis/anova_dk68_degree_edge_with_site.R

echo "=========================================="
echo "DK-68 ANOVA 分析完成"
echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================="
