#!/bin/bash
#SBATCH --job-name=DK68_ANOVA
#SBATCH --output=analysis/logs/DK68_ANOVA%j.out
#SBATCH --error=analysis/logs/DK68_ANOVA%j.err
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

# 进入代码目录（按你的实际路径改）
cd /data/home/tqi/data1/share/after_freesurfer/CODE/analysis

# 设置 R 临时目录，避免节点默认临时目录不可写
mkdir -p tmp
export TMPDIR="$(pwd)/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

echo "TMPDIR: ${TMPDIR}"

# 跑R脚本
Rscript anova_dk68_degree_edge.R

echo "End: $(date)"