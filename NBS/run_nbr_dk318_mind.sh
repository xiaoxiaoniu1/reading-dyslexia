#!/bin/bash
#SBATCH --job-name=DK318_NBR
#SBATCH --output=NBS/logs/DK318_NBR_%j.out
#SBATCH --error=NBS/logs/DK318_NBR_%j.err
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
cd /data/home/tqi/data1/share/after_freesurfer/CODE/NBS

# 设置 R 临时目录，避免节点默认临时目录不可写
mkdir -p tmp
export TMPDIR="$(pwd)/tmp"

echo "TMPDIR: ${TMPDIR}"

# 控制分析模式：debug 或 formal
export NBR_ANALYSIS_MODE="formal"
echo "NBR_ANALYSIS_MODE: ${NBR_ANALYSIS_MODE}"

# 跑 NBR R 脚本
Rscript run_nbr_dk318_mind.R

echo "NBR R analysis finished. Submitting Python plotting shell..."
sbatch /data/home/tqi/data1/share/after_freesurfer/CODE/NBS/run_plot_DK318_NBS_components.sh

echo "End: $(date)"
