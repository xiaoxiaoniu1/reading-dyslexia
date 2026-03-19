#!/bin/bash
#SBATCH --job-name=DK318_STRAT_DD_TD
#SBATCH --output=logs/DK318_STRAT_DD_TD_%j.out
#SBATCH --error=logs/DK318_STRAT_DD_TD_%j.err
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

# 方式A：直接用默认 python（前提是这个环境里装好了 pandas/numpy/statsmodels 等）
#python stratified_analysis_dd_td.py

# 如果你希望指定 conda 环境，可以改成下面这种：
conda run -n rd_env python stratified_analysis_dd_td.py

echo "End: $(date)"