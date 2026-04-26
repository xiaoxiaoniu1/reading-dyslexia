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

# ========================================
# Default plotting options
# 直接修改下面两个变量即可控制绘图内容
# PLOT_MODE: overall / sex / both / all
# PLOT_VARIANT: auto / all / default / main / sensitivity
# ========================================
PLOT_MODE="overall"
PLOT_VARIANT="auto"

cd /data/home/tqi/data1/share/after_freesurfer/CODE/degree

echo "Start: $(date)"

# 激活 Python 环境
set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env
set -u

# 运行 Python 脚本
MODE="${DK318_DGLM_PLOT_MODE:-${PLOT_MODE}}"
VARIANT="${DK318_DGLM_VARIANT:-${PLOT_VARIANT}}"

CMD=(python degree_analysis_DK318_DGLM.py --mode "${MODE}" --variant "${VARIANT}")

if [[ -n "${DK318_DGLM_STRATIFIED_DIR:-}" ]]; then
  CMD+=(--stratified-dir "${DK318_DGLM_STRATIFIED_DIR}")
fi

echo "Plot mode: ${MODE}"
echo "Plot variant: ${VARIANT}"
if [[ -n "${DK318_DGLM_STRATIFIED_DIR:-}" ]]; then
  echo "Custom stratified dir: ${DK318_DGLM_STRATIFIED_DIR}"
fi

echo "Running command: ${CMD[*]}"
"${CMD[@]}"

echo "End: $(date)"
