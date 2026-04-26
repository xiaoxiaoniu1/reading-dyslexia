#!/bin/bash
#SBATCH --job-name=DK318_VIOLIN_strat
#SBATCH --output=logs/DK318_VIOLIN_stratified_%j.out
#SBATCH --error=logs/DK318_VIOLIN_stratified_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""

set -euo pipefail

mkdir -p logs
cd /data/home/tqi/data1/share/after_freesurfer/CODE/analysis

echo "Host: $(hostname)"
echo "Start: $(date)"

set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

mkdir -p tmp
export TMPDIR="$(pwd)/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

echo "TMPDIR: ${TMPDIR}"

RESULT_DIR="${1:-/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_DGLM_stratified}"
LABEL="${2:-DK318 stratified}"

echo "Result dir: ${RESULT_DIR}"
echo "Label: ${LABEL}"

Rscript violin_dk318_stratified_degree.R \
  --result-dir "${RESULT_DIR}" \
  --analysis-label "${LABEL}"

echo "End: $(date)"
