#!/bin/bash
#SBATCH --job-name=DK68_VIOLIN
#SBATCH --output=analysis/logs/DK68_VIOLIN_%j.out
#SBATCH --error=analysis/logs/DK68_VIOLIN_%j.err
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

set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

cd /data/home/tqi/data1/share/after_freesurfer/CODE/analysis
mkdir -p tmp
export TMPDIR="$(pwd)/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

echo "TMPDIR: ${TMPDIR}"

VARIANT="${1:-overall}"
case "${VARIANT}" in
  overall)
    RESULT_DIR="/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK68_DGLM"
    LABEL="DK68 overall"
    ;;
  main)
    RESULT_DIR="/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK68_DGLM_main"
    LABEL="DK68 main"
    ;;
  sensitivity)
    RESULT_DIR="/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK68_DGLM_sensitivity"
    LABEL="DK68 sensitivity"
    ;;
  main_sex)
    RESULT_DIR="/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK68_DGLM_main_sex"
    LABEL="DK68 main sex-stratified"
    ;;
  sensitivity_sex)
    RESULT_DIR="/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK68_DGLM_sensitivity_sex"
    LABEL="DK68 sensitivity sex-stratified"
    ;;
  *)
    echo "Unknown variant: ${VARIANT}"
    echo "Use one of: overall | main | sensitivity | main_sex | sensitivity_sex"
    exit 1
    ;;
esac

echo "Variant: ${VARIANT}"
echo "Result dir: ${RESULT_DIR}"

Rscript violin_dk68_degree_edge.R \
  --result-dir "${RESULT_DIR}" \
  --analysis-label "${LABEL}"

echo "End: $(date)"
