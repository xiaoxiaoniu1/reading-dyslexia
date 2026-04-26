#!/bin/bash
#SBATCH --job-name=DK68_DGLM_all
#SBATCH --output=analysis/logs/DK68_DGLM_all_%j.out
#SBATCH --error=analysis/logs/DK68_DGLM_all_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=120:00:00
#SBATCH --account=""
set -euo pipefail

mkdir -p logs

echo "Host: $(hostname)"
echo "Start: $(date)"

cd /data/home/tqi/data1/share/after_freesurfer/CODE/analysis

VARIANT="${1:-all}"
ANALYSIS_TYPE="${2:-both}"

echo "Variant: ${VARIANT}"
echo "Analysis type: ${ANALYSIS_TYPE}"

run_one() {
  local label="$1"
  local script_path="$2"
  echo ""
  echo "=================================================="
  echo "Running ${label}"
  echo "Script: ${script_path}"
  echo "Analysis type: ${ANALYSIS_TYPE}"
  echo "=================================================="
  bash "${script_path}" "${ANALYSIS_TYPE}"
}

case "${VARIANT}" in
  overall)
    run_one "DK68 overall" "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_dglm_dk68.sh"
    ;;
  main)
    run_one "DK68 main" "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_dglm_dk68_main.sh"
    ;;
  sensitivity)
    run_one "DK68 sensitivity" "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_dglm_dk68_sensitivity.sh"
    ;;
  main_sex)
    run_one "DK68 main sex" "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_dglm_dk68_main_sex.sh"
    ;;
  sensitivity_sex)
    run_one "DK68 sensitivity sex" "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_dglm_dk68_sensitivity_sex.sh"
    ;;
  nonsex)
    run_one "DK68 overall" "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_dglm_dk68.sh"
    run_one "DK68 main" "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_dglm_dk68_main.sh"
    run_one "DK68 sensitivity" "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_dglm_dk68_sensitivity.sh"
    ;;
  sex)
    run_one "DK68 main sex" "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_dglm_dk68_main_sex.sh"
    run_one "DK68 sensitivity sex" "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_dglm_dk68_sensitivity_sex.sh"
    ;;
  all)
    run_one "DK68 overall" "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_dglm_dk68.sh"
    run_one "DK68 main" "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_dglm_dk68_main.sh"
    run_one "DK68 sensitivity" "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_dglm_dk68_sensitivity.sh"
    run_one "DK68 main sex" "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_dglm_dk68_main_sex.sh"
    run_one "DK68 sensitivity sex" "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_dglm_dk68_sensitivity_sex.sh"
    ;;
  *)
    echo "Unknown variant: ${VARIANT}"
    echo "Allowed: overall | main | sensitivity | main_sex | sensitivity_sex | nonsex | sex | all"
    exit 1
    ;;
esac

echo "End: $(date)"
