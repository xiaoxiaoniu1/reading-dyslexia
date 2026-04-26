#!/bin/bash
#SBATCH --job-name=DK318_DGLM_strat
#SBATCH --output=logs/DK318_DGLM_stratified_%j.out
#SBATCH --error=logs/DK318_DGLM_stratified_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --account=""
set -euo pipefail

mkdir -p logs
cd /data/home/tqi/data1/share/after_freesurfer/CODE/analysis

echo "Host: $(hostname)"
echo "Start: $(date)"
echo "Working directory: $(pwd)"

set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

ANALYSIS_TYPE="${1:-both}"
if [[ "${ANALYSIS_TYPE}" != "both" && "${ANALYSIS_TYPE}" != "degree" && "${ANALYSIS_TYPE}" != "edge" ]]; then
  echo "Invalid analysis type: ${ANALYSIS_TYPE}. Use one of: both, degree, edge." >&2
  exit 1
fi

RESULT_DIR="/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_DGLM_stratified"
BRAINMAP_DIR="${RESULT_DIR}/stratified_brainmaps"
VIOLIN_DIR="${RESULT_DIR}/VIOLIN_Followup_stratified"
R_SCRIPT="/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/dglm_stratified_dk318.R"

echo "Analysis type: ${ANALYSIS_TYPE}"
echo "DGLM result dir: ${RESULT_DIR}"
echo "Brainmap output dir: ${BRAINMAP_DIR}"
echo "VIOLIN output dir: ${VIOLIN_DIR}"
echo "R script: ${R_SCRIPT}"

echo "===== STEP 1/3: Run stratified DGLM ====="
echo "Output dir (STEP 1/3): ${RESULT_DIR}"
Rscript --vanilla "${R_SCRIPT}" "${ANALYSIS_TYPE}" --out-dir "${RESULT_DIR}" --run-plots false
echo "===== STEP 1/3 DONE ====="
echo "Output dir (STEP 1/3): ${RESULT_DIR}"

if [[ "${ANALYSIS_TYPE}" == "both" || "${ANALYSIS_TYPE}" == "degree" ]]; then
  echo "===== STEP 2/3: Plot stratified significant degree brainmaps ====="
  echo "Output dir (STEP 2/3): ${BRAINMAP_DIR}"
  python /data/home/tqi/data1/share/after_freesurfer/CODE/degree/plot_significant_DK318_stratified_brainmaps.py \
    --result-dir "${RESULT_DIR}"
  echo "===== STEP 2/3 DONE ====="
  echo "Output dir (STEP 2/3): ${BRAINMAP_DIR}"

  echo "===== STEP 3/3: Run stratified VIOLIN follow-up for significant interaction ROIs ====="
  echo "Output dir (STEP 3/3): ${VIOLIN_DIR}"
  bash /data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_violin_dk318_stratified.sh \
    "${RESULT_DIR}" \
    "DK318 stratified"
  echo "===== STEP 3/3 DONE ====="
  echo "Output dir (STEP 3/3): ${VIOLIN_DIR}"
else
  echo "===== STEP 2/3 SKIPPED: Brainmaps require degree analysis ====="
  echo "===== STEP 3/3 SKIPPED: VIOLIN follow-up requires degree interaction ROIs ====="
fi

echo "End: $(date)"

