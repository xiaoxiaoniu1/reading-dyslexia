#!/bin/bash
#SBATCH --job-name=DK68_DGLM_main
#SBATCH --output=analysis/logs/DK68_DGLM_main_%j.out
#SBATCH --error=analysis/logs/DK68_DGLM_main_%j.err
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
ANALYSIS_TYPE="${1:-both}"
RESULT_DIR="/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK68_DGLM_main"
VIOLIN_DIR="${RESULT_DIR}/VIOLIN_Followup"
echo "Analysis type: ${ANALYSIS_TYPE}"
echo "DGLM result dir: ${RESULT_DIR}"
echo "VIOLIN follow-up dir: ${VIOLIN_DIR}"
echo "===== STEP 1/3: Run DGLM model ====="
echo "Output dir (STEP 1/3): ${RESULT_DIR}"
Rscript dglm_dk68_degree_edge_main.R "${ANALYSIS_TYPE}"
echo "===== STEP 1/3 DONE ====="
if [[ "${ANALYSIS_TYPE}" == "both" || "${ANALYSIS_TYPE}" == "degree" ]]; then
  echo "===== STEP 2/3: Generate DGLM degree brainmaps ====="
  echo "Output dir (STEP 2/3): ${RESULT_DIR}"
  python /data/home/tqi/data1/share/after_freesurfer/CODE/degree/degree_analysis_DK68_DGLM.py --result-dir "${RESULT_DIR}"
  echo "===== STEP 2/3 DONE ====="
  echo "===== STEP 3/3: Run VIOLIN follow-up plots ====="
  echo "Output dir (STEP 3/3): ${VIOLIN_DIR}"
  bash /data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_violin_dk68.sh main
  echo "===== STEP 3/3 DONE ====="
fi
echo "End: $(date)"
