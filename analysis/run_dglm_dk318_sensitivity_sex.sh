#!/bin/bash
#SBATCH --job-name=DK318_DGLM_sens_sex
#SBATCH --output=analysis/logs/DK318_DGLM_sensitivity_sex_%j.out
#SBATCH --error=analysis/logs/DK318_DGLM_sensitivity_sex_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""


set -euo pipefail

mkdir -p logs
echo "==== TEMP DEBUG INFO ===="
echo "Host: $(hostname)"
echo "PWD: $(pwd)"
echo "TMPDIR(before): ${TMPDIR:-<unset>}"

echo "--- env tmp vars ---"
env | grep -E 'TMP|TEMP'

echo "--- system tmp status ---"
ls -ld /tmp || true
df -h /tmp || true
df -i /tmp || true

echo "--- current dir tmp status ---"
mkdir -p tmp
ls -ld tmp
df -h .
df -i .

echo "--- write test /tmp ---"
touch /tmp/test_tmp_write_$$ && echo "write /tmp ok" || echo "write /tmp failed"
rm -f /tmp/test_tmp_write_$$ || true

echo "--- write test local tmp ---"
touch tmp/test_tmp_write_$$ && echo "write local tmp ok" || echo "write local tmp failed"
rm -f tmp/test_tmp_write_$$ || true

echo "========================="
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

# analysis_type: "both" (default), "degree", or "edge"
ANALYSIS_TYPE="${1:-both}"
RESULT_DIR="/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_DGLM_sensitivity_sex"
VIOLIN_DIR="${RESULT_DIR}/VIOLIN_Followup"
echo "Analysis type: ${ANALYSIS_TYPE}"
echo "DGLM result dir: ${RESULT_DIR}"
echo "VIOLIN follow-up dir: ${VIOLIN_DIR}"

echo "===== STEP 1/3: Run DGLM model ====="
echo "Output dir (STEP 1/3): ${RESULT_DIR}"
Rscript dglm_dk318_degree_edge_sensitivity_sex.R "${ANALYSIS_TYPE}"

echo "===== STEP 1/3 DONE ====="
echo "Output dir (STEP 1/3): ${RESULT_DIR}"

if [[ "${ANALYSIS_TYPE}" == "both" || "${ANALYSIS_TYPE}" == "degree" ]]; then
  echo "===== STEP 2/3: Generate DGLM degree brainmaps ====="
  echo "Output dir (STEP 2/3): ${RESULT_DIR}"
  python /data/home/tqi/data1/share/after_freesurfer/CODE/degree/degree_analysis_DK318_DGLM.py \
    --result-dir "${RESULT_DIR}"

  echo "===== STEP 2/3 DONE ====="
  echo "Output dir (STEP 2/3): ${RESULT_DIR}"
  echo "===== STEP 3/3: Run VIOLIN follow-up plots ====="
  echo "Output dir (STEP 3/3): ${VIOLIN_DIR}"
  bash /data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_violin_dk318.sh sensitivity_sex
  echo "===== STEP 3/3 DONE ====="
  echo "Output dir (STEP 3/3): ${VIOLIN_DIR}"
fi

echo "End: $(date)"
