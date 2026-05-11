#!/bin/bash
# ============================================================
# DK318 3-Network Degree Analysis - SLURM Submission Script
# Purpose: Run statistical analysis on DM/MD/Reading network degrees
# ============================================================

#SBATCH --job-name=DK318_3NET_ANALYSIS
#SBATCH --output=logs/DK318_3NET_ANALYSIS_%j.out
#SBATCH --error=logs/DK318_3NET_ANALYSIS_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# ============================================================
# 1. Print Job Information
# ============================================================

echo "========================================"
echo "DK318 3-Network Degree Analysis"
echo "========================================"
echo "Host: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Start time: $(date)"
echo ""

# ============================================================
# 2. Activate Conda Environment
# ============================================================

echo "Activating conda environment..."
set +u  # Temporarily disable undefined variable check for conda
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u
echo "Conda environment: ${CONDA_DEFAULT_ENV}"
echo ""

# ============================================================
# 3. Set Working Directory and Paths
# ============================================================

CODE_DIR="/data/home/tqi/data1/share/after_freesurfer/CODE/Network"
ANALYSIS_SCRIPT="${CODE_DIR}/01_DK318_3network_degree_analysis.R"

cd "${CODE_DIR}"
echo "Working directory: $(pwd)"

# Create necessary directories
mkdir -p logs tmp
export TMPDIR="$(pwd)/tmp"
echo "TMPDIR: ${TMPDIR}"
echo ""

# ============================================================
# 4. Verify Script Exists
# ============================================================

if [[ ! -f "${ANALYSIS_SCRIPT}" ]]; then
  echo "ERROR: Analysis script not found: ${ANALYSIS_SCRIPT}"
  exit 1
fi

echo "Analysis script: ${ANALYSIS_SCRIPT}"
echo ""

# ============================================================
# 5. Run Analysis Script
# ============================================================

echo "Starting R analysis..."
echo "----------------------------------------"

Rscript "${ANALYSIS_SCRIPT}"

EXIT_CODE=$?

echo "----------------------------------------"
echo ""

# ============================================================
# 6. Check Exit Status
# ============================================================

if [[ ${EXIT_CODE} -eq 0 ]]; then
  echo "✓ Analysis completed successfully"
else
  echo "✗ Analysis failed with exit code: ${EXIT_CODE}"
  exit ${EXIT_CODE}
fi

# ============================================================
# 7. Summary
# ============================================================

echo ""
echo "========================================"
echo "Job Summary"
echo "========================================"
echo "End time: $(date)"
echo "Output directory: /data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_3network_analysis/analysis"
echo ""
echo "Next step: Run visualization script"
echo "  sbatch ${CODE_DIR}/run_DK318_3network_degree_visualization.sh"
echo "========================================"
