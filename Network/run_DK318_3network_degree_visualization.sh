#!/bin/bash
# ============================================================
# DK318 3-Network Degree Visualization - SLURM Submission Script
# Purpose: Generate figures from analysis results
# ============================================================

#SBATCH --job-name=DK318_3NET_VIZ
#SBATCH --output=logs/DK318_3NET_VIZ_%j.out
#SBATCH --error=logs/DK318_3NET_VIZ_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --account=""

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# ============================================================
# 1. Print Job Information
# ============================================================

echo "========================================"
echo "DK318 3-Network Degree Visualization"
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
VIZ_SCRIPT="${CODE_DIR}/02_DK318_3network_degree_visualization.R"
ANALYSIS_DIR="/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_3network_analysis/analysis"

cd "${CODE_DIR}"
echo "Working directory: $(pwd)"

# Create necessary directories
mkdir -p logs tmp
export TMPDIR="$(pwd)/tmp"
echo "TMPDIR: ${TMPDIR}"
echo ""

# ============================================================
# 4. Verify Script and Input Files Exist
# ============================================================

if [[ ! -f "${VIZ_SCRIPT}" ]]; then
  echo "ERROR: Visualization script not found: ${VIZ_SCRIPT}"
  exit 1
fi

if [[ ! -d "${ANALYSIS_DIR}" ]]; then
  echo "ERROR: Analysis directory not found: ${ANALYSIS_DIR}"
  echo "Please run analysis script first."
  exit 1
fi

echo "Visualization script: ${VIZ_SCRIPT}"
echo "Analysis directory: ${ANALYSIS_DIR}"
echo ""

# ============================================================
# 5. Run Visualization Script
# ============================================================

echo "Starting R visualization..."
echo "----------------------------------------"

Rscript "${VIZ_SCRIPT}"

EXIT_CODE=$?

echo "----------------------------------------"
echo ""

# ============================================================
# 6. Check Exit Status
# ============================================================

if [[ ${EXIT_CODE} -eq 0 ]]; then
  echo "✓ Visualization completed successfully"
else
  echo "✗ Visualization failed with exit code: ${EXIT_CODE}"
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
echo "Figure directory: /data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_3network_analysis/figures"
echo ""
echo "Generated figures:"
echo "  - fig_01_network_overlap_membership_bar"
echo "  - fig_02_mixed_model_adjusted_means_by_network"
echo "  - fig_03_raw_network_degree_distribution"
echo "  - fig_04_mixed_model_TD_minus_DD_t_heatmap"
echo "  - fig_05_mixed_model_TD_minus_DD_contrast_profile"
echo "  - fig_06_per_network_TD_minus_DD_t_heatmap"
echo "========================================"
