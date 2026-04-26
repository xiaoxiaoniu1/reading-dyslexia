#!/bin/bash
#SBATCH --job-name=zscore_viz
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs/zscore_viz_%j.log
#SBATCH --error=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs/zscore_viz_%j.err

source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r

SCRIPT_DIR="/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model"
R_SCRIPT="Visualize_Zscore_by_type.R"

echo "========================================" 
echo "Starting Z-score Visualization (FDR by type)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "========================================" 

cd "${SCRIPT_DIR}"
Rscript "${R_SCRIPT}"

EXIT_CODE=$?

echo "========================================" 
echo "End time: $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "========================================" 

exit ${EXIT_CODE}
