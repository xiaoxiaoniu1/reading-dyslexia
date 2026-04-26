#!/bin/bash
#SBATCH --job-name=disease_app
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs/disease_application_%j.log
#SBATCH --error=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs/disease_application_%j.err

# Shell script to run Disease-application-normative-model.R with SLURM
# This script applies calibrated normative models to disease patients

# Activate conda environment
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r

# Set variables
SCRIPT_DIR="/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model"
R_SCRIPT="Disease-application-normative-model.R"
SETUP_FLAG="${HOME}/.R/packages_installed_rd_env.flag"

# Print start information
echo "========================================" 
echo "Starting Disease Application Analysis"
echo "Job ID: ${SLURM_JOB_ID}"
echo "R Environment: rd_env_r"
echo "R Path: $(which R)"
echo "Start time: $(date)"
echo "Script: ${R_SCRIPT}"
echo "========================================" 

# Install required R packages only if not already done
if [ ! -f "$SETUP_FLAG" ]; then
  echo "First run detected. Installing R packages..."
  Rscript - <<'EOF'
required_packages <- c('gamlss', 'readxl', 'ggplot2', 'dplyr', 'tidyr', 'stringr')
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse=", "), "\n")
  install.packages(missing_packages, dependencies=TRUE)
}

# Verify all packages are installed and loadable
cat("\nVerifying packages...\n")
for(pkg in required_packages) {
  if(require(pkg, character.only=TRUE, quietly=TRUE)) {
    cat(pkg, "- OK\n")
  } else {
    cat(pkg, "- FAILED TO LOAD\n")
    stop(paste("Failed to load package:", pkg))
  }
}
EOF
  
  # Create flag file to indicate setup is complete
  mkdir -p "$(dirname "$SETUP_FLAG")"
  touch "$SETUP_FLAG"
  echo "Package installation complete."
else
  echo "Packages already installed. Skipping installation step."
fi

# Run the R script
cd "${SCRIPT_DIR}"
Rscript "${R_SCRIPT}"

# Capture exit code
EXIT_CODE=$?

# Print end information
echo "========================================" 
echo "End time: $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "========================================" 

exit ${EXIT_CODE}
