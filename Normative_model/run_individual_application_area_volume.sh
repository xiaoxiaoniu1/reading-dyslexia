#!/bin/bash
#SBATCH --job-name=individual_app_area_volume
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs/individual_application_area_volume_%j.log
#SBATCH --error=/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs/individual_application_area_volume_%j.err

set -euo pipefail

mkdir -p /data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs

echo "========================================"
echo "Starting Individual Application Analysis"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "R Environment: rd_env_r"
echo "R Path: $(which R || true)"
echo "Start time: $(date)"
echo "Script: Individual-application-normative-model.R"
echo "========================================"

set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

SCRIPT_DIR="/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model"
R_SCRIPT="Individual-application-normative-model.R"
SETUP_FLAG="${HOME}/.R/packages_installed_rd_env.flag"

if [ ! -f "$SETUP_FLAG" ]; then
  echo "First run detected. Installing R packages..."
  Rscript - <<'EOF'
required_packages <- c('gamlss', 'readxl', 'ggplot2', 'dplyr', 'tidyr', 'stringr')
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse=", "), "\n")
  install.packages(missing_packages, dependencies=TRUE)
}

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

  mkdir -p "$(dirname "$SETUP_FLAG")"
  touch "$SETUP_FLAG"
  echo "Package installation complete."
else
  echo "Packages already installed. Skipping installation step."
fi

cd "${SCRIPT_DIR}"
Rscript "${R_SCRIPT}"

EXIT_CODE=$?

echo "========================================"
echo "End time: $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "========================================"

exit ${EXIT_CODE}
