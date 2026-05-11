#!/bin/bash
#SBATCH --job-name=DK318_NMF_DD
#SBATCH --output=logs/DK318_NMF_DD_%j.out
#SBATCH --error=logs/DK318_NMF_DD_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""

set -euo pipefail

OUTPUT_DIR=${DD_NMF_OUT_DIR:-"/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NMF_DD_only"}
NMF_K_FINAL=${NMF_K_FINAL:-4}
NMF_NRUN=${NMF_NRUN:-100}
NMF_SEED=${NMF_SEED:-2026}
NMF_METHOD=${NMF_METHOD:-brunet}

mkdir -p logs

echo "=========================================="
echo "DK318 MIND Degree NMF - DD Only"
echo "=========================================="
echo "Host: $(hostname)"
echo "Start: $(date)"
echo "Output directory: ${OUTPUT_DIR}"
echo "NMF_K_FINAL: ${NMF_K_FINAL}"
echo "NMF_NRUN: ${NMF_NRUN}"
echo "NMF_SEED: ${NMF_SEED}"
echo "NMF_METHOD: ${NMF_METHOD}"
echo ""

set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

echo "Conda environment: rd_env_r"
echo "R version: $(Rscript --version 2>&1 | head -n 1)"
echo ""

cd /data/home/tqi/data1/share/after_freesurfer/CODE/NMF
mkdir -p tmp
export TMPDIR="$(pwd)/tmp"
export DD_NMF_OUT_DIR="${OUTPUT_DIR}"
export NMF_K_FINAL
export NMF_NRUN
export NMF_SEED
export NMF_METHOD

echo "Working directory: $(pwd)"
echo "TMPDIR: ${TMPDIR}"
echo ""

echo "Checking required R packages..."
Rscript - <<'EOF'
required_packages <- c("readxl", "dplyr", "NMF", "car", "ggplot2", "reshape2", "readr", "tidyr", "pheatmap", "emmeans", "effectsize")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  cat("Missing packages:", paste(missing_packages, collapse = ", "), "\n")
  if ("NMF" %in% missing_packages) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cran.r-project.org")
    if (!requireNamespace("Biobase", quietly = TRUE)) BiocManager::install("Biobase", update = FALSE, ask = FALSE)
  }
  for (pkg in missing_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cran.r-project.org", dependencies = TRUE)
  }
}
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package still missing: ", pkg)
}
cat("All required packages are available.\n")
EOF

echo ""
echo "Running DD-only NMF analysis..."
echo "Start: $(date)"
Rscript 05_nmf_mind_degree_DD_only.R
echo "End: $(date)"
echo ""

echo "Running DD-only post-hoc/descriptive analysis..."
echo "Start: $(date)"
if Rscript 08_nmf_posthoc_DD_only.R; then
  echo "DD-only post-hoc/descriptive analysis completed successfully."
else
  echo "WARNING: DD-only post-hoc/descriptive analysis failed. Continuing..."
fi
echo "End: $(date)"
echo ""

echo "Running DD-only visualization..."
echo "Start: $(date)"
if Rscript 07_nmf_visualization_DD_only.R; then
  echo "DD-only visualization completed successfully."
else
  echo "WARNING: DD-only visualization failed. Continuing..."
fi
echo "End: $(date)"
echo ""

echo "Running DD-only network enrichment analysis..."
echo "Start: $(date)"
if Rscript 06_nmf_network_enrichment_DD_only.R; then
  echo "DD-only network enrichment completed successfully."
else
  echo "WARNING: DD-only network enrichment failed. Main NMF outputs are still available."
fi
echo "End: $(date)"
echo ""

echo "Checking key output files:"
KEY_FILES=(
  "NMF_DD_only_K_selection_metrics.csv"
  "NMF_DD_only_W_subject_loadings_and_subtypes.csv"
  "NMF_DD_only_subtype_summary.csv"
  "NMF_DD_only_H_component_roi_weights.csv"
  "NMF_DD_only_component_top_ROIs.csv"
  "NMF_DD_only_ROI_hard_assignment.csv"
  "NMF_DD_only_posthoc_age_sex_contrasts.csv"
  "NMF_DD_only_descriptive_statistics_age_sex.csv"
  "NMF_DD_only_descriptive_statistics_by_subtype.csv"
  "plots/DD_only_all_components_overview.png"
  "plots/DD_only_component_scores_heatmap.png"
  "plots/DD_only_dominant_subtype_counts_by_age.png"
  "NMF_DD_only_network_enrichment_hard_assignment.csv"
  "NMF_DD_only_network_enrichment_top_rois.csv"
  "NMF_DD_only_degree_QC_summary.csv"
)

for file in "${KEY_FILES[@]}"; do
  if [ -f "${OUTPUT_DIR}/${file}" ]; then
    echo "  ✓ ${file}"
  else
    echo "  ✗ ${file} (missing)"
  fi
done

echo ""
if [ -d "${OUTPUT_DIR}" ]; then
  shopt -s nullglob
  CSV_FILES=("${OUTPUT_DIR}"/*.csv)
  PNG_FILES=("${OUTPUT_DIR}"/*.png)
  CSV_COUNT=${#CSV_FILES[@]}
  PNG_COUNT=${#PNG_FILES[@]}
  shopt -u nullglob
  echo "Output files generated:"
  echo "  CSV files: ${CSV_COUNT}"
  echo "  PNG files: ${PNG_COUNT}"
fi

echo ""
echo "=========================================="
echo "DD-only NMF Pipeline Complete"
echo "=========================================="
echo "End: $(date)"

rm -rf tmp
