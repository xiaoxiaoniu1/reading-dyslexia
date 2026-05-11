#!/bin/bash
#SBATCH --job-name=DK318_NMF
#SBATCH --output=logs/DK318_NMF_%j.out
#SBATCH --error=logs/DK318_NMF_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""

set -euo pipefail

# ============================================================
# Configuration Options
# ============================================================

# Set to "true" to skip steps that have already completed
# Set to "false" to rerun all steps
SKIP_COMPLETED=${SKIP_COMPLETED:-true}

# Output directory
OUTPUT_DIR="/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NMF_3"

# ============================================================

# 创建日志目录
mkdir -p logs

echo "=========================================="
echo "DK318 MIND Degree NMF Analysis Pipeline"
echo "=========================================="
echo "Host: $(hostname)"
echo "Start: $(date)"
echo ""
echo "Configuration:"
echo "  SKIP_COMPLETED: $SKIP_COMPLETED"
echo ""

# 激活 conda 环境（临时关闭 -u 以避免 conda 激活脚本的未定义变量错误）
set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

echo "Conda environment: rd_env_r"
echo "R version: $(Rscript --version 2>&1 | head -n 1)"
echo ""

# 进入代码目录
cd /data/home/tqi/data1/share/after_freesurfer/CODE/NMF

# 设置 R 临时目录，避免节点默认临时目录不可写
mkdir -p tmp
export TMPDIR="$(pwd)/tmp"

echo "Working directory: $(pwd)"
echo "TMPDIR: ${TMPDIR}"
echo ""

# ==========================================
# Step 0: Check and Install R Packages
# ==========================================

echo "=========================================="
echo "Step 0: Check and Install R Packages"
echo "=========================================="
echo "Start: $(date)"

Rscript - <<'EOF'
# Required packages
required_packages <- list(
  core = c("readxl", "dplyr", "stringr"),
  nmf = c("NMF"),
  stats = c("broom", "car", "emmeans", "effectsize"),
  viz = c("ggplot2", "reshape2", "RColorBrewer", "pheatmap", "gridExtra")
)

all_packages <- unlist(required_packages)

cat("\nChecking R packages...\n")
cat("----------------------------------------\n")

missing_packages <- c()

for (pkg in all_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  ✗", pkg, "(missing)\n")
    missing_packages <- c(missing_packages, pkg)
  } else {
    cat("  ✓", pkg, "\n")
  }
}

if (length(missing_packages) > 0) {
  cat("\n")
  cat("Found", length(missing_packages), "missing package(s):", paste(missing_packages, collapse = ", "), "\n")
  cat("Installing missing packages...\n")
  cat("This may take a few minutes...\n\n")
  
  # Check if NMF is in missing packages - it needs Bioconductor dependency
  if ("NMF" %in% missing_packages) {
    cat("----------------------------------------\n")
    cat("NMF requires Bioconductor packages.\n")
    
    # Install BiocManager if not available
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      cat("Installing BiocManager...\n")
      install.packages("BiocManager", repos = "https://cran.r-project.org", quiet = FALSE)
      cat("  ✓ BiocManager installed\n")
    } else {
      cat("  ✓ BiocManager already installed\n")
    }
    
    # Check and install Biobase from Bioconductor
    if (!requireNamespace("Biobase", quietly = TRUE)) {
      cat("Installing Biobase from Bioconductor...\n")
      BiocManager::install("Biobase", update = FALSE, ask = FALSE)
      cat("  ✓ Biobase installed\n")
    } else {
      cat("  ✓ Biobase already installed\n")
    }
    cat("----------------------------------------\n\n")
  }
  
  # Install missing packages one by one
  for (pkg in missing_packages) {
    # Double check if package is still missing (might have been installed as dependency)
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Installing", pkg, "...\n")
      tryCatch({
        install.packages(pkg, repos = "https://cran.r-project.org", dependencies = TRUE, quiet = FALSE)
        
        # Verify installation
        if (requireNamespace(pkg, quietly = TRUE)) {
          cat("  ✓", pkg, "installed successfully\n")
        } else {
          cat("  ✗", pkg, "installation failed\n")
        }
      }, error = function(e) {
        cat("  ✗ ERROR installing", pkg, ":", conditionMessage(e), "\n")
      })
    } else {
      cat("  ✓", pkg, "already installed (installed as dependency)\n")
    }
  }
  
  cat("\n")
  cat("Package installation complete.\n")
  
  # Final verification
  cat("\nFinal verification:\n")
  cat("----------------------------------------\n")
  all_ok <- TRUE
  for (pkg in all_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat("  ✓", pkg, "\n")
    } else {
      cat("  ✗", pkg, "(STILL MISSING)\n")
      all_ok <- FALSE
    }
  }
  
  if (all_ok) {
    cat("\n✓ All packages are ready!\n")
  } else {
    cat("\n✗ Some packages are still missing. Analysis may fail.\n")
  }
  
} else {
  cat("\n")
  cat("✓ All required packages are already installed.\n")
}

cat("----------------------------------------\n")
EOF

echo "End: $(date)"
echo ""

# ==========================================
# Step 1: Main NMF Analysis
# ==========================================

echo "=========================================="
echo "Step 1: Main NMF Analysis"
echo "=========================================="

# Check if main analysis is already completed
MAIN_DONE=false
if [ "$SKIP_COMPLETED" = "true" ] && [ -f "${OUTPUT_DIR}/NMF_component_lm_anova_results.csv" ]; then
    echo "✓ Main NMF analysis already completed. Skipping..."
    MAIN_DONE=true
else
    echo "Start: $(date)"
    
    Rscript 01_nmf_mind_degree.R
    
    if [ $? -eq 0 ]; then
        echo "Main NMF analysis completed successfully."
        MAIN_DONE=true
    else
        echo "ERROR: Main NMF analysis failed!"
        exit 1
    fi
    
    echo "End: $(date)"
fi

echo ""

# ==========================================
# Step 2: Network Enrichment Analysis
# ==========================================

echo "=========================================="
echo "Step 2: Network Enrichment Analysis"
echo "=========================================="

# Check if network enrichment is already completed
if [ "$SKIP_COMPLETED" = "true" ] && [ -f "${OUTPUT_DIR}/NMF_network_enrichment_hard_assignment.csv" ]; then
    echo "✓ Network enrichment analysis already completed. Skipping..."
else
    echo "Start: $(date)"
    
    Rscript 02_nmf_network_enrichment.R
    
    if [ $? -eq 0 ]; then
        echo "Network enrichment analysis completed successfully."
    else
        echo "WARNING: Network enrichment analysis failed. Continuing..."
    fi
    
    echo "End: $(date)"
fi

echo ""

# ==========================================
# Step 3: Enhanced Visualization
# ==========================================

echo "=========================================="
echo "Step 3: Enhanced Visualization"
echo "=========================================="

# Check if visualization is already completed
if [ "$SKIP_COMPLETED" = "true" ] && [ -d "${OUTPUT_DIR}/plots" ] && [ -f "${OUTPUT_DIR}/plots/all_components_overview.png" ]; then
    echo "✓ Visualization already completed. Skipping..."
else
    echo "Start: $(date)"
    
    Rscript 03_nmf_visualization.R
    
    if [ $? -eq 0 ]; then
        echo "Visualization completed successfully."
    else
        echo "WARNING: Visualization failed. Continuing..."
    fi
    
    echo "End: $(date)"
fi

echo ""

# ==========================================
# Step 4: Post-hoc Analysis
# ==========================================

echo "=========================================="
echo "Step 4: Post-hoc Analysis"
echo "=========================================="

# Check if post-hoc analysis is already completed
if [ "$SKIP_COMPLETED" = "true" ] && [ -f "${OUTPUT_DIR}/NMF_posthoc_pairwise_contrasts.csv" ]; then
    echo "✓ Post-hoc analysis already completed. Skipping..."
else
    echo "Start: $(date)"
    
    Rscript 04_nmf_posthoc.R
    
    if [ $? -eq 0 ]; then
        echo "Post-hoc analysis completed successfully."
    else
        echo "WARNING: Post-hoc analysis failed. Continuing..."
    fi
    
    echo "End: $(date)"
fi

echo ""

# ==========================================
# Summary
# ==========================================

echo "=========================================="
echo "Pipeline Summary"
echo "=========================================="

OUTPUT_DIR="/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NMF"

echo "Output directory: ${OUTPUT_DIR}"
echo ""

# 检查关键输出文件
echo "Checking key output files:"

KEY_FILES=(
    "NMF_W_subject_loadings.csv"
    "NMF_H_component_roi_weights.csv"
    "NMF_component_lm_anova_results.csv"
    "NMF_component_top_ROIs.csv"
    "NMF_ROI_hard_assignment.csv"
)

for file in "${KEY_FILES[@]}"; do
    if [ -f "${OUTPUT_DIR}/${file}" ]; then
        echo "  ✓ ${file}"
    else
        echo "  ✗ ${file} (missing)"
    fi
done

echo ""

# 统计输出文件数量
if [ -d "${OUTPUT_DIR}" ]; then
    CSV_COUNT=$(find "${OUTPUT_DIR}" -maxdepth 1 -name "*.csv" | wc -l)
    PNG_COUNT=$(find "${OUTPUT_DIR}" -maxdepth 1 -name "*.png" | wc -l)
    PLOT_COUNT=$(find "${OUTPUT_DIR}/plots" -name "*.png" 2>/dev/null | wc -l)
    
    echo "Output files generated:"
    echo "  CSV files: ${CSV_COUNT}"
    echo "  PNG files (main): ${PNG_COUNT}"
    echo "  PNG files (plots): ${PLOT_COUNT}"
    echo "  Total: $((CSV_COUNT + PNG_COUNT + PLOT_COUNT))"
else
    echo "WARNING: Output directory not found!"
fi

echo ""
echo "=========================================="
echo "Pipeline Complete"
echo "=========================================="
echo "End: $(date)"
echo ""

# 清理临时目录
if [ -d "tmp" ]; then
    echo "Cleaning up temporary directory..."
    rm -rf tmp
fi

echo "All done!"
