#!/bin/bash
#SBATCH --job-name=plot_DK318_NBS
#SBATCH --output=logs/plot_DK318_NBS_%j.out
#SBATCH --error=logs/plot_DK318_NBS_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""
set -euo pipefail

cd /data/home/tqi/data1/share/after_freesurfer/CODE/NBS

echo "Start: $(date)"

NBS_OUT_DIR="/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NBS_100"
TABLE_DIR="${NBS_OUT_DIR}/tables"
ROI_META_FILE="${NBS_OUT_DIR}/roi_metadata_with_mni_coordinates.csv"
EDGE_TABLE_1="${TABLE_DIR}/NBR_ALL_edges_in_FWE_significant_components.csv"
EDGE_TABLE_2="${TABLE_DIR}/NBR_ALL_suprathreshold_component_edges.csv"

if [[ ! -d "$NBS_OUT_DIR" ]]; then
  echo "ERROR: NBS output directory does not exist: $NBS_OUT_DIR" >&2
  exit 1
fi

if [[ ! -d "$TABLE_DIR" ]]; then
  echo "ERROR: NBS tables directory does not exist: $TABLE_DIR" >&2
  exit 1
fi

if [[ ! -f "$ROI_META_FILE" ]]; then
  echo "ERROR: ROI metadata file not found: $ROI_META_FILE" >&2
  exit 1
fi

if [[ ! -f "$EDGE_TABLE_1" && ! -f "$EDGE_TABLE_2" ]]; then
  echo "ERROR: No NBS edge table found. Checked:" >&2
  echo "  $EDGE_TABLE_1" >&2
  echo "  $EDGE_TABLE_2" >&2
  exit 1
fi

if [[ -f "$EDGE_TABLE_1" ]]; then
  echo "Using edge table: $EDGE_TABLE_1"
else
  echo "Using edge table: $EDGE_TABLE_2"
fi

echo "Using ROI metadata: $ROI_META_FILE"

set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env
set -u

python connect_MIND_DK318_NBS_components.py --nbs-dir "$NBS_OUT_DIR"

echo "End: $(date)"
