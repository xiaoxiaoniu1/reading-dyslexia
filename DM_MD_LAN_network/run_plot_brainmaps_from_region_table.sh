#!/bin/bash
#SBATCH --job-name=DK318_region_table_brainmaps
#SBATCH --output=logs/DK318_region_table_brainmaps_%j.out
#SBATCH --error=logs/DK318_region_table_brainmaps_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""
set -euo pipefail

cd "/data/home/tqi/data1/share/after_freesurfer/CODE/DM_MD_LAN_network"

echo "Start: $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"

set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env
set -u

python "/data/home/tqi/data1/share/after_freesurfer/CODE/DM_MD_LAN_network/plot_brainmaps_from_region_table.py" \
  --csvs \
  "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/brain_regions_at_least_one_task_mean_gt_2.csv" \
  "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/brain_regions_at_least_two_tasks_mean_gt_2.csv" \
  "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/alphabetic_structural.csv" \
  "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/morphosyllabic_functional.csv" \
  "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/morphosyllabic_structural.csv" \
  "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/alphabetic_functional.csv" \
  --output-dir "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/brain_region_table_brainmaps"

echo "End: $(date)"
