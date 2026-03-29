#!/bin/bash
#SBATCH --job-name=DK68_STRAT_DD_TD
#SBATCH --output=analysis/logs/t_test_stratified_dk68_%j.out
#SBATCH --error=analysis/logs/t_test_stratified_dk68_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""

set -euo pipefail

mkdir -p logs
echo "Host: $(hostname)"
echo "Start: $(date)"

source ~/.bashrc

cd /data/home/tqi/data1/share/after_freesurfer/CODE/analysis

conda run -n rd_env python t_test_stratified_dk68.py

echo "End: $(date)"
