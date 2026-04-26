#!/bin/bash
#SBATCH --job-name=DK318_VIOLIN_main
#SBATCH --output=analysis/logs/DK318_VIOLIN_main_%j.out
#SBATCH --error=analysis/logs/DK318_VIOLIN_main_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""

set -euo pipefail
bash /data/home/tqi/data1/share/after_freesurfer/CODE/analysis/run_violin_dk318.sh main
