#!/bin/bash
#SBATCH --job-name=DK68_Interaction_Trajectories
#SBATCH --output=logs/DK68_Interaction_Trajectories_%j.out
#SBATCH --error=logs/DK68_Interaction_Trajectories_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""
set -euo pipefail
mkdir -p logs
cd /data/home/tqi/data1/share/after_freesurfer/CODE/degree

echo "Start: $(date)"


# 激活 Python 环境
set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

# 运行 Python 脚本
Rscript DK68_Interaction_Trajectories.R

echo "End: $(date)"
