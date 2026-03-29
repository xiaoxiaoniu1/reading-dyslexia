#!/bin/bash
#SBATCH --job-name=DK_normative_validation_newsite
#SBATCH --output=Normative_model/logs/DK_normative_validation_newsite_%j.out
#SBATCH --error=Normative_model/logs/DK_normative_validation_newsite_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --account=""

set -euo pipefail

mkdir -p /data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/logs
echo "Host: $(hostname)"
echo "Start: $(date)"

set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

cd /data/home/tqi/data1/share/after_freesurfer/CODE

# Set a writable temp directory to avoid 'cannot create R_TempDir' error
MY_TMPDIR=/data/home/tqi/data1/share/after_freesurfer/tmp_R_$$
mkdir -p "$MY_TMPDIR"
export TMPDIR="$MY_TMPDIR"
export TMP="$MY_TMPDIR"
export TEMP="$MY_TMPDIR"

Rscript /data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/Validate-calibration-newsite.R

# Clean up temp dir
rm -rf "$MY_TMPDIR"

echo "End: $(date)"
