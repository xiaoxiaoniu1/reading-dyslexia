#!/bin/bash
#SBATCH --job-name=DK318_DGLM
#SBATCH --output=logs/DK318_DGLM_%j.out
#SBATCH --error=logs/DK318_DGLM_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""


set -euo pipefail

mkdir -p logs
echo "==== TEMP DEBUG INFO ===="
echo "Host: $(hostname)"
echo "PWD: $(pwd)"
echo "TMPDIR(before): ${TMPDIR:-<unset>}"

echo "--- env tmp vars ---"
env | grep -E 'TMP|TEMP'

echo "--- system tmp status ---"
ls -ld /tmp || true
df -h /tmp || true
df -i /tmp || true

echo "--- current dir tmp status ---"
mkdir -p tmp
ls -ld tmp
df -h .
df -i .

echo "--- write test /tmp ---"
touch /tmp/test_tmp_write_$$ && echo "write /tmp ok" || echo "write /tmp failed"
rm -f /tmp/test_tmp_write_$$ || true

echo "--- write test local tmp ---"
touch tmp/test_tmp_write_$$ && echo "write local tmp ok" || echo "write local tmp failed"
rm -f tmp/test_tmp_write_$$ || true

echo "========================="
echo "Host: $(hostname)"
echo "Start: $(date)"

# 激活 conda 环境（临时关闭 -u 以避免 conda 激活脚本的未定义变量错误）
set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

# 进入代码目录
cd /data/home/tqi/data1/share/after_freesurfer/CODE/analysis

# 设置 R 临时目录，避免节点默认临时目录不可写
mkdir -p tmp
export TMPDIR="$(pwd)/tmp"

echo "TMPDIR: ${TMPDIR}"

# 跑 DGLM R 脚本
Rscript dglm_dk318_degree_edge.R

echo "End: $(date)"
