#!/bin/bash
#SBATCH -J MIND_all_subjects         # 作业名
#SBATCH -p partition_1              # 队列名：按你们集群实际改
#SBATCH -N 1
#SBATCH --cpus-per-task=8   # 用几个核
#SBATCH --mem=32G           # 内存
#SBATCH -t 72:00:00    # 最长运行时间
#SBATCH -o logs/mind_all_%j.out
#SBATCH -e logs/mind_all_%j.err

# 确保有日志目录
mkdir -p logs

# 1) 加载 freesurfer（按你们讲义）：
export FREESURFER_HOME=/data/software/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# 2) 激活 conda 环境
source /data/software/miniconda/bin/activate rd_env

# 3) 切到代码目录
cd /data1/tqi/share/after_freesurfer/CODE

# 4) 跑你刚才那个 python 脚本（不加 --sub，默认全体被试）
python run_mind_all_subjects.py
