#!/bin/bash
#SBATCH --job-name=DK318_Tmap
#SBATCH --output=logs/DK318_Tmap_%j.out
#SBATCH --error=logs/DK318_Tmap_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""

# ============================================================
# DK318 MIND Degree T-map 完整分析流程
# 依次运行: DGLM → 可视化 → 聚类
# ============================================================

set -euo pipefail

echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                              ║"
echo "║              DK318 MIND Degree T-map 完整分析流程                            ║"
echo "║                                                                              ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"
echo ""

# 创建日志目录
mkdir -p logs

echo "Host: $(hostname)"
echo "Start: $(date)"
echo ""

# 激活 conda 环境（临时关闭 -u 以避免 conda 激活脚本的未定义变量错误）
set +u
source /data/software/miniconda/etc/profile.d/conda.sh
conda activate rd_env_r
set -u

# 设置工作目录
SCRIPT_DIR="/data/home/tqi/data1/share/after_freesurfer/CODE/Clustering"
cd "$SCRIPT_DIR" || exit 1

echo "工作目录: $SCRIPT_DIR"
echo ""

# 设置 R 临时目录，避免节点默认临时目录不可写
mkdir -p tmp
export TMPDIR="$(pwd)/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

echo "TMPDIR: ${TMPDIR}"
echo ""

# 检查 R 是否可用
if ! command -v Rscript &> /dev/null; then
    echo "❌ 错误: Rscript 未找到"
    echo "   请确保 R 已安装并在 PATH 中"
    exit 1
fi

echo "✓ R 版本:"
Rscript --version
echo ""

# ============================================================
# 步骤 0: 可选测试（仅在交互模式下运行）
# ============================================================
if [ -t 0 ]; then
    read -p "是否运行逻辑测试? (y/n, 默认 n): " run_test
    run_test=${run_test:-n}
else
    # 非交互模式，跳过测试
    run_test="n"
fi

if [[ "$run_test" == "y" || "$run_test" == "Y" ]]; then
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "步骤 0: 逻辑测试"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    
    Rscript test_logic.R
    
    if [ $? -ne 0 ]; then
        echo ""
        echo "❌ 逻辑测试失败"
        exit 1
    else
        echo ""
        echo "✓ 逻辑测试通过"
    fi
    echo ""
fi

# ============================================================
# 步骤 1: DGLM T-map 提取
# ============================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "步骤 1: DGLM T-map 提取"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "运行: 01_dglm_degree_tmap_TD_minus_DD.R"
echo "预计时间: 5-15 分钟"
echo ""

START_TIME=$(date +%s)

Rscript 01_dglm_degree_tmap_TD_minus_DD.R

if [ $? -ne 0 ]; then
    echo ""
    echo "❌ DGLM 分析失败"
    exit 1
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo ""
echo "✓ DGLM 分析完成 (用时: ${ELAPSED}秒)"
echo ""

# 检查输出文件
TMAP_DIR="/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Clustering_5"
FULL_CSV="$TMAP_DIR/DGLM_DK318_degree_Tmap_TD_minus_DD_full.csv"
TMAP_CSV="$TMAP_DIR/DGLM_DK318_degree_Tmap_TD_minus_DD_318x2.csv"
SUMMARY_CSV="$TMAP_DIR/DGLM_DK318_degree_Tmap_TD_minus_DD_summary.csv"

echo "检查输出文件:"
if [ -f "$FULL_CSV" ]; then
    echo "  ✓ 完整结果: $(wc -l < "$FULL_CSV") 行"
else
    echo "  ❌ 完整结果文件未找到"
    exit 1
fi

if [ -f "$TMAP_CSV" ]; then
    echo "  ✓ 318×2 矩阵: $(wc -l < "$TMAP_CSV") 行"
else
    echo "  ❌ 318×2 矩阵文件未找到"
    exit 1
fi

if [ -f "$SUMMARY_CSV" ]; then
    echo "  ✓ 汇总统计"
else
    echo "  ❌ 汇总统计文件未找到"
    exit 1
fi

echo ""

# ============================================================
# 步骤 2: 可视化（自动运行）
# ============================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "步骤 2: T-map 可视化"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "运行: 02_visualize_tmap.R"
echo ""

Rscript 02_visualize_tmap.R

if [ $? -ne 0 ]; then
    echo ""
    echo "⚠ 可视化失败，但继续流程"
else
    echo ""
    echo "✓ 可视化完成"
fi
echo ""

# ============================================================
# 步骤 3: 聚类分析（自动运行）
# ============================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "步骤 3: 聚类分析"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "运行: 03_cluster_tmap.R"
echo ""

Rscript 03_cluster_tmap.R

if [ $? -ne 0 ]; then
    echo ""
    echo "⚠ 聚类分析失败"
else
    echo ""
    echo "✓ 聚类分析完成"
fi
echo ""

# ============================================================
# 步骤 4: 网络富集分析与可视化（自动运行）
# ============================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "步骤 4: 网络富集分析与可视化"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "运行: 04_cluster_network_enrichment_and_visualization.R"
echo ""

Rscript 04_cluster_network_enrichment_and_visualization.R

if [ $? -ne 0 ]; then
    echo ""
    echo "⚠ 网络富集分析失败"
else
    echo ""
    echo "✓ 网络富集分析完成"
fi
echo ""

# ============================================================
# 完成总结
# ============================================================
echo ""
echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                              ║"
echo "║                          分析流程完成                                        ║"
echo "║                                                                              ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"
echo ""

echo "输出文件位置:"
echo ""
echo "1. T-map 结果:"
echo "   $TMAP_DIR/"
echo "   ├── DGLM_DK318_degree_Tmap_TD_minus_DD_full.csv"
echo "   ├── DGLM_DK318_degree_Tmap_TD_minus_DD_318x2.csv"
echo "   └── DGLM_DK318_degree_Tmap_TD_minus_DD_summary.csv"
echo ""
echo "2. 可视化图表:"
echo "   $TMAP_DIR/figures/"
echo ""
echo "3. 聚类结果:"
echo "   $TMAP_DIR/clustering/"
echo "   ├── degree_Tmap_cluster_K3.csv"
echo "   ├── degree_Tmap_cluster_K3_centroids.csv"
echo "   ├── degree_Tmap_cluster_K3_summary.csv"
echo "   └── k_selection_metrics.csv"
echo ""
echo "4. 网络富集分析:"
echo "   $TMAP_DIR/clustering/"
echo "   ├── cluster_network_enrichment.csv"
echo "   ├── cluster_network_enrichment_with_pattern.csv"
echo "   ├── cluster_network_overlap_details.csv"
echo "   ├── enrichment_summary_table.csv"
echo "   ├── enrichment_summary_report.txt"
echo "   ├── enrichment_heatmap_enhanced.png"
echo "   ├── enrichment_barplot_by_cluster.png"
echo "   └── network_distribution_pie_charts.png"
echo ""

echo "查看汇总统计:"
echo "  cat $SUMMARY_CSV"
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "End: $(date)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
