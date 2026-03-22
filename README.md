# Reading Dyslexia — MIND Network Analysis Pipeline

A full analysis pipeline for constructing **Morphometric INverse Distance (MIND)** cortical similarity networks, performing site-harmonization, and running group-level statistical analyses (ANOVA / ANCOVA) on children and adults with and without reading dyslexia.

---

## Table of Contents

1. [Background](#background)
2. [Repository Structure](#repository-structure)
3. [Pipeline Overview](#pipeline-overview)
4. [Step-by-step Usage](#step-by-step-usage)
   - [Step 1 — Compute MIND matrices](#step-1--compute-mind-matrices)
   - [Step 2 — ComBat harmonization](#step-2--combat-harmonization)
   - [Step 3 — Post-ComBat degree computation](#step-3--post-combat-degree-computation)
   - [Step 4 — Degree analysis (ANOVA / ANCOVA)](#step-4--degree-analysis-anova--ancova)
   - [Step 5 — Edge-level analysis](#step-5--edge-level-analysis)
   - [Step 6 — Brain map visualization](#step-6--brain-map-visualization)
5. [Parcellations](#parcellations)
6. [Dependencies](#dependencies)
7. [Data Format](#data-format)
8. [Notes](#notes)

---

## Background

MIND networks capture structural brain organization by computing, for each pair of cortical parcels, the similarity of their vertex-wise morphometric profiles across five features: cortical thickness (CT), myelination proxy (MC), volume (Vol), sulcal depth (SD), and surface area (SA). This repository applies MIND to a multi-site developmental cohort to investigate cortical network differences between typically developing (TD) readers and individuals with developmental dyslexia (DD), across children and adults.

---

## Repository Structure

```
CODE/
├── MIND/                          # MIND network construction
│   ├── MIND.py                    # Core compute_MIND() function
│   ├── get_vertex_df.py           # Vertex-level feature extraction
│   ├── register_and_vol2surf.py   # Volume-to-surface registration (BBRegister + mri_vol2surf)
│   └── requirements.txt           # Python dependencies
│
├── combat/                        # Site harmonization
│   ├── combat_dk318.py            # ComBat on DK-318 MIND edges
│   └── combat_dk68.py             # ComBat on DK-68 MIND edges
│
├── degree/                        # Node-level (degree) analyses
│   ├── degree_analysis_DK318.py                # ANOVA: Diagnosis × AgeGroup (DK-318)
│   ├── degree_analysis_DK318_site.py           # Site-effect check (DK-318)
│   ├── degree_analysis_DK318_t_test.py         # Post-hoc t-tests (DK-318)
│   ├── degree_analysis_DK318_adult_with_IQ.py  # ANCOVA with IQ covariate (adults, DK-318)
│   ├── degree_analysis_DK68.py                 # ANOVA (DK-68)
│   ├── degree_analysis_DK68_site.py            # Site-effect check (DK-68)
│   ├── degree_MIND_heatmap.py                  # Group-mean MIND matrix heatmap
│   ├── DK318_brain_region_check.py             # Sanity check: region coverage (DK-318)
│   ├── DK68_brain_region_check.py              # Sanity check: region coverage (DK-68)
│   ├── DK318_Diagnosis_degree_barplot.py       # Bar plots: degree by diagnosis (DK-318)
│   ├── DK68_Diagnosis_degree_barplot.py        # Bar plots: degree by diagnosis (DK-68)
│   ├── DK318_Interaction_Trajectories.R        # Developmental trajectory plots (DK-318)
│   ├── DK68_Interaction_Trajectories.R         # Developmental trajectory plots (DK-68)
│   └── run_*.sh                                # SLURM submission scripts
│
├── edge/                          # Edge-level analyses
│   ├── connect_ANOVA.py           # ANOVA on every edge (DK-318)
│   ├── connect_ANOVA_site.py      # Site-effect check at edge level (DK-318)
│   ├── connect_ANOVA_dk68.py      # ANOVA on every edge (DK-68)
│   ├── connect_ANOVA_dk68_site.py # Site-effect check at edge level (DK-68)
│   └── connect_t_test.py          # Post-hoc edge t-tests
│
├── analysis/                      # Downstream visualization & analysis
│
├── aftercombat_MIND_degree.py     # Post-ComBat: recompute degree from harmonized matrices
├── run_mind_all_subjects.py       # Batch MIND computation for all subjects
├── run_mind_all_subjects.sh       # SLURM wrapper for batch MIND computation
├── get_MNI.R                      # Extract MNI coordinates for parcels
└── README.md
```

---

## Pipeline Overview

```
FreeSurfer output
       │
       ▼
 [Step 1] Compute MIND matrices per subject (DK-318 or DK-68)
       │  run_mind_all_subjects.py
       │  → <subject>_MIND_DK318.npy  (318×318)
       │  → <subject>_MIND_DK318_edges.npy  (50403-dim upper triangle)
       │  → <subject>_MIND_DK318_degree.npy (318-dim)
       │
       ▼
 [Step 2] ComBat site harmonization on edge vectors
       │  combat/combat_dk318.py
       │  → <subject>_MIND_DK318_combat.npy  (318×318, harmonized)
       │
       ▼
 [Step 3] Recompute degree from harmonized matrices
       │  aftercombat_MIND_degree.py
       │  → <subject>_degree.csv
       │  → <subject>_MIND_DK318_combat_labeled.csv
       │
       ▼
 [Step 4] Node-level statistics (degree)
       │  degree/degree_analysis_DK318.py
       │  → ANOVA_DK318_degree_results.csv  (p-values, η², FDR)
       │
       ▼
 [Step 5] Edge-level statistics
       │  edge/connect_ANOVA.py
       │  → ANOVA_DK318_edge_results.csv
       │
       ▼
 [Step 6] Brain map visualization
          degree/degree_analysis_DK318.py  (calls plot_degree_brainmaps)
          edge/connect_ANOVA.py            (calls plot_connectome)
          → PNG brain surface maps and connectome plots
```

---

## Step-by-step Usage

### Step 1 — Compute MIND matrices

Edit the path variables at the top of `run_mind_all_subjects.py`:

```python
MIND_CODE_DIR = "/path/to/CODE/MIND"
BASE_FS_DIR   = "/path/to/freesurfer/subjects"   # root containing all subject folders
MIND_OUT_DIR  = "/path/to/output/MIND_out"
PARCELLATION  = "DK318"   # must match lh.<name>.annot / rh.<name>.annot in label/
```

Run directly or via the provided SLURM script:

```bash
# Interactive
python run_mind_all_subjects.py

# SLURM
sbatch run_mind_all_subjects.sh
```

Outputs per subject (in `MIND_OUT_DIR/`):

| File | Description |
|---|---|
| `<sub>_MIND_DK318.npy` | 318×318 MIND matrix (float32) |
| `<sub>_MIND_DK318.csv` | 318×318 matrix with ROI labels |
| `<sub>_MIND_DK318_edges.npy` | 50 403-dim upper-triangle vector |
| `<sub>_MIND_DK318_degree.npy` | 318-dim degree vector |
| `<sub>_MIND_DK318_degree.csv` | Degree with ROI names |

Subjects whose MIND matrix already exists are skipped automatically (checkpoint restart). Failed subjects are logged to `bad_subjects.log`.

---

### Step 2 — ComBat harmonization

Edit path variables in `combat/combat_dk318.py`:

```python
EXCEL_PATH = "/path/to/all_data_cqt.xlsx"   # covariate table (site, age_month, sex, group_d_or_c)
MIND_DIR   = "/path/to/MIND_out"
OUT_DIR    = "/path/to/MIND_out_combat_group"
```

Run:

```bash
python combat/combat_dk318.py
```

The script uses **neuroCombat** with `batch_col="site"` and covariates `age`, `sex`, `group`. Sites that contain only one diagnostic group are automatically excluded before harmonization.

Outputs per subject (in `OUT_DIR/`):

| File | Description |
|---|---|
| `<sub>_MIND_DK318_combat.npy` | 318×318 harmonized MIND matrix |

---

### Step 3 — Post-ComBat degree computation

Edit paths in `aftercombat_MIND_degree.py`, then run:

```bash
python aftercombat_MIND_degree.py
```

Outputs per subject:

| File | Description |
|---|---|
| `<sub>_degree.csv` | 318-dim degree with ROI names |
| `<sub>_MIND_DK318_combat_labeled.csv` | 318×318 matrix with ROI labels |

Degree is defined as the mean similarity to all other parcels:

```
degree_i = (sum_j≠i MIND_ij) / (N - 1)
```

---

### Step 4 — Degree analysis (ANOVA / ANCOVA)

**Full cohort — two-way ANOVA (Diagnosis × AgeGroup):**

```bash
sbatch degree/run_DK318.sh
# or: python degree/degree_analysis_DK318.py
```

**Adults only — ANCOVA with IQ covariate:**

```bash
sbatch degree/run_DK318_adult_with_IQ.sh
# or: python degree/degree_analysis_DK318_adult_with_IQ.py
```

**Site-effect verification:**

```bash
python degree/degree_analysis_DK318_site.py
```

All scripts output a CSV with per-ROI statistics including F-statistics, p-values, effect sizes (η²), and Benjamini–Hochberg FDR-corrected p-values. Significant results (p < 0.05 and FDR < 0.05) are saved as separate tables, and brain surface maps are rendered automatically.

---

### Step 5 — Edge-level analysis

```bash
python edge/connect_ANOVA.py
```

Runs a two-way ANOVA on every edge of the 318×318 MIND network (50 403 edges). Outputs:
- `ANOVA_DK318_edge_results.csv` — per-edge F/p/η² for Diagnosis, AgeGroup, and their interaction
- Connectome plots (`.png`) colored by signed −log₁₀(p) with effect direction (TD > DD = red, DD > TD = blue)

---

### Step 6 — Brain map visualization

Brain surface maps are generated automatically at the end of the degree and edge analysis scripts using **nilearn** `plot_surf_stat_map` and `plot_connectome`. Maps show:

- Signed −log₁₀(p) on the inflated fsaverage surface (lateral & medial views, both hemispheres)
- **coolwarm** colormap: red = group with higher degree / connectivity, blue = lower
- **viridis** colormap for unsigned effects (interaction, IQ)
- Parcel boundaries drawn as contour lines

Separate output directories are created per statistical column (e.g., `p_Diagnosis_FDR/`, `p_AgeGroup_FDR/`).

---

## Parcellations

| Name | ROIs | Annotation files | Notes |
|---|---|---|---|
| **DK-318** | 318 | `lh.DK318.annot`, `rh.DK318.annot` | High-resolution parcellation; primary analysis |
| **DK-68** | 68 | fsaverage `lh.aparc.annot`, `rh.aparc.annot` | Standard Desikan–Killiany atlas; replication |

Annotation files for DK-318 must be placed in each subject's FreeSurfer `label/` directory before running Step 1. The fsaverage-space `.annot` files should also be placed in `FILE/DK-318/` for visualization.

MNI coordinates for DK-318 parcels are extracted with `get_MNI.R` and stored in `FILE/DK-318/DK318_MNI_Coordinates.csv`.

---

## Dependencies

### Python

Install with:

```bash
pip install -r MIND/requirements.txt
pip install neuroCombat nilearn mne nibabel
```

| Package | Purpose |
|---|---|
| `numpy`, `pandas` | Array and data-frame operations |
| `scipy` | Statistical tests |
| `nibabel` | Reading FreeSurfer surface/annotation files |
| `nipype` | BBRegister / mri_vol2surf wrappers (Step 1 preprocessing) |
| `neuroCombat` | ComBat site harmonization |
| `nilearn` | Brain surface and connectome visualization |
| `mne` | Fetching fsaverage surface data |
| `matplotlib` | General plotting |
| `statsmodels` | FDR correction (Benjamini–Hochberg) |

### R

```r
install.packages(c("ggplot2", "dplyr", "lme4", "emmeans"))
```

Used for developmental trajectory plots (`DK318_Interaction_Trajectories.R`, `DK68_Interaction_Trajectories.R`) and MNI coordinate extraction (`get_MNI.R`).

### System

- **FreeSurfer** (≥ 6.0) — must be installed and `SUBJECTS_DIR` set correctly
- **FSL** — required by BBRegister for initialization
- **AFNI** — used in T1/T2 ratio computation (`3dcalc`)
- A SLURM-based HPC cluster is assumed for batch jobs (`.sh` scripts use `sbatch`); all Python scripts can also be run interactively

---

## Data Format

### Covariate table (`all_data_cqt.xlsx`)

Each row represents one subject. Required columns:

| Column | Description |
|---|---|
| `original-project` | Site / project identifier (string) |
| `id_old` | Subject ID within the project |
| `site` | Scanner site label used as ComBat batch variable |
| `age_month` | Age in months (continuous covariate) |
| `sex` | Biological sex (categorical, e.g. 0/1) |
| `group_d_or_c` | Diagnostic group: `0` = TD, `1` = DD |
| `group_age` | Age group: `1` = Child, `2` = Adult |

Subject directories in `BASE_FS_DIR` must be named `{original-project}_{id_old}` and contain a standard FreeSurfer `surf/` subdirectory.

### ROI name file (`DK318_roi_names.csv`)

A single-column CSV with header `region`, listing all 318 parcel names in the same order as the rows/columns of the MIND matrix.

---

## Notes

- **Checkpoint restart**: `run_mind_all_subjects.py` skips subjects whose output `.npy` already exists. Failed subjects are appended to `bad_subjects.log`.
- **ComBat site filtering**: Sites containing only TD or only DD subjects are automatically excluded from ComBat to ensure stable harmonization.
- **Parcellation naming**: The MIND code looks for annotation files named `lh.<PARCELLATION>.annot` and `rh.<PARCELLATION>.annot` inside each subject's FreeSurfer `label/` directory. The `PARCELLATION` variable must match exactly (e.g., `"DK318"` → `lh.DK318.annot`).
- **DK-68 vs DK-318**: All scripts have DK-68 counterparts (suffix `_dk68` or `_DK68`). The DK-68 pipeline is a replication analysis using the standard Desikan–Killiany parcellation bundled with fsaverage.
- **Effect direction in plots**: For Diagnosis comparisons, red = TD > DD; for AgeGroup, red = Adult > Child; for IQ correlations, red = positive correlation. All directions are computed from group mean matrices after ComBat.

---

---

# 阅读障碍 — MIND 网络分析流程（中文说明）

本仓库提供完整的分析流程，用于构建 **形态学逆距离（MIND）** 皮层相似性网络，进行多站点数据协调（ComBat），并对阅读障碍儿童及成人开展组间统计分析（ANOVA / ANCOVA）。

---

## 目录

1. [研究背景](#研究背景)
2. [仓库结构](#仓库结构)
3. [流程概览](#流程概览)
4. [逐步使用说明](#逐步使用说明)
   - [第一步 — 计算 MIND 矩阵](#第一步--计算-mind-矩阵)
   - [第二步 — ComBat 站点校正](#第二步--combat-站点校正)
   - [第三步 — 校正后 degree 计算](#第三步--校正后-degree-计算)
   - [第四步 — Degree 统计分析](#第四步--degree-统计分析)
   - [第五步 — 边水平统计分析](#第五步--边水平统计分析)
   - [第六步 — 脑图可视化](#第六步--脑图可视化)
5. [脑区分区方案](#脑区分区方案)
6. [依赖环境](#依赖环境)
7. [数据格式](#数据格式)
8. [注意事项](#注意事项)

---

## 研究背景

MIND 网络通过计算每对皮层脑区在顶点水平形态特征（皮层厚度 CT、髓鞘化指标 MC、体积 Vol、脑沟深度 SD、表面积 SA）上的相似性，来刻画大脑的结构组织方式。本仓库将 MIND 方法应用于多站点发育队列，探究典型发育（TD）读者与发育性阅读障碍（DD）个体在儿童期和成年期的皮层网络差异。

---

## 仓库结构

```
CODE/
├── MIND/                          # MIND 网络构建
│   ├── MIND.py                    # 核心函数 compute_MIND()
│   ├── get_vertex_df.py           # 顶点水平特征提取
│   ├── register_and_vol2surf.py   # 体积到表面的配准（BBRegister + mri_vol2surf）
│   └── requirements.txt           # Python 依赖包
│
├── combat/                        # 站点协调
│   ├── combat_dk318.py            # 对 DK-318 MIND 边向量做 ComBat
│   └── combat_dk68.py             # 对 DK-68 MIND 边向量做 ComBat
│
├── degree/                        # 节点水平（degree）分析
│   ├── degree_analysis_DK318.py                # 双因素 ANOVA：诊断 × 年龄组（DK-318）
│   ├── degree_analysis_DK318_site.py           # 站点效应检验（DK-318）
│   ├── degree_analysis_DK318_t_test.py         # 事后 t 检验（DK-318）
│   ├── degree_analysis_DK318_adult_with_IQ.py  # 含 IQ 协变量的 ANCOVA（成人，DK-318）
│   ├── degree_analysis_DK68.py                 # ANOVA（DK-68）
│   ├── degree_analysis_DK68_site.py            # 站点效应检验（DK-68）
│   ├── degree_MIND_heatmap.py                  # 组均值 MIND 矩阵热图
│   ├── DK318_brain_region_check.py             # 脑区覆盖率检查（DK-318）
│   ├── DK68_brain_region_check.py              # 脑区覆盖率检查（DK-68）
│   ├── DK318_Diagnosis_degree_barplot.py       # 按诊断组绘制 degree 柱状图（DK-318）
│   ├── DK68_Diagnosis_degree_barplot.py        # 按诊断组绘制 degree 柱状图（DK-68）
│   ├── DK318_Interaction_Trajectories.R        # 发育轨迹图（DK-318）
│   ├── DK68_Interaction_Trajectories.R         # 发育轨迹图（DK-68）
│   └── run_*.sh                                # SLURM 提交脚本
│
├── edge/                          # 边水平分析
│   ├── connect_ANOVA.py           # 对每条边做 ANOVA（DK-318）
│   ├── connect_ANOVA_site.py      # 边水平站点效应检验（DK-318）
│   ├── connect_ANOVA_dk68.py      # 对每条边做 ANOVA（DK-68）
│   ├── connect_ANOVA_dk68_site.py # 边水平站点效应检验（DK-68）
│   └── connect_t_test.py          # 事后边水平 t 检验
│
├── analysis/                      # 下游可视化与分析
│
├── aftercombat_MIND_degree.py     # ComBat 后：从校正矩阵重新计算 degree
├── run_mind_all_subjects.py       # 批量计算所有被试的 MIND 矩阵
├── run_mind_all_subjects.sh       # SLURM 批量提交脚本
├── get_MNI.R                      # 提取脑区 MNI 坐标
└── README.md
```

---

## 流程概览

```
FreeSurfer 输出
       │
       ▼
 [第一步] 逐被试计算 MIND 矩阵（DK-318 或 DK-68）
       │  run_mind_all_subjects.py
       │  → <被试>_MIND_DK318.npy        （318×318）
       │  → <被试>_MIND_DK318_edges.npy  （50403 维上三角向量）
       │  → <被试>_MIND_DK318_degree.npy （318 维 degree）
       │
       ▼
 [第二步] ComBat 站点校正（作用于边向量）
       │  combat/combat_dk318.py
       │  → <被试>_MIND_DK318_combat.npy  （318×318，已校正）
       │
       ▼
 [第三步] 从校正矩阵重新计算 degree
       │  aftercombat_MIND_degree.py
       │  → <被试>_degree.csv
       │  → <被试>_MIND_DK318_combat_labeled.csv
       │
       ▼
 [第四步] 节点水平统计（degree）
       │  degree/degree_analysis_DK318.py
       │  → ANOVA_DK318_degree_results.csv  （p 值、η²、FDR）
       │
       ▼
 [第五步] 边水平统计
       │  edge/connect_ANOVA.py
       │  → ANOVA_DK318_edge_results.csv
       │
       ▼
 [第六步] 脑图可视化
          degree/degree_analysis_DK318.py  （调用 plot_degree_brainmaps）
          edge/connect_ANOVA.py            （调用 plot_connectome）
          → PNG 脑表面图和连接组图
```

---

## 逐步使用说明

### 第一步 — 计算 MIND 矩阵

修改 `run_mind_all_subjects.py` 顶部的路径变量：

```python
MIND_CODE_DIR = "/path/to/CODE/MIND"
BASE_FS_DIR   = "/path/to/freesurfer/subjects"   # 包含所有被试文件夹的根目录
MIND_OUT_DIR  = "/path/to/output/MIND_out"
PARCELLATION  = "DK318"   # 需与 label/ 中 lh.<name>.annot / rh.<name>.annot 一致
```

直接运行或通过 SLURM 提交：

```bash
# 交互式运行
python run_mind_all_subjects.py

# SLURM 提交
sbatch run_mind_all_subjects.sh
```

每个被试的输出文件（保存在 `MIND_OUT_DIR/`）：

| 文件 | 说明 |
|---|---|
| `<被试>_MIND_DK318.npy` | 318×318 MIND 矩阵（float32）|
| `<被试>_MIND_DK318.csv` | 带 ROI 标签的 318×318 矩阵 |
| `<被试>_MIND_DK318_edges.npy` | 50403 维上三角边向量 |
| `<被试>_MIND_DK318_degree.npy` | 318 维 degree 向量 |
| `<被试>_MIND_DK318_degree.csv` | 带 ROI 名称的 degree 表 |

已存在输出文件的被试会自动跳过（断点续算）。计算失败的被试记录在 `bad_subjects.log`。

---

### 第二步 — ComBat 站点校正

修改 `combat/combat_dk318.py` 中的路径变量：

```python
EXCEL_PATH = "/path/to/all_data_cqt.xlsx"   # 协变量表（含 site、age_month、sex、group_d_or_c）
MIND_DIR   = "/path/to/MIND_out"
OUT_DIR    = "/path/to/MIND_out_combat_group"
```

运行：

```bash
python combat/combat_dk318.py
```

脚本使用 **neuroCombat**，以 `site` 为批次变量，协变量为 `age`、`sex`、`group`。仅含单一诊断组的站点会在校正前自动剔除。

每个被试的输出文件（保存在 `OUT_DIR/`）：

| 文件 | 说明 |
|---|---|
| `<被试>_MIND_DK318_combat.npy` | 318×318 站点校正后 MIND 矩阵 |

---

### 第三步 — 校正后 degree 计算

修改 `aftercombat_MIND_degree.py` 中的路径后运行：

```bash
python aftercombat_MIND_degree.py
```

每个被试的输出文件：

| 文件 | 说明 |
|---|---|
| `<被试>_degree.csv` | 带 ROI 名称的 318 维 degree |
| `<被试>_MIND_DK318_combat_labeled.csv` | 带 ROI 标签的 318×318 矩阵 |

Degree 定义为与所有其他脑区的平均相似度：

```
degree_i = (∑_{j≠i} MIND_ij) / (N - 1)
```

---

### 第四步 — Degree 统计分析

**全样本 — 双因素 ANOVA（诊断 × 年龄组）：**

```bash
sbatch degree/run_DK318.sh
# 或：python degree/degree_analysis_DK318.py
```

**仅成人 — 含 IQ 协变量的 ANCOVA：**

```bash
sbatch degree/run_DK318_adult_with_IQ.sh
# 或：python degree/degree_analysis_DK318_adult_with_IQ.py
```

**站点效应验证：**

```bash
python degree/degree_analysis_DK318_site.py
```

所有脚本均输出包含逐 ROI F 统计量、p 值、效应量（η²）及 Benjamini–Hochberg FDR 校正 p 值的 CSV 文件。显著结果（p < 0.05 及 FDR < 0.05）单独保存，并自动生成脑表面图。

---

### 第五步 — 边水平统计分析

```bash
python edge/connect_ANOVA.py
```

对 318×318 MIND 网络中的每条边（共 50 403 条）进行双因素 ANOVA。输出：
- `ANOVA_DK318_edge_results.csv` — 逐边的诊断、年龄组及其交互效应的 F/p/η²
- 连接组图（`.png`），以有符号 −log₁₀(p) 着色，效应方向：TD > DD = 红色，DD > TD = 蓝色

---

### 第六步 — 脑图可视化

脑表面图在 degree 和边分析脚本末尾自动生成，使用 **nilearn** 的 `plot_surf_stat_map` 和 
 