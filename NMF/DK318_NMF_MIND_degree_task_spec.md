# DK318 MIND Degree NMF Analysis Task Specification

## 1. 分析目标

请新写一个独立的 R 脚本，用于完成 DK318 MIND degree 的 NMF 分析。这个脚本不要继续加长前面的 DGLM 脚本。

建议脚本名称为：

```text
01_nmf_mind_degree.R
```

本分析的核心目标是：

```text
基于所有被试的 DK318 MIND degree 矩阵，进行非负矩阵分解，得到数据驱动的 MIND degree 脑区模式。
随后用每个被试在这些模式上的 loading，检验这些模式是否与 AgeGroup、Diagnosis 以及 AgeGroup × Diagnosis 有关。
```

这个分析和前一个 Tmap 聚类分析不同。

前一个 Tmap 分析是：

```text
先对每个 ROI 做 DGLM，得到每个 ROI 的 Child 和 Adult TD-DD t 值，然后对 ROI 聚类。
```

本 NMF 分析是：

```text
先对 subject × ROI 的原始 MIND degree 矩阵做 NMF，得到 K 个数据驱动的脑区模式，再对每个模式的 subject loading 做统计检验。
```

---

## 2. 输入矩阵

输入矩阵为：

```text
X = n × 318
```

其中：

```text
n = 被试数量
318 = DK318 ROI 数量
X[i, j] = 第 i 个被试在第 j 个 ROI 上的 MIND degree
```

矩阵的行是被试，列是 ROI。

只需要读取每个被试的 degree 文件，不需要读取完整 MIND matrix。

---

## 3. NMF 分解目标

对输入矩阵 X 做非负矩阵分解：

```text
X ≈ W × H
```

其中：

```text
W = n × K
H = K × 318
```

含义如下：

```text
W 矩阵：每个被试在每个 NMF component 上的 loading，也可以叫 score。
H 矩阵：每个 NMF component 在 318 个 ROI 上的权重分布。
K：NMF 分解出的 latent components 数量。
```

后续统计检验使用 W 矩阵。

每一列 `score_k` 表示一个被试在第 k 个 NMF component 上的表达强度。

H 矩阵用于解释每个 component 对应哪些 ROI 或网络。

---

## 4. 关于 K 的解释

老师笔记中写到 K clusters，但在 NMF 中更准确的说法是：

```text
K latent components
```

NMF 默认得到的是 K 个 component，不是严格的 hard clustering。

如果后续需要把 ROI 分成 K 类，可以根据 H 矩阵把每个 ROI 分配给权重最大的 component：

```text
ROI_j belongs to argmax_k H[k, j]
```

这样可以得到 K 个 ROI clusters。

但本脚本中主要输出 NMF components 和 subject loadings。

---

## 5. 输入文件和路径

沿用之前脚本的路径。

### 5.1 demo 文件

```r
demo_file <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/all_data_cqt_mean_1.5.xlsx"
```

### 5.2 degree 文件所在目录

```r
mind_combat_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_combat"
```

每个被试的 degree 文件命名方式沿用旧脚本：

```r
subj_prefix = paste0(original_project, "_", id_old)
file_base = paste0(subj_prefix, "_MIND_DK318_combat")
degree_file = file.path(mind_combat_dir, paste0(file_base, "_degree.csv"))
```

### 5.3 输出目录

建议输出目录为：

```r
out_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NMF"
```

如果目录不存在，自动创建。

---

## 6. 被试信息和因子编码

从 demo 文件中读取被试信息，并沿用之前的编码方式。

```r
Diagnosis = ifelse(group_d_or_c == 0, "TD", "DD")
AgeGroup  = ifelse(group_age == 1, "Adult", "Child")
Sex       = ifelse(sex == 1, "Male", "Female")
```

因子水平必须固定为：

```r
Diagnosis = factor(Diagnosis, levels = c("TD", "DD"))
AgeGroup  = factor(AgeGroup,  levels = c("Child", "Adult"))
Sex       = factor(Sex,       levels = c("Female", "Male"))
```

只保留存在 degree 文件的被试。

如果某个被试缺失 Diagnosis、AgeGroup 或 Sex，则在统计模型中排除。

---

## 7. 是否分年龄组单独做 NMF

不要分别在 Child 和 Adult 中单独做 NMF。

主分析应当使用所有被试共同构建一个统一的 NMF basis。

原因是：

```text
如果 Child 和 Adult 分开做 NMF，两个年龄组得到的 component 不一定一一对应。
这样后面无法合理地做 score_k ~ AgeGroup * Diagnosis + Sex。
```

因此，本脚本应当：

```text
把 Child、Adult、TD、DD 的所有有效被试放在一起。
构建统一的 n × 318 degree 矩阵。
在全体被试上做 NMF。
得到共同的 H 矩阵。
再用 W 矩阵中的 score_k 做 AgeGroup × Diagnosis 检验。
```

---

## 8. NMF 输入数据的非负处理

NMF 要求输入矩阵非负。

请在做 NMF 前检查：

```r
min(X, na.rm = TRUE)
```

如果 X 中已经没有负值，可以直接使用 X。

如果 X 中有负值，需要做非负化处理。

推荐处理方式：对每个 ROI 做 min-max scaling 到 `[0, 1]`。

对第 j 个 ROI：

```text
X_scaled[, j] = (X[, j] - min(X[, j])) / (max(X[, j]) - min(X[, j]))
```

如果某个 ROI 的方差接近 0，记录并处理为 0 或排除该 ROI。由于 DK318 需要保留 318 个 ROI，优先保留该 ROI，并将该列 scaling 后设为 0。

请保存 scaling 信息，便于复现：

```text
ROI
min_value
max_value
range
was_constant
```

输出文件建议为：

```text
NMF_degree_scaling_info.csv
```

---

## 9. K 的选择

不要只固定一个 K。

需要尝试：

```r
K = 2:10
```

每个 K 需要进行多次 NMF 初始化，建议：

```r
nrun = 100
```

如果运行太慢，可以先用：

```r
nrun = 50
```

K 选择至少输出以下指标：

```text
reconstruction error
RSS 或 residual sum of squares
cophenetic coefficient，如果所用 R 包支持
silhouette score，如果可以基于 W 或 ROI assignment 计算
```

最终选择 K 时，需要综合：

```text
reconstruction error 是否继续明显下降
cophenetic coefficient 或 stability 是否稳定
component 是否具有可解释性
每个 component 是否有足够 ROI 权重或被试 loading 变化
```

如果没有明确数据驱动的最佳 K，可以保留一个假设驱动 K，比如 K=4，并同时输出 K=2:10 的评估结果。

输出文件建议为：

```text
NMF_K_selection_metrics.csv
NMF_K_selection_reconstruction_error.png
NMF_K_selection_cophenetic.png
```

---

## 10. 推荐使用的 R 包

推荐使用：

```r
NMF
```

如果没有安装，请自动安装或提示安装。

需要的 R 包包括：

```r
readxl
dplyr
stringr
NMF
broom
car
ggplot2
```

如果使用 `car::Anova` 做 type III ANOVA，请设置 contrasts：

```r
options(contrasts = c("contr.sum", "contr.poly"))
```

也可以使用线性模型 `lm()` 后提取系数和 ANOVA 结果。

---

## 11. 最终 NMF 模型

确定主分析 K 后，拟合最终 NMF：

```r
fit <- nmf(X_scaled, rank = K_final, method = "brunet", nrun = 100, seed = 2026)
```

需要提取：

```r
W <- basis(fit)
H <- coef(fit)
```

注意不同包里 W 和 H 的函数名称可能不同。请检查维度，确保：

```text
W = n × K
H = K × 318
```

如果 R 包返回的是相反方向，请转置到这个格式。

---

## 12. 对 subject loadings 做统计检验

对 W 矩阵中的每个 component score 做统计模型。

模型为：

```r
score_k ~ AgeGroup * Diagnosis + Sex
```

其中：

```text
score_k 是 W 矩阵中第 k 个 component 的 subject loading。
AgeGroup 和 Diagnosis 有交互。
Sex 只作为协变量。
```

不要加入：

```text
age_within_group
IQ
Diagnosis × Sex
AgeGroup × Sex
Diagnosis × AgeGroup × Sex
```

每个 component 都要输出以下 term 的统计结果：

```text
Diagnosis
AgeGroup
AgeGroup:Diagnosis
Sex
```

最重要的 term 是：

```text
AgeGroup:Diagnosis
```

它回答：

```text
该 NMF component 的 TD-DD 差异是否在 Child 和 Adult 中不同。
```

---

## 13. 统计输出要求

统计结果输出文件建议为：

```text
NMF_component_lm_anova_results.csv
```

每一行是一个 component 的一个 term。

建议列包括：

```text
component
term
df
statistic
p_value
q_value
```

对多个 component 和多个 term 做 FDR 校正。

建议分别对每类 term 做 FDR：

```text
Diagnosis across K components
AgeGroup across K components
AgeGroup:Diagnosis across K components
Sex across K components
```

也可以额外输出一个 across all tests 的 q 值。

列名建议为：

```text
q_value_by_term
q_value_all_terms
```

---

## 14. Subject loading 输出

输出 W 矩阵，并合并被试信息。

文件名建议为：

```text
NMF_W_subject_loadings.csv
```

每行一个被试，列包括：

```text
subject_id
original_project
id_old
file_base
Diagnosis
AgeGroup
Sex
score_1
score_2
...
score_K
```

如果后续要画箱线图，这个文件是主要输入。

---

## 15. Component ROI weight 输出

输出 H 矩阵。

文件名建议为：

```text
NMF_H_component_roi_weights.csv
```

格式建议为 wide format：

```text
component
ROI_1
ROI_2
...
ROI_318
```

也可以额外输出 long format：

```text
NMF_H_component_roi_weights_long.csv
```

long format 列包括：

```text
component
ROI
weight
rank_within_component
```

`rank_within_component` 按每个 component 内的 weight 从高到低排序。

---

## 16. Top ROI 输出

对每个 component 输出权重最高的 ROI。

建议输出 top 10、top 20 和 top 10% ROI。

文件名建议为：

```text
NMF_component_top_ROIs.csv
```

列包括：

```text
component
ROI
weight
rank_within_component
top_type
```

其中 `top_type` 可以是：

```text
top10
top20
top10percent
```

---

## 17. ROI hard assignment

为了把 NMF components 解释成 K 个 ROI clusters，可以根据 H 矩阵给每个 ROI 分配一个主 component。

规则：

```text
ROI_j assigned_component = argmax_k H[k, j]
```

输出文件建议为：

```text
NMF_ROI_hard_assignment.csv
```

列包括：

```text
ROI
assigned_component
max_weight
second_max_weight
assignment_margin
```

其中：

```text
assignment_margin = max_weight - second_max_weight
```

assignment_margin 越大，说明该 ROI 更明确地属于某个 component。

---

## 18. 可视化输出

请至少输出以下图：

### 18.1 K 选择图

```text
NMF_K_selection_reconstruction_error.png
NMF_K_selection_cophenetic.png
```

如果无法计算 cophenetic，则至少输出 reconstruction error。

### 18.2 Subject loading 组间箱线图

对每个 component 的 score_k 画图。

x 轴：

```text
AgeGroup
```

颜色或分组：

```text
Diagnosis
```

建议文件名：

```text
NMF_component_scores_boxplot.png
```

如果 component 很多，可以每个 component 单独保存，或用 facet_wrap。

### 18.3 Component ROI weight 图

每个 component 画 318 个 ROI 的权重柱状图或排序图。

建议文件名：

```text
NMF_component_roi_weights.png
```

---

## 19. 显著 component 的解释

对于统计检验中显著的 component，尤其是：

```text
AgeGroup:Diagnosis q_value_by_term < 0.05
```

需要进一步查看：

```text
该 component 的 H 权重最高的 ROI 是哪些
这些 ROI 是否集中在 DMN、LAN、MD 或其他预定义网络中
不同组别的 score 分布如何
```

如果没有 FDR 显著的 component，也需要输出 raw p 值最小的 component，并作为探索性结果报告。

---

## 20. 网络富集分析预留接口

本脚本可以先不实现网络富集，但要为后续网络富集保存足够文件。

后续网络富集会使用：

```text
NMF_component_top_ROIs.csv
NMF_ROI_hard_assignment.csv
```

富集分析的目标是：

```text
检查每个 component 的 top ROI 或 hard-assigned ROI 是否富集在 DMN、LAN、MD 等网络中。
```

可用方法包括：

```text
Fisher exact test
hypergeometric test
permutation test
```

---

## 21. 质量控制要求

需要输出基本 QC 信息。

文件名建议为：

```text
NMF_degree_QC_summary.csv
```

内容包括：

```text
n_subjects_total_in_demo
n_subjects_with_degree_file
n_subjects_used
n_roi
min_raw_degree
max_raw_degree
min_scaled_degree
max_scaled_degree
n_missing_values
n_negative_raw_values
K_range_tested
K_final
NMF_method
nrun
seed
```

如果 degree matrix 中有 NA，需要处理。

建议处理方式：

```text
如果某个被试或 ROI 有 NA，优先报错并输出位置。
不要静默填补，除非明确指定填补策略。
```

---

## 22. 主要输出文件清单

最终至少应输出这些文件：

```text
NMF_degree_QC_summary.csv
NMF_degree_scaling_info.csv
NMF_K_selection_metrics.csv
NMF_K_selection_reconstruction_error.png
NMF_K_selection_cophenetic.png
NMF_W_subject_loadings.csv
NMF_H_component_roi_weights.csv
NMF_H_component_roi_weights_long.csv
NMF_component_lm_anova_results.csv
NMF_component_top_ROIs.csv
NMF_ROI_hard_assignment.csv
NMF_component_scores_boxplot.png
NMF_component_roi_weights.png
```

如果某些图由于包或数据问题不能生成，需要在日志中说明原因，但不要影响核心 CSV 输出。

---

## 23. 推荐脚本结构

R 脚本建议分成以下模块：

```text
1. Load packages
2. Define paths and options
3. Load demo file and recode variables
4. Read degree files and build X matrix
5. Check non-negativity and scale X if needed
6. Evaluate K from 2 to 10
7. Fit final NMF model with K_final
8. Extract W and H
9. Run lm or ANOVA for each component score
10. Apply FDR correction
11. Save W, H, top ROI, ROI assignment
12. Generate plots
13. Save QC summary
```

---

## 24. 后续如何解释结果

统计模型：

```r
score_k ~ AgeGroup * Diagnosis + Sex
```

解释方式：

```text
Diagnosis 显著：该 component 的 subject loading 在 TD 和 DD 之间不同。
AgeGroup 显著：该 component 的 subject loading 在 Child 和 Adult 之间不同。
AgeGroup × Diagnosis 显著：TD-DD 差异随年龄组变化，提示该 MIND degree pattern 具有发育阶段特异性。
Sex 显著：该 component 的 loading 与性别有关。
```

最重要的结果是：

```text
AgeGroup × Diagnosis
```

如果某个 component 的 `AgeGroup × Diagnosis` 显著，就查看该 component 的 H 权重最高 ROI，并进一步做网络富集。

---

## 25. 最终分析逻辑总结

本任务需要完成：

```text
第一步，构建所有被试的 DK318 MIND degree 矩阵 X，维度为 n × 318。
第二步，对 X 做非负处理，满足 NMF 输入要求。
第三步，在 K = 2:10 范围内评估 NMF 分解质量。
第四步，选择 K_final，拟合最终 NMF，得到 W 和 H。
第五步，对 W 中每个 component 的 subject loading 做 score_k ~ AgeGroup * Diagnosis + Sex。
第六步，输出每个 component 的统计检验结果。
第七步，解释显著 component 的 H 权重分布和 top ROI。
```

这个分析的核心问题是：

```text
数据驱动得到的 MIND degree 脑区模式，是否存在 TD-DD 差异，是否存在年龄组差异，以及 TD-DD 差异是否随儿童和成人阶段而改变。
```
