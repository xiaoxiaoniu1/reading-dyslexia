# DK318 三网络聚合分析任务说明

## 0. 任务概述

请新写一套独立脚本，用于完成 DK318 MIND degree 的三网络聚合分析。

本分析不再使用 17 个 Yeo 网络作为主分析，而是直接使用已有的三个理论网络：

```text
DM / Default network
MD / Multiple Demand network
Reading network
```

三个网络允许有重叠。不要把 ROI 强行互斥分配，不要按优先级去掉重叠 ROI。每个网络都按照它自己的原始 ROI 列表独立计算 network-level degree。

需要完成两类统计模型：

```r
Degree ~ Network * Diagnosis * AgeGroup + Sex + (1 | SubjectID)
```

以及每个网络单独模型：

```r
network_degree ~ Diagnosis * AgeGroup + Sex
```

注意，分析脚本和可视化脚本必须分开。

建议脚本命名：

```text
01_DK318_3network_degree_analysis.R
02_DK318_3network_degree_visualization.R
```

---

## 1. 分析目标

本分析的目标是从预定义网络层面研究 TD、DD 在儿童和成人中的 MIND degree 差异模式。

原始输入是：

```text
X = n × 318
```

其中：

```text
n = 被试数量
318 = DK318 ROI 数量
X[i, j] = 第 i 个被试在第 j 个 DK318 ROI 上的 MIND degree
```

通过三个网络标签，将 `n × 318` 聚合为：

```text
X_net = n × 3
```

三列分别是：

```text
DM_degree
MD_degree
Reading_degree
```

每个网络 degree 推荐用该网络内所有 ROI degree 的均值计算：

```text
network_degree_i = mean(subject i 在该网络所有 ROI 上的 degree)
```

不要用 sum，因为不同网络包含的 ROI 数不同。

---

## 2. 和前两个分析的区别

这个分析是预定义网络层面的聚合统计分析。

它不同于前面的 Tmap 聚类：

```text
Tmap 聚类：先对每个 ROI 做 TD-DD 统计检验，再基于 Child 和 Adult 的 t 值对 318 个 ROI 聚类。
```

它也不同于 NMF：

```text
NMF：先从 n × 318 原始 degree 矩阵中数据驱动分解出 K 个 component，再检验 subject loading 与 Diagnosis、AgeGroup 的关系。
```

当前三网络聚合分析是：

```text
先按已有的 DM、MD、Reading 网络标签聚合 ROI，再直接在网络层面检验 Diagnosis、AgeGroup 及其交互。
```

这个分析回答的问题是：

```text
TD-DD 的发育差异是否体现在 DM、MD、Reading 这些理论网络层面？
不同网络是否具有不同的 Diagnosis × AgeGroup 差异模式？
```

---

## 3. 输入文件

### 3.1 被试信息文件

沿用之前 DGLM 脚本中的 demo 文件：

```r
demo_file <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/all_data_cqt_mean_1.5.xlsx"
```

读取：

```r
sheet = "Sheet1"
```

### 3.2 DK318 degree 文件目录

沿用之前 DGLM 脚本中的 degree 文件目录：

```r
mind_combat_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_combat"
```

每个被试的 degree 文件命名方式沿用旧脚本：

```r
subj_prefix = paste0(original_project, "_", id_old)
file_base = paste0(subj_prefix, "_MIND_DK318_combat")
degree_file = file.path(mind_combat_dir, paste0(file_base, "_degree.csv"))
```

每个 degree 文件需要包含：

```text
ROI
degree
```

### 3.3 三个网络标签文件

请在脚本开头设置以下路径变量，方便后续修改：

```r
reading_file <- "mean_gt1p0_union_reading_v4_v5_regions_used.csv"
md_file      <- "multiple_demand_regions_used.csv"
dm_file      <- "default_regions_used.csv"
```
这三个网络的路径为
reading:/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/brain_region_table_brainmaps/mean_gt1p0_union_reading_v4_v5/mean_gt1p0_union_reading_v4_v5_regions_used.csv
md:/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/brain_region_table_brainmaps/multiple_demand/multiple_demand_regions_used.csv
dm:/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/brain_region_table_brainmaps/default/default_regions_used.csv
每个网络文件中至少应包含：

```text
feature
```

`feature` 列与 degree 文件中的 `ROI` 名称对应。

三个网络名称统一为：

```text
Reading
MD
DM
```

### 3.4 输出目录

建议输出目录：

```r
out_root <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_3network_analysis"
```

分析表格输出到：

```r
analysis_dir <- file.path(out_root, "analysis")
```

图像输出到：

```r
fig_dir <- file.path(out_root, "figures")
```

如果目录不存在，自动创建。

---

## 4. 因子编码要求

沿用之前脚本的分组规则。

```r
Diagnosis = ifelse(group_d_or_c == 0, "TD", "DD")
AgeGroup  = ifelse(group_age == 1, "Adult", "Child")
Sex       = ifelse(sex == 1, "Male", "Female")
```

因子水平必须设置为：

```r
Diagnosis = factor(Diagnosis, levels = c("TD", "DD"))
AgeGroup  = factor(AgeGroup,  levels = c("Child", "Adult"))
Sex       = factor(Sex,       levels = c("Female", "Male"))
Network   = factor(Network,   levels = c("DM", "MD", "Reading"))
```

TD-DD contrast 的方向必须是：

```text
TD - DD
```

所以解释方向是：

```text
estimate > 0 或 t > 0 表示 TD > DD
estimate < 0 或 t < 0 表示 TD < DD，也就是 DD > TD
```

---

## 5. 网络重叠处理

三个网络允许有重叠，不要做互斥处理。

当前已知情况大致是：

```text
Reading 与 MD 有较多重叠
Reading 与 DM 有少量重叠
MD 与 DM 无明显重叠
三者共同重叠为 0
```

脚本需要重新从输入文件中自动计算并输出完整重叠结果，不要硬编码这些数量。

需要输出：

```text
network_membership_3net.csv
network_overlap_summary_3net.csv
```

### 5.1 network_membership_3net.csv

每个唯一 ROI 一行，至少包含：

```text
feature
in_DM
in_MD
in_Reading
n_networks
membership_pattern
```

其中：

```text
in_DM, in_MD, in_Reading 为 TRUE/FALSE
n_networks 表示该 ROI 属于几个网络
membership_pattern 例如 DM_only, MD_only, Reading_only, Reading+MD, Reading+DM, MD+DM, DM+MD+Reading
```

### 5.2 network_overlap_summary_3net.csv

输出每个网络大小和两两重叠：

```text
comparison
n_overlap
n_network_1
n_network_2
percent_of_network_1
percent_of_network_2
```

例如：

```text
Reading_vs_MD
Reading_vs_DM
MD_vs_DM
All_three
```

---

## 6. 数据构建要求

### 6.1 读取所有被试 degree

只读取 degree 文件，不需要读取完整 MIND matrix。

构造：

```text
degree_mat = n × 318
```

行是被试，列是 ROI。

同时保留被试信息：

```text
SubjectID
subj_prefix
file_base
Diagnosis
AgeGroup
Sex
degree_file
```

建议：

```r
SubjectID = subj_prefix
```

或者：

```r
SubjectID = file_base
```

但整个脚本中必须保持一致。

### 6.2 计算三网络 degree

对每个被试和每个网络，计算网络内 ROI degree 的均值：

```text
DM_degree_i = mean(degree_mat[i, DM_ROIs], na.rm = TRUE)
MD_degree_i = mean(degree_mat[i, MD_ROIs], na.rm = TRUE)
Reading_degree_i = mean(degree_mat[i, Reading_ROIs], na.rm = TRUE)
```

因为允许重叠，所以同一个 ROI 可以同时参与 MD_degree 和 Reading_degree 的计算。

### 6.3 ROI 标准化选项

请在脚本开头设置一个选项：

```r
standardize_roi_before_network_mean <- FALSE
```

默认使用原始 degree 求网络平均。

如果该选项为 TRUE，则先对每个 ROI 在所有被试中 z-score，再计算网络均值：

```r
degree_z = scale(degree_mat)
```

然后使用 z 值求网络均值。

无论是否标准化，都需要在输出 summary 中记录该设置。

---

## 7. 输出聚合矩阵

分析脚本需要输出两个数据文件。

### 7.1 宽表

文件名：

```text
DK318_3network_degree_wide.csv
```

格式：

```text
SubjectID
subj_prefix
file_base
Diagnosis
AgeGroup
Sex
DM_degree
MD_degree
Reading_degree
```

### 7.2 长表

文件名：

```text
DK318_3network_degree_long.csv
```

格式：

```text
SubjectID
subj_prefix
file_base
Diagnosis
AgeGroup
Sex
Network
Degree
```

每个被试有三行，分别对应：

```text
DM
MD
Reading
```

---

## 8. 模型一：三网络混合效应模型

### 8.1 模型公式

使用长表数据拟合：

```r
Degree ~ Network * Diagnosis * AgeGroup + Sex + (1 | SubjectID)
```

建议使用：

```r
lme4::lmer()
lmerTest::anova(type = 3)
emmeans
```

为了 Type III ANOVA 稳定，建议在模型前设置：

```r
options(contrasts = c("contr.sum", "contr.poly"))
```

但 emmeans 里的 TD-DD contrast 必须按 factor levels 显式指定。

### 8.2 需要检验的固定效应

输出 Type III ANOVA 结果，至少包含：

```text
Network
Diagnosis
AgeGroup
Sex
Network:Diagnosis
Network:AgeGroup
Diagnosis:AgeGroup
Network:Diagnosis:AgeGroup
```

最关键的是：

```text
Network:Diagnosis:AgeGroup
```

它回答：

```text
TD-DD 的儿童/成人差异是否因网络不同而不同？
```

输出文件：

```text
mixed_model_3network_type3_anova.csv
mixed_model_3network_summary.txt
```

### 8.3 mixed model 的 adjusted means

用 `emmeans` 输出四组调整均值：

```r
emmeans(fit_mixed, ~ Diagnosis * AgeGroup | Network)
```

输出文件：

```text
mixed_model_3network_emmeans_adjusted_means.csv
```

每行至少包含：

```text
Network
Diagnosis
AgeGroup
emmean
SE
df
lower.CL
upper.CL
```

### 8.4 mixed model 的 TD-DD follow-up contrast

对每个 Network 和 AgeGroup，提取：

```text
Child: TD - DD
Adult: TD - DD
```

建议：

```r
emm_diag <- emmeans(fit_mixed, ~ Diagnosis | Network * AgeGroup)
contrast(emm_diag, method = list(TD_minus_DD = c(1, -1)), by = c("Network", "AgeGroup"))
```

输出文件：

```text
mixed_model_3network_TD_minus_DD_by_network_age.csv
```

每行至少包含：

```text
Network
AgeGroup
contrast
estimate
SE
df
t.ratio
p.value
q_by_age
q_all
```

FDR 处理要求：

```text
q_by_age: 分别在 Child 和 Adult 内，对 3 个网络做 FDR
q_all: 对所有 6 个 Network × AgeGroup contrast 一起做 FDR
```

### 8.5 mixed model 的发育差异 contrast

需要对每个网络提取：

```text
Adult(TD-DD) - Child(TD-DD)
```

这个 contrast 不要用两个 t 值相减。必须从模型的 emmeans cell means 中构造 contrast。

权重逻辑：

```text
Adult TD = +1
Adult DD = -1
Child TD = -1
Child DD = +1
```

输出文件：

```text
mixed_model_3network_developmental_contrast.csv
```

每行至少包含：

```text
Network
contrast
estimate
SE
df
t.ratio
p.value
q_value
```

这个结果回答：

```text
某个网络的 TD-DD 差异是否随 AgeGroup 改变。
```

---

## 9. 模型二：每个网络单独分析

除了 mixed model，还需要每个网络单独做线性模型。

### 9.1 模型公式

对宽表中的每个网络分别拟合：

```r
DM_degree ~ Diagnosis * AgeGroup + Sex
MD_degree ~ Diagnosis * AgeGroup + Sex
Reading_degree ~ Diagnosis * AgeGroup + Sex
```

也可以使用长表筛选每个网络后拟合：

```r
Degree ~ Diagnosis * AgeGroup + Sex
```

建议使用：

```r
lm()
car::Anova(type = 3)
emmeans
```

### 9.2 输出 per-network ANOVA

输出文件：

```text
per_network_lm_type3_anova.csv
```

每行至少包含：

```text
Network
term
statistic
df
p.value
q_value
```

term 包括：

```text
Diagnosis
AgeGroup
Sex
Diagnosis:AgeGroup
```

FDR 需要至少对三个网络的 `Diagnosis:AgeGroup` p 值做校正。也可以对所有网络和所有 term 一起输出一个 `q_all`。

### 9.3 per-network adjusted means

输出：

```text
per_network_emmeans_adjusted_means.csv
```

每行至少包含：

```text
Network
Diagnosis
AgeGroup
emmean
SE
df
lower.CL
upper.CL
```

### 9.4 per-network TD-DD follow-up

每个网络内提取：

```text
Child: TD - DD
Adult: TD - DD
Adult(TD-DD) - Child(TD-DD)
```

输出文件：

```text
per_network_TD_minus_DD_by_age.csv
per_network_developmental_contrast.csv
```

列名和 mixed model 对应，至少包含：

```text
Network
AgeGroup
estimate
SE
df
t.ratio
p.value
q_by_age
q_all
```

发育差异 contrast 文件至少包含：

```text
Network
estimate
SE
df
t.ratio
p.value
q_value
```

---

## 10. 质量控制和日志

分析脚本需要输出以下 QC 文件。

### 10.1 被试数量和分组数量

文件：

```text
QC_subject_cell_counts.csv
```

至少包含：

```text
Diagnosis
AgeGroup
Sex
n
```

并在 console 中打印：

```text
总被试数
有 degree 文件的被试数
进入最终分析的被试数
Diagnosis × AgeGroup × Sex 表
```

### 10.2 网络 ROI 检查

文件：

```text
QC_network_roi_coverage.csv
```

每个网络一行：

```text
Network
n_roi_in_network_file
n_roi_found_in_degree_matrix
n_roi_missing_from_degree_matrix
missing_roi_list
```

如果某个网络在 degree 矩阵中找不到任何 ROI，脚本应报错停止。

如果只有少量 ROI 缺失，允许继续，但要记录 warning。

### 10.3 模型状态

文件：

```text
QC_model_status.csv
```

包含 mixed model 和每个 per-network model 是否成功：

```text
model_name
status
message
```

如果某个模型失败，不要让整个脚本静默失败。需要记录错误，并尽可能继续输出其他结果。

---

## 11. 可视化脚本要求

第二个脚本只负责可视化，不重新拟合模型。

脚本名称：

```text
02_DK318_3network_degree_visualization.R
```

它需要读取第一个脚本输出的 CSV 文件，并在 `figures` 目录下生成图像。

建议每张图同时保存：

```text
.png
.pdf
```

如果只保存一种，优先保存 `.png`，分辨率 300 dpi。

---

## 12. 可视化图一：网络重叠图

### 12.1 membership pattern bar plot

读取：

```text
network_membership_3net.csv
```

画：

```text
x = membership_pattern
y = ROI count
```

输出：

```text
fig_01_network_overlap_membership_bar.png
fig_01_network_overlap_membership_bar.pdf
```

目的：

```text
说明三个网络允许重叠，并展示 Reading、MD、DM 的实际重叠情况。
```

---

## 13. 可视化图二：三网络四组 adjusted mean 图

读取：

```text
mixed_model_3network_emmeans_adjusted_means.csv
```

画 adjusted mean profile。

建议格式：

```text
x = AgeGroup
y = emmean
color = Diagnosis
facet = Network
error bar = ± SE 或 95% CI
```

输出：

```text
fig_02_mixed_model_adjusted_means_by_network.png
fig_02_mixed_model_adjusted_means_by_network.pdf
```

这张图回答：

```text
DM、MD、Reading 三个网络中，TD 和 DD 在儿童和成人中的调整后 network degree 如何变化。
```

---

## 14. 可视化图三：原始数据分布图

读取：

```text
DK318_3network_degree_long.csv
```

画 raw network degree 的四组分布。

建议格式：

```text
x = AgeGroup
y = Degree
color = Diagnosis
facet = Network
使用 boxplot 或 violin plot
叠加 jitter points
```

输出：

```text
fig_03_raw_network_degree_distribution.png
fig_03_raw_network_degree_distribution.pdf
```

目的：

```text
展示原始数据分布，避免只看模型均值。
```

---

## 15. 可视化图四：TD-DD t 值 heatmap

读取：

```text
mixed_model_3network_TD_minus_DD_by_network_age.csv
mixed_model_3network_developmental_contrast.csv
```

构造矩阵：

```text
rows = DM, MD, Reading
columns = Child TD-DD, Adult TD-DD, Adult(TD-DD)-Child(TD-DD)
value = t.ratio
```

画 heatmap。

要求：

```text
颜色表示 t.ratio
t > 0 表示 TD > DD
t < 0 表示 DD > TD
每个格子标注 t 值
如果 q < 0.05，在格子中加 *
如果 q < 0.01，加 **
如果 q < 0.001，加 ***
```

输出：

```text
fig_04_mixed_model_TD_minus_DD_t_heatmap.png
fig_04_mixed_model_TD_minus_DD_t_heatmap.pdf
```

这张图是该分析的核心结果图之一。

---

## 16. 可视化图五：TD-DD contrast profile

读取：

```text
mixed_model_3network_TD_minus_DD_by_network_age.csv
```

画每个网络在 Child 和 Adult 中的 TD-DD contrast。

建议格式：

```text
x = AgeGroup
y = estimate
group/color = Network
error bar = ± SE
加 y = 0 参考线
```

或者：

```text
x = Network
y = estimate
color = AgeGroup
error bar = ± SE
```

输出：

```text
fig_05_mixed_model_TD_minus_DD_contrast_profile.png
fig_05_mixed_model_TD_minus_DD_contrast_profile.pdf
```

目的：

```text
直接展示每个网络中 TD-DD 差异在 Child 和 Adult 之间是否变化。
```

---

## 17. 可视化图六：per-network 模型结果对照图

读取：

```text
per_network_TD_minus_DD_by_age.csv
per_network_developmental_contrast.csv
```

生成与 mixed model 类似的 heatmap 或 contrast profile。

输出：

```text
fig_06_per_network_TD_minus_DD_t_heatmap.png
fig_06_per_network_TD_minus_DD_t_heatmap.pdf
```

目的：

```text
展示 per-network 独立模型的结果，作为 mixed model 的补充。
```

---

## 18. 推荐最终汇报逻辑

这个分析最后可以这样报告：

```text
We aggregated DK318 MIND degree values into three theory-driven network summaries, including DM, MD, and Reading networks. Because the networks were not mutually exclusive, each network summary was computed independently using its original ROI definition. We then tested whether network-level degree showed diagnosis-by-age-group effects using both a mixed-effects model across networks and separate per-network linear models.
```

中文表述：

```text
我们将 DK318 MIND degree 聚合到 DM、MD 和 Reading 三个理论网络中。由于三个网络不是完全互斥的，每个网络按照其原始 ROI 定义独立计算 network-level degree。随后，我们使用跨网络混合效应模型和每个网络单独线性模型，检验网络层面的 MIND degree 是否存在 Diagnosis × AgeGroup 差异。
```

---

## 19. 结果解释重点

### 19.1 mixed model

重点看：

```text
Network:Diagnosis:AgeGroup
```

如果显著，说明：

```text
TD-DD 的发育差异模式在 DM、MD、Reading 三个网络之间不同。
```

然后看每个网络的 follow-up：

```text
Child TD-DD
Adult TD-DD
Adult(TD-DD)-Child(TD-DD)
```

### 19.2 per-network model

重点看每个网络的：

```text
Diagnosis:AgeGroup
```

如果某个网络显著，说明：

```text
该网络内 TD-DD 差异随 AgeGroup 改变。
```

### 19.3 方向解释

所有 follow-up contrast 方向都是：

```text
TD - DD
```

所以：

```text
t > 0 表示 TD > DD
t < 0 表示 DD > TD
```

### 19.4 重叠解释

由于允许网络重叠，结果解释时不要说三个网络是完全独立或互斥的 ROI 集合。

应表述为：

```text
DM、MD 和 Reading 是三个理论定义的 network summaries。部分 ROI 可同时参与多个 network summary，尤其 Reading 与 MD 存在部分重叠。
```

---

## 20. 最终需要生成的文件清单

### 20.1 分析脚本

```text
01_DK318_3network_degree_analysis.R
```

### 20.2 可视化脚本

```text
02_DK318_3network_degree_visualization.R
```

### 20.3 分析输出

```text
DK318_3network_degree_wide.csv
DK318_3network_degree_long.csv
network_membership_3net.csv
network_overlap_summary_3net.csv
QC_subject_cell_counts.csv
QC_network_roi_coverage.csv
QC_model_status.csv

mixed_model_3network_type3_anova.csv
mixed_model_3network_summary.txt
mixed_model_3network_emmeans_adjusted_means.csv
mixed_model_3network_TD_minus_DD_by_network_age.csv
mixed_model_3network_developmental_contrast.csv

per_network_lm_type3_anova.csv
per_network_emmeans_adjusted_means.csv
per_network_TD_minus_DD_by_age.csv
per_network_developmental_contrast.csv
```

### 20.4 图像输出

```text
fig_01_network_overlap_membership_bar.png
fig_02_mixed_model_adjusted_means_by_network.png
fig_03_raw_network_degree_distribution.png
fig_04_mixed_model_TD_minus_DD_t_heatmap.png
fig_05_mixed_model_TD_minus_DD_contrast_profile.png
fig_06_per_network_TD_minus_DD_t_heatmap.png
```

如有条件，每张图同时保存 `.pdf` 版本。

---

## 21. 不要做的事情

本任务中不要做以下内容：

```text
不要使用 17 个 Yeo 网络作为主分析
不要把三个网络强行互斥化
不要删除 Reading 和 MD 的重叠 ROI
不要重新做 Tmap 聚类
不要重新做 NMF
不要做 edge 分析
不要加入 age_within_group
不要加入 IQ
不要加入 Sex interaction
不要把 Sex 写进 Diagnosis 或 AgeGroup 的交互项
```

本任务只做：

```text
三网络 degree 聚合
mixed model
per-network model
对应的可视化
```
