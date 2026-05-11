# DK318 MIND Degree DGLM Tmap 分析任务说明

## 1. 当前任务目标

请新写一个独立的 R 脚本，不要继续加长之前的全量 DGLM 脚本。

脚本名称建议为：

```text
01_dglm_degree_tmap_TD_minus_DD.R
```

这个脚本只做 DK318 MIND degree 的 Tmap 提取，为后续 ROI 聚类准备输入矩阵。

当前分析目标是：基于 DK318 MIND degree，计算每个 ROI 在儿童组和成人组内的 TD-DD 差异 t 值。

最终输出一个：

```text
318 × 2
```

的 Tmap 矩阵。

每一行是一个 DK318 ROI。

两列分别是：

```text
t_child_TD_minus_DD
t_adult_TD_minus_DD
```

方向必须是：

```text
TD - DD
```

解释方向为：

```text
t > 0 表示 TD > DD
t < 0 表示 TD < DD，也就是 DD > TD
```

这个 318 × 2 的矩阵后面会用于 ROI 聚类，研究 TD-DD 差异在儿童和成人中的发育模式。

---

## 2. 模型要求

对每个 ROI 的 degree 单独拟合一个 DGLM。

mean model 使用：

```r
y ~ Diagnosis * AgeGroup + Sex
```

dispersion model 使用同样的公式：

```r
~ Diagnosis * AgeGroup + Sex
```

注意事项：

```text
Sex 只作为协变量使用
不要加入 Diagnosis × Sex
不要加入 AgeGroup × Sex
不要加入 Diagnosis × AgeGroup × Sex
不要使用 age_within_group
不要使用 IQ
```

所以最终模型只包含：

```text
Diagnosis 主效应
AgeGroup 主效应
Diagnosis × AgeGroup 交互
Sex 协变量
```

当前主分析只做 degree mean Tmap。dispersion model 可以保留在 DGLM 中，但主输出和聚类只使用 mean model 的 TD-DD t 值。

---

## 3. 因子编码要求

Diagnosis 的水平必须是：

```r
Diagnosis = factor(Diagnosis, levels = c("TD", "DD"))
```

AgeGroup 的水平必须是：

```r
AgeGroup = factor(AgeGroup, levels = c("Child", "Adult"))
```

Sex 的水平必须是：

```r
Sex = factor(Sex, levels = c("Female", "Male"))
```

原始数据中的映射规则沿用旧脚本：

```r
Diagnosis = ifelse(group_d_or_c == 0, "TD", "DD")
AgeGroup  = ifelse(group_age == 1, "Adult", "Child")
Sex       = ifelse(sex == 1, "Male", "Female")
```

---

## 4. 输入文件和路径

沿用之前脚本的路径。

demo 文件：

```r
demo_file <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/all_data_cqt_mean_1.5.xlsx"
```

degree 文件所在目录：

```r
mind_combat_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_combat"
```

每个被试的 degree 文件命名方式沿用旧脚本：

```r
subj_prefix = paste0(original_project, "_", id_old)
file_base = paste0(subj_prefix, "_MIND_DK318_combat")
degree_file = file.path(mind_combat_dir, paste0(file_base, "_degree.csv"))
```

只需要读取 degree 文件，不需要读取完整 MIND matrix。

输出目录建议为：

```r
out_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_Tmap"
```

如果目录不存在，需要自动创建。

---

## 5. 每个 ROI 需要提取的统计量

每个 ROI 拟合 DGLM 后，用 `emmeans` 从 mean model 中提取：

```r
emmeans(model, ~ Diagnosis | AgeGroup, weights = "equal")
```

然后在每个 AgeGroup 内做 TD-DD contrast。

由于 Diagnosis 水平顺序是：

```text
TD, DD
```

所以 contrast 权重必须是：

```r
TD_minus_DD = c(1, -1)
```

也就是：

```r
contrast(
  emm,
  method = list(TD_minus_DD = c(1, -1)),
  by = "AgeGroup"
)
```

需要分别提取 Child 和 Adult 的：

```text
estimate
SE
t value
p value
```

输出列名建议为：

```text
estimate_child_TD_minus_DD
se_child_TD_minus_DD
t_child_TD_minus_DD
p_child_TD_minus_DD

estimate_adult_TD_minus_DD
se_adult_TD_minus_DD
t_adult_TD_minus_DD
p_adult_TD_minus_DD
```

然后对 318 个 ROI 分别做 FDR 校正：

```text
q_child_TD_minus_DD
q_adult_TD_minus_DD
```

---

## 6. 还需要输出一个发育差异检验

除了 Child 和 Adult 各自的 TD-DD t 值，还需要输出一个正式的发育阶段差异检验：

```text
Adult(TD-DD) - Child(TD-DD)
```

这个结果不是聚类输入，但用于解释聚类结果时验证某个 ROI 的 TD-DD 差异是否真的随年龄组改变。

输出列名建议为：

```text
estimate_dev_adult_minus_child
se_dev_adult_minus_child
t_dev_adult_minus_child
p_dev_adult_minus_child
q_dev_adult_minus_child
```

这个 contrast 可以通过 `emmeans(model, ~ Diagnosis * AgeGroup, weights = "equal")` 构造。

权重逻辑是：

```text
Adult TD = +1
Adult DD = -1
Child TD = -1
Child DD = +1
```

也就是：

```text
Adult(TD-DD) - Child(TD-DD)
```

不要直接从两个 t 值相减得到这个统计量，要从模型 contrast 中计算。

---

## 7. 质量控制要求

每个 ROI 需要输出一个 `fit_status`。

至少包含这些状态：

```text
ok
low_n
near_zero_variance
empty_diagnosis_age_cell
dglm_error
```

建模前需要检查：

```text
y 是否有限
Diagnosis 是否非缺失
AgeGroup 是否非缺失
Sex 是否非缺失
每个 Diagnosis × AgeGroup cell 是否至少有样本
Sex 是否至少包含 Female 和 Male 两个水平
degree 是否接近零方差
```

因为 Sex 只是协变量，所以不要求每个 Diagnosis × AgeGroup × Sex 小格子都有样本。

但至少四个 Diagnosis × AgeGroup 组合都要有样本：

```text
TD Child
DD Child
TD Adult
DD Adult
```

如果某个 ROI 模型失败，不要中断整个脚本。该 ROI 对应的统计量填 NA，并记录错误信息到 `fit_status`。

---

## 8. 需要保存的输出文件

### 8.1 完整统计结果

文件名：

```text
DGLM_DK318_degree_Tmap_TD_minus_DD_full.csv
```

这个文件每个 ROI 一行，包含：

```text
feature
n
fit_status

raw_n_child_TD
raw_mean_child_TD
raw_sd_child_TD

raw_n_child_DD
raw_mean_child_DD
raw_sd_child_DD

raw_n_adult_TD
raw_mean_adult_TD
raw_sd_adult_TD

raw_n_adult_DD
raw_mean_adult_DD
raw_sd_adult_DD

estimate_child_TD_minus_DD
se_child_TD_minus_DD
t_child_TD_minus_DD
p_child_TD_minus_DD
q_child_TD_minus_DD

estimate_adult_TD_minus_DD
se_adult_TD_minus_DD
t_adult_TD_minus_DD
p_adult_TD_minus_DD
q_adult_TD_minus_DD

estimate_dev_adult_minus_child
se_dev_adult_minus_child
t_dev_adult_minus_child
p_dev_adult_minus_child
q_dev_adult_minus_child
```

### 8.2 聚类输入矩阵

文件名：

```text
DGLM_DK318_degree_Tmap_TD_minus_DD_318x2.csv
```

这个文件只保留：

```text
feature
t_child_TD_minus_DD
t_adult_TD_minus_DD
```

这是后续聚类脚本的核心输入。

### 8.3 简单 summary

文件名：

```text
DGLM_DK318_degree_Tmap_TD_minus_DD_summary.csv
```

里面统计：

```text
contrast
n_tested
n_sig_FDR_0.05
min_p
min_q
```

contrast 包括：

```text
Child_TD_minus_DD
Adult_TD_minus_DD
Adult_minus_Child_change
```

---

## 9. 后续聚类脚本接口

请再预留一个后续脚本，名称建议为：

```text
02_cluster_degree_tmap.R
```

这个脚本读取：

```text
DGLM_DK318_degree_Tmap_TD_minus_DD_318x2.csv
```

然后用：

```r
X <- as.matrix(tmap[, c("t_child_TD_minus_DD", "t_adult_TD_minus_DD")])
rownames(X) <- tmap$feature
```

作为聚类输入。

聚类对象是 ROI，不是被试。

每个 ROI 的二维特征是：

```text
儿童组 TD-DD 的 t 值
成人组 TD-DD 的 t 值
```

后续聚类目标是识别不同的发育差异模式，比如：

```text
儿童和成人中都 TD > DD
儿童和成人中都 DD > TD
儿童 TD > DD 但成人 DD > TD
儿童 DD > TD 但成人 TD > DD
儿童差异强而成人差异弱
成人差异强而儿童差异弱
```

---

## 10. 聚类分析要求

主分析可以先使用 K=4，因为这是一个假设驱动的方案。但同时要做 K 选择。

K 选择范围：

```r
K = 2:8
```

需要输出：

```text
elbow curve
silhouette score
cluster size
cluster centroid
```

K=4 的主结果需要输出：

```text
degree_Tmap_cluster_K4.csv
degree_Tmap_cluster_K4_centroids.csv
degree_Tmap_cluster_K4_scatter.png
degree_Tmap_k_selection.csv
degree_Tmap_elbow.png
degree_Tmap_silhouette.png
```

scatter 图要求：

```text
x 轴是 t_child_TD_minus_DD
y 轴是 t_adult_TD_minus_DD
每个点是一个 ROI
颜色表示 cluster
画 x = 0 参考线
画 y = 0 参考线
画 y = x 参考线
```

解释 cluster 时，需要输出每个 cluster 的：

```text
n
mean_t_child
mean_t_adult

median_t_child
median_t_adult

mean_overall_TD_minus_DD
mean_developmental_change
```

其中：

```text
overall_TD_minus_DD = (t_child_TD_minus_DD + t_adult_TD_minus_DD) / 2

developmental_change = t_adult_TD_minus_DD - t_child_TD_minus_DD
```

注意：`developmental_change` 只是描述性指标。正式的发育差异统计仍然看前一个 DGLM 脚本输出的：

```text
t_dev_adult_minus_child
p_dev_adult_minus_child
q_dev_adult_minus_child
```

---

## 11. 最终分析逻辑

这个项目现在要完成的是：

```text
第一步，基于 DK318 MIND degree，对每个 ROI 拟合 DGLM。
第二步，在 Child 和 Adult 内分别提取 TD-DD 的模型调整后 t 值。
第三步，得到 318 × 2 Tmap 矩阵。
第四步，基于这个 318 × 2 矩阵对 ROI 聚类。
第五步，用聚类结果解释 TD-DD 差异在儿童和成人中的发育模式。
```

模型必须是：

```r
y ~ Diagnosis * AgeGroup + Sex
```

不要使用：

```text
age_within_group
IQ
Sex interaction
edge analysis
dispersion Tmap 聚类
```

当前主分析只做 degree mean Tmap。
