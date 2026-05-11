# Yeo 7网络标签对应关系说明

## 来源和依据

这个映射关系基于以下权威来源：

### 1. 原始论文
**Yeo BT, Krienen FM, Sepulcre J, et al. (2011)**  
*The organization of the human cerebral cortex estimated by intrinsic functional connectivity.*  
Journal of Neurophysiology, 106(3):1125-1165.  
DOI: 10.1152/jn.00338.2011

### 2. FreeSurfer官方文档
FreeSurfer Wiki: CorticalParcellation_Yeo2011  
https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011

### 3. 7个网络的标准定义

根据Yeo et al. 2011论文，7个大尺度脑网络按照以下顺序定义：

| 索引 | FreeSurfer标签 | 标准网络名称 | 中文名称 | 主要功能 |
|------|---------------|-------------|---------|---------|
| 1 | 7Networks_1 | Visual | 视觉网络 | 视觉信息处理 |
| 2 | 7Networks_2 | Somatomotor | 躯体运动网络 | 感觉运动功能 |
| 3 | 7Networks_3 | Dorsal Attention | 背侧注意网络 | 自上而下的注意控制 |
| 4 | 7Networks_4 | Ventral Attention | 腹侧注意网络 | 刺激驱动的注意 |
| 5 | 7Networks_5 | Limbic | 边缘系统网络 | 情绪和记忆 |
| 6 | 7Networks_6 | Frontoparietal | 额顶网络 | 执行控制和认知灵活性 |
| 7 | 7Networks_7 | Default | 默认模式网络 | 自我参照和内部思维 |

### 4. 网络的解剖位置

- **Visual (视觉网络)**: 主要位于枕叶，包括初级和次级视觉皮层
- **Somatomotor (躯体运动网络)**: 中央前回和中央后回，包括运动和感觉皮层
- **Dorsal Attention (背侧注意网络)**: 额叶眼动区和顶内沟
- **Ventral Attention (腹侧注意网络)**: 腹侧额叶和颞顶联合区
- **Limbic (边缘系统网络)**: 眶额皮层、颞极和海马旁回
- **Frontoparietal (额顶网络)**: 背外侧前额叶和后顶叶皮层
- **Default (默认模式网络)**: 内侧前额叶、后扣带回和角回

### 5. 验证方法

可以通过以下方式验证这个映射关系：

1. **查看FreeSurfer源代码和文档**
2. **使用nilearn库** (Python神经影像分析库)
   ```python
   from nilearn import datasets
   yeo = datasets.fetch_atlas_yeo_2011()
   # 查看网络标签
   ```
3. **参考Thomas Yeo实验室的GitHub仓库**
   https://github.com/ThomasYeoLab/CBIG

### 6. 注意事项

⚠️ **重要提醒**：
- 这个顺序是FreeSurfer中Yeo2011_7Networks_N1000.annot文件的标准顺序
- 不同的软件包或工具可能使用不同的编号或颜色方案
- 建议在使用前通过可视化验证映射是否正确
- 可以通过FreeView或其他可视化工具加载annot文件来确认

### 7. 如何验证

运行以下命令来检查实际的标签：
```bash
cd /data/home/tqi/data1/share/after_freesurfer/CODE/DM_MD_LAN_network
python3 check_yeo7_labels.py
```

或者在FreeView中可视化：
```bash
freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.inflated:annot=lh.Yeo2011_7Networks_N1000.annot
```

### 8. 参考资源

- **Yeo Lab官网**: https://sites.google.com/view/yeolab
- **论文PDF**: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3174820/
- **Nilearn文档**: https://nilearn.github.io/modules/description/yeo_2011.html

---

**结论**: 这个映射关系是基于Yeo et al. 2011年发表的权威论文和FreeSurfer官方实现，不是随意编造的。但建议在实际使用前通过可视化或读取annot文件来验证具体的索引对应关系。
