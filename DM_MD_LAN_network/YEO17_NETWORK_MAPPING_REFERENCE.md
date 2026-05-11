# Yeo 7 和 Yeo 17 网络映射关系

## 官方对应关系

根据 Yeo et al. (2011) J Neurophysiol 论文和官方文档

### Yeo 7 网络

| 编号 | 网络名称 |
|------|---------|
| 1 | Visual |
| 2 | Somatomotor |
| 3 | Dorsal Attention |
| 4 | Salience / Ventral Attention |
| 5 | Limbic |
| 6 | Control |
| 7 | Default |

### Yeo 17 网络

| 编号 | 网络名称 | 简称 | 对应的Yeo7网络 |
|------|---------|------|---------------|
| 1 | Visual Central, Visual A | VisCent | Visual (1) |
| 2 | Visual Peripheral, Visual B | VisPeri | Visual (1) |
| 3 | Somatomotor A | SomMotA | Somatomotor (2) |
| 4 | Somatomotor B | SomMotB | Somatomotor (2) |
| 5 | Dorsal Attention A | DorsAttnA | Dorsal Attention (3) |
| 6 | Dorsal Attention B | DorsAttnB | Dorsal Attention (3) |
| 7 | Salience / Ventral Attention A | SalVentAttnA | Salience/Ventral Attention (4) |
| 8 | Salience / Ventral Attention B | SalVentAttnB | Salience/Ventral Attention (4) |
| 9 | Limbic A | LimbicA | Limbic (5) |
| 10 | Limbic B | LimbicB | Limbic (5) |
| 11 | Control C | ContC | Control (6) |
| 12 | Control A | ContA | Control (6) |
| 13 | Control B | ContB | Control (6) |
| 14 | Temporal Parietal | TempPar | Temporal Parietal (独立) |
| 15 | Default C | DefaultC | Default (7) |
| 16 | Default A | DefaultA | Default (7) |
| 17 | Default B | DefaultB | Default (7) |

## 映射总结

### Yeo7 → Yeo17 的细分

1. **Visual (1)** → 2个子网络
   - 17Networks_1: VisCent (Visual Central)
   - 17Networks_2: VisPeri (Visual Peripheral)

2. **Somatomotor (2)** → 2个子网络
   - 17Networks_3: SomMotA
   - 17Networks_4: SomMotB

3. **Dorsal Attention (3)** → 2个子网络
   - 17Networks_5: DorsAttnA
   - 17Networks_6: DorsAttnB

4. **Salience/Ventral Attention (4)** → 2个子网络
   - 17Networks_7: SalVentAttnA
   - 17Networks_8: SalVentAttnB

5. **Limbic (5)** → 2个子网络
   - 17Networks_9: LimbicA
   - 17Networks_10: LimbicB

6. **Control (6)** → 3个子网络
   - 17Networks_11: ContC
   - 17Networks_12: ContA
   - 17Networks_13: ContB

7. **Default (7)** → 3个子网络
   - 17Networks_15: DefaultC
   - 17Networks_16: DefaultA
   - 17Networks_17: DefaultB

8. **Temporal Parietal** → 独立网络
   - 17Networks_14: TempPar

## 注意事项

1. **Temporal Parietal (17Networks_14)** 在Yeo17中是独立的网络，不属于Yeo7的任何一个网络的细分
2. **Control** 在某些文献中也称为 **Frontoparietal** 网络
3. **Salience/Ventral Attention** 在某些文献中简称为 **Ventral Attention**
4. FreeSurfer的annot文件中只包含编号（如17Networks_1），不包含语义名称
5. 语义名称（如VisCent）是根据官方文档添加的

## 参考文献

Yeo BT, Krienen FM, Sepulcre J, et al. The organization of the human cerebral cortex estimated by intrinsic functional connectivity. J Neurophysiol. 2011;106(3):1125-1165. DOI: 10.1152/jn.00338.2011

## 验证方法

可以通过以下脚本验证映射关系：
- `check_yeo7_labels.py` - 查看Yeo7网络的实际标签
- `check_yeo17_labels.py` - 查看Yeo17网络的实际标签
- `map_yeo7_to_dk318.py` - 将Yeo7映射到DK318
- `map_yeo17_to_dk318.py` - 将Yeo17映射到DK318

## 更新日期

2024年 - 根据官方对应表修正了映射关系
