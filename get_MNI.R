library(freesurferformats)
library(tidyverse)

# ==========================================
# 1. 设置路径
# ==========================================
# 你的 DK-318.annot 文件所在的文件夹路径
your_annot_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/DK-318"

# FreeSurfer 标准模板路径
fs_surf_dir <- "/data/software/freesurfer/subjects/fsaverage/surf/"

# ==========================================
# 2. 定义计算函数 (使用更稳健的直接读取法)
# ==========================================
get_centroids <- function(hemi, annot_file, surf_file) {
  
  message(paste("正在处理:", hemi, "半球..."))
  
  # 1. 读取大脑表面几何文件
  if(!file.exists(surf_file)) stop(paste("找不到文件:", surf_file))
  surf <- read.fs.surface(surf_file)
  coords <- as.data.frame(surf$vertices)
  colnames(coords) <- c("x", "y", "z")
  
  # 2. 读取分区注释文件
  if(!file.exists(annot_file)) stop(paste("找不到文件:", annot_file))
  annot <- read.fs.annot(annot_file)
  
  # =======================================================
  # 核心修正点：直接使用包提供的 label_names
  # =======================================================
  # 官方文档显示 read.fs.annot 返回的列表中直接包含 "label_names"
  # 这是一个已经匹配好的字符向量，长度与顶点数一致
  
  if("label_names" %in% names(annot)) {
    vertex_labels <- annot$label_names
  } else {
    # 备用方案：如果万一没有 label_names，尝试从 colortable 提取
    # 这里同时尝试 struct_names (复数) 和 struct_name (单数)
    message("警告: 未找到直接的 label_names，尝试手动匹配...")
    codes <- annot$label_codes
    ctab <- annot$colortable
    
    # 自动寻找正确的列名
    name_col <- grep("struct", names(ctab), value = TRUE)[1] 
    if(is.na(name_col)) stop("无法在颜色表中找到结构名称列！")
    
    vertex_labels <- ctab[[name_col]][match(codes, ctab$code)]
  }

  # 3. 合并并计算质心
  df <- coords %>%
    mutate(label = vertex_labels) %>%
    # 移除无效的标签 (NA 或 unknown)
    filter(!is.na(label)) %>%
    filter(label != "unknown") %>% 
    filter(label != "???") %>%
    group_by(label) %>%
    summarise(
      x = mean(x),
      y = mean(y),
      z = mean(z),
      n_vertices = n()
    ) %>%
    mutate(hemi = hemi)
    
  return(df)
}

# ==========================================
# 3. 执行提取
# ==========================================

# 左脑
lh_data <- get_centroids(
  hemi = "left",
  annot_file = file.path(your_annot_dir, "lh.DK318.annot"), 
  surf_file  = file.path(fs_surf_dir, "lh.pial")
)

# 右脑
rh_data <- get_centroids(
  hemi = "right",
  annot_file = file.path(your_annot_dir, "rh.DK318.annot"), 
  surf_file  = file.path(fs_surf_dir, "rh.pial")
)

# ==========================================
# 4. 合并并保存
# ==========================================
final_coords <- bind_rows(lh_data, rh_data)

print("------------------------------------------------")
print(head(final_coords))
print(paste("共提取脑区数量:", nrow(final_coords)))
print("------------------------------------------------")

if(nrow(final_coords) > 0) {
  save_path <- "/data/home/tqi/data1/share/after_freesurfer/FILE/DK-318/DK318_MNI_Coordinates.csv"
  write.csv(final_coords, save_path, row.names = FALSE)
  message(paste("成功！文件已保存至:", save_path))
} else {
  message("错误：未能提取到任何脑区，请检查 annot 文件是否为空或格式异常。")
}