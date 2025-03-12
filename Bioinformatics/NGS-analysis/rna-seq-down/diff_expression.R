library(edgeR)
library(ggplot2)
library(pheatmap)

run_deg = function(normal_file, tumor_file, output_name) {
  # 读取数据
  normal = read.csv(normal_file, row.names = 1)
  tumor = read.csv(tumor_file, row.names = 1)
  
  # 合并数据
  merged_data <- merge(normal, tumor, by = 'row.names')
  rownames(merged_data) <- merged_data$Row.names
  merged_data <- merged_data[, -1]
  
  # 设置分组信息
  normal_g = rep(1, dim(normal)[2])
  tumor_g = rep(2, dim(tumor)[2])
  group = factor(c(normal_g, tumor_g))
  
  # 创建 DGEList 对象
  y <- DGEList(counts = merged_data, group = group)
  
  # 过滤低表达基因
  keep_exprs <- filterByExpr(y, group = group)
  y <- y[keep_exprs, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y)
  
  # 差异分析
  design <- model.matrix(~group)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  
  # 获取结果
  diff = topTags(qlf, n = nrow(y), sort.by = 'logFC')
  diff = as.data.frame(diff)
  result = merge(diff, merged_data, by = 'row.names')
  
  # 整理结果
  rownames(result) <- result$Row.names
  result <- result[, -1]
  result = result[row.names(diff),]
  
  # 绘制火山图
  volcano_plot <- ggplot(result, aes(x = logFC, y = -log10(FDR))) +
    geom_point(aes(color = abs(logFC) > 1 & FDR < 0.05)) +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    labs(title = paste(output_name, "Volcano Plot"),
         x = "log2 Fold Change",
         y = "-log10 FDR")
  
  # 保存火山图
  ggsave(paste0("D:/JiaoYuan/DEG/", output_name, "_volcano.png"), 
         volcano_plot, 
         width = 8, 
         height = 6)
  
  # 选择显著差异的前50个基因进行热图绘制
  top_genes <- head(rownames(result[order(result$FDR), ]), 50)
  heatmap_data <- merged_data[top_genes, ]
  
  # 数据标准化
  heatmap_data <- t(scale(t(heatmap_data)))
  
  # 准备注释信息
  annotation_col <- data.frame(
    Group = c(rep("Normal", dim(normal)[2]), rep("Tumor", dim(tumor)[2]))
  )
  rownames(annotation_col) <- colnames(heatmap_data)
  
  # 绘制热图
  png(paste0("D:/JiaoYuan/DEG/", output_name, "_heatmap.png"), 
      width = 2000, 
      height = 3000, 
      res = 300)
  pheatmap(heatmap_data,
           annotation_col = annotation_col,
           show_colnames = FALSE,
           main = paste(output_name, "Top 50 DEGs"),
           fontsize = 8)
  dev.off()
  
  # 保存结果
  write.csv(result, paste0("D:/JiaoYuan/DEG/", output_name, "_deg.csv"))
}

# 定义要分析的项目
projects <- c("CHOL", "HNSC")

# 遍历每个项目进行差异分析
for (project in projects) {
  print(paste("Processing", project))
  
  normal_file <- paste0("D:/JiaoYuan/TCGA_NORMAL/", project, "_normal.csv")
  tumor_file <- paste0("D:/JiaoYuan/TCGA_TUMOR/", project, "_tumor.csv")
  
  # 确保输出目录存在
  dir.create("D:/JiaoYuan/DEG", recursive = TRUE, showWarnings = FALSE)
  
  # 运行差异分析
  run_deg(normal_file, tumor_file, project)
}