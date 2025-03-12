# 加载必要的包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("recount3", quietly = TRUE))
    BiocManager::install("recount3")

# 设置代理
Sys.setenv(http_proxy = "http://127.0.0.1:7897")
Sys.setenv(https_proxy = "http://127.0.0.1:7897")

# 设置 SSL 验证选项
options(download.file.method = "libcurl")
# 如果上面的方法不行，可以尝试禁用 SSL 验证（不推荐，但可以临时使用）
options(download.file.method = "libcurl", ssl.verifypeer = FALSE)

library(recount3)
library(dplyr)

# 确保缓存目录存在
dir.create("D:/JiaoYuan/cache", recursive = TRUE, showWarnings = FALSE)

# 确保输出目录存在
dir.create("D:/JiaoYuan/TCGA_NORMAL", recursive = TRUE, showWarnings = FALSE)
dir.create("D:/JiaoYuan/TCGA_TUMOR", recursive = TRUE, showWarnings = FALSE)

# 设置 recount3 URL
url <- "https://recount-opendata.s3.amazonaws.com/recount3/release"
options(
    recount3_url = url,
    recount3_organism = "human"
)

# 验证 URL 设置
message("Current recount3 URL: ", getOption("recount3_url"))
message("Current organism: ", getOption("recount3_organism"))

# 设置缓存路径
options(recount3_cache = "D:/JiaoYuan/cache")
message("Cache directory: ", getOption("recount3_cache"))

# 获取可用的项目
print("Fetching available projects...")
human_projects <- available_projects()
print(human_projects)

# 筛选出 TCGA 的项目
print("Filtering TCGA projects...")
sub_proj <- subset(human_projects, file_source == "tcga")
print(sub_proj)

# 筛选出 HNSC 和 CHOL 项目
target_projects <- c("HNSC", "CHOL")
target_sub_proj <- sub_proj[sub_proj$project %in% target_projects, ]
print(target_sub_proj)

# 定义处理计数数据的函数，用于筛选正常组织和肿瘤组织样本
process_counts <- function(count_trans, sub_proj_rse) {
  # 转换 barcode，方便挑选样本
  barcode_to_file <- sub_proj_rse@colData@listData$tcga.tcga_barcode
  names(barcode_to_file) <- sub_proj_rse@colData@listData$tcga.gdc_file_id
  # 重命名 data 数据框的列名
  new_names <- lapply(colnames(count_trans), function(x) barcode_to_file[x])
  colnames(count_trans) <- unlist(new_names)
  
  # 创建空的正常组织样本列表和肿瘤组织样本列表
  normal_samples <- character(0)
  tumor_samples <- character(0)
  
  # 遍历列名
  # <10 是肿瘤组织，>=10 是正常组织
  for (col in colnames(count_trans)) {
    if (as.numeric(substring(col, 14, 15)) >= 10) {
      normal_samples <- c(normal_samples, col)
    } else {
      tumor_samples <- c(tumor_samples, col)
    }
  }
  
  # 根据正常组织样本列表筛选出正常组织数据
  count_trans_normal <- count_trans[, normal_samples]
  # 根据肿瘤组织样本列表筛选出肿瘤组织数据
  count_trans_tumor <- count_trans[, tumor_samples]
  
  return(list(normal = count_trans_normal, tumor = count_trans_tumor))
}

# 定义一个函数来处理行名（去除 _PAR_Y 和版本号）
process_rownames <- function(data) {
  # 找出行名中包含 "_PAR_Y" 的行索引
  rows_to_drop <- grep("_PAR_Y", rownames(data))
  # 删除这些行
  data <- data[-rows_to_drop, ]
  # 再去除版本号
  rownames(data) <- gsub("\\.\\d+$", "", rownames(data))
  return(data)
}

# 遍历目标项目
for (i in 1:nrow(target_sub_proj)) {
  print(paste("Processing project:", target_sub_proj[i, ]$project))
  
  # 创建 RangedSummarizedExperiment 对象
  print("Creating RSE object...")
  sub_proj_rse <- create_rse(target_sub_proj[i, ], type = c("gene"))
  print(sub_proj_rse)
  
  # 根据文献，read count 要缩放，用 auc 比较精确
  print("Transforming counts...")
  count_trans <- transform_counts(sub_proj_rse, by = "auc")
  print(head(count_trans))
  
  # 处理计数数据，筛选出正常组织和肿瘤组织样本
  print("Processing counts...")
  processed_counts <- process_counts(count_trans, sub_proj_rse)
  count_trans_normal <- processed_counts$normal
  count_trans_tumor <- processed_counts$tumor
  
  # 处理行名
  print("Processing rownames...")
  count_trans_normal <- process_rownames(count_trans_normal)
  count_trans_tumor <- process_rownames(count_trans_tumor)
  
  # 保存数据
  print("Saving data...")
  write.csv(count_trans_normal, file = paste0("D:/JiaoYuan/TCGA_NORMAL/", target_sub_proj[i, ]$project, "_normal.csv"))
  write.csv(count_trans_tumor, file = paste0("D:/JiaoYuan/TCGA_TUMOR/", target_sub_proj[i, ]$project, "_tumor.csv"))
}