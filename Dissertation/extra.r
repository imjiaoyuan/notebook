# 加载必要的库
library(raster)
library(sp)
library(readr)  # 用于读取CSV文件
library(dplyr)  # 用于数据处理

# 设定文件夹路径和CSV文件路径
tif_folder <- "/home/yuanj/work/dissertation/模拟用数据"
csv_file <- "/home/yuanj/work/dissertation/BD.csv"

# 读取细菌.csv文件，假设Lon列和Lat列存在
bacteria_data <- read.csv(csv_file)

# 获取文件夹内所有TIF文件
tif_files <- list.files(tif_folder, pattern = "\\.tif$", full.names = TRUE)

# 对每个TIF文件进行处理，并提取值
for (tif_file in tif_files) {
  
  # 读取TIF文件
  raster_layer <- raster(tif_file)
  
  # 提取经纬度数据的栅格值
  # 创建一个SpatialPoints对象，用来存放Lon和Lat
  coords <- data.frame(lon = bacteria_data$Lon, lat = bacteria_data$Lat)
  coordinates(coords) <- ~lon+lat
  proj4string(coords) <- CRS(proj4string(raster_layer))  # 确保坐标系统一致
  
  # 提取这些点对应的栅格值
  extracted_values <- extract(raster_layer, coords)
  
  # 将提取的值添加到原始数据中
  bacteria_data[[basename(tif_file)]] <- extracted_values
}

# 将所有TIF文件提取的数据保存为一个新的CSV文件
write.csv(bacteria_data, "/home/yuanj/work/dissertation/提取结果.csv", row.names = FALSE, fileEncoding = "UTF-8")

print("数据提取完成，已保存为提取结果.csv")