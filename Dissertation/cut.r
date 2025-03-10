# 加载必要的库
library(raster)
library(rgdal)

# 设定文件夹路径和shapefile路径
tif_folder <- "/home/yuanj/work/dissertation"
shp_file <- "/home/yuanj/work/dissertation/huaxi3.0.shp"

# 读取掩膜Shapefile
mask_shapefile <- readOGR(shp_file)

# 获取文件夹内所有TIF文件
tif_files <- list.files(tif_folder, pattern = "\\.tif$", full.names = TRUE)

# 对每个TIF文件进行掩膜操作
for (tif_file in tif_files) {
  
  # 读取TIF文件
  raster_layer <- raster(tif_file)
  
  # 使用Shapefile进行掩膜
  masked_raster <- mask(raster_layer, mask_shapefile)
  
  # 将掩膜后的图像覆盖原始文件
  writeRaster(masked_raster, tif_file, format = "GTiff", overwrite = TRUE)
}

print("掩膜操作完成！")
