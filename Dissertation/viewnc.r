library(ncdf4)
library(terra)

# 打开 NetCDF 文件
nc_file <- nc_open("/home/yuanj/work/dissertation/pre_2023.nc")

# 打开一个文本文件用于输出
sink("nc_file_info.txt")

# 输出文件的基本信息
print(nc_file)

# 输出文件的维度信息
dimensions <- nc_file$dim
print(dimensions)

# 输出文件的变量信息
variables <- nc_file$var
print(variables)

# 输出全局属性（如文件描述、作者等）
global_attrs <- ncatt_get(nc_file, 0)
print(global_attrs)

# 关闭文本文件输出
sink()

# 关闭 NetCDF 文件
nc_close(nc_file)