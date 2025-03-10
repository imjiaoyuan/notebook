library(raster)
library(caret)
library(rgdal)
library(magrittr)
library(xgboost)
library(Matrix)

aspect <- raster("/home/yuanj/work/dissertation/模拟用数据/aspect.tif")
PET <- raster("/home/yuanj/work/dissertation/模拟用数据/PET.tif")
S1 <- raster("/home/yuanj/work/dissertation/模拟用数据/S1.tif")
Dc <- raster("/home/yuanj/work/dissertation/模拟用数据/Dc.tif")
NDVI <- raster("/home/yuanj/work/dissertation/模拟用数据/NDVI.tif")
GRAV <- raster("/home/yuanj/work/dissertation/模拟用数据/GRAV.tif")
CW1 <- raster("/home/yuanj/work/dissertation/模拟用数据/CW1.tif")
BD <- raster("/home/yuanj/work/dissertation/模拟用数据/BD.tif")
MG <- raster("/home/yuanj/work/dissertation/模拟用数据/MG.tif")
slope <- raster("/home/yuanj/work/dissertation/模拟用数据/slope.tif")

# h <- raster("/home/yuanj/work/dissertation/模拟用数据/Dh.tif")
# H <- raster("/home/yuanj/work/dissertation/模拟用数据/H.tif")
# AL <- raster("/home/yuanj/work/dissertation/模拟用数据/AL.tif")
# CL <- raster("/home/yuanj/work/dissertation/模拟用数据/CL.tif")
# LNUM <- raster("/home/yuanj/work/dissertation/模拟用数据/LNUM.tif")
# LDEP <- raster("/home/yuanj/work/dissertation/模拟用数据/LDEP.tif")
# AN <- raster("/home/yuanj/work/dissertation/模拟用数据/AN.tif")
# CA <- raster("/home/yuanj/work/dissertation/模拟用数据/CA.tif")
# DEM <- raster("/home/yuanj/work/dissertation/模拟用数据/DEM.tif")
# PH <- raster("/home/yuanj/work/dissertation/模拟用数据/PH.tif")
# K <- raster("/home/yuanj/work/dissertation/模拟用数据/K.tif")
# N <- raster("/home/yuanj/work/dissertation/模拟用数据/NA.tif")
# PDEP <- raster("/home/yuanj/work/dissertation/模拟用数据/PDEP.tif")
# Dc <- raster("/home/yuanj/work/dissertation/模拟用数据/Dc.tif")
# AK <- raster("/home/yuanj/work/dissertation/模拟用数据/AK.tif")
# MG <- raster("/home/yuanj/work/dissertation/模拟用数据/MG.tif")
# CEC <- raster("/home/yuanj/work/dissertation/模拟用数据/CEC.tif")
# CW1 <- raster("/home/yuanj/work/dissertation/模拟用数据/CW1.tif")
# AP <- raster("/home/yuanj/work/dissertation/模拟用数据/AP.tif")
# PRE <- raster("/home/yuanj/work/dissertation/模拟用数据/PRE.tif")
# TN <- raster("/home/yuanj/work/dissertation/模拟用数据/TN.tif")
# TMP <- raster("/home/yuanj/work/dissertation/模拟用数据/TMP.tif")

#xgbDART

data <- read.csv("/home/yuanj/work/dissertation/RESULTS2.csv")

nrow(data)

names(data)

### remove outliers by mean +- 2*sd
mean(data$BD);min(data$BD);max(data$BD)
mean <- mean(data$BD)
sd <- sd(data$BD)

up.outlier <- mean +2*sd
low.outlier <- mean - 2*sd

# Songwon Seo, 2006, Kyunghee University

data <- subset(data, data$BD < up.outlier)
names(data)
dat <- data[,c(3,4:36)]
names(dat)

# 将自变量转化为矩阵
data1 <- data.matrix(dat[,c(2:ncol(dat))]) 
# 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
data2 <- Matrix(data1,sparse=T)  
data3 <- as.numeric(as.character(dat[,1]))
# 将自变量和因变量拼接为list
data4 <- list(data=data2,label=data3) 
# 构造模型需要的xgb.DMatrix对象，处理对象为稀疏矩阵
dtrain <- xgb.DMatrix(data = data4$data, label = data4$label) 

param <- list(max_depth = 2, eta = 0.2, silent = 1, objective = 'reg:squarederror')
nround = 150
bst = xgb.train(params = param, data = dtrain, nrounds = nround, nthread = 2)
names <- dimnames(data.matrix(dat[,c(2:ncol(dat))]))[[2]]
importance_matrix <- xgb.importance(names,model=bst)
xgb.plot.importance(importance_matrix[,],xlab = "Importance",cex=.9,main="XGBoost")

#选择6个变量参与建模

ecowater.dat <- dat[,c("BD",importance_matrix$Feature[1:8])]#XGB
names(ecowater.dat)

set.seed(65)#65#181
rf_grid <- expand.grid(mtry=seq(1:8))
rf_ctrl <- trainControl(method = "cv",number = 10,savePredictions = 'final')
rf_fit <- caret::train(BD ~ ., data = ecowater.dat,
                       method = "xgbDART",
                       # tuneGrid = rf_grid,
                       verbosity = 0,
                       trControl = rf_ctrl
)       
rf_fit
rf.dat <- data.frame(rf_fit$pred)
head(rf.dat)
#k <- data.frame(rf_fit$results)

#rf_res <- k[order(k[,"RMSE"]),][1,]

Water.lm.rf.train1<-lm(rf.dat$obs~rf.dat$pred)
Water.lm.rf.train<-lm(rf.dat$pred~rf.dat$obs)
summary(Water.lm.rf.train)


#####读取变量
# [1] "BD"     "aspect" "S1"     "PET"    "Dc"     "CW1"    "GRAV"   "slope" 
# [9] "MG"
variables <- brick(PET,MG,S1,CW1,GRAV,aspect,slope,BD,Dc)
names(variables)[1:8] <- c("PET", "MG", "S1", "CW1", "GRAV","aspect","slope", "BD","Dc")
names(variables)
rf.variables <- variables
rf.BD.predscale <- raster::predict(rf.variables,rf_fit)
rf.BD.predscale[rf.BD.predscale[] <=0 ] <- NA
plot(rf.BD.predscale)
writeRaster(rf.BD.predscale,"/home/yuanj/work/dissertation/XGBBD0.5687.tif",overwrite=T)