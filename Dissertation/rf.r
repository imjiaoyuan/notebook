library(ggplot2)
library(caret)
library(raster)
library(httpgd)
Ra <- read.csv("/home/yuanj/work/dissertation/RESULTS.csv")
### remove outliers by mean +- 2*sd
mean(Ra$BD);min(Ra$BD);max(Ra$BD)
mean <- mean(Ra$BD)
sd <- sd(Ra$BD)

up.outlier <- mean + 1.5*sd
low.outlier <- mean - 1.5*sd

# Songwon Seo, 2006, Kyunghee University

Ra.subset <- subset(Ra, BD < up.outlier)
nrow(Ra.subset)
names(Ra.subset)

subet<-Ra.subset[,c(3,4:36)]

names(subet)

set.seed(123)

#333，123(best)， 

subsets <- c(1:ncol(subet))
subsets
ctrl.rfe <- rfeControl(functions = rfFuncs,
                       method = "cv",
                       repeats = 10,
                       verbose = FALSE)

lmProfile <- rfe(x = subet[,-1], y = subet[,1],
                 sizes = 10,
                 rfeControl = ctrl.rfe)

lmProfile

lmProfile$optVariables
train.dat <- subset(subet,select = c("BD",c(paste0(lmProfile$optVariables,sep = ''))))
library(caretEnsemble)
library(randomForest)
set.seed(123)
# Stacking Algorithms - Run multiple algos in one call.
tr_control <- trainControl(method = "cv", 
                           number = 10, 
                           search = "random",
                           savePredictions = "final"
)

rf_grid <- expand.grid(mtry = c(2, 3, 4, 5,6)#,
                       #splitrule = c("gini", "extratrees"),
                       #min.node.size = c(1, 3, 5)
)
rf_grid

#algorithmList <- c('rf', 'treebag','bstTree', 'gbm')

set.seed(123)

models <- train(BD~ ., data = train.dat, metric ="Rsquared",
                maximize = TRUE,
                tuneLength = 15,
                #preProcess = c("center","scale"),
                trControl = tr_control, 
                tunegrid = rf_grid,
                method = "rf", 
                importance = T)

models$dots$importance

varImp(models)

names(models)

plot(as.data.frame(importance(models$finalModel))[,1])

imp.df <- data.frame(c(paste0(lmProfile$optVariables,sep = '')),
                     importance(models$finalModel)[,1])

imp.df.sort <- imp.df[order(imp.df[,2],decreasing=TRUE),] 

imp.df.sort <- imp.df[with(imp.df, order(imp.df[,2], decreasing = T)),]


# 调整图形边距
par(mar = c(9, 5,4, 2) + 0.1)

#pdf("E:/BigData/CarbonAllocation/Results/variable_importance.pdf", height = 4, width = 4)

barplot(imp.df.sort[,2],names.arg = imp.df.sort[,1], ylim = c(0, 5), 
        ylab = "%IncMSE", xlab = "", las = 2)
imp.df.sort[,2]

train.dat <- subset(subet,select = c("BD",c(paste0(lmProfile$optVariables,sep = ''))))

set.seed(123)


# train the model using rf = random forest.
tr_control <- trainControl(method = "cv", 
                           number = 10, 
                           search = "random",
                           savePredictions = "final"
)

rf_grid <- expand.grid(mtry = c(2, 3, 4, 5,6)#,
                       #splitrule = c("gini", "extratrees"),
                       #min.node.size = c(1, 3, 5)
)
rf_grid

set.seed(123)
models <- train(BD ~ ., data = train.dat, metric ="Rsquared",
                maximize = TRUE,
                tuneLength = 15,
                preProcess = c("center","scale"),
                trControl = tr_control, 
                tunegrid = rf_grid,
                method = "rf", 
                importance = T)
models

# plot correlation between observed and predicted
cv.dat.rf <- models$pred # models$glm_ensemble$pred

BD.lm.rf.train <- lm(cv.dat.rf$obs ~ cv.dat.rf$pred)
summary(BD.lm.rf.train)

rmse_value <- RMSE(cv.dat.rf$pred, cv.dat.rf$obs)
print(rmse_value)


# 假设前面的代码已经运行，cv.dat.rf 等相关数据已经存在

# 加载 ggplot2包
library(ggplot2)

# 创建数据框用于绘图
plot_df <- data.frame(obs = cv.dat.rf$obs, pred = cv.dat.rf$pred)

# 拟合一条线性回归直线（用于展示拟合情况）
lm_fit <- lm(pred ~ obs, data = plot_df)

# 计算R²值
r_squared <- summary(lm_fit)$r.squared

# 绘制散点图展示观测值和预测值的关系
p <- ggplot(plot_df, aes(x = obs, y = pred)) +
  geom_point() +
  labs(x = "观测值", y = "预测值", title = "物种丰富度模型精度") +
  theme(plot.title = element_text(hjust = 0.5))

# 添加1:1线
p <- p + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")

# 添加拟合直线
p <- p + geom_abline(intercept = coef(lm_fit)[1], slope = coef(lm_fit)[2], linetype = "solid", color = "blue")

# 在右下角添加R²和RMSE值的标注
p <- p + annotate("text", x = max(plot_df$obs), y = min(plot_df$pred), 
                  label = paste("R²:", round(r_squared, 2), "\nRMSE:", round(rmse_value, 2)), hjust = 1, vjust = 0)

# 显示图形
print(p)

#####读取变量
aspect <- raster("/home/yuanj/work/dissertation/模拟用数据/aspect.tif")
S1 <- raster("/home/yuanj/work/dissertation/模拟用数据/S1.tif")
H <- raster("/home/yuanj/work/dissertation/模拟用数据/H.tif")
SW1 <- raster("/home/yuanj/work/dissertation/模拟用数据/SW1.tif")
PET <- raster("/home/yuanj/work/dissertation/模拟用数据/PET.tif")
Dc <- raster("/home/yuanj/work/dissertation/模拟用数据/Dc.tif")
GRAV <- raster("/home/yuanj/work/dissertation/模拟用数据/GRAV.tif")
DEM <- raster("/home/yuanj/work/dissertation/模拟用数据/DEM.tif")
AK <- raster("/home/yuanj/work/dissertation/模拟用数据/AK.tif")
N <- raster("/home/yuanj/work/dissertation/模拟用数据/NA.tif")

# CL <- raster("/home/yuanj/work/dissertation/模拟用数据/CL.tif")
# LNUM <- raster("/home/yuanj/work/dissertation/模拟用数据/LNUM.tif")
# BD <- raster("/home/yuanj/work/dissertation/模拟用数据/BD.tif")
# LDEP <- raster("/home/yuanj/work/dissertation/模拟用数据/LDEP.tif")
# AN <- raster("/home/yuanj/work/dissertation/模拟用数据/AN.tif")
# CA <- raster("/home/yuanj/work/dissertation/模拟用数据/CA.tif")
# PH <- raster("/home/yuanj/work/dissertation/模拟用数据/PH.tif")
# GRAV <- raster("/home/yuanj/work/dissertation/模拟用数据/GRAV.tif")
# K <- raster("/home/yuanj/work/dissertation/模拟用数据/K.tif")
# PDEP <- raster("/home/yuanj/work/dissertation/模拟用数据/PDEP.tif")
# NDVI <- raster("/home/yuanj/work/dissertation/模拟用数据/NDVI.tif")
# MG <- raster("/home/yuanj/work/dissertation/模拟用数据/MG.tif")
# CEC <- raster("/home/yuanj/work/dissertation/模拟用数据/CEC.tif")
# CW1 <- raster("/home/yuanj/work/dissertation/模拟用数据/CW1.tif")
# AP <- raster("/home/yuanj/work/dissertation/模拟用数据/AP.tif")
qw<-brick(aspect,S1,H,SW1,PET,Dc,GRAV,DEM,AK,N)

names(qw)<-lmProfile$optVariables
qw[is.na(qw[])==T]<-0
pred<-predict(qw,models)
par(mar = c(4, 4, 2, 2))
plot(pred, cex.axis = 0.8, cex.lab = 0.8)
writeRaster(pred,"/home/yuanj/work/dissertation/BD.tif")
