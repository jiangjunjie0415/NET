library(ggplot2)
library(lattice)
library(caret)
setwd("D:/3.1 科研课题/5.0 WHY-胃癌NET分析/修回2-Heliyon/Figure S_Cross_validation")
data(iris)

# 创建模型
irisData <- iris[1:60,]
irisData$Species <- as.factor(as.character(irisData$Species))
irisData2 <- read.csv("NET_Model4.csv",sep=',',header=TRUE)
irisData2$Response <- as.factor(as.character(irisData2$Response))

library('caret')

fitControl <- trainControl(
  method = 'LOOCV',                
  number = 1,                     
  savePredictions = 'final',        
  classProbs = T ,
  seed = as.list(rep(1,46)),                
  summaryFunction=twoClassSummary 
) 
names(getModelInfo())
Model <- train(Response ~ ACTA2 + AKT1 + HIST2H2AC + ITGB3 + PRKCG, data=irisData2 ,method='rf',   
      tuneGrid=data.frame(mtry=2)  ,trControl = fitControl)
print(Model)


fitControl2 <- trainControl(
  method = 'boot',                
  number = 10,                     
  savePredictions = 'final', 
  p = 0.75,
  classProbs = T ,
  seed = as.list(rep(1,46)),                
  summaryFunction=twoClassSummary 
) 
names(getModelInfo())
Model2 <- train(Response ~ ACTA2 + AKT1 + HIST2H2AC + ITGB3 + PRKCG, data=irisData2 ,method='rf',   
               tuneGrid=data.frame(mtry=2)  ,trControl = fitControl2)
print(Model2)
