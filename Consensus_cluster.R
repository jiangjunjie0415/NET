library(readxl)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(ggsci)
library(ggplot2)
library(pheatmap) 
library(RColorBrewer)
library(ConsensusClusterPlus)
library(survival)
library(survminer)
library(export)

input_path <- "D:/科研/5.0 WHY-胃癌NET分析/2.0 一致性聚类/"
fig_path <- "D:/科研/5.0 WHY-胃癌NET分析/2.0 一致性聚类/"

DRG <- read.csv(file = paste0(input_path,"Model1_Matrix.csv"), 
                header =TRUE,
                sep = ",",
                row.names = 1,
                check.names = F)
DRG <- as.matrix(DRG)

ConsensusClusterPlus(DRG, 
                     maxK = 10, 
                     reps = 1000, 
                     pItem = 0.8, 
                     pFeature = 1,
                     clusterAlg = "km",
                     corUse = "pairwise.complete.obs",
                     distance='euclidean',
                     title="DRG_euclidean",
                     seed=123, 
                     plot="pdf", 
                     writeTable=T)



