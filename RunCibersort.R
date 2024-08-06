setwd("G:/ZJU/实验室/18.0 Wu文章/5.3 免疫细胞浸润分析/CIBERSORT")
source('Cibersort.R')

result1 <- CIBERSORT('LM22.txt','ComBat_data.txt', perm = 1000, QN = T)
