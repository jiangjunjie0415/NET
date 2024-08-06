library(ggplot2)
require(tidyverse)
require(ggplot2)
require(cowplot)
library(ggpubr)
library(reshape2)
Expr <- read.table("ESTIMATE.txt",header = T)
A= melt(Expr,
        id.vars = c('ID','Group'),
        variable.name='score',
        value.name='abundance')
p <-ggplot(data = A, aes(x=score,y=abundance))+geom_boxplot(aes(fill=Group))
q <-p+ facet_wrap(~ score, scales="free")
q+ stat_compare_means(aes(group = Group), label = "p.signif",method="wilcox.test")


p <- ggboxplot(A, x = "Immunecell", y = "abundance",
               color = "cluster")
p + stat_compare_means(aes(group = cluster), label = "p.signif",method="wilcox.test")


