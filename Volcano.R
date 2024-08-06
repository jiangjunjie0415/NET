library(RColorBrewer)
library(ggplot2)
diff <- read.csv("Volcano_Input.csv",header = T)
logFC <-diff$logFC
padj <- diff$adj.P.Val
data <- data.frame(logFC=logFC,padj=padj)

data$sig[(data$padj > 0.05|data$padj=="NA")|(data$logFC < 0.584963)& data$logFC > -0.584963] <- "no"
data$sig[data$padj <= 0.05 & data$logFC >= 0.584963] <- "up"
data$sig[data$padj <= 0.05 & data$logFC <= -0.584963] <- "down"
x_lim <- max(logFC,-logFC)
theme_set(theme_bw())
p <- ggplot(data,aes(logFC,-1*log10(padj),
                     color = sig))+geom_point()+xlim(-4,4)+ labs(x="log2(FoldChange)",y="-log10(FDR)")
p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-0.5,0.5),linetype=4)
p <- p +theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0))
p <- p  +guides(colour = FALSE)
p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
p
