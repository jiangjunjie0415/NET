library(limma)
exp.rma <- read.table("Exp-STAD.txt",sep = "\t",header = T)
group <- read.table("Group.txt",sep = "\t",header = T)

class <- as.character(group$Group)
design <- model.matrix(~0+factor(class)) 
colnames(design) <- c("GC","Normal")

cont.matrix<-makeContrasts(GC-Normal,levels=design)

fit <- lmFit(exp.rma,design) 
fit <- contrasts.fit(fit, cont.matrix) 

fit <- eBayes(fit)

analysis_result <- topTable(fit,adjust='BH',number = 50000)

diff_result <- analysis_result[(abs(analysis_result$logFC) >= 1 & analysis_result$adj.P.Val < 0.05),]

write.table(analysis_result,"limma_analysis_result.txt",quote = F,sep = "\t")
write.table(diff_result,"limma_analysis_diff_result.txt",quote = F,sep = "\t")