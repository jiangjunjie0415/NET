library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
library(stringr)

gene = read.table("Genelist.txt", header = TRUE)
geneList<-gene$logFC 
names(geneList)=gene$ID 
geneList=sort(geneList,decreasing = T) 

GMT<-read.gmt("KEGG+REACTOME.gmt") 
B<-GSEA(geneList,TERM2GENE = GMT)
library(enrichplot)
B@result[["ID"]]
paths <- c("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION")
gseaplot2(B,paths,pvalue_table = T) 

