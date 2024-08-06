library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
library(stringr)
library(enrichplot)

gene = read.table("Genelist-GO_KEGG.txt", header = TRUE)
geneList<-gene$logFC 
names(geneList)=gene$ID
geneList=sort(geneList,decreasing = T)

A=gene[,1]

ego_BP <- enrichGO(gene = gene[,1],
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
BP <- barplot(ego_BP, showCategory=10,title="BP")

ego_CC <- enrichGO(gene = gene[,1],
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
CC <- barplot(ego_CC, showCategory=10,title="CC")

ego_MF <- enrichGO(gene = gene[,1],
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
MF <- barplot(ego_MF, showCategory=10,title="MF")

H = BP+CC+MF
H

B = read.table("KEGG.txt", header = TRUE)
A=B[,1]
CaseGeneSet=bitr(A, 'SYMBOL', "ENTREZID", "org.Hs.eg.db",drop =F)[, "ENTREZID"]
CaseGeneSet=as.data.frame(CaseGeneSet)
D = CaseGeneSet[,1]
kk <- enrichKEGG(gene = D,
                 organism = 'hsa', 
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05)
dotplot(kk,color = "p.adjust",showCategory=10) 

GMT<-read.gmt("REACTOME.gmt") 
B<-GSEA(geneList,TERM2GENE = GMT) 

paths <- c("REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM","REACTOME_PROGRAMMED_CELL_DEATH","REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING",
           "REACTOME_TCR_SIGNALING","REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
           "REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING","REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
           "REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY","REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION","REACTOME_DOWNSTREAM_TCR_SIGNALING")#ѡȡ????Ҫչʾ??ͨ·ID
gseaplot2(B,paths,pvalue_table = T) 