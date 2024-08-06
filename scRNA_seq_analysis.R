library(dplyr)
library(Seurat)
library(data.table)
library("R.utils")
library(SingleR)
library(matrixStats)
library(celldex)
library(AUCell) 
library(clusterProfiler) 
library(ggplot2)
library(monocle)

load("pbmc.Rdata")
m1_26=merge(x = m1_16, y = m17_26, merge.data = TRUE)
save(m1_26,file = "m1_26.Rdata")

sample39 <- fread("GSM5573504_sample39.csv.gz",
                  data.table = F)
d1=sample39[,-1]
rownames(d1)=sample39[,1]
pbmc1 <- CreateSeuratObject(counts = d1,
                            min.cells = 3, 
                            min.features = 200)

sample40 <- fread("GSM5573505_sample40.csv.gz",
                  data.table = F)
d2=sample40[,-1]
rownames(d2)=sample40[,1]
pbmc2 <- CreateSeuratObject(counts = d2,
                            min.cells = 3, 
                            min.features = 200)
m = merge(x = pbmc1, y = pbmc2, add.cell.ids = c("sample39", "sample40"), merge.data = TRUE)
save(m,file = "m25_26.Rdata")
##AnalyzeData
pbmc = m1_26
rm(m1_26)
gc()
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") 
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)   
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)      
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10) 
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt") 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5) 
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") 
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) 
DimHeatmap(pbmc, dims = 1:5, cells = 500, balanced = TRUE) 
pbmc <- JackStraw(pbmc, num.replicate = 100) 
pbmc <- ScoreJackStraw(pbmc, dims = 1:20) 
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc,ndims = 50)
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.5) 
head(Idents(pbmc), 5) 
pbmc <- RunUMAP(pbmc, dims = 1:19) 
DimPlot(pbmc, reduction = "umap", label = TRUE)
saveRDS(pbmc,file = 'ClusterResult.rds')
save(pbmc,file = "pbmc.Rdata")
save(pbmc.markers,file = "pbmc-markers.Rdata")
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
library(dplyr)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10,file = "top10.txt",sep = "\t")
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
VlnPlot(pbmc, features = top10$gene[1:12],pt.size=0)
saveRDS(pbmc, file = "pbmc.rds")
save.image(file = 'pbmc.RData')
bfreaname.pbmc <- pbmc
pbmc <- bfreaname.pbmc 
new.cluster.ids <- c("Tcell","Tcell", "Plasma_cell", "Plasma_cell", "Epithelial_cell", "Macrophage","Plasma_cell", "DC", 
                     "Macrophage", "Bcell","Fibroblast", "Endothelial_cell","Mast_cell","Fibroblast","Epithelial_cell","Fibroblast","Epithelial_cell","Plasma_cell","Endocrine_cell")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) 

write.table(pbmc_tumor@meta.data,"metadata_Primary_GC_tumor.txt")
tissuetype = read.table("Input_Metadata_tissuetype.txt")
tissuetype = as.vector(tissuetype)
pbmc <- AddMetaData(object = pbmc,               
                    metadata = tissuetype,           
                    col.name = "tissue") 
DimPlot(pbmc, reduction = "umap", group.by = "tissue",label = TRUE)

pbmc_tumor <- subset(x = pbmc, subset = tissue == "Primary_GC_tumor")
Lauren <-  read.table("Input_Metadata_Lauren.txt")
Lauren = as.vector(Lauren)
pbmc_tumor <- AddMetaData(object = pbmc_tumor,               
                          metadata = Lauren,             
                          col.name = "Lauren") 
DimPlot(pbmc_tumor, reduction = "umap", group.by = "Lauren",label = TRUE)
VlnPlot(pbmc_tumor, features = c("PRF1"),pt.size=0)
pbmc_tumor$tissue
{
  ##AUcell
  cells_rankings <- AUCell_buildRankings(pbmc_tumor@assays$RNA@data, splitByBlocks = T) 
  Hallmarker <- read.gmt("SAM.gmt") 
  geneSets <- lapply(unique(Hallmarker$term), function(x){print(x);Hallmarker$gene[Hallmarker$term == x]})
  names(geneSets) <- unique(Hallmarker$term)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
  
  geneSet <- "SAM"
  aucs <- as.numeric(cells_AUC@assays@data@listData$AUC)
  write.table(aucs,"Result_auc.txt")
  pbmc_tumor$AUC  <- aucs
  
  VlnPlot(pbmc_tumor, features = c("AUC"),pt.size=0)
  library(ggraph)
  ggplot(data.frame(pbmc_tumor@meta.data, pbmc_tumor@reductions$umap@cell.embeddings), aes(UMAP_1, UMAP_2, color=AUC)
  ) + geom_point(size=0.5
  ) + scale_color_viridis(option="A")  + theme_light(base_size = 15)+labs(title = "SAM")+
    theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
    theme(plot.title = element_text(hjust = 0.5))
  min(pbmc_tumor@meta.data$AUC)
  table(Idents(pbmc_tumor))
  pbmc_tumor$seurat_annoation = pbmc_tumor@active.ident
  pbmc_tumor$seurat_annoation
  pbmc_tumor$AUC_group <- cut(pbmc_tumor@meta.data$AUC, breaks = c(0, 0.0499432603118216, 0.21388395500296), labels = c("Low","high"), right=FALSE)
  ggplot(pbmc_tumor@meta.data,aes(pbmc_tumor$AUC_group,fill = pbmc_tumor$seurat_annoation))+
    geom_bar(position = "fill", alpha = 0.9)+
    scale_fill_brewer(palette = "Set1")+
    theme_classic()
} ##AUC

pbmc_cancercell <- subset(x = pbmc_tumor, subset =  seurat_annoation == "Epithelial_cell")
save(pbmc_cancercell,file = "pbmc_cancercell.RData")
expr_matrix <- as(as.matrix(pbmc_cancercell@assays$RNA@counts),'sparseMatrix')
p_data <- pbmc_cancercell@meta.data
p_data$seurat_annoation <- pbmc_cancercell@active.ident
f_data <- data.frame(gene_short_name = row.names(pbmc_cancercell),row.names = row.names(pbmc_cancercell))
pd <- new('AnnotatedDataFrame',data = p_data)
fd <- new('AnnotatedDataFrame',data = f_data)
cds <- newCellDataSet(expr_matrix, phenoData = pd,featureData = fd, lowerDetectionLimit = 0.1,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
cds <- orderCells(cds)
