library(maftools)
laml = read.maf(maf = 'TCGA.STAD.somatic.maf',isTCGA = TRUE)
laml

#Shows sample summry.
getSampleSummary(laml)
SampleSummary = getSampleSummary(laml)
write.table (SampleSummary, file ="Result1_SampleSummary.txt")

#Shows gene summary.
getGeneSummary(laml)
GeneSummary = getGeneSummary(laml)
write.table (GeneSummary, file ="Result2_GeneSummary.txt")

#Shows all fields in MAF
getFields(laml)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'mean', dashboard = TRUE, titvRaw = FALSE)

#Shows sample summry.
getSampleSummary(laml)
MatrixSample2 = getSampleSummary(laml)
write.table (MatrixSample2, file ="Result3_MatrixSample2.txt")

#oncoplot for top ten mutated genes.
Gene = read.table(file = "genelist.txt")
Gene1 = as.matrix(Gene)
oncoplot(maf = laml, top = 10, fontSize = 0.6,genes = Gene1,cBioPortal = F,removeNonMutated = T )
A = oncoplot(maf = laml, fontSize = 0.2, cBioPortal = F,top = 10, genes = Gene1,removeNonMutated = T )

#plot titv summary
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
write.table (laml.titv, file ="TiTv.txt")