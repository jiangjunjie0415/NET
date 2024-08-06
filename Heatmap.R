library(circlize)
library(ComplexHeatmap)
Data= read.table("Input.txt",header = T, row.names = 1)
Data = as.matrix(Data)
z=t(Data)
z2 = scale(z, center = T, scale = T)
Matrix = t(z2)
Matrix = as.matrix(Matrix)
range(Matrix)

colgroup=read.table("Annoation.txt",sep="\t",header=TRUE,row.names=1)
Group2= subset(colgroup, select = Group)
Group2= t(Group2)
Group2= as.numeric(Group2)


col2 = list(Group2 = c("0"="#006dbb","1"="#e60012"))
class(col2)
df2 = data.frame(Group2)
top_annotation =  HeatmapAnnotation(df = df2, col = col2)

m = Heatmap(Data,name = " ",
            top_annotation = top_annotation,
            show_heatmap_legend = T,
            border = F,
            show_column_names = F,
            show_row_names = F,
            column_title = NULL,
            cluster_columns = F,
            cluster_rows = T)
m

