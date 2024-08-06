library(ggplot2)
library(ggrepel)
library(ggplot2)
iris = read.table(file = "Input.txt",header = T, sep = "",row.names = 1)
iris_input <- iris 
head(iris_input) 
pca1 <- prcomp(iris_input[,-ncol(iris_input)],center = TRUE,scale. = TRUE)

df1 <- pca1$x 
df1 <- as.data.frame(df1) 
summ1 <- summary(pca1)
summ1
summ1 <- summary(pca1)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")


p.pca1 <- ggplot(data = df1,aes(x = PC1,y = PC2,color = iris_input$cluster))+
  stat_ellipse(aes(fill = iris_input$cluster),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
  geom_point(size = 1)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c("blue","red"))+
  scale_colour_manual(values = c("blue","red"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
p.pca1
