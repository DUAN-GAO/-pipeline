library("FactoMineR")
library("factoextra")
Args <- commandArgs()
result_file<- Args[6] #获取输入文件名字信息
result_file1<- Args[7] #获取输出文件名字信息
a <- read.table(file=result_file,header=T,row.names=1,sep=",")
res.pca <- PCA(a, graph = FALSE)
png(result_file1)
fviz_pca_ind(res.pca, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
             )
dev.off()             