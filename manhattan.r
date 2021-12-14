library("qqman")
Args <- commandArgs()
result_file<- Args[6] #获取输入文件名字信息
result_file1<- Args[7] #获取输出文件名字信息
print(result_file)
a <- read.table(file=result_file,header=T)
png(result_file1)
manhattan(a,highlight = snpsOfInterest,col=c("#A4D3EE","#DDA0DD"))
dev.off()             