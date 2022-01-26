Args <- commandArgs()
dat<-read.table(Args[6],sep = "\t",row.names = 1,header = T)
pheatmap::pheatmap(dat[1:50,],filename = paste(Args[7],"/heatmap.png",sep=""))