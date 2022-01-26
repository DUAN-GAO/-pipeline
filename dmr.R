.libPaths()
library(DSS)
library(bsseq)
Args <- commandArgs()
numOfArgs<-length(commandArgs())
file_names<- unlist(strsplit(Args[6],split=",")) #获取文件名字
bseq_names<- unlist(strsplit(Args[7],split=",")) #获取bseq obj命名
bseq_names<-as.vector(bseq_names)
group<- unlist(strsplit(Args[8],split="DUAN")) #获取分组信息
if(grepl(",", group[1], fixed = TRUE)){
    group1 = as.vector(unlist(strsplit(group[1],split=",")))
} else {
   group1 = group[1]
}
if(grepl(",", group[2], fixed = TRUE)){
    group2 = as.vector(unlist(strsplit(group[2],split=",")))
} else {
   group2 = group[2]
}
result_file<- Args[9] #获取结果文件名字信息
result_file1<- Args[10] #获取结果文件名字信息
list_names<-''
lis<-list()
for(i in 1:length(file_names)){ #生成dat1，dat2变量存储甲基化cov信息
    name<-paste("dat",as.character(i),sep = "")
    if(list_names==''){
        list_names<-name
    }
    else{
        list_names<-paste(list_names,name,sep = ",")
    }
    assign(name,read.table(file_names[i], header = T,col.names=c("chr","pos","N","X")))
    lis[[i]]<-get(paste("dat",as.character(i),sep = ""))
}
BSobj <- makeBSseqData(lis,bseq_names) #主分析流程
dmlResult <- DMLtest(BSobj, group1 = group1, group2 =group2, smoothing = T)
dmls = callDML(dmlResult,p.threshold = 1) #,p.threshold = 0.001,delta =  0.1
write.table(dmls,sep="\t",row.names =FALSE,file=result_file1)
dmrs = callDMR(dmlResult ,p.threshold = 0.05,delta = 0.1) 
write.table(dmrs,sep="\t",row.names =FALSE,file=result_file)
