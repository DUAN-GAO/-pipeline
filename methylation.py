# -*- coding:utf-8 -*-
#示例程序 python /public/DUAN/methylation/methylation.py /home/DUAN/methyl_rawdata/yinle/293FT_1.fq.gz,/home/DUAN/methyl_rawdata/yinle/293FT_2.fq.gz~/home/DUAN/methyl_rawdata/yinle/293FT-HOSK_1.fq.gz,/home/DUAN/methyl_rawdata/yinle/293FT-HOSK_2.fq.gz~/home/DUAN/methyl_rawdata/yinle/J-H-0_1.fq.gz,/home/DUAN/methyl_rawdata/yinle/J-H-0_2.fq.gz~/home/DUAN/methyl_rawdata/yinle/J-H-5_1.fq.gz,/home/DUAN/methyl_rawdata/yinle/J-H-5_2.fq.gz~/home/DUAN/methyl_rawdata/yinle/J-H-10_1.fq.gz,/home/DUAN/methyl_rawdata/yinle/J-H-10_2.fq.gz 293FT,293FT-HOSKDUANJ-H-0,J-H-5,J-H-10 /public/mlrna_seq/sample_result/16/ 16
import os
import sys
import time
import subprocess
import requests
from pathlib import Path
import pandas as pd
import numpy as np
import shelve
import math
files = sys.argv[1].split('~')  #第一个命令行变量，样本名称，用,分隔，如sample1,sample2
group = sys.argv[2]  #第二个命令行变量，样本的分组信息，用;分隔，如sample1,sample2DUANsample3,sample4
index = sys.argv[3] #结果路径
reid = sys.argv[4] #任务id
"""生成中间结果文件夹和最终结果文件夹"""
os.system("mkdir -p "+index)
os.system("export LD_LIBRARY_PATH=/usr/local/lib:/home/wangxiaoqi/miniconda3/lib")
start = time.time() 
path="/home/DUAN/methylation/"
os.system("mkdir -p "+path+"analysis")
os.system("mkdir -p "+path+"analysis/pre_data")
# --------------------------------------------------------------------------------
# 开始Bismark主流程,只支持双端测序
# --------------------------------------------------------------------------------
cov = []
for fil in files:
    result = fil.split(",")
    fil = result[0].split("/")[-1]
    os.system("mkdir -p "+path+"analysis/result/"+fil) #为每个样本名在result路径下生成一个文件夹
    if len(result) == 1:    #结束不符合要求的单端测序
        print("请确保是双端测序")
        exit(0)
    else:
        file1 = path+"analysis/pre_data/"+fil+"_R2_pre"
        file2 = path+"analysis/result/"+fil+"/"+fil+"_R1_pre_bismark_bt2_pe.sam"
        file3 = path+"analysis/result/"+fil+"/"+fil+"_R1_pre_bismark_bt2_pe.deduplicated.sam"
        file4 = path+"analysis/result/"+fil+".cov"
        if not Path(file1).exists(): #开始fastp流程
            os.system("fastp -w 8 -i "+result[0]+" -o "+path+"analysis/pre_data/"+fil+"_R1_pre "+"-I "+result[1]+" -O "+path+"analysis/pre_data/"+fil+"_R2_pre -h "+path+fil+".html") 
        os.system("cp "+path+fil+".html "+index)
        if not Path(file2).exists(): #开始bismark比对流程
            os.system("/public/DUAN/methylation/bismark_ref/bismark --temp_dir /home/DUAN --bowtie2 -p 4 --path_to_bowtie /public/DUAN/methylation/bismark_ref/bowtie2-2.3.4.2-linux-x86_64 -N 0 -L 20 --quiet --un --ambiguous --sam -o "+path+"analysis/result/"+fil+" /home/lvhongyi/lvhy/reference/human/ensembl_hg19/bismark_ref -1 "+path+"analysis/pre_data/"+fil+"_R1_pre"+" -2 "+path+"analysis/pre_data/"+fil+"_R2_pre")
        os.system("echo "+fil+" >> "+index+"mapping.txt") #开始将比对结果写入文件，包含覆盖度信息
        os.system("grep -w Mapping "+path+"analysis/result/"+fil+"/*_PE_report.txt >> "+index+"mapping.txt")
        if not Path(file3).exists(): #开始bismark去重流程
            os.system("/public/DUAN/methylation/deduplicate_bismark -p "+path+"analysis/result/"+fil+"/"+fil+"_R1_pre_bismark_bt2_pe.sam"+" --output_dir "+path+"analysis/result/"+fil)
        if not Path(file4).exists(): #开始从bismark结果提取甲基化信息生成cov文件流程
            os.system("python /public/DUAN/DMR/extract_CpG_data.py -i "+path+"analysis/result/"+fil+"/"+fil+"_R1_pre_bismark_bt2_pe.deduplicated.sam"+ " -o "+path+"analysis/result/"+fil+".cov")
        cov.append(fil+".cov")
# --------------------------------------------------------------------------------
# 开始PCA检测流程，只支持三个样本以上
# --------------------------------------------------------------------------------
count = 1
names = globals()
all_set = set()
pca_file = index+group+".pca.txt"
if len(cov)!=2:
    if not Path(pca_file).exists():
        for covs in cov: #
            names['n' + str(count)] = {}
            names['a' + str(count)] = set()
            with open(path+"analysis/result/"+covs) as f:
                for i in f:
                    pos = i.split("\t")[0]+i.split("\t")[1]
                    globals()['a' + str(count)].add(pos) #给set赋值，找寻公共染色体位置
                    globals()['n' + str(count)][pos] =  int(i.split("\t")[3]) / int(i.split("\t")[2]) #按照染色体位置存储每个样本的单位点甲基化值
            if count>=2:
                all_set &= globals()['a' + str(count)]
            else:
                all_set = globals()['a' + str(count)]
            count+=1
        dataframe_dic = {}
        for sett in all_set:
            dataframe_dic[sett] = [globals()['n' + str(t)][sett] for t in range(1,count)]
        df = pd.DataFrame(dataframe_dic,index=[i.split("/")[-1].split(".")[0] for i in cov])
        df = df.iloc[:,np.random.randint(1,df.shape[-1],20000)] #只选取了20000个位点
        df.to_csv(pca_file)
        os.system("Rscript /public/DUAN/methylation/pca.R "+pca_file+" "+index+"pca.png")
# --------------------------------------------------------------------------------
# 开始manhattan和按照基因十等分统计甲基化位点检测流程
# --------------------------------------------------------------------------------
cov_files = ",".join([path+"analysis/result/"+i for i in cov])
result_file_name = index+group+".txt"
result_file1_name = index+group+".dml.txt"
manhattan_file = index+group+".manhattan.txt"
slice_file = index+"slice.png"
sample_name = []
for i in group.split("DUAN"): #从用DUAN分割的组间信息名字提取样本名
    if "," in i:
        sample_name+=i.split(",")
    else:
        sample_name.append(i)
sample_name = ",".join(sample_name)
os.system("python /public/DUAN/methylation/manhattan.py "+result_file1_name+" "+manhattan_file+" "+slice_file) #test
def sigmoid_function(num): #sigmoid函数，归一化用途
    return (1-1/(1 + math.exp(-abs(num))))*0.05
if not Path(result_file_name).exists(): #开始DSS检验流程，整合cov文件，开始wald检验
    os.system("Rscript /public/DUAN/methylation/dmr.R "+cov_files+" "+sample_name+" "+group+" "+result_file_name+" "+result_file1_name)
if not Path(index+"manhattan.png").exists():
    os.system("python /public/DUAN/methylation/manhattan.py "+result_file1_name+" "+manhattan_file+" "+slice_file) #开始绘制Manhattan和10等份基因坐标同时处理
    os.system("Rscript /public/DUAN/methylation/manhattan.r "+manhattan_file+" "+index+"manhattan.png")
os.system("sed -i '"+'s/"/'+"/g' "+result_file_name)
os.system("/public/wangxiaoqi/WES_PIPLINE/softwares/bedtools intersect -a "+result_file_name+" -b /public/DUAN/DMRfinder/Annotations/hg19_RefSeq_Genes.bed -wb > "+result_file_name+".bed")
os.system("awk"+" -F "+ "' ' '"+ '{ print $1"\\t"$2"\\t"$3"\\t"$6"\t"$7"\\t"$6/$7"\\t"$9"\\t"$13}'+"' "+result_file_name+".bed"+" > "+result_file_name+".bed.genes")
os.system("cat "+result_file_name+".bed.genes |uniq > "+result_file_name+".genes")
os.system("awk '!a[$8]++' "+result_file_name+".genes > "+index+"table.txt") #结果图表保存入txt
with open(index+"table.txt") as f,open(result_file_name+".volcano.txt","w") as w:
    w.write("gene"+"\t"+"log2(Fold_change)"+"\t"+"p-value"+"\n")
    for i in f:
        genes = i.split("\t")[-1].strip("\n")
        fold = str(math.log2(float(i.split("\t")[3])/float(i.split("\t")[4])))
        p = str(sigmoid_function(float(i.split("\t")[-2])))
        w.write(genes+"\t"+fold+"\t"+p+"\n")
# --------------------------------------------------------------------------------
# 开始绘制甲基化基因的volcano plot
# --------------------------------------------------------------------------------
os.system("/usr/bin/Rscript /public/DUAN/methylation/volcano.r "+result_file_name+".volcano.txt "+index+"volcano.png")
# --------------------------------------------------------------------------------
# 开始绘制甲基化基因的heatmap plot
# --------------------------------------------------------------------------------
os.system("awk"+" -F "+ "' ' '"+ '{ print $8"\\t"$4"\\t"$5}'+"' "+index+"table.txt"+" > "+result_file_name+".heatmap")
os.system("sed -i '"+'s/基因名称/'+"/g' "+result_file_name+".heatmap")
os.system("sed -i '1i\基因\\t实验组该区域平均甲基化水平\\t对照组该区域平均甲基化水平' "+result_file_name+".heatmap") #处理后的矩阵文件
res1 = subprocess.call("/usr/bin/Rscript /public/DUAN/methylation/plot.r "+result_file_name+".heatmap "+index,shell=True) #绘制heatmap图片
# --------------------------------------------------------------------------------
# 程序结束后发送状态代码到API
# --------------------------------------------------------------------------------
if res1 == 0:
    r = requests.get("http://192.168.4.203:81/api/methylcell/methylcellresult?id="+reid+"&status=3") #成功
else:
    r = requests.get("http://192.168.4.203:81/api/methylcell/methylcellresult?id="+reid+"&status=4") #失败
end = time.time()


