# -*- coding:utf-8 -*-
import os
import sys
import shelve
import matplotlib.pyplot as plt
import numpy as np
in_file = sys.argv[1]
out_file = sys.argv[2]
plot_file = sys.argv[3]
help_dict = {"-2kb":0,"10%":0,"20%":0,"30%":0,"40%":0,"50%":0,"60%":0,"70%":0,"80%":0,"90%":0,"100%":0,"+2kb":0}
with shelve.open("/home/DUAN/methylation/data_files/dic") as dbs:
    new_dic = dbs["dic"]

with open(in_file) as f,open(out_file,"w") as w:
    w.write("SNP"+"\t"+"CHR"+"\t"+"BP"+"\t"+"P"+"\n")
    counter = 0
    dic = {"chr1":1,"chr2":2,"chr3":3,"chr4":4,"chr5":5,"chr6":6,"chr7":7,"chr8":8,"chr9":9,"chr10":10,"chr11":11,"chr12":12,"chr13":13,"chr14":14,"chr15":15,"chr16":16,"chr17":17,"chr18":18,"chr19":19,"chr20":20,"chr21":21,"chr22":22,"chrX":23,"chrY":24}
    for i in f:
        if i.split("\t")[0].strip('"') == 'chr':
            pass
        else:
            if i.split("\t")[0].strip('"') not in list(dic.keys()):
                continue
            else:
                SNP = str(counter)
                chr_ori = i.split("\t")[0].strip('"')

                CHR = str(dic[i.split("\t")[0].strip('"')])
                BP = i.split("\t")[1]
                P = i.split("\t")[9]
                for pos in new_dic[chr_ori]:
                    if int(BP) in pos:
                        v = (int(BP) - list(pos)[0]) / len(list(pos))
                        if 0<=v<0.1:
                            help_dict["10%"] +=1
                        elif 0.1<=v<0.2:
                            help_dict["20%"] +=1
                        elif 0.2<=v<0.3:
                            help_dict["30%"] +=1
                        elif 0.3<=v<0.4:
                            help_dict["40%"] +=1
                        elif 0.4<=v<0.5:
                            help_dict["50%"] +=1
                        elif 0.5<=v<0.6:
                            help_dict["60%"] +=1
                        elif 0.6<=v<0.7:
                            help_dict["70%"] +=1
                        elif 0.7<=v<0.8:
                            help_dict["80%"] +=1
                        elif 0.8<=v<0.9:
                            help_dict["90%"] +=1
                        else:
                            help_dict["100%"] +=1
                    break


                if float(P) < 1e-30:
                    P = str(1e-30)
                counter+=1
                w.write(SNP+"\t"+CHR+"\t"+BP+"\t"+P+"\n")

y = list(help_dict.values())[1:11]
x = list(help_dict.keys())[1:11]
plt.plot(x,y)
plt.xlabel('gene position sliced to ten equal parts')
plt.ylabel('methylation counts in certain area')
my_y_ticks = np.arange(0, 10000, 50)
plt.yticks(my_y_ticks)
plt.show()
plt.savefig(plot_file)