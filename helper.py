import sys
import pymysql
import re
in_file = sys.argv[1]
out_file = sys.argv[2]
with open(in_file) as f,open(out_file,"w") as w: #refseq gene id转换到gene name
    conn = pymysql.connect(host="genome-mysql.cse.ucsc.edu", user="genome",database="hg19")
    cursor = conn.cursor()
    for i in f:
        en_id = i.split("\t")[0]
        left = i.split("\t")[1:]
        sql = "select distinct G.gene,N.value from ensGtp as G, ensemblToGeneName as N where G.transcript=N.name and G.gene = '"+en_id+"' ;"
        cursor.execute(sql)
        ret = cursor.fetchone()
        if len(en_id)==0:
            w.write(i)
        elif ret == None:
            pass
        else:
            w.write(ret[1]+"\t"+"\t".join(left))



        
