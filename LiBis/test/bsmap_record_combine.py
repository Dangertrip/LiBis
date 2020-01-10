import sys,os
import re
import pandas as pd

def record_process(files):
    ans=[]
    title=['fastq_file_name','output_file','total_read_pairs','aligned_pairs','unique_pairs','non-unique_pairs','#1_unpaired_read','#1_unique_reads','#1_non-unique_reads','#2_unpaired_read','#2_unique_reads','#2_non-unique_reads']
    for file in files:
        if not '.record' in file: continue
        with open(file) as f:
            lines = f.readlines()
        #print(file)
        a=[]
        a.append(lines[4][lines[4].find(':')+2:lines[4].rfind('(')-1].strip())
        a.append(lines[6][lines[6].find(':')+2:lines[6].rfind('(')-1].strip())
        a.append(lines[7].split()[-6].strip())
        for i in range(8,11):
            temp=re.split(':|,',lines[i])
            a.append(temp[1].strip())
            a.append(temp[3].strip())
            a.append(temp[5].strip())
        ans.append(a)
    data = pd.DataFrame(ans,columns=title)
    return data


if __name__=="__main__":
    import os
    files = os.listdir()
    print(record_process(files))

