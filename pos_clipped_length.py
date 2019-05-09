def comb(bam_file,fq_start,fq_end):
    fq={}
    with open(fq_start) as f:
        lines = f.readlines()
    for i in range(1,len(lines),4):
        l = len(lines[i].strip())
        name = lines[i-1].strip()[1:]
        fq[name]=[str(l)]
    with open(fq_end) as f:
        lines = f.readlines()
    for i in range(1,len(lines),4):
        l = len(lines[i].strip())
        name = lines[i-1].strip()[1:]
        fq[name].append(str(l))
    with open(bam_file) as f:
        lines = f.readlines()
    result=[]
    for line in lines:
        temp = line.strip().split()
        if temp[0] in fq:
            result.append(temp[0]+'\t'+temp[1]+'\t'+temp[2]+'\t'+fq[temp[0]][0]+'\t'+fq[temp[0]][1]+'\n')
    with open(bam_file+'.comb','w') as f:
        f.writelines(result)




if __name__=="__main__":
    comb("M6G_split.bam.pos","M6G_clipped_start.fastq","M6G_clipped_end.fastq")
