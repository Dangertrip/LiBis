def mates(dic,file):
    #chr start front_del end_del
    result=[]
    for name in dic:
        if name[-2:]=='_0':
            _name = name[:-2]
            if _name[-1]=='1':
                continue
            if (_name[:-1]+'1'+'_0' in dic) or (_name[:-1]+'1'+'_1' in dic):
                continue
            mate = _name +'_1'
            if mate in dic:
                a = '\t'.join(dic[name])
                b = '\t'.join(dic[mate])
                s = _name+'\t'+a+'\t'+b+'\n'
                result.append(s)
    with open(file+'.mate','w') as f:
        f.writelines(result)

def duplicates(dic,file):
    result=[]
    for name in dic:
        if name[-3]=='0':
            dupname = name[:-3]+'1'+name[-2:]
            if dupname in dic:
                a = '\t'.join(dic[name])
                b = '\t'.join(dic[dupname])
                s = name+'\t'+a+'\t'+b+'\n'
                result.append(s)
    with open(file+'.dup','w') as f:
        f.writelines(result)

def single(dic,file):
    result=[]
    for name in dic:
        _name = name[:-4]
        tail=['_0_0','_0_1','_1_1','_1_0']
        sum=0
        for t in tail:
            if _name+t in dic:
                sum+=1
        if sum!=1: continue
        s = name+'\t'+'\t'.join(dic[name])+'\n'
        result.append(s)
    with open(file+'.single','w') as f:
        f.writelines(result)

def readfile(file):
    dic={}
    with open(file) as f:
        for line in f:
            temp = line.strip().split()
            dic[temp[0]] = [temp[1],temp[2],temp[3],temp[4]]
    return dic

if __name__=='__main__':
    ff = 'M6G_split.bam.pos.comb'
    dic=readfile(ff)
    mates(dic,ff)
    duplicates(dic,ff)
    single(dic,ff)
