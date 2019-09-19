import numpy as np
import math
import random
import sys

#==============================Variable setting===================================
readsname = sys.argv[1]
conversion_ratio=1
methylation_ratio=0.75
readsnum=5000000
readlength=80
#================================================================================


#==============================Get Genome fasta dictory=========================
chr=''
ans={}
pos=0
mark=False
s=''
length=0
l=[]
chr_order=[]
with open('/data/dsun/ref/humanigenome/hg19.fa') as f:
    for line in f:
        temp = line.strip()
        if temp[0]=='>':
            if chr!='':
                ans[chr]=s
                ll=len(s)
                length+=ll
                l.append(ll)
                s=''
            chr= temp[1:]
            chr_order.append(chr)
            continue
        s+=temp
ans[chr]=s
ll=len(s)
l.append(ll)
length+=ll



#print(ans['chr3'][43483166:43483246].upper())
#print(ans['chr3'][43483379:43483459].upper())
#sys.exit()
#=========================================================================

'''
def getcg(read):
    num=0
    if random.random()>0.5:
        for i in range(1,len(read)):
            if ('G'==read[i] or 'g'==read[i]) and ('C'==read[i-1] or 'c'==read[i-1]):
                num+=1
    else:
        for i in range(1,len(read)):
            if ('C'==read[i] or 'c'==read[i]) and ('G'==read[i-1] or 'g'==read[i-1]):
                num+=1

    return num
'''

def fake(fl,read):
    base=['A','T','C','G']
    l = len(read)
    if l>=fl:
        return read,0,0
    front = int((fl - l)/2)
    end = fl-l-front
    r=read
    for i in range(front):
        pos = random.randint(0,3)
        r=base[pos]+r
    for i in range(end):
        pos = random.randint(0,3)
        r=r+base[pos]
    #print(r,front,end,l)
    return r,front,end

def reverse(read):
    dic={'A':'T','T':'A','C':'G','G':'C','N':'N'}
    r=''
#    print('r')
    for rr in read:
        r=dic[rr.upper()]+r
    return r

def bisulfite(conv_r,meth_r,read):
    r=''
    l = len(read)
    for i in range(1,l-1):
        base=read[i].upper()
        if base=='C':
            c = random.random()
            if c<conv_r:
                if read[i+1].upper()=='G':
                    m = random.random()
                    if m>meth_r:
                        base='T'
                else:
                    base='T'
        r=r+base
    return r

#time1=time.time()
num=[]
zero=0
finalbed=[]
for i in range(readsnum):
    pos = random.randint(0,length)
    chr=0
    f=0
    t=0
    while pos>l[chr]:
        pos-=l[chr]
        chr+=1
    #readlength=length
    start=pos-1
    end=pos+readlength+1
    if pos<1:continue
    if pos>l[chr]:continue
    read=ans[chr_order[chr]][start:end]
    r=read
    a1=random.random()
    a2=random.random()
    if a1>0.5:
        r=reverse(read)#Get reads from +/- strand
    r = bisulfite(conversion_ratio,methylation_ratio,r)
    if a2>0.5:
        r=reverse(r)#PCR +/-
    #r,f,t = fake(fake_length,r)
    quality = 'E'*readlength
#    print(readlength)
    #print(a1,a2)
    #print(read)
    print('@'+str(i)+'_'+readsname)
    print(r)
    print('+')
    print(quality)
    finalbed.append(chr_order[chr]+'\t'+str(start+1)+'\t'+str(end-1)+'\t'+str(i)+'_'+readsname+'\t'+str(f)+'\t'+str(t)+'\n')
with open(readsname+'_simulation.bed','w') as f:
    f.writelines(finalbed)

#print(time.time()-time1)
    #g=getcg(read)
    #num.append(g)
    #if g==0: zero+=1
#print(zero/float(100000000))
#print((100000000-zero)/float(100000000))

#import seaborn as sns
#from matplotlib import pyplot as plt
#import numpy as np
#n=np.array_split(num,10)
#for nn in n:
#    sns.distplot(nn)
#plt.savefig('125.pdf')


##
